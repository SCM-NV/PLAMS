import itertools
import os
import re
import shutil
import threading
from collections import defaultdict
from functools import cached_property
from genericpath import isdir, isfile
from os.path import join as opj
from typing import TYPE_CHECKING, Callable, Dict, List, Optional

from scm.plams.core.basejob import MultiJob
from scm.plams.core.enums import JobStatus
from scm.plams.core.errors import FileError, PlamsError
from scm.plams.core.formatters import JobCSVFormatter
from scm.plams.core.functions import config, get_logger, log
from scm.plams.core.logging import Logger

if TYPE_CHECKING:
    from scm.plams.core.basejob import Job
    from scm.plams.core.settings import Settings

__all__ = ["JobManager"]


class JobManager:
    """Class responsible for jobs and files management.

    Every instance has the following attributes:

    *   ``foldername`` -- the working folder name.
    *   ``workdir`` -- the absolute path to the working folder.
    *   ``logfile`` -- the absolute path to the text logfile.
    *   ``job_logger`` -- the logger used to write job summaries.
    *   ``input`` -- the absolute path to the copy of the input file in the working folder.
    *   ``settings`` -- a |Settings| instance for this job manager (see below).
    *   ``jobs`` -- a list of all jobs managed with this instance (in order of |run| calls).
    *   ``names`` -- a dictionary with names of jobs. For each name an integer value is stored indicating how many jobs with that basename have already been run.
    *   ``hashes`` -- a dictionary working as a hash-table for jobs.

    The *path* argument should be a path to a directory inside which the main working folder will be created. If ``None``, the directory from where the whole script was executed is used.

    The ``foldername`` attribute is initially set to the *folder* argument. If such a folder already exists (and ``use_existing_folder`` is False), the suffix ``.002`` is appended to *folder* and the number is increased (``.003``, ``.004``...) until a non-existsing name is found. If *folder* is ``None``, the name ``plams_workdir`` is used, followed by the same procedure to find a unique ``foldername``.

    The ``settings`` attribute is directly set to the value of *settings* argument (unlike in other classes where they are copied) and it should be a |Settings| instance with the following keys:

    *   ``hashing`` -- chosen hashing method (see |RPM|).
    *   ``counter_len`` -- length of number appended to the job name in case of a name conflict.
    *   ``remove_empty_directories`` -- if ``True``, all empty subdirectories of the working folder are removed on |finish|.

    """

    def __init__(
        self,
        settings: "Settings",
        path: Optional[str] = None,
        folder: Optional[str] = None,
        use_existing_folder: bool = False,
        job_logger: Optional[Logger] = None,
        load_all: bool = False,
        default_job_loader: Optional[Callable[[str], "Job"]] = None,
    ):

        self.settings = settings
        self.jobs: List[Job] = []
        self.names: Dict[str, int] = {}
        self.hashes: Dict[str, Job] = {}
        self.multijob_hashes: Dict[str, Job] = {}
        self._register_lock = threading.RLock()
        self._lazy_lock = threading.Lock()

        if path is None:
            ams_resultsdir = os.getenv("AMS_RESULTSDIR")
            if ams_resultsdir is not None and os.path.isdir(ams_resultsdir):
                self.path = ams_resultsdir
            else:
                self.path = os.getcwd()
        elif os.path.isdir(path):
            self.path = os.path.abspath(path)
        else:
            raise PlamsError("Invalid path: {}".format(path))

        basename = os.path.normpath(folder) if folder else "plams_workdir"
        self.foldername = basename

        if not use_existing_folder:
            n = 2
            while os.path.exists(opj(self.path, self.foldername)):
                self.foldername = basename + "." + str(n).zfill(3)
                n += 1

        self._workdir = opj(self.path, self.foldername)
        self.logfile = os.environ["SCM_LOGFILE"] if ("SCM_LOGFILE" in os.environ) else opj(self._workdir, "logfile")
        self.input = opj(self._workdir, "input")
        self._create_workdir = not (use_existing_folder and os.path.exists(self._workdir))
        self._job_logger = job_logger
        if use_existing_folder and load_all:
            self.load_all(self.workdir, register=True, default_job_loader=default_job_loader)

    @property
    def workdir(self) -> str:
        """
        Absolute path to the working directory
        """
        # Create the working directory only when first required
        # Avoids creating working directory only for e.g. load_job
        with self._lazy_lock:
            if self._create_workdir:
                os.mkdir(self._workdir)
                self._create_workdir = False
        return self._workdir

    @property
    def job_logger(self) -> Logger:
        """
        Logger used to write job summaries.
        If not specified on initialization, defaults to a csv logger with file ``job_logfile.csv``.
        """
        if self._job_logger is None:
            self._job_logger = get_logger(os.path.basename(self.workdir), fmt="csv")
            self._job_logger.configure(
                logfile_level=config.log.csv,
                logfile_path=opj(self.workdir, "job_logfile.csv"),
                csv_formatter=JobCSVFormatter,
                include_date=True,
                include_time=True,
            )
        return self._job_logger

    def load_all(self, path, register=False, default_job_loader: Optional[Callable[[str], "Job"]] = None):
        """Load all jobs from *path*.

        This function works as multiple executions of |load_job|. It searches for ``.dill`` files inside the directory given by *path*, yet not directly in it, but one level deeper. In other words, all files matching ``path/*/*.dill`` are used. That way a path to the main working folder of a previously run script can be used to import all the jobs run by that script.

        In case of partially failed |MultiJob| instances (some children jobs finished successfully, but not all) the function will search for ``.dill`` files in children folders. That means, if ``path/[multijobname]/`` contains some subfolders (for children jobs) but does not contail a ``.dill`` file (the |MultiJob| was not fully successful), it will look into these subfolders. This behavior is recursive up to any folder tree depth.

        The purpose of this function is to provide a quick way of restarting a script. Loading all successful jobs from the previous run prevents double work and allows the new execution of the script to proceed directly to the place where the previous execution failed.

        Jobs are loaded using default job manager stored in ``config.default_jobmanager``. If you wish to use a different one you can pass it as *jobmanager* argument of this function.

        Returned value is a dictionary containing all loaded jobs as values and absolute paths to ``.dill`` files as keys.
        """
        loaded_jobs = {}
        for foldername in filter(lambda x: isdir(opj(path, x)), os.listdir(path)):
            maybedill = opj(path, foldername, foldername + ".dill")
            maybedefault = opj(path, foldername, foldername + ".in")
            job = None
            if isfile(maybedill):
                job = self.load_job(maybedill, default_job_loader=default_job_loader, register=register)
                if job is not None:
                    loaded_jobs[os.path.abspath(maybedill)] = job
            elif isfile(maybedefault):
                job = self.load_job(opj(path, foldername), default_job_loader=default_job_loader, register=register)
                if job is not None:
                    loaded_jobs[os.path.abspath(opj(path, foldername))] = job
            if job is None:
                loaded_jobs.update(
                    self.load_all(path=opj(path, foldername), register=register, default_job_loader=default_job_loader)
                )
        if register:
            # self._register_jobs_in_hashes()
            self.jobs = list(self.jobs_hashed)
        return loaded_jobs

    def load_job(self, filename, register: bool = False, default_job_loader: Callable[[str], "Job"] = None):
        """Load previously saved job from *filename*.

        if *Filename* is a folder it will try to find in the folder the ``.dill`` file
        *Filename* should be a path to a ``.dill`` file in some job folder. A |Job| instance stored there is loaded and returned. All attributes of this instance removed before pickling are restored. That includes ``jobmanager``, ``path`` (the absolute path to the folder containing *filename* is used) and ``default_settings`` (a list containing only ``config.job``).

        If unpickling does not work and and default_job_loader is not None filename.parent is taken as possible path , it will try to load the job by calling default_job_loader(filename.parent)
        See |pickling| for details.
        """
        try:
            import dill as pickle
        except ImportError:
            import pickle

        filename = os.path.abspath(filename)
        path = None
        if os.path.isfile(filename):  # if it is a dill file
            path = os.path.dirname(filename)
        elif os.path.isdir(filename):  # if is a job folder
            path = filename
            filename = opj(filename, os.path.basename(filename) + ".dill")

        path_exists = path is not None and os.path.exists(path)
        dill_exists_or_try_default_loader = os.path.isfile(filename) or default_job_loader is not None
        if not path_exists or not dill_exists_or_try_default_loader:
            if not path_exists:
                raise FileError("File {} not present".format(filename))
            else:
                raise FileError("Default job loader {} is None".format(default_job_loader))

        path = os.path.dirname(filename)
        # if os.path.isfile(filename) or
        if path in set([j.path for j in self.jobs]):
            return

        job = None
        if os.path.exists(filename):
            with open(filename, "rb") as f:
                try:
                    job = pickle.load(f)
                except Exception as e:
                    log("Unpickling of {} failed. Caught the following Exception:\n{}".format(filename, e), 1)

        if job is None and default_job_loader is not None:
            try:
                job = default_job_loader(path)
            except Exception as e:
                log(f"Default job loader of {filename} failed. Caught the following Exception:\n{e}", 1)
                job = None
        if job is None:
            return None
        self._setstate(job, path, register=register)
        return job

    def _setstate(self, job, path, parent=None, register: bool = False):
        job.parent = parent
        job.jobmanager = self
        job.default_settings = [config.job]
        job.path = path
        if isinstance(job, MultiJob):
            job._lock = threading.Lock()
            for child in job:
                self._setstate(child, opj(path, child.name), job, register=register)
            for otherjob in job.other_jobs():
                self._setstate(otherjob, opj(path, otherjob.name), job, register=register)

        job.results.refresh()
        h = job.hash()
        multi_h = None
        if isinstance(job, MultiJob):
            multi_h = job.children_hash()

        if h is not None:
            self.hashes[h] = job
        if multi_h is not None:
            self.multijob_hashes[multi_h] = job

        for key in job._dont_pickle:
            job.__dict__[key] = None

        if register:
            fname = re.sub(r"(\.\d{%i})+$" % (self.settings.counter_len), "", job._full_name())
            if fname in self.names:
                self.names[fname] += 1
            else:
                self.names[fname] = 1

    @property
    def jobs_hashed(self):
        return itertools.chain(self.hashes.values(), self.multijob_hashes.values())

    @cached_property
    def names_count(self):
        names = defaultdict(lambda: 0)
        for j in self.jobs_hashed:
            fname = re.sub(r"(\.\d{%i})+$" % (self.settings.counter_len), "", j._full_name())
            names[fname] += 1
        return dict(names)

    def _register_jobs_in_hashes(self):
        """It would have worked if the jobs_hashed would contains also the copied jobs!"""
        # register the job_name taking into account it might be in a MultiJob
        for job in self.jobs_hashed:
            orgfname, fname = self._register_job_name(job)
            job.jobmanager = self
        # update the names and jobs accordingly with the jobs saved in self.jobs_hashed
        self.names = self.names_count.copy()
        self.jobs = list(self.jobs_hashed)

    def remove_job(self, job):
        """Remove *job* from the job manager. Forget its hash."""
        with self._register_lock:
            if job in self.jobs:
                self.jobs.remove(job)
                job.jobmanager = None
            h = job.hash()
            if h in self.hashes and self.hashes[h] == job:
                del self.hashes[h]
            if isinstance(job, MultiJob):
                for child in job:
                    self.remove_job(child)
                for otherjob in job.other_jobs():
                    self.remove_job(otherjob)
            shutil.rmtree(job.path)

    def _register(self, job: "Job"):
        """Register the *job*. Register job's name (rename if needed) and create the job folder.

        If a job with the same name was already registered, *job* is renamed by appending consecutive integers. The number of digits in the appended number is defined by the ``counter_len`` value in ``settings``.
        Note that jobs whose name already contains a counting suffix, e.g. ``myjob.002`` will have the suffix stripped as the very first step.
        """
        with self._register_lock:

            log("Registering job {}".format(job.name), 7)
            job.jobmanager = self

            # If the name ends with the counting suffix, e.g. ".002", remove it.
            # The suffix is just not part of a legitimate job name and users will have to live with it potentially changing.
            orgfname, fname = self._register_job_name(job)
            if fname != orgfname:
                log("Renaming job {} to {}".format(orgfname, fname), 3)

            if job.path is None:
                if job.parent:
                    job.path = opj(job.parent.path, job.name)
                else:
                    job.path = opj(self.workdir, job.name)
            os.mkdir(job.path)

            self.jobs.append(job)
            job.status = JobStatus.REGISTERED
            log("Job {} registered".format(job.name), 7)

    def _register_job_name(self, job: "Job"):
        orgfname = job._full_name()
        job.name = re.sub(r"(\.\d{%i})+$" % (self.settings.counter_len), "", job.name)
        fname = job._full_name()
        if fname in self.names:
            self.names[fname] += 1
            job.name += "." + str(self.names[fname]).zfill(self.settings.counter_len)
            fname = job._full_name()
        else:
            self.names[fname] = 1
        return orgfname, fname

    def _check_hash(self, job):
        """Calculate the hash of *job* and, if it is not ``None``, search previously run jobs for the same hash. If such a job is found, return it. Otherwise, return ``None``"""
        h = job.hash()
        if h is not None:
            with self._register_lock:
                if h in self.hashes:
                    prev = self.hashes[h]
                    log("Job {} previously run as {}, using old results".format(job.name, prev.name), 1)
                    return prev
                else:
                    self.hashes[h] = job
        return None

    def _clean(self):
        """Clean all registered jobs according to the ``save`` parameter in their ``settings``. If ``remove_empty_directories`` is ``True``,  traverse the working directory and delete all empty subdirectories."""
        log("Cleaning job manager", 7)

        for job in self.jobs:
            job.results._clean(job.settings.save)

        if self.settings.remove_empty_directories:
            for root, dirs, files in os.walk(self._workdir, topdown=False):
                for dirname in dirs:
                    fullname = opj(root, dirname)
                    if not os.listdir(fullname):
                        os.rmdir(fullname)

        if self._job_logger is not None:
            self._job_logger.close()

        log("Job manager cleaned", 7)

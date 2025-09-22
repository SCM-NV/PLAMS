import os
import re
import threading
from os.path import join as opj
from pathlib import Path
from typing import TYPE_CHECKING, Optional, List, Dict, Tuple

from scm.plams.core.basejob import MultiJob
from scm.plams.core.errors import FileError, PlamsError
from scm.plams.core.functions import get_logger, log, config, _get_dir_for_jobs
from scm.plams.core.logging import Logger
from scm.plams.core.formatters import JobCSVFormatter

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
    ):

        self.settings = settings
        self.jobs: List[Job] = []
        self.names: Dict[str, int] = {}
        self._job_full_name_map: Dict[Job, Tuple[str, int]] = {}
        self.hashes: Dict[str, Job] = {}

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
            raise PlamsError(f"Invalid path: {path}")

        basename = os.path.normpath(folder) if folder else "plams_workdir"
        self.foldername = basename

        if not use_existing_folder:
            n = 2
            while os.path.exists(opj(self.path, self.foldername)):
                self.foldername = basename + "." + str(n).zfill(3)
                n += 1

        self._workdir = Path(self.path, self.foldername)
        self.logfile = os.environ["SCM_LOGFILE"] if ("SCM_LOGFILE" in os.environ) else opj(self._workdir, "logfile")
        self.input = opj(self._workdir, "input")
        self._create_workdir = not (use_existing_folder and self._workdir.exists())
        self._job_logger = job_logger

    @property
    def workdir(self) -> str:
        """
        Absolute path to the |JobManager| working directory.

        This is the top-level directory which contains subdirectories and job directories.
        """
        # Create the working directory only when first required
        # Avoids creating working directory only for e.g. load_job
        with self._lazy_lock:
            if self._create_workdir:
                os.mkdir(self._workdir)
                self._create_workdir = False
        return str(self._workdir.resolve())

    @property
    def current_dir_for_jobs(self) -> Path:
        """
        Absolute path to the current directory where new jobs will be run.

        This is the directory which will directly contain the job directories for any newly run jobs.
        It is located within the ``workdir``.
        """
        rel_dir = _get_dir_for_jobs()
        path = self._workdir / rel_dir if rel_dir else self._workdir
        return path.resolve()

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

    def load_job(self, filename: str) -> Optional["Job"]:
        """Load previously saved job from *filename*.

        *Filename* should be a path to a ``.dill`` file in some job folder. A |Job| instance stored there is loaded and returned. All attributes of this instance removed before pickling are restored. That includes ``jobmanager``, ``path`` (the absolute path to the folder containing *filename* is used) and ``default_settings`` (a list containing only ``config.job``).

        See |pickling| for details.
        """
        try:
            import dill as pickle
        except ImportError:
            import pickle

        def setstate(job: "Job", path: str, parent: Optional["MultiJob"] = None) -> None:
            job.parent = parent
            job.jobmanager = self
            job.default_settings = []
            job.path = path
            if isinstance(job, MultiJob):
                job._lock = threading.Lock()
                for child in job:
                    setstate(child, opj(path, child.name), job)
                for otherjob in job.other_jobs():
                    setstate(otherjob, opj(path, otherjob.name), job)

            job.results.refresh()
            h = job.hash()
            if h is not None:
                self.hashes[h] = job
            for key in job._dont_pickle:
                job.__dict__[key] = None

        if os.path.isfile(filename):
            filename = os.path.abspath(filename)
        else:
            raise FileError(f"File {filename} not present")
        path = os.path.dirname(filename)
        with open(filename, "rb") as f:

            def resolve_missing_attributes(j: "Job") -> None:
                # For backwards compatibility (before attributes added/converted to properties)
                if not hasattr(j, "_status"):
                    j._status = j.__dict__["status"]
                    j._status_log = []
                if not hasattr(j, "_error_msg"):
                    j._error_msg = None
                if isinstance(j, MultiJob):
                    for child in j:
                        resolve_missing_attributes(child)
                    for otherjob in j.other_jobs():
                        resolve_missing_attributes(otherjob)

            try:
                job = pickle.load(f)
                resolve_missing_attributes(job)
            except Exception as e:
                log(f"Unpickling of {filename} failed. Caught the following Exception:\n{e}", 1)
                return None

        setstate(job, path)
        return job

    def rename_job(self, job: "Job", name: str) -> None:
        """
        Rename the job in the job manager, de-registering the old name and re-registering it under the new name.

        .. note::
            The job manager is not responsible for moving the job files and directories

        :param job: job to rename
        :param name: new name for the job
        """
        with self._register_lock:
            if name == job.name:
                return

            log(f"Renaming job {job.name} to {name}", 7)

            if job.jobmanager is not self:
                raise PlamsError(f"Cannot rename job '{job.name}' as it does not belong to this job manager")

            orig_name = job.name

            # renaming a job always preserves the parent directory it was run in, so determine this
            try:
                if job.parent is not None or job.path is None:
                    rel_dir = Path(".")
                else:
                    rel_dir = Path(job.path).parent.resolve().relative_to(self._workdir)
            except ValueError:
                raise PlamsError(
                    f"Cannot rename job '{job.name}' as it does not reside in the working directory of this job manager: '{self._workdir}'."
                )

            # deregister original job from the job manager then change the name and register it again
            self.remove_job(job)
            job.name = name
            job.path = None
            self._register(job, rel_dir)

            # child jobs have been removed, so re-register these to update the paths
            def reregister_child_job(child_job: "Job") -> None:
                child_job.path = None
                self._register(child_job, rel_dir_for_jobs=Path("."), auto_rename=False)

            MultiJob.apply_to_children(job, reregister_child_job, recursive=True)

            log(f"Job {orig_name} renamed to {job.name}", 7)

    def remove_job(self, job: "Job") -> None:
        """
        Remove *job* from the job manager.
        This removes its hash and resets the name count to the last remaining job with the same name.

        .. note::
            The job manager is not responsible for removing the job files and directories

        :param job: job to remove
        """
        with self._register_lock:
            log(f"Removing job {job.name}", 7)

            if job in self.jobs:
                self.jobs.remove(job)
                job.jobmanager = None
            if job in self._job_full_name_map:
                name, cnt = self._job_full_name_map.pop(job)
                if name in self.names:
                    if cnt == self.names[name]:
                        if cnt == 1:
                            self.names.pop(name)
                        else:
                            remaining = [c for n, c in self._job_full_name_map.values() if n == name]
                            if remaining:
                                self.names[name] = max(remaining)
                            else:
                                self.names.pop(name)
            h = job.hash()
            if h in self.hashes and self.hashes[h] == job:
                del self.hashes[h]
            MultiJob.apply_to_children(job, self.remove_job)

            log(f"Job {job.name} removed", 7)

    def _register(self, job: "Job", rel_dir_for_jobs: Optional[Path] = None, auto_rename: bool = True) -> None:
        """Register the *job*. Register job's name. Rename if needed and ``auto_rename=True``.

        If a job with the same name was already registered, *job* is renamed by appending consecutive integers.
        The number of digits in the appended number is defined by the ``counter_len`` value in ``settings``.
        Note that jobs whose name already contains a counting suffix, e.g. ``myjob.002`` will have the suffix stripped as the very first step.

        If a relative directory for the jobs is provided, jobs will be registered in this directory, relative to the job manager working directory.
        Otherwise, the value from :func:`~scm.plams.core.functions._get_dir_for_jobs` will be used.
        An exception is for registering a child of a |MultiJob|, when the parent directory is always used.

        If a job with the same name in the same directory is already registered and ``auto_rename=False``, raises a ``PlamsError``.

        .. note::
            The job manager is not responsible for creating the job files and directories

        :param job: job to register
        :param rel_dir_for_jobs: relative directory to the job manager working directory, in which to register jobs
        :param auto_rename: whether to automatically rename job if a job with the same name already exists
        """
        with self._register_lock:

            log(f"Registering job {job.name}", 7)
            job.jobmanager = self

            # get current directory for jobs
            # this directory should be used unless job is within a multi-job (as then the parent job directory should be used)
            if rel_dir_for_jobs:
                dir_for_jobs = (self._workdir / rel_dir_for_jobs).resolve()
            else:
                dir_for_jobs = self.current_dir_for_jobs
                # N.B. resolve does not remove the Windows extended prefix if the directory does not yet exist
                # so sanitise the directory for jobs here when finding the relative path
                dir_for_jobs_str = str(dir_for_jobs)
                clean_dir_for_jobs = Path(
                    dir_for_jobs_str[4:] if dir_for_jobs_str.startswith(("\\\\?\\", "//?/")) else dir_for_jobs_str
                ).resolve()
                rel_dir_for_jobs = clean_dir_for_jobs.relative_to(self.workdir)
            rel_dir_for_jobs = rel_dir_for_jobs if rel_dir_for_jobs != Path(".") else None

            # check the name of the current job for collisions with existing jobs, and rename if required
            orig_job_full_name = job._full_name(rel_dir_for_jobs)
            pattern = rf"(?:\.(\d{{{self.settings.counter_len}}}))+?$"
            match = re.search(pattern, job.name)
            job_full_name_counter = int(match.group(1)) if match else 1
            job_full_name_no_counter = re.sub(pattern, "", job._full_name(rel_dir_for_jobs))

            # if the name is not yet registered there are no collisions
            if job_full_name_no_counter not in self.names:
                self.names[job_full_name_no_counter] = job_full_name_counter
            else:
                # otherwise rename the job with an incrementing counter suffix
                if auto_rename:
                    self.names[job_full_name_no_counter] += 1
                    job_full_name_counter = self.names[job_full_name_no_counter]
                    job.name = (
                        f"{re.sub(pattern, '', job.name)}.{str(job_full_name_counter).zfill(self.settings.counter_len)}"
                    )
                    # self.names[job_full_name_no_counter] = job_full_name_counter
                    if orig_job_full_name != job._full_name(rel_dir_for_jobs):
                        log(f"Renaming job {orig_job_full_name} to {job._full_name(rel_dir_for_jobs)}", 3)
                # alternatively do a strict check that the job with this suffix is not already registered
                else:
                    counts_in_use = [c for n, c in self._job_full_name_map.values() if n == job_full_name_no_counter]
                    if job_full_name_counter not in counts_in_use:
                        max_count_in_use = max(counts_in_use) if counts_in_use else 1
                        if job_full_name_counter > max_count_in_use:
                            self.names[job_full_name_no_counter] = job_full_name_counter
                    else:
                        raise PlamsError(f"Job {job.name} already registered and cannot automatically be renamed", 1)

            if job.path is None:
                if job.parent and job.parent.path:
                    job.path = opj(job.parent.path, job.name)
                else:
                    job.path = opj(dir_for_jobs, job.name)

            self.jobs.append(job)
            self._job_full_name_map[job] = (job_full_name_no_counter, job_full_name_counter)

            log(f"Job {job.name} registered", 7)

    def _check_hash(self, job: "Job") -> Optional["Job"]:
        """Calculate the hash of *job* and, if it is not ``None``, search previously run jobs for the same hash. If such a job is found, return it. Otherwise, return ``None``"""
        h = job.hash()
        if h is not None:
            with self._register_lock:
                if h in self.hashes:
                    prev = self.hashes[h]
                    log(f"Job {job.name} previously run as {prev.name}, using old results", 1)
                    return prev
                else:
                    self.hashes[h] = job
        return None

    def _clean(self) -> None:
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

import os
import re
import shutil
import subprocess
import sys
import threading
import types
from os.path import dirname, expandvars, isdir, isfile
from os.path import join as opj
from typing import Dict, Iterable, Optional, TYPE_CHECKING
import atexit
from importlib.util import find_spec
import functools

from scm.plams.core.logging import get_logger
from scm.plams.core.errors import FileError, MissingOptionalPackageError
from scm.plams.core.private import retry
from scm.plams.core.settings import Settings, ConfigSettings
from scm.plams.core.enums import JobStatus

if TYPE_CHECKING:
    from scm.plams.core.jobmanager import JobManager
    from scm.plams.core.basejob import Job

__all__ = [
    "init",
    "finish",
    "log",
    "load",
    "load_all",
    "delete_job",
    "add_to_class",
    "add_to_instance",
    "requires_optional_package",
    "config",
    "read_molecules",
    "read_all_molecules_in_xyz_file",
]

# Initialise config with default values for standard settings
# An initial call to '_init' and a final call to '_finish' are made later in this module
config = ConfigSettings()
# ===========================================================================

_logger = get_logger("plams")


def log(message: str, level: int = 0) -> None:
    """
    Log a *message* with verbosity *level*.

    Logs are printed independently to the text logfile (a file called ``logfile`` in the main working folder) and to the standard output.
    If *level* is equal or lower than verbosity (defined by ``config.log.file`` or ``config.log.stdout``) the message is printed.
    By convention in PLAMS, level should be between 0-7, with 0 indicating no loggin and 7 indicating the most verbose logging.
    Date and/or time can be added based on ``config.log.date`` and ``config.log.time``.
    All logging activity is thread safe.
    """
    if config.init and "log" in config:
        logfile = config.default_jobmanager.logfile if config["default_jobmanager"] is not None else None
        _logger.configure(config.log.stdout, config.log.file, logfile, config.log.date, config.log.time)
    else:
        # By default write to stdout with level 3
        _logger.configure(3)
    _logger.log(message, level)


# ===========================================================================


def _init() -> None:
    """
    Internal initialisation of the config settings, which should only be called once.
    This finishes setup of the config settings, adding any settings from a plams_defaults file, any slurm settings, and marking the settings as initialised.
    It also attempts to import the dill package and logs a warning if this was unsuccessful.
    """

    try:
        import dill  # noqa F401
    except ImportError:
        log(
            "WARNING: importing dill package failed. Falling back to the default pickle module. Expect problems with pickling",
            1,
        )

    _load_defaults_file()

    config.slurm = _init_slurm() if "SLURM_JOB_ID" in os.environ else None

    config.init = True


def _load_defaults_file() -> Optional[str]:
    """
    Searches for a raw python defaults file and executes it to apply settings to the global config.
    The defaults file is specified either by the ``$PLAMSDEFAULTS`` environment variable,
    the path ``$AMSHOME/scripting/scm/plams/plams_defaults``, or the path ``../plams_defaults`` relative to the current file.

    Note that this is a legacy functionality which is no longer encouraged but still used.

    :return: path of the defaults file used to populate the global config
    """
    if "PLAMSDEFAULTS" in os.environ and isfile(expandvars("$PLAMSDEFAULTS")):
        defaults = expandvars("$PLAMSDEFAULTS")
    elif "AMSHOME" in os.environ and isfile(opj(expandvars("$AMSHOME"), "scripting", "scm", "plams", "plams_defaults")):
        defaults = opj(expandvars("$AMSHOME"), "scripting", "scm", "plams", "plams_defaults")
    else:
        defaults = opj(dirname(dirname(__file__)), "plams_defaults")
    if isfile(defaults):
        with open(defaults, "r") as f:
            exec(compile(f.read(), defaults, "exec"))
        return defaults
    else:
        return None


def _init_slurm() -> Optional[Settings]:
    """
    If PLAMS is running under Slurm, query some information about the batch system and set up the environment for the PLAMS/Slurm integration
    """
    ret = Settings()
    if "SLURM_JOB_ID" not in os.environ:
        log("Slurm setup aborted: SLURM_JOB_ID is not set")
        return None
    try:
        srun = subprocess.run(["srun", "--version"], stdout=subprocess.PIPE, timeout=60)
    except subprocess.TimeoutExpired:
        log("Slurm setup failed: timeout for srun --version")
        return None
    if srun.returncode != 0:
        log("Slurm setup failed: srun --version exited with non-zero return code")
        return None
    try:
        ret.slurm_version = srun.stdout.decode().split()[1].split(".")
        int(ret.slurm_version[0])
        int(ret.slurm_version[1])
    except Exception:
        log("Slurm setup failed: could not determine Slurm version")
        return None
    if int(ret.slurm_version[0]) < 15:
        log("Slurm setup failed: Slurm version >=15 is required for Slurm/PLAMS integration")
        return None
    if "SLURM_TASKS_PER_NODE" not in os.environ:
        log("Slurm setup failed: SLURM_TASKS_PER_NODE is not set")
        return None
    ret.tasks_per_node = []
    for tasks in os.environ["SLURM_TASKS_PER_NODE"].split(","):
        try:
            if "(" in tasks:
                tasks, _, num_nodes = tasks.partition("(")
                num_nodes = num_nodes[1:-1]
                ret.tasks_per_node += [int(tasks)] * int(num_nodes)
            else:
                ret.tasks_per_node.append(int(tasks))
        except Exception:
            log("Slurm setup failed: can not determine the number of tasks per node")
            return None
    ret.tasks_per_node.sort(reverse=True)

    # General setup of the environment when running under SLURM.
    os.environ["SCM_SRUN_OPTIONS"] = "-m block:block:block,NoPack --use-min-nodes"
    # There was a change in the behaviour of the srun --exclusive flag in Slurm 20.11.
    # Newer versions should use --exact instead, see: https://www.nsc.liu.se/support/batch-jobs/slurm/20.11/
    if int(ret.slurm_version[0]) >= 21 or (int(ret.slurm_version[0]) == 20 and int(ret.slurm_version[1]) >= 11):
        os.environ["SCM_SRUN_OPTIONS"] += " --exact"
    else:
        os.environ["SCM_SRUN_OPTIONS"] += " --exclusive"

    return ret


# Initial call to finish set up of config
_init()


def init(
    path: Optional[str] = None,
    folder: Optional[str] = None,
    config_settings: Optional[Dict] = None,
    quiet: bool = False,
    use_existing_folder: bool = False,
) -> None:
    """
    Initialize PLAMS environment. Set up the global ``config`` and its default |JobManager|.

    This function only needs to be called to override default values on the global ``config``, without interacting with it directly.
    For example, both of the following approaches are equivalent:

    .. code-block:: python

        # Call to init
        init(path="my/path", folder="my_folder", config_settings=my_settings)

        # Direct modification of config
        config.update(my_settings) # direct modification of config
        config.default_jobmanager = JobManager(config.jobmanager, "my/path", "my_folder")

    .. note::
        A call to ``init()`` with no arguments is therefore a no-op, and is not required. It will however force initialisation of the lazy default |JobManager|.
        This is in contrast to AMS2024 and earlier versions, where this function had to be called before any other PLAMS command could be executed.

    :param path: for the default |JobManager|
    :param folder: for the default |JobManager|
    :param use_existing_folder: for the default |JobManager|
    :param config_settings: |Settings| to update config with - these will overwrite any existing items
    :param quiet: do not log header with information about the PLAMS environment
    """

    if config.init and config._explicit_init:
        return

    # Build the config settings via:
    # * (Re)Initialise config with default values
    # * Update with any values from a PLAMS defaults file
    # * Update with any settings passed to the method
    config.update(ConfigSettings())
    defaults_file = _load_defaults_file()
    config.update(config_settings or {})

    from scm.plams.core.jobmanager import JobManager

    config.default_jobmanager = JobManager(config.jobmanager, path, folder, use_existing_folder)

    config.slurm = _init_slurm() if "SLURM_JOB_ID" in os.environ else None

    if not quiet:
        log("Running PLAMS located in {}".format(dirname(dirname(__file__))), 5)
        log("Using Python {}.{}.{} located in {}".format(*sys.version_info[:3], sys.executable), 5)
        if defaults_file is not None:
            log("PLAMS defaults were loaded from {}".format(defaults_file), 5)
        log("PLAMS environment initialized", 5)
        log("PLAMS working folder: {}".format(config.default_jobmanager.workdir), 1)

    config.init = True
    config._explicit_init = True


# ===========================================================================


def _finish():
    """
    Internal clean up of the PLAMS environment, which will be called at the end of the script.
    """

    if not config.init:
        return

    for thread in threading.enumerate():
        if thread.name == "plamsthread":
            thread.join()

    # Only clean the default lazy job manager if it has been initialised
    # as otherwise accessing it will create an empty workdir
    if config["default_jobmanager"] is not None:
        config.default_jobmanager._clean()

        if config.erase_workdir is True:
            shutil.rmtree(config.default_jobmanager.workdir)

    config.init = False


def finish(otherJM: Optional[Iterable["JobManager"]] = None):
    """
    Clean up the PLAMS environment. This can be explicitly called for |cleaning| to take place.
    If you used some other job managers than just the default one, they need to be passed as *otherJM*.

    A call is registered to be performed at exit, which will perform cleaning for the default job manager.

    .. note::
        A call to ``finish()`` with no arguments once at the end of a script is therefore a no-op, and is not required.
        This is in contrast to AMS2024 and earlier versions, where this function had to be called at the end of a script for |cleaning| to take place.

    :param otherJM: additional job managers used in a workflow to be cleaned
    """
    if not config.init:
        return

    _finish()

    if otherJM:
        for jm in otherJM:
            jm._clean()

    log("PLAMS environment cleaned up successfully", 5)
    log("PLAMS run finished. Goodbye", 3)


# Register call to _finish on workflow end
atexit.register(_finish)


# ===========================================================================


def load(filename):
    """Load previously saved job from ``.dill`` file. This is just a shortcut for |load_job| method of the default |JobManager| ``config.default_jobmanager``."""
    return config.default_jobmanager.load_job(filename)


# ===========================================================================


def load_all(path, jobmanager=None):
    """Load all jobs from *path*.

    This function works as multiple executions of |load_job|. It searches for ``.dill`` files inside the directory given by *path*, yet not directly in it, but one level deeper. In other words, all files matching ``path/*/*.dill`` are used. That way a path to the main working folder of a previously run script can be used to import all the jobs run by that script.

    In case of partially failed |MultiJob| instances (some children jobs finished successfully, but not all) the function will search for ``.dill`` files in children folders. That means, if ``path/[multijobname]/`` contains some subfolders (for children jobs) but does not contail a ``.dill`` file (the |MultiJob| was not fully successful), it will look into these subfolders. This behavior is recursive up to any folder tree depth.

    The purpose of this function is to provide a quick way of restarting a script. Loading all successful jobs from the previous run prevents double work and allows the new execution of the script to proceed directly to the place where the previous execution failed.

    Jobs are loaded using default job manager stored in ``config.default_jobmanager``. If you wish to use a different one you can pass it as *jobmanager* argument of this function.

    Returned value is a dictionary containing all loaded jobs as values and absolute paths to ``.dill`` files as keys.
    """
    jm = jobmanager or config.default_jobmanager
    loaded_jobs = {}
    for foldername in filter(lambda x: isdir(opj(path, x)), os.listdir(path)):
        maybedill = opj(path, foldername, foldername + ".dill")
        if isfile(maybedill):
            job = jm.load_job(maybedill)
            if job:
                loaded_jobs[os.path.abspath(maybedill)] = job
        else:
            loaded_jobs.update(load_all(path=opj(path, foldername), jobmanager=jm))
    return loaded_jobs


# ===========================================================================


@retry()
def delete_job(job: "Job"):
    """Remove *job* from its corresponding |JobManager| and delete the job folder from the disk. Mark *job* as 'deleted'."""

    if job.status != JobStatus.CREATED:
        job.results.wait()

    # In case job.jobmanager is None, run() method was not called yet, so no JobManager knows about this job and no folder exists.
    if job.jobmanager is not None:
        job.jobmanager.remove_job(job)

    if job.parent is not None:
        job.parent.remove_child(job)

    job.status = JobStatus.DELETED
    job._log_status(5)


# ===========================================================================


def read_molecules(folder, formats=None):
    """Read all molecules from *folder*.

    Read all the files present in *folder* with extensions compatible with :meth:`Molecule.read<scm.plams.mol.molecule.Molecule.read>`. Returned value is a dictionary with keys being molecule names (filename without extension) and values being |Molecule| instances.

    The optional argument *formats* can be used to narrow down the search to files with specified extensions::

        molecules = read_molecules('mymols', formats=['xyz', 'pdb'])
    """
    from scm.plams.mol.molecule import Molecule

    extensions = formats or list(Molecule._readformat.keys())
    is_valid = lambda x: isfile(opj(folder, x)) and any([x.endswith("." + ext) for ext in extensions])
    filenames = filter(is_valid, os.listdir(folder))
    ret = {}
    for f in filenames:
        m = Molecule(opj(folder, f))
        ret[m.properties.name] = m
    return ret


# ===========================================================================


def read_all_molecules_in_xyz_file(filename):
    """The .xyz format allows to store multiple geometries on a single file (such file is essentially a concatenated series of xyz files)

    This function returns a *list* of all molecules found in the file *filename*

    *filename*: path (absolute or relative to the current working directory) to the xyz file.
    """
    from scm.plams.mol.molecule import Molecule

    mols = []
    with open(filename, "r") as f:
        while True:
            try:
                mol = Molecule()
                mol.readxyz(f)
                mols.append(mol)
            except FileError:
                break
    return mols


# ===========================================================================


def add_to_class(classname):
    """Add decorated function as a method to the whole class *classname*.

    The decorated function should follow a method-like syntax, with the first argument ``self`` that references the class instance.
    Example usage::

        @add_to_class(ADFResults)
        def get_energy(self):
            return self.readkf('Energy', 'Bond Energy')

    After executing the above code all instances of ``ADFResults`` in the current script (even the ones created beforehand) are enriched with ``get_energy`` method that can be invoked by::

        someadfresults.get_energy()

    The added method is accessible also from subclasses of *classname* so ``@add_to_class(Results)`` in the above example will work too.

    If *classname* is |Results| or any of its subclasses, the added method will be wrapped with the thread safety guard (see |parallel|).
    """
    from scm.plams.core.results import ApplyRestrict, _restrict

    def decorator(func):
        if isinstance(classname, ApplyRestrict):
            func = _restrict(func)
        setattr(classname, func.__name__, func)

    return decorator


# ===========================================================================


def add_to_instance(instance):
    """Add decorated function as a method to one particular *instance*.

    The decorated function should follow a method-like syntax, with the first argument ``self`` that references the class instance.
    Example usage::

        results = myjob.run()

        @add_to_instance(results)
        def get_energy(self):
            return self.readkf('Energy', 'Bond Energy')

        results.get_energy()

    The added method is accessible only for that one particular instance and it overrides any methods with the same name defined on a class level (in original class' source) or added with :func:`add_to_class` decorator.

    If *instance* is an instance of |Results| or any of its subclasses, the added method will be wrapped with the thread safety guard (see |parallel|).
    """
    from scm.plams.core.results import Results, _restrict

    def decorator(func):
        if isinstance(instance, Results):
            func = _restrict(func)
        func = types.MethodType(func, instance)
        setattr(instance, func.__func__.__name__, func)

    return decorator


# ===========================================================================


def requires_optional_package(package_name: str):
    """
    Ensures a given package is available before running a function, otherwise raises an ImportError.
    This can be used to check for optional dependencies which are required for specific functionality.
    :param package_name: name of the required package
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if find_spec(package_name) is None:
                raise MissingOptionalPackageError(package_name)
            return func(*args, **kwargs)

        return wrapper

    return decorator


# ===========================================================================


def parse_heredoc(bash_input: str, heredoc_delimit: str = "eor") -> str:
    """Take a string and isolate the content of a bash-style `Here Document`_.

    The input string, *bash_input*, is returned unaltered if no heredoc block is found.
    If multiple heredoc blocks are present only the first one is returned.

    An example bash input file for ADF:

    .. code:: bash

        #!/bin/bash

        $AMSBIN/adf << eor
        ATOMS
            1.H  0.0  0.0  0.0
            2.H  1.0  0.0  0.0
        END

        BASIS
            type TZ2P
        END

        XC
            GGA BP86
        END
        eor

        echo "Job finished"

    The matching :func:`parse_heredoc` output:

    .. code:: python

        >>> filename: str = ...  # The bash input file
        >>> with open(filename, 'r') as f:
        ...     output = parse_heredoc(f.read())

        >>> print(output)
        ATOMS
            1.H  0.0  0.0  0.0
            2.H  1.0  0.0  0.0
        END

        BASIS
            type TZ2P
        END

        XC
            GGA BP86
        END

    .. _`Here Document`: https://en.wikipedia.org/wiki/Here_document

    """
    # Find the start of the heredoc block
    start_pattern = r"<<(-)?(\s+)?{}".format(heredoc_delimit)
    start_heredoc = re.search(start_pattern, bash_input)
    if not start_heredoc:
        return bash_input

    # Find the end of the heredoc block
    end_pattern = r"\n(\s+)?{}(\s+)?\n".format(heredoc_delimit)
    end_heredoc = re.search(end_pattern, bash_input)

    # Prepare the slices
    if end_heredoc:
        i, j = start_heredoc.end(), end_heredoc.start()
    else:
        err = f"parse_heredoc: failed to find the final '{heredoc_delimit}' delimiter"
        raise ValueError(err)

    # Grab heredoced block and parse it
    _, ret = bash_input[i:j].split("\n", maxsplit=1)
    return ret

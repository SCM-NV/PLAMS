import os
import re
import shutil
import sys
import threading
import time
import types
import warnings
from typing import Callable, Dict, NoReturn

from os.path import join as opj
from os.path import isfile, isdir, expandvars, dirname

from .errors import PlamsError, FileError
from .settings import Settings
from .private import retry

__all__ = ['init', 'finish', 'log', 'load', 'load_all', 'delete_job', 'add_to_class', 'add_to_instance', 'config', 'read_molecules', 'read_all_molecules_in_xyz_file']

config = Settings()
config.init = False

#===========================================================================


def init(path=None, folder=None, config_settings:Dict=None):
    """Initialize PLAMS environment. Create global ``config`` and the default |JobManager|.

    An empty |Settings| instance is created and populated with default settings by executing ``plams_defaults``. The following locations are used to search for the defaults file, in order of precedence:

    *   If ``$PLAMSDEFAULTS`` variable is in your environment and it points to a file, this file is used (executed as a Python script).
    *   If ``$AMSHOME`` variable is in your environment and ``$AMSHOME/scripting/scm/plams/plams_defaults`` exists, it is used.
    *   Otherwise, the path ``../plams_defaults`` relative to the current file (``functions.py``) is checked. If defaults file is not found there, an exception is raised.

    Then a |JobManager| instance is created as ``config.default_jobmanager`` using *path* and *folder* to determine the main working folder. Settings for this instance are taken from ``config.jobmanager``. If *path* is not supplied, the current directory is used. If *folder* is not supplied, ``plams_workdir`` is used.

    Optionally, an additional `dict` (or |Settings| instance) can be provided to the `config_settings` argument which will be used to update the values from the ``plams_defaults``.

    .. warning::
      This function **must** be called before any other PLAMS command can be executed. Trying to do anything without it results in a crash. See also |master-script|.
    """

    if config.init:
        return

    if 'PLAMSDEFAULTS' in os.environ and isfile(expandvars('$PLAMSDEFAULTS')):
        defaults = expandvars('$PLAMSDEFAULTS')
    elif 'AMSHOME' in os.environ and isfile(opj(expandvars('$AMSHOME'), 'scripting', 'scm', 'plams', 'plams_defaults')):
        defaults = opj(expandvars('$AMSHOME'), 'scripting', 'scm', 'plams', 'plams_defaults')
    else:
        defaults = opj(dirname(dirname(__file__)), 'plams_defaults')
        if not isfile(defaults):
            raise PlamsError('plams_defaults not found, please set PLAMSDEFAULTS or AMSHOME in your environment')
    with open(defaults, 'r') as f:
        exec(compile(f.read(), defaults, 'exec'))

    config.update(config_settings or {})

    from .jobmanager import JobManager
    config.default_jobmanager = JobManager(config.jobmanager, path, folder)

    log('Running PLAMS located in {}'.format(dirname(dirname(__file__))), 5)
    log('Using Python {}.{}.{} located in {}'.format(*sys.version_info[:3], sys.executable), 5)
    log('PLAMS defaults were loaded from {}'.format(defaults), 5)

    log('PLAMS environment initialized', 5)
    log('PLAMS working folder: {}'.format(config.default_jobmanager.workdir), 1)

    try:
        import dill
    except ImportError:
        log('WARNING: importing dill package failed. Falling back to the default pickle module. Expect problems with pickling', 1)

    config.init = True


#===========================================================================


def finish(otherJM=None):
    """Wait for all threads to finish and clean the environment.

    This function must be called at the end of your script for |cleaning| to take place. See |master-script| for details.

    If you used some other job managers than just the default one, they need to be passed as *otherJM* list.
    """
    if not config.init:
        return

    for thread in threading.enumerate():
        if thread.name == 'plamsthread':
            thread.join()

    config.default_jobmanager._clean()
    if otherJM:
        for jm in otherJM:
            jm._clean()
    log('PLAMS environment cleaned up successfully', 5)
    log('PLAMS run finished. Goodbye', 3)

    if config.erase_workdir is True:
        shutil.rmtree(config.default_jobmanager.workdir)

    config.init = False


#===========================================================================


def load(filename):
    """Load previously saved job from ``.dill`` file. This is just a shortcut for |load_job| method of the default |JobManager| ``config.default_jobmanager``."""
    return config.default_jobmanager.load_job(filename)


#===========================================================================


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
    for foldername in filter(lambda x: isdir(opj(path,x)), os.listdir(path)):
        maybedill = opj(path,foldername,foldername+'.dill')
        if isfile(maybedill):
            job = jm.load_job(maybedill)
            if job:
                loaded_jobs[os.path.abspath(maybedill)] = job
        else:
            loaded_jobs.update(load_all(path=opj(path,foldername), jobmanager=jm))
    return loaded_jobs


#===========================================================================


@retry()
def delete_job(job):
    """Remove *job* from its corresponding |JobManager| and delete the job folder from the disk. Mark *job* as 'deleted'."""

    if job.status != 'created':
        job.results.wait()

    #In case job.jobmanager is None, run() method was not called yet, so no JobManager knows about this job and no folder exists.
    if job.jobmanager is not None:
        job.jobmanager.remove_job(job)

    if job.parent is not None:
        job.parent.remove_child(job)

    job.status = 'deleted'
    job._log_status(5)


#===========================================================================


def read_molecules(folder, formats=None):
    """Read all molecules from *folder*.

    Read all the files present in *folder* with extensions compatible with :meth:`Molecule.read<scm.plams.mol.molecule.Molecule.read>`. Returned value is a dictionary with keys being molecule names (filename without extension) and values being |Molecule| instances.

    The optional argument *formats* can be used to narrow down the search to files with specified extensions::

        molecules = read_molecules('mymols', formats=['xyz', 'pdb'])
    """
    from ..mol.molecule import Molecule
    extensions = formats or list(Molecule._readformat.keys())
    is_valid = lambda x: isfile(opj(folder,x)) and any([x.endswith('.'+ext) for ext in extensions])
    filenames = filter(is_valid, os.listdir(folder))
    ret = {}
    for f in filenames:
        m = Molecule(opj(folder,f))
        ret[m.properties.name] = m
    return ret


#===========================================================================


def read_all_molecules_in_xyz_file(filename):
    """The .xyz format allows to store multiple geometries on a single file (such file is essentially a concatenated series of xyz files)

    This function returns a *list* of all molecules found in the file *filename*

    *filename*: path (absolute or relative to the current working directory) to the xyz file.
    """
    from ..mol.molecule import Molecule
    mols = []
    with open(filename, 'r') as f:
        while True:
            try:
                mol = Molecule()
                mol.readxyz(f)
                mols.append(mol)
            except FileError:
                break
    return mols


#===========================================================================


_stdlock = threading.Lock()
_filelock = threading.Lock()

def log(message, level=0):
    """Log *message* with verbosity *level*.

    Logs are printed independently to the text logfile (a file called ``logfile`` in the main working folder) and to the standard output. If *level* is equal or lower than verbosity (defined by ``config.log.file`` or ``config.log.stdout``) the message is printed. Date and/or time can be added based on ``config.log.date`` and ``config.log.time``. All logging activity is thread safe.
    """
    if config.init and 'log' in config:
        if level <= config.log.file or level <= config.log.stdout:
            message = str(message)
            prefix = ''
            if config.log.date:
                prefix += '%d.%m|'
            if config.log.time:
                prefix += '%H:%M:%S'
            if prefix:
                prefix = '[' + prefix.rstrip('|') + '] '
                message = time.strftime(prefix) + message
            if level <= config.log.stdout:
                with _stdlock:
                    print(message)
            if level <= config.log.file and 'default_jobmanager' in config:
                with _filelock, open(config.default_jobmanager.logfile, 'a') as f:
                    f.write(message + '\n')
    elif level <= 3:
        # log() is called before plams.init() was called ...
        with _stdlock: print(message)


#===========================================================================


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
    from .results import _restrict, _MetaResults
    def decorator(func):
        if isinstance(classname, _MetaResults):
            func = _restrict(func)
        setattr(classname, func.__name__, func)
    return decorator


#===========================================================================


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
    from .results import _restrict, Results
    def decorator(func):
        if isinstance(instance, Results):
            func = _restrict(func)
        func = types.MethodType(func, instance)
        setattr(instance, func.__func__.__name__, func)
    return decorator


#===========================================================================


def parse_heredoc(bash_input: str, heredoc_delimit: str = 'eor') -> str:
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
    start_pattern = r'<<(-)?(\s+)?{}'.format(heredoc_delimit)
    start_heredoc = re.search(start_pattern, bash_input)
    if not start_heredoc:
        return bash_input

    # Find the end of the heredoc block
    end_pattern = r'\n(\s+)?{}(\s+)?\n'.format(heredoc_delimit)
    end_heredoc = re.search(end_pattern, bash_input)

    # Prepare the slices
    try:
        i, j = start_heredoc.end(), end_heredoc.start()
    except AttributeError as ex:
        err = f"parse_heredoc: failed to find the final '{heredoc_delimit}' delimiter"
        raise ValueError(err).with_traceback(ex.__traceback__)

    # Grab heredoced block and parse it
    _, ret = bash_input[i:j].split('\n', maxsplit=1)
    return ret

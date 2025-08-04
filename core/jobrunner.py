import functools
import os
import threading
import time
from os.path import join as opj
from subprocess import DEVNULL, PIPE

from scm.plams.core.errors import PlamsError
from scm.plams.core.functions import config, log
from scm.plams.core.private import saferun
from scm.plams.core.settings import Settings
from scm.plams.core.threading_utils import LimitedSemaphore

__all__ = ["JobRunner", "GridRunner"]


def _in_thread(func):
    """Decorator for an instance method. If ``parallel`` attribute of given instance is ``True``, run decorated method in a separate :class:`~threading.Thread`. This thread is usually a daemon thread, the decision is based on ``config.daemon_threads`` entry."""

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self.parallel:
            t = threading.Thread(name="plamsthread", target=func, args=(self,) + args, kwargs=kwargs)
            t.daemon = config.daemon_threads
            t.start()
        else:
            func(self, *args, **kwargs)

    return wrapper


def _in_limited_thread(func):
    """More careful version of the ``in_thread`` decorator: a new job thread will only be launched if the currently active number of threads is below the ``maxthreads`` limit set in :meth:`~scm.plams.core.jobrunner.JobRunner.__init__`."""

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self.parallel:
            t = threading.Thread(name="plamsthread", target=func, args=(self,) + args, kwargs=kwargs)
            t.daemon = config.daemon_threads
            if self._jobthread_limit:
                self._jobthread_limit.acquire()
            try:
                t.start()
            except Exception as e:
                if self._jobthread_limit:
                    self._jobthread_limit.release()
                raise e
        else:
            func(self, *args, **kwargs)

    return wrapper


def _limit(func):
    """Decorator for an instance method. If ``_job_limit`` attribute of given instance is not ``None``, use this attribute to wrap decorated method via :ref:`with<with-locks>` statement."""

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self._job_limit:
            with self._job_limit:
                return func(self, *args, **kwargs)
        else:
            return func(self, *args, **kwargs)

    return wrapper


class _MetaRunner(type):
    """Metaclass for |JobRunner|. During an instance creation wrap the :meth:`~scm.plams.core.jobrunner.JobRunner.call` method with :func:`_limit` decorator which enforces a limit on the number of simultaneous :meth:`~scm.plams.core.jobrunner.JobRunner.call` calls."""

    def __new__(meta, name, bases, dct):
        dct["call"] = _limit(dct["call"])
        return type.__new__(meta, name, bases, dct)


# ===========================================================================
# ===========================================================================
# ===========================================================================


class JobRunner(metaclass=_MetaRunner):
    """Class defining the basic job runner interface. Instances of this class represent local job runners -- job runners that execute computational jobs on the current machine.

    The goal of the job runner is to take care of two important things -- parallelization and runscript execution:

    *   When the method |run| of any |Job| instance is executed, the control, after some initial preparations, is passed to a |JobRunner| instance. This |JobRunner| instance decides if a separate thread should be spawned for the job or if the execution should proceed in the current thread. This decision is based on the ``parallel`` attribute which can be set on |JobRunner| creation. There are no separate classes for serial and parallel job runner, both cases are covered by |JobRunner| depending on the ``parallel`` parameter.
    *   If the executed job is an instance of |SingleJob|, it creates a runscript which contains most of the actual computational work (usually it's just an execution of some external binary). The runscript is then submitted to a |JobRunner| instance using its :meth:`call` method. This method executes the runscript as a separate subprocess and takes care of output and error streams handling, setting a proper working directory etc.

    For a job runner with parallel execution enabled the number of simultaneously running jobs can be limited using the *maxjobs* parameter. If *maxjobs* is 0, no limit is enforced. If *parallel* is ``False``, *maxjobs* is ignored. If *parallel* is ``True`` and *maxjobs* is a positive integer, a :class:`BoundedSemaphore<threading.BoundedSemaphore>` of that size is used to limit the number of simultaneously running :meth:`call` methods.

    For a parallel |JobRunner| the maximum amount of threads started is limited by the *maxthreads* argument. As each job is running in a separate thread, this number necessarily acts as an upper limit for the number of jobs that can run in parallel. If the limit is exhausted, running further jobs with this |JobRunner| instance will block execution until another already running job thread terminates. The ``maxthreads`` limit should be set as large as possible, but not so large as to exceed any limits imposed by the operating system.

    A |JobRunner| instance can be passed to |run| with a keyword argument ``jobrunner``. If this argument is omitted, the instance stored in ``config.default_jobrunner`` is used.
    """

    def __init__(self, parallel: bool = False, maxjobs: int = 0, maxthreads: int = 256):
        # Initialize protected attributes
        self._parallel = parallel
        self._job_limit = None
        self._jobthread_limit = None

        # Set properties
        self.maxjobs = maxjobs
        self.maxthreads = maxthreads

    @property
    def parallel(self) -> bool:
        """
        Whether jobs can be run simultaneously or not. Defaults to ``False``.

        When set to ``True``, jobs can be run in parallel, with the level of parallelization determined by the properties :attr:`~scm.plams.core.jobrunner.JobRunner.maxjobs` and :attr:`~scm.plams.core.jobrunner.JobRunner.maxthreads`.
        """
        return self._parallel

    @parallel.setter
    def parallel(self, value: bool) -> None:
        self._parallel = value

    @property
    def maxjobs(self) -> int:
        """
        Maximum number of jobs which can be run in parallel. Defaults to ``0``.

        When set to ``0``, no job limit is enforced.

        Note that to take effect, this requires :attr:`~scm.plams.core.jobrunner.JobRunner.parallel` to be ``True``.
        The number of jobs may also be limited by :attr:`~scm.plams.core.jobrunner.JobRunner.maxthreads`.

        When jobs are running using this, updating this value will block until its usage is complete.
        """
        return self._job_limit.max_value if self._job_limit else 0

    @maxjobs.setter
    def maxjobs(self, value: int) -> None:
        if value < 0:
            raise ValueError(f"Value of 'maxjobs' must be greater or equal to zero, but was {value}")

        if self._job_limit:
            # safely update maximum value (this will block on decrease)
            self._job_limit.max_value = value
            if value == 0:
                self._job_limit = None
        else:
            self._job_limit = LimitedSemaphore(value) if value else None

    @property
    def maxthreads(self) -> int:
        """
        Maximum number of threads which can be started in parallel. Defaults to ``256``.

        When set to ``0``, no thread limit is enforced.

        As each job is runs in a separate thread, this number necessarily acts as an upper limit for the number of jobs that can run in parallel.
        If the limit is exhausted, running further jobs with this |JobRunner| instance will block execution until another already running job thread terminates.
        The ``maxthreads`` limit should be set as large as possible, but not so large as to exceed any limits imposed by the operating system.

        Note that to take effect, this requires :attr:`~scm.plams.core.jobrunner.JobRunner.parallel` to be ``True``.
        The number of jobs may also be limited by :attr:`~scm.plams.core.jobrunner.JobRunner.maxjobs`.

        When jobs are running using this, updating this value will block until its usage is complete.
        """
        return self._jobthread_limit.max_value if self._jobthread_limit else 0

    @maxthreads.setter
    def maxthreads(self, value: int) -> None:
        if value < 0:
            raise ValueError(f"Value of 'maxthreads' must be greater or equal to zero, but was {value}")

        if self._jobthread_limit:
            # safely update maximum value (this will block on decrease)
            self._jobthread_limit.max_value = value
            if value == 0:
                self._jobthread_limit = None
        else:
            self._jobthread_limit = LimitedSemaphore(value) if value else None

    def call(self, runscript, workdir, out, err, runflags):
        """call(runscript, workdir, out, err, runflags)
        Execute the *runscript* in the folder *workdir*. Redirect output and error streams to *out* and *err*, respectively.

        Arguments *runscript*, *workdir*, *out* and *err* should be strings with paths to corresponding files or folders.

        *runflags* is a |Settings| instance containing the `run` branch of running job's settings. The basic job runner defined here ignores them, but they can be useful in |JobRunner| subclasses (see :meth:`GridRunner.call`).

        Returns an integer with the exit code returned by the *runscript*.

        This method can be safely overridden in |JobRunner| subclasses. For example, in |GridRunner| it submits the runscript to a queueing system instead of executing it locally.

        .. note::
            This method is used automatically during |run| and should never be explicitly called in your script.
        """
        log("Executing {}".format(runscript), 5)
        command = ["./" + runscript] if os.name == "posix" else ["sh", runscript]
        if out is not None:
            with open(opj(workdir, err), "w") as e, open(opj(workdir, out), "w") as o:
                process = saferun(command, cwd=workdir, stderr=e, stdout=o)
        else:
            with open(opj(workdir, err), "w") as e:
                process = saferun(command, cwd=workdir, stderr=e)
        log("Execution of {} finished with returncode {}".format(runscript, process.returncode), 5)
        return process.returncode

    @_in_limited_thread
    def _run_job(self, job, jobmanager):
        """_run_job(job, jobmanager)
        This method aggregates the parts of :ref:`job-life-cycle` that are supposed to be run in a separate thread in case of parallel job execution. It is wrapped with :func:`_in_limited_thread` decorator.

        This method should not be overridden.
        """
        try:
            if job._prepare(jobmanager):
                job._execute(self)
                job._finalize()
        finally:
            # Log job summaries
            try:
                # Log any error messages to the standard logger
                if not job.ok(False) or not job.check():
                    err_msg = job.get_errormsg()
                    err_lines = err_msg.splitlines()
                    max_lines = 30
                    if len(err_lines) > max_lines:
                        err_lines = err_lines[:max_lines]
                        err_lines.append("... (see output for full error)")
                    err_msg = str.join("\n", [f"\t{l}" for l in err_lines])
                    log(f"""Error message for job {job.name} was:\n{err_msg}""", 3)

                # Log job summary to the csv logger
                jobmanager.job_logger.log(job, level=3)
            except:  # logging should never throw, but best to make sure
                pass

            if self.parallel and self._jobthread_limit:
                self._jobthread_limit.release()


# ===========================================================================
# ===========================================================================
# ===========================================================================


class GridRunner(JobRunner):
    """Subclass of |JobRunner| that submits the runscript to a queueing system instead of executing it locally. Besides two new keyword arguments (*grid* and *sleepstep*) it behaves and is meant to be used just like a regular |JobRunner|.

    .. note::

        The default value of the *parallel* argument is ``True``, unlike in the regular |JobRunner|.

    There are many different queueing systems that are popular nowadays (for example: TORQUE, SLURM, OGE). Usually they use different commands for submitting jobs or checking the queue status. |GridRunner| class tries to build a common interface to these systems. The commands used to communicate with the queueing system are not hard-coded, but rather taken from a |Settings| instance. Thanks to that the user has almost full control over the behavior of a |GridRunner| instance and the behavior can be ajdusted dynamically.

    The behavior of a |GridRunner| instance is determined by the contents of a |Settings| instance stored in its ``settings`` attribute. That |Settings| instance can be manually supplied by the user or taken from a collection of predefined instances stored as branches of ``GridRunner.config``. The adjustment is done with the *grid* parameter which should be a string or a |Settings| instance. If it's a string, it has to be a key occurring in ``GridRunner.config`` or ``'auto'`` for autodetection. For example, if ``grid='slurm'`` is passed, ``GridRunner.config.slurm`` is used as settings. If ``grid='auto'`` then entries present in ``GridRunner.config`` are tested and the first one that works (its submit command is present on your system) is chosen. When a |Settings| instance is passed as *grid*, it is directly used as ``settings``.

    Currently two predefined schemes are available (see :ref:`global-settings`): ``slurm`` for SLURM and ``pbs`` for queueing systems following PBS syntax (PBS, TORQUE, Oracle Grid Engine etc.).

    The settings of |GridRunner| should have the following structure:

    *   ``output`` -- flag for specifying the output file path.
    *   ``error`` -- flag for specifying the error file path.
    *   ``workdir`` -- flag for specifying path to the working directory.
    *   ``commands.submit`` -- submit command.
    *   ``commands.check`` -- queue status check command.
    *   ``commands.getid`` -- function extracting submitted job's ID from the output of the submit command.
    *   ``commands.running`` -- function extracting a list of all running jobs from the output of queue check command
    *   ``commands.special`` -- branch storing definitions of special |run| keyword arguments.

    See :meth:`call` for more details and examples.

    The *sleepstep* parameter defines how often the queue check is performed. It should be a numerical value telling how many seconds should the interval between two consecutive checks last.

    .. note::
        Usually queueing systems are configured in such a way that output of your calculation is captured somewhere else and copied to the location indicated by the output flag only when the job is finished. Because of that it is not possible to have a peek at your output while your job is running (for example, to see if your calculation is going well). This limitation can be circumvented with ``myjob.settings.runscript.stdout_redirect`` flag. If set to ``True``, the output redirection will not be handled by the queueing system, but rather placed in the runscript using the shell redirection ``>``. That forces the output file to be created directly in *workdir* and updated live as the job proceeds.
    """

    # GridRunner mechanism for testing if a job is finished:
    # if [...].commands.finished exists it is used to check if the job is finished. It should be a function that takes a single string (job_id) as an argument and returns True or False
    # otherwise [...].commands.check is combined with job_id, executed as a subprocess and returned exit code is tested (nonzero return code indicates that job has finished)

    # N.B. these SLURM/PBS functions are always used in a static way
    # but as stored and accessed via Settings, cannot be decorated as such
    def __slurm_get_jobid(output: str):  # type: ignore
        s = output.split()
        if len(s) > 0 and all([ch.isdigit() for ch in s[-1]]):
            return s[-1]
        return None

    def __slurm_running(output: str):  # type: ignore
        lines = output.splitlines()[1:]
        return [line.split()[0] for line in lines]

    def __pbs_get_jobid(output: str):  # type: ignore
        s = output.split(".")
        if len(s) > 0 and all([ch.isdigit() for ch in s[0]]):
            return s[0]
        return None

    def __pbs_running(output: str):  # type: ignore
        lines = output.splitlines()[2:]
        return [line.split()[0].split(".")[0] for line in lines]

    # Static config holding the preconfigured grids:
    config = Settings()
    # PBS
    config.pbs.workdir = "-d"
    config.pbs.output = "-o"
    config.pbs.error = "-e"
    config.pbs.special.nodes = "-l nodes="
    config.pbs.special.walltime = "-l walltime="
    config.pbs.special.memory = "-l mem="
    config.pbs.special.queue = "-q "
    config.pbs.commands.submit = "qsub"
    config.pbs.commands.check = "qstat"
    config.pbs.commands.getid = __pbs_get_jobid
    config.pbs.commands.running = __pbs_running
    # Slurm
    config.slurm.workdir = "-D"
    config.slurm.output = "-o"
    config.slurm.error = "-e"
    config.slurm.special.nodes = "-N "
    config.slurm.special.cores = "-n "
    config.slurm.special.walltime = "-t "
    config.slurm.special.memory = "--mem="
    config.slurm.special.queue = "-p "
    config.slurm.commands.submit = "sbatch"
    config.slurm.commands.check = "squeue"
    config.slurm.commands.getid = __slurm_get_jobid
    config.slurm.commands.running = __slurm_running

    def __init__(self, grid="auto", sleepstep=5, parallel=True, maxjobs=0):
        JobRunner.__init__(self, parallel=parallel, maxjobs=maxjobs)
        self.sleepstep = sleepstep
        self._active_jobs = {}
        self._active_lock = threading.Lock()
        self._mainlock = threading.Lock()

        if isinstance(grid, Settings):
            self.settings = grid
        elif grid == "auto":
            self.settings = self._autodetect()
        elif grid in GridRunner.config:
            self.settings = GridRunner.config[grid]
            try:
                saferun([self.settings.commands.submit, "--version"], stdout=DEVNULL, stderr=DEVNULL)
            except OSError:
                raise PlamsError("GridRunner: {} command not found".format(self.settings.commands.submit))
        else:
            raise PlamsError(
                "GridRunner: invalid 'grid' argument. 'grid' should be either a Settings instance (see documentations for details) or a string occurring in GridRunner.config or 'auto' for autodetection"
            )

    def call(self, runscript, workdir, out, err, runflags):
        """call(runscript, workdir, out, err, runflags)
        Submit *runscript* to the queueing system with *workdir* as the working directory. Redirect output and error streams to *out* and *err*, respectively. *runflags* stores varoius submit command options.

        The submit command has the following structure::

            <commands.submit>_<workdir>_{workdir}_<error>_{err}[_<output>_{out}][FLAGS]_{runscript}

        Underscores denote spaces, parts in pointy brackets correspond to ``settings`` entries, parts in curly brackets to :meth:`call` arguments, square brackets contain optional parts. Output part is added if *out* is not ``None``. This is handled automatically based on ``runscript.stdout_redirect`` value in job's ``settings``.

        ``FLAGS`` part is built based on *runflags* argument, which is a |Settings| instance storing |run| keyword arguments. For every *(key,value)* pair in *runflags* the string ``_-key_value`` is appended to ``FLAGS`` **unless** the *key* is a special key occurring in ``commands.special``. In that case ``_<commands.special.key>value`` is used (mind the lack of space in between). For example, a |Settings| instance defining interaction with SLURM has the following entries::

            workdir = '-D'
            output  = '-o'
            error   = '-e'
            special.nodes    = '-N '
            special.walltime = '-t '
            special.memory = '--mem='
            special.queue    = '-p '
            commands.submit  = 'sbatch'
            commands.check  = 'squeue'

        The submit command produced by::

            gr = GridRunner(parallel=True, maxjobs=4, grid='slurm')
            j.run(jobrunner=gr, queue='short', nodes=2, J='something', O='')

        will be:

        .. code-block:: none

            sbatch -D {workdir} -e {err} -o {out} -p short -N 2 -J something -O  {runscript}

        In certain queueing systems some flags don't have a short form with semantics ``-key value``. For example, in SLURM the flag ``--nodefile=value`` has a short form ``-F value``, but the flag ``--export=value`` does not. One can still use such a flag using the special keys logic::

            gr = GridRunner(parallel=True, maxjobs=4, grid='slurm')
            gr.settings.special.export = '--export='
            j.run(jobrunner=gr, queue='short', export='value')

        That results in the command::

            sbatch -D {workdir} -e {err} -o {out} -p short --export=value {runscript}

        The submit command is then executed and the output returned by it is used to determine the submitted job's ID. The value stored in ``commands.getid`` is used for that purpose. It should be a function taking a single string (the whole output of the submit command) and returning a string with job's ID.

        The submitted job's ID is then added to ``_active_jobs`` dictionary, with the key being job's ID and the value being an instance of :class:`threading.Lock`. This lock is used to singal the fact that the job is finished and the thread handling it can continue. Then the :meth:`_check_queue` method starts the thread querying the queue and unlocking finished jobs.

        Since it is difficult to automatically obtain job's exit code, the returned value is 0 (or 1, if the submit command failed). From |run| perspective it means that a job executed with |GridRunner| is *crashed* only if it never entered the queue (usually due to improper submit command).

        .. note::
            This method is used automatically during |run| and should never be explicitly called in your script.
        """
        s = self.settings
        cmd = " ".join([s.commands.submit, s.workdir, workdir, s.error, err])
        if out is not None:
            cmd += " " + s.output + " " + out
        for k, v in runflags.items():
            if k in s.special:
                cmd += " " + s.special[k] + str(v)
            else:
                cmd += " -" + k + " " + str(v)
        cmd += " " + opj(workdir, runscript)

        log("Submitting {} with command {}".format(runscript, cmd), 5)
        process = saferun(cmd.split(" "), stdout=PIPE, stderr=PIPE)
        subout = process.stdout.decode()
        log("Output of {} submit command: {}".format(runscript, subout), 5)

        jobid = s.commands.getid(subout)
        if jobid is None:
            log("Submitting of {} failed. Stderr of submit command:\n{}".format(runscript, process.stderr.decode()), 1)
            return 1
        log("{} submitted successfully as job {}".format(runscript, jobid), 3)

        event = threading.Event()
        with self._active_lock:
            self._active_jobs[jobid] = event
        self._check_queue()
        event.wait()

        log("Execution of {} finished".format(runscript), 5)
        return 0

    @_in_thread
    def _check_queue(self):
        """Query the queueing system to obtain a list of currently running jobs. Check for active jobs that are not any more in the queue and release their locks. Repeat this procedure every ``sleepstep`` seconds until there are no more active jobs. The ``_mainlock`` lock ensures that there is at most one thread executing the main loop of this method at the same time."""
        if self._mainlock.acquire(blocking=False):
            try:
                while True:
                    with self._active_lock:
                        active_jobs = set(self._active_jobs.keys())

                    process = saferun([self.settings.commands.check], stdout=PIPE)
                    output = process.stdout.decode()
                    running_jobs = set(self.settings.commands.running(output))

                    with self._active_lock:
                        for jobid in active_jobs - running_jobs:
                            self._active_jobs[jobid].set()
                            del self._active_jobs[jobid]
                        if len(self._active_jobs) == 0:
                            return
                    time.sleep(self.sleepstep)
            finally:
                self._mainlock.release()

    def _autodetect(self):
        """Try to autodetect the type of queueing system.

        The autodetection mechanism is very simple. For each entry in ``GridRunner.config`` the submit command followed by ``--version`` is executed (for example ``qsub --version``). If the execution was successful (which is indicated by the exit code 0), that queueing system is present and it is chosen. Thus if there are multiple queueing systems installed, only one of them is picked -- the one which "name" (indicated by a key in ``GridRunner.config``) is first in the lexicographical order.

        Returned value is one of ``GridRunner.config`` branches. If autodetection was not successful, an exception is raised.
        """
        for grid in GridRunner.config:
            try:
                process = saferun(
                    [GridRunner.config[grid].commands.submit, "--version"], stdout=DEVNULL, stderr=DEVNULL
                )
            except OSError:
                continue
            if process.returncode == 0:
                log("Grid type autodetected as '{}'".format(grid), 5)
                return GridRunner.config[grid]
        raise PlamsError("GridRunner: Failed to autodetect grid type")

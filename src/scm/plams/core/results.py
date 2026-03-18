import copy
import functools
import glob
import inspect
import operator
import os
import shutil
import threading
from os.path import join as opj
from subprocess import PIPE
from typing import List, Dict, TYPE_CHECKING, Optional, Callable, Union, TypeVar, Tuple, Any
from typing_extensions import ParamSpec, Concatenate
from types import FrameType

from scm.plams.core.errors import FileError, ResultsError
from scm.plams.core.functions import get_config, log
from scm.plams.core.private import saferun
from scm.plams.core.enums import JobStatus

if TYPE_CHECKING:
    from scm.plams.core.basejob import Job
    from scm.plams.core.settings import Settings
    from scm.plams.mol.molecule import Molecule

T = TypeVar("T")
P = ParamSpec("P")
TSelf = TypeVar("TSelf", bound="Results")

__all__ = ["Results"]


def _caller_name_and_arg(frame: Optional[FrameType]) -> Tuple[Optional[str], Optional[Any]]:
    """Extract information about name and arguments of a function call from a *frame* object"""
    if frame is None:
        return None, None
    caller_name = frame.f_code.co_name
    caller_varnames = frame.f_code.co_varnames
    caller_arg = None
    if len(caller_varnames) > 0:
        try:
            loc = frame.f_locals
        except:
            loc = {}
        if caller_varnames[0] in loc:
            caller_arg = loc[caller_varnames[0]]
    return caller_name, caller_arg


def _privileged_access() -> bool:
    """Analyze contents of the current stack to find out if privileged access to the |Results| methods should be granted.

    Privileged access is granted to two |Job| methods: |postrun| and :meth:`~scm.plams.core.basejob.Job.check`, but only if they are called from :meth:`~scm.plams.core.basejob.Job._finalize` of the same |Job| instance.
    """
    from scm.plams.core.basejob import Job

    for frame in inspect.getouterframes(inspect.currentframe()):
        cal, arg = _caller_name_and_arg(frame[0])
        prev_cal, prev_arg = _caller_name_and_arg(frame[0].f_back)
        if cal in ["postrun", "check"] and prev_cal == "_finalize" and arg == prev_arg and isinstance(arg, Job):
            return True
    return False


def _restrict(func: Callable[Concatenate["Results", P], T]) -> Callable[Concatenate["Results", P], Optional[T]]:
    """Decorator that wraps methods of |Results| instances.

    Whenever decorated method is called, the status of associated job is checked. Depending on its value access to the method is granted, refused or the calling thread is forced to wait for the right :ref:`event<event-objects>` to be set.
    """

    @functools.wraps(func)
    def guardian(self: "Results", /, *args: P.args, **kwargs: P.kwargs) -> Optional[T]:
        if not self.job:
            raise ResultsError("Using Results not associated with any Job")

        cfg = get_config()
        if self.job.status in [JobStatus.SUCCESSFUL, JobStatus.COPIED]:
            return func(self, *args, **kwargs)

        elif self.job.status == JobStatus.PREVIEW:
            if cfg.ignore_failure:
                log(
                    f"WARNING: Trying to obtain results of job {self.job.name} run in a preview mode. Returned value is None",
                    3,
                )
                return None
            else:
                raise ResultsError("Using Results associated with job run in a preview mode")

        elif self.job.status == JobStatus.DELETED:
            raise ResultsError("Using Results associated with deleted job")

        elif self.job.status in [JobStatus.CRASHED, JobStatus.FAILED]:
            # waiting for crashed of failed job should not trigger any warnings/exceptions
            # and neither should checking the status from this job with 'ok', 'check' or 'get_errormsg'
            suppress_errors = func.__name__ == "wait"
            if not suppress_errors:
                for frame in inspect.getouterframes(inspect.currentframe()):
                    cal, arg = _caller_name_and_arg(frame[0])
                    if arg == self.job and cal in ["ok", "check", "get_errormsg"]:
                        suppress_errors = True
                        break

            if suppress_errors:
                cal, arg = _caller_name_and_arg(inspect.currentframe())
                if isinstance(arg, Results):
                    return func(self, *args, **kwargs)

            if cfg.ignore_failure:
                log(f"WARNING: Trying to obtain results of crashed or failed job {self.job.name}", 3)
                try:
                    ret = func(self, *args, **kwargs)
                except:
                    log(f"Obtaining results of {self.job.name} failed. Returned value is None", 3)
                    return None
                log(
                    f"Obtaining results of {self.job.name} successful. However, no guarantee that they make sense",
                    3,
                )
                return ret
            else:
                raise ResultsError("Using Results associated with crashed or failed job")

        elif self.job.status in [JobStatus.CREATED, JobStatus.STARTED, JobStatus.REGISTERED, JobStatus.RUNNING]:
            log(f"Waiting for job {self.job.name} to finish", 1)
            if _privileged_access():
                self.finished.wait()
            else:
                self.done.wait()
            return func(self, *args, **kwargs)

        elif self.job.status == JobStatus.FINISHED:
            if _privileged_access():
                return func(self, *args, **kwargs)
            log(f"Waiting for job {self.job.name} to finish", 1)
            self.done.wait()
            return func(self, *args, **kwargs)

        return None

    return guardian


# ===========================================================================
# ===========================================================================
# ===========================================================================


class ApplyRestrict:
    """Parent class that wraps all methods with :func:`_restrict` decorator when a subclass of it is defined
    Used as parent class for |Results|, ensuring proper synchronization and thread safety.
    Methods listed in ``_dont_restrict`` as well as "magic methods" are not wrapped."""

    def __init_subclass__(cls) -> None:
        _dont_restrict = ["refresh", "collect", "_clean", "get_errormsg", "collect_rkfs"]
        for name, attr in cls.__dict__.items():
            # don't touch magic methods
            if name.startswith("__") and name.endswith("__"):
                continue
            # don't apply decorator to attributes that are classes (copied over, unsure why necessary)
            if type(attr) is type:
                continue
            # the _restrict decorator does not make sense for staticmethods, since they don't get `self` passed as the first argument
            if isinstance(attr, (staticmethod, classmethod)):
                continue
            # Note that in python 3.10 staticmethod objects are considered callable, hence the previous line to be futureproof.
            if callable(attr) and name not in _dont_restrict:
                setattr(cls, name, _restrict(attr))


# ===========================================================================
# ===========================================================================
# ===========================================================================


class Results(ApplyRestrict):
    """General concrete class for job results.

    ``job`` attribute stores a reference to associated job. ``files`` attribute is a list with contents of the job folder. ``_rename_map`` is a class attribute with the dictionary storing the default renaming scheme.

    Bracket notation (``myresults[filename]``) can be used to obtain full absolute paths to files in the job folder.

    Instance methods are automatically wrapped with the "access guardian" that ensures thread safety (see |parallel|).
    """

    _rename_map: Dict[str, str] = {}

    def __init__(self, job: "Job"):
        self.job = job
        self.files: List[str] = []
        self.finished = threading.Event()
        self.done = threading.Event()

    def refresh(self) -> None:
        """Refresh the contents of the ``files`` list. Traverse the job folder (and all its subfolders) and collect relative paths to all files found there, except files with ``.dill`` extension.

        This is a cheap and fast method that should be used every time there is a risk the contents of the job folder changed and ``files`` is no longer up-to-date. For proper working of various PLAMS elements it is crucial that ``files`` always contains up-to-date information about the contents of the job folder.

        All functions and methods defined in PLAMS that could change the state of the job folder refresh the ``files`` list, so there is no need to manually call :meth:`~Results.refresh` after, for example, :meth:`~Results.rename`. If you are implementing a new method of that kind, please don't forget about refreshing.
        """
        if self.job.path is None:
            return
        self.files = []
        for pth, dirs, files in os.walk(self.job.path):
            relpath = os.path.relpath(pth, self.job.path)
            self.files += [opj(relpath, x) if relpath != "." else x for x in files]
        self.files = [x for x in self.files if not x.endswith(".dill")]

    def collect(self) -> None:
        """Collect the files present in the job folder after execution of the job is finished. This method is simply :meth:`~Results.refresh` followed by renaming according to the ``_rename_map``.

        If you wish to override this function, you have to call the parent version at the beginning.
        """
        self.refresh()
        for old, new in self.__class__._rename_map.items():
            old = old.replace("$JN", self.job.name)
            new = new.replace("$JN", self.job.name)
            if old in self.files:
                os.rename(self.job.get_path() / old, self.job.get_path() / new)
                self.files[self.files.index(old)] = new
        self.refresh()

    def wait(self) -> None:
        """wait()
        Wait for associated job to finish.

        .. technical::

            This is **not** an abstract method. It does exactly what it should: nothing. All the work is done by :func:`_restrict` decorator that is wrapped around it.
        """

    def grep_file(self, filename: str, pattern: str = "", options: str = "") -> List[str]:
        """grep_file(filename, pattern='', options='')
        Execute ``grep`` on a file given by *filename* and search for *pattern*.

        Additional ``grep`` flags can be passed with *options*, which should be a single string containing all flags, space separated.

        Returned value is a list of lines (strings). See ``man grep`` for details.
        """
        cmd = ["grep"] + [pattern] + options.split()
        return self._process_file(filename, cmd)

    def grep_output(self, pattern: str = "", options: str = "") -> List[str]:
        """grep_output(pattern='', options='')
        Shortcut for :meth:`~Results.grep_file` on the output file."""
        try:
            output = self.job._filename("out")  # type: ignore[attr-defined]
        except AttributeError:
            raise ResultsError(
                f"Job {self.job.name} does not seem to be an instance of SingleJob, it does not have _filenames dictionary"
            )
        except KeyError:
            raise ResultsError(f"Job {self.job.name} does not have an output")
        return self.grep_file(output, pattern, options)

    def read_file(self, filename: str) -> str:
        """
        Returns the contents of the `filename`,
        where `filename` has to be in `self.files`.
        If `filename` contains the  ``$JN`` string, it will be replaced with
        `self.job.name`.

        .. note::

            For text files only. Reading binary files such as `*.rkf` will result in an error.
        """
        filename = filename.replace("$JN", self.job.name)
        if filename not in self.files:
            raise ResultsError(f"No `{filename}` associated with job `{self.job.name}`")
        with open(self.job.get_path() / filename) as f:
            return f.read()

    def regex_file(self, filename: str, regex: str) -> List:
        """
        Applies a regular expression pattern to the
        output of :meth:`read_file` such that the returned value
        is ``re.findall(regex, read_file(filename))``.
        """
        from re import findall

        txt = self.read_file(filename)
        return findall(regex, txt)

    def awk_file(self, filename: str, script: str = "", progfile: Optional[str] = None, **kwargs: Any) -> List[str]:
        """awk_file(filename, script='', progfile=None, **kwargs)
        Execute an AWK script on a file given by *filename*.

        The AWK script can be supplied in two ways: either by directly passing the contents of the script (should be a single string) as the *script* argument, or by providing the path (absolute or relative to *filename*) to a text file with an AWK script as the *progfile* argument. If *progfile* is not ``None``, *script* is ignored.

        Other keyword arguments (*\\*\\*kwargs*) can be used to pass additional variables to AWK (see ``-v`` flag in AWK manual)

        Returned value is a list of lines (strings). See ``man awk`` for details.
        """
        cmd = ["awk"]
        for k, v in kwargs.items():
            cmd += ["-v", f"{k}={v}"]
        if progfile:
            if os.path.isfile(progfile):
                cmd += ["-f", progfile]
            else:
                raise FileError(f"File {progfile} not present")
        else:
            cmd += [script]
        return self._process_file(filename, cmd)

    def awk_output(self, script: str = "", progfile: Optional[str] = None, **kwargs: Any) -> List[str]:
        """awk_output(script='', progfile=None, **kwargs)
        Shortcut for :meth:`~Results.awk_file` on the output file."""
        try:
            output = self.job._filename("out")  # type: ignore[attr-defined]
        except AttributeError:
            raise ResultsError(
                f"Job {self.job.name} does not seem to be an instance of SingleJob, it does not have _filenames dictionary"
            )
        except KeyError:
            raise ResultsError(f"Job {self.job.name} does not have an output")
        return self.awk_file(output, script, progfile, **kwargs)

    def rename(self, old: str, new: str) -> None:
        """rename(old, new)
        Rename a file from ``files``. In both *old* and *new* the shortcut ``$JN`` for job name can be used."""
        old = old.replace("$JN", self.job.name)
        new = new.replace("$JN", self.job.name)
        self.refresh()
        if old in self.files:
            os.rename(self.job.get_path() / old, self.job.get_path() / new)
            self.files[self.files.index(old)] = new
        else:
            raise FileError(f"File {old} not present in {str(self.job.get_path())}")

    def get_file_chunk(
        self,
        filename: str,
        begin: Optional[str] = None,
        end: Optional[str] = None,
        match: int = 0,
        inc_begin: bool = False,
        inc_end: bool = False,
        process: Optional[Callable[[str], T]] = None,
    ) -> Union[List[str], List[T]]:
        """get_file_chunk(filename, begin=None, end=None, match=0, inc_begin=False, inc_end=False, process=None)

        Extract a chunk of a text file given by *filename*, consisting of all the lines between a line containing *begin* and a line containing *end*.

        *begin* and *end* should be simple strings (no regular expressions allowed) or ``None`` (in that case matching is done from the beginning or until the end of the file). If multiple blocks delimited by *begin* end *end* are present in the file, *match* can be used to indicate which one should be printed (*match*=0 prints all of them). *inc_begin* and *inc_end* can be used to include the delimiting lines in the final result (by default they are excluded).

        The returned value is a list of strings. *process* can be used to provide a function executed on each element of this list before returning it.
        """
        current_match = 0
        ret: List[str] = []
        switch = begin is None

        def append(x: str) -> None:
            if match in [0, current_match]:
                ret.append(x.rstrip("\n"))

        with open(self[filename], "r") as f:
            for line in f:
                if switch and end and (end in line):
                    switch = False
                    if inc_end:
                        append(line)
                    if match == current_match:
                        break
                if switch:
                    append(line)
                if (not switch) and begin and (begin in line):
                    switch = True
                    current_match += 1
                    if inc_begin:
                        append(line)

        return list(map(process, ret)) if process else ret

    def get_output_chunk(
        self,
        begin: Optional[str] = None,
        end: Optional[str] = None,
        match: int = 0,
        inc_begin: bool = False,
        inc_end: bool = False,
        process: Optional[Callable[[str], T]] = None,
    ) -> Union[List[str], List[T]]:
        """get_output_chunk(begin=None, end=None, match=0, inc_begin=False, inc_end=False, process=None)
        Shortcut for :meth:`~Results.get_file_chunk` on the output file."""
        try:
            output = self.job._filename("out")  # type: ignore[attr-defined]
        except AttributeError:
            raise ResultsError(f"Job {self.job.name} is not an instance of SingleJob, it does not have an output")
        return self.get_file_chunk(output, begin, end, match, inc_begin, inc_end, process)

    def recreate_molecule(self) -> Union[None, "Molecule", Dict[str, "Molecule"]]:
        """Recreate the input molecule for the corresponding job based on files present in the job folder. This method is used by |load_external|.

        The definition here serves as a default fall-back template preventing |load_external| from crashing when a particular |Results| subclass does not define it's own :meth:`recreate_molecule`.
        """
        return None

    def recreate_settings(self) -> Optional["Settings"]:
        """Recreate the input |Settings| instance for the corresponding job based on files present in the job folder. This method is used by |load_external|.

        The definition here serves as a default fall-back template preventing |load_external| from crashing when a particular |Results| subclass does not define it's own :meth:`recreate_settings`.
        """
        return None

    # =======================================================================

    def _clean(self, arg: Union[str, None, List[str]]) -> None:
        """Clean the job folder. *arg* should be a string or a list of strings. See |cleaning| for details."""
        if arg == "all":
            return

        path = self.job.get_path()
        absfiles = [path / f for f in self.files]
        childnames = [child.name for child in self.job] if hasattr(self.job, "children") else []  # type: ignore[attr-defined]
        if arg in ["none", [], None]:
            for f in absfiles:
                if os.path.isfile(f):
                    os.remove(f)

        elif isinstance(arg, list):
            rev = False
            if arg[0] == "-":
                rev = True
                arg = arg[1:]

            absarg = []
            for i in arg:
                s = i.replace("$JN", self.job.name)
                if s.find("$CH") != -1:
                    absarg += [opj(path, s.replace("$CH", ch)) for ch in childnames]
                else:
                    absarg.append(opj(path, s))

            if absarg:
                absarg = functools.reduce(operator.iadd, map(glob.glob, absarg))

            for f in absfiles:
                if (f in absarg) == rev and os.path.isfile(f):  # type: ignore[comparison-overlap]
                    os.remove(f)
                    log(f"Deleting file {str(f)}", 5)

        else:
            log(f"WARNING: {arg} is not a valid keep/save argument", 3)
        self.refresh()

    def _copy_to(self: TSelf, newresults: TSelf) -> None:
        """_copy_to(newresults)
        Copy these results to *newresults*.

        This method is used when |RPM| discovers an attempt to run a job identical to a previously run job. Instead of the execution, results of the previous job are copied/linked to the new one.

        This method is called from |Results| of the old job and *newresults* should be |Results| of the new job. The goal is to faithfully recreate the state of this |Results| instance in ``newresults``. To achieve that, all the contents of the job folder are copied (or hardlinked, if your platform allows that and ``job.settings.link_files`` is ``True``) to other's job folder. Moreover, all attributes of this |Results| instance (other than ``job`` and ``files``) are exported to *newresults* using :meth:`~Results._export_attribute` method.
        """
        for name in self.files:
            newname = Results._replace_job_name(name, self.job.name, newresults.job.name)
            oldpath = self.job.get_path() / name
            newpath = newresults.job.get_path() / newname
            os.makedirs(os.path.dirname(newpath), exist_ok=True)
            if os.name == "posix" and self.job.settings.link_files is True:  # type: ignore[has-type]
                os.link(oldpath, newpath)
            else:
                shutil.copy(oldpath, newpath)
            newresults.files.append(newname)
        for k, v in self.__dict__.items():
            if k in ["job", "files", "done", "finished"]:
                continue
            newresults.__dict__[k] = self._export_attribute(v, newresults)

    def _export_attribute(self: TSelf, attr: T, other: TSelf) -> T:
        """_export_attribute(attr, other)
        Export this instance's attribute to *other*. This method should be overridden in your |Results| subclass if it has some attributes that are not properly handled by :func:`python3:copy.deepcopy`.

        *other* is the |Results| instance, *attr* is the **value** of the attribute to be copied. See :meth:`SCMJob._export_attribute<scm.plams.interfaces.adfsuite.scmjob.SCMResults._export_attribute>` for an example implementation.
        """
        return copy.deepcopy(attr)

    @staticmethod
    def _replace_job_name(string: str, oldname: str, newname: str) -> str:
        """If *string* starts with *oldname*, maybe followed by some extension, replace *oldname* with *newname*."""
        return string.replace(oldname, newname) if (os.path.splitext(string)[0] == oldname) else string

    # =======================================================================

    def __getitem__(self, name: str) -> str:
        """Magic method to enable bracket notation. Elements from ``files`` can be used to get absolute paths."""
        name = name.replace("$JN", self.job.name)
        if name in self.files:
            return str(self.job.get_path() / name)
        else:
            raise FileError(f"File {name} not present in {str(self.job.get_path())}")

    def __contains__(self, name: str) -> bool:
        """Magic method to enable the Python ``in`` operator notation for checking if a filename with a particular name is present."""
        name = name.replace("$JN", self.job.name)
        return name in self.files

    def _process_file(self, filename: str, command: List[str]) -> List[str]:
        """_process_file(filename, command)
        Skeleton for all file processing methods. Execute *command* (should be a list of strings) on *filename* and return output as a list of lines.
        """
        filename = filename.replace("$JN", self.job.name)
        if filename in self.files:
            process = saferun(command + [filename], cwd=str(self.job.get_path()), stdout=PIPE)
            if process.returncode != 0:
                return []
            ret: List[str] = process.stdout.decode().splitlines()
            return ret
        else:
            raise FileError(f"File {filename} not present in {str(self.job.get_path())}")

import contextlib
import textwrap
from functools import wraps
from typing import TYPE_CHECKING, TypeVar, Union, Tuple, Type

__all__ = [
    "Settings",
    "SafeRunSettings",
    "LogSettings",
    "RunScriptSettings",
    "JobSettings",
    "JobManagerSettings",
    "ConfigSettings",
]

if TYPE_CHECKING:
    from scm.plams.core.jobmanager import JobManager
    from scm.plams.core.jobrunner import JobRunner

TSelf = TypeVar("TSelf", bound="Settings")


class Settings(dict):
    """Automatic multi-level dictionary. Subclass of built-in class :class:`dict`.

    The shortcut dot notation (``s.basis`` instead of ``s['basis']``) can be used for keys that:

    *   are strings
    *   don't contain whitespaces
    *   begin with a letter or an underscore
    *   don't both begin and end with two or more underscores.

    Iteration follows lexicographical order (via :func:`sorted` function)

    Methods for displaying content (:meth:`~object.__str__` and :meth:`~object.__repr__`) are overridden to recursively show nested instances in easy-readable format.

    Regular dictionaries (also multi-level ones) used as values (or passed to the constructor) are automatically transformed to |Settings| instances::

        >>> s = Settings({'a': {1: 'a1', 2: 'a2'}, 'b': {1: 'b1', 2: 'b2'}})
        >>> s.a[3] = {'x': {12: 'q', 34: 'w'}, 'y': 7}
        >>> print(s)
        a:
          1:    a1
          2:    a2
          3:
            x:
              12:   q
              34:   w
            y:  7
        b:
          1:    b1
          2:    b2

    """

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

        for k, v in self.items():
            if isinstance(v, dict) and not isinstance(v, Settings):
                self[k] = Settings(v)
            if isinstance(v, list):
                self[k] = [Settings(i) if (isinstance(i, dict) and not isinstance(i, Settings)) else i for i in v]

    def copy(self: TSelf) -> TSelf:
        """Return a new instance that is a copy of this one. Nested |Settings| instances are copied recursively, not linked.

        In practice this method works as a shallow copy: all "proper values" (leaf nodes) in the returned copy point to the same objects as the original instance (unless they are immutable, like ``int`` or ``tuple``). However, nested |Settings| instances (internal nodes) are copied in a deep-copy fashion. In other words, copying a |Settings| instance creates a brand new "tree skeleton" and populates its leaf nodes with values taken directly from the original instance.

        This behavior is illustrated by the following example::

            >>> s = Settings()
            >>> s.a = 'string'
            >>> s.b = ['l','i','s','t']
            >>> s.x.y = 12
            >>> s.x.z = {'s','e','t'}
            >>> c = s.copy()
            >>> s.a += 'word'
            >>> s.b += [3]
            >>> s.x.u = 'new'
            >>> s.x.y += 10
            >>> s.x.z.add(1)
            >>> print(c)
            a:  string
            b:  ['l', 'i', 's', 't', 3]
            x:
              y:    12
              z:    set([1, 's', 'e', 't'])
            >>> print(s)
            a:  stringword
            b:  ['l', 'i', 's', 't', 3]
            x:
              u:    new
              y:    22
              z:    set([1, 's', 'e', 't'])

        This method is also used when :func:`python3:copy.copy` is called.
        """
        cls = type(self)
        ret = cls()
        for name in self:
            if isinstance(self[name], Settings):
                ret[name] = self[name].copy()
            else:
                ret[name] = self[name]
        return ret

    def soft_update(self: TSelf, other: "Settings") -> TSelf:
        """Update this instance with data from *other*, but do not overwrite existing keys. Nested |Settings| instances are soft-updated recursively.

        In the following example ``s`` and ``o`` are previously prepared |Settings| instances::

            >>> print(s)
            a:  AA
            b:  BB
            x:
              y1:   XY1
              y2:   XY2
            >>> print(o)
            a:  O_AA
            c:  O_CC
            x:
              y1:   O_XY1
              y3:   O_XY3
            >>> s.soft_update(o)
            >>> print(s)
            a:  AA        #original value s.a not overwritten by o.a
            b:  BB
            c:  O_CC
            x:
              y1:   XY1   #original value s.x.y1 not overwritten by o.x.y1
              y2:   XY2
              y3:   O_XY3

        *Other* can also be a regular dictionary. Of course in that case only top level keys are updated.

        Shortcut ``A += B`` can be used instead of ``A.soft_update(B)``.
        """
        for name in other:
            if isinstance(other[name], Settings):
                if name not in self:
                    self[name] = other[name].copy()
                elif isinstance(self[name], Settings):
                    self[name].soft_update(other[name])
            elif name not in self:
                self[name] = other[name]
        return self

    def update(self, other):  # type: ignore
        """Update this instance with data from *other*, overwriting existing keys. Nested |Settings| instances are updated recursively.

        In the following example ``s`` and ``o`` are previously prepared |Settings| instances::

            >>> print(s)
            a:  AA
            b:  BB
            x:
              y1:   XY1
              y2:   XY2
            >>> print(o)
            a:  O_AA
            c:  O_CC
            x:
              y1:   O_XY1
              y3:   O_XY3
            >>> s.update(o)
            >>> print(s)
            a:  O_AA        #original value s.a overwritten by o.a
            b:  BB
            c:  O_CC
            x:
              y1:   O_XY1   #original value s.x.y1 overwritten by o.x.y1
              y2:   XY2
              y3:   O_XY3

        *Other* can also be a regular dictionary. Of course in that case only top level keys are updated.
        """
        for name in other:
            if isinstance(other[name], Settings):
                if name not in self or not isinstance(self[name], Settings):
                    self[name] = other[name].copy()
                else:
                    self[name].update(other[name])
            else:
                self[name] = other[name]

    def merge(self: TSelf, other: "Settings") -> TSelf:
        """Return new instance of |Settings| that is a copy of this instance soft-updated with *other*.

        Shortcut ``A + B`` can be used instead of ``A.merge(B)``.
        """
        ret = self.copy()
        ret.soft_update(other)
        return ret

    def find_case(self, key):
        """Check if this instance contains a key consisting of the same letters as *key*, but possibly with different case. If found, return such a key. If not, return *key*."""
        if not isinstance(key, str):
            return key
        lowkey = key.lower()
        for k in self:
            try:
                if k.lower() == lowkey:
                    return k
            except (AttributeError, TypeError):
                pass
        return key

    def get(self, key, default=None):
        """Like regular ``get``, but ignore the case."""
        return dict.get(self, self.find_case(key), default)

    def pop(self, key, *args):
        """Like regular ``pop``, but ignore the case."""
        # A single positional argument can be supplied `*args`,
        # functioning as a default return value in case `key` is not present in this instance
        return dict.pop(self, self.find_case(key), *args)

    def popitem(self):
        """Like regular ``popitem``, but ignore the case."""
        return dict.popitem(self)

    def setdefault(self, key, default=None):
        """Like regular ``setdefault``, but ignore the case and if the value is a dict, convert it to |Settings|."""
        if isinstance(default, dict) and not isinstance(default, Settings):
            default = Settings(default)
        return dict.setdefault(self, self.find_case(key), default)

    def as_dict(self):
        """Return a copy of this instance with all |Settings| replaced by regular Python dictionaries."""
        d = {}
        for k, v in self.items():
            if isinstance(v, Settings):
                d[k] = v.as_dict()
            elif isinstance(v, list):
                d[k] = [i.as_dict() if isinstance(i, Settings) else i for i in v]
            else:
                d[k] = v

        return d

    @classmethod
    def suppress_missing(cls):
        """A context manager for temporary disabling the :meth:`.Settings.__missing__` magic method: all calls now raising a :exc:`KeyError`.

        As a results, attempting to access keys absent from an arbitrary |Settings| instance will raise a :exc:`KeyError`, thus reverting to the default dictionary behaviour.

        .. note::
            The :meth:`.Settings.__missing__` method is (temporary) suppressed at the class level to ensure consistent invocation by the Python interpreter.
            See also `special method lookup`_.

        Example:

        .. code:: python

            >>> s = Settings()

            >>> with s.suppress_missing():
            ...     s.a.b.c = True
            KeyError: 'a'

            >>> s.a.b.c = True
            >>> print(s.a.b.c)
            True

        .. _`special method lookup`: https://docs.python.org/3/reference/datamodel.html#special-method-lookup

        """
        return SuppressMissing(Settings)

    def get_nested(self, key_tuple, suppress_missing=False):
        """Retrieve a nested value by, recursively, iterating through this instance using the keys in *key_tuple*.

        The :meth:`.Settings.__getitem__` method is called recursively on this instance until all keys in key_tuple are exhausted.

        Setting *suppress_missing* to ``True`` will internally open the :meth:`.Settings.suppress_missing` context manager, thus raising a :exc:`KeyError` if a key in *key_tuple* is absent from this instance.

        .. code:: python

            >>> s = Settings()
            >>> s.a.b.c = True
            >>> value = s.get_nested(('a', 'b', 'c'))
            >>> print(value)
            True
        """
        s = self
        with contextlib.suppress() if not suppress_missing else s.suppress_missing():
            for k in key_tuple:
                s = s[k]
        return s

    def set_nested(self, key_tuple, value, suppress_missing=False):
        """Set a nested value by, recursively, iterating through this instance using the keys in *key_tuple*.

        The :meth:`.Settings.__getitem__` method is called recursively on this instance, followed by :meth:`.Settings.__setitem__`, until all keys in key_tuple are exhausted.


        Setting *suppress_missing* to ``True`` will internally open the :meth:`.Settings.suppress_missing` context manager, thus raising a :exc:`KeyError` if a key in *key_tuple* is absent from this instance.

        .. code:: python

            >>> s = Settings()
            >>> s.set_nested(('a', 'b', 'c'), True)
            >>> print(s)
            a:
              b:
                c: 	True
        """
        s = self
        with contextlib.suppress() if not suppress_missing else s.suppress_missing():
            for k in key_tuple[:-1]:
                s = s[k]
        s[key_tuple[-1]] = value

    def flatten(self, flatten_list=True) -> "Settings":
        """Return a flattened copy of this instance.

        New keys are constructed by concatenating the (nested) keys of this instance into tuples.

        Opposite of the :meth:`.Settings.unflatten` method.

        If *flatten_list* is ``True``, all nested lists will be flattened as well. Dictionary keys are replaced with list indices in such case.

        .. code-block:: python

            >>> s = Settings()
            >>> s.a.b.c = True
            >>> print(s)
            a:
              b:
                c: 	True

            >>> s_flat = s.flatten()
            >>> print(s_flat)
            ('a', 'b', 'c'): 	True
        """
        if flatten_list:
            nested_type: Union[Type, Tuple[Type, ...]] = (Settings, list)
            iter_type = lambda x: x.items() if isinstance(x, Settings) else enumerate(x)
        else:
            nested_type = Settings
            iter_type = Settings.items

        def _concatenate(key_ret, sequence):
            # Switch from Settings.items() to enumerate() if a list is encountered
            for k, v in iter_type(sequence):
                k = key_ret + (k,)
                if isinstance(v, nested_type) and v:  # Empty lists or Settings instances will return ``False``
                    _concatenate(k, v)
                else:
                    ret[k] = v

        # Changes keys into tuples
        ret = Settings()
        _concatenate((), self)
        return ret

    def unflatten(self, unflatten_list=True) -> "Settings":
        """Return a nested copy of this instance.

        New keys are constructed by expanding the keys of this instance (*e.g.* tuples) into new nested |Settings| instances.

        If *unflatten_list* is ``True``, integers will be interpretted as list indices and are used for creating nested lists.

        Opposite of the :meth:`.Settings.flatten` method.

        .. code-block:: python

            >>> s = Settings()
            >>> s[('a', 'b', 'c')] = True
            >>> print(s)
            ('a', 'b', 'c'): 	True

            >>> s_nested = s.unflatten()
            >>> print(s_nested)
            a:
              b:
                c: 	True
        """
        ret = Settings()
        for key, value in self.items():
            s = ret
            for k1, k2 in zip(key[:-1], key[1:]):
                if not unflatten_list:
                    s = s[k1]
                    continue

                if isinstance(k2, int):  # Apply padding to s
                    if not isinstance(s[k1], list):
                        s[k1] = []
                    s[k1] += [Settings()] * (k2 - len(s[k1]) + 1)
                s = s[k1]
            s[key[-1]] = value

        return ret

    # =======================================================================

    def __iter__(self):
        """Iteration through keys follows lexicographical order. All keys are sorted as if they were strings."""
        return iter(sorted(self.keys(), key=str))

    def __missing__(self, name):
        """When requested key is not present, add it with an empty |Settings| instance as a value.

        This method is essential for automatic insertions in deeper levels. Without it things like::

            >>> s = Settings()
            >>> s.a.b.c = 12

        will not work.

        The behaviour of this method can be suppressed by initializing the :class:`.Settings.suppress_missing` context manager.
        """
        self[name] = Settings()
        return self[name]

    def __contains__(self, name):
        """Like regular ``__contains`__``, but ignore the case."""
        return dict.__contains__(self, self.find_case(name))

    def __getitem__(self, name):
        """Like regular ``__getitem__``, but ignore the case."""
        return dict.__getitem__(self, self.find_case(name))

    def __setitem__(self, name, value):
        """Like regular ``__setitem__``, but ignore the case and if the value is a dict, convert it to |Settings|."""
        if isinstance(value, dict) and not isinstance(value, Settings):
            value = Settings(value)
        dict.__setitem__(self, self.find_case(name), value)

    def __delitem__(self, name):
        """Like regular ``__detitem__``, but ignore the case."""
        return dict.__delitem__(self, self.find_case(name))

    def __getattr__(self, name):
        """If name is not a magic method, redirect it to ``__getattribute__``."""
        if name.startswith("__") and name.endswith("__"):
            return dict.__getattribute__(self, name)
        return self[name]

    def __setattr__(self, name, value):
        """If name is not a magic method, redirect it to ``__setattr__``."""
        if name.startswith("__") and name.endswith("__"):
            dict.__setattr__(self, name, value)
        else:
            self[name] = value

    def __delattr__(self, name):
        """If name is not a magic method, redirect it to ``__delattr__``."""
        if name.startswith("__") and name.endswith("__"):
            dict.__delattr__(self, name)
        else:
            del self[name]

    def _str(self, indent):
        """Print contents with *indent* spaces of indentation. Recursively used for printing nested |Settings| instances with proper indentation."""
        ret = ""
        for key, value in self.items():
            ret += " " * indent + str(key) + ": \t"
            if isinstance(value, Settings):
                if len(value) == 0:
                    ret += "<empty Settings>\n"
                else:
                    ret += "\n" + value._str(indent + len(str(key)) + 1)
            else:  # Apply consistent indentation at every '\n' character
                indent_str = " " * (2 + indent + len(str(key))) + "\t"
                ret += textwrap.indent(str(value), indent_str)[len(indent_str) :] + "\n"
        return ret if ret else "<empty Settings>"

    def __str__(self):
        return self._str(0)

    def __dir__(self):
        """
        Return standard attributes, plus dynamically added keys which can be accessed via dot notation.
        """
        return [x for x in super().__dir__()] + [k for k in self.keys() if isinstance(k, str) and k.isidentifier()]

    __repr__ = __str__
    __iadd__ = soft_update
    __add__ = merge
    __copy__ = copy


class SuppressMissing(contextlib.AbstractContextManager):
    """A context manager for temporary disabling the :meth:`.Settings.__missing__` magic method. See :meth:`Settings.suppress_missing` for more details."""

    def __init__(self, obj: type):
        """Initialize the :class:`SuppressMissing` context manager."""
        # Ensure that obj is a class, not a class instance
        self.obj = obj if isinstance(obj, type) else type(obj)
        self.missing = obj.__missing__ if hasattr(obj, "__missing__") else None

    def __enter__(self):
        """Enter the :class:`SuppressMissing` context manager: delete :meth:`.Settings.__missing__` at the class level."""

        @wraps(self.missing)
        def __missing__(self, name):
            raise KeyError(name)

        # The __missing__ method is replaced for as long as the context manager is open
        setattr(self.obj, "__missing__", __missing__)

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the :class:`SuppressMissing` context manager: reenable :meth:`.Settings.__missing__` at the class level."""
        setattr(self.obj, "__missing__", self.missing)


class SafeRunSettings(Settings):
    """
    Safe run settings for global config.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.repeat = 10
        self.delay = 1

    @property
    def repeat(self) -> int:
        """
        Number of attempts for each run() call. Defaults to ``10``.
        """
        return self["repeat"]

    @repeat.setter
    def repeat(self, value: int) -> None:
        self["repeat"] = value

    @property
    def delay(self) -> int:
        """
        Delay between attempts for each run() call. Defaults to ``1``.
        """
        return self["delay"]

    @delay.setter
    def delay(self, value: int) -> None:
        self["delay"] = value


class LogSettings(Settings):
    """
    Log settings for global config.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.file = 5
        self.stdout = 3
        self.time = True
        self.date = True

    @property
    def file(self) -> int:
        """
        Verbosity of the log printed to .log file in the main working folder. Defaults to ``5``.
        """
        return self["file"]

    @file.setter
    def file(self, value: int) -> None:
        self["file"] = value

    @property
    def stdout(self) -> int:
        """
        Verbosity of the log printed to the standard output. Defaults to ``3``.
        """
        return self["stdout"]

    @stdout.setter
    def stdout(self, value: int) -> None:
        self["stdout"] = value

    @property
    def time(self) -> bool:
        """
        When enabled, include write time for each log event. Defaults to ``True``.
        """
        return self["time"]

    @time.setter
    def time(self, value: bool) -> None:
        self["time"] = value

    @property
    def date(self) -> bool:
        """
        When enabled, include write date for each log event. Defaults to ``True``.
        """
        return self["date"]

    @date.setter
    def date(self, value: bool) -> None:
        self["date"] = value


class RunScriptSettings(Settings):
    """
    Run script settings for global config.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.shebang = "#!/bin/sh"
        self.stdout_redirect = False

    @property
    def shebang(self) -> str:
        """
        First line of all produced runscripts. Defaults to ``#!/bin/sh``.
        """
        return self["shebang"]

    @shebang.setter
    def shebang(self, value: str) -> None:
        self["shebang"] = value

    @property
    def stdout_redirect(self) -> bool:
        """
        When enabled, the standard output redirection is handled by the operating system (by using '>[jobname].out' in the runscript), instead of being handled by native Python mechanism.
        Defaults to ``False``.
        """
        return self["stdout_redirect"]

    @stdout_redirect.setter
    def stdout_redirect(self, value: bool) -> None:
        self["stdout_redirect"] = value


class JobSettings(Settings):
    """
    Job settings for global config.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.pickle = True
        self.pickle_protocol = -1
        self.keep = "all"
        self.save = "all"
        self.runscript = RunScriptSettings()
        self.link_files = True

    @property
    def pickle(self) -> bool:
        """
        Enables pickle for the whole job object to [jobname].dill, after job execution is finished. Defaults to ``True``.
        """
        return self["pickle"]

    @pickle.setter
    def pickle(self, value: bool) -> None:
        self["pickle"] = value

    @property
    def pickle_protocol(self) -> int:
        """
        Protocol used for pickling. Defaults to ``-1``.
        """
        return self["pickle_protocol"]

    @pickle_protocol.setter
    def pickle_protocol(self, value: int) -> None:
        self["pickle_protocol"] = value

    @property
    def keep(self) -> str:
        """
        Defines which files should be kept on disk. Defaults to ``all``.
        """
        return self["keep"]

    @keep.setter
    def keep(self, value: str) -> None:
        self["keep"] = value

    @property
    def save(self) -> str:
        """
        Defines which files should be kept on disk. Defaults to ``all``.
        """
        return self["save"]

    @save.setter
    def save(self, value: str) -> None:
        self["save"] = value

    @property
    def runscript(self) -> RunScriptSettings:
        """
        See :class:`~scm.plams.core.settings.RunscriptSettings`.
        """
        return self["runscript"]

    @runscript.setter
    def runscript(self, value: RunScriptSettings) -> None:
        self["runscript"] = value

    @property
    def link_files(self) -> bool:
        """
        When enabled, re-run files will be hardlinked instead of copied, unless on Windows when files are always copied.
        Defaults to ``True``.
        """
        return self["link_files"]

    @link_files.setter
    def link_files(self, value: bool) -> None:
        self["link_files"] = value


class JobManagerSettings(Settings):
    """
    Job manager settings for global config.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.counter_len = 3
        self.hashing = "input"
        self.remove_empty_directories = True

    @property
    def counter_len(self) -> int:
        """
        Number of digits for the counter used when two or more jobs have the same name, when all the jobs apart from the first one are renamed to [jobname].002 ([jobname].003 etc.)
        Defaults to ``3``.
        """
        return self["counter_len"]

    @counter_len.setter
    def counter_len(self, value: int) -> None:
        self["counter_len"] = value

    @property
    def hashing(self) -> str:
        """
        Hashing method used for testing if some job was previously run. Defaults to ``input``.
        """
        return self["hashing"]

    @hashing.setter
    def hashing(self, value: str) -> None:
        self["hashing"] = value

    @property
    def remove_empty_directories(self) -> bool:
        """
        When enabled, removes of all empty subdirectories in the main working folder at the end of the script.
        Defaults to ``True``.
        """
        return self["remove_empty_directories"]

    @remove_empty_directories.setter
    def remove_empty_directories(self, value: bool) -> None:
        self["remove_empty_directories"] = value


class ConfigSettings(Settings):
    """
    Extends the default |Settings| with standard options which are required for global config.
    The values for these options are initialised to default values.
    The default |JobRunner| and |JobManager| are lazily initialised when first accessed.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.init = False
        self._explicit_init = False
        self.preview = False
        self.sleepstep = 5
        self.ignore_failure = True
        self.daemon_threads = True
        self.erase_workdir = False
        self.jobmanager = JobManagerSettings()
        self.job = JobSettings()
        self.log = LogSettings()
        self.saferun = SafeRunSettings()

        # Default job runner and job manager are lazily initialised on first access
        # This is to allow users to change their settings before initialisation (due to side effects in init)
        self.default_jobrunner = None
        self.default_jobmanager = None

    @property
    def init(self) -> bool:
        """
        Whether config has been marked as fully initialized and jobs are ready to be run. Defaults to ``False``.
        """
        return self["init"]

    @init.setter
    def init(self, value: bool):
        self["init"] = value

    @property
    def _explicit_init(self) -> bool:
        """
        Whether config has been explicitly initialized by the user by calling |init|. Defaults to ``False``.
        """
        return self["_explicit_init"]

    @_explicit_init.setter
    def _explicit_init(self, value: bool):
        self["_explicit_init"] = value

    @property
    def preview(self) -> bool:
        """
        When enabled, no actual calculations are run, only inputs and runscripts are prepared. Defaults to ``False``.
        """
        return self["preview"]

    @preview.setter
    def preview(self, value: bool):
        self["preview"] = value

    @property
    def sleepstep(self) -> int:
        """
        Unit of time which is used whenever some action needs to be repeated until a certain condition is met. Defaults to ``5``.
        """
        return self["sleepstep"]

    @sleepstep.setter
    def sleepstep(self, value: int) -> None:
        self["sleepstep"] = value

    @property
    def ignore_failure(self) -> bool:
        """
        When enabled, accessing a failed/crashed job gives a log message instead of an error. Defaults to ``True``.
        """
        return self["ignore_failure"]

    @ignore_failure.setter
    def ignore_failure(self, value: bool) -> None:
        self["ignore_failure"] = value

    @property
    def daemon_threads(self) -> bool:
        """
        When enabled, all threads started by JobRunner are daemon threads, which are terminated when the main thread finishes,
        and hence allow immediate end of the parallel script when Ctrl-C is pressed. Defaults to ``True``.
        """
        return self["daemon_threads"]

    @daemon_threads.setter
    def daemon_threads(self, value: bool) -> None:
        self["daemon_threads"] = value

    @property
    def erase_workdir(self) -> bool:
        """
        When enabled, the entire main working folder is deleted at the end of script. Defaults to ``False``.
        :return:
        """
        return self["erase_workdir"]

    @erase_workdir.setter
    def erase_workdir(self, value: bool) -> None:
        self["erase_workdir"] = value

    @property
    def jobmanager(self) -> JobManagerSettings:
        """
        See :class:`~scm.plams.core.settings.JobManagerSettings`.
        """
        return self["jobmanager"]

    @jobmanager.setter
    def jobmanager(self, value: JobManagerSettings) -> None:
        self["jobmanager"] = value

    @property
    def job(self) -> JobSettings:
        """
        See :class:`~scm.plams.core.settings.JobSettings`.
        """
        return self["job"]

    @job.setter
    def job(self, value: JobSettings) -> None:
        self["job"] = value

    @property
    def log(self) -> LogSettings:
        """
        See :class:`~scm.plams.core.settings.LogSettings`.
        """
        return self["log"]

    @log.setter
    def log(self, value: LogSettings) -> None:
        self["log"] = value

    @property
    def saferun(self) -> SafeRunSettings:
        """
        See :class:`~scm.plams.core.settings.SafeRunSettings`.
        """
        return self["saferun"]

    @saferun.setter
    def saferun(self, value: SafeRunSettings) -> None:
        self["saferun"] = value

    @property
    def default_jobrunner(self) -> "JobRunner":
        """
        Default |JobRunner| that will be used for running jobs, if an explicit runner is not provided.
        """
        from scm.plams.core.jobrunner import JobRunner

        if self["default_jobrunner"] is None:
            self["default_jobrunner"] = JobRunner()
        return self["default_jobrunner"]

    @default_jobrunner.setter
    def default_jobrunner(self, value: "JobRunner") -> None:
        self["default_jobrunner"] = value

    @property
    def default_jobmanager(self) -> "JobManager":
        """
        Default |JobManager| that will be used when running jobs, if an explicit manager is not provided.
        """
        from scm.plams.core.jobmanager import JobManager

        if self["default_jobmanager"] is None:
            self["default_jobmanager"] = JobManager(self.jobmanager)
        return self["default_jobmanager"]

    @default_jobmanager.setter
    def default_jobmanager(self, value: "JobManager") -> None:
        self["default_jobmanager"] = value

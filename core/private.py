import copy
import hashlib
import os
import subprocess
import sys
import time
import warnings
from contextlib import AbstractContextManager
from os.path import join as opj
from typing import Callable, Dict, NoReturn, List, Optional

__all__: List[str] = []


def smart_copy(obj, owncopy=[], without=[]):
    """Return a copy of *obj*. Attributes of *obj* listed in *without* are ignored. Attributes listed in *owncopy* are copied by calling their own ``copy()`` methods. All other attributes are copied using :func:`copy.deepcopy`."""

    ret = obj.__class__()
    for k in owncopy:
        ret.__dict__[k] = obj.__dict__[k].copy()
    for k in obj.__dict__:
        if k not in (without + owncopy):
            ret.__dict__[k] = copy.deepcopy(obj.__dict__[k])
    return ret


# ===========================================================================


def sha256(string):
    """A small utility wrapper around :ref:`hashlib.sha256<hash-algorithms>`."""
    if not isinstance(string, bytes):
        string = str(string).encode()
    h = hashlib.sha256()
    h.update(string)
    return h.hexdigest()


# ===========================================================================


def saferun(*args, **kwargs):
    """A wrapper around :func:`subprocess.run` repeating the call ``config.saferun.repeat`` times with ``config.saferun.delay`` interval in case of :exc:`BlockingIOError` being raised (any other exception is not caught and directly passed above). All arguments (*args* and *kwargs*) are passed directly to :func:`~subprocess.run`. If all attempts fail, the last raised :exc:`BlockingIOError` is reraised."""
    from scm.plams.core.functions import config, log

    attempt = 0
    (repeat, delay) = (config.saferun.repeat, config.saferun.delay) if ("saferun" in config) else (5, 1)
    while attempt <= repeat:
        try:
            return subprocess.run(*args, **kwargs)
        except OSError as e:
            attempt += 1
            log("subprocess.run({}) attempt {} failed with {}".format(args[0], attempt, e), 5)
            last_error = e
            time.sleep(delay)
    raise last_error


# ===========================================================================


class UpdateSysPath(AbstractContextManager):
    """A context manager for temporary adding ``$AMSHOME/scripting`` to :data:`sys.path`.

    While the context manager is opened modules can be imported (directly) from the ``$AMSHOME/scripting`` directory,
    allowing one access to the python modules distributed with ADF.
    An :exc:`EnvironmentError` is raised if the ``AMSHOME`` environment variable has not been set.

    Usage example:

    .. code:: python
        >>> import sys
        >>> import os

        >>> scripting = os.path.join(os.environ['AMSHOME'], 'scripting')
        >>> with UpdateSysPath():
        ...     print(scripting in sys.path)  # The context manager is opened
        True

        >>> print(scripting in sys.path)  # The context manager is closed
        False

    """

    def __init__(self, path: Optional[str] = None) -> None:
        try:
            parser_path = path if path is not None else opj(os.environ["AMSHOME"], "scripting")
        except KeyError as ex:
            raise EnvironmentError("The 'AMSHOME' environment variable has not been set").with_traceback(
                ex.__traceback__
            )

        if parser_path not in sys.path:
            sys.path.append(parser_path)
            self.path: Optional[str] = parser_path
        else:  # The specified path is already present in sys.path
            self.path = None

    def __enter__(self) -> None:
        """If not ``None``, add :attr:`UpdateSyspath.path` to :data:`sys.path`."""
        if self.path is not None:
            sys.path.append(self.path)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """If not ``None``, remove :attr:`UpdateSyspath.path` from :data:`sys.path`."""
        if self.path is not None:
            try:
                idx = sys.path.index(self.path)
            except ValueError:  # self.path has been (manually) removed by the user
                pass
            else:
                del sys.path[idx]


# ===========================================================================


def _raise(exc: Exception) -> NoReturn:
    """Raise *exc*; this tiny function is here because ``lambda`` expressions can't handle the ``raise`` statement."""
    raise exc


#: Map a :class:`str` representation of an action to an actual callable.
_ACTION_DICT: Dict[str, Callable[[Exception], None]] = {
    "ignore": lambda n: None,
    "raise": lambda n: _raise(n),
    "warn": lambda n: warnings.warn(str(n), stacklevel=2),
}


def parse_action(action: str) -> Callable[[Exception], None]:
    """Return a callable which takes an :exc:`Exception` as argument and performs some action base on the value of *action*.

    * ``action = 'ignore'``: Do nothing.
    * ``action = 'warn'``: Issue a warning.
    * ``action = 'raise'``: Raise the passed exception.

    """
    try:
        return _ACTION_DICT[action]
    except KeyError as ex:
        raise ValueError("`action` expected 'ignore', 'raise' or 'warn'; " f"observed value: {action!r}") from ex


# ===========================================================================


def retry(sleep=0.1, maxtries=10):
    # wrapper for sleep-retrying a function call. use with `@retry()`
    from time import sleep as _sleep

    def wrap1(f):
        def wrap2(*a, __count=0, **kw):
            try:
                return f(*a, **kw)
            except Exception as e:
                if __count > maxtries:
                    raise e from None  # ignores stack trace
                _sleep(sleep)
                wrap2(
                    *a,
                    __count=__count + 1,
                    **kw,
                )

        return wrap2

    return wrap1

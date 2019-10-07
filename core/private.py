import os
import sys
import copy
import hashlib
import subprocess
import time

from os.path import join as opj
from contextlib import AbstractContextManager


from .functions import config, log


__all__ = []



def smart_copy(obj, owncopy=[], without=[]):
    """Return a copy of *obj*. Attributes of *obj* listed in *without* are ignored. Attributes listed in *owncopy* are copied by calling their own ``copy()`` methods. All other attributes are copied using :func:`copy.deepcopy`."""

    ret = obj.__class__()
    for k in owncopy:
        ret.__dict__[k] = obj.__dict__[k].copy()
    for k in obj.__dict__:
        if k not in (without + owncopy):
            ret.__dict__[k] = copy.deepcopy(obj.__dict__[k])
    return ret


#===========================================================================


def sha256(string):
    """A small utility wrapper around :ref:`hashlib.sha256<hash-algorithms>`."""
    if not isinstance(string, bytes):
        string = str(string).encode()
    h = hashlib.sha256()
    h.update(string)
    return h.hexdigest()


#===========================================================================


def saferun(*args, **kwargs):
    """A wrapper around :func:`subprocess.run` repeating the call ``config.saferun.repeat`` times with ``config.saferun.delay`` interval in case of :exc:`BlockingIOError` being raised (any other exception is not caught and directly passed above). All arguments (*args* and *kwargs*) are passed directly to :func:`~subprocess.run`. If all attempts fail, the last raised :exc:`BlockingIOError` is reraised."""
    attempt = 0
    (repeat, delay) = (config.saferun.repeat, config.saferun.delay) if ('saferun' in config) else (5,1)
    while attempt <= repeat:
        try:
            return subprocess.run(*args, **kwargs)
        except BlockingIOError as e:
            attempt += 1
            log('subprocess.run({}) attempt {} failed with {}'.format(args[0], attempt, e), 5)
            last_error = e
            time.sleep(delay)
    raise last_error


#===========================================================================


class UpdateSysPath(AbstractContextManager):
    """A context manager for temporary adding ``$ADFHOME/scripting`` to :data:`sys.path`.

    While the context manager is opened modules can be imported (directly) from the ``$ADFHOME/scripting`` directory,
    allowing one access to the python modules distributed with ADF.
    An :exc:`EnvironmentError` is raised if the ``ADFHOME`` environment variable has not been set.

    Usage example:

    .. code:: python
        >>> import sys
        >>> import os

        >>> scripting = os.path.join(os.environ['ADFHOME'], 'scripting')
        >>> with UpdateSysPath():
        ...     print(scripting in sys.path)  # The context manager is opened
        True

        >>> print(scripting in sys.path)  # The context manager is closed
        False

    """

    def __init__(self, path: str = None) -> None:
        try:
            parser_path = path if path is not None else opj(os.environ['ADFHOME'], 'scripting')
        except KeyError as ex:
            raise EnvironmentError("The 'ADFHOME' environment variable has not been set").with_traceback(ex.__traceback__)

        if parser_path not in sys.path:
            sys.path.append(parser_path)
            self.path = parser_path
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

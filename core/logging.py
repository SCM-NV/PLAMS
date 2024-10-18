import logging
import sys
from typing import Optional, Dict
import threading


from scm.plams.core.errors import FileError

__all__ = ["get_logger", "Logger"]


def get_logger(name: str) -> "Logger":
    """
    Get a logger with the specified name.
    If there is an existing logger with the same name this is returned, otherwise a new logger is created.
    """
    return LogManager.get_logger(name)


class LogManager:
    """
    Manages PLAMS logger instances.
    The manager should not be instantiated directly, but loggers accessed through the ``get_logger`` method.
    """

    _loggers: Dict[str, "Logger"] = {}

    def __new__(cls, *args, **kwargs):
        raise TypeError("LoggerManager cannot be directly instantiated.")

    @classmethod
    def get_logger(cls, name: str) -> "Logger":
        """
        Get a logger with the specified name.
        If there is an existing logger with the same name this is returned, otherwise a new logger is created.
        """
        if name not in cls._loggers:
            cls._loggers[name] = Logger(name)
        return cls._loggers[name]


class Logger:
    """
    Wrapper around default logger, which handles logging to stdout and a text logfile depending on the options in the
    global logging config.
    """

    def __init__(self, name: str):
        """
        Get a logger with given name.
        """
        self._logger = logging.Logger(name)
        self._stdout_handler: Optional[logging.StreamHandler] = None
        self._file_handler: Optional[logging.FileHandler] = None
        self._formatter: Optional[logging.Formatter] = None
        self._lock = threading.Lock()

    def configure(
        self,
        stdout_level: int = 0,
        logfile_level: int = 0,
        logfile_path: Optional[str] = None,
        include_date: bool = False,
        include_time: bool = False,
    ) -> None:
        """
        Configure logging to stdout and the logfile, and its formatting.

        For backwards compatibility, the logging level is between 0-7, with 0 indicating no logging and 7 the most verbose logging.
        Note that this is only a PLAMS logging level, it will be mapped to a level between INFO and WARNING for the
        standard python logger.

        The logfile path must be unique across all loggers, otherwise a |FileError| is raised. If set to None this will disable file logging.

        :param stdout_level: value between 0-7, with 0 indicating no logging and 7 the most verbose logging to stdout
        :param logfile_level: value between 0-7, with 0 indicating no logging and 7 the most verbose logging to logfile
        :param logfile_path: path for the logfile, if set to None this will remove file logging
        :param include_date: whether to include date stamp at the start of a log line
        :param include_time: whether to include time stamp at the start of a log line
        """
        with self._lock:

            # Initialise the stdout handler once
            if self._stdout_handler is None:
                self._stdout_handler = logging.StreamHandler(sys.stdout)
                self._stdout_handler.setFormatter(self._formatter)
                self._logger.addHandler(self._stdout_handler)

            # Update the stdout handler level if required
            if stdout_level != self._stdout_handler.level:
                self._stdout_handler.setLevel(28 - stdout_level)

            # Remove and close existing file handler if present and required
            if self._file_handler is not None and (
                logfile_path is None or logfile_path != self._file_handler.baseFilename
            ):
                self._logger.removeHandler(self._file_handler)
                self._file_handler.flush()
                self._file_handler.close()

            # Add new file handler if required
            if logfile_path is not None and (
                self._file_handler is None or logfile_path != self._file_handler.baseFilename
            ):
                # Check logfile is not already in use by another logger
                # This will cause permission errors on Windows and generally is not a good idea
                for name, logger in LogManager._loggers.items():
                    if logger._file_handler is not None and logger._file_handler.baseFilename == logfile_path:
                        raise FileError(f"Logger '{name}' already exists with logfile path '{logfile_path}'")

                self._file_handler = logging.FileHandler(logfile_path)
                self._file_handler.setFormatter(self._formatter)
                self._logger.addHandler(self._file_handler)

            # Update the logfile handler level if required
            if self._file_handler is not None:
                self._file_handler.setLevel(28 - logfile_level)

            # Configure the formatter if required
            datefmt = None
            if include_date and include_time:
                datefmt = "[%d.%m|%H:%M:%S]"
            elif include_date:
                datefmt = "[%d.%m]"
            elif include_time:
                datefmt = "[%H:%M:%S]"

            if self._formatter is None or datefmt != self._formatter.datefmt:
                fmt = "%(asctime)s %(message)s" if datefmt is not None else None
                self._formatter = logging.Formatter(fmt, datefmt=datefmt)
                if self._stdout_handler is not None:
                    self._stdout_handler.setFormatter(self._formatter)
                if self._file_handler is not None:
                    self._file_handler.setFormatter(self._formatter)

    def log(self, message: str, level: int) -> None:
        """
        Log a message with the given level of verbosity.
        """
        # Shouldn't really have to take a lock here as logging itself is thread-safe
        # but in the PLAMS log function the logfile can be configured on a per-call basis
        # which could easily lead to dropped logs if this were multi-threaded.
        with self._lock:
            self._logger.log(28 - level, message)

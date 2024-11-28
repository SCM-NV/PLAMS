import logging
import sys
from typing import Optional, Dict, Literal, Any
import threading
from abc import ABC, abstractmethod


from scm.plams.core.errors import FileError, PlamsError

__all__ = ["get_logger", "TextLogger"]


def get_logger(name: str, fmt: Optional[Literal["txt"]] = None) -> "Logger":
    """
    Get a logger with the specified name.
    If there is an existing logger with the same name this is returned, otherwise a new logger is created.
    If no format is specified, the logger will be a simple ``TextLogger``.
    """
    return LogManager.get_logger(name, fmt)


class LogManager:
    """
    Manages PLAMS logger instances.
    The manager should not be instantiated directly, but loggers accessed through the ``get_logger`` method.
    """

    _loggers: Dict[str, "Logger"] = {}

    def __new__(cls, *args, **kwargs):
        raise TypeError("LoggerManager cannot be directly instantiated.")

    @classmethod
    def get_logger(cls, name: str, fmt: Optional[Literal["txt"]] = None) -> "Logger":
        """
        Get a logger with the specified name.
        If there is an existing logger with the same name this is returned, otherwise a new logger is created.
        """
        if name not in cls._loggers:
            if fmt == "txt" or fmt is None:
                logger = TextLogger(name)
            else:
                raise PlamsError(f"Cannot create logger '{name}' with format '{fmt}'")
            cls._loggers[name] = logger
        return cls._loggers[name]


class Logger(ABC):
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
        self._stdout_formatter: Optional[logging.Formatter] = None
        self._file_handler: Optional[logging.FileHandler] = None
        self._file_formatter: Optional[logging.Formatter] = None
        self._lock = threading.Lock()

    @abstractmethod
    def configure(
        self,
        stdout_level: int = 0,
        logfile_level: int = 0,
        logfile_path: Optional[str] = None,
        *args,
        **kwargs,
    ) -> None:
        """
        Configure logging to stdout and the logfile, and its formatting.

        For backwards compatibility, the logging level is between 0-7, with 0 indicating no logging and 7 the most verbose logging.
        Note that this is only a PLAMS logging level, it will be mapped to a level between INFO and WARNING for the
        standard python logger.

        The logfile path must be unique across all loggers, otherwise a |FileError| is raised. If set to None this will disable file logging.
        """
        pass

    def _configure_stdout_handler(self, level: int):
        """
        Configure the stdout handler, initializing and adjusting the level if required
        """

        # Initialise the stdout handler if not already done so
        if self._stdout_handler is None:
            self._stdout_handler = logging.StreamHandler(sys.stdout)
            self._stdout_handler.setFormatter(self._stdout_formatter)
            self._logger.addHandler(self._stdout_handler)

        # Update the logfile handler level if required
        if level != 28 - self._stdout_handler.level:
            self._stdout_handler.setLevel(28 - level)

    def _configure_file_handler(self, level: int, logfile_path: Optional[str]):
        """
        Configure the file handler, setting the logfile and adjusting the level if required
        """

        # Remove and close existing file handler if present and required
        if self._file_handler is not None and (logfile_path is None or logfile_path != self._file_handler.baseFilename):
            self._remove_handler(self._file_handler)
            self._file_handler = None

        # Add new file handler if required
        if logfile_path is not None and (self._file_handler is None or logfile_path != self._file_handler.baseFilename):
            # Check logfile is not already in use by another logger
            # This will cause permission errors on Windows and generally is not a good idea
            for name, logger in LogManager._loggers.items():
                if logger._file_handler is not None and logger._file_handler.baseFilename == logfile_path:
                    raise FileError(f"Logger '{name}' already exists with logfile path '{logfile_path}'")

            self._file_handler = logging.FileHandler(logfile_path)
            self._file_handler.setFormatter(self._file_formatter)
            self._logger.addHandler(self._file_handler)

        # Update the logfile handler level if required
        if self._file_handler is not None:
            self._file_handler.setLevel(28 - level)

    def _configure_stdout_formatter(self, formatter: logging.Formatter) -> None:
        """
        Configure the formatter for the stdout handler
        """
        self._stdout_formatter = formatter
        if self._stdout_handler is not None:
            self._stdout_handler.setFormatter(self._stdout_formatter)

    def _configure_file_formatter(self, formatter: logging.Formatter) -> None:
        """
        Configure the formatter for the file handler
        """
        self._file_formatter = formatter
        if self._file_handler is not None:
            self._file_handler.setFormatter(self._file_formatter)

    def _remove_handler(self, handler: Optional[logging.Handler]) -> None:
        """
        Flush logs, close the handler and remove it from the logger.
        """
        if handler is not None:
            self._logger.removeHandler(handler)
            try:
                handler.flush()
            except ValueError:
                pass  # Already closed
            handler.close()

    def close(self) -> None:
        """
        Flush logs to stdout and logfile, then close the resources.
        No further logs will then be written.
        """
        with self._lock:
            self._remove_handler(self._file_handler)
            self._remove_handler(self._stdout_handler)
            self._file_handler = None
            self._stdout_handler = None

    @abstractmethod
    def log(self, message: Any, level: int) -> None:
        """
        Log a message with the given level of verbosity.
        """
        pass


class TextLogger(Logger):
    """
    Wrapper around default logger, which handles simple text logging to stdout and a text logfile depending on the options in the
    global logging config.
    """

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

            self._configure_stdout_handler(stdout_level)
            self._configure_file_handler(logfile_level, logfile_path)

            # Configure formatter if required
            datefmt = None
            if include_date and include_time:
                datefmt = "[%d.%m|%H:%M:%S]"
            elif include_date:
                datefmt = "[%d.%m]"
            elif include_time:
                datefmt = "[%H:%M:%S]"
            fmt = "%(asctime)s %(message)s" if datefmt is not None else None

            if self._stdout_formatter is None or datefmt != self._stdout_formatter.datefmt:
                self._configure_stdout_formatter(logging.Formatter(fmt, datefmt=datefmt))

            if self._file_formatter is None or datefmt != self._file_formatter.datefmt:
                self._configure_file_formatter(logging.Formatter(fmt, datefmt=datefmt))

    def log(self, message: str, level: int) -> None:
        """
        Log a message with the given level of verbosity.
        """
        # Shouldn't really have to take a lock here as logging itself is thread-safe
        # but in the PLAMS log function the logfile can be configured on a per-call basis
        # which could easily lead to dropped logs if this were multi-threaded.
        with self._lock:
            self._logger.log(28 - level, message)

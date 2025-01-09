import csv
import logging
import os
import sys
from io import StringIO
from typing import Any, Type, Dict, Literal, Optional, overload
import threading
from abc import ABC, abstractmethod

from scm.plams.core.errors import FileError, PlamsError

__all__ = ["get_logger", "TextLogger", "CSVLogger", "CSVFormatter"]


@overload
def get_logger(name: str, fmt: Optional[Literal["txt"]] = None) -> "TextLogger": ...
@overload
def get_logger(name: str, fmt: Literal["csv"]) -> "CSVLogger": ...
def get_logger(name: str, fmt: Optional[Literal["txt", "csv"]] = None) -> "Logger":
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
    def get_logger(cls, name: str, fmt: Optional[Literal["txt", "csv"]] = None) -> "Logger":
        """
        Get a logger with the specified name.
        If there is an existing logger with the same name this is returned, otherwise a new logger is created.
        """
        if name not in cls._loggers:
            if fmt == "csv":
                logger: Logger = CSVLogger(name)
            elif fmt == "txt" or fmt is None:
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
        logfile_path = (
            os.path.abspath(logfile_path) if logfile_path is not None else None
        )  # convert any relative paths to abs paths
        if self._file_handler is not None and (logfile_path is None or logfile_path != self._file_handler.baseFilename):
            self._remove_handler(self._file_handler)
            self._file_handler = None
            self._file_formatter = None

        # Add new file handler if required
        if logfile_path is not None and (self._file_handler is None or logfile_path != self._file_handler.baseFilename):
            # Check logfile is not already in use by another logger
            # This will cause permission errors on Windows and generally is not a good idea
            for name, logger in LogManager._loggers.items():
                if logger._file_handler is not None and logger._file_handler.baseFilename == logfile_path:
                    raise FileError(f"Logger '{name}' already exists with logfile path '{logfile_path}'")

            try:
                self._file_handler = logging.FileHandler(logfile_path)
                self._file_handler.setFormatter(self._file_formatter)
                self._logger.addHandler(self._file_handler)
            except FileNotFoundError:
                # Logger should not error if logfile directory does not exist
                pass

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

    @property
    def logfile(self) -> Optional[str]:
        """
        Path of the logfile currently used for the logger.
        """
        return self._file_handler.baseFilename if self._file_handler is not None else None

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
            self._file_formatter = None
            self._stdout_formatter = None

    def log(self, message: Any, level: int) -> None:
        """
        Log a message with the given level of verbosity.

        :param message: The message to log.
        :param level: Verbosity level (1=important, 3=normal, 5=verbose, 7=debug).
        """
        # Shouldn't really have to take a lock here as logging itself is thread-safe
        # but in the PLAMS log function the logfile can be configured on a per-call basis
        # which could easily lead to dropped logs if this were multi-threaded.
        with self._lock:
            self._logger.log(28 - level, message)


class TextLogger(Logger):
    """
    Logger which performs simple text logging to stdout and a logfile.
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


class CSVFormatter(logging.Formatter):
    """
    Formatter which creates comma-separated log lines from a log record.
    """

    def __init__(
        self,
        datefmt: Optional[str] = None,
        include_level: bool = False,
        include_name: bool = False,
    ) -> None:
        """
        Initialize a new formatter with given logging options.

        :param datefmt: format of datetime for ``asctime`` field, by default none is included
        :param include_level: whether to include the logging level as a ``level`` field
        :param include_name: whether to include the logger name as the ``logger_name`` field
        """

        super().__init__(datefmt=datefmt)

        self.log_time = datefmt is not None
        self.include_level = include_level
        self.include_name = include_name
        self.headers: Optional[list[str]] = None
        self._write_headers = True

    @property
    def write_headers(self) -> bool:
        """
        Whether the formatter needs to write the csv headers for the next log line.
        """
        return self._write_headers

    @write_headers.setter
    def write_headers(self, value: bool):
        self._write_headers = value

    def format(self, record: logging.LogRecord) -> str:
        """
        Format the log record into a CSV row.

        :param record: The log record to format.
        :return: A CSV-formatted string.
        """
        # Extract core log fields
        log_record: Dict[str, Any] = {}
        if self.include_name:
            log_record["logger_name"] = record.name
        if self.include_level:
            log_record["level"] = 28 - record.levelno
        if self.log_time:
            log_record["logged_at"] = self.formatTime(record, self.datefmt)

        if isinstance(record.msg, dict):
            log_record.update(record.msg)
        else:
            log_record["message"] = record.getMessage()

        # Create headers if they are not set
        if self.headers is None:
            self.headers = list(log_record.keys())

        # Write headers if this is the first entry
        row = StringIO()
        csv_writer = csv.DictWriter(row, fieldnames=self.headers, lineterminator="\n")
        if self.write_headers:
            csv_writer.writeheader()
            self.write_headers = False

        csv_writer.writerow(log_record)
        return row.getvalue().strip()

    def __eq__(self, other):
        if type(other) is not type(self):
            return False
        return (
            self.datefmt == other.datefmt
            and self.include_level == other.include_level
            and self.include_name == other.include_name
        )


class CSVLogger(Logger):
    """
    Logger that logs dictionary messages in CSV format with optional date and time stamps.
    """

    _file_formatter: CSVFormatter

    def configure(
        self,
        stdout_level: int = 0,
        logfile_level: int = 0,
        logfile_path: Optional[str] = None,
        include_date: bool = False,
        include_time: bool = False,
        include_level: bool = False,
        include_name: bool = False,
        csv_formatter: Type[CSVFormatter] = CSVFormatter,
    ):
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
        :param include_level: whether to include the logging level in the log lines
        :param include_name: whether to include the logger name in the log lines
        :param csv_formatter: type of CSV formatter to be used for the logger
        """
        with self._lock:

            self._configure_stdout_handler(stdout_level)
            self._configure_file_handler(logfile_level, logfile_path)

            # Configure formatter if required
            datefmt = None
            if include_date and include_time:
                # Format chosen because it supports datetime and pandas.to_datetime conversion
                datefmt = "%Y-%m-%d %H:%M:%S"
            elif include_date:
                datefmt = "%Y-%m-%d"
            elif include_time:
                datefmt = "%H:%M:%S"

            stdout_formatter = csv_formatter(datefmt, include_level, include_name)
            if self._stdout_formatter is None or stdout_formatter != self._stdout_formatter:
                self._configure_stdout_formatter(stdout_formatter)

            file_formatter = csv_formatter(datefmt, include_level, include_name)
            if self._file_formatter is None or file_formatter != self._file_formatter:
                self._configure_file_formatter(file_formatter)

            # For files, check if the logfile is empty and therefore if the headers need to be written
            if self._file_handler is not None and self._file_formatter is not None:
                logfile = self._file_handler.baseFilename
                has_logs = os.path.exists(logfile) and os.path.getsize(logfile) > 0
                self._file_formatter.write_headers = not has_logs

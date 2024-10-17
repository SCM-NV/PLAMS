import csv
import re
import threading
import time
from pathlib import Path
from typing import TYPE_CHECKING, Union, Mapping

_stdlock = threading.Lock()
_filelock = threading.Lock()


if TYPE_CHECKING:
    from scm.plams.core.basejob import Job


class LoggerCSV:
    """basic logger that mimic the behavior of plams.log with adding option of personalize the file logging location"""

    def __init__(self, log_path="logfile.csv", stdout_level=3, date=True, time=True, init: bool = True):
        self.log_path = Path(log_path)
        self.stdout_level = stdout_level
        self.date = date
        self.time = time
        self.init = init

    def log(self, message: Mapping[str, Union[str, int, float]], level=7):
        """Log a dictionary message in CSV format with verbosity level.

        Logs are printed independently to the CSV logfile and to the standard output based on configuration levels.
        Date and/or time can be added based on configuration. All logging activity is/should thread safe.

        :param message: The message to log as a dictionary.
        :type message: dict
        :param level: 1: important; 3: normal; 5: verbose; 7: debug, defaults to 0
        :type level: int, optional
        """
        if not self.init:
            if level <= 3:
                with _stdlock:
                    print(str(message))
            return

        current_headers = get_current_headers(self.log_path)
        new_headers_needed = False

        # Prepare row and check headers
        row = []
        headers = current_headers.copy()
        for key in message:
            if key not in headers:
                headers.append(key)
                new_headers_needed = True
            row.append(message.get(key))
        extra_headers_needed = new_headers_needed and (current_headers != [])

        if extra_headers_needed:
            raise NotImplementedError(
                f"found new headers! {current_headers=} BUT {headers=}, this option is not currently supported"
            )

        # Add date and time if configured
        if self.date:
            row.insert(0, time.strftime("%Y-%m-%d"))
            if "Date" not in headers:
                headers.insert(0, "Date")
        if self.time:
            row.insert(1 if self.date else 0, time.strftime("%H:%M:%S"))
            if "Time" not in headers:
                headers.insert(1 if self.date else 0, "Time")

        # Log to file
        with _filelock:
            with open(self.log_path, mode="a", newline="") as file:
                csv_writer = csv.writer(file)
                if not current_headers:
                    csv_writer.writerow(headers)
                csv_writer.writerow(row)

        # Log to stdout if needed
        if level <= self.stdout_level:
            with _stdlock:
                print(", ".join([str(x) for x in row]))

    def job_log(self, job: "Job", level=7):
        message = get_message_to_log(job)
        self.log(message=message, level=level)


def get_message_to_log(job: "Job"):

    pattern = r"(\.\d{%i})+$" % 3  # (jm.settings.counter_len)
    base_name = re.sub(pattern, "", job.name)
    message = {
        "job_base_name": base_name,
        "job_path": job.path,
        "job_name": job.name,
        "job_ok": job.ok(),
        "job_check": "",
        "job_get_errormsg": "",
        "job_parent_name": "",
        "job_parent_path": "",
    }
    try:
        message.update({"job_check": job.check()})
    except TypeError:
        pass
    try:
        # this one it is not supported by the Job class but on many jobs they have it implemented!
        message.update({"job_get_errormsg": job.get_errormsg()})
    except (AttributeError, TypeError):
        pass
    if job.parent:
        add = {
            "job_parent_name": job.parent.name,
            "job_parent_path": job.parent.path,
        }
        message.update(add)
    return message


def get_current_headers(filepath: Path):
    """Retrieve current headers from the CSV file if it exists."""
    if filepath.exists():
        with open(filepath, mode="r", newline="") as file:
            reader = csv.reader(file)
            return next(reader, [])
    return []

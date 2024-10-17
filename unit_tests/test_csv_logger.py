import os
from scm.plams.core.csv_logger import LoggerCSV
import pytest


def test_csv_logger():
    csv_logger = LoggerCSV(log_path=os.getcwd() + "/logfile.csv")
    print(csv_logger.log_path)

    message = {"aaa": 5, "bbb": 33}
    csv_logger.log(message=message)
    message = {"aaa": 5, "bbb": 777}
    csv_logger.log(message=message)
    with pytest.raises(NotImplementedError):
        message = {"aaa": 5, "ff": 777}
        csv_logger.log(message=message)

    with open(csv_logger.log_path, "r") as f:
        lines = f.readlines()
    assert lines[0] == "Date,Time,aaa,bbb\n"
    assert len(lines) == 3
    os.remove(csv_logger.log_path)

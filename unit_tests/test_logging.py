import re
import uuid
from unittest.mock import patch
from io import StringIO
import threading
import pytest
import time
import random

from scm.plams.core.errors import FileError
from scm.plams.core.logging import get_logger
from scm.plams.unit_tests.test_helpers import temp_file_path


class TestGetLogger:

    def test_get_logger_returns_existing_or_creates_new(self):
        name1 = str(uuid.uuid4())
        name2 = str(uuid.uuid4())
        logger1 = get_logger(name1)
        logger2 = get_logger(name2)
        logger3 = get_logger(name1)

        assert logger1 == logger3 != logger2


class TestLogger:

    def test_no_logging_to_stdout_by_default(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            logger = get_logger(str(uuid.uuid4()))
            logger.log("hello", 1)

            assert mock_stdout.getvalue() == ""

    def test_configure_writes_to_stdout_up_to_and_including_level(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            logger = get_logger(str(uuid.uuid4()))
            logger.configure(3)
            for i in range(10):
                logger.log(f"log line {i}", i)

            assert (
                mock_stdout.getvalue()
                == """log line 0
log line 1
log line 2
log line 3
"""
            )

    def test_configure_writes_to_logfile_up_to_and_including_level(self):
        with temp_file_path(suffix=".log") as temp_log_file:
            logger = get_logger(str(uuid.uuid4()))
            logger.configure(logfile_path=temp_log_file, logfile_level=3)
            for i in range(10):
                logger.log(f"log line {i}", i)

            with open(temp_log_file) as tf:
                assert (
                    tf.read()
                    == """log line 0
log line 1
log line 2
log line 3
"""
                )
            logger.configure()  # close logfile

    def test_multiple_loggers_cannot_write_to_same_file(self):
        with temp_file_path(suffix=".log") as temp_log_file:
            logger1 = get_logger(str(uuid.uuid4()))
            logger2 = get_logger(str(uuid.uuid4()))
            logger1.configure(logfile_path=temp_log_file, logfile_level=2)
            with pytest.raises(FileError):
                logger2.configure(logfile_path=temp_log_file, logfile_level=3)
            logger1.configure()  # close logfile
            logger2.configure()  # close logfile

    def test_multiple_loggers_can_write_to_stdout_and_different_files(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file1, temp_file_path(suffix=".log") as temp_log_file2:
                logger1 = get_logger(str(uuid.uuid4()))
                logger2 = get_logger(str(uuid.uuid4()))
                logger1.configure(2, 1, temp_log_file1)
                logger2.configure(3, 2, temp_log_file2)

                for i in range(5):
                    logger1.log(f"From 1, level {i}", i)
                    logger2.log(f"From 2, level {i}", i)

                assert (
                    mock_stdout.getvalue()
                    == """From 1, level 0
From 2, level 0
From 1, level 1
From 2, level 1
From 1, level 2
From 2, level 2
From 2, level 3
"""
                )

                with open(temp_log_file1) as tf1:
                    assert (
                        tf1.read()
                        == """From 1, level 0
From 1, level 1
"""
                    )

                with open(temp_log_file2) as tf2:
                    assert (
                        tf2.read()
                        == """From 2, level 0
From 2, level 1
From 2, level 2
"""
                    )
                logger1.configure()  # close logfile
                logger2.configure()  # close logfile

    def test_same_logger_can_switch_write_files(self):
        with temp_file_path(suffix=".log") as temp_log_file1, temp_file_path(suffix=".log") as temp_log_file2:
            logger = get_logger(str(uuid.uuid4()))
            logger.configure(logfile_path=temp_log_file1, logfile_level=2)

            for i in range(5):
                logger.log(f"To 1, level {i}", i)

            logger.configure(logfile_path=temp_log_file2, logfile_level=1)

            for i in range(5):
                logger.log(f"To 2, level {i}", i)

            logger.configure()

            for i in range(5):
                logger.log(f"To None, level {i}", i)

            logger.configure(logfile_path=temp_log_file1, logfile_level=2)
            for i in range(5):
                logger.log(f"To 1 again, level {i}", i)

            with open(temp_log_file1) as tf1:
                assert (
                    tf1.read()
                    == """To 1, level 0
To 1, level 1
To 1, level 2
To 1 again, level 0
To 1 again, level 1
To 1 again, level 2
"""
                )

            with open(temp_log_file2) as tf2:

                assert (
                    tf2.read()
                    == """To 2, level 0
To 2, level 1
"""
                )
            logger.configure()  # close logfile

    def test_configure_prefixes_date_and_or_time_for_stdout_and_file(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file:
                logger = get_logger(str(uuid.uuid4()))

                dts = [[tf1, tf2] for tf1 in [True, False] for tf2 in [True, False]]
                for d, t in dts:
                    logger.configure(4, 4, temp_log_file, d, t)
                    for i in range(3, 8):
                        logger.log(f"d={d}, t={t}, level={i}", i)

                logger.configure(4, 4, temp_log_file, True, False)

                with open(temp_log_file) as tf:
                    for i, (l1, l2) in enumerate(
                        zip(
                            [l for l in mock_stdout.getvalue().replace("\r\n", "\n").split("\n") if l],
                            [l for l in tf.read().replace("\r\n", "\n").split("\n") if l],
                        )
                    ):
                        date, time = dts[i // 2]
                        pattern = f"d={date}, t={time}, level={i % 2 + 3}"
                        if date and time:
                            pattern = "\[\d{2}\.\d{2}\|\d{2}:\d{2}:\d{2}\] " + pattern
                        elif date:
                            pattern = "\[\d{2}\.\d{2}\] " + pattern
                        elif time:
                            pattern = "\[\d{2}:\d{2}:\d{2}\] " + pattern
                        assert re.fullmatch(pattern, l1) is not None and re.fullmatch(pattern, l2) is not None

                logger.configure()  # close logfile

    def test_thread_safe(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file1, temp_file_path(suffix=".log") as temp_log_file2:
                name = str(uuid.uuid4())
                num_threads = 10
                num_msgs = 1000

                def log(id):
                    # Introduce random variation into when threads start
                    time.sleep(random.uniform(0.0, 0.05))
                    logger = get_logger(name)
                    logger.configure(5, 5, temp_log_file1)
                    for i in range(num_msgs):
                        logger.configure(
                            5, 5, temp_log_file1 if i % 2 == 0 else temp_log_file2, i % 5 == 0, i % 11 == 0
                        )
                        logger.log(f"id {id} msg {i}", 5)

                threads = [threading.Thread(target=log, args=(i,)) for i in range(num_threads)]

                for thread in threads:
                    thread.start()
                for thread in threads:
                    thread.join()

                get_logger(name).configure()  # close logfile

                assert len(mock_stdout.getvalue().replace("\r\n", "\n").split("\n")) == num_threads * num_msgs + 1

                with open(temp_log_file1) as tf1, open(temp_log_file2) as tf2:
                    assert (
                        len(tf1.read().replace("\r\n", "\n").split("\n"))
                        + len(tf2.read().replace("\r\n", "\n").split("\n"))
                        == num_threads * num_msgs + 2
                    )

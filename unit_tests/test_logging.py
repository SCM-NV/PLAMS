import os.path
import re
import uuid
from io import StringIO
from unittest.mock import patch
import threading
import pytest
import time
import random

from scm.plams.core.errors import FileError, PlamsError
from scm.plams.core.logging import get_logger, TextLogger, CSVLogger
from scm.plams.core.formatters import JobCSVFormatter
from scm.plams.unit_tests.test_helpers import temp_file_path
from scm.plams.unit_tests.test_basejob import DummySingleJob


class TestGetLogger:

    def test_get_logger_returns_existing_or_creates_new(self):
        name1 = str(uuid.uuid4())
        name2 = str(uuid.uuid4())
        logger1 = get_logger(name1)
        logger2 = get_logger(name2)
        logger3 = get_logger(name1)

        assert logger1 == logger3 != logger2

    def test_get_logger_returns_logger_of_correct_type(self):
        logger1 = get_logger(str(uuid.uuid4()))
        logger2 = get_logger(str(uuid.uuid4()), "txt")
        logger3 = get_logger(str(uuid.uuid4()), "csv")

        assert isinstance(logger1, TextLogger)
        assert isinstance(logger2, TextLogger)
        assert isinstance(logger3, CSVLogger)

    def test_get_logger_errors_with_unsupported_format(self):
        with pytest.raises(PlamsError):
            get_logger("foo", "bar")


class TestTextLogger:

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
            logger.close()

    def test_configure_does_not_error_with_invalid_logfile_location(self):
        log_file = "not/a/file"
        logger = get_logger(str(uuid.uuid4()))
        logger.configure(logfile_path=log_file, logfile_level=3)
        logger.log("log line", 3)

        logger.close()

    def test_configure_only_updates_logfile_handler_when_abspath_changes(self):
        with temp_file_path(suffix=".log") as temp_log_file1:
            with temp_file_path(suffix=".log") as temp_log_file2:
                rel_path1 = temp_log_file1.split("/")[-1]
                rel_path2 = temp_log_file2.split("/")[-1]

                logger = get_logger(str(uuid.uuid4()))
                logger.configure(logfile_path=rel_path1, logfile_level=3)
                fh1 = logger._file_handler
                logger.configure(logfile_path=rel_path1, logfile_level=3)
                fh2 = logger._file_handler
                logger.configure(logfile_path=rel_path2, logfile_level=3)
                fh3 = logger._file_handler
                logger.configure(logfile_path=rel_path2, logfile_level=3)
                fh4 = logger._file_handler

                assert fh1 == fh2 != fh3 == fh4

                logger.close()

    def test_close_removes_stdout_and_logfile_handlers(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file:
                logger = get_logger(str(uuid.uuid4()))
                logger.configure(logfile_path=temp_log_file, logfile_level=3, stdout_level=3)
                for i in range(10):
                    logger.log(f"log line {i}", i)

                logger.close()

                for i in range(10):
                    logger.log(f"log line {i}", i)

                with open(temp_log_file) as tf:
                    assert (
                        tf.read()
                        == mock_stdout.getvalue()
                        == """log line 0
log line 1
log line 2
log line 3
"""
                    )

    def test_multiple_loggers_cannot_write_to_same_file(self):
        with temp_file_path(suffix=".log") as temp_log_file:
            logger1 = get_logger(str(uuid.uuid4()))
            logger2 = get_logger(str(uuid.uuid4()))
            logger1.configure(logfile_path=temp_log_file, logfile_level=2)
            with pytest.raises(FileError):
                logger2.configure(logfile_path=temp_log_file, logfile_level=3)
            logger1.close()
            logger2.close()

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
                logger1.close()
                logger2.close()

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
            logger.close()

    def test_configure_prefixes_date_and_or_time_for_stdout_and_file(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file:
                logger = get_logger(str(uuid.uuid4()))

                dts = [[tf1, tf2] for tf1 in [True, False] for tf2 in [True, False]]
                for d, t in dts:
                    logger.configure(4, 4, temp_log_file, d, t)
                    for i in range(3, 8):
                        logger.log(f"d={d}, t={t}, level={i}", i)

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

                logger.close()

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

                get_logger(name).close()

                assert len(mock_stdout.getvalue().replace("\r\n", "\n").split("\n")) == num_threads * num_msgs + 1

                with open(temp_log_file1) as tf1, open(temp_log_file2) as tf2:
                    assert (
                        len(tf1.read().replace("\r\n", "\n").split("\n"))
                        + len(tf2.read().replace("\r\n", "\n").split("\n"))
                        == num_threads * num_msgs + 2
                    )


class TestCSVLogger:

    def test_no_logging_to_stdout_by_default(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            logger = get_logger(str(uuid.uuid4()), "csv")
            logger.log("hello", 1)

            assert mock_stdout.getvalue() == ""

    def test_configure_writes_to_stdout_up_to_and_including_level(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            logger = get_logger(str(uuid.uuid4()), "csv")
            logger.configure(3)
            for i in range(10):
                logger.log(f"log line {i}", i)

            assert (
                mock_stdout.getvalue()
                == """message
log line 0
log line 1
log line 2
log line 3
"""
            )

    def test_configure_writes_to_logfile_up_to_and_including_level(self):
        with temp_file_path(suffix=".log") as temp_log_file:
            logger = get_logger(str(uuid.uuid4()), "csv")
            logger.configure(logfile_path=temp_log_file, logfile_level=3)
            for i in range(10):
                logger.log(f"log line {i}", i)

            with open(temp_log_file) as tf:
                assert (
                    tf.read()
                    == """message
log line 0
log line 1
log line 2
log line 3
"""
                )
            logger.close()

    def test_close_removes_stdout_and_logfile_handlers(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file:
                logger = get_logger(str(uuid.uuid4()), "csv")
                logger.configure(logfile_path=temp_log_file, logfile_level=3, stdout_level=3)
                for i in range(10):
                    logger.log(f"log line {i}", i)

                logger.close()

                for i in range(10):
                    logger.log(f"log line {i}", i)

                with open(temp_log_file) as tf:
                    assert (
                        tf.read()
                        == mock_stdout.getvalue()
                        == """message
log line 0
log line 1
log line 2
log line 3
"""
                    )

    def test_multiple_loggers_cannot_write_to_same_file(self):
        with temp_file_path(suffix=".log") as temp_log_file:
            logger1 = get_logger(str(uuid.uuid4()), "csv")
            logger2 = get_logger(str(uuid.uuid4()), "csv")
            logger1.configure(logfile_path=temp_log_file, logfile_level=2)
            with pytest.raises(FileError):
                logger2.configure(logfile_path=temp_log_file, logfile_level=3)
            logger1.close()
            logger2.close()

    def test_multiple_loggers_can_write_to_stdout_and_different_files(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file1, temp_file_path(suffix=".log") as temp_log_file2:
                logger1 = get_logger(str(uuid.uuid4()), "csv")
                logger2 = get_logger(str(uuid.uuid4()), "csv")
                logger1.configure(2, 1, temp_log_file1)
                logger2.configure(3, 2, temp_log_file2)

                for i in range(5):
                    logger1.log(f"From 1, level {i}", i)
                    logger2.log(f"From 2, level {i}", i)

                assert (
                    mock_stdout.getvalue()
                    == """message
"From 1, level 0"
message
"From 2, level 0"
"From 1, level 1"
"From 2, level 1"
"From 1, level 2"
"From 2, level 2"
"From 2, level 3"
"""
                )

                with open(temp_log_file1) as tf1:
                    assert (
                        tf1.read()
                        == """message
"From 1, level 0"
"From 1, level 1"
"""
                    )

                with open(temp_log_file2) as tf2:
                    assert (
                        tf2.read()
                        == """message
"From 2, level 0"
"From 2, level 1"
"From 2, level 2"
"""
                    )
                logger1.close()
                logger2.close()

    def test_same_logger_can_switch_write_files_and_only_writes_headers_once(self):
        with temp_file_path(suffix=".log") as temp_log_file1, temp_file_path(suffix=".log") as temp_log_file2:
            logger = get_logger(str(uuid.uuid4()), "csv")
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
                    == """message
"To 1, level 0"
"To 1, level 1"
"To 1, level 2"
"To 1 again, level 0"
"To 1 again, level 1"
"To 1 again, level 2"
"""
                )

            with open(temp_log_file2) as tf2:

                assert (
                    tf2.read()
                    == """message
"To 2, level 0"
"To 2, level 1"
"""
                )
            logger.close()

    def test_configure_adds_date_and_or_time_fields_to_logging(self):
        with temp_file_path(suffix=".log") as temp_log_file:
            logger = get_logger(str(uuid.uuid4()), "csv")

            dts = [[tf1, tf2] for tf1 in [True, False] for tf2 in [True, False]]
            for d, t in dts:
                logger.configure(4, 4, temp_log_file, d, t, include_level=True)
                for i in range(3, 8):
                    logger.log({"d": d, "t": t}, i)

            with open(temp_log_file) as tf:
                for i, l in enumerate([l for l in tf.read().replace("\r\n", "\n").split("\n") if l][1:]):
                    date, time = dts[i // 2]
                    pattern = f"{date},{time}"
                    if date and time:
                        pattern = "[3-4],\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}," + pattern
                    elif date:
                        pattern = "[3-4],\d{4}-\d{2}-\d{2}," + pattern
                    elif time:
                        pattern = "[3-4],\d{2}:\d{2}:\d{2}," + pattern
                    else:
                        pattern = "[3-4]," + pattern
                    assert re.fullmatch(pattern, l)

                logger.close()

    def test_configure_adds_fields_to_logging(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file:
                name = str(uuid.uuid4())
                logger = get_logger(name, "csv")
                logger.configure(3, 3, temp_log_file, include_level=True, include_name=True)

                for i in range(10):
                    logger.log(f"log line {i}", i)

                with open(temp_log_file) as tf:
                    assert (
                        tf.read()
                        == mock_stdout.getvalue()
                        == f"""logger_name,level,message
{name},0,log line 0
{name},1,log line 1
{name},2,log line 2
{name},3,log line 3
"""
                    )
                logger.close()

    def test_log_gets_field_names_from_dictionary(self):
        with temp_file_path(suffix=".log") as temp_log_file:
            name = str(uuid.uuid4())
            logger = get_logger(name, "csv")
            logger.configure(3, 3, temp_log_file)

            for i in range(3):
                logger.log({"foo": f"bar{i}", "fizz": f"buzz{i}"}, 3)

            with open(temp_log_file) as tf:
                assert (
                    tf.read()
                    == """foo,fizz
bar0,buzz0
bar1,buzz1
bar2,buzz2
"""
                )
            logger.close()

    def test_log_cannot_change_field_names_or_add_fields(self):
        with patch("sys.stderr", new_callable=StringIO) as mock_stderr:
            with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                name = str(uuid.uuid4())
                logger = get_logger(name, "csv")
                logger.configure(3)

                logger.log({"foo": "bar", "fizz": "buzz"}, 3)
                logger.log({"abc": "bar2", "def": "buzz2"}, 3)
                logger.log({"foo": "bar3", "fizz": "buzz3", "abc": "def"}, 3)

            assert "ValueError: dict contains fields not in fieldnames" in mock_stderr.getvalue()
            assert mock_stdout.getvalue() == (
                """foo,fizz
bar,buzz
"""
            )

            logger.close()

    def test_log_correctly_escapes_commas_and_multiline_strings(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            name = str(uuid.uuid4())
            logger = get_logger(name, "csv")
            logger.configure(3)

            logger.log(
                {
                    "commas": "a,b,c,'(d,e,f)'",
                    "multi-line": """See the following:
    a,b,c
    d,e,f
""",
                },
                3,
            )
            logger.log(
                {
                    "commas": 'm,"n,",,,o\'',
                    "multi-line": """See the following:
    m,n,o

    p,q,r
""",
                },
                3,
            )

        assert mock_stdout.getvalue() == (
            """commas,multi-line
"a,b,c,\'(d,e,f)\'","See the following:
    a,b,c
    d,e,f
"
"m,""n,"",,,o\'","See the following:
    m,n,o

    p,q,r
"
"""
        )

        logger.close()


class TestJobCSVFormatter:

    def test_formatter_populates_job_fields(self):
        with temp_file_path(suffix=".log") as temp_log_file:

            job1 = DummySingleJob(name="test_job_csv_formatter")
            job2 = DummySingleJob(name="test_job_csv_formatter.002", cmd="err")

            def get_errormsg():
                return "some error"

            setattr(job2, "get_errormsg", get_errormsg)

            logger = get_logger(str(uuid.uuid4()), "csv")
            logger.configure(logfile_level=7, csv_formatter=JobCSVFormatter, logfile_path=temp_log_file)

            logger.log(job1, 3)
            logger.log(job2, 3)

            job1.run()
            job2.run()

            logger.log(job1, 3)
            logger.log(job2, 3)

            dir = os.getcwd()
            path1 = os.path.join(dir, "plams_workdir", "test_job_csv_formatter")
            path2 = os.path.join(dir, "plams_workdir", "test_job_csv_formatter.002")

            with open(temp_log_file) as tf:
                assert (
                    tf.read()
                    == f"""job_base_name,job_name,job_status,job_path,job_ok,job_check,job_get_errormsg,job_parent_name,job_parent_path
test_job_csv_formatter,test_job_csv_formatter,created,,,,,,
test_job_csv_formatter,test_job_csv_formatter.002,created,,,,,,
test_job_csv_formatter,test_job_csv_formatter,successful,{path1},True,True,,,
test_job_csv_formatter,test_job_csv_formatter.002,crashed,{path2},False,False,some error,,
"""
                )

            logger.close()

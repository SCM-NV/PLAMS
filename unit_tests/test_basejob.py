import pytest
import uuid
from unittest.mock import patch
from datetime import datetime
from collections import namedtuple
import shutil
import re
from io import StringIO
import csv
from functools import wraps
import time
from pathlib import Path

from scm.plams.core.settings import Settings
from scm.plams.core.basejob import SingleJob, MultiJob
from scm.plams.core.errors import PlamsError, FileError, ResultsError
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.jobmanager import JobManager
from scm.plams.core.functions import add_to_instance, jobs_in_directory
from scm.plams.core.enums import JobStatus

LogEntry = namedtuple("LogEntry", ["method", "args", "kwargs", "start", "end"])


def log_call(method):
    """
    Decorator to log calls to instance.
    """

    @wraps(method)
    def wrapper(self, *args, **kwargs):
        start = datetime.now()
        try:
            return method(self, *args, **kwargs)
        finally:
            end = datetime.now()
            entry = LogEntry(method=method.__name__, args=args, kwargs=kwargs, start=start, end=end)
            self.add_call_log_entry(entry)

    return wrapper


class DummySingleJob(SingleJob):
    """
    Dummy Single Job for testing PLAMS components.
    Calls to methods are logged in order, with passed arguments and a timestamp.
    """

    def __init__(self, inp: str = None, cmd: str = None, wait: float = 0.0, **kwargs):
        """
        Initialize new dummy single job instance. Each job will have a unique id.
        By default, a job will run a simple sed command on an input string containing the unique id,
        which will transform it to the results in the output file.

        :param inp: input string for input file, where any %ID% value will be replaced by a unique job id
        :param cmd: command to execute on the input file
        :param wait: delay before executing command
        :param kwargs: kwargs for base single job
        """
        super().__init__(**kwargs)
        self._call_log = {}
        self.id = uuid.uuid4()
        self.input = inp.replace("%ID%", str(self.id)) if inp is not None else f"Dummy input {self.id}"
        self.command = cmd if cmd is not None else "sed 's/input/output/g'"
        self.wait = wait

    def add_call_log_entry(self, entry):
        if entry.method in self._call_log:
            if isinstance(self._call_log[entry.method], list):
                self._call_log[entry.method].append(entry)
            else:
                self._call_log[entry.method] = [self._call_log[entry.method], entry]
        else:
            self._call_log[entry.method] = entry

    def get_call_log_entries(self, method):
        # Try with delay to avoid race condition in entry being added
        for _ in range(4):
            try:
                return self._call_log[method]
            except KeyError:
                time.sleep(0.1)
        return self._call_log[method]

    @log_call
    def prerun(self) -> None:
        pass

    @log_call
    def postrun(self) -> None:
        pass

    @log_call
    def get_input(self) -> str:
        return self.input

    @log_call
    def get_runscript(self) -> str:
        return f"sleep {self.wait} && {self.command} {self._filename('inp')}"

    def _prepare(self, jobmanager) -> bool:
        # Decorator log_call causes pickling issues here for parent override so just do it manually here...
        start = datetime.now()
        try:
            return super()._prepare(jobmanager)
        finally:
            end = datetime.now()
            entry = LogEntry(method="_prepare", args=[], kwargs={}, start=start, end=end)
            self.add_call_log_entry(entry)

    @log_call
    def _execute(self, jobrunner) -> None:
        super()._execute(jobrunner)

    @log_call
    def _finalize(self) -> None:
        super()._finalize()

    def check(self) -> bool:
        try:
            return self.results.read_file(self._filename("err")) == ""
        except ResultsError:
            return True

    def get_errormsg(self) -> str:
        if self._error_msg:
            return self._error_msg

        msg = self.results.read_file(self._filename("err"))
        return msg if msg else None


class TestSingleJob:
    """
    Test suite for the Single Job.
    Not truly independent as relies upon the job runner/manager and results components.
    But this suite focuses on testing the methods on the job class itself.
    """

    def test_full_runscript_applies_pre_and_post_runscript_settings(self):
        # Given runscript settings
        s = Settings()
        s.runscript.shebang = "#!/bin/sh"
        s.runscript.pre = "# Add pre line"
        s.runscript.post = "\n# Add post line"

        # When get the runscript from the job
        job = DummySingleJob(settings=s)
        runscript = job.full_runscript()

        # Then the script is composed as expected
        assert (
            runscript
            == """#!/bin/sh

# Add pre line

sleep 0.0 && sed 's/input/output/g' plamsjob.in
# Add post line

""".replace(
                "\r\n", "\n"
            )
        )

    def test_run_wrapped_in_pre_and_postrun(self):
        # When run a job
        job = DummySingleJob()
        job.run()

        # Then makes calls to pre- and post-run
        assert (
            job.get_call_log_entries("prerun").end
            <= job.get_call_log_entries("_execute").start
            <= job.get_call_log_entries("_execute").end
            <= job.get_call_log_entries("postrun").start
        )

    def test_run_writes_output_file_on_success(self):
        # When run a successful job
        job = DummySingleJob()
        results = job.run()

        # Then job check passes and output file written
        assert job.check()
        assert job.ok()
        assert job.status == JobStatus.SUCCESSFUL
        assert results.read_file("$JN.in") == f"Dummy input {job.id}"
        assert results.read_file("$JN.out") == f"Dummy output {job.id}"
        assert results.read_file("$JN.err") == ""

    def test_run_writes_error_file_on_failure(self):
        # When run an erroring job
        job = DummySingleJob(cmd="not_a_cmd")
        results = job.run()

        # Then job check fails and error file written
        assert not job.check()
        assert not job.ok()
        assert job.status == JobStatus.CRASHED
        assert results.read_file("$JN.in") == f"Dummy input {job.id}"
        assert results.read_file("$JN.out") == ""
        assert results.read_file("$JN.err") != ""

    def test_run_marks_jobs_as_failed_and_stores_exception_on_prerun_or_postrun_error(self):
        # Given job which error in pre- or post-run
        job1 = DummySingleJob()
        job2 = DummySingleJob()

        @add_to_instance(job1)
        def prerun(s):
            raise RuntimeError("something went wrong")

        @add_to_instance(job2)
        def postrun(s):
            raise RuntimeError("something went wrong")

        # When run jobs
        job1.run()
        job2.run()

        # Then status is marked as failed
        assert job1.status == JobStatus.FAILED
        assert job2.status == JobStatus.FAILED
        assert not job1.ok()
        assert not job2.ok()
        assert "self.prerun()" in job1.get_errormsg()
        assert "RuntimeError: something went wrong" in job1.get_errormsg()
        assert "self.postrun()" in job2.get_errormsg()
        assert "RuntimeError: something went wrong" in job2.get_errormsg()

    def test_run_marks_jobs_as_failed_and_stores_exception_on_execution_error(self):
        # Given job which errors in execution
        job1 = DummySingleJob()

        def filename_errors(t):
            import inspect

            s = inspect.stack()
            caller = s[1].function
            if caller == "_execute":
                raise RuntimeError("something went wrong")
            else:
                return job1._filenames[t].replace("$JN", job1.name)

        job1._filename = filename_errors

        # When run job
        job1.run()

        # Then status is marked as failed
        assert job1.status == JobStatus.FAILED
        assert not job1.ok()
        assert "_execute" in job1.get_errormsg()
        assert "RuntimeError: something went wrong" in job1.get_errormsg()

    @pytest.mark.parametrize(
        "mode,expected",
        [
            [None, {(1, 1): None, (1, 2): None, (1, 3): None, (2, 2): None, (2, 3): None, (3, 3): None}],
            ["input", {(1, 1): True, (1, 2): False, (1, 3): True, (2, 2): True, (2, 3): False, (3, 3): True}],
            ["runscript", {(1, 1): True, (1, 2): True, (1, 3): False, (2, 2): True, (2, 3): False, (3, 3): True}],
            [
                "input+runscript",
                {(1, 1): True, (1, 2): False, (1, 3): False, (2, 2): True, (2, 3): False, (3, 3): True},
            ],
            ["not_a_mode", None],
        ],
        ids=["no_hashing", "input_hashing", "runscript_hashing", "both_hashing", "invalid_hashing"],
    )
    def test_hash_respects_mode(self, mode, expected, config):
        # Given jobs with different inputs and/or runscripts
        config.jobmanager.hashing = mode
        s = Settings()
        s.runscript.shebang = "#!/bin/sh"
        job1 = DummySingleJob(settings=s)
        job2 = DummySingleJob(inp=job1.input.replace("input", "inputx"), settings=s)
        job3 = DummySingleJob(inp=job1.input, cmd="echo 'foo' && sed 's/input/output/g'", settings=s)
        jobs = [job1, job2, job3]

        # When call hash with different modes
        # Then hashes match as expected
        if expected is None:
            with pytest.raises(PlamsError):
                job1.hash()
        else:
            hashes = [
                (i + 1, j + 1, job_i.hash(), job_j.hash())
                for i, job_i in enumerate(jobs)
                for j, job_j in enumerate(jobs)
                if j >= i
            ]
            matches = {(i, j): None if h_i is None and h_j is None else h_i == h_j for i, j, h_i, h_j in hashes}
            assert matches == expected

    def test_run_multiple_independent_jobs_in_parallel(self):
        # Given parallel job runner
        runner = JobRunner(parallel=True, maxjobs=2)

        # When set up two jobs with no dependencies
        job1 = DummySingleJob(wait=0.5)
        job2 = DummySingleJob(wait=0.01)
        results = [job1.run(runner), job2.run(runner)]

        # Then both run in parallel
        # Shorter job finishes first even though started second
        # Execute call of second job is made before finalize call of first job
        results[1].wait()
        assert job2.status == JobStatus.SUCCESSFUL
        assert job1.status == JobStatus.RUNNING

        results[0].wait()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.get_call_log_entries("_execute").start < job1.get_call_log_entries("_finalize").start

    def test_run_multiple_dependent_jobs_in_serial(self):
        # Given parallel job runner
        runner = JobRunner(parallel=True, maxjobs=2)

        # When set up two jobs with dependency
        job1 = DummySingleJob(wait=0.2)
        job2 = DummySingleJob(wait=0.01, depend=[job1])
        results = [job1.run(runner), job2.run(runner)]

        # Then run in serial
        # Second job finishes second even though shorter
        # Pre-run call of second job is made after post-run call of first job
        results[0].wait()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.status in [JobStatus.REGISTERED, JobStatus.STARTED, JobStatus.RUNNING]

        results[1].wait()
        assert job2.status == JobStatus.SUCCESSFUL
        assert job2.get_call_log_entries("prerun").start >= job1.get_call_log_entries("postrun").end

    def test_run_multiple_prerun_dependent_jobs_in_serial(self):
        # Given parallel job runner
        runner = JobRunner(parallel=True, maxjobs=2)

        # When set up two jobs with dependency via prerun
        job1 = DummySingleJob(wait=0.2)
        job2 = DummySingleJob(wait=0.01)

        @add_to_instance(job2)
        def prerun(s):
            job1.results.wait()

        results = [job1.run(runner), job2.run(runner)]

        # Then run in serial
        # Second job finishes second even though shorter
        # Execute call of second job is made after finalize call of first job
        results[0].wait()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.status in [JobStatus.REGISTERED, JobStatus.STARTED, JobStatus.RUNNING]

        results[1].wait()
        assert job2.status == JobStatus.SUCCESSFUL
        assert job2.get_call_log_entries("_execute").start >= job1.get_call_log_entries("_finalize").end

    def test_run_multiple_dependent_jobs_first_job_fails_in_prerun(self):
        # Given two dependent jobs where first errors in prerun
        job1 = DummySingleJob()
        job2 = DummySingleJob(depend=[job1])

        @add_to_instance(job1)
        def prerun(s):
            raise RuntimeError("something went wrong")

        # When run jobs
        job1.run()
        job2.run()

        # Then first job is marked as failed, the dependent job as successful
        assert job1.status == JobStatus.FAILED
        assert job2.status == JobStatus.SUCCESSFUL
        assert not job1.ok()
        assert job2.ok()
        assert "RuntimeError" in job1.get_errormsg()
        assert job2.get_errormsg() is None

    def test_run_multiple_independent_jobs_with_use_subdir_in_parallel(self, config):
        # Given parallel job runner
        runner = JobRunner(parallel=True, maxjobs=4)

        # When run jobs in subdir
        jobs = []
        with jobs_in_directory("results"):
            for outer in ["A", "B", "C"]:
                with jobs_in_directory(outer):
                    for inner1 in ["X", "Y", "Z"]:
                        with jobs_in_directory(inner1):
                            for inner2 in range(2):
                                job = DummySingleJob(name=f"{outer}_{inner1}_{inner2}")
                                job.run(runner)
                                jobs.append(job)

        # Then jobs are located in the correct subdirectory
        assert len(jobs) == 18
        assert all(j.ok() for j in jobs)

        def normalise_path(p):
            s = str(p)
            if s.startswith(("\\\\?\\", "//?/")):
                s = s[4:]
            return Path(s).resolve()

        for j in jobs:
            o, i1, i2 = j.name.split("_")
            assert normalise_path(j.path) == normalise_path(
                Path(config.default_jobmanager.workdir, "results", o, i1, j.name)
            )

    def test_ok_waits_on_results_and_checks_status(self):
        # Given job and a copy
        job1 = DummySingleJob(wait=0.1)
        job2 = DummySingleJob(inp=job1.input)
        job3 = DummySingleJob(cmd="not_a_cmd")

        # Then not ok before being started
        assert not job1.ok()

        # When call run and check ok
        job1.run()
        job2.run()
        job3.run()

        # Then waits for job to finish and checks status
        assert job1.ok()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.ok()
        assert job2.status == JobStatus.COPIED
        assert not job3.ok()
        assert job3.status == JobStatus.CRASHED

    def test_load_non_existent_job_errors(self):
        # Given path that does not point to a job
        # When call load, then gives error
        with pytest.raises(FileError):
            DummySingleJob.load("not_a_path")

    def test_load_previously_run_job_succeeds(self):
        # Given successful previously run job
        job1 = DummySingleJob()
        job1.run()
        job1.results.wait()

        # When call load
        job2 = DummySingleJob.load(job1.path)

        # Then job loaded successfully
        assert job1.name == job2.name
        assert job1.id == job2.id
        assert job1.path == job2.path
        assert job1.settings == job2.settings
        assert job1._filenames == job2._filenames

    def test_load_legacy_job_succeeds(self):
        # Given job run before additional properties were added
        job1 = DummySingleJob()
        job1.run()
        job1.results.wait()
        status = job1.status
        delattr(job1, "_status")
        delattr(job1, "_status_log")
        job1.__dict__["status"] = status
        job1.pickle()

        # When call load
        job2 = DummySingleJob.load(job1.path)

        # Then job loaded successfully
        assert job1.name == job2.name
        assert job1.id == job2.id
        assert job1.path == job2.path
        assert job1.settings == job2.settings
        assert job1._filenames == job2._filenames
        assert job2.status == "successful"
        assert job2.status_log == []
        assert job2.get_errormsg() is None

    def test_job_summaries_logged(self, config):
        job1 = DummySingleJob()
        job2 = DummySingleJob(inp=job1.input)
        job3 = DummySingleJob(cmd="not_a_cmd")

        job_manager = JobManager(config.jobmanager, folder="test_logging")

        job1.run(jobmanager=job_manager)
        job2.run(jobmanager=job_manager)
        job3.run(jobmanager=job_manager)

        dt_fmt = r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}"

        def assert_csv_entry(entry, statuses, postfix="", ok="True", check="True", error_msg=""):
            assert re.match(dt_fmt, entry[0])
            assert entry[1] == "plamsjob"
            assert entry[2] == f"plamsjob{postfix}"
            assert entry[3] == statuses[-1]
            assert entry[4].endswith(f"plamsjob{postfix}")
            assert entry[5] == ok
            assert entry[6] == check
            assert error_msg in entry[7]
            status_pattern = str.join(" -> ", [rf"{dt_fmt} {s}" for s in statuses])
            assert re.match(status_pattern, entry[8])

        with open(job_manager.job_logger.logfile) as f:
            reader = csv.reader(f)

            assert next(reader) == [
                "logged_at",
                "job_base_name",
                "job_name",
                "job_status",
                "job_path",
                "job_ok",
                "job_check",
                "job_get_errormsg",
                "job_timeline",
                "job_parent_name",
                "job_parent_path",
            ]

            assert_csv_entry(next(reader), ["created", "started", "registered", "running", "finished", "successful"])
            assert_csv_entry(next(reader), ["created", "started", "registered", "copied"], postfix=".002")
            assert_csv_entry(
                next(reader),
                ["created", "started", "registered", "running", "crashed"],
                postfix=".003",
                error_msg="not_a_cmd",
                ok="False",
                check="False",
            )

        job_manager.job_logger.close()
        shutil.rmtree(job_manager.workdir)

    def test_job_errors_logged_to_stdout(self, config):
        from scm.plams.core.logging import get_logger

        logger_name = f"plams-{uuid.uuid4()}"
        logger = get_logger(logger_name)

        with patch("scm.plams.core.functions._logger", logger):
            with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                job1 = DummySingleJob(cmd="not_a_cmd")
                job2 = DummySingleJob(cmd="x\n" * 50)

                job1.run()
                job2.run()

                stdout = mock_stdout.getvalue()
                assert re.match(
                    f".*Error message for job {job1.name} was:.* 3: not_a_cmd: (command ){{0,1}}not found",
                    stdout,
                    re.DOTALL,
                )
                assert re.match(
                    f".*Error message for job {job2.name} was:.* 3: x: (command ){{0,1}}not found.* 32: x: (command ){{0,1}}not found.*(see output for full error)",
                    stdout,
                    re.DOTALL,
                )
        logger.close()

    def test_full_name(self):
        job = DummySingleJob(name="dummy_job")
        job.run()

        assert job.ok()
        assert job._full_name() == "dummy_job"
        assert job._full_name("some/rundir") == "some/rundir/dummy_job"

    def test_delete_created_job(self, config):
        # Given job
        job = DummySingleJob()

        # When deleted
        job.delete()

        # Then status set to deleted
        assert job.status == JobStatus.DELETED
        assert job.path is None
        with pytest.raises(ResultsError):
            _ = job.results.grep_output("")

    def test_delete_running_job(self, config):
        # Given job
        job = DummySingleJob(wait=0.5)
        job.run()
        path = job.path

        # When deleted while running
        job.delete()

        # Then waits, job files removed and job path removed
        assert job.status == JobStatus.DELETED
        assert job.name not in config.default_jobmanager.names
        assert job not in config.default_jobmanager.jobs
        assert not Path(path).exists()
        assert job.path is None
        with pytest.raises(ResultsError):
            _ = job.results.grep_output("")

    def test_delete_job_then_rerun_with_same_name(self, config):
        # Given job
        name = f"to-be-deleted-{uuid.uuid4()}"
        job1 = DummySingleJob(name=name)
        job1.run()
        path1 = job1.path

        # When deleted
        job1.delete()

        # Then can rerun job with the same name, and results cannot be accessed from first job
        job2 = DummySingleJob(name=name)
        job2.run()

        assert job1.name == job2.name
        assert job1.path is None
        with pytest.raises(ResultsError):
            job1.results.grep_output("")
        assert job2.path == path1
        assert job2.results.grep_output("")

    def test_delete_many_jobs_then_rerun_with_same_name(self, config):
        # Given jobs with same name
        name = f"to-be-deleted-{uuid.uuid4()}"
        jobs = [DummySingleJob(name=name) for _ in range(10)]
        for job in jobs:
            job.run()
        for job in jobs:
            job.results.wait()
        path = jobs[0].path

        # When jobs which are not the final job are deleted
        jobs[0].delete()
        jobs[4].delete()
        jobs[7].delete()

        # Then next job follows on from last job
        jobs.append(DummySingleJob(name=name))
        jobs[-1].run().wait()
        assert int(jobs[-1].name[-3:]) == 11

        # When last jobs are deleted
        jobs[-2].delete()
        jobs[-1].delete()

        # Then next job follows on from highest remaining job
        jobs.append(DummySingleJob(name=name))
        jobs[-1].run().wait()
        assert int(jobs[-1].name[-3:]) == 10

        # When all other jobs are deleted
        jobs[1].delete()
        jobs[2].delete()
        jobs[3].delete()
        jobs[5].delete()
        jobs[-1].delete()
        jobs[6].delete()
        jobs[8].delete()

        # Then count is fully reset and next job has the base name
        jobs.append(DummySingleJob(name=name))
        jobs[-1].run().wait()
        last_job = jobs.pop()
        assert last_job.name == name
        assert last_job.results.grep_output("")
        assert jobs[0].path is None
        with pytest.raises(ResultsError):
            jobs[0].results.grep_output("")
        assert all(job.status == JobStatus.DELETED for job in jobs)
        assert len(list(Path(path).parent.glob(f"{name}*"))) == 1

    def test_delete_many_jobs_with_same_name_in_different_sub_dirs(self, config):
        # Given jobs with same name in different sub dirs
        name = f"to-be-deleted-{uuid.uuid4()}"
        with jobs_in_directory("dir1") as dir1:
            jobs1 = [DummySingleJob(name=name) for _ in range(5)]
            for j in jobs1:
                j.run()
        with jobs_in_directory("dir2") as dir2:
            jobs2 = [DummySingleJob(name=name) for _ in range(5)]
            for j in jobs2:
                j.run()
            with jobs_in_directory("dir3") as dir3:
                jobs3 = [DummySingleJob(name=name) for _ in range(5)]
                for j in jobs3:
                    j.run()
        for j in jobs1 + jobs2 + jobs3:
            j.results.wait()

        # When jobs deleted from a subdir
        for job in jobs2[::-1]:
            job.delete()

        # Then jobs in other dirs unaffected
        assert len(config.default_jobmanager.names) == 2
        assert all(job.status == JobStatus.SUCCESSFUL for job in jobs1)
        assert all(job.status == JobStatus.DELETED for job in jobs2)
        assert all(job.status == JobStatus.SUCCESSFUL for job in jobs3)
        with pytest.raises(ResultsError):
            jobs2[0].results.grep_output("")
        assert jobs1[0].results.grep_output("")
        assert jobs3[0].results.grep_output("")
        assert len(list(dir2.glob(f"{name}*"))) == 0
        assert len(list(dir1.glob(f"{name}*"))) == 5
        assert len(list(dir3.glob(f"{name}*"))) == 5

        # When jobs deleted from a subdir
        for job in jobs1:
            job.delete()

        # Then jobs in other dirs unaffected
        assert len(config.default_jobmanager.names) == 1
        assert all(job.status == JobStatus.DELETED for job in jobs1)
        assert all(job.status == JobStatus.SUCCESSFUL for job in jobs3)

        # When jobs deleted from a subdir
        jobs3[1].delete()
        jobs3[3].delete()
        jobs3[4].delete()
        jobs3[2].delete()
        jobs3[0].delete()

        # Then jobs in other dirs unaffected
        assert len(config.default_jobmanager.names) == 0

    def test_rename_created_job(self, config):
        # Given job
        id = uuid.uuid4()
        name1 = f"to-be-renamed-{id}"
        name2 = f"renamed-{id}"
        job = DummySingleJob(name=name1)

        # When renamed
        job.rename(name2)

        # Then name changed
        assert job.status == JobStatus.CREATED
        assert job.path is None
        assert job.name == name2

    def test_rename_running_job(self, config):
        # Given job
        id = uuid.uuid4()
        name1 = f"to-be-renamed-{id}"
        name2 = f"renamed-{id}"
        job = DummySingleJob(name=name1, wait=0.5)
        job.run()
        path1 = job.path

        # When renamed while running
        job.rename(name2)
        path2 = job.path

        # Then waits, job files moved and renamed, re-registered in job manager
        assert job.status == JobStatus.SUCCESSFUL
        assert not Path(path1).exists()
        assert Path(path2).exists()
        assert Path(path2) == Path(path1).parent / name2
        assert name1 not in config.default_jobmanager.names
        assert name2 in config.default_jobmanager.names

        # And results can be read
        assert job.results.read_file("$JN.in")
        assert job.results.read_file("$JN.out")
        assert job.results.read_file("$JN.run")

        # And job can be loaded from .dill file
        loaded_job = DummySingleJob.load(str(Path(job.path) / f"{job.name}.dill"))
        assert loaded_job.name == name2

    def test_rename_job_with_same_name(self, config):
        # Given two jobs
        name1 = f"to-be-renamed-{uuid.uuid4()}"
        name2 = f"renamed-{uuid.uuid4()}"
        job1 = DummySingleJob(name=name1)
        job2 = DummySingleJob(name=name2)
        job1.run()
        job2.run()
        path1 = job1.path

        # When job renamed with same name as itself
        job1.rename(name1)

        # Then job and files unchanged
        assert job1.name == name1
        assert job1.path == path1

        # When job renamed with same name as another job
        job1.rename(name2)
        path2 = job1.path

        # Then job files moved and renamed, re-registered in job manager
        name3 = f"{name2}.002"
        assert not Path(path1).exists()
        assert Path(path2).exists()
        assert Path(path2) == Path(path1).parent / name3
        assert name1 not in config.default_jobmanager.names
        assert name2 in config.default_jobmanager.names
        assert config.default_jobmanager.names[name2] == 2

    def test_rename_many_jobs_in_different_sub_dirs(self, config):
        # Given jobs with same name in different sub dirs
        name1 = f"name1-{uuid.uuid4()}"
        name2 = f"name2-{uuid.uuid4()}"
        name3 = f"name3-{uuid.uuid4()}"
        with jobs_in_directory("dir1") as dir1:
            jobs1 = [DummySingleJob(name=name1) for _ in range(3)]
            for j in jobs1:
                j.run()
        with jobs_in_directory("dir2") as dir2:
            jobs2 = [DummySingleJob(name=name2) for _ in range(3)]
            for j in jobs2:
                j.run()
            with jobs_in_directory("dir3") as dir3:
                jobs3 = [DummySingleJob(name=name3) for _ in range(3)]
                for j in jobs3:
                    j.run()
        jobs = jobs1 + jobs2 + jobs3
        for j in jobs:
            j.results.wait()

        # When jobs renamed
        for i in range(3):
            for j, n in enumerate([name1, name2, name3]):
                jobs[i * 3 + j].rename(n)

        # Then jobs remain in the same subdirectory but are renamed
        for i, d in enumerate([dir1, dir2, dir3]):
            for j, n in enumerate([name1, name2, name3]):
                assert jobs[i * 3 + j].path.startswith(str(d / n))


class TestMultiJob:
    """
    Test suite for the Multi Job.
    Not truly independent as relies upon the job runner/manager and results components.
    But this suite focuses on testing the methods on the job class itself.
    """

    def test_run_multiple_independent_single_jobs_all_succeed(self, config):
        runner = JobRunner(parallel=True, maxjobs=3)
        config.sleepstep = 0.1

        # Given 3 jobs which are independent
        jobs = [DummySingleJob() for _ in range(3)]
        multi_job = MultiJob(children=jobs)

        # When run multi-job
        multi_job.run(jobrunner=runner).wait()

        # Then multi-job ran ok
        assert multi_job.check()
        assert multi_job.ok()
        assert multi_job.status == JobStatus.SUCCESSFUL
        assert all([j.status == JobStatus.SUCCESSFUL for j in jobs])

    def test_run_multiple_independent_single_jobs_one_fails(self, config):
        runner = JobRunner(parallel=True, maxjobs=3)
        config.sleepstep = 0.1

        # Given 3 jobs which are independent, one of which fails
        jobs = [DummySingleJob(), DummySingleJob(cmd="not_a_cmd"), DummySingleJob()]
        multi_job = MultiJob(children=jobs)

        # When run multi-job
        multi_job.run(jobrunner=runner).wait()

        # Then multi-job fails
        assert not multi_job.check()
        assert not multi_job.ok()
        assert multi_job.status == JobStatus.FAILED
        assert [j.status for j in jobs] == [JobStatus.SUCCESSFUL, JobStatus.CRASHED, JobStatus.SUCCESSFUL]

    def test_run_multiple_dependent_single_jobs_all_succeed(self, config):
        runner = JobRunner(parallel=True, maxjobs=3)
        config.sleepstep = 0.1

        # Given 3 jobs which are dependent
        jobs = [DummySingleJob() for _ in range(3)]

        @add_to_instance(jobs[1])
        def prerun(s):
            jobs[0].results.wait()

        @add_to_instance(jobs[2])
        def postrun(s):
            jobs[1].results.wait()

        multi_job = MultiJob(children=jobs)

        # When run multi-job
        multi_job.run(jobrunner=runner).wait()

        # Then multi-job ran ok
        assert multi_job.check()
        assert multi_job.ok()
        assert multi_job.status == JobStatus.SUCCESSFUL
        assert all([j.status == JobStatus.SUCCESSFUL for j in jobs])

    def test_run_multiple_independent_single_jobs_error_in_prerun_or_postrun(self, config):
        runner = JobRunner(parallel=True, maxjobs=3)
        config.sleepstep = 0.1

        # Given 3 jobs which are dependent
        jobs = [DummySingleJob() for _ in range(3)]

        @add_to_instance(jobs[1])
        def prerun(s):
            raise RuntimeError("something went wrong")

        @add_to_instance(jobs[2])
        def postrun(s):
            raise RuntimeError("something went wrong")

        multi_job = MultiJob(children=jobs)

        # When run multi-job
        multi_job.run(jobrunner=runner).wait()

        # Then multi-job failed
        assert not multi_job.check()
        assert not multi_job.ok()
        assert multi_job.status == JobStatus.FAILED
        assert [j.status for j in jobs] == [JobStatus.SUCCESSFUL, JobStatus.FAILED, JobStatus.FAILED]

    def test_run_multiple_independent_single_jobs_error_in_execute(self, config):
        runner = JobRunner(parallel=True, maxjobs=3)
        config.sleepstep = 0.1

        # Given 3 jobs which are dependent
        jobs = [DummySingleJob() for _ in range(3)]

        def filename_errors(t):
            import inspect

            s = inspect.stack()
            caller = s[1].function
            if caller == "_execute":
                raise RuntimeError("something went wrong")
            else:
                return jobs[1]._filenames[t].replace("$JN", jobs[1].name)

        jobs[1]._filename = filename_errors

        multi_job = MultiJob(children=jobs)

        # When run multi-job
        multi_job.run(jobrunner=runner).wait()

        # Then multi-job failed
        assert not multi_job.check()
        assert not multi_job.ok()
        assert multi_job.status == JobStatus.FAILED
        assert [j.status for j in jobs] == [JobStatus.SUCCESSFUL, JobStatus.FAILED, JobStatus.SUCCESSFUL]

    def test_run_multiple_independent_multijobs_all_succeed(self, config):
        runner = JobRunner(parallel=True, maxjobs=3)
        config.sleepstep = 0.1

        # Given multi-job with multiple multi-jobs
        jobs = [[DummySingleJob() for _ in range(3)] for _ in range(3)]
        multi_jobs = [MultiJob(children=js) for js in jobs]
        multi_job = MultiJob(children=multi_jobs)

        # When run top level job
        multi_job.run(runner=runner).wait()

        # Then multi-job ran ok
        assert multi_job.check()
        assert multi_job.ok()
        assert multi_job.status == JobStatus.SUCCESSFUL
        assert all([mj.status == JobStatus.SUCCESSFUL for mj in multi_jobs])
        assert all([j.status == JobStatus.SUCCESSFUL for js in jobs for j in js])

    def test_run_multiple_dependent_multijobs_all_succeed(self, config):
        runner = JobRunner(parallel=True, maxjobs=5)
        config.sleepstep = 0.1

        # Given multi-job with multiple dependent multi-jobs
        jobs = [[DummySingleJob() for _ in range(3)] for _ in range(3)]
        multi_jobs = [MultiJob(children=js) for js in jobs]
        multi_job = MultiJob(children=multi_jobs)

        @add_to_instance(jobs[1][0])
        def prerun(s):
            jobs[0][0].results.wait()

        @add_to_instance(multi_jobs[1])
        def postrun(s):
            multi_jobs[0].results.wait()

        # When run top level job
        multi_job.run(runner=runner).wait()

        # Then multi-job ran ok
        assert multi_job.check()
        assert multi_job.ok()
        assert multi_job.status == JobStatus.SUCCESSFUL
        assert all([mj.status == JobStatus.SUCCESSFUL for mj in multi_jobs])
        assert all([j.status == JobStatus.SUCCESSFUL for js in jobs for j in js])

    def test_run_multiple_dependent_multijobs_one_fails(self, config):
        runner = JobRunner(parallel=True, maxjobs=5)
        config.sleepstep = 0.1

        # Given multi-job with multiple dependent multi-jobs and one which fails
        jobs = [[DummySingleJob() for _ in range(3)] for _ in range(3)]
        jobs[1][1].command = "not a cmd"
        multi_jobs = [MultiJob(children=js) for js in jobs]
        multi_job = MultiJob(children=multi_jobs)

        # When run top level job
        multi_job.run(runner=runner).wait()

        # Then overall multi-job failed
        assert not multi_job.check()
        assert not multi_job.ok()
        assert multi_job.status == JobStatus.FAILED
        assert [mj.status for mj in multi_jobs] == [
            JobStatus.SUCCESSFUL,
            JobStatus.FAILED,
            JobStatus.SUCCESSFUL,
        ]
        assert all([j.status == JobStatus.SUCCESSFUL for j in jobs[0]])
        assert all([j.status == JobStatus.SUCCESSFUL for j in jobs[2]])
        assert [j.status for j in jobs[1]] == [
            JobStatus.SUCCESSFUL,
            JobStatus.CRASHED,
            JobStatus.SUCCESSFUL,
        ]

    def test_run_multijobs_with_use_subdir_in_parallel(self, config):
        # Given parallel job runner
        runner = JobRunner(parallel=True, maxjobs=4)

        # When run multi jobs in subdir
        outer_multi_jobs = []
        with jobs_in_directory("results"):
            for mj in ["A", "B", "C"]:
                with jobs_in_directory(mj):
                    jobs = [[DummySingleJob(name=f"dummy_{i}_{j}") for i in range(3)] for j in range(3)]
                    inner_multi_jobs = [MultiJob(children=js, name=f"multi_inner_{i}") for i, js in enumerate(jobs)]
                    multi_job = MultiJob(children=inner_multi_jobs, name=f"multi_outer_{mj}")
                    multi_job.run(runner=runner)
                    outer_multi_jobs.append(multi_job)

        # Then top-level multi-jobs are located in the correct subdirectory
        workdir = config.default_jobmanager.workdir
        for mj_outer in outer_multi_jobs:
            assert mj_outer.ok()

            outer = mj_outer.name.split("_")[-1]
            assert Path(mj_outer.path) == Path(workdir, "results", outer, mj_outer.name)

            for mj_inner in mj_outer.children:
                assert mj_inner.ok()
                assert Path(mj_inner.path) == Path(
                    workdir,
                    "results",
                    outer,
                    mj_outer.name,
                    mj_inner.name,
                )

                for j in mj_inner.children:
                    assert j.ok()
                    assert Path(j.path) == Path(workdir, "results", outer, mj_outer.name, mj_inner.name, j.name)

    def test_full_name(self):
        job = DummySingleJob(name="dummy_job")
        inner_multi_job = MultiJob(children=[job], name="multi_inner_job")
        multi_job = MultiJob(children=[inner_multi_job], name="multi_outer")
        multi_job.run()

        assert multi_job.ok()
        assert multi_job._full_name() == "multi_outer"
        assert multi_job._full_name("some/rundir") == "some/rundir/multi_outer"
        assert inner_multi_job._full_name() == "multi_outer/multi_inner_job"
        assert inner_multi_job._full_name("some/rundir") == "some/rundir/multi_outer/multi_inner_job"
        assert job._full_name() == "multi_outer/multi_inner_job/dummy_job"
        assert job._full_name("some/rundir") == "some/rundir/multi_outer/multi_inner_job/dummy_job"

    def test_apply_to_children(self):
        def add_tag(j):
            j.tag = True

        # Given nested multi-jobs
        job1 = DummySingleJob(name="dummy_job")
        job2 = DummySingleJob(name="dummy_job")
        inner_multi_job = MultiJob(children=[job1, job2], name="multi_inner_job")
        multi_job = MultiJob(children=[inner_multi_job], name="multi_outer")

        # When apply tagging function to non multi-job
        MultiJob.apply_to_children(job1, add_tag)

        # Then apply is a no-op
        assert not hasattr(multi_job, "tag")
        assert not hasattr(inner_multi_job, "tag")
        assert not hasattr(job1, "tag")
        assert not hasattr(job2, "tag")

        # When apply tagging function to children
        MultiJob.apply_to_children(multi_job, add_tag)

        # Then applies to children non-recursively
        assert not hasattr(multi_job, "tag")
        assert hasattr(inner_multi_job, "tag")
        assert not hasattr(job1, "tag")
        assert not hasattr(job2, "tag")

        # When apply tagging function recursively
        MultiJob.apply_to_children(multi_job, add_tag, recursive=True)
        assert not hasattr(multi_job, "tag")
        assert hasattr(inner_multi_job, "tag")
        assert hasattr(job1, "tag")
        assert hasattr(job2, "tag")

    def test_delete_created_multijob(self, config):
        # Given multi job
        jobs = [DummySingleJob() for _ in range(3)]
        multi_job = MultiJob(children=[j for j in jobs])

        # When deleted
        multi_job.delete()

        # Then status set to deleted for parent and child jobs
        assert multi_job.status == JobStatus.DELETED
        assert all(j.status == JobStatus.DELETED for j in jobs)
        assert multi_job.path is None
        assert all(j.path is None for j in jobs)
        assert multi_job.children == []
        assert all(j.parent is None for j in jobs)
        with pytest.raises(ResultsError):
            _ = multi_job.results.grep_output("")
        with pytest.raises(ResultsError):
            _ = jobs[0].results.grep_output("")

    def test_delete_running_multijob(self, config):
        # Given multi job
        jobs = [DummySingleJob() for _ in range(3)]
        multi_job = MultiJob(children=[j for j in jobs])
        multi_job.run()
        path = multi_job.path

        # When deleted while running
        multi_job.delete()

        # Then waits, job files removed and job path removed for parent and child jobs
        assert multi_job.status == JobStatus.DELETED
        assert all(j.status == JobStatus.DELETED for j in jobs)
        assert multi_job.name not in config.default_jobmanager.names
        assert multi_job not in config.default_jobmanager.jobs
        assert not Path(path).exists()
        assert multi_job.path is None
        assert all(j.path is None for j in jobs)
        assert multi_job.children == []
        assert all(j.parent is None for j in jobs)
        with pytest.raises(ResultsError):
            _ = multi_job.results.grep_output("")
        with pytest.raises(ResultsError):
            _ = jobs[0].results.grep_output("")

    def test_delete_nested_multijob(self, config):
        # Given multi job
        jobs = [DummySingleJob() for _ in range(3)]
        multi_job = MultiJob(children=[j for j in jobs])
        top_multi_job = MultiJob(children=[multi_job])
        top_multi_job.run()
        path = top_multi_job.path

        # When deleted
        top_multi_job.delete()

        # Then waits, job files removed and job path removed for parent and child jobs
        assert top_multi_job.status == JobStatus.DELETED
        assert multi_job.status == JobStatus.DELETED
        assert all(j.status == JobStatus.DELETED for j in jobs)
        assert top_multi_job.name not in config.default_jobmanager.names
        assert multi_job.name not in config.default_jobmanager.names
        assert top_multi_job not in config.default_jobmanager.jobs
        assert multi_job not in config.default_jobmanager.jobs
        assert not Path(path).exists()
        assert top_multi_job.path is None
        assert multi_job.path is None
        assert all(j.path is None for j in jobs)
        assert top_multi_job.children == []
        assert multi_job.children == []
        assert all(j.parent is None for j in jobs)
        with pytest.raises(ResultsError):
            _ = top_multi_job.results.grep_output("")
        with pytest.raises(ResultsError):
            _ = multi_job.results.grep_output("")
        with pytest.raises(ResultsError):
            _ = jobs[0].results.grep_output("")

    def test_rename_created_multijob(self, config):
        # Given multi job
        id = str(uuid.uuid4())[:8]
        name1 = f"to-be-renamed-{id}"
        name2 = f"renamed-{id}"
        jobs = [DummySingleJob(name=name1) for _ in range(3)]
        multi_job = MultiJob(children=[j for j in jobs], name=name1)

        # When renamed
        multi_job.rename(name2)
        multi_job.children[0].rename(name2)

        # Then name changed
        assert multi_job.status == JobStatus.CREATED
        assert multi_job.path is None
        assert multi_job.name == name2
        assert jobs[0].status == JobStatus.CREATED
        assert jobs[0].path is None
        assert jobs[0].name == name2

    def test_rename_running_multijob(self, config):
        # Given multi job
        id = str(uuid.uuid4())[:8]
        name1 = f"to-be-renamed-{id}"
        name2 = f"renamed-{id}"
        jobs = [DummySingleJob(name=name1, wait=0.2) for _ in range(3)]
        multi_job = MultiJob(children=[j for j in jobs], name=name1)
        multi_job.run()
        path1 = multi_job.path

        # When renamed while running
        multi_job.rename(name2)
        jobs[0].rename(name2)
        path2 = multi_job.path

        # Then waits, multi job files moved and renamed, re-registered in job manager
        assert multi_job.status == JobStatus.SUCCESSFUL
        assert not Path(path1).exists()
        assert Path(path2).exists()
        assert Path(path2) == Path(path1).parent / name2
        assert name1 not in config.default_jobmanager.names
        assert name2 in config.default_jobmanager.names

        # And child jobs are also moved and renamed and re-registered in job manager
        for i, job in enumerate(jobs):
            assert Path(job.path).exists()
            assert Path(job.path) == Path(path2) / job.name
            assert f"{name1}/{name1}" not in config.default_jobmanager.names
            if i == 0:
                assert f"{name2}/{name2}" in config.default_jobmanager.names
            assert f"{name2}/{name1}" in config.default_jobmanager.names

            # And results can be read
            assert job.results.read_file("$JN.in")
            assert job.results.read_file("$JN.out")
            assert job.results.read_file("$JN.run")

    def test_rename_nested_multijob(self, config):
        # Given multi job
        id = str(uuid.uuid4())[:8]
        name1 = f"top-to-be-renamed-{id}"
        name2 = f"middle-to-be-renamed-{id}"
        name3 = f"bottom-to-be-renamed-{id}"
        name4 = f"top-renamed-{id}"
        name5 = f"middle-renamed-{id}"
        name6 = f"bottom-renamed-{id}"
        jobs = [DummySingleJob(name=name3) for _ in range(3)]
        multi_job = MultiJob(children=[j for j in jobs], name=name2)
        top_multi_job = MultiJob(children=[multi_job], name=name1)
        top_multi_job.run()
        orig_path = top_multi_job.path

        # When renamed
        jobs[0].rename(name6)
        top_multi_job.rename(name4)
        multi_job.rename(name5)

        # Then waits, multi job files moved and renamed, re-registered in job manager
        assert not Path(orig_path).exists()
        assert Path(top_multi_job.path) == Path(orig_path).parent / name4
        assert Path(multi_job.path) == Path(orig_path).parent / name4 / name5
        assert Path(jobs[0].path) == Path(orig_path).parent / name4 / name5 / name6
        assert (Path(orig_path).parent / name4).exists()
        assert (Path(orig_path).parent / name4 / name5).exists()
        assert (Path(orig_path).parent / name4 / name5 / name6).exists()
        assert name1 not in config.default_jobmanager.names
        assert f"{name1}/{name2}" not in config.default_jobmanager.names
        assert f"{name1}/{name2}/{name3}" not in config.default_jobmanager.names
        assert name4 in config.default_jobmanager.names
        assert f"{name4}/{name5}" in config.default_jobmanager.names
        assert f"{name4}/{name5}/{name6}" in config.default_jobmanager.names

        # And results can be read
        for job in jobs:
            assert job.results.read_file("$JN.in")
            assert job.results.read_file("$JN.out")
            assert job.results.read_file("$JN.run")

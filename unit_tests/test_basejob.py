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

from scm.plams.core.settings import Settings
from scm.plams.core.basejob import SingleJob, MultiJob
from scm.plams.core.errors import PlamsError, FileError, ResultsError
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.jobmanager import JobManager
from scm.plams.core.functions import add_to_instance
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
        with patch("scm.plams.core.basejob.config", config):
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

        logger = get_logger(f"plams-{uuid.uuid4()}")

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


class TestMultiJob:
    """
    Test suite for the Multi Job.
    Not truly independent as relies upon the job runner/manager and results components.
    But this suite focuses on testing the methods on the job class itself.
    """

    def test_run_multiple_independent_single_jobs_all_succeed(self, config):
        runner = JobRunner(parallel=True, maxjobs=3)
        config.sleepstep = 0.1

        with patch("scm.plams.core.basejob.config", config):
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

        with patch("scm.plams.core.basejob.config", config):
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

        with patch("scm.plams.core.basejob.config", config):
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

        with patch("scm.plams.core.basejob.config", config):
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

    def test_run_multiple_independent_multijobs_all_succeed(self, config):
        runner = JobRunner(parallel=True, maxjobs=3)
        config.sleepstep = 0.1

        with patch("scm.plams.core.basejob.config", config):
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

        with patch("scm.plams.core.basejob.config", config):
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

        with patch("scm.plams.core.basejob.config", config):
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

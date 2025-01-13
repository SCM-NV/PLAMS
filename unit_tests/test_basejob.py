import pytest
import uuid
from unittest.mock import patch
from datetime import datetime
from collections import namedtuple
import shutil
import re
from io import StringIO

from scm.plams.core.settings import Settings
from scm.plams.core.basejob import SingleJob
from scm.plams.core.errors import PlamsError, FileError
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.jobmanager import JobManager
from scm.plams.core.functions import add_to_instance
from scm.plams.core.enums import JobStatus


def log_call(method):
    """
    Decorator to log calls to instance.
    """

    def wrapper(self, *args, **kwargs):
        log_entry = namedtuple("LogEntry", ["method", "args", "kwargs", "timestamp"])
        self.call_log.append(log_entry(method=method.__name__, args=args, kwargs=kwargs, timestamp=datetime.now()))
        return method(self, *args, **kwargs)

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
        self.call_log = []
        self.id = uuid.uuid4()
        self.input = inp.replace("%ID%", str(self.id)) if inp is not None else f"Dummy input {self.id}"
        self.command = cmd if cmd is not None else "sed 's/input/output/g'"
        self.wait = wait

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

    def check(self) -> bool:
        return self.results.read_file(self._filename("err")) == ""

    def get_errormsg(self) -> str:
        return self.results.read_file(self._filename("err"))


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
        assert [c.method for c in job.call_log] == ["prerun", "get_input", "get_input", "get_runscript", "postrun"]

    def test_run_writes_output_file_on_success(self):
        # When run a successful job
        job = DummySingleJob()
        results = job.run()

        # Then job check passes and output file written
        assert job.check()
        assert results.read_file("$JN.in") == f"Dummy input {job.id}"
        assert results.read_file("$JN.out") == f"Dummy output {job.id}"
        assert results.read_file("$JN.err") == ""

    def test_run_writes_error_file_on_failure(self):
        # When run an erroring job
        job = DummySingleJob(cmd="not_a_cmd")
        results = job.run()

        # Then job check fails and error file written
        assert not job.check()
        assert results.read_file("$JN.in") == f"Dummy input {job.id}"
        assert results.read_file("$JN.out") == ""
        assert results.read_file("$JN.err") != ""

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
        job1 = DummySingleJob(wait=0.2)
        job2 = DummySingleJob(wait=0.01)
        results = [job1.run(runner), job2.run(runner)]

        # Then both run in parallel
        # Shorter job finishes first even though started second
        # Pre-run call of second job is made before post-run call of first job
        results[1].wait()
        assert job2.status == JobStatus.SUCCESSFUL
        assert job1.status == JobStatus.RUNNING

        results[0].wait()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.call_log[0].timestamp < job1.call_log[-1].timestamp

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
        assert job2.call_log[0].timestamp >= job1.call_log[-1].timestamp

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
        # Get runscript call of second job is made after post-run call of first job
        results[0].wait()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.status in [JobStatus.REGISTERED, JobStatus.STARTED, JobStatus.RUNNING]

        results[1].wait()
        assert job2.status == JobStatus.SUCCESSFUL
        assert job2.call_log[2].timestamp >= job1.call_log[-1].timestamp

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

        # When call load
        job2 = DummySingleJob.load(job1.path)

        # Then job loaded successfully
        assert job1.name == job2.name
        assert job1.id == job2.id
        assert job1.path == job2.path
        assert job1.settings == job2.settings
        assert job1._filenames == job2._filenames

    def test_job_summaries_logged(self, config):
        job1 = DummySingleJob()
        job2 = DummySingleJob(inp=job1.input)
        job3 = DummySingleJob(cmd="not_a_cmd")

        job_manager = JobManager(config.jobmanager, folder="test_logging")

        job1.run(jobmanager=job_manager)
        job2.run(jobmanager=job_manager)
        job3.run(jobmanager=job_manager)

        with open(job_manager.job_logger.logfile) as f:
            assert (
                f.readline()
                == """logged_at,job_base_name,job_name,job_status,job_path,job_ok,job_check,job_get_errormsg,job_parent_name,job_parent_path
"""
            )

            assert (
                re.match(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},plamsjob,plamsjob,successful", f.readline()) is not None
            )
            assert (
                re.match(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},plamsjob,plamsjob\.002,copied", f.readline())
                is not None
            )
            assert (
                re.match(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},plamsjob,plamsjob\.003,crashed", f.readline())
                is not None
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

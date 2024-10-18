import pytest
from unittest.mock import patch, MagicMock
import time
import threading
import os
from pathlib import Path
import inspect
from io import StringIO
import re
import uuid

from scm.plams.core.functions import (
    _init,
    init,
    _finish,
    finish,
    requires_optional_package,
    add_to_class,
    add_to_instance,
    log,
)
from scm.plams.core.settings import Settings
from scm.plams.core.errors import MissingOptionalPackageError
from scm.plams.unit_tests.test_helpers import assert_config_as_expected, get_mock_open_function, temp_file_path


class TestInitAndFinish:
    """
    Test suite for the init() and finish() public function
    """

    @pytest.fixture(autouse=True)
    def mock_jobmanager(self):
        """
        Mock out the job manager to avoid creating run directories etc.
        """
        with patch("scm.plams.core.jobmanager.JobManager") as mock_jobmanager:
            yield mock_jobmanager

    @pytest.fixture(autouse=True)
    def mock_file(self):
        """
        Mock out file for logging.
        """
        with get_mock_open_function(lambda f: isinstance(f, MagicMock), None):
            yield

    @pytest.mark.parametrize("explicit_init", [False, True])
    def test_init_sets_config_settings_default_values(self, config, explicit_init):
        # When call _init and optionally init
        if explicit_init:
            init()

        # Then config settings set to default values
        assert_config_as_expected(config, explicit_init=explicit_init)

    def test_init_updates_config_settings_from_arg(self, config):
        settings = Settings()
        settings.preview = True
        settings.job.runscript.stdout_redirect = True
        settings.log.stdout = 5

        # When call init with additional settings
        init(config_settings=settings)

        # Then config settings set to default values with overrides
        assert_config_as_expected(config, preview=True, stdout_redirect=True, stdout=5)

    @pytest.mark.parametrize(
        "plams_defaults,ams_home,file1_exists,file2_exists,file3_exists,explicit_init,expected_source",
        [
            ["my/defaults/file", None, True, True, True, True, 1],
            [None, "my/home", True, True, True, True, 2],
            [None, None, True, True, True, True, 3],
            ["my/defaults/file", "my/home", True, True, True, True, 1],
            ["my/defaults/file", "my/home", False, True, True, True, 2],
            ["my/defaults/file", "my/home", False, False, True, True, 3],
            ["my/defaults/file", None, True, True, True, False, 1],
        ],
        ids=[
            "plams_defaults_only_set",
            "amshome_only_set",
            "no_env_vars_but_default_exists",
            "both_env_var_set_and_files_exist_uses_first",
            "both_env_var_set_but_only_second_file_exists_uses_second",
            "both_env_var_set_but_neither_file_exists_uses_default",
            "plams_defaults_only_set_with_implicit_init",
        ],
    )
    def test_init_updates_config_settings_from_defaults_file(
        self,
        config,
        plams_defaults,
        ams_home,
        file1_exists,
        file2_exists,
        file3_exists,
        explicit_init,
        expected_source,
        monkeypatch,
    ):
        # Patch env vars for PLAMSDEFAULTS and AMSHOME
        if plams_defaults is not None:
            plams_defaults = Path(plams_defaults)
            monkeypatch.setenv("PLAMSDEFAULTS", str(plams_defaults))
        else:
            monkeypatch.delenv("PLAMSDEFAULTS", raising=False)

        if ams_home is not None:
            ams_home = Path(ams_home)
            monkeypatch.setenv("AMSHOME", str(ams_home))
        else:
            monkeypatch.delenv("AMSHOME", raising=False)

        # Set up the content for the defaults file
        mock_open = get_mock_open_function(
            lambda p: Path(p) == plams_defaults or Path(p).name == "plams_defaults",
            content="""
config.preview = True
config.job.runscript.stdout_redirect = True
config.log.stdout = 5
        """,
        )

        # Set up portable paths for the expected location of the defaults file
        path1 = plams_defaults
        path2 = Path(f"{ams_home}/scripting/scm/plams/plams_defaults") if ams_home is not None else None
        path3 = Path(inspect.getfile(init)).parent.parent / "plams_defaults"

        # Set up which files are expected to exist
        def file_exists(path):
            p = Path(path)
            return (
                (path1 is not None and p == path1 and file1_exists)
                or (path2 is not None and p == path2 and file2_exists)
                or (p == path3 and file3_exists)
            )

        with mock_open, patch("scm.plams.core.functions.isfile") as mock_isfile:
            mock_isfile.side_effect = file_exists

            # When call _init and optionally init
            _init()  # need to re-call after env vars and file patched
            if explicit_init:
                init()

            # Then config settings set to default values with overrides
            assert_config_as_expected(config, explicit_init=explicit_init, preview=True, stdout_redirect=True, stdout=5)

            # Then the source used is the expected file
            source = Path(mock_open.new.mock_calls[0].args[0])
            if expected_source == 1:
                assert source == path1
            elif expected_source == 2:
                assert source == path2
            else:
                assert source == path3

    @pytest.mark.parametrize("explicit_init", [False, True])
    def test_init_can_set_default_jobrunner_in_defaults_file(self, explicit_init, config):
        content = r"""
from .settings import Settings
from .jobrunner import GridRunner
import re

def __psubmit_get_jobid(output):
    match = re.search(r'> Job ID:\s*([0-9]+)', output)
    if match is not None:
        return match[1]
    else:
        return None

def __pbs_running(output):
    lines = output.splitlines()[2:]
    return [line.split()[0].split('.')[0] for line in lines]

grid_config = Settings()
grid_config.workdir = '-d'
grid_config.output  = '-o'
grid_config.error   = '-e'
grid_config.commands.submit  = '/home/psubmit-plams-wrapper'
grid_config.commands.check  = 'qstat'
grid_config.commands.getid   = __psubmit_get_jobid
grid_config.commands.running = __pbs_running

config.default_jobrunner = GridRunner(grid=grid_config, sleepstep=30)
    """

        # Set up the content for the defaults file
        mock_open = get_mock_open_function(lambda p: p.endswith("plams_defaults"), content=content)

        # Set up which files are expected to exist
        with mock_open, patch("scm.plams.core.functions.isfile") as mock_isfile:
            mock_isfile.side_effect = lambda p: p.endswith("plams_defaults")

            # When call _init and optionally init
            _init()  # need to re-call after file patched
            if explicit_init:
                init()

            # Then default jobrunner is set to the gridrunner
            jobrunner = config.default_jobrunner
            from scm.plams.core.jobrunner import GridRunner

            assert isinstance(jobrunner, GridRunner)
            assert jobrunner.settings.workdir == "-d"
            assert jobrunner.settings.commands.submit == "/home/psubmit-plams-wrapper"
            assert jobrunner.settings.commands.getid("> Job ID: 1234") == "1234"

    def test_init_passes_args_to_default_jobmanager(self, config, mock_jobmanager):
        # When call init with job manager args
        init(path="foo/bar", folder="test_folder", use_existing_folder=True)

        # Then default jobmanager initialised with the passed args
        assert mock_jobmanager.call_args_list[0].args == (
            {"counter_len": 3, "hashing": "input", "remove_empty_directories": True},
            "foo/bar",
            "test_folder",
            True,
        )

    @pytest.mark.parametrize("explicit_init", [False, True])
    def test_init_updates_init_flag_on_config_settings(self, explicit_init, config):
        # When call _init and optionally init
        if explicit_init:
            init()

        # Then config is marked as initialised
        assert config.init

    def test_init_idempotent(self, config, mock_jobmanager):
        # When call init twice, the second time with updated settings
        settings = Settings()
        settings.preview = True
        settings.job.runscript.stdout_redirect = True
        settings.log.stdout = 5
        init()
        init(config_settings=settings)

        # Then the second call is a no-op
        assert_config_as_expected(config)
        assert mock_jobmanager.call_count == 1

    def test_init_initialises_again_when_config_init_flag_reset(self, config, mock_jobmanager):
        # When call init twice, the second time with updated settings
        settings = Settings()
        settings.preview = True
        settings.job.runscript.stdout_redirect = True
        settings.log.stdout = 5
        init()
        config.init = False
        init(config_settings=settings)

        # Then the second call is effective
        assert_config_as_expected(config, preview=True, stdout_redirect=True, stdout=5)
        assert mock_jobmanager.call_count == 2

    @pytest.mark.parametrize(
        "code,version,tasks,explicit_init,should_error",
        [
            (0, "20.10", "16", True, False),
            (0, "21.10", "16", True, False),
            (0, "21.10", "16", False, False),
            (1, "20.11", "16", True, True),
            (1, "20.11", "16", False, True),
            (0, "10.11", "16", True, True),
            (0, "21_11", "16", True, True),
            (0, "21.10", None, True, True),
            (0, "20.10", "xyz", True, True),
        ],
        ids=[
            "happy_v20.10",
            "happy_v21.10",
            "happy_v21.10_with_implicit_init",
            "unhappy_code",
            "unhappy_code_with_implicit_init",
            "unhappy_v10.11",
            "unhappy_v20_11",
            "unhappy_no_tasks",
            "unhappy_tasks",
        ],
    )
    def test_init_sets_slurm_settings(self, code, version, tasks, explicit_init, should_error, monkeypatch, config):
        # When running under slurm and init() called
        monkeypatch.setenv("SLURM_JOB_ID", "123456")
        if tasks is not None:
            monkeypatch.setenv("SLURM_TASKS_PER_NODE", tasks)
        else:
            monkeypatch.delenv("SLURM_TASKS_PER_NODE", raising=False)
        monkeypatch.setenv("SCM_SRUN_OPTIONS", "")

        mock_result = MagicMock()
        mock_stdout = MagicMock()
        mock_stdout.decode = lambda: f"x {version}"
        mock_result.returncode = code
        mock_result.stdout = mock_stdout

        with patch("scm.plams.core.functions.subprocess.run", return_value=mock_result):
            # When call _init and optionally init
            _init()  # need to re-call after env vars patched
            if explicit_init:
                init()

            # Then config should have defaults set
            assert_config_as_expected(config, explicit_init=explicit_init)

            # And slurm settings set when happy, otherwise None
            if not should_error:
                assert config.slurm == Settings(
                    {
                        "slurm_version": version.split("."),
                        "tasks_per_node": [int(tasks)],
                    }
                )
                assert os.environ["SCM_SRUN_OPTIONS"].startswith("-m block:block:block,NoPack --use-min-nodes")
            else:
                assert config.slurm is None

    def test_finish_awaits_plams_threads(self, config, mock_jobmanager):
        # When set flag on settings object on a different thread
        config.init = True
        config.default_jobmanager = mock_jobmanager

        def set_flag():
            time.sleep(0.1)
            config.flag = True

        thread = threading.Thread(target=set_flag, name="plamsthread")
        thread.start()

        # Then thread awaited and flag therefore set
        finish()
        assert not thread.is_alive()
        assert config.flag

    def test_finish_calls_clean_on_job_managers(self, config):
        # When finish called with job manager arguments
        config.init = True
        jobmanagers = [MagicMock() for _ in range(3)]
        config.default_jobmanager = jobmanagers[0]
        finish(jobmanagers[1:])

        # Then clean called on both default and passed job managers
        for mock_jobmanager in jobmanagers:
            assert mock_jobmanager._clean.call_count == 1

    @pytest.mark.parametrize("erase_workdir", [True, False])
    def test_finish_respects_default_job_manager_erase_workdir_flag(self, erase_workdir, config, mock_jobmanager):
        with patch("shutil.rmtree") as mock_rmtree:
            # When erase_workdir configured
            config.init = True
            config.default_jobmanager = mock_jobmanager
            config.erase_workdir = erase_workdir
            finish()
            # Then rmtree only called if flag enabled
            assert mock_rmtree.call_count == (1 if erase_workdir else 0)

    def test_finish_updates_init_flag_on_config(self, config, mock_jobmanager):
        # When call finish
        config.init = True
        config.default_jobmanager = mock_jobmanager
        finish()

        # Then init flag always set to False
        assert not config.init

    def test_finish_idempotent(self, config, mock_jobmanager):
        # When call finish twice
        config.init = True
        config.default_jobmanager = mock_jobmanager
        finish()
        finish()

        # Then second call a no-op
        assert config.default_jobmanager._clean.call_count == 1

    def test_finish_cleans_again_when_config_flag_reset(self, config, mock_jobmanager):
        # When call finish twice
        config.init = True
        config.default_jobmanager = mock_jobmanager
        finish()
        config.init = True
        finish()

        # Then second call effective
        assert config.default_jobmanager._clean.call_count == 2

    @pytest.mark.parametrize("explicit_finish", [True, False])
    def test_finish_does_not_clean_uninitialised_default_jobmanager(self, explicit_finish, config):
        # When call _finish or finish but without initialising the lazy default jobmanager
        if explicit_finish:
            finish()
        else:
            _finish()

        # Then the default job manager is still uninitialised and so not cleaned
        assert config["default_jobmanager"] is None

    @pytest.mark.parametrize("explicit_init", [True, False])
    def test_init_then_finish_as_expected(self, explicit_init, config, mock_jobmanager):
        # When call _init and optionally init then finish
        if explicit_init:
            init()
        finish()

        # Then config defaults initialised and job manager created and cleaned (on explicit init which forces job manager materialisation)
        assert_config_as_expected(config, init=False, explicit_init=explicit_init)
        assert config.default_jobmanager._clean.call_count == (1 if explicit_init else 0)
        assert mock_jobmanager.call_args_list[0].args[0] == {
            "counter_len": 3,
            "hashing": "input",
            "remove_empty_directories": True,
        }
        if explicit_init:
            assert mock_jobmanager.call_args_list[0].args[1:] == (None, None, False)

    @pytest.mark.parametrize("explicit_init", [True, False])
    def test_init_then_finish_successive_calls_as_expected(self, explicit_init, config, mock_jobmanager):
        # When call _init and optionally init then finish successively
        if explicit_init:
            init()
        config.preview = True
        config.log.stdout = 3
        config.job.runscript.stdout_redirect = True
        config.foo = "bar"
        finish()
        init()
        finish()

        # Then on the second call the default settings should be reset, but the custom ones remain
        assert_config_as_expected(config, init=False)
        assert config.foo == "bar"
        assert config.default_jobmanager._clean.call_count == (2 if explicit_init else 1)
        assert mock_jobmanager.call_args_list[0].args[0] == {
            "counter_len": 3,
            "hashing": "input",
            "remove_empty_directories": True,
        }
        if explicit_init:
            assert mock_jobmanager.call_args_list[0].args[1:] == (None, None, False)

    def test_init_then_finish_in_loop_as_expected(self, config, mock_jobmanager):
        # When call init then finish in loop
        for i in range(3):
            init(folder=f"folder{i}")
            finish()

            # Then config defaults initialised and job manager created and cleaned
            assert_config_as_expected(config, init=False)
            assert config.default_jobmanager._clean.call_count == i + 1
            assert mock_jobmanager.call_args_list[i].args == (
                {
                    "counter_len": 3,
                    "hashing": "input",
                    "remove_empty_directories": True,
                },
                None,
                f"folder{i}",
                False,
            )


class TestDecorators:
    """
    Test suite for PLAMS decorators
    """

    class EmptyClass:
        pass

    def test_add_to_class(self):
        # Given initially empty class
        empty_class = self.EmptyClass()

        with pytest.raises(AttributeError):
            empty_class.is_added_to_class()

        # When add function to class
        @add_to_class(self.EmptyClass)
        def is_added_to_class(self):
            return True

        # Then is available on all instances
        another_empty_class = self.EmptyClass()
        assert empty_class.is_added_to_class()
        assert another_empty_class.is_added_to_class()

    def test_add_to_instance(self):
        # Given initially empty class
        empty_class = self.EmptyClass()

        with pytest.raises(AttributeError):
            empty_class.is_added_to_instance()

        # When add function to instance
        @add_to_instance(empty_class)
        def is_added_to_instance(self):
            return True

        # Then is available on that instance only
        another_empty_class = self.EmptyClass()
        assert empty_class.is_added_to_instance()
        with pytest.raises(AttributeError):
            assert another_empty_class.is_added_to_instance()

    class OptionalRequirementsClass:

        def no_requirements(self):
            return True

        @requires_optional_package("numpy")
        def requires_numpy_package(self):
            import numpy  # noqa F401

            return True

        @requires_optional_package("__this_is_definitely_not_an_available_package__")
        def requires_unavailable_package(self):
            import __this_is_definitely_not_an_available_package__  # noqa F401

            return True

    def test_requires_optional_package(self):
        # Given class which has methods, some of which require optional packages
        req_class = self.OptionalRequirementsClass()

        def maybe_calls_requires_unavailable_package(does_call):
            if does_call:
                req_class.requires_unavailable_package()
            return True

        # When call methods where package requirements are satisfied
        # Then no errors are raised
        assert req_class.no_requirements()
        assert req_class.requires_numpy_package()
        assert maybe_calls_requires_unavailable_package(False)

        # When call methods where package is missing (or method which calls said method)
        # Then error is raised
        with pytest.raises(MissingOptionalPackageError):
            req_class.requires_unavailable_package()
        with pytest.raises(MissingOptionalPackageError):
            maybe_calls_requires_unavailable_package(True)

    def test_requires_optional_package_with_add_to_class_and_instance(self):
        # Given initially empty class
        empty_class = self.EmptyClass()

        @add_to_class(self.EmptyClass)
        @requires_optional_package("__this_is_definitely_not_an_available_package__")
        def requires_unavailable_package_and_is_added_to_class(self):
            return True

        @add_to_instance(empty_class)
        @requires_optional_package("__this_is_definitely_not_an_available_package__")
        def requires_unavailable_package_and_is_added_to_instance(self):
            return True

        # When call methods that require non-existing package which were added to the class
        # Then error raised
        with pytest.raises(MissingOptionalPackageError):
            empty_class.requires_unavailable_package_and_is_added_to_class()

        with pytest.raises(MissingOptionalPackageError):
            empty_class.requires_unavailable_package_and_is_added_to_instance()


class TestLog:

    @pytest.fixture(autouse=True)
    def logger(self):
        """
        Instead of re-using the same global logger object, patch with a fresh logger instance.
        """
        from scm.plams.core.logging import get_logger

        logger = get_logger(f"plams-{uuid.uuid4()}")

        with patch("scm.plams.core.functions._logger", logger):
            yield logger

    def test_log_no_init_writes_message_to_stdout(self, config):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            config.init = False

            for i in range(10):
                log(f"log with level {i}", level=i)

            self.assert_logs(mock_stdout.getvalue(), expected_lines=4)
            assert (
                mock_stdout.getvalue()
                == """log with level 0
log with level 1
log with level 2
log with level 3
"""
            )

    def test_log_with_init_writes_message_to_stdout_and_default_jobmanger_file(self, config):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file1, temp_file_path(suffix=".log") as temp_log_file2:

                # Log to both stdout and file with date and time
                config.init = True
                config.default_jobmanager = MagicMock()
                config.default_jobmanager.logfile = temp_log_file1
                config.log.stdout = 3
                config.log.file = 5
                config.log.time = True
                config.log.date = True

                for i in range(1, 10):
                    log(f"date and time log with level {i}", level=i)

                # Log to both stdout and file without date and time
                config.log.time = False
                config.log.date = False
                for i in range(1, 10):
                    log(f"log with level {i}", level=i)

                # Log to stdout and file switching logfile location
                config.default_jobmanager = MagicMock()
                config.default_jobmanager.logfile = temp_log_file2
                config.log.file = 4
                config.log.time = True
                for i in range(1, 10):
                    log(f"time log with level {i}", level=i)

                # Log only to stdout
                config.default_jobmanager = None
                config.log.date = True
                config.log.time = False
                for i in range(1, 10):
                    log(f"date log with level {i}", level=i)

                # Stdout has a log line for each level up to and including 3
                stdout_logs = mock_stdout.getvalue()
                self.assert_logs(stdout_logs, line_end=3, date_expected=True, time_expected=True)
                self.assert_logs(stdout_logs, line_start=3, line_end=6)
                self.assert_logs(stdout_logs, line_start=6, line_end=9, time_expected=True)
                self.assert_logs(stdout_logs, line_start=9, expected_lines=3, date_expected=True)

                # File1 has a log for each level up to and including 5
                with open(temp_log_file1) as tf1:
                    file1_logs = tf1.read()
                self.assert_logs(file1_logs, line_end=5, date_expected=True, time_expected=True)
                self.assert_logs(file1_logs, line_start=5, expected_lines=5)

                # File2 has a log for each level up to and including 4
                with open(temp_log_file2) as tf2:
                    file2_logs = tf2.read()
                self.assert_logs(file2_logs, expected_lines=4, time_expected=True)

    def assert_logs(
        self, logs, line_start=0, line_end=None, expected_lines=None, date_expected=False, time_expected=False
    ):
        # Convert string text to individual lines
        lines = [l for l in logs.replace("\r\n", "\n").split("\n") if l]

        # Select the required logs section
        line_end = len(lines) if line_end is None else line_end
        lines = lines[line_start:line_end]
        expected_lines = line_end - line_start if expected_lines is None else expected_lines
        assert len(lines) == expected_lines

        # Check all log lines match the expected pattern
        if date_expected and time_expected:
            pattern = "\[\d{2}\.\d{2}\|\d{2}:\d{2}:\d{2}\] date and time log with level \d"
        elif date_expected:
            pattern = "\[\d{2}\.\d{2}] date log with level \d"
        elif time_expected:
            pattern = "\[\d{2}:\d{2}:\d{2}] time log with level \d"
        else:
            pattern = "log with level \d"

        assert all([re.fullmatch(pattern, line) is not None for line in lines])

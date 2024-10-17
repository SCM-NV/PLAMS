import builtins
from unittest.mock import patch, mock_open
import pytest
import os
from importlib.util import find_spec
from pathlib import Path
import uuid
from contextlib import contextmanager

from scm.plams.core.settings import (
    SafeRunSettings,
    LogSettings,
    RunScriptSettings,
    JobSettings,
    JobManagerSettings,
    ConfigSettings,
)

real_import = builtins.__import__
real_open = builtins.open


def get_mock_import_function(import_error_name):
    """
    Gets a mock version of the import function, which will raise an ImportError for the specified module,
    but import all other modules as usual.
    """

    def import_with_error(name, globals=None, locals=None, fromlist=(), level=0):
        if name == import_error_name:
            raise ImportError(f"Mock import error for '{import_error_name}'")
        return real_import(name, globals=globals, locals=locals, fromlist=fromlist, level=level)

    return import_with_error


def get_mock_open_function(predicate, content):
    """
    Gets a patched version of the open function, which will return the supplied content for a file which matches the given predicate only,
    but use the 'real' read for all other files.
    """
    mock = mock_open(read_data=content)

    def read(file, mode="r", *args, **kwargs):
        if predicate(file):
            new_mock = mock_open(read_data=content)  # return file with reset file pointer
            return new_mock.return_value
        else:
            return real_open(file, mode, *args, **kwargs)

    mock.side_effect = read

    return patch("builtins.open", new=mock)


def get_mock_find_spec(patch_module: str, package_name: str):
    """
    Gets a patched version of the find_spec function, which will return None for the specified package,
    but use the 'real' function for all other packages
    """

    def mock_find_spec(name):
        return None if name == package_name else find_spec(name)

    return patch(f"{patch_module}.find_spec", side_effect=mock_find_spec)


@contextmanager
def temp_file_path(suffix=".txt"):
    """
    Get the path to a temporary file in the current directory, which will be deleted after the context.
    An alternative to tempfile, as the file is not opened on entry into the context.
    """
    temp_path = str(Path(__file__).parent.absolute() / f"{uuid.uuid4()}{suffix}")
    try:
        yield str(temp_path)
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


def assert_config_as_expected(
    config, init=True, explicit_init=True, preview=False, stdout_redirect=False, stdout=3, verify_derived_types=True
):
    """
    Assert the passed config object has the expected values.
    These are the default values, unless overridden via the args.
    """
    assert config.init == init
    assert config._explicit_init == explicit_init
    assert config.preview == preview
    assert config.sleepstep == 5
    assert config.ignore_failure
    assert config.daemon_threads
    assert not config.erase_workdir
    assert config.jobmanager.counter_len == 3
    assert config.jobmanager.hashing == "input"
    assert config.jobmanager.remove_empty_directories
    assert config.job.pickle
    assert config.job.pickle_protocol == -1
    assert config.job.keep == "all"
    assert config.job.save == "all"
    assert config.job.runscript.shebang == "#!/bin/sh"
    assert config.job.runscript.stdout_redirect == stdout_redirect
    assert config.job.link_files
    assert config.log.file == 5
    assert config.log.stdout == stdout
    assert config.log.time
    assert config.log.date
    assert config.saferun.repeat == 10
    assert config.saferun.delay == 1

    if verify_derived_types:
        assert isinstance(config.saferun, SafeRunSettings)
        assert isinstance(config.log, LogSettings)
        assert isinstance(config.job, JobSettings)
        assert isinstance(config.job.runscript, RunScriptSettings)
        assert isinstance(config.jobmanager, JobManagerSettings)
        assert isinstance(config, ConfigSettings)


def skip_if_no_ams_installation():
    """
    Check whether the AMSBIN environment variable is set, and therefore if there is an AMS installation present.
    If there is no installation, skip the test with a warning.
    """
    if os.getenv("AMSBIN") is None:
        pytest.skip("Skipping test as cannot find AMS installation. '$AMSBIN' environment variable is not set.")


def skip_if_no_scm_pisa():
    """
    Check whether SCM PISA is available, and skip the test with a warning if it is not available.
    """
    try:
        import scm.pisa  # noqa F401
    except ImportError:
        pytest.skip("Skipping test as cannot find scm.pisa package.")


def skip_if_no_scm_libbase():
    """
    Check whether SCM libbase is available, and skip the test with a warning if it is not available.
    """
    try:
        import scm.libbase  # noqa F401
    except ImportError:
        pytest.skip("Skipping test as cannot find scm.libbase package.")

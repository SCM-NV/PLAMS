import builtins
import pytest
import os

real_import = builtins.__import__


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


def skip_if_no_ams_installation():
    """
    Check whether the AMSBIN environment variable is set, and therefore if there is an AMS installation present.
    If there is no installation, skip the test with a warning.
    """
    if os.getenv("AMSBIN") is None:
        pytest.skip("Skipping test as cannot find AMS installation. '$AMSBIN' environment variable is not set.")

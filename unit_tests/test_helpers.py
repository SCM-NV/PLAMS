import builtins

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

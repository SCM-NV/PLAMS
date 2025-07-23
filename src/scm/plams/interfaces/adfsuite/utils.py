from typing import Optional
import functools
import os
import subprocess
import re

from scm.plams.interfaces.adfsuite.errors import AMSBINEnvVarNotSetError, AMSExecutionError, AMSVersionError, AMSError


def requires_ams(minimum_version: Optional[str] = None):
    """
    Ensures the AMS executable is available, with the required version constraints, before running a function.
    Otherwise, raises an :class:`~scm.plams.interfaces.adfsuite.errors.AMSError`.

    :param minimum_version: minimum version of AMS required to run the function.
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            amsbin = os.environ.get("AMSBIN")
            if not amsbin:
                raise AMSBINEnvVarNotSetError()

            ams = os.path.join(amsbin, "ams")
            try:
                # N.B. return code is 1 for version check
                result = subprocess.run([ams, "--version"], capture_output=True, text=True)

                if result.stderr:
                    raise AMSExecutionError(command="$AMSBIN/ams --version", error=result.stderr)

                match = re.search(r"release=(\S+)", result.stdout)
                version = match.group(1) if match else None
            except Exception as e:
                raise AMSExecutionError(command="$AMSBIN/ams --version", error=e)

            if version and (not minimum_version or version >= minimum_version):
                return func(*args, **kwargs)
            else:
                raise AMSVersionError(version=version, minimum_version=minimum_version)

        return wrapper

    return decorator

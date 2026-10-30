import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


# Keep Matplotlib's font cache private to this test process.
# A stale/broken cache must not affect later test runs.
_original_mplconfigdir = os.environ.get("MPLCONFIGDIR")
_matplotlib_config_dir = tempfile.TemporaryDirectory(prefix="plams-matplotlib-")
os.environ["MPLCONFIGDIR"] = _matplotlib_config_dir.name


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    _matplotlib_config_dir.cleanup()
    if _original_mplconfigdir is None:
        os.environ.pop("MPLCONFIGDIR", None)
    else:
        os.environ["MPLCONFIGDIR"] = _original_mplconfigdir


UNIT_TESTS_DIR = Path(__file__).resolve().parent
if str(UNIT_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(UNIT_TESTS_DIR))


@pytest.fixture(autouse=True)
def config():
    """
    Instead of re-using the same global config object, patch with a fresh config settings instance.
    Wrap in implicitly called _init and _finish functions, as occurs in PLAMS.
    """
    from scm.plams.core.functions import _init, _finish, config
    from scm.plams.core.settings import ConfigSettings

    # Set the workdir created by the "real" config (triggered by the loading of the module) to be deleted
    config.erase_workdir = True

    config_settings = ConfigSettings()

    with patch("scm.plams.core.functions.config", config_settings):
        _init()

        yield config_settings

        # Erase any scratch test directory, if the job manager is not mocked
        config_settings.erase_workdir = config_settings["default_jobmanager"] is not None and not isinstance(
            config_settings.default_jobmanager, MagicMock
        )

        _finish()


@pytest.fixture(scope="session")
def xyz_folder():
    """
    Returns the path to the XYZ folder
    """
    p = Path(__file__).parent.absolute() / "xyz"
    assert p.exists()
    return p


@pytest.fixture(scope="session")
def inp_folder() -> Path:
    """Directory containing standalone AMS input fixtures."""
    path = Path(__file__).resolve().parent / "inp"
    assert path.is_dir()
    return path


@pytest.fixture
def pdb_folder():
    """
    Returns the path to the PDB folder
    """
    p = Path(__file__).parent.absolute() / "pdb"
    assert p.exists()
    return p


@pytest.fixture
def rkf_folder():
    """
    Returns the path to the RKF folder
    """
    p = Path(__file__).parent.absolute() / "rkf"
    assert p.exists()
    return p


@pytest.fixture
def coskf_folder():
    """
    Returns the path to the COSKF folder
    """
    p = Path(__file__).parent.absolute() / "coskf"
    assert p.exists()
    return p


@pytest.fixture(scope="session")
def dill_folder():
    """
    Returns the path to the dill fixture folder.
    """
    p = Path(__file__).parent.absolute() / "dill"
    assert p.exists()
    return p


@pytest.fixture
def vasp_folder():
    "Return the path to VASP"
    p = Path(__file__).parent.absolute() / "VASP"
    assert p.exists()
    return p

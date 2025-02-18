import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path


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


@pytest.fixture
def xyz_folder():
    """
    Returns the path to the XYZ folder
    """
    p = Path(__file__).parent.absolute() / "xyz"
    assert p.exists()
    return p


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

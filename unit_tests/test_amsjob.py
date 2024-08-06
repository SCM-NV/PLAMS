import dill as pickle
import pytest

from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.core.settings import Settings


def test_pickle():
    job = AMSJob()
    pickle_bytes = pickle.dumps(job)
    job2 = pickle.loads(pickle_bytes)
    assert isinstance(job2, AMSJob)


def test_pickle_settings():
    s = Settings()
    s.input.ams.Task = "GeometryOptimization"
    s.input.ams.Properties.NormalModes = "Yes"
    s.input.DFTB.Model = "GFN1-xTB"
    job = AMSJob(settings=s)
    pickle_bytes = pickle.dumps(job)
    job2 = pickle.loads(pickle_bytes)
    assert isinstance(job2, AMSJob)
    assert job2.settings == job.settings


def test_pickle_pisa():
    try:
        from scm.input_classes.drivers import AMS
        from scm.input_classes.engines import DFTB

    except ImportError:
        pytest.skip("Skipping test because optional 'scm.pisa' package is not available")
    driver = AMS()
    driver.Task = "GeometryOptimization"
    driver.Properties.NormalModes = True
    driver.Engine = DFTB()
    driver.Engine.Model = "GFN1-xTB"
    job = AMSJob(settings=driver)
    pickle_bytes = pickle.dumps(job)
    job2 = pickle.loads(pickle_bytes)
    assert isinstance(job2, AMSJob)
    assert job2.settings == job.settings

import dill as pickle
import pytest
from unittest.mock import MagicMock

from scm.plams.interfaces.adfsuite.ams import AMSJob, AMSResults
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


@pytest.mark.parametrize(
    "status,expected",
    [
        ["NORMAL TERMINATION", True],
        ["NORMAL TERMINATION with warnings", True],
        ["NORMAL TERMINATION with errors", False],
        ["Input error", False],
        [None, False],
    ],
)
def test_check_returns_true_for_normal_termination_with_no_errors_otherwise_false(status, expected):
    # Given job with results of certain status
    job = AMSJob()
    job.results = MagicMock(spec=AMSResults)
    job.results.readrkf.return_value = status

    # When check the job
    # Then job check is ok only for normal termination with no errors
    assert job.check() == expected

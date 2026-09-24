from copy import deepcopy
from typing import Any, Dict

import pytest

from scm.plams import CRSJob, Settings
from scm.plams.interfaces.adfsuite import crs_method_parameters as helpers


@pytest.fixture
def presets(monkeypatch: pytest.MonkeyPatch) -> Dict[str, Any]:
    data = {
        "adf-combi2005": {
            "method": "COSMO-RS",
            "parameter_blocks": {
                "CRSParameters": {"lambda": 0.13, "chb": 8850.0},
                "Dispersion": {"H": -0.034},
            },
        },
        "2013-adf-xiong": {
            "method": "COSMOSAC2013",
            "parameter_blocks": {"SACParameters": {"aeff": 6.4813}, "Epsilon": {"C.sp3": 29160.92}},
        },
    }
    monkeypatch.setattr(helpers, "get_crs_method_parameters_metadata", lambda: data)
    return data


@pytest.mark.parametrize("typed", [False, True], ids=["settings", "inputs"])
def test_apply_and_switch_presets(presets: Dict[str, Any], typed: bool) -> None:
    original_metadata = deepcopy(presets)
    if typed:
        CRS = pytest.importorskip("scm.inputs").CRS
        target = CRS(PROPERTY=CRS.PROPERTYBlock(header="LLE"), COMPOUND=[CRS.COMPOUNDBlock(header="water.coskf")])
        target.TECHNICAL.LLE.eps_g = 1.0e-5
        apply = CRSJob.apply_parameter_set_to_inputs
    else:
        target = Settings()
        target.input.property._h = "LLE"
        target.input.compound = [Settings({"_h": "water.coskf"})]
        target.input.technical.lle.eps_g = 1.0e-5
        target.runscript.nproc = 2
        apply = CRSJob.apply_parameter_set_to_settings

    assert apply(target, " adf-combi2005 ") is target
    text = CRSJob(settings=target).get_input().lower()
    assert "lambda 0.13" in text
    assert "h -0.034" in text
    if typed:
        target.CRSPARAMETERS.chb = 9000.0
    else:
        target.input.crsparameters.chb = 9000.0
    assert presets == original_metadata

    assert apply(target, "2013-adf-xiong") is target
    text = CRSJob(settings=target).get_input().lower()
    assert "method cosmosac2013" in text
    assert "sacparameters" in text
    assert "aeff 6.4813" in text
    assert "epsilon" in text
    assert "c.sp3 29160.92" in text
    assert "crsparameters" not in text
    assert "dispersion" not in text
    assert "property lle" in text
    assert "compound water.coskf" in text
    assert "eps_g" in text
    if typed:
        assert target.TECHNICAL.LLE.eps_g == pytest.approx(1.0e-5)
        assert "CRSPARAMETERS" not in target.model_fields_set
    else:
        assert target.input.technical.lle.eps_g == pytest.approx(1.0e-5)
        assert target.runscript.nproc == 2

    with pytest.raises(ValueError, match="Unknown parameter_set"):
        apply(target, "missing")
    assert CRSJob(settings=target).get_input().lower() == text
    apply(target, "adf-combi2005")
    text = CRSJob(settings=target).get_input().lower()
    assert "sacparameters" not in text
    assert "epsilon" not in text


def test_invalid_typed_parameters_leave_model_unchanged(presets: Dict[str, Any]) -> None:
    CRS = pytest.importorskip("scm.inputs").CRS
    target = CRS(METHOD="COSMOSAC2013")
    before = target.to_input()
    presets["adf-combi2005"]["parameter_blocks"]["CRSParameters"]["chb"] = "invalid"
    with pytest.raises(ValueError):
        CRSJob.apply_parameter_set_to_inputs(target, "adf-combi2005")
    assert target.to_input() == before


def test_parameter_set_options(presets: Dict[str, Any]) -> None:
    presets["adf-combi2005-alternative"] = {"method": "COSMO-RS", "parameter_blocks": {}}
    assert CRSJob.get_parameter_set_options() == {
        "COSMO-RS": ("adf-combi2005", "adf-combi2005-alternative"),
        "COSMOSAC2013": ("2013-adf-xiong",),
    }

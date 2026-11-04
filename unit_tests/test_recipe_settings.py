import subprocess
import sys
from typing import Any

import pytest

from scm.plams import AMSJob, Settings
from scm.plams.recipes._utils import to_plams_settings


@pytest.fixture
def typed_inputs() -> Any:
    return pytest.importorskip("scm.inputs")


def test_import_does_not_load_generated_ams_model(typed_inputs) -> None:
    subprocess.run(
        [sys.executable, "-c", "import scm.plams, scm.inputs; assert 'AMS' not in vars(scm.inputs)"],
        check=True,
        capture_output=True,
        text=True,
    )


def test_legacy_recipes_import_without_optional_scm_packages() -> None:
    code = """
import sys
sys.modules['scm.base'] = None
sys.modules['scm.inputs'] = None
from scm.plams import AMSNVEJob, Settings
settings = Settings()
settings.input.ADF.Basis.Type = 'DZ'
assert AMSNVEJob(settings=settings).settings.input.ams.Task == 'MolecularDynamics'
"""
    subprocess.run([sys.executable, "-c", code], check=True, capture_output=True, text=True)


@pytest.mark.parametrize("wrapped", [False, True])
def test_settings_conversion_preserves_options_and_does_not_mutate_input(typed_inputs, wrapped: bool) -> None:
    model = typed_inputs.AMS(Task="SinglePoint", Engine=typed_inputs.ADF())
    settings = Settings()
    settings.input = model
    settings.runscript.nproc = 2
    settings.adfnbo = ["write", "fock"]
    source = AMSJob(settings=settings) if wrapped else settings

    converted = to_plams_settings(source, "test")
    converted.input.ADF.Basis.Type = "TZP"
    converted.runscript.nproc = 4

    assert converted.adfnbo == ["write", "fock"]
    assert converted.input.ams.Task == "SinglePoint"
    original_settings = source.settings if wrapped else settings
    assert original_settings.runscript.nproc == 2
    assert not model.Engine.is_present("Basis")


def test_legacy_settings_work_without_optional_inputs(monkeypatch) -> None:
    monkeypatch.setattr("scm.plams.recipes._utils._has_scm_inputs", False)
    settings = Settings()
    settings.input.ADF.Basis.Type = "DZ"

    assert to_plams_settings(settings, "test") is settings


def test_legacy_job_settings_are_unwrapped() -> None:
    settings = Settings()
    settings.input.ADF.Basis.Type = "DZ"
    source = AMSJob(settings=settings)

    assert to_plams_settings(source, "test") is source.settings


def test_bare_ams_model_is_converted(typed_inputs) -> None:
    model = typed_inputs.AMS(Task="SinglePoint", Engine=typed_inputs.ADF())

    converted = to_plams_settings(model, "test")

    assert isinstance(converted, Settings)
    assert converted.input.ams.Task == "SinglePoint"
    assert "ADF" in converted.input


def test_bare_engine_model_is_rejected(typed_inputs) -> None:
    with pytest.raises(TypeError, match="scm.inputs.AMS"):
        to_plams_settings(typed_inputs.ADF(), "test")


@pytest.mark.parametrize("names", [("",), ("reactant",), ("", "final")])
def test_system_block_conversion_preserves_headers_and_charge(typed_inputs, names) -> None:
    base = pytest.importorskip("scm.base")
    model = typed_inputs.AMS()
    for name in names:
        system = base.ChemicalSystem()
        system.charge = 1
        model.Systems[name] = system

    converted = to_plams_settings(model, "test")
    systems = converted.input.ams.System

    if names == ("",):
        assert isinstance(systems, Settings)
        assert float(systems.Charge) == 1
    else:
        assert isinstance(systems, list)
        assert {system.get("_h", "") for system in systems} == set(names)
    text = AMSJob(settings=converted).get_input()
    assert text.count("Charge 1") == len(names)
    assert set(model.Systems) == set(names)

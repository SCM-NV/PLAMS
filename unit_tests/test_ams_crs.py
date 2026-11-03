from unittest.mock import MagicMock

import pytest

from scm.plams import Atom, Molecule, Settings
from scm.plams.recipes import ams_crs


@pytest.mark.parametrize("use_system", [False, True])
@pytest.mark.parametrize("typed", [False, True])
@pytest.mark.parametrize("solute_mode", ["single", "list", "none", "empty"])
def test_run_crs_ams_inputs(monkeypatch, use_system, typed, solute_mode) -> None:
    if use_system:
        structure = pytest.importorskip("scm.base").ChemicalSystem.from_smiles("O")
    else:
        structure = Molecule()
        structure.add_atom(Atom(symbol="O"))
    if typed:
        inputs = pytest.importorskip("scm.inputs")
        settings = inputs.AMS(Engine=inputs.ADF())
        settings.Engine.Basis.Type = "DZP"
        original = settings.to_input()
    else:
        settings = Settings()
        settings.input.ADF.Basis.Type = "DZP"
        original = settings.copy()
    solvent_results = MagicMock()
    solute_results = MagicMock()
    ams = MagicMock(side_effect=[solvent_results, solute_results])
    crs_results = MagicMock()
    crs_results.job.name = "CRSJob.test"
    crs = MagicMock(return_value=crs_results)
    monkeypatch.setattr(ams_crs, "run_amsjob", ams)
    monkeypatch.setattr(ams_crs, "run_crsjob", crs)
    crs_settings = Settings()
    solutes = {"single": structure, "list": [structure], "none": None, "empty": []}[solute_mode]

    result, solvents, solutes_out = ams_crs.run_crs_ams(
        settings, crs_settings, structure, solutes, return_amsresults=True, watch=True
    )

    has_solute = solute_mode in ("single", "list")
    assert result == {"CRSJob.test": crs_results}
    assert solvents == [solvent_results]
    assert solutes_out == ([solute_results] if has_solute else [None])
    assert ams.call_count == (2 if has_solute else 1)
    for call in ams.call_args_list:
        assert call.args[0] is structure
        converted = call.args[1]
        assert converted.input.ADF.Basis.Type == "DZP"
        assert converted.input.ADF.XC.GGA == "BP86"
        assert "name=CRS" in converted.input.ADF.Solvation.Solv
        assert call.kwargs == {"watch": True}
    crs.assert_called_once_with(
        solvent_results, crs_settings, solute=solute_results if has_solute else None, watch=True
    )
    assert (settings.to_input() if typed else settings) == original


@pytest.mark.parametrize("use_system", [False, True])
def test_run_amsjob_structure_and_restart(monkeypatch, use_system) -> None:
    if use_system:
        structure = pytest.importorskip("scm.base").ChemicalSystem.from_smiles("O")
        name = "AMSJob.mol"
    else:
        structure = Molecule()
        structure.properties.name = "water"
        name = "AMSJob.water"
    settings = Settings()
    settings.input.ADF.Solvation.Solv = "name=CRS"
    original = settings.copy()
    gas = MagicMock()
    gas.run.return_value.__getitem__.side_effect = lambda key: "/gas/" + key
    solv = MagicMock()
    solv.name = name
    constructor = MagicMock(side_effect=[gas, solv])
    convert = MagicMock()
    monkeypatch.setattr(ams_crs, "AMSJob", constructor)
    monkeypatch.setattr(ams_crs, "convert_to_coskf", convert)

    result = ams_crs.run_amsjob(structure, settings, watch=True)

    first, second = constructor.call_args_list
    assert first.kwargs["molecule"] is structure
    assert first.kwargs["name"] == name + ".gas"
    assert second.kwargs["molecule"] is None
    assert second.kwargs["depend"] == [gas]
    assert second.kwargs["settings"].input.ams.EngineRestart == "/gas/adf.rkf"
    assert second.kwargs["settings"].input.ams.LoadSystem.File == "/gas/ams.rkf"
    assert result is solv.run.return_value
    convert.assert_called_once_with(result, name.split(".")[1])
    gas.run.assert_called_once_with(watch=True)
    solv.run.assert_called_once_with(watch=True)
    assert settings == original

from unittest.mock import MagicMock

import pytest

from scm.plams import AMSJob, Atom, Molecule, Settings
from scm.plams.recipes.reorganization_energy import ReorganizationEnergyJob


@pytest.mark.parametrize("use_system", [False, True])
@pytest.mark.parametrize("typed", [False, True])
def test_reorganization_energy_preserves_structure_types_and_state_charges(use_system: bool, typed: bool) -> None:
    structure = Molecule()
    structure.add_atom(Atom(symbol="He"))
    structure.properties.charge = 1
    if use_system:
        base = pytest.importorskip("scm.base")
        structure = base.ChemicalSystem(AMSJob(molecule=structure).get_input())
    if typed:
        inputs = pytest.importorskip("scm.inputs")
        base = pytest.importorskip("scm.base")
        common = inputs.AMS(Engine=inputs.ADF())
        state_a = inputs.AMS()
        state_a.System = base.ChemicalSystem()
        state_a.System.charge = 1
        state_b = inputs.AMS()
    else:
        common = Settings()
        common.input.ADF.Basis.Type = "DZ"
        state_a = Settings()
        state_a.input.ams.System.Charge = 1
        state_b = Settings()
    job = ReorganizationEnergyJob(structure, common, state_a, state_b)

    assert len(job.children) == 4
    getter_name = "get_main_system" if use_system else "get_main_molecule"
    for source_name, target_name, charge in (
        ("go_A", "sp_B_for_geo_A", 0),
        ("go_B", "sp_A_for_geo_B", 1),
    ):
        source = job.children[source_name]
        target = job.children[target_name]
        getter = MagicMock(side_effect=structure.copy)
        setattr(source.results, getter_name, getter)

        target.prerun()

        assert isinstance(target.molecule, type(structure))
        assert source.settings.input.ams.Task == "GeometryOptimization"
        assert target.settings.input.ams.Task == "SinglePoint"
        assert float(target.settings.input.ams.System.Charge) == charge
        getter.assert_called_once_with()
    if typed:
        assert common.Task is None
        assert state_a.System.charge == 1
        assert not state_b.is_present("Systems")
    else:
        assert "Task" not in common.input.ams
        assert "System" not in state_b.input.ams

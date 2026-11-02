from unittest.mock import MagicMock

import pytest

from scm.plams import AMSJob, Atom, Molecule, Settings
from scm.plams.recipes.counterpoise import CounterpoiseEnergyJob


@pytest.mark.parametrize("use_system", [False, True])
@pytest.mark.parametrize("typed", [False, True])
@pytest.mark.parametrize("from_job", [False, True])
def test_counterpoise_preserves_structure_type_and_ghost_atoms(use_system: bool, typed: bool, from_job: bool) -> None:
    structure = Molecule()
    for symbol, coords in (("O", (0, 0, 0)), ("H", (0, 0, 1)), ("H", (0, 1, 0))):
        structure.add_atom(Atom(symbol=symbol, coords=coords))
    if use_system:
        base = pytest.importorskip("scm.base")
        structure = base.ChemicalSystem(AMSJob(molecule=structure).get_input())
    if typed:
        inputs = pytest.importorskip("scm.inputs")
        settings = inputs.AMS(Engine=inputs.ADF())
    else:
        settings = Settings()
        settings.input.ADF.Basis.Type = "DZ"
    source = AMSJob(molecule=structure, settings=settings) if from_job else structure
    job = CounterpoiseEnergyJob(source, ids_state_A=[1], settings_state_AB=settings)
    main = source if from_job else job.children["AB_geo_AB_basis_AB"]
    getter_name = "get_main_system" if use_system else "get_main_molecule"
    getter = MagicMock(side_effect=structure.copy)
    setattr(main.results, getter_name, getter)

    assert len(job.children) == (6 if from_job else 7)
    for name, child in job.children.items():
        if child is not main:
            child.prerun()
            assert isinstance(child.molecule, type(structure))
            assert child.settings.input.ams.System.Charge == 0
            if name.endswith("basis_AB"):
                ghosts = [atom.is_ghost if use_system else atom.properties.ghost for atom in child.molecule]
                assert ghosts == ([False, True, True] if name.startswith("A_") else [True, False, False])
            else:
                assert len(child.molecule) == (1 if name.startswith("A_") else 2)
    assert getter.call_count == 6
    if typed:
        assert not settings.Engine.is_present("Basis")

from unittest.mock import MagicMock

import pytest

from scm.plams import Atom, Molecule, Settings
from scm.plams.core.enums import JobStatus
from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob
from scm.plams.recipes._utils import get_structure_charge


@pytest.mark.parametrize("use_system", [False, True])
@pytest.mark.parametrize("typed", [False, True])
def test_compound_inputs_and_preoptimization(use_system, typed) -> None:
    if use_system:
        structure = pytest.importorskip("scm.base").ChemicalSystem.from_smiles("O")
        structure.charge = -1
    else:
        structure = Molecule()
        structure.add_atom(Atom(symbol="O"))
        structure.add_atom(Atom(symbol="H", coords=(1, 0, 0)))
        structure.add_atom(Atom(symbol="H", coords=(0, 1, 0)))
        structure.properties.charge = -1
    if typed:
        inputs = pytest.importorskip("scm.inputs")
        settings = inputs.AMS(Engine=inputs.ADF())
        original = settings.to_input()
    else:
        settings = Settings()
        settings.input.ADF.Basis.Type = "DZP"
        original = settings.copy()

    job = ADFCOSMORSCompoundJob(structure, settings=settings, preoptimization="UFF")
    info, atomic_ion = job.get_compound_properties(structure)
    assert info["Formula"] == "H2O"
    assert info["Molar Mass"] > 18
    assert not atomic_ion
    assert get_structure_charge(structure) == -1
    assert isinstance(job.settings, Settings)

    getter = "get_main_system" if use_system else "get_main_molecule"
    optimized = structure.copy()
    preoptimization = job.children["preoptimization"]
    setattr(preoptimization.results, getter, MagicMock(return_value=optimized))
    gas = job.children["gas"]
    gas.parent = job
    gas.prerun()
    assert gas.molecule is optimized

    job.children["solv"].results.get_main_system = MagicMock(return_value=optimized)
    preoptimization.results.get_input_system = MagicMock(return_value=structure)
    job.status = JobStatus.SUCCESSFUL
    assert job.results.get_main_system() is optimized
    assert job.results.get_input_system() is structure
    if typed:
        assert settings.to_input() == original
    else:
        assert settings.input == original.input

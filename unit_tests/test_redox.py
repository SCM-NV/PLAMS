from unittest.mock import MagicMock

import pytest

from scm.plams import Atom, Molecule, Settings
from scm.plams.recipes.redox import (
    AMSRedoxDirectJob,
    AMSRedoxScreeningJob,
    AMSRedoxThermodynamicCycleJob,
)


@pytest.fixture(params=[False, True], ids=["molecule", "system"])
def structure(request):
    if request.param:
        result = pytest.importorskip("scm.base").ChemicalSystem.from_smiles("O")
        result.charge = -1
    else:
        result = Molecule()
        for symbol, coords in (("O", (0, 0, 0)), ("H", (1, 0, 0)), ("H", (0, 1, 0))):
            result.add_atom(Atom(symbol=symbol, coords=coords))
        result.properties.charge = -1
    return result


@pytest.fixture(params=[False, True], ids=["settings", "typed"])
def settings(request):
    if request.param:
        inputs = pytest.importorskip("scm.inputs")
        result = inputs.AMS(Engine=inputs.ADF())
        result.Engine.Basis.Type = "DZP"
    else:
        result = Settings()
        result.input.ADF.Basis.Type = "DZP"
    return result


@pytest.mark.parametrize("workflow", ["direct", "screening", "cycle"])
def test_redox_charges_spin_and_geometry_transfers(structure, settings, workflow: str) -> None:
    original = settings.input.copy() if isinstance(settings, Settings) else settings.to_input()
    if workflow == "direct":
        job = AMSRedoxDirectJob(structure, "Water", settings=settings, oxidation=True, reduction=True)
        transfers = [("job_0", "job_ox"), ("job_0", "job_red")]
        charge_jobs = [(f"job_{state}", charge) for state, charge in (("0", -1), ("ox", 0), ("red", -2))]
    elif workflow == "screening":
        job = AMSRedoxScreeningJob(structure, "water.coskf", settings=settings, oxidation=True, reduction=True)
        transfers = [("dftb_go_0", "dftb_go_ox"), ("dftb_go_0", "dftb_go_red")]
        charge_jobs = [(f"dftb_go_{state}", charge) for state, charge in (("0", -1), ("ox", 0), ("red", -2))]
    else:
        job = AMSRedoxThermodynamicCycleJob(structure, settings=settings, oxidation=True, reduction=True)
        transfers = [("go_0_vacuum", "go_ox_vacuum"), ("go_0_vacuum", "go_red_vacuum")]
        charge_jobs = []
        for state, charge in (("0", -1), ("ox", 0), ("red", -2)):
            transfers.extend(
                [
                    (f"go_{state}_vacuum", f"go_{state}_vacuum_sp_solvated"),
                    (f"go_{state}_vacuum", f"go_{state}_solvated"),
                    (f"go_{state}_solvated", f"go_{state}_solvated_sp_vacuum"),
                ]
            )
            charge_jobs.extend(
                (f"go_{state}_{suffix}", charge)
                for suffix in ("vacuum", "vacuum_sp_solvated", "solvated", "solvated_sp_vacuum")
            )

    use_system = not isinstance(structure, Molecule)
    assert job._use_chemical_system == use_system
    assert job.orig_charge == -1
    assert job.ox_charge == 0
    assert job.red_charge == -2
    assert job.name == "plamsjob"
    assert isinstance(job.settings, Settings)
    assert job.settings.input.ADF.Basis.Type == "DZP"

    getter_name = "get_main_system" if use_system else "get_main_molecule"
    for source_name, target_name in transfers:
        optimized = structure.copy()
        getter = MagicMock(return_value=optimized)
        source = job.children[source_name]
        target = job.children[target_name]
        setattr(source.results, getter_name, getter)
        target.parent = job
        target.prerun()
        assert target.molecule is optimized
        getter.assert_called_once_with()

    for name, charge in charge_jobs:
        child = job.children[name]
        assert child.settings.input.ams.System.Charge == charge
        if workflow != "screening":
            assert child.settings.input.ADF.SpinPolarization == (10 - charge) % 2
            assert child.settings.input.ADF.Unrestricted == ("Yes" if (10 - charge) % 2 else "No")

    if workflow == "screening":
        for state, charge in (("0", -1), ("ox", 0), ("red", -2)):
            optimized = structure.copy()
            if use_system:
                optimized.charge = charge
            else:
                optimized.properties.charge = charge
            source = job.children[f"dftb_go_{state}"]
            target = job.children[f"gencoskf_{state}"]
            getter = MagicMock(return_value=optimized)
            setattr(source.results, getter_name, getter)
            target.parent = job
            target.prerun()
            assert target.input_molecule is optimized
            assert target._use_chemical_system == use_system
            assert not target.atomic_ion
            assert target.settings.input.ADF.SpinPolarization == (10 - charge) % 2
            assert target.settings.input.ams.System.Charge == charge
            gas = target.children["gas"]
            assert gas.settings.input.ams.System.Charge == charge
            gas.results.wait = MagicMock()
            setattr(gas.results, getter_name, MagicMock(return_value=optimized))
            solv = target.children["solv"]
            solv.parent = target
            solv.prerun()
            assert "System" not in solv.settings.input.ams
            assert solv.settings.input.ams.LoadSystem.File == "../gas/ams.rkf"
            getter.assert_called_once_with()

    assert (settings.input if isinstance(settings, Settings) else settings.to_input()) == original

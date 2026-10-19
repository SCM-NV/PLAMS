from pathlib import Path
from typing import Type
from unittest.mock import MagicMock, patch

import pytest

from scm.plams import AMSJob, AMSMDJob, AMSNPTJob, AMSNVEJob, AMSNVTJob, Atom, Molecule, Settings
from scm.plams.core.enums import JobStatus
from scm.plams.recipes.md.movingrestraints import AMSMovingRestraintBondSwitchJob, AMSMovingRestraintRMSDJob
from scm.plams.recipes.md.amsmdjob import AMSNPTResults
from scm.plams.recipes.md.equilibratedensity import EquilibrateDensityResults
from scm.plams.recipes.md.scandensity import AMSMDScanDensityJob


def _periodic_molecule() -> Molecule:
    molecule = Molecule()
    molecule.add_atom(Atom(symbol="H", coords=(0, 0, 0)))
    molecule.lattice = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    return molecule


@pytest.mark.parametrize("job_type", [AMSMDJob, AMSNVEJob, AMSNVTJob, AMSNPTJob])
@pytest.mark.parametrize("wrapped", [False, True])
def test_md_jobs_convert_ams_model_to_settings(job_type: Type[AMSMDJob], wrapped: bool) -> None:
    inputs = pytest.importorskip("scm.inputs")
    model = inputs.AMS()
    model.MolecularDynamics.NSteps = 10
    model.Engine = inputs.DFTB(Model="GFN1-xTB")
    settings = Settings()
    settings.input = model
    settings.runscript.nproc = 2

    job = job_type(settings=settings if wrapped else model, nsteps=20)

    assert isinstance(job.settings.input, Settings)
    assert job.settings.input.ams.Task == "MolecularDynamics"
    assert job.settings.input.ams.MolecularDynamics.NSteps == 20
    assert job.settings.input.DFTB.Model == "GFN1-xTB"
    assert model.MolecularDynamics.NSteps == 10
    if wrapped:
        assert settings.input is model
        assert job.settings.runscript.nproc == 2


@pytest.mark.parametrize("job_type", [AMSMDJob, AMSNVEJob, AMSNVTJob, AMSNPTJob])
def test_md_jobs_keep_plams_settings(job_type: Type[AMSMDJob]) -> None:
    settings = Settings()
    settings.input.ams.MolecularDynamics.NSteps = 12

    job = job_type(settings=settings, timestep=0.5)

    assert job.settings.input.ams.MolecularDynamics.NSteps == 12
    assert job.settings.input.ams.MolecularDynamics.TimeStep == 0.5


def test_md_job_converts_input_from_settings_job() -> None:
    inputs = pytest.importorskip("scm.inputs")
    model = inputs.AMS()
    model.MolecularDynamics.NSteps = 10
    source = AMSJob(settings=model)

    job = AMSNVEJob(settings=source, nsteps=20)

    assert isinstance(job.settings.input, Settings)
    assert job.settings.input.ams.MolecularDynamics.NSteps == 20
    assert isinstance(source.settings.input, inputs.AMS)
    assert source.settings.input.MolecularDynamics.NSteps == 10


@pytest.mark.parametrize("wrapped", [False, True])
def test_md_job_rejects_other_input_model(wrapped: bool) -> None:
    inputs = pytest.importorskip("scm.inputs")
    model = inputs.DFTB()
    settings = Settings()
    settings.input = model
    with pytest.raises(TypeError, match="scm.inputs.AMS"):
        AMSMDJob(settings=settings if wrapped else model)


def test_nve_removes_ensemble_blocks_from_converted_settings() -> None:
    inputs = pytest.importorskip("scm.inputs")
    model = inputs.AMS()
    model.MolecularDynamics.Thermostat.append({"Type": "NHC"})
    model.MolecularDynamics.Barostat.Type = "MTK"

    job = AMSNVEJob(settings=model)

    assert "Thermostat" not in job.settings.input.ams.MolecularDynamics
    assert "Barostat" not in job.settings.input.ams.MolecularDynamics
    assert model.MolecularDynamics.is_present("Thermostat")


def test_converted_md_job_restart_hash_does_not_resolve_rkf_path() -> None:
    inputs = pytest.importorskip("scm.inputs")
    previous = AMSNVTJob(settings=inputs.AMS(), nsteps=10)
    job = AMSNVEJob(settings=inputs.AMS(), velocities=previous)

    with patch.object(previous.results, "rkfpath", return_value="/previous/ams.rkf") as rkfpath:
        initial_hash = job.hash_input()
        rkfpath.assert_not_called()
        assert "File /previous/ams.rkf" in job.get_input()
        assert job.hash_input() == initial_hash


def test_converted_md_job_restart_from_typed_ams_job() -> None:
    inputs = pytest.importorskip("scm.inputs")
    model = inputs.AMS()
    model.MolecularDynamics.NSteps = 8
    previous = AMSJob(settings=model)
    restarted = AMSNVEJob.restart_from(previous, use_prerun=True, nsteps=12)

    assert isinstance(restarted.settings.input, Settings)
    assert restarted.settings.input.ams.MolecularDynamics.NSteps == 12
    assert isinstance(previous.settings.input, inputs.AMS)


def test_npt_equilibrated_molecule_and_system_select_same_history_step() -> None:
    job = MagicMock()
    job.status = JobStatus.SUCCESSFUL
    results = AMSNPTResults(job=job)
    results.get_history_property = MagicMock(return_value=[100.0, 4.0, 8.0, 6.0])
    results.get_history_molecule = MagicMock(return_value="molecule")
    results.get_history_system = MagicMock(return_value="system")

    molecule, molecule_step = results.get_equilibrated_molecule(equilibration_fraction=0.25, return_index=True)
    system, system_step = results.get_equilibrated_system(equilibration_fraction=0.25, return_index=True)

    assert (molecule, molecule_step) == ("molecule", 4)
    assert (system, system_step) == ("system", 4)
    results.get_history_molecule.assert_called_once_with(4)
    results.get_history_system.assert_called_once_with(4)


def test_equilibrate_density_results_forwards_get_equilibrated_system() -> None:
    npt_results = MagicMock()
    npt_results.get_equilibrated_system.return_value = ("system", 12)
    job = MagicMock()
    job.status = JobStatus.SUCCESSFUL
    job.children = {"npt": MagicMock(results=npt_results)}
    results = EquilibrateDensityResults(job=job)

    output = results.get_equilibrated_system(equilibration_fraction=0.75, return_index=True)

    assert output == ("system", 12)
    npt_results.get_equilibrated_system.assert_called_once_with(equilibration_fraction=0.75, return_index=True)


@pytest.mark.parametrize("typed", [False, True])
def test_bond_switch_job_accepts_input(typed: bool) -> None:
    inputs = pytest.importorskip("scm.inputs") if typed else None
    molecule = _periodic_molecule()
    settings = inputs.AMS() if typed else Settings()

    with patch.object(AMSMovingRestraintBondSwitchJob, "_get_bond_switches", return_value=([(1, 1, 1.0)], [], [])):
        job = AMSMovingRestraintBondSwitchJob(
            molecule=molecule, visit_molecules=[molecule.copy()], kappa=[1000], settings=settings
        )

    assert "MOVINGRESTRAINT ARG=bond1_1_1" in job.settings.input.ams.MolecularDynamics.Plumed.Input


@pytest.mark.parametrize("typed", [False, True])
def test_rmsd_job_accepts_input(typed: bool, tmp_path: Path) -> None:
    inputs = pytest.importorskip("scm.inputs") if typed else None
    settings = inputs.AMS() if typed else Settings()
    job = AMSMovingRestraintRMSDJob(visit_molecules=[_periodic_molecule()], kappa=[1000], settings=settings)
    job.path = str(tmp_path)

    with patch("scm.plams.recipes.md.movingrestraints.molecules_to_rkf"):
        with patch.object(job, "_write_plumed_compatible_pdb"):
            job.prerun()

    assert "MOVINGRESTRAINT ARG=rmsd_mol0" in job.settings.input.ams.MolecularDynamics.Plumed.Input


@pytest.mark.parametrize("typed", [False, True])
def test_scan_density_job_accepts_input(typed: bool) -> None:
    inputs = pytest.importorskip("scm.inputs") if typed else None
    settings = inputs.AMS() if typed else Settings()
    job = AMSMDScanDensityJob(molecule=_periodic_molecule(), settings=settings)

    deformation = job.settings.input.ams.MolecularDynamics.Deformation
    assert len(deformation.TargetLength.split()) == 3
    assert deformation.StartStep == 1

import pytest

pytest.importorskip("scm.base")
pytest.importorskip("scm.conformers")

from scm.base import ChemicalSystem
from scm.conformers import ConformersJob
from scm.plams import Atom, Molecule, Settings
from scm.plams.core.enums import JobStatus
from scm.plams.recipes.adfcosmorsconformers import ADFCOSMORSConfFilter, ADFCOSMORSConfJob


@pytest.mark.parametrize("use_system", [False, True])
@pytest.mark.parametrize("input_kind", ["legacy", "typed", "wrapped"])
def test_refinement_inputs(use_system, input_kind, tmp_path) -> None:
    molecule = Molecule()
    molecule.add_atom(Atom(symbol="He"))
    structure = ChemicalSystem.from_smiles("[He]", optimize_with_uff=False) if use_system else molecule
    if input_kind == "legacy":
        settings = Settings()
        settings.input.AMS.Task = "Optimize"
        settings.input.DFTB.Model = "GFN1-xTB"
    else:
        inputs = pytest.importorskip("scm.inputs")
        model = inputs.Conformers(Task="Optimize", Engine=inputs.DFTB(Model="GFN1-xTB"))
        original = model.to_input()
        settings = model
        if input_kind == "wrapped":
            settings = Settings()
            settings.input = model
            settings.runscript.nproc = 2

    job = ADFCOSMORSConfJob(
        structure,
        first_filter=ADFCOSMORSConfFilter(20, 22),
        additional=[(settings, ADFCOSMORSConfFilter(10, 12))],
    )
    job.prerun()
    job.children["job_0"].path = str(tmp_path / "generator")
    job.children["job_0"].status = JobStatus.SUCCESSFUL
    refinement = job._make_intermediate_job(1)
    assert isinstance(refinement, ConformersJob)
    assert refinement.settings.input.AMS.Task == "Optimize"
    assert refinement.settings.input.DFTB.Model == "GFN1-xTB"
    assert refinement.settings.input.AMS.InputMaxConfs == 20
    assert refinement.settings.input.AMS.InputMaxEnergy == 22
    assert str(tmp_path / "generator" / "conformers.rkf") in refinement.get_input()
    assert "InputMaxConfs" not in job.job_settings[1].input.AMS
    if input_kind != "legacy":
        assert model.to_input() == original
    if input_kind == "wrapped":
        assert refinement.settings.runscript.nproc == 2

    job.children["job_1"] = refinement
    adf = job._make_adf_job(2)
    assert adf.settings.input.AMS.InputMaxConfs == 10
    assert adf.settings.input.AMS.InputMaxEnergy == 12


@pytest.mark.parametrize("wrapped", [False, True])
def test_reject_wrong_typed_program(wrapped) -> None:
    inputs = pytest.importorskip("scm.inputs")
    settings = inputs.AMS(Engine=inputs.DFTB())
    if wrapped:
        outer = Settings()
        outer.input = settings
        settings = outer
    molecule = Molecule()
    molecule.add_atom(Atom(symbol="He"))
    with pytest.raises(TypeError, match="scm.inputs.Conformers"):
        ADFCOSMORSConfJob(molecule, additional=[(settings, None)])


def test_supplied_typed_generator() -> None:
    inputs = pytest.importorskip("scm.inputs")
    structure = ChemicalSystem.from_smiles("CC(=O)O", optimize_with_uff=True)
    settings = inputs.Conformers(Task="Generate")
    settings.Generator.RDKit.InitialNConformers = 50
    generator = ConformersJob(molecule=structure, settings=settings)
    original = settings.to_input()
    job = ADFCOSMORSConfJob(
        structure,
        conf_gen=generator,
        additional=[(inputs.Conformers(Task="Optimize", Engine=inputs.DFTB()), None)],
    )
    job.prerun()
    assert job.children["job_0"] is generator
    assert job.mol_info["Formula"] == structure.formula()
    text = generator.get_input()
    assert "50" in text
    assert "Atoms" in text
    assert settings.to_input() == original

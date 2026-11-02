from pathlib import Path

import pytest

from scm.plams import Atom, Molecule, Settings
from scm.plams.recipes.bandfragment import BANDFragmentJob, NOCVBandFragmentJob


@pytest.fixture
def typed_inputs():
    return pytest.importorskip("scm.inputs")


@pytest.fixture
def chemical_system():
    base = pytest.importorskip("scm.base")
    return base.ChemicalSystem.from_smiles("O")


@pytest.mark.parametrize("custom_mapping", [False, True])
def test_band_fragment_creates_full_job_and_applies_fragment_settings(
    typed_inputs, chemical_system, custom_mapping: bool
) -> None:
    full_settings = Settings()
    full_settings.input.BAND.XC.GGA = "PBE"
    if custom_mapping:
        full_settings.input.BAND.Fragment = [Settings({"FileName": "custom.rkf", "AtomMapping": {"1": "1"}})]
    job = BANDFragmentJob(
        fragment1=chemical_system.extract_atoms([0]),
        fragment2=chemical_system.extract_atoms([1, 2]),
        full_settings=full_settings,
        frag1_settings=typed_inputs.AMS(Engine=typed_inputs.BAND(Basis={"Type": "DZP"})),
        frag2_settings=typed_inputs.AMS(Engine=typed_inputs.BAND(Basis={"Type": "TZP"})),
    )

    job.prerun()
    children = job.new_children()

    assert children == [job.full]
    assert job.f1.settings.input.BAND.Basis.Type == "DZP"
    assert job.f2.settings.input.BAND.Basis.Type == "TZP"
    assert job.full.settings.input.BAND.XC.GGA == "PBE"
    assert isinstance(job.full.molecule, type(chemical_system))
    if custom_mapping:
        assert job.full.settings.input.BAND.Fragment[0].FileName == "custom.rkf"
    else:
        assert len(job.full.settings.input.BAND.Fragment) == 2
    assert job.new_children() is None


def test_band_fragment_accepts_typed_input(typed_inputs, chemical_system) -> None:
    job = BANDFragmentJob(
        fragment1=chemical_system.extract_atoms([0]),
        fragment2=chemical_system.extract_atoms([1, 2]),
        full_settings=typed_inputs.AMS(Engine=typed_inputs.BAND()),
    )

    job.prerun()

    assert len(job.children) == 2
    assert isinstance(job.full_settings, Settings)


def test_nocv_band_fragment_accepts_typed_settings(typed_inputs, tmp_path: Path) -> None:
    fragment = Molecule()
    fragment.add_atom(Atom(symbol="He"))
    model = typed_inputs.AMS(Engine=typed_inputs.BAND(Basis={"Type": "TZP"}))
    job = NOCVBandFragmentJob(fragment1=fragment, fragment2=fragment, nocv_settings=model)

    job.prerun()
    assert job.new_children() == [job.full]
    for child in (job.f1, job.f2, job.full):
        child.path = str(tmp_path / child.name)
    assert job.new_children() == [job.nocv]

    assert job.nocv.settings.input.BAND.Basis.Type == "TZP"
    assert job.nocv.settings.input.BAND.Restart.File == "full.rkf"
    assert len(job.nocv.settings.input.BAND.Fragment) == 2
    assert job.new_children() is None
    assert model.Engine.Basis.Type == "TZP"

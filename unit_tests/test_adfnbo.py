import pytest

from scm.plams import Atom, Molecule, Settings
from scm.plams.recipes.adfnbo import ADFNBOJob


@pytest.mark.parametrize("wrapped", [False, True])
def test_adfnbo_accepts_complete_ams_model(tmp_path, wrapped: bool) -> None:
    inputs = pytest.importorskip("scm.inputs")
    model = inputs.AMS(Engine=inputs.ADF())
    settings = Settings()
    settings.input = model
    settings.runscript.nproc = 2
    job = ADFNBOJob(settings=settings if wrapped else model)
    job.path = str(tmp_path)

    job.prerun()

    assert job.settings.input.ADF.FullFock
    assert job.settings.input.ADF.Basis.Core == "None"
    assert not model.Engine.is_present("FullFock")
    if wrapped:
        assert job.settings.runscript.nproc == 2


def test_adfnbo_rejects_bare_engine_model() -> None:
    inputs = pytest.importorskip("scm.inputs")
    with pytest.raises(TypeError, match="scm.inputs.AMS"):
        ADFNBOJob(settings=inputs.ADF())


@pytest.mark.parametrize("typed", [False, True])
def test_adfnbo_preserves_positional_molecule(typed: bool) -> None:
    molecule = Molecule()
    molecule.add_atom(Atom(symbol="He"))
    if typed:
        inputs = pytest.importorskip("scm.inputs")
        settings = inputs.AMS(Engine=inputs.ADF())
    else:
        settings = Settings()
        settings.input.ADF.Basis.Type = "DZ"

    job = ADFNBOJob(molecule, settings=settings)

    assert len(job.molecule) == 1
    assert job.molecule[1].symbol == "He"
    assert isinstance(job.settings.input, Settings)


def test_adfnbo_accepts_molecule_without_settings() -> None:
    molecule = Molecule()
    molecule.add_atom(Atom(symbol="He"))
    job = ADFNBOJob(molecule, name="nbo")

    assert job.name == "nbo"
    assert len(job.molecule) == 1

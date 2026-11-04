from unittest.mock import MagicMock

import pytest

from scm.plams import Atom, Molecule, Settings
from scm.plams.core.enums import JobStatus
from scm.plams.recipes.adffragment import ADFFragmentJob, ADFFragmentResults


@pytest.mark.parametrize("use_system", [False, True])
@pytest.mark.parametrize("typed", [False, True])
def test_adf_fragment_preserves_structures_and_input_settings(use_system: bool, typed: bool) -> None:
    fragment1 = Molecule()
    fragment1.add_atom(Atom(symbol="O"))
    fragment2 = Molecule()
    fragment2.add_atom(Atom(symbol="H", coords=(0, 0, 1)))
    fragment2.add_atom(Atom(symbol="H", coords=(0, 1, 0)))
    if use_system:
        base = pytest.importorskip("scm.base")
        from scm.plams import AMSJob

        fragment1 = base.ChemicalSystem(AMSJob(molecule=fragment1).get_input())
        fragment2 = base.ChemicalSystem(AMSJob(molecule=fragment2).get_input())
    if typed:
        inputs = pytest.importorskip("scm.inputs")
        settings = inputs.AMS(Engine=inputs.ADF())
    else:
        settings = Settings()
        settings.input.ADF.Basis.Type = "DZ"

    original_fragments = [str(fragment1), str(fragment2)]
    job = ADFFragmentJob(
        fragment1=fragment1,
        fragment2=fragment2,
        fragment1_opt=fragment1,
        fragment2_opt=fragment2,
        settings=settings,
        full_settings=settings,
        frag1_settings=settings,
        frag2_settings=settings,
    )

    job.prerun()

    assert isinstance(job.full.molecule, type(fragment1))
    assert len(job.full.molecule) == 3
    assert isinstance(job.f1_opt.molecule, type(fragment1))
    assert isinstance(job.f2_opt.molecule, type(fragment2))
    if use_system:
        assert all(atom.adf.f == "subsystem1" for atom in job.fragment1.atoms)
        assert all(atom.adf.f == "subsystem2" for atom in job.fragment2.atoms)
    else:
        assert all(atom.properties.suffix == "adf.f=subsystem1" for atom in job.fragment1)
        assert all(atom.properties.suffix == "adf.f=subsystem2" for atom in job.fragment2)
    assert "adf.f=" not in job.f1.get_input()
    assert "adf.f=" not in job.f2.get_input()
    assert [str(fragment1), str(fragment2)] == original_fragments
    if typed:
        assert not settings.Engine.is_present("Fragments")
    else:
        assert "Fragments" not in settings.input.ADF


def test_fragment_results_forward_system_methods() -> None:
    full_results = MagicMock()
    full_results.get_main_system.return_value = "main system"
    full_results.get_input_system.return_value = "input system"
    job = MagicMock()
    job.status = JobStatus.SUCCESSFUL
    job.full.results = full_results
    results = ADFFragmentResults(job=job)

    assert results.get_main_system() == "main system"
    assert results.get_input_system() == "input system"

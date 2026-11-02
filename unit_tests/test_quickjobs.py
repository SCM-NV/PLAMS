from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from scm.plams import Atom, Molecule
from scm.plams.core.errors import PlamsError
from scm.plams.interfaces.adfsuite import quickjobs


@pytest.mark.parametrize("container", [None, list, tuple])
def test_preoptimize_preserves_input_structure_types(container) -> None:
    base = pytest.importorskip("scm.base")
    system = base.ChemicalSystem.from_smiles("O")
    molecule = Molecule()
    molecule.add_atom(Atom(symbol="He"))
    structures = [system] if container is None else [system, molecule]
    results = []
    for structure in structures:
        result = MagicMock()
        result.get_main_system.return_value = structure.copy()
        result.get_main_molecule.return_value = structure.copy()
        results.append(result)

    with patch.object(quickjobs, "AMSWorker") as worker_type:
        worker = worker_type.return_value.__enter__.return_value
        worker.GeometryOptimization.side_effect = results
        output = quickjobs.preoptimize(system if container is None else container(structures))

    output_structures = [output] if container is None else output
    assert [type(structure) for structure in output_structures] == [type(structure) for structure in structures]
    results[0].get_main_system.assert_called_once_with()
    results[0].get_main_molecule.assert_not_called()
    if container is not None:
        assert isinstance(output, list)
        results[1].get_main_molecule.assert_called_once_with()
        results[1].get_main_system.assert_not_called()


def test_refined_system_preserves_metadata_and_leaves_input_unchanged() -> None:
    base = pytest.importorskip("scm.base")
    system = base.ChemicalSystem.from_smiles("[NH4+]")
    system.charge = 1
    original_coords = np.asarray(system.coords).copy()
    frame = Molecule()
    for atom in system.atoms:
        frame.add_atom(Atom(symbol=atom.symbol, coords=atom.coords + [1, 2, 3]))
    frame.lattice = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]

    output = quickjobs._with_updated_geometry(system, frame)

    assert isinstance(output, base.ChemicalSystem)
    assert output.charge == 1
    assert len(output.bonds) == len(system.bonds)
    assert np.allclose(output.coords, frame.as_array())
    assert np.allclose(output.lattice.vectors, frame.lattice)
    assert np.allclose(system.coords, original_coords)


def test_shakemd_reports_missing_final_structure() -> None:
    molecule = Molecule()
    with patch.object(quickjobs, "_ensure_init", return_value=False), patch.object(quickjobs, "AMSJob") as job_type:
        with patch.object(quickjobs, "delete_job") as delete_job:
            job_type.return_value.results.get_main_molecule.return_value = None
            with pytest.raises(PlamsError, match="SHAKE MD failed"):
                quickjobs.shakemd(molecule)
            delete_job.assert_not_called()

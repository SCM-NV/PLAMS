import pytest
from pathlib import Path
import numpy as np

from scm.plams.mol.molecule import Molecule
from scm.plams.trajectories.rkfhistoryfile import RKFHistoryFile


@pytest.fixture
def conformers_rkf(rkf_folder):
    return str(Path(rkf_folder) / "conformers" / "conformers.rkf")


class TestRKFHistoryFile:

    def test_frames(self, conformers_rkf):
        # Given rkf file with multiple conformers
        history_file = RKFHistoryFile(conformers_rkf)

        num_frames = history_file.get_length()
        assert num_frames == 12

        input_mol = history_file.get_plamsmol()
        assert input_mol.get_formula() == "C2H7NO"

        # When read successive frames
        for i in range(1, num_frames):
            mol = Molecule()
            crds, _ = history_file.read_frame(i, mol)

            # Then coordinates are different to the original molecule
            # But the labels still match considering bond connectivity
            assert not np.allclose(crds, input_mol.as_array())
            assert mol.label(3) == input_mol.label(3)

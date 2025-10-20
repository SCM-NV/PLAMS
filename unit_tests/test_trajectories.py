import pytest
from pathlib import Path
import numpy as np

from scm.plams.mol.molecule import Molecule
from scm.plams.trajectories.rkfhistoryfile import RKFHistoryFile
from scm.plams.trajectories.xyzfile import XYZTrajectoryFile
from scm.plams.trajectories.xyzhistoryfile import XYZHistoryFile

@pytest.fixture
def conformers_rkf(rkf_folder):
    return str(Path(rkf_folder) / "conformers" / "conformers.rkf")

@pytest.fixture
def trajectory_xyz(xyz_folder):
    return str(Path(xyz_folder) / "water_box_traj.xyz")


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


class TestXYZFiles:

    # def test_historyfile(self, trajectory_xyz):
    #
    #     xyz = XYZHistoryFile(trajectory_xyz)
    #     mol = xyz.get_plamsmol()
    #
    #     xyzout = XYZHistoryFile('new.xyz',mode='w')
    #     for i in range(xyz.get_length()) :
    #         crd,cell = xyz.read_frame(i,molecule=mol)
    #         xyzout.write_next(molecule=mol)
    #     xyzout.close()


    def test_trajectoryfile(self, trajectory_xyz):

        trajectory_file = XYZTrajectoryFile(trajectory_xyz)
        mol = trajectory_file.get_plamsmol()

        xyzout = XYZTrajectoryFile('new.xyz',mode='w')
        xyzout.style = "scm"
        for i in range(trajectory_file.get_length()) :
            crd,cell = trajectory_file.read_frame(i,molecule=mol)
            xyzout.write_next(molecule=mol)
        xyzout.close()

        with open('new.xyz') as f: output = f.read()
        expectedOutput = """\
3

       O        -0.1995727721        -0.2577814972        -0.0000000000 
       H        -0.4018214824        -0.2097333578        -0.9372229380 
       H         0.7422877088        -0.4019259155         0.1171376168 
VEC0         6.0000000000         0.0000000000         0.0000000000 
VEC1         3.0000000000         5.2000000000         0.0000000000 
VEC2         0.0000000000         0.0000000000         6.0000000000 \
"""

        assert output[0:428] == expectedOutput

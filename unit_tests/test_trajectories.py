import pytest
from pathlib import Path
import numpy as np

from scm.plams.mol.molecule import Molecule
from scm.plams.trajectories.rkfhistoryfile import RKFHistoryFile
from scm.plams.trajectories.rkffile import RKFTrajectoryFile
from scm.plams.trajectories.xyzfile import XYZTrajectoryFile
from scm.plams.trajectories.sdffile import SDFTrajectoryFile


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


class TestTrajectoryFileFormats:

    def test_trajectoryfile(self, conformers_rkf):
        """
        Check that XYZ trajectory files in SCM format can be read and written
        """
        from io import StringIO

        classes = [XYZTrajectoryFile, SDFTrajectoryFile]

        rkf = RKFTrajectoryFile(conformers_rkf)
        mol = rkf.get_plamsmol()
        nats = len(mol)
        nframes = len(rkf)

        # Add periodic boundary
        cell = np.zeros((3, 3))
        np.fill_diagonal(cell, 10)
        mol.lattice = cell.tolist()

        for trajcls in classes:
            # Write periodic SCM style trajectory
            fileobj = StringIO()
            outfile = trajcls(fileobject=fileobj, mode="w")
            outfile.style = "scm"

            for i in range(nframes):
                rkf.read_frame(i, molecule=mol)
                outfile.write_next(molecule=mol)

            # Read the XYZ file into a string
            fileobj.seek(0)
            lines = fileobj.readlines()
            fileobj.seek(0)

            if trajcls is XYZTrajectoryFile:
                assert len(lines) == nframes * (nats + 5)
                for i in range(3):
                    iline = i + nats + 2
                    assert lines[iline].split()[0] == "VEC%i" % (i + 1)
            else:
                assert "VEC1" in "\n".join(lines)

            # Read the XYZ file, and check that the periodicity is retained
            infile = trajcls(fileobject=fileobj)
            mol_pbc = infile.get_plamsmol()
            assert hasattr(mol_pbc, "lattice")

            outfile.close()

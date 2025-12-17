import pytest
import os
import random
from pathlib import Path
import numpy as np

from scm.plams.mol.molecule import Molecule
from scm.plams.trajectories.rkfhistoryfile import RKFHistoryFile
from scm.plams.trajectories.rkffile import RKFTrajectoryFile
from scm.plams.trajectories.xyzfile import XYZTrajectoryFile
from scm.plams.trajectories.sdffile import SDFTrajectoryFile

from test_helpers import skip_if_no_ams_installation


@pytest.fixture
def conformers_rkf(rkf_folder):
    return str(Path(rkf_folder) / "conformers" / "conformers.rkf")


@pytest.fixture
def trajectory_xyz(xyz_folder):
    return str(Path(xyz_folder) / "water_box_traj.xyz")


@pytest.fixture
def molecules(xyz_folder):
    """
    Return set of single molecules
    """
    path = str(Path(xyz_folder))
    filenames = [os.path.join(path, fn) for fn in os.listdir(path)]

    molecules = []
    for i, fn in enumerate(filenames):
        mol = Molecule(fn)
        mol.guess_bonds()
        if len(mol.lattice) > 0:
            continue
        nmols = len(mol.separate())
        if nmols > 1:
            continue
        molecules.append(mol)
    return molecules


class TestRKFHistoryFile:

    @pytest.fixture
    def rkffilename(self, tmp_path_factory, molecules):
        """
        Write RKFHistoryFile with properties NOT stored as blocks
        """
        skip_if_no_ams_installation()
        rkfname = (tmp_path_factory.mktemp("data") / "molecules.rkf").as_posix()
        rkf = RKFHistoryFile(rkfname, mode="wb")
        rkf.store_historydata()
        rkf.store_mddata()

        # Store iframe as a list of integers, but not every step
        for iframe, mol in enumerate(molecules):
            historydata = {"Step": iframe, "Energy": 0.0}
            mddata = {}
            if iframe % 2 != 0:
                mddata = {"ListOfInts": [iframe]}
            rkf.write_next(molecule=mol, historydata=historydata, mddata=mddata)
        rkf.close()
        return rkfname

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

    def test_property_reading(self, rkffilename):
        """
        Test reading of properties from the RKFHistoryFile
        """
        assert os.path.isfile(rkffilename)

        rkf = RKFHistoryFile(rkffilename)
        rkf.store_mddata()
        assert len(rkf) >= 23

        indices = [i for i in range(len(rkf))]
        indices = random.sample(indices, len(rkf))
        results = []
        for iframe in indices:
            crd, cell = rkf.read_frame(iframe)
            if "ListOfInts" in rkf.mddata:
                results.append(rkf.mddata["ListOfInts"])
                assert rkf.mddata["ListOfInts"] == iframe
        assert len(results) > 0
        assert len(results) < len(rkf)


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

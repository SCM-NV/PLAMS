import os
import numpy as np
import pytest

from scm.plams import Settings
from scm.plams import Molecule
from scm.plams import AMSJob
from scm.plams import AMSMSDJob
from scm.plams import AMSRDFJob

from test_helpers import skip_if_no_ams_installation
from test_helpers import skip_if_no_scm_pisa


@pytest.fixture(scope="session")
def mdjob(tmp_path_factory, xyz_folder):
    """
    Create AMSJob for an MD simulation (hopefully once per session)
    """
    skip_if_no_ams_installation()
    pwd = os.getcwd()
    mol = Molecule(xyz_folder / "water_box.xyz")

    dirname = tmp_path_factory.mktemp("data")
    os.chdir(dirname.as_posix())

    s = Settings()
    s.input.ams.Task = "MolecularDynamics"
    s.input.ReaxFF.ForceField = "Water2017.ff"
    s.input.ams.RNGSeed = "1 2 3 4 5 6 7 8 9"
    s.input.ams.MolecularDynamics.CalcPressure = "Yes"
    s.input.ams.MolecularDynamics.InitialVelocities.Temperature = 300
    s.input.ams.MolecularDynamics.Trajectory.SamplingFreq = 1
    s.input.ams.MolecularDynamics.TimeStep = 0.5
    s.input.ams.MolecularDynamics.NSteps = 200

    s.runscript.nproc = 1
    os.environ["OMP_NUM_THREADS"] = "1"
    job = AMSJob(settings=s, molecule=mol, name="md")
    job.run()
    os.chdir(pwd)

    return job


class TestMSDJob:
    """
    Test of the AMS Job representing a MeanSquareDisplacement calculation
    """

    def test_default_use(self, mdjob):
        """
        Test plainest use of AMSMSDJob
        """
        mol = mdjob.results.get_main_molecule()
        oxygens = [i + 1 for i, at in enumerate(mol.atoms) if at.symbol == "O"]
        msd_job = AMSMSDJob(mdjob, atom_indices=oxygens, start_time_fit_fs=20)

        # Check the input
        txt = msd_job.get_input()
        for iat in oxygens:
            assert "Atom %i" % (iat) in txt

        # Run the job and check output
        msd_job.run()
        D = msd_job.results.get_diffusion_coefficient()
        assert abs(D - 2.7430808506924798e-08) < 1e-18

    def test_settings_use(self, mdjob):
        """
        Test AMSMSDJob with settings object
        """
        s = Settings()
        s.input.MeanSquareDisplacement.Atoms.Element = "O"
        msd_job = AMSMSDJob(mdjob, settings=s, start_time_fit_fs=20)

        # Check the input
        txt = msd_job.get_input()
        assert "Element O" in txt

        # Run the job and check output
        msd_job.run()
        D = msd_job.results.get_diffusion_coefficient()
        assert abs(D - 2.7430808506924798e-08) < 1e-18

    def test_pisa_use(self, mdjob):
        """
        Test AMSMSDJob with Pisa object
        """
        skip_if_no_scm_pisa()
        from scm.input_classes import Analysis

        sets = Analysis()
        sets.MeanSquareDisplacement.Atoms.Element = "O"
        msd_job = AMSMSDJob(mdjob, settings=sets, start_time_fit_fs=20)

        # Check the input
        txt = msd_job.get_input()
        assert "Element O" in txt

        # Run the job and check output
        msd_job.run()
        D = msd_job.results.get_diffusion_coefficient()
        assert abs(D - 2.7430808506924798e-08) < 1e-18


class TestRDFJob:
    """
    Test of the AMS Job representing a RadialDistribution calculation
    """

    def test_default_use(self, mdjob):
        """
        Test plainest use of AMSRDFJob
        """
        mol = mdjob.results.get_main_molecule()
        nats = len(mol)
        rdf_job = AMSRDFJob(mdjob)

        # Check the input
        txt = rdf_job.get_input()
        for iat in range(nats):
            assert "Atom %i" % (iat + 1) in txt

        # Run the job and check output
        rdf_job.run()
        x, rdf = rdf_job.results.get_rdf()
        ind = rdf.argmax()
        peak = x[ind]
        assert abs(peak - 0.9074074074074073) < 1e-8

    def test_indices_use(self, mdjob):
        """
        Test AMSRDFJob supplying atom indices as arguments
        """
        mol = mdjob.results.get_main_molecule()
        oxygens = [i + 1 for i, at in enumerate(mol.atoms) if at.symbol == "O"]
        rdf_job = AMSRDFJob(mdjob, atom_indices=oxygens, atom_indices_to=oxygens)

        # Check the input
        txt = rdf_job.get_input()
        for iat in oxygens:
            assert "Atom %i" % (iat) in txt
        assert not "Atom 2\n" in txt

        # Run the job and check output
        rdf_job.run()
        x, rdf = rdf_job.results.get_rdf()
        ind = rdf.argmax()
        peak = x[ind]
        assert abs(peak - 2.7407407407407405) < 1e-8

    def test_settings_use(self, mdjob):
        """
        Test AMSRDFJob with a Settings object
        """
        mol = mdjob.results.get_main_molecule()
        oxygens = [i + 1 for i, at in enumerate(mol.atoms) if at.symbol == "O"]

        s = Settings()
        s.input.RadialDistribution.AtomsFrom.Atom = oxygens
        s.input.RadialDistribution.AtomsTo.Atom = oxygens

        rdf_job = AMSRDFJob(mdjob, atom_indices=oxygens, atom_indices_to=oxygens, settings=s)

        # Check the input
        failed = False
        try:
            txt = rdf_job.get_input()
        except Exception:
            failed = True
        assert failed

        # Now do this correctly
        rdf_job = AMSRDFJob(mdjob, settings=s)
        txt = rdf_job.get_input()
        for iat in oxygens:
            assert "Atom %i" % (iat) in txt
        assert not "Atom 2\n" in txt

        # Run the job and check output
        rdf_job.run()
        x, rdf = rdf_job.results.get_rdf()
        ind = rdf.argmax()
        peak = x[ind]
        assert abs(peak - 2.7407407407407405) < 1e-8

    def test_pisa_use(self, mdjob):
        """
        Test AMSRDFJob with Pisa settings
        """
        skip_if_no_scm_pisa()
        from scm.input_classes import Analysis

        sets = Analysis()
        sets.RadialDistribution.AtomsFrom.Element = "O"
        sets.RadialDistribution.AtomsTo.Element = "O"

        rdf_job = AMSRDFJob(mdjob, settings=sets)
        txt = rdf_job.get_input()
        assert "Element O" in txt

        # Run the job and check output
        rdf_job.run()
        x, rdf = rdf_job.results.get_rdf()
        ind = rdf.argmax()
        peak = x[ind]
        assert abs(peak - 2.7407407407407405) < 1e-8

    def test_multiple_rdfs(self, mdjob):
        """
        Test using multiple RadialDistribution blocks in the settings
        """
        mol = mdjob.results.get_main_molecule()
        oxygens = [i + 1 for i, at in enumerate(mol.atoms) if at.symbol == "O"]

        # Create the settings with multiply RadialDistribution blocks
        s = Settings()
        subsettings = [Settings(), Settings()]
        subsettings[0].AtomsFrom.Atom = oxygens
        subsettings[0].AtomsTo.Atom = oxygens
        s.input.RadialDistribution = subsettings

        # Create the RDFJob, test supplying the job as a list, and then check the input
        rdf_job = AMSRDFJob([mdjob], settings=s)
        txt = rdf_job.get_input()
        for iat in oxygens:
            assert "Atom %i" % (iat) in txt

        # Run the job and check output
        rdf_job.run()

        # Get the peak from the first RDF
        x, rdf = rdf_job.results.get_rdf()
        ind = rdf.argmax()
        peak = x[ind]
        assert abs(peak - 2.7407407407407405) < 1e-8

        # Get the peak from the second RDF
        xy = rdf_job.results.get_xy(i=2)
        x = np.array(xy.x[0])
        rdf = np.array(xy.y)
        peak = x[rdf.argmax()]
        assert abs(peak - 0.9074074074074073) < 1e-8

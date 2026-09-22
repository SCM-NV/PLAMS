import os
import numpy as np
import pytest

from scm.plams import Settings
from scm.plams import Molecule
from scm.plams import AMSJob
from scm.plams import AMSMSDJob
from scm.plams import AMSRDFJob
from scm.plams import AMSVACFJob
from scm.plams.recipes.md.trajectoryanalysis import AMSViscosityFromBinLogJob

from test_helpers import skip_if_no_ams_installation
from test_helpers import skip_if_no_scm_inputs


@pytest.fixture(scope="session")
def mdjob(xyz_folder):
    """
    Create AMSJob for an MD simulation once per session
    """
    skip_if_no_ams_installation()
    mol = Molecule(xyz_folder / "water_box.xyz")

    s = Settings()
    s.input.ams.Task = "MolecularDynamics"
    s.input.ReaxFF.ForceField = "Water2017.ff"
    s.input.ams.RNGSeed = "1 2 3 4 5 6 7 8 9"
    s.input.ams.MolecularDynamics.CalcPressure = "Yes"
    s.input.ams.MolecularDynamics.InitialVelocities.Temperature = 300
    s.input.ams.MolecularDynamics.Trajectory.SamplingFreq = 1
    s.input.ams.MolecularDynamics.TimeStep = 0.5
    s.input.ams.MolecularDynamics.NSteps = 200
    s.input.ams.MolecularDynamics.BinLog.Step = "Yes"
    s.input.ams.MolecularDynamics.BinLog.PressureTensor = "Yes"

    s.runscript.nproc = 1
    os.environ["OMP_NUM_THREADS"] = "1"
    job = AMSJob(settings=s, molecule=mol, name="md")
    job.run()

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
        assert np.isclose(D, 2.7430808506924798e-08, atol=1e-13)

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
        assert np.isclose(D, 2.7430808506924798e-08, atol=1e-13)

    def test_scm_inputs_use(self, mdjob):
        """
        Test AMSMSDJob with an scm.inputs model
        """
        skip_if_no_scm_inputs()
        from scm.inputs import Analysis

        sets = Analysis()
        msd = Analysis.MeanSquareDisplacementBlock()
        msd.Atoms.Element = ["O"]
        sets.MeanSquareDisplacement = [msd]
        msd_job = AMSMSDJob(mdjob, settings=sets, start_time_fit_fs=20)

        # Check the input
        txt = msd_job.get_input()
        assert "Element O" in txt

        # Run the job and check output
        msd_job.run()
        D = msd_job.results.get_diffusion_coefficient()
        assert np.isclose(D, 2.7430808506924798e-08, atol=1e-13)


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
        assert np.isclose(peak, 0.9074074074074073, atol=1e-8)

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
        assert np.isclose(peak, 2.7407407407407405, atol=1e-8)

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
        assert np.isclose(peak, 2.7407407407407405, atol=1e-8)

    def test_scm_inputs_use(self, mdjob):
        """
        Test AMSRDFJob with an scm.inputs model
        """
        skip_if_no_scm_inputs()
        from scm.inputs import Analysis

        sets = Analysis()
        rdf = Analysis.RadialDistributionBlock()
        rdf.AtomsFrom.Element = ["O"]
        rdf.AtomsTo.Element = ["O"]
        sets.RadialDistribution = [rdf]

        rdf_job = AMSRDFJob(mdjob, settings=sets)
        txt = rdf_job.get_input()
        assert "Element O" in txt

        # Run the job and check output
        rdf_job.run()
        x, rdf = rdf_job.results.get_rdf()
        ind = rdf.argmax()
        peak = x[ind]
        assert np.isclose(peak, 2.7407407407407405, atol=1e-8)

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
        assert np.isclose(peak, 2.7407407407407405, atol=1e-8)

        # Get the peak from the second RDF
        xy = rdf_job.results.get_xy(i=2)
        x = np.array(xy.x[0])
        rdf = np.array(xy.y)
        peak = x[rdf.argmax()]
        assert np.isclose(peak, 0.9074074074074073, atol=1e-8)


class TestVACFJob:
    """
    Test of the AMS Job representing a velocity autocorrelation function calculation
    """

    def test_default_use(self, mdjob):
        """
        Test plainest use of AMSVACFJob
        """
        mol = mdjob.results.get_main_molecule()
        oxygens = [i + 1 for i, at in enumerate(mol.atoms) if at.symbol == "O"]
        vacf_job = AMSVACFJob(mdjob, atom_indices=oxygens)

        # Check the input
        txt = vacf_job.get_input()
        for iat in oxygens:
            assert "Atom %i" % (iat) in txt

        # Run the job and check output
        vacf_job.run()
        D, units = vacf_job.results.get_D()
        assert np.isclose(D, 3.2882389727568384e-08, atol=1e-13)


class TestViscosityFromBinLogJob:
    """
    Test of the AMS Job representing a velocity autocorrelation function calculation
    """

    def test_default_use(self, mdjob):
        """
        Test plainest use of AMSVACFJob
        """
        visc_job = AMSViscosityFromBinLogJob(mdjob)

        # Run the job and check output
        visc_job.run()
        visc = visc_job.results.get_double_exponential_fit()[-1][-1]
        assert np.isclose(visc, 8.855124346230889e-05, atol=1e-08)

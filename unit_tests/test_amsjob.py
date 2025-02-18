import re
import dill as pickle
import pytest
from unittest.mock import MagicMock, patch
from collections import namedtuple
from io import StringIO
import os
import time
import threading

from scm.plams.interfaces.adfsuite.ams import AMSJob, AMSResults
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Atom, Molecule
from scm.plams.unit_tests.test_helpers import skip_if_no_scm_pisa, skip_if_no_scm_libbase


class TestAMSJob:
    """
    Test suite for AMSJob without using PISA / CS for input.
    Sets up a geometry optimization of water.
    """

    @pytest.fixture
    def job_input(self):
        JobInput = namedtuple("JobInput", "molecule settings input")
        molecule = self.get_input_molecule()
        settings = self.get_input_settings()
        expected_input = self.get_expected_input()
        return JobInput(molecule, settings, expected_input)

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        molecule = Molecule()
        molecule.add_atom(Atom(symbol="O", coords=(0, 0, 0)))
        molecule.add_atom(Atom(symbol="H", coords=(1, 0, 0)))
        molecule.add_atom(Atom(symbol="H", coords=(0, 1, 0)))
        return molecule

    @staticmethod
    def get_input_settings():
        """
        Get instance of the Settings class passed to the AMSJob
        """
        settings = Settings()
        settings.input.ams.Task = "GeometryOptimization"
        settings.input.ams.Properties.NormalModes = "Yes"
        settings.input.DFTB.Model = "GFN1-xTB"
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
Properties
  NormalModes Yes
End

Task GeometryOptimization

System
  Atoms
              O       0.0000000000       0.0000000000       0.0000000000
              H       1.0000000000       0.0000000000       0.0000000000
              H       0.0000000000       1.0000000000       0.0000000000
  End
End

Engine DFTB
  Model GFN1-xTB
EndEngine

"""

    def test_init_deep_copies_molecule(self, job_input):
        # Given job with molecule
        job = AMSJob(molecule=job_input.molecule)

        # When get molecule from job
        # Then job molecule is a deep copy
        if job_input.molecule is not None:
            assert job.molecule is not job_input.molecule
            if isinstance(job_input.molecule, dict):
                for name, mol in job.molecule.items():
                    assert mol is not job_input.molecule[name]
                    assert mol.atoms is not job_input.molecule[name].atoms
            else:
                assert job.molecule.atoms is not job_input.molecule.atoms
        else:
            assert job.molecule is None

    def test_pickle_dumps_and_loads_job_successfully(self, job_input):
        # Given job with molecule and settings
        job = AMSJob(molecule=job_input.molecule, settings=job_input.settings)

        # When round trip via pickling
        pickle_bytes = pickle.dumps(job)
        job2 = pickle.loads(pickle_bytes)

        # Then job is still of correct type
        assert isinstance(job2, AMSJob)

    def test_get_input_generates_expected_input_string(self, job_input):
        # Given job with molecule and settings
        job = AMSJob(molecule=job_input.molecule, settings=job_input.settings)

        # When get the job input for the input file
        # Then the input matches the expected input
        assert job.get_input() == job_input.input

    def test_get_runscript_generates_expected_string(self, job_input):
        # Given job
        job = AMSJob(molecule=job_input.molecule, settings=job_input.settings)

        # When get the runscript
        runscript = job.get_runscript()

        # Then standard runscript returned
        assert (
            runscript
            == """unset AMS_SWITCH_LOGFILE_AND_STDOUT
unset SCM_LOGFILE
AMS_JOBNAME="plamsjob" AMS_RESULTSDIR=. $AMSBIN/ams --input="plamsjob.in" < /dev/null

"""
        )

    def test_get_runscript_with_runscript_settings_generates_expected_string(self, job_input, config):
        # Given job with additional runscript settings
        job_input.settings.runscript.preamble_lines = ["# Start"]
        job_input.settings.runscript.postamble_lines = ["# End"]
        job_input.settings.runscript.nproc = 8
        job_input.settings.runscript.stdout_redirect = True
        job = AMSJob(molecule=job_input.molecule, settings=job_input.settings)

        # When get the runscript
        with patch("scm.plams.interfaces.adfsuite.ams.config", config):
            config.slurm = None  # Remove any specific settings when running under slurm
            runscript = job.get_runscript()

        # Then runscript with additional lines returned
        assert (
            runscript
            == """unset AMS_SWITCH_LOGFILE_AND_STDOUT
unset SCM_LOGFILE
# Start
AMS_JOBNAME="plamsjob" AMS_RESULTSDIR=. $AMSBIN/ams -n 8 --input="plamsjob.in" < /dev/null >"plamsjob.out"

# End

"""
        )

    @pytest.mark.parametrize(
        "status,expected",
        [
            ["NORMAL TERMINATION", True],
            ["NORMAL TERMINATION with warnings", True],
            ["NORMAL TERMINATION with errors", False],
            ["Input error", False],
            [None, False],
        ],
    )
    def test_check_returns_true_for_normal_termination_with_no_errors_otherwise_false(self, status, expected):
        # Given job with results of certain status
        job = AMSJob()
        job.results = MagicMock(spec=AMSResults)
        job.results.readrkf.return_value = status

        # When check the job
        # Then job check is ok only for normal termination with no errors
        assert job.check() == expected

    @pytest.mark.parametrize(
        "status,expected", [["NORMAL TERMINATION", None], ["NORMAL TERMINATION with errors", "something bad"]]
    )
    def test_get_errormsg_returns_message_from_logfile_on_error_otherwise_none(self, status, expected):
        # Given job with results of certain status
        job = AMSJob()
        job.results = MagicMock(spec=AMSResults)
        job.results.readrkf.return_value = status
        job.results.grep_file.return_value = ["ERROR: something bad"]

        # When get the error message
        # Then the message is none only for normal termination with no errors
        assert job.get_errormsg() == expected


class TestAMSJobWithPisa(TestAMSJob):
    """
    Test suite for AMSJob using PISA for settings input.
    Sets up a geometry optimization of water.
    """

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS
        from scm.input_classes.engines import DFTB

        settings = Settings()
        driver = AMS()
        driver.Task = "GeometryOptimization"
        driver.Properties.NormalModes = "True"
        driver.Engine = DFTB()
        driver.Engine.Model = "GFN1-xTB"
        settings.input = driver
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
Properties
  NormalModes True
End
Task GeometryOptimization

Engine DFTB
  Model GFN1-xTB
EndEngine

System
  Atoms
              O       0.0000000000       0.0000000000       0.0000000000
              H       1.0000000000       0.0000000000       0.0000000000
              H       0.0000000000       1.0000000000       0.0000000000
  End
End
"""


class TestAMSJobWithPisaOnly(TestAMSJob):
    """
    Test suite for AMSJob using PISA for settings and molecule input.
    Sets up a geometry optimization of water.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        return None

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS
        from scm.input_classes.engines import DFTB

        settings = Settings()
        driver = AMS()
        driver.Task = "GeometryOptimization"
        driver.Properties.NormalModes = "True"
        driver.Engine = DFTB()
        driver.Engine.Model = "GFN1-xTB"
        driver.System.Atoms = [
            "O       0.0000000000       0.0000000000       0.0000000000",
            "H       1.0000000000       0.0000000000       0.0000000000",
            "H       0.0000000000       1.0000000000       0.0000000000",
        ]
        settings.input = driver
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
Properties
  NormalModes True
End
System
  Atoms
    O       0.0000000000       0.0000000000       0.0000000000
    H       1.0000000000       0.0000000000       0.0000000000
    H       0.0000000000       1.0000000000       0.0000000000
  End
End
Task GeometryOptimization

Engine DFTB
  Model GFN1-xTB
EndEngine
"""


class TestAMSJobWithChemicalSystem(TestAMSJob):
    """
    Test suite for AMSJob using ChemicalSystem for molecule input
    """

    @staticmethod
    def get_input_molecule():
        """
        Instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

        molecule = ChemicalSystem()
        molecule.add_atom("O", coords=[0, 0, 0])
        molecule.add_atom("H", coords=[1, 0, 0])
        molecule.add_atom("H", coords=[0, 1, 0])
        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
Properties
  NormalModes Yes
End

Task GeometryOptimization

System
   Atoms
      O 0 0 0
      H 1 0 0
      H 0 1 0
   End
End

Engine DFTB
  Model GFN1-xTB
EndEngine

"""


class TestAMSJobWithChemicalSystemAndPisa(TestAMSJobWithPisa):
    """
    Test suite for AMSJob using ChemicalSystem for molecule input and PISA for settings input
    """

    @staticmethod
    def get_input_molecule():
        """
        Instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

        molecule = ChemicalSystem()
        molecule.add_atom("O", coords=[0, 0, 0])
        molecule.add_atom("H", coords=[1, 0, 0])
        molecule.add_atom("H", coords=[0, 1, 0])
        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
Properties
  NormalModes True
End
Task GeometryOptimization

Engine DFTB
  Model GFN1-xTB
EndEngine

System
   Atoms
      O 0 0 0
      H 1 0 0
      H 0 1 0
   End
End"""


class TestAMSJobWithMultipleMolecules(TestAMSJob):
    """
    Test suite for AMSJob using multiple molecules.
    Sets up a NEB calculation for the isomerisation of HCN.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        main_molecule = Molecule()
        main_molecule.add_atom(Atom(symbol="C", coords=(0, 0, 0)))
        main_molecule.add_atom(Atom(symbol="N", coords=(1.18, 0, 0)))
        main_molecule.add_atom(Atom(symbol="H", coords=(2.196, 0, 0)))
        final_molecule = main_molecule.copy()
        final_molecule.atoms[1].x = 1.163
        final_molecule.atoms[2].x = -1.078

        molecule = {"": main_molecule, "final": final_molecule}

        return molecule

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        settings = Settings()
        settings.input.ams.Task = "NEB"
        settings.input.ams.NEB.Images = 9
        settings.input.ams.NEB.Iterations = 100
        settings.input.DFTB.Model = "DFTB3"
        settings.input.DFTB.ResourcesDir = "3ob-3-1"
        settings.input.DFTB.DispersionCorrection = "D3-BJ"
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
NEB
  Images 9
  Iterations 100
End

Task NEB

System
  Atoms
              C       0.0000000000       0.0000000000       0.0000000000
              N       1.1800000000       0.0000000000       0.0000000000
              H       2.1960000000       0.0000000000       0.0000000000
  End
End
System final
  Atoms
              C       0.0000000000       0.0000000000       0.0000000000
              N       1.1630000000       0.0000000000       0.0000000000
              H      -1.0780000000       0.0000000000       0.0000000000
  End
End

Engine DFTB
  DispersionCorrection D3-BJ
  Model DFTB3
  ResourcesDir 3ob-3-1
EndEngine

"""


class TestAMSJobWithMultipleMoleculesAndPisa(TestAMSJobWithMultipleMolecules):
    """
    Test suite for AMSJob using multiple molecules and PISA for settings input.
    Sets up a NEB calculation for the isomerisation of HCN.
    """

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS
        from scm.input_classes.engines import DFTB

        settings = Settings()
        driver = AMS()
        driver.Task = "NEB"
        driver.NEB.Images = 9
        driver.NEB.Iterations = 100
        driver.Engine = DFTB()
        driver.Engine.Model = "DFTB3"
        driver.Engine.ResourcesDir = "3ob-3-1"
        driver.Engine.DispersionCorrection = "D3-BJ"
        settings.input = driver
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
NEB
  Images 9
  Iterations 100
End
Task NEB

Engine DFTB
  DispersionCorrection D3-BJ
  Model DFTB3
  ResourcesDir 3ob-3-1
EndEngine

System
  Atoms
              C       0.0000000000       0.0000000000       0.0000000000
              N       1.1800000000       0.0000000000       0.0000000000
              H       2.1960000000       0.0000000000       0.0000000000
  End
End

System final
  Atoms
              C       0.0000000000       0.0000000000       0.0000000000
              N       1.1630000000       0.0000000000       0.0000000000
              H      -1.0780000000       0.0000000000       0.0000000000
  End
End
"""


class TestAMSJobWithMultipleMoleculesAndPisaOnly(TestAMSJobWithMultipleMolecules):
    """
    Test suite for AMSJob using multiple molecules and PISA for settings and molecule input.
    Sets up a NEB calculation for the isomerisation of HCN.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        return None

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS
        from scm.input_classes.engines import DFTB

        settings = Settings()
        driver = AMS()
        driver.Task = "NEB"
        driver.NEB.Images = 9
        driver.NEB.Iterations = 100
        driver.Engine = DFTB()
        driver.Engine.Model = "DFTB3"
        driver.Engine.ResourcesDir = "3ob-3-1"
        driver.Engine.DispersionCorrection = "D3-BJ"
        driver.System[0].Atoms = [
            "C       0.0000000000       0.0000000000       0.0000000000",
            "N       1.1800000000       0.0000000000       0.0000000000",
            "H       2.1960000000       0.0000000000       0.0000000000",
        ]
        driver.System[1].header = "final"
        driver.System[1].Atoms = [
            "C       0.0000000000       0.0000000000       0.0000000000",
            "N       1.1630000000       0.0000000000       0.0000000000",
            "H      -1.0780000000       0.0000000000       0.0000000000",
        ]
        settings.input = driver
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
NEB
  Images 9
  Iterations 100
End
System
  Atoms
    C       0.0000000000       0.0000000000       0.0000000000
    N       1.1800000000       0.0000000000       0.0000000000
    H       2.1960000000       0.0000000000       0.0000000000
  End
End
System final
  Atoms
    C       0.0000000000       0.0000000000       0.0000000000
    N       1.1630000000       0.0000000000       0.0000000000
    H      -1.0780000000       0.0000000000       0.0000000000
  End
End
Task NEB

Engine DFTB
  DispersionCorrection D3-BJ
  Model DFTB3
  ResourcesDir 3ob-3-1
EndEngine
"""


class TestAMSJobWithMultipleChemicalSystems(TestAMSJobWithMultipleMolecules):
    """
    Test suite for AMSJob using multiple Chemical Systems for molecule input.
    Sets up a NEB calculation for the isomerisation of HCN.
    """

    @staticmethod
    def get_input_molecule():
        """
        Instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

        main_molecule = ChemicalSystem()
        main_molecule.add_atom("C", coords=(0, 0, 0))
        main_molecule.add_atom("N", coords=(1, 0, 0))
        main_molecule.add_atom("H", coords=(2, 0, 0))
        final_molecule = main_molecule.copy()
        final_molecule.atoms[2].coords[0] = -1
        molecule = {"": main_molecule, "final": final_molecule}

        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
NEB
  Images 9
  Iterations 100
End

Task NEB

System
   Atoms
      C 0 0 0
      N 1 0 0
      H 2 0 0
   End
End
System final
   Atoms
      C 0 0 0
      N 1 0 0
      H -1 0 0
   End
End

Engine DFTB
  DispersionCorrection D3-BJ
  Model DFTB3
  ResourcesDir 3ob-3-1
EndEngine

"""


class TestAMSJobWithMultipleChemicalSystemsAndPisa(TestAMSJobWithMultipleMoleculesAndPisa):
    """
    Test suite for AMSJob using multiple Chemical Systems for molecule input and PISA for settings input.
    Sets up a NEB calculation for the isomerisation of HCN.
    """

    @staticmethod
    def get_input_molecule():
        """
        Instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

        main_molecule = ChemicalSystem()
        main_molecule.add_atom("C", coords=(0, 0, 0))
        main_molecule.add_atom("N", coords=(1, 0, 0))
        main_molecule.add_atom("H", coords=(2, 0, 0))
        final_molecule = main_molecule.copy()
        final_molecule.atoms[2].coords[0] = -1
        molecule = {"": main_molecule, "final": final_molecule}

        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
NEB
  Images 9
  Iterations 100
End
Task NEB

Engine DFTB
  DispersionCorrection D3-BJ
  Model DFTB3
  ResourcesDir 3ob-3-1
EndEngine

System
   Atoms
      C 0 0 0
      N 1 0 0
      H 2 0 0
   End
End
System final
   Atoms
      C 0 0 0
      N 1 0 0
      H -1 0 0
   End
End"""


class TestAMSJobWithSystemBlockSettings(TestAMSJob):
    """
    Test suite for AMSJob with system block overrides/settings in the settings object.
    Sets up MD of Lennard-Jones system.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        molecule = Molecule()
        molecule.add_atom(Atom(symbol="Ar", coords=(0, 0, 0)))
        molecule.add_atom(Atom(symbol="Ar", coords=(1.605, 0.9266471820493496, 2.605)))
        molecule.lattice = [[3.21, 0.0, 0.0], [1.605, 2.779941546148048, 0.0], [0.0, 0.0, 5.21]]
        molecule.properties.charge = 42  # value to be overridden
        return molecule

    @staticmethod
    def get_input_settings():
        """
        Get instance of the Settings class passed to the AMSJob
        """
        settings = Settings()
        settings.input.ams.Task = "MolecularDynamics"
        settings.input.ams.MolecularDynamics.nSteps = 200
        settings.input.ams.MolecularDynamics.TimeStep = 5.0
        settings.input.ams.MolecularDynamics.Thermostat.Type = "NHC"
        settings.input.ams.MolecularDynamics.Thermostat.Temperature = 298.15
        settings.input.ams.MolecularDynamics.Thermostat.Tau = 100
        settings.input.ams.MolecularDynamics.Trajectory.SamplingFreq = 20
        settings.input.ams.MolecularDynamics.InitialVelocities.Type = "random"
        settings.input.ams.MolecularDynamics.InitialVelocities.Temperature = 200
        settings.input.ams.System.SuperCell = "4 4 4"
        settings.input.ams.System.PerturbCoordinates = 0.1
        settings.input.ams.System.Charge = 0

        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
MolecularDynamics
  InitialVelocities
    Temperature 200
    Type random
  End
  Thermostat
    Tau 100
    Temperature 298.15
    Type NHC
  End
  TimeStep 5.0
  Trajectory
    SamplingFreq 20
  End
  nSteps 200
End

Task MolecularDynamics

System
  Atoms
             Ar       0.0000000000       0.0000000000       0.0000000000
             Ar       1.6050000000       0.9266471820       2.6050000000
  End
  Charge 0
  Lattice
         3.2100000000     0.0000000000     0.0000000000
         1.6050000000     2.7799415461     0.0000000000
         0.0000000000     0.0000000000     5.2100000000
  End
  PerturbCoordinates 0.1
  SuperCell 4 4 4
End
"""


class TestAMSJobWithSystemBlockSettingsAndPisa(TestAMSJobWithSystemBlockSettings):
    """
    Test suite for AMSJob with system block overrides/settings in the PISA settings object.
    Sets up MD of Lennard-Jones system.
    """

    @staticmethod
    def get_input_settings():
        """
        Get instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS

        settings = Settings()
        driver = AMS()
        driver.Task = "MolecularDynamics"
        driver.MolecularDynamics.NSteps = 200
        driver.MolecularDynamics.TimeStep = 5.0
        driver.MolecularDynamics.Thermostat.Type = "NHC"
        driver.MolecularDynamics.Thermostat.Temperature = [298.15]
        driver.MolecularDynamics.Thermostat.Tau = 100
        driver.MolecularDynamics.Trajectory.SamplingFreq = 20
        driver.MolecularDynamics.InitialVelocities.Type = "random"
        driver.MolecularDynamics.InitialVelocities.Temperature = 200
        driver.System.SuperCell = [4, 4, 4]
        driver.System.PerturbCoordinates = 0.1
        driver.System.Charge = 0
        settings.input = driver

        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
MolecularDynamics
  InitialVelocities
    Temperature 200.0
    Type Random
  End
  NSteps 200
  Thermostat
    Tau 100.0
    Temperature 298.15
    Type NHC
  End
  TimeStep 5.0
  Trajectory
    SamplingFreq 20
  End
End

System
  Atoms
             Ar       0.0000000000       0.0000000000       0.0000000000
             Ar       1.6050000000       0.9266471820       2.6050000000
  End
  Charge 0.0
  Lattice
         3.2100000000     0.0000000000     0.0000000000
         1.6050000000     2.7799415461     0.0000000000
         0.0000000000     0.0000000000     5.2100000000
  End
  PerturbCoordinates 0.1
  SuperCell 4 4 4
End
Task MolecularDynamics
"""


class TestAMSJobWithSystemBlockSettingsAndPisaOnly(TestAMSJobWithSystemBlockSettings):
    """
    Test suite for AMSJob with system block overrides/settings in the PISA settings.
    Sets up MD of Lennard-Jones system.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        return None

    @staticmethod
    def get_input_settings():
        """
        Get instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS

        settings = Settings()
        driver = AMS()
        driver.Task = "MolecularDynamics"
        driver.MolecularDynamics.NSteps = 200
        driver.MolecularDynamics.TimeStep = 5.0
        driver.MolecularDynamics.Thermostat.Type = "NHC"
        driver.MolecularDynamics.Thermostat.Temperature = [298.15]
        driver.MolecularDynamics.Thermostat.Tau = 100
        driver.MolecularDynamics.Trajectory.SamplingFreq = 20
        driver.MolecularDynamics.InitialVelocities.Type = "random"
        driver.MolecularDynamics.InitialVelocities.Temperature = 200
        driver.System.Atoms = [
            "Ar 0.0000000000       0.0000000000       0.0000000000",
            "Ar 1.6050000000       0.9266471820       2.6050000000",
        ]
        driver.System.Lattice = [
            "3.2100000000     0.0000000000     0.0000000000",
            "1.6050000000     2.7799415461     0.0000000000",
            "0.0000000000     0.0000000000     5.2100000000",
        ]
        driver.System.SuperCell = [4, 4, 4]
        driver.System.PerturbCoordinates = 0.1
        driver.System.Charge = 0
        settings.input = driver

        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
MolecularDynamics
  InitialVelocities
    Temperature 200.0
    Type Random
  End
  NSteps 200
  Thermostat
    Tau 100.0
    Temperature 298.15
    Type NHC
  End
  TimeStep 5.0
  Trajectory
    SamplingFreq 20
  End
End
System
  Atoms
    Ar 0.0000000000       0.0000000000       0.0000000000
    Ar 1.6050000000       0.9266471820       2.6050000000
  End
  Charge 0.0
  Lattice
    3.2100000000     0.0000000000     0.0000000000
    1.6050000000     2.7799415461     0.0000000000
    0.0000000000     0.0000000000     5.2100000000
  End
  PerturbCoordinates 0.1
  SuperCell 4 4 4
End
Task MolecularDynamics
"""


class TestAMSJobWithSystemBlockSettingsAndChemicalSystem(TestAMSJobWithSystemBlockSettings):
    """
    Test suite for AMSJob with system block overrides/settings in the settings object and a chemical system as the molecule input.
    Sets up MD of Lennard-Jones system.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem, UnifiedLattice as Lattice

        molecule = ChemicalSystem()
        molecule.add_atom("Ar", coords=(0, 0, 0))
        molecule.add_atom("Ar", coords=(1.605, 0.9266471820493496, 2.605))
        molecule.lattice = Lattice([[3.21, 0.0, 0.0], [1.605, 2.779941546148048, 0.0], [0.0, 0.0, 5.21]])
        molecule.charge = 42  # value to be overridden
        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
MolecularDynamics
  InitialVelocities
    Temperature 200
    Type random
  End
  Thermostat
    Tau 100
    Temperature 298.15
    Type NHC
  End
  TimeStep 5.0
  Trajectory
    SamplingFreq 20
  End
  nSteps 200
End

Task MolecularDynamics

System
  Atoms
     Ar 0 0 0
     Ar 1.605 0.9266471820493496 2.605
  End
  Charge 0
  Lattice
     3.21 0 0
     1.605 2.7799415461480477 0
     0 0 5.21
  End
  PerturbCoordinates 0.1
  SuperCell 4 4 4
End
"""


class TestAMSJobWithSystemBlockSettingsAndChemicalSystemAndPisa(TestAMSJobWithSystemBlockSettingsAndPisa):
    """
    Test suite for AMSJob with system block overrides/settings in the Pisa settings object and a chemical system as the molecule input.
    Sets up MD of Lennard-Jones system.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem, UnifiedLattice as Lattice

        molecule = ChemicalSystem()
        molecule.add_atom("Ar", coords=(0, 0, 0))
        molecule.add_atom("Ar", coords=(1.605, 0.9266471820493496, 2.605))
        molecule.lattice = Lattice([[3.21, 0.0, 0.0], [1.605, 2.779941546148048, 0.0], [0.0, 0.0, 5.21]])
        molecule.charge = 42  # value to be overridden
        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
MolecularDynamics
  InitialVelocities
    Temperature 200.0
    Type Random
  End
  NSteps 200
  Thermostat
    Tau 100.0
    Temperature 298.15
    Type NHC
  End
  TimeStep 5.0
  Trajectory
    SamplingFreq 20
  End
End

System
  Atoms
     Ar 0 0 0
     Ar 1.605 0.9266471820493496 2.605
  End
  Charge 0.0
  Lattice
     3.21 0 0
     1.605 2.7799415461480477 0
     0 0 5.21
  End
  PerturbCoordinates 0.1
  SuperCell 4 4 4
End
Task MolecularDynamics
"""


class TestAMSJobWithSystemBlockSettingsAndMultipleMolecules(TestAMSJob):
    """
    Test suite for AMSJob using multiple molecules with system block overrides/settings.
    Sets up a PES scan.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        main_molecule = Molecule()
        main_molecule.add_atom(Atom(symbol="C", coords=(-0.13949723, -0.08053322, -0.12698191)))
        main_molecule.add_atom(Atom(symbol="C", coords=(0.14718597, -0.17422444, 1.25723648)))
        main_molecule.add_atom(Atom(symbol="O", coords=(1.12643603, 0.42644528, 1.83497366)))
        main_molecule.add_atom(Atom(symbol="C", coords=(-1.17969400, -0.73436273, -0.69479265)))
        main_molecule.add_atom(Atom(symbol="H", coords=(1.73538143, 1.00860600, 1.24790832)))
        main_molecule.add_atom(Atom(symbol="H", coords=(0.51558278, 0.54974352, -0.75249667)))
        main_molecule.add_atom(Atom(symbol="H", coords=(-0.42536492, -0.76321853, 2.00541072)))
        main_molecule.add_atom(Atom(symbol="H", coords=(-1.40444698, -0.66207335, -1.76629803)))
        main_molecule.add_atom(Atom(symbol="H", coords=(-1.87423948, -1.37900373, -0.15129934)))

        second_molecule = main_molecule.copy()
        second_molecule.properties.charge = 42

        molecule = {"": main_molecule, "state2": second_molecule}

        return molecule

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        settings = Settings()
        settings.input.ams.Task = "PESExploration"
        settings.input.ams.PESExploration.Job = "LandscapeRefinement"
        settings.input.ams.PESExploration.LoadEnergyLandscape.Path = "foo.results"
        settings.input.ams.PESExploration.LoadEnergyLandscape.Remove = "2 5 6 7"
        settings.input.ams.PESExploration.Optimizer.ConvergedForce = 0.01
        settings.input.ams.PESExploration.SaddleSearch.RelaxFromSaddlePoint = "T"
        settings.input.ams.GeometryOptimization.InitialHessian.Type = "Calculate"
        settings.input.ams.system = []
        settings.input.ams.system.append(Settings({"Charge": 1.0}))
        settings.input.ams.system.append(Settings({"_h": "state2", "Charge": 1.0, "PerturbCoordinates": 0.1}))
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
GeometryOptimization
  InitialHessian
    Type Calculate
  End
End

PESExploration
  Job LandscapeRefinement
  LoadEnergyLandscape
    Path foo.results
    Remove 2 5 6 7
  End
  Optimizer
    ConvergedForce 0.01
  End
  SaddleSearch
    RelaxFromSaddlePoint T
  End
End

Task PESExploration

System
  Atoms
              C      -0.1394972300      -0.0805332200      -0.1269819100
              C       0.1471859700      -0.1742244400       1.2572364800
              O       1.1264360300       0.4264452800       1.8349736600
              C      -1.1796940000      -0.7343627300      -0.6947926500
              H       1.7353814300       1.0086060000       1.2479083200
              H       0.5155827800       0.5497435200      -0.7524966700
              H      -0.4253649200      -0.7632185300       2.0054107200
              H      -1.4044469800      -0.6620733500      -1.7662980300
              H      -1.8742394800      -1.3790037300      -0.1512993400
  End
  Charge 1.0
End
System state2
  Atoms
              C      -0.1394972300      -0.0805332200      -0.1269819100
              C       0.1471859700      -0.1742244400       1.2572364800
              O       1.1264360300       0.4264452800       1.8349736600
              C      -1.1796940000      -0.7343627300      -0.6947926500
              H       1.7353814300       1.0086060000       1.2479083200
              H       0.5155827800       0.5497435200      -0.7524966700
              H      -0.4253649200      -0.7632185300       2.0054107200
              H      -1.4044469800      -0.6620733500      -1.7662980300
              H      -1.8742394800      -1.3790037300      -0.1512993400
  End
  Charge 1.0
  PerturbCoordinates 0.1
End
"""


class TestAMSJobWithSystemBlockSettingsAndMultipleMoleculesAndPisa(
    TestAMSJobWithSystemBlockSettingsAndMultipleMolecules
):
    """
    Test suite for AMSJob using multiple molecules with system block overrides/settings and PISA for settings input.
    Sets up a PES scan.
    """

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS

        settings = Settings()
        driver = AMS()
        driver.Task = "PESExploration"
        driver.PESExploration.Job = "LandscapeRefinement"
        driver.PESExploration.LoadEnergyLandscape.Path = "foo.results"
        driver.PESExploration.LoadEnergyLandscape.Remove = [2, 5, 6, 7]
        driver.PESExploration.Optimizer.ConvergedForce = 0.01
        driver.PESExploration.SaddleSearch.RelaxFromSaddlePoint = "T"
        driver.GeometryOptimization.InitialHessian.Type = "Calculate"
        driver.System[0].Charge = 1.0
        driver.System[1].Charge = 1.0
        driver.System[1].PerturbCoordinates = 0.1
        driver.System[1].header = "state2"
        settings.input = driver
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
GeometryOptimization
  InitialHessian
    Type Calculate
  End
End
PESExploration
  Job LandscapeRefinement
  LoadEnergyLandscape
    Path foo.results
    Remove 2 5 6 7
  End
  Optimizer
    ConvergedForce 0.01
  End
  SaddleSearch
    RelaxFromSaddlePoint True
  End
End

System
  Atoms
              C      -0.1394972300      -0.0805332200      -0.1269819100
              C       0.1471859700      -0.1742244400       1.2572364800
              O       1.1264360300       0.4264452800       1.8349736600
              C      -1.1796940000      -0.7343627300      -0.6947926500
              H       1.7353814300       1.0086060000       1.2479083200
              H       0.5155827800       0.5497435200      -0.7524966700
              H      -0.4253649200      -0.7632185300       2.0054107200
              H      -1.4044469800      -0.6620733500      -1.7662980300
              H      -1.8742394800      -1.3790037300      -0.1512993400
  End
  Charge 1.0
End

System state2
  Atoms
              C      -0.1394972300      -0.0805332200      -0.1269819100
              C       0.1471859700      -0.1742244400       1.2572364800
              O       1.1264360300       0.4264452800       1.8349736600
              C      -1.1796940000      -0.7343627300      -0.6947926500
              H       1.7353814300       1.0086060000       1.2479083200
              H       0.5155827800       0.5497435200      -0.7524966700
              H      -0.4253649200      -0.7632185300       2.0054107200
              H      -1.4044469800      -0.6620733500      -1.7662980300
              H      -1.8742394800      -1.3790037300      -0.1512993400
  End
  Charge 1.0
  PerturbCoordinates 0.1
End
Task PESExploration
"""


class TestAMSJobWithSystemBlockSettingsAndMultipleMoleculesAndPisaOnly(
    TestAMSJobWithSystemBlockSettingsAndMultipleMolecules
):
    """
    Test suite for AMSJob using multiple molecules with system block overrides/settings and PISA for settings and molecule input.
    Sets up a PES scan.
    """

    @staticmethod
    def get_input_molecule():
        return None

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS

        settings = Settings()
        driver = AMS()
        driver.Task = "PESExploration"
        driver.PESExploration.Job = "LandscapeRefinement"
        driver.PESExploration.LoadEnergyLandscape.Path = "foo.results"
        driver.PESExploration.LoadEnergyLandscape.Remove = [2, 5, 6, 7]
        driver.PESExploration.Optimizer.ConvergedForce = 0.01
        driver.PESExploration.SaddleSearch.RelaxFromSaddlePoint = "T"
        driver.GeometryOptimization.InitialHessian.Type = "Calculate"
        driver.System[0].Atoms = [
            "C      -0.1394972300      -0.0805332200      -0.1269819100",
            "C       0.1471859700      -0.1742244400       1.2572364800",
            "O       1.1264360300       0.4264452800       1.8349736600",
            "C      -1.1796940000      -0.7343627300      -0.6947926500",
            "H       1.7353814300       1.0086060000       1.2479083200",
            "H       0.5155827800       0.5497435200      -0.7524966700",
            "H      -0.4253649200      -0.7632185300       2.0054107200",
            "H      -1.4044469800      -0.6620733500      -1.7662980300",
            "H      -1.8742394800      -1.3790037300      -0.1512993400",
        ]
        driver.System[1].Atoms = [
            "C      -0.1394972300      -0.0805332200      -0.1269819100",
            "C       0.1471859700      -0.1742244400       1.2572364800",
            "O       1.1264360300       0.4264452800       1.8349736600",
            "C      -1.1796940000      -0.7343627300      -0.6947926500",
            "H       1.7353814300       1.0086060000       1.2479083200",
            "H       0.5155827800       0.5497435200      -0.7524966700",
            "H      -0.4253649200      -0.7632185300       2.0054107200",
            "H      -1.4044469800      -0.6620733500      -1.7662980300",
            "H      -1.8742394800      -1.3790037300      -0.1512993400",
        ]
        driver.System[0].Charge = 1.0
        driver.System[1].Charge = 1.0
        driver.System[1].PerturbCoordinates = 0.1
        driver.System[1].header = "state2"
        settings.input = driver
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
GeometryOptimization
  InitialHessian
    Type Calculate
  End
End
PESExploration
  Job LandscapeRefinement
  LoadEnergyLandscape
    Path foo.results
    Remove 2 5 6 7
  End
  Optimizer
    ConvergedForce 0.01
  End
  SaddleSearch
    RelaxFromSaddlePoint True
  End
End
System
  Atoms
    C      -0.1394972300      -0.0805332200      -0.1269819100
    C       0.1471859700      -0.1742244400       1.2572364800
    O       1.1264360300       0.4264452800       1.8349736600
    C      -1.1796940000      -0.7343627300      -0.6947926500
    H       1.7353814300       1.0086060000       1.2479083200
    H       0.5155827800       0.5497435200      -0.7524966700
    H      -0.4253649200      -0.7632185300       2.0054107200
    H      -1.4044469800      -0.6620733500      -1.7662980300
    H      -1.8742394800      -1.3790037300      -0.1512993400
  End
  Charge 1.0
End
System state2
  Atoms
    C      -0.1394972300      -0.0805332200      -0.1269819100
    C       0.1471859700      -0.1742244400       1.2572364800
    O       1.1264360300       0.4264452800       1.8349736600
    C      -1.1796940000      -0.7343627300      -0.6947926500
    H       1.7353814300       1.0086060000       1.2479083200
    H       0.5155827800       0.5497435200      -0.7524966700
    H      -0.4253649200      -0.7632185300       2.0054107200
    H      -1.4044469800      -0.6620733500      -1.7662980300
    H      -1.8742394800      -1.3790037300      -0.1512993400
  End
  Charge 1.0
  PerturbCoordinates 0.1
End
Task PESExploration
"""


class TestAMSJobWithSystemBlockSettingsAndMultipleChemicalSystems(
    TestAMSJobWithSystemBlockSettingsAndMultipleMolecules
):
    """
    Test suite for AMSJob using multiple chemical systems with system block overrides/settings.
    Sets up a PES scan.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

        main_molecule = ChemicalSystem()
        main_molecule.add_atom("C", coords=(-0.13949723, -0.08053322, -0.12698191))
        main_molecule.add_atom("C", coords=(0.14718597, -0.17422444, 1.25723648))
        main_molecule.add_atom("O", coords=(1.12643603, 0.42644528, 1.83497366))
        main_molecule.add_atom("C", coords=(-1.17969400, -0.73436273, -0.69479265))
        main_molecule.add_atom("H", coords=(1.73538143, 1.00860600, 1.24790832))
        main_molecule.add_atom("H", coords=(0.51558278, 0.54974352, -0.75249667))
        main_molecule.add_atom("H", coords=(-0.42536492, -0.76321853, 2.00541072))
        main_molecule.add_atom("H", coords=(-1.40444698, -0.66207335, -1.76629803))
        main_molecule.add_atom("H", coords=(-1.87423948, -1.37900373, -0.15129934))

        second_molecule = main_molecule.copy()
        second_molecule.charge = 42

        molecule = {"": main_molecule, "state2": second_molecule}

        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
GeometryOptimization
  InitialHessian
    Type Calculate
  End
End

PESExploration
  Job LandscapeRefinement
  LoadEnergyLandscape
    Path foo.results
    Remove 2 5 6 7
  End
  Optimizer
    ConvergedForce 0.01
  End
  SaddleSearch
    RelaxFromSaddlePoint T
  End
End

Task PESExploration

System
  Atoms
     C -0.13949722999999997 -0.08053322 -0.12698191
     C 0.14718597 -0.17422443999999998 1.25723648
     O 1.1264360299999998 0.42644528 1.83497366
     C -1.1796939999999998 -0.73436273 -0.69479265
     H 1.73538143 1.008606 1.2479083199999998
     H 0.51558278 0.54974352 -0.75249667
     H -0.42536492 -0.7632185299999998 2.00541072
     H -1.4044469799999997 -0.66207335 -1.7662980299999997
     H -1.87423948 -1.37900373 -0.15129933999999998
  End
  Charge 1.0
End
System state2
  Atoms
     C -0.13949722999999997 -0.08053322 -0.12698191
     C 0.14718597 -0.17422443999999998 1.25723648
     O 1.1264360299999998 0.42644528 1.83497366
     C -1.1796939999999998 -0.73436273 -0.69479265
     H 1.73538143 1.008606 1.2479083199999998
     H 0.51558278 0.54974352 -0.75249667
     H -0.42536492 -0.7632185299999998 2.00541072
     H -1.4044469799999997 -0.66207335 -1.7662980299999997
     H -1.87423948 -1.37900373 -0.15129933999999998
  End
  Charge 1.0
  PerturbCoordinates 0.1
End
"""


class TestAMSJobWithSystemBlockSettingsAndMultipleChemicalSystemsAndPisa(
    TestAMSJobWithSystemBlockSettingsAndMultipleMolecules
):
    """
    Test suite for AMSJob using multiple chemical systems with system block overrides/settings and PISA for settings input.
    Sets up a PES scan.
    """

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS

        settings = Settings()
        driver = AMS()
        driver.Task = "PESExploration"
        driver.PESExploration.Job = "LandscapeRefinement"
        driver.PESExploration.LoadEnergyLandscape.Path = "foo.results"
        driver.PESExploration.LoadEnergyLandscape.Remove = [2, 5, 6, 7]
        driver.PESExploration.Optimizer.ConvergedForce = 0.01
        driver.PESExploration.SaddleSearch.RelaxFromSaddlePoint = "T"
        driver.GeometryOptimization.InitialHessian.Type = "Calculate"
        driver.System[0].Charge = 1.0
        driver.System[1].Charge = 1.0
        driver.System[1].PerturbCoordinates = 0.1
        driver.System[1].header = "state2"
        settings.input = driver
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
GeometryOptimization
  InitialHessian
    Type Calculate
  End
End
PESExploration
  Job LandscapeRefinement
  LoadEnergyLandscape
    Path foo.results
    Remove 2 5 6 7
  End
  Optimizer
    ConvergedForce 0.01
  End
  SaddleSearch
    RelaxFromSaddlePoint True
  End
End

System
  Atoms
              C      -0.1394972300      -0.0805332200      -0.1269819100
              C       0.1471859700      -0.1742244400       1.2572364800
              O       1.1264360300       0.4264452800       1.8349736600
              C      -1.1796940000      -0.7343627300      -0.6947926500
              H       1.7353814300       1.0086060000       1.2479083200
              H       0.5155827800       0.5497435200      -0.7524966700
              H      -0.4253649200      -0.7632185300       2.0054107200
              H      -1.4044469800      -0.6620733500      -1.7662980300
              H      -1.8742394800      -1.3790037300      -0.1512993400
  End
  Charge 1.0
End

System state2
  Atoms
              C      -0.1394972300      -0.0805332200      -0.1269819100
              C       0.1471859700      -0.1742244400       1.2572364800
              O       1.1264360300       0.4264452800       1.8349736600
              C      -1.1796940000      -0.7343627300      -0.6947926500
              H       1.7353814300       1.0086060000       1.2479083200
              H       0.5155827800       0.5497435200      -0.7524966700
              H      -0.4253649200      -0.7632185300       2.0054107200
              H      -1.4044469800      -0.6620733500      -1.7662980300
              H      -1.8742394800      -1.3790037300      -0.1512993400
  End
  Charge 1.0
  PerturbCoordinates 0.1
End
Task PESExploration
"""


class TestAMSJobWithChainOfMolecules(TestAMSJob):
    """
    Test suite for AMSJob with a chain of water molecules.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        mol = TestAMSJob.get_input_molecule()
        mol.lattice = [[3, 0, 0]]
        return mol.supercell(4)

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        settings = Settings()
        settings.input.Mopac.SCF.ConvergenceThreshold = 1.0e-8
        settings.input.Mopac.model = "pm6"
        settings.input.AMS.Task = "SinglePoint"
        settings.input.AMS.Properties.Gradients = "Yes"
        settings.input.AMS.NumericalDifferentiation.Parallel.nCoresPerGroup = 1
        settings.input.AMS.NumericalDifferentiation.NuclearStepSize = 0.0001
        settings.input.AMS.EngineDebugging.IgnoreGradientsRequest = "No"
        settings.input.AMS.System.ElectrostaticEmbedding.ElectricField = "0.0 0.0 0.0"
        settings.input.AMS.Task = "SinglePoint"
        settings.input.AMS.Properties.Gradients = "Yes"
        settings.input.AMS.NumericalDifferentiation.Parallel.nCoresPerGroup = 1
        settings.input.AMS.NumericalDifferentiation.NuclearStepSize = 0.0001
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
EngineDebugging
  IgnoreGradientsRequest No
End

NumericalDifferentiation
  NuclearStepSize 0.0001
  Parallel
    nCoresPerGroup 1
  End
End

Properties
  Gradients Yes
End

Task SinglePoint

System
  Atoms
              O       0.0000000000       0.0000000000       0.0000000000
              H       1.0000000000       0.0000000000       0.0000000000
              H       0.0000000000       1.0000000000       0.0000000000
              O       3.0000000000       0.0000000000       0.0000000000
              H       4.0000000000       0.0000000000       0.0000000000
              H       3.0000000000       1.0000000000       0.0000000000
              O       6.0000000000       0.0000000000       0.0000000000
              H       7.0000000000       0.0000000000       0.0000000000
              H       6.0000000000       1.0000000000       0.0000000000
              O       9.0000000000       0.0000000000       0.0000000000
              H      10.0000000000       0.0000000000       0.0000000000
              H       9.0000000000       1.0000000000       0.0000000000
  End
  ElectrostaticEmbedding
    ElectricField 0.0 0.0 0.0
  End
  Lattice
        12.0000000000     0.0000000000     0.0000000000
  End
End

Engine Mopac
  SCF
    ConvergenceThreshold 1e-08
  End
  model pm6
EndEngine

"""


class TestAMSJobRun:

    def test_run_with_watch_forwards_ams_logs_to_stdout(self, config):
        # Patch the config and the stdout
        from scm.plams.core.logging import get_logger

        with patch("scm.plams.interfaces.adfsuite.ams.config", config), patch(
            "sys.stdout", new_callable=StringIO
        ) as mock_stdout:
            logger = get_logger("test_run_with_watch_forwards_ams_logs_to_stdout")
            with patch("scm.plams.core.functions._logger", logger):
                config.log.date = False
                config.log.time = False

                # Given a dummy job
                job = AMSJob()

                # Which writes logs to logfile periodically on background thread
                logfile = os.path.join(config.default_jobmanager.workdir, job.name, "ams.log")

                def write_logs():
                    for i in range(10):
                        time.sleep(0.01)
                        try:
                            with open(logfile, "a") as f:
                                f.write(f"<Oct21-2024> <09:34:54>  line {i}\n")
                                f.flush()
                        except FileNotFoundError:
                            pass

                background_thread = threading.Thread(target=write_logs)

                def get_runscript() -> str:
                    background_thread.start()
                    return "sleep 1"

                job.get_runscript = get_runscript

                # When run job and watching output
                job.run(watch=True)

                # Then ams logs are also forwarded to the standard output
                stdout = mock_stdout.getvalue().replace("\r\n", "\n").split("\n")
                postrun_lines = [l for l in stdout if re.fullmatch("plamsjob: line \d", l)]
                status_lines = [l for l in stdout if re.fullmatch("JOB plamsjob (STARTED|RUNNING|FINISHED|FAILED)", l)]

                assert len(postrun_lines) == 10
                assert len(status_lines) == 4

                logger.close()

    def test_get_errormsg_populates_correctly_for_different_scenarios(self):
        from scm.plams.core.errors import FileError

        # Invalid license
        job = AMSJob()
        job._error_msg = None
        results = MagicMock(spec=AMSResults)
        results.readrkf.side_effect = FileError()
        results.grep_file.side_effect = FileError()
        results.get_output_chunk.return_value = [
            "LICENSE INVALID",
            "---------------",
            "Your license does not include module AMS version 2024.206 on this machine.",
            "Module AMS",
            "Version 2024.206",
            "Machine: Foo",
            "License file: ./license.txt",
        ]
        job.results = results

        assert (
            job.get_errormsg()
            == """LICENSE INVALID
---------------
Your license does not include module AMS version 2024.206 on this machine.
Module AMS
Version 2024.206
Machine: Foo
License file: ./license.txt"""
        )

        # Invalid input
        job._error_msg = None
        results.grep_file.side_effect = None
        results.grep_file.return_value = [
            '<Dec05-2024> <12:03:49>  ERROR: Input error: value "foo" found in line 13 for multiple choice key "Task" is not an allowed choice'
        ]
        assert (
            job.get_errormsg()
            == 'Input error: value "foo" found in line 13 for multiple choice key "Task" is not an allowed choice'
        )

        # Error in calculation
        job._error_msg = None
        results.readrkf.return_value = "NORMAL TERMINATION with errors"
        results.readrkf.side_effect = None
        results.grep_file.return_value = [
            "<Dec05-2024> <13:44:55>  ERROR: Geometry optimization failed! (Did not converge.)"
        ]
        assert job.get_errormsg() == "Geometry optimization failed! (Did not converge.)"

        # Error in prerun
        job._error_msg = "RuntimeError: something went wrong"
        assert job.get_errormsg() == "RuntimeError: something went wrong"

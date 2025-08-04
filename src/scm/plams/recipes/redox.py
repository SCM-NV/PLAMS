import os
import shutil
from collections import OrderedDict

from scm.plams.core.basejob import MultiJob
from scm.plams.core.functions import add_to_instance
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.crs import CRSJob, CRSResults
from scm.plams.mol.molecule import Molecule
from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob
from scm.plams.tools.units import Units

"""
    Classes for calculating oxidation potentials and reduction potentials in implicit solvent.

    To get the results, call

    job.results.get_oxidation_potential(vibrations=False)

    or

    job.results.get_reduction_potential(vibrations=False)

    Set ``vibrations=True`` to include vibrational effects if the job was run
    with ``vibrations=True``.

    The potential is returned in V. It is on an absolute scale. To get
    potential relative to SHE, subtract 4.47 V.

    AMSRedoxDirectJob
    --------------------

    Calculates reaction energies for 
    A + e^- --> A^-     if ``reduction=True``,
    and/or
    A --> A^+ + e^-     if ``oxidation=True``.

    using direct geometry optimizations in implicit solvent. The solvent must be one supported by ADF.

    If ``vibrations`` = True, then calculate and use Gibbs free energies at room temperature.

    ``settings`` should be set to the engine settings, excluding spin polarization and implicit solvation.

    Requires an ADF license.

    AMSRedoxThermodynamicCycleJob
    -----------------------------

    Sometimes more efficient (and less accurate) alternative to AMSRedoxDirectJob. Can be
    useful if the molecule is large and if ``vibrations=True``. For small molecules just use AMSRedoxDirectJob.

    Requires an ADF license.

    AMSRedoxScreeningJob
    -----------------------

    The fastest and least accurate option. Geometries are optimized at the DFTB
    level, ADF is then used to generate a .coskf file and the solvation free
    energies are evaluated with COSMO-RS.

    Note: you must supply ``solvent_coskf``, a path to the solvent .coskf file.
    If you do not have one, you can generate one with the
    ``ADFCOSMORSCompoundJob`` recipe.

    Requires DFTB, ADF, and COSMO-RS licenses.

    Note: The functional used is different from Belic et al. This class uses the
    settings that are meant to be used for generating .coskf files.

"""

__all__ = [
    "AMSRedoxScreeningJob",
    "AMSRedoxScreeningResults",
    "AMSRedoxDirectJob",
    "AMSRedoxDirectResults",
    "AMSRedoxThermodynamicCycleJob",
    "AMSRedoxThermodynamicCycleResults",
]


def _spinpol_settings(molecule, charge=0):
    sett = Settings()
    num_electrons = sum(atom.atnum for atom in molecule) - charge
    spinpol = num_electrons % 2
    if spinpol == 0:
        unrestricted = "No"
    else:
        unrestricted = "Yes"
        # sett.input.adf.OCCUPATIONS = "ElectronicTemperature=100"
    sett.input.adf.SpinPolarization = spinpol
    sett.input.adf.Unrestricted = unrestricted
    return sett


class AMSRedoxParentJob(MultiJob):
    def __init__(
        self,
        molecule: Molecule,  # Molecule
        name: str = None,
        settings: Settings = None,
        oxidation: bool = True,
        reduction: bool = False,
    ):

        MultiJob.__init__(self, children=OrderedDict(), name=name)

        self.oxidation = oxidation
        self.reduction = reduction
        self.input_molecule = molecule
        self.settings = settings or Settings()
        self.orig_charge = molecule.properties.get("charge", 0)
        self.ox_charge = self.orig_charge + 1
        self.red_charge = self.orig_charge - 1

    @staticmethod
    def solvation_settings(solvent):
        sett = Settings()
        if solvent is None:
            solvent = "vacuum"
        if isinstance(solvent, str):
            sett.input.adf.Solvation.Solv = f"name={solvent}"
        elif isinstance(solvent, tuple):
            sett.input.adf.Solvation.Solv = f"eps={solvent[0]} rad={solvent[1]}"
        else:
            raise TypeError(f"solvent is of type {type(solvent)}, expected str or 2-tuple")

        return sett

    def get_dft_0_settings(self, vibrations=False, perturbcoordinates=False):
        """all settings except task and solvation"""
        s = Settings()
        s += self.settings
        s.update(_spinpol_settings(self.input_molecule, self.orig_charge))
        s.input.ams.System.Charge = self.orig_charge
        s.input.ams.System.PerturbCoordinates = 0.01 if perturbcoordinates else 0.00
        s.input.ams.Properties.NormalModes = str(vibrations)
        return s

    def get_dft_ox_settings(self, vibrations=False, perturbcoordinates=False):
        sox = self.get_dft_0_settings(vibrations=vibrations, perturbcoordinates=perturbcoordinates)
        sox.update(_spinpol_settings(self.input_molecule, self.ox_charge))
        sox.input.ams.System.Charge = self.ox_charge
        return sox

    def get_dft_red_settings(self, vibrations=False, perturbcoordinates=False):
        sred = self.get_dft_0_settings(vibrations=vibrations, perturbcoordinates=perturbcoordinates)
        sred.update(_spinpol_settings(self.input_molecule, self.red_charge))
        sred.input.ams.System.Charge = self.red_charge
        return sred


class AMSRedoxParentResults(Results):
    Gelectron = -0.0375 / 27.211  # in hartree

    @staticmethod
    def _get_energy(job, vibrations: bool):
        if not vibrations:
            return job.results.get_energy()
        e = job.results.readrkf("Thermodynamics", "Gibbs free Energy", file="engine")
        if isinstance(e, list):
            return e[0]
        return e

    @staticmethod
    def _get_solvation_energy(job):
        dG_solvation = job.results.readrkf("Energy", "Solvation Energy (el)", "adf") + job.results.readrkf(
            "Energy", "Solvation Energy (cd)", "adf"
        )
        return dG_solvation


class CRSActivityCoefficientResults(CRSResults):
    def get_gibbs_energy(self):
        """Gibbs energy in hartree"""
        return float(self.get_results()["G solute"][1]) * 0.00159376  # kcal/mol to hartree


class CRSActivityCoefficientJob(CRSJob):
    _result_type = CRSActivityCoefficientResults

    def __init__(self, solvent_coskf, solute_coskf, name=None, temperature=298.15, copy_coskf=False):
        CRSJob.__init__(self, name=name)
        self.solvent_coskf = solvent_coskf
        self.solute_coskf = solute_coskf
        self.temperature = temperature
        self.copy_coskf = copy_coskf

    def prerun(self):  # noqa F811
        self._prerun()

    def _prerun(self):
        if self.copy_coskf:
            new_solvent = f"solvent_{os.path.basename(self.solvent_coskf)}"
            new_solute = f"solute_{os.path.basename(self.solute_coskf)}"
            shutil.copy(self.solvent_coskf, os.path.join(self.path, new_solvent))
            shutil.copy(self.solute_coskf, os.path.join(self.path, new_solute))
            self.solvent_coskf = new_solvent
            self.solute_coskf = new_solute
        self.settings.input.property._h = "ACTIVITYCOEF"
        compounds = [Settings(), Settings()]
        compounds[0]._h = self.solvent_coskf
        compounds[1]._h = self.solute_coskf
        compounds[0].frac1 = 1
        compounds[1].frac1 = 0

        self.settings.input.temperature = str(self.temperature)
        self.settings.input.compound = compounds


class AMSRedoxDirectResults(AMSRedoxParentResults):
    def get_oxidation_potential(self, vibrations=True, unit="eV"):
        Greact = self._get_energy(self.job.children["job_0"], vibrations=vibrations)
        Gprod = self._get_energy(self.job.children["job_ox"], vibrations=vibrations) + self.Gelectron
        ret = (Gprod - Greact) * Units.convert(1.0, "hartree", unit)
        return ret

    def get_reduction_potential(self, vibrations=True, unit="eV"):
        Greact = self._get_energy(self.job.children["job_0"], vibrations=vibrations) + self.Gelectron
        Gprod = self._get_energy(self.job.children["job_red"], vibrations=vibrations)
        ret = (Gprod - Greact) * Units.convert(1.0, "hartree", unit)
        ret *= -1  # deltaG = -nFE
        return ret


class AMSRedoxDirectJob(AMSRedoxParentJob):
    _result_type = AMSRedoxDirectResults

    def __init__(
        self,
        molecule: Molecule,  # Molecule
        solvent: str,  # ADF pre-defined solvent
        name: str = None,
        settings: Settings = None,
        oxidation: bool = True,
        reduction: bool = False,
        vibrations: bool = True,
    ):

        AMSRedoxParentJob.__init__(
            self, name=name, settings=settings, molecule=molecule, oxidation=oxidation, reduction=reduction
        )

        self.solvent = solvent
        self.vibrations = vibrations

        s = Settings()
        s += self.settings
        s.input.ams.Task = "GeometryOptimization"
        s += self.solvation_settings(self.solvent)
        s.update(_spinpol_settings(self.input_molecule, self.orig_charge))
        s.input.ams.System.Charge = self.orig_charge
        s.input.ams.System.PerturbCoordinates = 0.01
        s.input.ams.Properties.NormalModes = str(self.vibrations)
        job_0 = AMSJob(settings=s, name="job_0", molecule=self.input_molecule)
        self.children["job_0"] = job_0

        if oxidation:
            sox = s.copy()
            sox.update(_spinpol_settings(self.input_molecule, self.ox_charge))
            sox.input.ams.System.Charge = self.ox_charge
            sox.input.ams.System.PerturbCoordinates = 0.01
            job_ox = AMSJob(settings=sox, name="job_ox")

            @add_to_instance(job_ox)
            def prerun(self):  # noqa F811
                self.molecule = job_0.results.get_main_molecule()

            self.children["job_ox"] = job_ox

        if reduction:
            sred = s.copy()
            sred.update(_spinpol_settings(self.input_molecule, self.red_charge))
            sred.input.ams.System.Charge = self.red_charge
            sred.input.ams.System.PerturbCoordinates = 0.01
            job_red = AMSJob(settings=sred, name="job_red")

            @add_to_instance(job_red)
            def prerun(self):  # noqa F811
                self.molecule = job_0.results.get_main_molecule()

            self.children["job_red"] = job_red


class AMSRedoxScreeningResults(Results):
    Gelectron = -0.0375 / 27.211  # in hartree

    def get_oxidation_potential(self, vibrations=True, unit="eV"):
        """Note: vibrations cannot be disabled for AMSRedoxScreening"""
        Greact = self.job.children["activitycoef_0"].results.get_gibbs_energy()
        Gprod = self.job.children["activitycoef_ox"].results.get_gibbs_energy() + self.Gelectron
        ret = (Gprod - Greact) * Units.convert(1.0, "hartree", unit)
        return ret

    def get_reduction_potential(self, vibrations=True, unit="eV"):
        """Note: vibrations cannot be disabled for AMSRedoxScreening"""
        Greact = self.job.children["activitycoef_0"].results.get_gibbs_energy() + self.Gelectron
        Gprod = self.job.children["activitycoef_red"].results.get_gibbs_energy()
        ret = (Gprod - Greact) * Units.convert(1.0, "hartree", unit)
        ret *= -1  # deltaG = -nFE
        return ret


class AMSRedoxScreeningJob(AMSRedoxParentJob):
    _result_type = AMSRedoxScreeningResults

    def __init__(
        self,
        molecule: Molecule,  # Molecule
        solvent_coskf: str,  # path to solvent .coskf file
        name: str = None,
        settings: Settings = None,
        oxidation: bool = True,
        reduction: bool = False,
        copy_coskf=False,
    ):
        """
        copy_coskf: bool
            Will copy the coskf files into the job's own directory before running.
        """
        AMSRedoxParentJob.__init__(
            self, name=name, settings=settings, molecule=molecule, oxidation=oxidation, reduction=reduction
        )

        self.solvent_coskf = solvent_coskf

        # dftb_go_0
        s = Settings()
        s.input.ams.Task = "GeometryOptimization"
        s.input.DFTB.Model = "GFN1-xTB"
        s.input.ams.System.Charge = self.orig_charge
        s.input.ams.System.PerturbCoordinates = 0.01
        dftb_go_0 = AMSJob(settings=s, molecule=self.input_molecule, name="dftb_go_0")
        self.children["dftb_go_0"] = dftb_go_0

        if oxidation:
            sox = s.copy()
            sox.input.ams.System.PerturbCoordinates = 0.01
            sox.input.ams.System.Charge = self.ox_charge
            dftb_go_ox = AMSJob(settings=sox, name="dftb_go_ox")

            @add_to_instance(dftb_go_ox)
            def prerun(self):  # noqa F811
                self.molecule = dftb_go_0.results.get_main_molecule()

            self.children["dftb_go_ox"] = dftb_go_ox

        if reduction:
            sred = s.copy()
            sred.input.ams.System.PerturbCoordinates = 0.01
            sred.input.ams.System.Charge = self.red_charge
            dftb_go_red = AMSJob(settings=sred, name="dftb_go_red")

            @add_to_instance(dftb_go_red)
            def prerun(self):  # noqa F811
                self.molecule = dftb_go_0.results.get_main_molecule()

            self.children["dftb_go_red"] = dftb_go_red

        # gencoskf_0
        s = self.settings.copy()
        s.update(_spinpol_settings(self.input_molecule, self.orig_charge))
        # s.input.ams.System.Charge = self.orig_charge
        gencoskf_0 = ADFCOSMORSCompoundJob(molecule=None, name="gencoskf_0", singlepoint=True, settings=s)

        @add_to_instance(gencoskf_0)
        def prerun(self):  # noqa F811
            self.input_molecule = dftb_go_0.results.get_main_molecule()  # will inherit the charge
            assert isinstance(self.input_molecule, Molecule)
            self.atomic_ion = len(self.input_molecule) == 1

        self.children["gencoskf_0"] = gencoskf_0

        if oxidation:
            sox = s.copy()
            sox.update(_spinpol_settings(self.input_molecule, self.ox_charge))
            # sox.input.ams.System.Charge = self.ox_charge
            gencoskf_ox = ADFCOSMORSCompoundJob(molecule=None, name="gencoskf_ox", singlepoint=True, settings=sox)

            @add_to_instance(gencoskf_ox)
            def prerun(self):  # noqa F811
                self.input_molecule = dftb_go_ox.results.get_main_molecule()
                self.atomic_ion = len(self.input_molecule) == 1

            self.children["gencoskf_ox"] = gencoskf_ox

        if reduction:
            sred = s.copy()
            sred.update(_spinpol_settings(self.input_molecule, self.red_charge))
            # sred.input.ams.System.Charge = self.red_charge
            gencoskf_red = ADFCOSMORSCompoundJob(molecule=None, name="gencoskf_red", singlepoint=True, settings=sred)

            @add_to_instance(gencoskf_red)
            def prerun(self):  # noqa F811
                self.input_molecule = dftb_go_red.results.get_main_molecule()
                self.atomic_ion = len(self.input_molecule) == 1

            self.children["gencoskf_red"] = gencoskf_red

        # activitycoef_0
        activitycoef_0 = CRSActivityCoefficientJob(
            name="activitycoef_0", solvent_coskf=self.solvent_coskf, solute_coskf=None, copy_coskf=copy_coskf
        )

        @add_to_instance(activitycoef_0)
        def prerun(self):  # noqa F811
            self.solute_coskf = gencoskf_0.results.coskfpath()
            self._prerun()

        self.children["activitycoef_0"] = activitycoef_0

        if oxidation:
            activitycoef_ox = CRSActivityCoefficientJob(
                name="activitycoef_ox", solvent_coskf=self.solvent_coskf, solute_coskf=None, copy_coskf=copy_coskf
            )

            @add_to_instance(activitycoef_ox)
            def prerun(self):  # noqa F811
                self.solute_coskf = gencoskf_ox.results.coskfpath()
                self._prerun()

            self.children["activitycoef_ox"] = activitycoef_ox

        if reduction:
            activitycoef_red = CRSActivityCoefficientJob(
                name="activitycoef_red", solvent_coskf=self.solvent_coskf, solute_coskf=None, copy_coskf=copy_coskf
            )

            @add_to_instance(activitycoef_red)
            def prerun(self):  # noqa F811
                self.solute_coskf = gencoskf_red.results.coskfpath()
                self._prerun()

            self.children["activitycoef_red"] = activitycoef_red


class AMSRedoxThermodynamicCycleResults(AMSRedoxParentResults):
    def get_oxidation_potential(self, vibrations=True, unit="eV"):
        Greact = (
            self._get_energy(self.job.children["go_0_vacuum"], vibrations=vibrations)
            + self._get_energy(self.job.children["go_0_vacuum_sp_solvated"], vibrations=False)
            + self._get_energy(self.job.children["go_0_solvated_sp_vacuum"], vibrations=False)
            - 2 * self._get_energy(self.job.children["go_0_vacuum"], vibrations=False)
        )

        Gprod = (
            self._get_energy(self.job.children["go_ox_vacuum"], vibrations=vibrations)
            + self._get_energy(self.job.children["go_ox_vacuum_sp_solvated"], vibrations=False)
            + self._get_energy(self.job.children["go_ox_solvated_sp_vacuum"], vibrations=False)
            - 2 * self._get_energy(self.job.children["go_ox_vacuum"], vibrations=False)
            + self.Gelectron
        )

        ret = (Gprod - Greact) * Units.convert(1.0, "hartree", unit)
        return ret

    def get_reduction_potential(self, vibrations=True, unit="eV"):
        Greact = (
            self._get_energy(self.job.children["go_0_vacuum"], vibrations=vibrations)
            + self._get_energy(self.job.children["go_0_vacuum_sp_solvated"], vibrations=False)
            + self._get_energy(self.job.children["go_0_solvated_sp_vacuum"], vibrations=False)
            - 2 * self._get_energy(self.job.children["go_0_vacuum"], vibrations=False)
            + self.Gelectron
        )

        Gprod = (
            self._get_energy(self.job.children["go_red_vacuum"], vibrations=vibrations)
            + self._get_energy(self.job.children["go_red_vacuum_sp_solvated"], vibrations=False)
            + self._get_energy(self.job.children["go_red_solvated_sp_vacuum"], vibrations=False)
            - 2 * self._get_energy(self.job.children["go_red_vacuum"], vibrations=False)
        )

        ret = (Gprod - Greact) * Units.convert(1.0, "hartree", unit)
        ret *= -1  # deltaG = -nFE
        return ret


class AMSRedoxThermodynamicCycleJob(AMSRedoxParentJob):
    _result_type = AMSRedoxThermodynamicCycleResults

    def __init__(
        self,
        molecule: Molecule,  # Molecule
        name: str = None,
        solvent: str = "Water",
        settings: Settings = None,
        oxidation: bool = True,
        reduction: bool = False,
        vibrations: bool = False,
    ):

        AMSRedoxParentJob.__init__(
            self, name=name, settings=settings, molecule=molecule, oxidation=oxidation, reduction=reduction
        )

        self.vibrations = vibrations
        self.solvent = solvent

        # go_0_vacuum
        s = self.get_dft_0_settings(vibrations=self.vibrations, perturbcoordinates=True)
        s.input.ams.Task = "GeometryOptimization"
        go_0_vacuum = AMSJob(settings=s, molecule=self.input_molecule, name="go_0_vacuum")

        @add_to_instance(go_0_vacuum)
        def prerun(self):  # noqa F811
            self.molecule = self.parent.input_molecule

        self.children["go_0_vacuum"] = go_0_vacuum

        # go_0_vacuum_sp_solvated
        ss = self.get_dft_0_settings(vibrations=False, perturbcoordinates=False)
        ss.input.ams.Task = "SinglePoint"
        ss += self.solvation_settings(self.solvent)
        go_0_vacuum_sp_solvated = AMSJob(settings=ss, molecule=None, name="go_0_vacuum_sp_solvated")

        @add_to_instance(go_0_vacuum_sp_solvated)
        def prerun(self):  # noqa F811
            self.molecule = go_0_vacuum.results.get_main_molecule()

        self.children["go_0_vacuum_sp_solvated"] = go_0_vacuum_sp_solvated

        # go_0_solvated
        s = self.get_dft_0_settings(vibrations=self.vibrations, perturbcoordinates=True)
        s.input.ams.Task = "GeometryOptimization"
        s += self.solvation_settings(self.solvent)
        go_0_solvated = AMSJob(settings=s, molecule=None, name="go_0_solvated")

        @add_to_instance(go_0_solvated)
        def prerun(self):  # noqa F811
            self.molecule = go_0_vacuum.results.get_main_molecule()

        self.children["go_0_solvated"] = go_0_solvated

        # go_0_solvated_sp_vacuum
        ss = self.get_dft_0_settings(vibrations=False, perturbcoordinates=False)
        ss.input.ams.Task = "SinglePoint"
        go_0_solvated_sp_vacuum = AMSJob(settings=ss, molecule=None, name="go_0_solvated_sp_vacuum")

        @add_to_instance(go_0_solvated_sp_vacuum)
        def prerun(self):  # noqa F811
            self.molecule = go_0_solvated.results.get_main_molecule()

        self.children["go_0_solvated_sp_vacuum"] = go_0_solvated_sp_vacuum

        if oxidation:
            # go_ox_vacuum
            s = self.get_dft_ox_settings(vibrations=self.vibrations, perturbcoordinates=True)
            s.input.ams.Task = "GeometryOptimization"
            go_ox_vacuum = AMSJob(settings=s, molecule=self.input_molecule, name="go_ox_vacuum")

            @add_to_instance(go_ox_vacuum)
            def prerun(self):  # noqa F811
                self.molecule = go_0_vacuum.results.get_main_molecule()

            self.children["go_ox_vacuum"] = go_ox_vacuum

            # go_ox_vacuum_sp_solvated
            ss = self.get_dft_ox_settings(vibrations=False, perturbcoordinates=False)
            ss.input.ams.Task = "SinglePoint"
            ss += self.solvation_settings(self.solvent)
            go_ox_vacuum_sp_solvated = AMSJob(settings=ss, molecule=None, name="go_ox_vacuum_sp_solvated")

            @add_to_instance(go_ox_vacuum_sp_solvated)
            def prerun(self):  # noqa F811
                self.molecule = go_ox_vacuum.results.get_main_molecule()

            self.children["go_ox_vacuum_sp_solvated"] = go_ox_vacuum_sp_solvated

            # go_ox_solvated
            s = self.get_dft_ox_settings(vibrations=self.vibrations, perturbcoordinates=True)
            s.input.ams.Task = "GeometryOptimization"
            s += self.solvation_settings(self.solvent)
            go_ox_solvated = AMSJob(settings=s, molecule=None, name="go_ox_solvated")

            @add_to_instance(go_ox_solvated)
            def prerun(self):  # noqa F811
                self.molecule = go_ox_vacuum.results.get_main_molecule()

            self.children["go_ox_solvated"] = go_ox_solvated

            # go_ox_solvated_sp_vacuum
            ss = self.get_dft_ox_settings(vibrations=False, perturbcoordinates=False)
            ss.input.ams.Task = "SinglePoint"
            go_ox_solvated_sp_vacuum = AMSJob(settings=ss, molecule=None, name="go_ox_solvated_sp_vacuum")

            @add_to_instance(go_ox_solvated_sp_vacuum)
            def prerun(self):  # noqa F811
                self.molecule = go_ox_solvated.results.get_main_molecule()

            self.children["go_ox_solvated_sp_vacuum"] = go_ox_solvated_sp_vacuum

        if reduction:
            # go_red_vacuum
            s = self.get_dft_red_settings(vibrations=self.vibrations, perturbcoordinates=True)
            s.input.ams.Task = "GeometryOptimization"
            go_red_vacuum = AMSJob(settings=s, molecule=self.input_molecule, name="go_red_vacuum")

            @add_to_instance(go_red_vacuum)
            def prerun(self):  # noqa F811
                self.molecule = go_0_vacuum.results.get_main_molecule()

            self.children["go_red_vacuum"] = go_red_vacuum

            # go_red_vacuum_sp_solvated
            ss = self.get_dft_red_settings(vibrations=False, perturbcoordinates=False)
            ss.input.ams.Task = "SinglePoint"
            ss += self.solvation_settings(self.solvent)
            go_red_vacuum_sp_solvated = AMSJob(settings=ss, molecule=None, name="go_red_vacuum_sp_solvated")

            @add_to_instance(go_red_vacuum_sp_solvated)
            def prerun(self):  # noqa F811
                self.molecule = go_red_vacuum.results.get_main_molecule()

            self.children["go_red_vacuum_sp_solvated"] = go_red_vacuum_sp_solvated

            # go_red_solvated
            s = self.get_dft_red_settings(vibrations=self.vibrations, perturbcoordinates=True)
            s.input.ams.Task = "GeometryOptimization"
            s += self.solvation_settings(self.solvent)
            go_red_solvated = AMSJob(settings=s, molecule=None, name="go_red_solvated")

            @add_to_instance(go_red_solvated)
            def prerun(self):  # noqa F811
                self.molecule = go_red_vacuum.results.get_main_molecule()

            self.children["go_red_solvated"] = go_red_solvated

            # go_red_solvated_sp_vacuum
            ss = self.get_dft_red_settings(vibrations=False, perturbcoordinates=False)
            ss.input.ams.Task = "SinglePoint"
            go_red_solvated_sp_vacuum = AMSJob(settings=ss, molecule=None, name="go_red_solvated_sp_vacuum")

            @add_to_instance(go_red_solvated_sp_vacuum)
            def prerun(self):  # noqa F811
                self.molecule = go_red_solvated.results.get_main_molecule()

            self.children["go_red_solvated_sp_vacuum"] = go_red_solvated_sp_vacuum

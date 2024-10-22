from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.tools.units import Units
from scm.plams.core.errors import FileError
from scm.plams.mol.molecule import Molecule


from os.path import join as opj
from os.path import abspath, isdir, basename
from typing import Dict

__all__ = ["ADFFragmentJob", "ADFFragmentResults"]


class ADFFragmentResults(Results):

    def check(self):
        """Check if the calculation finished successfully.

        Overwriting the method of |MultiJob| because it only checks if the job
        finished, not how it finished.
        """
        for job in self.job.children:
            if not job.check():
                return False
        return True

    def get_properties(self):
        """Redirect to |ADFResults| of the full calculation."""
        return self.job.full.results.get_properties()

    def get_main_molecule(self):
        """Redirect to |ADFResults| of the full calculation."""
        return self.job.full.results.get_main_molecule()

    def get_input_molecule(self):
        """Redirect to |ADFResults| of the full calculation."""
        return self.job.full.results.get_input_molecule()

    def get_energy(self, unit="au"):
        """Redirect to |ADFResults| of the full calculation."""
        return self.job.full.results.get_energy(unit)

    def get_dipole_vector(self, unit="au"):
        """Redirect to |ADFResults| of the full calculation."""
        return self.job.full.results.get_dipole_vector(unit)

    def get_energy_decomposition(self, unit="kJ/mol") -> Dict[str, float]:
        """Get the energy decomposition of the fragment calculation.

        Args:
            unit (str, optional): The unit of the energy. Defaults to 'kJ/mol'.

        Returns:
            Dict[str, float]: The energy decomposition.
        """
        energy_section = self.job.full.results.read_rkf_section("Energy", file="adf")
        ret = {}
        for k in [
            "Electrostatic Energy",
            "Kinetic Energy",
            "Elstat Interaction",
            "XC Energy",
        ]:
            ret[k] = Units.convert(energy_section[k], "au", unit)

        # most information not available from the KF file
        res = self.job.full.results
        res1 = self.job.f1.results
        res2 = self.job.f2.results
        pos = -4  # position of the energy in au the output
        # E_int appears in a comment below the PEDA Table of the output
        ret["E_int"] = Units.convert(float(res.grep_output("Total Bonding Energy:")[-2].split()[pos]), "au", unit)
        ret["E_int_disp"] = Units.convert(float(res.grep_output("Dispersion Energy:")[-1].split()[pos]), "au", unit)
        ret["E_Pauli"] = Units.convert(
            float(res.grep_output("Pauli Repulsion (Delta")[-1].split()[pos]),
            "au",
            unit,
        )
        ret["E_elstat"] = Units.convert(
            float(res.grep_output("Electrostatic Interaction:")[-1].split()[pos]),
            "au",
            unit,
        )
        ret["E_orb"] = Units.convert(
            float(res.grep_output("Total Orbital Interactions:")[-1].split()[pos]),
            "au",
            unit,
        )

        ret["E_1"] = res1.get_energy(unit=unit)
        ret["E_2"] = res2.get_energy(unit=unit)

        if hasattr(self.job, "f1_opt"):
            ret["E_1_opt"] = self.job.f1_opt.results.get_energy(unit=unit)
            ret["E_prep_1"] = ret["E_1"] - ret["E_1_opt"]
        if hasattr(self.job, "f2_opt"):
            ret["E_2_opt"] = self.job.f2_opt.results.get_energy(unit=unit)
            ret["E_prep_2"] = ret["E_2"] - ret["E_2_opt"]

        if ("E_1_opt" in ret) and ("E_2_opt" in ret):
            ret["E_prep"] = ret["E_prep_1"] + ret["E_prep_2"]
            ret["E_bond"] = ret["E_int"] + ret["E_prep"]

        return ret


class ADFFragmentJob(MultiJob):
    """Subclass of |MultiJob| for ADF fragment calculations."""

    _result_type = ADFFragmentResults

    def __init__(
        self,
        fragment1: Molecule = None,
        fragment2: Molecule = None,
        fragment1_opt: Molecule = None,
        fragment2_opt: Molecule = None,
        full_settings: Settings = None,
        frag1_settings: Settings = None,
        frag2_settings: Settings = None,
        **kwargs,
    ):
        """Create an ADFFragmentJob with the given fragments.

        The optimized fragment structures can be given as arguments.
        Single point calculations on these structures will then be performed
        to obtain the preparation energy.

        Subclass of |MultiJob|.

        Args:
            fragment1 (Molecule): The first fragment.
            fragment2 (Molecule): The second fragment.
            fragment1_opt (Molecule, optional): The optimized first fragment.
            fragment2_opt (Molecule, optional): The optimized second fragment.
            full_settings (Settings): The settings for the full calculation.
            frag1_settings (Settings): Specific settings for the first fragment calculation.
            frag2_settings (Settings): Specific settings for the second fragment calculation.
            **kwargs: Further keyword arguments for |MultiJob|.
        """
        MultiJob.__init__(self, **kwargs)
        self.fragment1 = fragment1.copy() if isinstance(fragment1, Molecule) else fragment1
        self.frag1_settings = frag1_settings or Settings()
        self.fragment2 = fragment2.copy() if isinstance(fragment2, Molecule) else fragment2
        self.frag2_settings = frag2_settings or Settings()
        self.full_settings = full_settings or Settings()

        if isinstance(fragment1_opt, Molecule):
            self.fragment1_opt = fragment1_opt.copy()
        elif fragment1_opt is not None:
            raise TypeError("fragment1_opt should be a Molecule object or None")
        if isinstance(fragment2_opt, Molecule):
            self.fragment2_opt = fragment2_opt.copy()
        elif fragment2_opt is not None:
            raise TypeError("fragment2_opt should be a Molecule object or None")

    def prerun(self):  # noqa F811
        """Prepare the fragments and the full calculation."""
        self.f1 = AMSJob(
            name="frag1",
            molecule=self.fragment1,
            settings=self.settings + self.frag1_settings,
        )
        self.f2 = AMSJob(
            name="frag2",
            molecule=self.fragment2,
            settings=self.settings + self.frag2_settings,
        )

        for at in self.fragment1:
            at.properties.suffix = "adf.f=subsystem1"
        for at in self.fragment2:
            at.properties.suffix = "adf.f=subsystem2"

        self.full = AMSJob(
            name="full",
            molecule=self.fragment1 + self.fragment2,
            settings=self.settings + self.full_settings,
        )

        self.full.settings.input.adf.fragments.subsystem1 = (self.f1, "adf")
        self.full.settings.input.adf.fragments.subsystem2 = (self.f2, "adf")
        self.children = [self.f1, self.f2, self.full]

        if hasattr(self, "fragment1_opt"):
            self.f1_opt = AMSJob(
                name="frag1_opt",
                molecule=self.fragment1_opt,
                settings=self.settings + self.frag1_settings,
            )
            self.children.append(self.f1_opt)
        if hasattr(self, "fragment2_opt"):
            self.f2_opt = AMSJob(
                name="frag2_opt",
                molecule=self.fragment2_opt,
                settings=self.settings + self.frag2_settings,
            )
            self.children.append(self.f2_opt)

    @classmethod
    def load_external(cls, path: str, jobname: str = None) -> "ADFFragmentJob":
        """Load the results of the ADFFragmentJob job from an external path.

        Args:
            path (str): The path to the job. It should at least have the
                        subfolders 'frag1', 'frag2' and 'full'.
            jobname (str, optional): The name of the job. Defaults to None.

        Returns:
            ADFFragmentJob: The job with the loaded results.
        """
        if not isdir(path):
            raise FileError("Path {} does not exist, cannot load from it.".format(path))
        path = abspath(path)
        jobname = basename(path) if jobname is None else str(jobname)

        job = cls(name=jobname)
        job.path = path
        job.status = "copied"

        job.f1 = AMSJob.load_external(opj(path, "frag1"))
        job.f2 = AMSJob.load_external(opj(path, "frag2"))
        job.full = AMSJob.load_external(opj(path, "full"))
        job.children = [job.f1, job.f2, job.full]

        if isdir(opj(path, "frag1_opt")):
            job.f1_opt = AMSJob.load_external(opj(path, "frag1_opt"))
            job.children.append(job.f1_opt)
        if isdir(opj(path, "frag2_opt")):
            job.f2_opt = AMSJob.load_external(opj(path, "frag2_opt"))
            job.children.append(job.f2_opt)

        return job

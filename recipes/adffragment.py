from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.mol.molecule import Molecule

__all__ = ["ADFFragmentJob", "ADFFragmentResults"]


class ADFFragmentResults(Results):

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

    def get_energy_decomposition(self):
        """Get the energy decomposition of the fragment calculation."""
        energy_section = self.job.full.results.read_rkf_section("Energy", file="adf")
        ret = {}
        for k in ["Electrostatic Energy", "Kinetic Energy", "Elstat Interaction", "XC Energy"]:
            ret[k] = energy_section[k]
        return ret


class ADFFragmentJob(MultiJob):
    """Subclass of |MultiJob| for ADF fragment calculations."""
    _result_type = ADFFragmentResults

    def __init__(self, fragment1=None, fragment2=None, full_settings=None, **kwargs):
        """
        Args:
            fragment1 (Molecule): The first fragment.
            fragment2 (Molecule): The second fragment.
            full_settings (Settings): The settings for the full calculation.
            **kwargs: Further keyword arguments for |MultiJob|.
        """
        MultiJob.__init__(self, **kwargs)
        self.fragment1 = fragment1.copy() if isinstance(fragment1, Molecule) else fragment1
        self.fragment2 = fragment2.copy() if isinstance(fragment2, Molecule) else fragment2
        self.full_settings = full_settings or Settings()

    def prerun(self):  # noqa F811
        """Prepare the fragments and the full calculation."""
        self.f1 = AMSJob(name="frag1", molecule=self.fragment1, settings=self.settings)
        self.f2 = AMSJob(name="frag2", molecule=self.fragment2, settings=self.settings)

        for at in self.fragment1:
            at.properties.suffix = "adf.f=subsystem1"
        for at in self.fragment2:
            at.properties.suffix = "adf.f=subsystem2"

        self.full = AMSJob(
            name="full", molecule=self.fragment1 + self.fragment2, settings=self.settings + self.full_settings
        )

        self.full.settings.input.adf.fragments.subsystem1 = (self.f1, "adf")
        self.full.settings.input.adf.fragments.subsystem2 = (self.f2, "adf")

        self.children = [self.f1, self.f2, self.full]

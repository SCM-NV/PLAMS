from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.mol.molecule import Molecule

__all__ = ["ADFFragmentJob", "ADFFragmentResults"]


class ADFFragmentResults(Results):

    def get_properties(self):
        return self.job.full.results.get_properties()

    def get_main_molecule(self):
        return self.job.full.results.get_main_molecule()

    def get_input_molecule(self):
        return self.job.full.results.get_input_molecule()

    def get_energy(self, unit="au"):
        return self.job.full.results.get_energy(unit)

    def get_dipole_vector(self, unit="au"):
        return self.job.full.results.get_dipole_vector(unit)

    def get_energy_decomposition(self):
        energy_section = self.job.full.results.read_rkf_section("Energy", file="adf")
        ret = {}
        for k in ["Electrostatic Energy", "Kinetic Energy", "Elstat Interaction", "XC Energy"]:
            ret[k] = energy_section[k]
        return ret


class ADFFragmentJob(MultiJob):
    _result_type = ADFFragmentResults

    def __init__(self, fragment1=None, fragment2=None, full_settings=None, **kwargs):
        MultiJob.__init__(self, **kwargs)
        self.fragment1 = fragment1.copy() if isinstance(fragment1, Molecule) else fragment1
        self.fragment2 = fragment2.copy() if isinstance(fragment2, Molecule) else fragment2
        self.full_settings = full_settings or Settings()

    def prerun(self):  # noqa F811
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

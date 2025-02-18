from ase.calculators.calculator import Calculator, all_changes
from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob

from scm.plams import finish, Units, config
import os


class SerenityCalculator(Calculator):
    implemented_properties = ["energy"]

    def __init__(self, settings: SerenitySettings):
        super().__init__()
        self.settings = settings
        self.results = dict()

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        config.default_jobmanager.hashing = None
        config.jobmanager.hashing = None
        atoms.write("my_system.xyz")
        for sys in self.settings.input.system:
            self.settings.input.system[sys].geometry = os.path.abspath("my_system.xyz")
        job = SerenityJob(settings=self.settings, name="serenity_job")
        job.run()
        self.results = dict()
        self.results["energy"] = job.results.get_energy() * Units.convert(1.0, "hartree", "eV")

    def clean_exit(self):
        finish()

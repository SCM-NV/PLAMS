from collections import OrderedDict

from scm.plams.core.basejob import MultiJob
from scm.plams.core.functions import add_to_instance
from scm.plams.core.results import Results
from scm.plams.interfaces.adfsuite.ams import AMSJob

__all__ = ["ReorganizationEnergyJob", "ReorganizationEnergyResults"]


# using this function to pass a molecule inside a MultiJob ensures proper parallel execution
def pass_molecule(source, target):
    @add_to_instance(target)
    def prerun(self):  # noqa F811
        self.molecule = source.results.get_main_molecule()


class ReorganizationEnergyResults(Results):
    """Results class for reorganization energy."""

    def reorganization_energy(self, unit="au"):
        energies = self.get_all_energies(unit)

        reorganization_energy = (
            energies["state B geo A"]
            - energies["state A geo A"]
            + energies["state A geo B"]
            - energies["state B geo B"]
        )
        return reorganization_energy

    def get_all_energies(self, unit="au"):
        energies = {}
        energies["state A geo A"] = self.job.children["go_A"].results.get_energy(unit=unit)
        energies["state B geo B"] = self.job.children["go_B"].results.get_energy(unit=unit)
        energies["state A geo B"] = self.job.children["sp_A_for_geo_B"].results.get_energy(unit=unit)
        energies["state B geo A"] = self.job.children["sp_B_for_geo_A"].results.get_energy(unit=unit)

        return energies


class ReorganizationEnergyJob(MultiJob):
    """A class for calculating the reorganization energy using AMS.
    Given two states, A and B, the reorganization energy is defined as follows:

    reorganization energy =
    E(state B at optimal geometry for state A) -
    E(state A at optimal geometry for state A) +
    E(state A at optimal geometry for state B) -
    E(state B at optimal geometry for state B)

    This job will run two geometry optimizations and two single point calculations.
    """

    _result_type = ReorganizationEnergyResults

    def __init__(self, molecule, common_settings, settings_state_A, settings_state_B, **kwargs):
        """
        molecule: the molecule
        common_settings: a setting object for an AMSJob containing the shared settings for all the calculations
        settings_state_A: Setting object for an AMSJob containing exclusivelt the options defining the state A (e.g. charge and spin)
        settings_state_B: Setting object for an AMSJob containing exclusivelt the options defining the state B (e.g. charge and spin)
        kwargs: other options to be passed to the MultiJob constructor
        """
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)

        go_settings = common_settings.copy()
        go_settings.input.ams.task = "GeometryOptimization"

        sp_settings = common_settings.copy()
        sp_settings.input.ams.task = "SinglePoint"

        # copy the settings so that we wont modify the ones provided as input by the user
        settings_state_A = settings_state_A.copy()
        settings_state_B = settings_state_B.copy()
        # In case the charge key is not specified, excplicitely set the value to 0.
        # This is to prevent the charge in molecule.properties.charge (set by get_main_molecule())
        # to be used in case of neutral systems
        for s in [settings_state_A, settings_state_B]:
            if not "charge" in s.input.ams.system:
                s.input.ams.system.charge = 0

        self.children["go_A"] = AMSJob(molecule=molecule, settings=go_settings + settings_state_A, name="go_A")
        self.children["go_B"] = AMSJob(molecule=molecule, settings=go_settings + settings_state_B, name="go_B")
        self.children["sp_A_for_geo_B"] = AMSJob(settings=sp_settings + settings_state_A, name="sp_A_geo_B")
        self.children["sp_B_for_geo_A"] = AMSJob(settings=sp_settings + settings_state_B, name="sp_B_geo_A")
        pass_molecule(self.children["go_A"], self.children["sp_B_for_geo_A"])
        pass_molecule(self.children["go_B"], self.children["sp_A_for_geo_B"])

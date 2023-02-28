from collections import OrderedDict
from ...core.functions import add_to_instance
from ...core.basejob import MultiJob
from ...core.results import Results
from ...core.settings import Settings
from ...mol.molecule import Molecule
from ...mol.atom import Atom
from ...interfaces.adfsuite.ams import AMSJob, AMSResults
from ...tools.units import Units
from .amsmdjob import AMSNVTJob
import numpy as np
from scipy.optimize import curve_fit

__all__ = ['AMSMDScanDensityJob', 'AMSMDScanDensityResults']

class AMSMDScanDensityResults(AMSResults):
    """Results class for AMSMDScanDensityJob
    """
    def get_lowest_energy_index(self, variable='Energy', history_section='History'):
        """
        Returns the 1-based index of the lowest energy molecule
        """
        energies = self.get_history_property(variable, history_section)
        minindex = np.argmin(energies) + 1
        return minindex

    def get_lowest_energy_molecule(self, variable='TotalEnergy'):
        return self.get_history_molecule(self.get_lowest_energy_index(variable, 'MDHistory'))
        

class AMSMDScanDensityJob(AMSNVTJob):
    """A class for scanning the density using MD Deformations
    """

    _result_type = AMSMDScanDensityResults

    def __init__(self, molecule, scan_density_upper=1.4, startstep=None, **kwargs):
        AMSNVTJob.__init__(self, molecule=molecule, **kwargs)

        from_density = self.molecule.get_density() * 1e-3
        orig_length = self.molecule.cell_lengths()
        density_ratio = from_density / scan_density_upper
        new_length = [x *  density_ratio**0.333333 for x in orig_length]

        nsteps = self.settings.input.ams.MolecularDynamics.NSteps

        self.scan_density_upper = scan_density_upper
        self.startstep = startstep or 1

        s = Settings()
        s.input.ams.MolecularDynamics.Deformation.TargetLength = ' '.join([str(x) for x in new_length]) 
        s.input.ams.MolecularDynamics.Deformation.StartStep = self.startstep

        self.settings += s




from collections import OrderedDict
from ...core.functions import add_to_instance
from ...core.basejob import MultiJob
from ...core.results import Results
from ...core.settings import Settings
from ...mol.molecule import Molecule
from ...mol.atom import Atom
from ...interfaces.adfsuite.ams import AMSJob
from ...tools.units import Units
from ...interfaces.molecule.packmol import packmol_mixture
from .scandensity import AMSMDScanDensityJob
from .amsmdjob import AMSNVTJob, AMSNPTJob, AMSNVEJob
import numpy as np
from scipy.optimize import curve_fit

__all__ = ['EquilibrateDensityJob', 'EquilibrateDensityResults']

class EquilibrateDensityResults(Results):
    """Results class for EquilibrateDensityJob
    """
        
    def get_equilibrated_molecule(self, equilibration_fraction=0.667, return_index=False):
        return self.job.children['npt'].results.get_equilibrated_molecule(equilibration_fraction=equilibration_fraction, return_index=return_index)

    def rkfpath(self):
        """ Returns the path to ams.rkf from the npt job """
        return self.job.children['npt'].rkfpath()

class EquilibrateDensityJob(MultiJob):
    """A class for equilibrating the density at a certain temperature and pressure
    """

    _result_type = EquilibrateDensityResults

    def _default_settings(self):
        s = Settings()
        s.input.ForceField.Type = 'GAFF'
        s.input.ForceField.AnteChamberIntegration = 'Yes'
        return s

    def _create_scan_density_job(self, initial_molecule):
        name = 'scan_density'
        self.children[name] = AMSMDScanDensityJob(
            name=name, 
            nsteps=self.nsteps[name],
            settings=self.settings,
            scan_density_upper=self.scan_density_upper, 
            molecule = initial_molecule,
            temperature = self.temperature,
        )

        return self.children[name]

    def _create_nvt_pre_eq_job(self, scan_density_job):
        name = 'nvt_pre_eq'
        job = AMSNVTJob(
            name=name,
            settings=self.settings, 
            nsteps=self.nsteps[name],
            timestep=self.timestep,
            tau = 10,
            thermostat = 'Berendsen',
            temperature = self.temperature,
            velocities = self.temperature,
            writevelocities = False,
            writebonds = True,
            writemolecules = False,
        )

        if scan_density_job is not None:
            @add_to_instance(job)
            def prerun(self):
                self.molecule = scan_density_job.results.get_lowest_energy_molecule()
        else:
            job.molecule = self.initial_molecule

        self.children[name] = job
        return self.children[name]

    def _create_npt_job(self, nvt_pre_eq_job):
        name = 'npt'
        job = AMSNPTJob(
            name=name,
            settings=self.settings,
            nsteps=self.nsteps[name],
            timestep=self.timestep,
            tau=100,
            thermostat='NHC',
            temperature=self.temperature,
            barostat='MTK',
            pressure=self.pressure,
            barostat_tau=1000,
            equal='XYZ',
            writevelocities=False,
            writebonds=True,
            writecharges=False
        )

        @add_to_instance(job)
        def prerun(self):
            self.get_velocities_from(nvt_pre_eq_job, update_molecule=True)

        self.children[name] = job
        return self.children[name]
        
    def __init__(self, 
                 molecule, 
                 settings=None, 
                 name='equilibrate_density',
                 nsteps=None, 
                 temperature=300, 
                 pressure=1.0, 
                 scan_density=True, 
                 scan_density_upper=1.5, 
                 **kwargs):
        """
        molecule: Molecule
            3D molecule (liquid/gas with multiple molecules). 

        settings: Settings
            All non-AMS-Driver settings, for example (``s.input.forcefield.type = 'GAFF'``, ``s.runscript.nproc = 1``)

        nsteps: dict
            Dictionary where the default key-values pairs are. Any keys present in the dictionary will override the default values.

            .. code-block:: python

                nsteps = {
                    'scan_density': 5000,
                    'nvt_pre_eq': 1000,
                    'npt': 100000,
                }

        temperature: float
            Temperature in K.

        pressure: float
            The pressure in bar.

        kwargs: other options to be passed to the MultiJob constructor (for example the name)
        """
        MultiJob.__init__(self, children=OrderedDict(), name=name, **kwargs)

        self.scan_density_upper = scan_density_upper
        self.timestep = 1.0
        self.temperature = temperature
        self.pressure = pressure
        self.nsteps = {
            'scan_density': 5000,
            'nvt_pre_eq': 1000,
            'npt': 100000,
        }
        if nsteps:
            self.nsteps.update(nsteps)

        self.settings = settings.copy() if settings is not None else self._default_settings()

        if scan_density:
            scan_density_job = self._create_scan_density_job(molecule)
        else:
            scan_density_job = None

        nvt_pre_eq_job = self._create_nvt_pre_eq_job(scan_density_job)

        npt_job = self._create_npt_job(nvt_pre_eq_job)


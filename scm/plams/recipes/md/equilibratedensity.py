from collections import OrderedDict

from scm.plams.core.basejob import MultiJob
from scm.plams.core.functions import add_to_instance
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.recipes.md.amsmdjob import AMSNPTJob, AMSNVTJob
from scm.plams.recipes.md.scandensity import AMSMDScanDensityJob

__all__ = ["EquilibrateDensityJob", "EquilibrateDensityResults"]


class EquilibrateDensityResults(Results):
    """Results class for EquilibrateDensityJob"""

    def get_equilibrated_molecule(self, equilibration_fraction=0.667, return_index=False):
        return self.job.children["npt"].results.get_equilibrated_molecule(
            equilibration_fraction=equilibration_fraction, return_index=return_index
        )

    def rkfpath(self):
        """Returns the path to ams.rkf from the npt job"""
        return self.job.children["npt"].rkfpath()


class EquilibrateDensityJob(MultiJob):
    """A class for equilibrating the density at a certain temperature and pressure"""

    _result_type = EquilibrateDensityResults

    def _default_settings(self):
        s = Settings()
        s.input.ForceField.Type = "GAFF"
        s.input.ForceField.AnteChamberIntegration = "Yes"
        return s

    def _create_scan_density_job(self, initial_molecule):
        name = "scan_density"
        self.children[name] = AMSMDScanDensityJob(
            name=name,
            nsteps=self.nsteps[name],
            settings=self.settings,
            scan_density_upper=self.scan_density_upper,
            molecule=initial_molecule,
            **self.kwargs,
        )

        return self.children[name]

    def _create_nvt_pre_eq_job(self, scan_density_job):
        name = "nvt_pre_eq"
        job = AMSNVTJob(name=name, settings=self.settings, nsteps=self.nsteps[name], **self.kwargs)

        if scan_density_job is not None:

            @add_to_instance(job)
            def prerun(self):  # noqa F811
                self.molecule = scan_density_job.results.get_lowest_energy_molecule()

        else:
            job.molecule = self.initial_molecule

        self.children[name] = job
        return self.children[name]

    def _create_npt_job(self, nvt_pre_eq_job):
        name = "npt"
        job = AMSNPTJob.restart_from(
            nvt_pre_eq_job,
            name=name,
            use_prerun=True,
            settings=self.settings,
            nsteps=self.nsteps[name],
            **self.kwargs,
        )

        # @add_to_instance(job)
        # def prerun(self):  # noqa F811
        # self.get_velocities_from(nvt_pre_eq_job, update_molecule=True)

        self.children[name] = job
        return self.children[name]

    def __init__(
        self,
        molecule,
        settings=None,
        name="equilibrate_density",
        nsteps=None,
        scan_density=True,
        scan_density_upper=1.5,
        **kwargs,
    ):
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

        kwargs: various options
            Other options for AMSMDJob (e.g. temperature, pressure, timestep)

        """
        MultiJob.__init__(self, children=OrderedDict(), name=name)

        self.scan_density_upper = scan_density_upper
        self.timestep = 1.0
        self.nsteps = {
            "scan_density": 5000,
            "nvt_pre_eq": 1000,
            "npt": 100000,
        }
        if nsteps:
            self.nsteps.update(nsteps)

        self.settings = settings.copy() if settings is not None else self._default_settings()

        self.kwargs = kwargs

        if scan_density:
            scan_density_job = self._create_scan_density_job(molecule)
        else:
            scan_density_job = None

        nvt_pre_eq_job = self._create_nvt_pre_eq_job(scan_density_job)

        self._create_npt_job(nvt_pre_eq_job)

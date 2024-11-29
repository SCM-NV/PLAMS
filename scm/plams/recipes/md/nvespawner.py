from collections import OrderedDict

from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.core.enums import JobStatus

import numpy as np
from scm.plams.recipes.md.amsmdjob import AMSNVEJob

__all__ = ["AMSNVESpawnerJob", "AMSNVESpawnerResults"]


class AMSNVESpawnerResults(Results):
    """Results class for AMSNVESpawnerJob"""

    def get_all_dipole_derivatives_acf(
        self, start_fs=0, end_fs=None, every_fs=None, max_dt_fs=None, x=True, y=True, z=True, normalize=False
    ):
        ret_x, ret_y = [], []
        for job in self.job.children.values():
            xv, yv = job.results.get_dipole_derivatives_acf(
                start_fs=start_fs,
                end_fs=end_fs,
                every_fs=every_fs,
                max_dt_fs=max_dt_fs,
                x=x,
                y=y,
                z=z,
                normalize=normalize,
            )
            ret_x.append(xv)
            ret_y.append(yv)

        return ret_x, ret_y

    def get_dipole_derivatives_acf(
        self, start_fs=0, end_fs=None, every_fs=None, max_dt_fs=None, x=True, y=True, z=True, normalize=False
    ):
        all_x, all_y = self.get_all_dipole_derivatives_acf(
            start_fs=start_fs, end_fs=end_fs, every_fs=every_fs, max_dt_fs=max_dt_fs, x=x, y=y, z=z, normalize=normalize
        )
        return np.mean(all_x, axis=0), np.mean(all_y, axis=0)

    def get_mean_temperature(self) -> float:
        avg_T = np.mean([job.results.readrkf("MDResults", "MeanTemperature") for job in self.job.children.values()])

        return avg_T


class AMSNVESpawnerJob(MultiJob):
    """A class for running multiple NVE simulations with initial structures/velocities taken from an NVT trajectory. The NVT trajectory must contain the velocities!"""

    _result_type = AMSNVESpawnerResults

    def _default_settings(self):
        s = Settings()
        s.input.ForceField.Type = "GAFF"
        s.input.ForceField.AnteChamberIntegration = "Yes"
        return s

    def __init__(self, previous_job, n_nve=1, name="nvespawnerjob", **kwargs):  # needs to be finished
        """
        previous_job: AMSJob
            An AMSJob with an MD trajectory. Must contain velocities (WriteVelocities Yes). Note that the trajectory should have been equilibrated before it starts.

        n_nve : int
            The number of NVE simulations to spawn

        All other settings can be set as for an AMSNVEJob (e.g. ``nsteps``).

        """
        MultiJob.__init__(self, children=OrderedDict(), name=name)

        self.previous_job = previous_job
        self.n_nve = n_nve
        self.nve_constructor_settings = kwargs
        self.nve_jobs = []

    def prerun(self):  # noqa F811
        """
        Constructs the children jobs
        """

        # use previously run previous_job
        assert self.previous_job.status != JobStatus.CREATED, "You can only pass in a finished AMSJob"
        try:
            self.previous_job.results.readrkf("MDHistory", "Velocities(1)")
        except KeyError:
            raise KeyError("Couldn't read velocities from {}".format(self.previous_job.results.rkfpath()))

        nframes_in_history = self.previous_job.results.readrkf("History", "nEntries")

        if self.n_nve > 0:
            interval = nframes_in_history // self.n_nve
            frames = np.linspace(interval, nframes_in_history, self.n_nve, dtype=int)
            for i, frame in enumerate(frames):
                # the nve jobs only get different inside prerun, so need to do something to make them different to prevent the PLAMS rerun prevention
                name = f"nve{i+1}"
                self.settings.input.ams["#"] = name  # add a comment line

                self.children[name] = AMSNVEJob.restart_from(
                    self.previous_job, frame=frame, name=name, use_prerun=True, **self.nve_constructor_settings
                )
                self.children[name].from_frame = frame

                self.nve_jobs.append(self.children[name])

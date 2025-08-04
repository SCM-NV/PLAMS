import os
from collections import OrderedDict
import numpy as np

from scm.plams.core.basejob import MultiJob
from scm.plams.core.functions import add_to_instance, requires_optional_package
from scm.plams.core.results import Results
from scm.plams.interfaces.adfsuite.ams import AMSJob, AMSResults
from scm.plams.mol.molecule import Molecule
from scm.plams.recipes.md.amsmdjob import AMSNVTJob
from scm.plams.tools.kftools import KFFile

__all__ = ["AMSGreenKuboViscosityJob", "AMSGreenKuboViscosityResults", "AMSPartialGreenKuboViscosityJob"]


@requires_optional_package("scipy")
def get_viscosity(job, max_dt_fs=50000, reuse=False):
    """
    job: AMSGreenKuboViscosityJob
        The job

    Returns the viscosity in mPa*s, and writes

    * viscosity.txt: viscosity in mPa*s
    * green_kubo_viscosity.txt: viscosity vs. time (curve that should converge to the viscosity as t goes to inifinity)
    * fit.txt : A fit of the form A*(1-exp(-x/tau)) to green_kubo_viscosity.txt; A is the final viscosity.
    * moving_avg_visc.txt : viscosity vs. time (moving average of green_kubo_viscosity.txt)
    * info.txt : Contains the upper time-limit of the fit and the viscosity

    The upper limit for the fit is determined by when the moving average starts to decrease.

    max_dt_fs : float
        Maximum correlation time in femtoseconds

    reuse : bool
        Will return the value of viscosity.txt if it exists (and not recalculate the autocorrelation function/integral)
    """
    from scipy.optimize import curve_fit

    def f(x, A, tau):
        return A * (1 - np.exp(-x / tau))  # A is the final converged viscosity

    if reuse and os.path.exists(job.path + "/viscosity.txt"):
        return float(np.loadtxt(job.path + "/viscosity.txt"))

    AMSGreenKuboViscosityResults._accumulate_results(job.path, previous_pressuretensors=job.previous_pressuretensors)
    pressuretensor_npy = job.path + "/concatenated_pressuretensor.npy"

    pressuretensor = np.load(pressuretensor_npy)

    for child in job.children.values():
        volume = child.results.get_main_molecule().unit_cell_volume()
        break

    temperature = job.temperature
    max_dt = int(max_dt_fs / (job.samplingfreq * job.timestep))

    t, visc = AMSResults._get_green_kubo_viscosity(
        pressuretensor=pressuretensor,
        max_dt=max_dt,
        time_step=job.samplingfreq * job.timestep,
        volume=volume,
        temperature=temperature,
    )

    A = np.stack((t, visc), axis=1)
    np.savetxt(job.path + "/green_kubo_viscosity.txt", A, header="Time(fs) Viscosity(mPa*s)")

    window = min(5000, int(t[-1] / 3))
    window = max(window, 1)  # window must be greater than 0
    moving_avg_t = np.convolve(t, np.ones(window) / window, mode="valid")
    moving_avg_visc = np.convolve(visc, np.ones(window) / window, mode="valid")

    flat_region_from = min(1000, int(t[-1] / 10))  # femtoseconds
    flat_region_from_indices = moving_avg_t > flat_region_from
    allowed_t = moving_avg_t[flat_region_from_indices]
    allowed_visc = moving_avg_visc[flat_region_from_indices]
    d_moving_avg_visc = np.diff(allowed_visc)
    zero_indices = np.argwhere(np.diff(np.sign(d_moving_avg_visc))).flatten()
    if len(zero_indices) == 0:
        max_index = len(allowed_visc) - 1
    else:
        max_index = zero_indices[0]

    fit_until_t = allowed_t[max_index]

    fit_x_indices = t < fit_until_t
    fit_x = t[fit_x_indices]
    fit_y = visc[fit_x_indices]

    popt, _ = curve_fit(f, fit_x, fit_y, p0=(1.0, 1000))
    prediction = f(t, popt[0], popt[1])
    with open(job.path + "/viscosity.txt", "w") as f:
        f.write("{}\n".format(popt[0]))

    with open(job.path + "/info.txt", "w") as f:
        f.write("FITTED UNTIL t = {} fs \n".format(fit_until_t))
        f.write("VISCOSITY = {} mPa*s\n".format(popt[0]))

    A = np.stack((t, prediction), axis=1)
    np.savetxt(job.path + "/fit.txt", A, header="Time(fs) Viscosity(mPa*s)")

    A = np.stack((moving_avg_t, moving_avg_visc), axis=1)
    np.savetxt(job.path + "/moving_avg_visc.txt", A, header="Time(fs) Viscosity(mPa*s)")

    return popt[0]


class AMSPartialGreenKuboViscosityResults(AMSResults):
    def get_pressuretensor_history(self):
        pt = self.get_history_property("PressureTensor", "MDHistory")
        return np.array(pt)


class AMSPartialGreenKuboViscosityJob(AMSNVTJob):
    _result_type = AMSPartialGreenKuboViscosityResults

    def __init__(self, keep_trajectory: bool = False, max_dt_fs: float = 50000, **kwargs):
        self.keep_trajectory = keep_trajectory
        self.max_dt_fs = max_dt_fs

        kwargs.pop("calcpressure", True)
        AMSNVTJob.__init__(
            self,
            writevelocities=kwargs.get("writevelocities", False),
            writemolecules=kwargs.get("writemolecules", False),
            writebonds=kwargs.get("writebonds", False),
            writecharges=kwargs.get("writecharges", False),
            calcpressure=True,
            **kwargs,
        )

    def postrun(self):
        pt = self.results.get_pressuretensor_history()
        np.save(self.path + "/pressuretensor.npy", pt)

        if not self.keep_trajectory:
            kf = KFFile(os.path.join(self.path, "ams.rkf"))
            kf.delete_section("History")

        self.results.collect()

        get_viscosity(self.parent, max_dt_fs=self.max_dt_fs, reuse=False)


class AMSGreenKuboViscosityResults(Results):
    """Results class for AMSGreenKuboViscosityJob"""

    def get_viscosity(self, max_dt_fs=50000, reuse=False):
        """
        Returns the viscosity in mPa*s, and writes

        * viscosity.txt: viscosity in mPa*s
        * green_kubo_viscosity.txt: viscosity vs. time (curve that should converge to the viscosity as t goes to inifinity)
        * fit.txt : A fit of the form A*(1-exp(-x/tau)) to green_kubo_viscosity.txt; A is the final viscosity.
        * moving_avg_visc.txt : viscosity vs. time (moving average of green_kubo_viscosity.txt)
        * info.txt : Contains the upper time-limit of the fit and the viscosity

        The upper limit for the fit is determined by when the moving average starts to decrease.

        max_dt_fs : float
            Maximum correlation time in femtoseconds

        reuse : bool
            Will return the value of viscosity.txt if it exists (and not recalculate the autocorrelation function/integral)
        """
        return get_viscosity(self.job, max_dt_fs=max_dt_fs, reuse=reuse)

    @staticmethod
    @requires_optional_package("natsort")
    def _accumulate_results(path=None, previous_pressuretensors=None):
        import glob
        from natsort import natsorted

        path = path or os.getcwd()
        complete_list = []

        if previous_pressuretensors:
            complete_list.extend(previous_pressuretensors)

        for npy in natsorted(glob.glob(f"{path}/*/pressuretensor.npy")):
            complete_list.append(np.load(npy)[:-1])

        if len(complete_list) == 0:
            raise ValueError("Couldn't find any pressuretensor.npy files!")

        complete_list = np.concatenate(complete_list).reshape(-1, 6)
        np.save(path + "/concatenated_pressuretensor.npy", complete_list)
        return complete_list

    def accumulate_results(self):
        return self._accumulate_results(path=self.job.path, previous_pressuretensors=self.job.previous_pressuretensors)


class AMSGreenKuboViscosityJob(MultiJob):
    """A class for calculating the Green-Kubo viscosity of a liquid"""

    _result_type = AMSGreenKuboViscosityResults

    def __init__(
        self,
        molecule: Molecule = None,
        name: str = "greenkuboviscosity",
        nsteps: int = 100000,
        fragment_length: int = 10000,
        samplingfreq: int = 5,
        timestep: float = 1.0,
        temperature: float = 300,
        keep_trajectory: bool = False,
        restart_from=None,
        previous_pressuretensors=None,
        **kwargs,
    ):
        """
        molecule: Molecule
            3D molecule (liquid) with a low density.

        nsteps: int
            Total number of MD steps in the simulation

        fragment_length: int
            Length per "fragment simulation". Should normally not be greater than 10000, as that will slow down the writing of the trajectory.

        samplingfreq: int
            How often to sample

        timestep: float
            The timestep in femtoseconds.

        temperature: float
            The temperature in K

        keep_trajectory: bool
            Whether to keep the structures in the trajectories. This takes up a lot of disk space.

        restart_from: AMSJob
            A previous job to restart from (structure and velocities). If given, ``molecule`` is ignored. Example: ``restart_from=AMSJob.load_external('plams_workdir/greenkubo/step59')

        previous_pressuretensors: str or list of np.array
            Load previous pressure tensors. Should normally be the path to "concatenated_pressuretensors.npy" when you use ``restart_from``.

        kwargs: other options to be passed to the AMSNVTJob constructor (engine settings, ...)
        """
        MultiJob.__init__(self, children=OrderedDict(), name=name)

        self.timestep = timestep
        self.temperature = temperature
        self.samplingfreq = samplingfreq
        self.nsteps = nsteps
        self.fragment_length = fragment_length
        self.keep_trajectory = keep_trajectory
        self.restart_from = restart_from

        self.n_fragments = self.nsteps // self.fragment_length
        self.n_fragments = max(1, self.n_fragments)  # at least 1 fragment

        if self.n_fragments > 100:
            # this check is here because you otherwise run into some maximum recursion error around n_fragments = 150
            raise ValueError(
                f"AMSGreenKuboViscosityJob can be used with most 100 fragments/steps, current value: {self.n_fragments}. Either decrease nsteps (current: {self.nsteps}) or increase fragment_length (current: {self.fragment_length})."
            )

        self.previous_pressuretensors = previous_pressuretensors or []
        if isinstance(self.previous_pressuretensors, str):
            self.previous_pressuretensors = [np.load(self.previous_pressuretensors)]

        if restart_from:
            if not isinstance(restart_from, AMSJob):
                raise ValueError(f"restart_from must be AMSJob, got {type(restart_from)}")
            self.children["step1"] = AMSPartialGreenKuboViscosityJob(
                name="step1",
                keep_trajectory=self.keep_trajectory,
                timestep=self.timestep,
                temperature=self.temperature,
                nsteps=self.fragment_length,
                samplingfreq=self.samplingfreq,
                **kwargs,
            )

            @add_to_instance(self.children["step1"])
            def prerun(self):  # noqa F811
                self.get_velocities_from(restart_from, update_molecule=True)

        else:
            self.children["step1"] = AMSPartialGreenKuboViscosityJob(
                name="step1",
                molecule=molecule,
                keep_trajectory=self.keep_trajectory,
                timestep=self.timestep,
                temperature=self.temperature,
                nsteps=self.fragment_length,
                samplingfreq=self.samplingfreq,
                **kwargs,
            )

        previous_name = "step1"
        previous_names = dict()
        kwargs.pop("velocities", None)
        for i in range(2, self.n_fragments + 1):
            name = f"step{i}"
            previous_names[name] = previous_name

            job = AMSPartialGreenKuboViscosityJob(
                name=name,
                keep_trajectory=keep_trajectory,
                timestep=self.timestep,
                temperature=self.temperature,
                nsteps=fragment_length,
                samplingfreq=self.samplingfreq,
                **kwargs,
            )

            @add_to_instance(job)
            def prerun(self):  # noqa F811
                prev = previous_names[self.name]
                self.get_velocities_from(self.parent.children[prev], update_molecule=True)

            self.children[name] = job
            previous_name = name

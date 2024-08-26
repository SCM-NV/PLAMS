from scm.plams.core.functions import add_to_instance
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob, AMSResults
import numpy as np
from scm.plams.tools.kftools import KFFile
from scm.plams.tools.units import Units
from typing import Union
import scm.plams as plams

__all__ = ["AMSMDJob", "AMSNVEJob", "AMSNVTJob", "AMSNPTJob"]


class AMSMDJob(AMSJob):
    """

    molecule: Molecule
        The initial structure

    name: str
        The name of the job

    settings: Settings
        Settings for the job. You should normally not populate neither settings.input.ams.MolecularDynamics nor settings.input.ams.Task

    velocities: float or AMSJob or str (path/to/ams.rkf) or 2-tuple (path/to/ams.rkf, frame-number)
        If float, it is taken as the temperature. If AMSJob or str, the velocities are taken from the EndVelocities section of the corresponding ams.rkf file. If 2-tuple, the first item must be a path to an ams.rkf, and the second item an integer specifying the frame number - the velocities are then read from History%Velocities(framenumber).

    timestep: float
        Timestep

    samplingfreq: int
        Sampling frequency

    nsteps: int
        Length of simulation

    **Trajectory options**:

    checkpointfrequency: int
        How frequently to write MDStep*.rkf checkpoint files

    writevelocities : bool
        Whether to save velocities to ams.rkf (needed for example to restart from individual frames or to calculate velocity autocorrelation functions)

    writebonds: bool
        Whether to save bonds to ams.rkf

    writemolecules: bool
        Whether to write molecules to ams.rkf

    writeenginegradients: bool
        Whether to save engine gradients (negative of atomic forces) to ams.rkf

    **Thermostat (NVT, NPT) options**:

    thermostat: str
        'Berendsen' or 'NHC'

    tau: float
        Thermostat time constant (fs)

    temperature: float or tuple of floats
        Temperature (K). If a tuple/list of floats, the Thermostat.Duration option will be set to evenly divided intervals.

    thermostat_region: str
        Region for which to apply the thermostat

    **Barostat (NPT) options**:

    barostat: str
        'Berendsen' or 'MTK'

    barostat_tau: float
        Barostat time constant (fs)

    pressure: float
        Barostat pressure (pascal)

    equal: str
        'XYZ' etc.

    scale: str
        'XYZ' etc.

    constantvolume: bool
        Constant volume?

    **Other options**:

    calcpressure: bool
        Whether to calculate pressure for each frame.

    binlog_time: bool
        Whether to log the time at every timestep in the BinLog section on ams.rkf

    binlog_pressuretensor: bool
        Whether to log the pressure tensor at every timestep in the BinLog section on ams.rkf



    """

    default_nsteps = 1000
    default_timestep = 0.25
    default_samplingfreq = 100

    default_thermostat = "NHC"
    default_temperature = 300
    default_tau_multiplier = 400  # get tau by multiplying the timestep with this number
    default_thermostat_region = None

    default_barostat = "MTK"
    default_pressure = 1e5
    default_barostat_tau_multiplier = 4000  # get barostat_tau by multiplying the timestep with this number
    default_scale = "XYZ"
    default_equal = "None"
    default_constantvolume = "False"

    default_checkpointfrequency = 1000000
    default_writevelocities = "True"
    default_writebonds = "True"
    default_writecharges = "True"
    default_writemolecules = "True"
    default_writeenginegradients = "False"

    default_calcpressure = "False"
    default_binlog_time = "False"
    default_binlog_pressuretensor = "False"
    default_binlog_dipolemoment = "False"

    def __init__(
        self,
        velocities=None,
        timestep=None,
        samplingfreq: int = None,
        nsteps: int = None,
        checkpointfrequency=None,
        writevelocities=None,
        writebonds=None,
        writemolecules=None,
        writecharges=None,
        writeenginegradients=None,
        calcpressure=None,
        molecule: Union[plams.Molecule, AMSJob, AMSResults] = None,
        temperature: float = None,
        thermostat=None,
        tau=None,
        thermostat_region=None,
        pressure=None,
        barostat=None,
        barostat_tau=None,
        scale=None,
        equal=None,
        constantvolume=None,
        binlog_time=None,
        binlog_pressuretensor=None,
        binlog_dipolemoment=None,
        _enforce_thermostat=False,
        _enforce_barostat=False,
        **kwargs,
    ):
        """ """

        if isinstance(molecule, AMSJob):
            molecule = molecule.results.get_main_molecule()
        if isinstance(molecule, AMSResults):
            molecule = molecule.get_main_molecule()
        AMSJob.__init__(self, molecule=molecule, **kwargs)

        self.settings.input.ams.Task = "MolecularDynamics"
        mdsett = self.settings.input.ams.MolecularDynamics

        mdsett.TimeStep = timestep or mdsett.TimeStep or self.default_timestep
        mdsett.Trajectory.SamplingFreq = samplingfreq or mdsett.Trajectory.SamplingFreq or self.default_samplingfreq
        mdsett.NSteps = nsteps or mdsett.NSteps or self.default_nsteps

        mdsett.Trajectory.WriteVelocities = (
            str(writevelocities)
            if writevelocities is not None
            else mdsett.Trajectory.WriteVelocities or self.default_writevelocities
        )
        mdsett.Trajectory.WriteBonds = (
            str(writebonds) if writebonds is not None else mdsett.Trajectory.WriteBonds or self.default_writebonds
        )
        mdsett.Trajectory.WriteMolecules = (
            str(writemolecules)
            if writemolecules is not None
            else mdsett.Trajectory.WriteMolecules or self.default_writemolecules
        )
        mdsett.Trajectory.WriteCharges = (
            str(writecharges)
            if writecharges is not None
            else mdsett.Trajectory.WriteCharges or self.default_writecharges
        )
        mdsett.Trajectory.WriteEngineGradients = (
            str(writeenginegradients)
            if writeenginegradients is not None
            else mdsett.Trajectory.WriteEngineGradients or self.default_writeenginegradients
        )
        mdsett.CalcPressure = (
            str(calcpressure) if calcpressure is not None else mdsett.CalcPressure or self.default_calcpressure
        )
        mdsett.Checkpoint.Frequency = (
            checkpointfrequency or mdsett.Checkpoint.Frequency or self.default_checkpointfrequency
        )

        mdsett.BinLog.Time = (
            str(binlog_time) if binlog_time is not None else mdsett.BinLog.Time or self.default_binlog_time
        )
        mdsett.BinLog.PressureTensor = (
            str(binlog_pressuretensor)
            if binlog_pressuretensor is not None
            else mdsett.BinLog.PressureTensor or self.default_binlog_pressuretensor
        )
        mdsett.BinLog.DipoleMoment = (
            str(binlog_dipolemoment)
            if binlog_dipolemoment is not None
            else mdsett.BinLog.DipoleMoment or self.default_binlog_dipolemoment
        )

        if velocities is None and temperature is not None:
            velocities = self._get_initial_temperature(temperature)
        self.settings += self._velocities2settings(velocities)

        if temperature or thermostat or _enforce_thermostat:
            self.settings.update(
                self._get_thermostat_settings(
                    thermostat=thermostat,
                    temperature=temperature,
                    tau=tau,
                    thermostat_region=thermostat_region,
                    nsteps=int(mdsett.NSteps),
                )
            )

        if pressure or barostat or _enforce_barostat:
            self.settings.update(
                self._get_barostat_settings(
                    pressure=pressure,
                    barostat=barostat,
                    barostat_tau=barostat_tau,
                    scale=scale,
                    equal=equal,
                    constantvolume=constantvolume,
                )
            )

    def remove_blocks(self, blocks=None):
        if blocks:
            for block in blocks:
                if block in self.settings.input.ams.MolecularDynamics:
                    del self.settings.input.ams.MolecularDynamics[block]

    @staticmethod
    def _get_initial_temperature(temperature):
        try:
            return temperature[0]
        except TypeError:
            return temperature

    @staticmethod
    def _velocities2settings(velocities):
        s = Settings()
        if isinstance(velocities, int) or isinstance(velocities, float) or velocities is None or velocities is False:
            s.input.ams.MolecularDynamics.InitialVelocities.Type = "Random"
            s.input.ams.MolecularDynamics.InitialVelocities.Temperature = velocities or AMSMDJob.default_temperature
        elif isinstance(velocities, tuple):
            # file and frame number
            f = velocities[0]
            frame = velocities[1]
            vels = KFFile(f).read("MDHistory", f"Velocities({frame})")
            vels = np.array(vels).reshape(-1, 3) * Units.convert(1.0, "bohr", "angstrom")  # angstrom/fs
            s.input.ams.MolecularDynamics.InitialVelocities.Type = "Input"
            values = ""
            for x in vels:
                values += 6 * " " + " ".join(str(y) for y in x) + "\n"
            s.input.ams.MolecularDynamics.InitialVelocities.Values._h = "   # From {} frame {}".format(f, frame)
            s.input.ams.MolecularDynamics.InitialVelocities.Values._1 = values
        else:
            s.input.ams.MolecularDynamics.InitialVelocities.Type = "FromFile"
            s.input.ams.MolecularDynamics.InitialVelocities.File = velocities
        return s

    def get_velocities_from(self, other_job, frame=None, update_molecule=True):
        """
        Function to update the InitialVelocities block in self. It is normally not needed, instead use the e.g. AMSNVEJob.restart_from() function.

        This function can be called in prerun() methods for MultiJobs
        """
        _, velocities, molecule, _ = self._get_restart_job_velocities_molecule(other_job, frame=frame)
        del self.settings.input.ams.MolecularDynamics.InitialVelocities
        self.settings += self._velocities2settings(velocities)

        if update_molecule:
            self.molecule = molecule

    @staticmethod
    def _get_restart_job_velocities_molecule(other_job, frame=None, settings=None, get_velocities_molecule=True):
        """
        other_job: str or some AMSMdJob

        get_velocities_molecule : bool
            Whether to get the velocities and molecule right now. Set to False to not access other_job.results (for use with MultiJob prerun() methods)

        Returns: (other_job [AMSJob], velocities, molecule, extra_settings [Settings])
        """
        if isinstance(other_job, str):
            other_job = AMSJob.load_external(other_job)
        if get_velocities_molecule:
            if frame:
                velocities = (other_job.results.rkfpath(), frame)
                molecule = other_job.results.get_history_molecule(frame)
            else:
                velocities = other_job
                molecule = other_job.results.get_main_molecule()
        else:
            velocities = None
            molecule = None

        extra_settings = other_job.settings.copy()
        if settings:
            extra_settings.update(settings)

        if "InitialVelocities" in extra_settings.input.ams.MolecularDynamics:
            del extra_settings.input.ams.MolecularDynamics.InitialVelocities

        if "System" in extra_settings.input.ams:
            del extra_settings.input.ams.System

        return other_job, velocities, molecule, extra_settings

    def _get_thermostat_settings(self, thermostat, temperature, tau, thermostat_region, nsteps: int):
        s = Settings()
        prev_thermostat_settings = self.settings.input.ams.MolecularDynamics.Thermostat
        if isinstance(prev_thermostat_settings, list):
            prev_thermostat_settings = prev_thermostat_settings[0]
        s.input.ams.MolecularDynamics.Thermostat.Type = (
            thermostat or prev_thermostat_settings.Type or self.default_thermostat
        )
        try:
            n_temperatures = len(temperature)
            my_temperature = " ".join(str(x) for x in temperature)
            if n_temperatures > 1:
                nsteps_per_temperature = nsteps // (len(temperature) - 1)
                s.input.ams.MolecularDynamics.Thermostat.Duration = " ".join(
                    [str(nsteps_per_temperature)] * (n_temperatures - 1)
                )
            else:
                s.input.ams.MolecularDynamics.Thermostat.Duration = None
        except TypeError:
            my_temperature = temperature
            s.input.ams.MolecularDynamics.Thermostat.Duration = None

        s.input.ams.MolecularDynamics.Thermostat.Region = (
            thermostat_region or prev_thermostat_settings.get("Region", None) or self.default_thermostat_region
        )
        s.input.ams.MolecularDynamics.Thermostat.Temperature = (
            my_temperature or prev_thermostat_settings.Temperature or self.default_temperature
        )
        s.input.ams.MolecularDynamics.Thermostat.Tau = (
            tau
            or prev_thermostat_settings.Tau
            or float(self.settings.input.ams.MolecularDynamics.TimeStep) * AMSMDJob.default_tau_multiplier
        )
        return s

    def _get_barostat_settings(self, pressure, barostat, barostat_tau, scale, equal, constantvolume):
        s = Settings()
        self.settings.input.ams.MolecularDynamics.Barostat.Type = (
            barostat or self.settings.input.ams.MolecularDynamics.Barostat.Type or self.default_barostat
        )
        self.settings.input.ams.MolecularDynamics.Barostat.Pressure = (
            pressure
            if pressure is not None
            else self.settings.input.ams.MolecularDynamics.Barostat.Pressure or self.default_pressure
        )
        self.settings.input.ams.MolecularDynamics.Barostat.Tau = (
            barostat_tau
            or self.settings.input.ams.MolecularDynamics.Barostat.Tau
            or float(self.settings.input.ams.MolecularDynamics.TimeStep) * AMSMDJob.default_barostat_tau_multiplier
        )
        self.settings.input.ams.MolecularDynamics.Barostat.Scale = (
            scale or self.settings.input.ams.MolecularDynamics.Barostat.Scale or self.default_scale
        )
        self.settings.input.ams.MolecularDynamics.Barostat.Equal = (
            equal or self.settings.input.ams.MolecularDynamics.Barostat.Equal or self.default_equal
        )
        self.settings.input.ams.MolecularDynamics.Barostat.ConstantVolume = (
            str(constantvolume)
            if constantvolume is not None
            else self.settings.input.ams.MolecularDynamics.Barostat.ConstantVolume or self.default_constantvolume
        )
        return s

    @classmethod
    def restart_from(cls, other_job, frame=None, settings=None, use_prerun=False, **kwargs):
        """

        other_job: AMSJob
            The job to restart from.

        frame: int
            Which frame to read the structure and velocities from. If None, the final structure and end velocities will be used (section Molecule and MDResults%EndVelocities)

        settings: Settings
            Settings that override any other settings. All settings from other_job (e.g. the engine settings) are inherited by default but they can be overridden here.

        use_prerun: bool
            If True, the molecule and velocities will only be read from other_job inside the prerun() method. Set this to True to prevent PLAMS from waiting for other_job to finish as soon as the new job is defined.

        kwargs: many options
            See the docstring for AMSMDJob.

        """
        other_job, velocities, molecule, extra_settings = cls._get_restart_job_velocities_molecule(
            other_job, frame, settings, get_velocities_molecule=not use_prerun
        )
        job = cls(settings=extra_settings, velocities=velocities, molecule=molecule, **kwargs)

        if use_prerun:

            @add_to_instance(job)
            def prerun(self):  # noqa F811
                self.get_velocities_from(other_job, frame=frame, update_molecule=True)

        return job


class AMSNVEJob(AMSMDJob):
    """
    A class for running NVE MD simulations
    """

    def __init__(self, **kwargs):
        AMSMDJob.__init__(self, **kwargs)
        self.remove_blocks(["thermostat", "barostat", "deformation", "nanoreactor"])


class AMSNVTJob(AMSMDJob):
    """
    A class for running NVT MD simulations
    """

    def __init__(self, **kwargs):
        AMSMDJob.__init__(self, _enforce_thermostat=True, **kwargs)
        self.remove_blocks(["barostat", "deformation", "nanoreactor"])


class AMSNPTResults(AMSResults):
    def get_equilibrated_molecule(self, equilibration_fraction=0.667, return_index=False):
        """

        Discards the first equilibration_fraction of the trajectory.
        Calculates the average density of the rest. Returns the molecule
        with the closest density to the average density among the remaining
        trajectory.

        """
        densities = self.job.results.get_history_property("Density", "MDHistory")
        analyze_from = int(len(densities) * equilibration_fraction)
        # take structure closest to the target density
        avg_density = np.mean(densities[analyze_from:])  # amu/bohr^3
        delta = np.array(densities[analyze_from:]) - avg_density
        delta = np.abs(delta)
        min_index = np.argmin(delta)
        min_index += analyze_from
        mol = self.job.results.get_history_molecule(min_index + 1)

        if return_index:
            return mol, min_index + 1

        return mol


class AMSNPTJob(AMSMDJob):
    """
    A class for running NPT MD simulations
    """

    _result_type = AMSNPTResults

    def __init__(self, **kwargs):
        AMSMDJob.__init__(self, _enforce_thermostat=True, _enforce_barostat=True, **kwargs)
        self.settings.input.ams.MolecularDynamics.CalcPressure = "True"
        self.remove_blocks(["deformation", "nanoreactor"])

from ...interfaces.adfsuite.ams import AMSJob, AMSResults
from ...core.settings import Settings
from ...tools.kftools import KFFile
from ...tools.units import Units
from typing import Union
import numpy as np

__all__ = ['NVEJob', 'NVTJob', 'NPTJob']

class AMSMDJob(AMSJob):
    default_nsteps = 1000
    default_timestep = 0.25
    default_samplingfreq = 100

    default_thermostat = 'NHC'
    default_temperature = 300
    default_tau_multiplier = 400

    default_barostat = 'MTK'
    default_pressure = '1.0 [bar]'
    default_barostat_tau_multiplier = 4000
    default_scale = 'XYZ'
    default_equal = 'None'
    default_constantvolume = 'False'

    default_writevelocities = 'True'
    default_writebonds = 'True'
    default_writemolecules = 'True'
    default_writeenginegradients = 'False'
    default_calcpressure = 'False'

    def __init__(
        self, 
        velocities=None,
        timestep=None,
        samplingfreq=None,
        nsteps=None,
        checkpointfrequency=None,
        writevelocities=None,
        writebonds=None,
        writemolecules=None,
        writecharges=None,
        writeenginegradients=None,
        calcpressure=None,
        molecule=None,
        **kwargs
    ):
        if isinstance(molecule, AMSJob):
            molecule = molecule.results.get_main_molecule()
        if isinstance(molecule, AMSResults):
            molecule = molecule.get_main_molecule()
        AMSJob.__init__(self, molecule=molecule, **kwargs)

        self.settings.input.ams.Task = 'MolecularDynamics'
        mdsett = self.settings.input.ams.MolecularDynamics

        mdsett.TimeStep = timestep or mdsett.TimeStep or self.default_timestep
        mdsett.Trajectory.SamplingFreq = samplingfreq or mdsett.Trajectory.SamplingFreq or self.default_samplingfreq
        mdsett.NSteps = nsteps or mdsett.NSteps or self.default_nsteps

        mdsett.Trajectory.WriteVelocities = str(writevelocities) if writevelocities is not None else mdsett.Trajectory.WriteVelocities or self.default_writevelocities
        mdsett.Trajectory.WriteBonds = str(writebonds) if writebonds is not None else mdsett.Trajectory.WriteBonds or self.default_writebonds
        mdsett.Trajectory.WriteMolecules = str(writemolecules) if writemolecules is not None else mdsett.Trajectory.WriteMolecules or self.default_writemolecules
        mdsett.Trajectory.WriteEngineGradients = str(writeenginegradients) if writeenginegradients is not None else mdsett.Trajectory.WriteEngineGradients or self.default_writeenginegradients
        mdsett.CalcPressure = str(calcpressure) if calcpressure is not None else mdsett.CalcPressure or self.default_calcpressure
        self.settings += self._velocities2settings(velocities)

    def remove_blocks(self, blocks=None):
        if blocks:
            for block in blocks:
                if block in self.settings.input.ams.MolecularDynamics:
                    del self.settings.input.ams.MolecularDynamics[block]

    @staticmethod
    def _velocities2settings(velocities):
        s = Settings()
        if isinstance(velocities, int) or isinstance(velocities, float) or velocities is None or velocities is False:
            s.input.ams.MolecularDynamics.InitialVelocities.Type = 'Random'
            s.input.ams.MolecularDynamics.InitialVelocities.Temperature = velocities or AMSMDJob.default_temperature
        elif isinstance(velocities, tuple):
            # file and frame number
            f = velocities[0]
            frame = velocities[1]
            vels = KFFile(f).read('MDHistory', f'Velocities({frame})')
            vels = np.array(vels).reshape(-1,3) * Units.convert(1.0, 'bohr', 'angstrom') # angstrom/fs
            s.input.ams.MolecularDynamics.InitialVelocities.Type = 'Input'
            values = ""
            for x in vels:
                values += 6*" " + " ".join(str(y) for y in x) + "\n"
            s.input.ams.MolecularDynamics.InitialVelocities.Values._h = '   # From {} frame {}'.format(f, frame)
            s.input.ams.MolecularDynamics.InitialVelocities.Values._1 = values
        else:
            s.input.ams.MolecularDynamics.InitialVelocities.Type = 'FromFile'
            s.input.ams.MolecularDynamics.InitialVelocities.File = velocities
        return s


    @staticmethod
    def _get_restart_job_velocities_molecule(other_job, frame=None, settings=None):
        """
            other_job: str or some AMSMdJob

            Returns: (other_job [AMSJob], velocities, molecule, extra_settings [Settings])
        """
        if isinstance(other_job, str):
            other_job = AMSJob.load_external(other_job)
        if frame:
            velocities = (other_job.results.rkfpath(), frame)
            molecule = other_job.results.get_history_molecule(frame)
        else:
            velocities = other_job
            molecule = other_job.results.get_main_molecule()

        extra_settings = other_job.settings.copy()
        if settings:
            extra_settings.update(settings)

        if 'InitialVelocities' in extra_settings.input.ams.MolecularDynamics:
            del extra_settings.input.ams.MolecularDynamics.InitialVelocities

        if 'System' in extra_settings.input.ams:
            del extra_settings.input.ams.System

        return other_job, velocities, molecule, extra_settings

    def _get_thermostat_settings(self, thermostat, temperature, tau):
        s= Settings()
        s.input.ams.MolecularDynamics.Thermostat.Type = thermostat or self.settings.input.ams.MolecularDynamics.Thermostat.Type or self.default_thermostat
        s.input.ams.MolecularDynamics.Thermostat.Temperature = temperature or self.settings.input.ams.MolecularDynamics.Thermostat.Temperature or self.default_temperature
        s.input.ams.MolecularDynamics.Thermostat.Tau = tau or self.settings.input.ams.MolecularDynamics.Thermostat.Tau or float(self.settings.input.ams.MolecularDynamics.TimeStep) * AMSMDJob.default_tau_multiplier
        return s

    def _get_barostat_settings(self, pressure, barostat, barostat_tau, scale, equal, constantvolume):
        s = Settings()
        self.settings.input.ams.MolecularDynamics.Barostat.Type = barostat or self.settings.input.ams.MolecularDynamics.Barostat.Type or self.default_barostat
        self.settings.input.ams.MolecularDynamics.Barostat.Pressure = str(pressure) + ' [bar]' if pressure is not None else self.settings.input.ams.MolecularDynamics.Barostat.Pressure or self.default_pressure
        self.settings.input.ams.MolecularDynamics.Barostat.Tau = barostat_tau or self.settings.input.ams.MolecularDynamics.Barostat.Tau or float(self.settings.input.ams.MolecularDynamics.TimeStep) * AMSMDJob.default_barostat_tau_multiplier
        self.settings.input.ams.MolecularDynamics.Barostat.Scale = scale or self.settings.input.ams.MolecularDynamics.Barostat.Scale or self.default_scale
        self.settings.input.ams.MolecularDynamics.Barostat.Equal = equal or self.settings.input.ams.MolecularDynamics.Barostat.Equal or self.default_equal 
        self.settings.input.ams.MolecularDynamics.Barostat.ConstantVolume = str(constantvolume) if constantvolume is not None else self.settings.input.ams.MolecularDynamics.Barostat.ConstantVolume or self.default_constantvolume
        return s

class NVEJob(AMSMDJob):
    def __init__( self, **kwargs):
        AMSMDJob.__init__(self, **kwargs)
        self.remove_blocks(['thermostat', 'barostat', 'deformation'])

    @classmethod
    def restart_from(cls, 
        other_job, 
        frame=None, 
        settings=None, 
        timestep=None,
        samplingfreq=None,
        nsteps=None,
        checkpointfrequency=None,
        writevelocities=None,
        writebonds=None,
        writemolecules=None,
        writecharges=None,
        writeenginegradients=None,
        calcpressure=None,
        molecule=None,
        **kwargs
    ):
        other_job, velocities, molecule, extra_settings = cls._get_restart_job_velocities_molecule(other_job, frame, settings)
        return cls(
            settings=extra_settings, 
            velocities=velocities, 
            molecule=molecule, 
            timestep=timestep,
            samplingfreq=samplingfreq,
            nsteps=nsteps,
            checkpointfrequency=checkpointfrequency,
            writevelocities=writevelocities,
            writebonds=writebonds,
            writemolecules=writemolecules,
            writecharges=writecharges,
            writeenginegradients=writeenginegradients,
            calcpressure=calcpressure,
            **kwargs)


class NVTJob(AMSMDJob):

    def __init__(self, 
        temperature=None,
        velocities=None,
        thermostat=None,
        tau=None,
        **kwargs):
        AMSMDJob.__init__(self, velocities = velocities or temperature or AMSMDJob.default_temperature, **kwargs)
        self.settings.update(self._get_thermostat_settings(thermostat, temperature, tau))
        self.remove_blocks(['barostat', 'deformation'])


    @classmethod
    def restart_from(cls, other_job, 
        settings=None,
        temperature=None,
        thermostat=None,
        tau=None,
        frame=None,
        **kwargs):

        other_job, velocities, molecule, extra_settings = cls._get_restart_job_velocities_molecule(other_job, frame, settings)
        return cls(molecule=molecule, settings=extra_settings, velocities=velocities, thermostat=thermostat, temperature=temperature, tau=tau, **kwargs)


class NPTJob(NVTJob):
    def __init__(self,
        pressure=None,
        barostat=None,
        barostat_tau=None,
        scale=None,
        equal=None,
        constantvolume=None,
        velocities=None,
        thermostat=None,
        temperature=None,
        tau=None,
        **kwargs
    ):
        AMSMDJob.__init__(self, velocities=velocities or temperature or AMSMDJob.default_temperature, **kwargs)
        self.settings.update(self._get_thermostat_settings(thermostat=thermostat, temperature=temperature, tau=tau))
        self.settings.update(self._get_barostat_settings(
            pressure=pressure,
            barostat=barostat,
            barostat_tau=barostat_tau,
            scale=scale,
            equal=equal,
            constantvolume=constantvolume
        ))
        self.settings.input.ams.MolecularDynamics.CalcPressure = 'True'

        self.remove_blocks(['deformation'])

    @classmethod
    def restart_from(cls,
        other_job, 
        frame=None, 
        settings=None,
        pressure=None,
        barostat=None,
        barostat_tau=None,
        scale=None,
        equal=None,
        constantvolume=None,
        temperature=None,
        thermostat=None,
        tau=None,
        **kwargs
    ):
        
        other_job, velocities, molecule, extra_settings = cls._get_restart_job_velocities_molecule(other_job, frame, settings)

        return cls(
            molecule=molecule,
            velocities=velocities,
            settings=extra_settings,
            thermostat=thermostat,
            temperature=temperature,
            tau=tau,
            pressure=pressure,
            barostat=barostat,
            barostat_tau=barostat_tau,
            scale=scale,
            equal=equal,
            constantvolume=constantvolume,
            **kwargs
        )



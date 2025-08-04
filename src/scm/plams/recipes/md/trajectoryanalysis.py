import copy
from collections import OrderedDict, defaultdict
from collections.abc import Iterable
from typing import List

import numpy as np
from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.core.functions import requires_optional_package
from scm.plams.interfaces.adfsuite.amsanalysis import AMSAnalysisJob, AMSAnalysisResults
from scm.plams.tools.units import Units
from scm.plams.core.enums import JobStatus

__all__ = ["AMSRDFJob", "AMSMSDJob", "AMSMSDResults", "AMSVACFJob", "AMSVACFResults"]


class AMSConvenientAnalysisJob(AMSAnalysisJob):
    _task = "None"

    def __init__(self, previous_job=None, atom_indices=None, **kwargs):  # needs to be finished
        """
        previous_job: AMSJob
            An AMSJob with an MD trajectory. Note that the trajectory should have been equilibrated before it starts.

        All other settings can be set as for AMS

        """
        AMSAnalysisJob.__init__(self, **kwargs)
        if not self._has_settings_entry(self.settings.input, "Task"):
            self.settings.input.Task = self._task

        self.previous_job = previous_job
        self.atom_indices = atom_indices

    def _get_max_dt_frames(self, max_correlation_time_fs):
        if max_correlation_time_fs is None:
            return None

        # Adjust for multiple previous_jobs
        prevjobs = self.previous_job
        if not isinstance(prevjobs, list):
            prevjobs = [prevjobs]

        # Read history and time step
        # FIXME: Does not account for the presence of Range of StepSize (which the user may have provided)
        historylength = 0
        for prevjob in prevjobs:
            historylength += prevjob.results.readrkf("History", "nEntries")
        timestep = prevjobs[0].results.get_time_step()
        max_dt_frames = int(max_correlation_time_fs / timestep)
        max_dt_frames = min(max_dt_frames, historylength // 2)
        return max_dt_frames

    def get_input(self):
        """
        Generate the input file
        """
        self._settings_to_list(self.settings.input, self._task)
        if self.atom_indices and self._parent_write_atoms:
            section = getattr(self.settings.input, self._task)
            for entry in section:
                if not (self._has_settings_entry(entry, "Atoms") and self._has_settings_entry(entry.Atoms, "Atom")):
                    self._add_nonunique_settings_entries(entry.Atoms.Atom, self.atom_indices)
        return super().get_input()

    @staticmethod
    def _settings_to_list(settings, entry, nentries=1):
        """
        Convert the main settings object for this Task to a list of settings objects
        """
        if isinstance(settings, Settings):
            # This is not a PISA object
            if AMSConvenientAnalysisJob._has_settings_entry(settings, entry):
                # If this is not a list-type object, make it into one
                if isinstance(settings[entry], dict) or not (isinstance(settings[entry], Iterable)):
                    settings[entry] = [settings[entry]]
            else:
                settings[entry] = [Settings() for i in range(nentries)]
        else:
            # Workaround for a PISA bug, where the object can have updated values but length 0
            block = getattr(settings, entry)
            if len(block) == 0:
                if block.value_changed:
                    s = copy.deepcopy(block)
                    block[0] = s
                else:
                    # I need there to be a single entry
                    for i in range(nentries):
                        s = block[i]

    @staticmethod
    def _has_settings_entry(settings, entry):
        """
        Check if this Settings or PISA object contains this entry already
        """
        if isinstance(settings, Settings):
            return entry in settings
        else:
            return getattr(settings, entry).value_changed

    @staticmethod
    def _add_nonunique_settings_entries(settings, entries):
        """
        For non-default entries, multiple entries can be supplied

        * ``settings`` -- Settings or PISA object to which the entries should be added
        * ``entries`` -- Iterator for multiple settings entries
        """
        if not isinstance(settings, Settings):
            for i, entry in enumerate(entries):
                settings[i] = entry
        else:
            settings = entries

    def _parent_prerun(self):
        """
        Possibly wait for previous job to finish before adding path
        """
        prevjobs = self.previous_job
        if not isinstance(prevjobs, list):
            prevjobs = [prevjobs]

        # Turn trajectory settings into a list,
        # If no entries yet exist, it will be a list of empty entries
        self._settings_to_list(self.settings.input.TrajectoryInfo, "Trajectory", nentries=len(prevjobs))

        nblocks = len(self.settings.input.TrajectoryInfo.Trajectory)
        msg = "Number of Trajectory blocks must match number of previous_jobs"
        assert nblocks == len(prevjobs) or nblocks == 0, msg

        # Overwrite KFFilename with that of previous_job
        for i, prevjob in enumerate(prevjobs):
            assert prevjob.status != JobStatus.CREATED, "You can only pass in a finished AMSJob"
            self.settings.input.TrajectoryInfo.Trajectory[i].KFFilename = prevjob.results.rkfpath()


class AMSMSDResults(AMSAnalysisResults):
    """Results class for AMSMSDJob"""

    def get_msd(self):
        """
        returns time [fs], msd [ang^2]
        """
        msd_xy = self.get_xy()
        time = np.array(msd_xy.x[0])  # fs
        y = np.array(msd_xy.y) * Units.convert(1.0, "bohr", "angstrom") ** 2

        return time, y

    @requires_optional_package("scipy")
    def get_linear_fit(self, start_time_fit_fs=None, end_time_fit_fs=None):
        """
        Fits the MSD between start_time_fit_fs and end_time_fit_fs

        Returns a 3-tuple LinRegress.result, fit_x (fs), fit_y (ang^2)

        result.slope is given in ang^2/fs

        """
        from scipy.stats import linregress

        time, y = self.get_msd()
        end_time_fit_fs = end_time_fit_fs or max(time)
        start_time_fit_fs = start_time_fit_fs or self.job.start_time_fit_fs

        if start_time_fit_fs >= end_time_fit_fs:
            start_time_fit_fs = end_time_fit_fs / 2

        y = y[time >= start_time_fit_fs]
        time = time[time >= start_time_fit_fs]
        y = y[time <= end_time_fit_fs]
        time = time[time <= end_time_fit_fs]

        result = linregress(time, y)
        fit_x = time
        fit_y = result.slope * fit_x + result.intercept

        return result, fit_x, fit_y

    def get_diffusion_coefficient(self, start_time_fit_fs=None, end_time_fit_fs=None):
        """
        Returns D in m^2/s
        """
        result, _, _ = self.get_linear_fit(start_time_fit_fs=start_time_fit_fs, end_time_fit_fs=end_time_fit_fs)
        D = result.slope * 1e-20 / (6 * 1e-15)  # convert from ang^2/fs to m^2/s, divide by 6 because 3-dimensional (2d)
        return D


class AMSMSDJob(AMSConvenientAnalysisJob):
    """A convenient class wrapping around the trajectory analysis MSD tool"""

    _result_type = AMSMSDResults
    _parent_write_atoms = True
    _task = "MeanSquareDisplacement"

    def __init__(
        self,
        previous_job=None,  # needs to be finished
        max_correlation_time_fs: float = 10000,
        start_time_fit_fs: float = 2000,
        atom_indices: List[int] = None,
        **kwargs,
    ):
        """
        previous_job: AMSJob
            An AMSJob with an MD trajectory. Note that the trajectory should have been equilibrated before it starts.

        max_correlation_time_fs: float
            Maximum correlation time in femtoseconds

        start_time_fit_fs : float
            Smallest correlation time for the linear fit

        atom_indices : List[int]
            If None, use all atoms. Otherwise use the specified atom indices (starting with 1)

        kwargs: dict
            Other options to AMSAnalysisJob

        """
        AMSConvenientAnalysisJob.__init__(self, previous_job=previous_job, atom_indices=atom_indices, **kwargs)

        self.max_correlation_time_fs = max_correlation_time_fs
        self.start_time_fit_fs = start_time_fit_fs

    def get_input(self):
        """
        Generate the input file
        """
        self._settings_to_list(self.settings.input, self._task)

        for settings in self.settings.input.MeanSquareDisplacement:
            if not self._has_settings_entry(settings, "Property"):
                settings.Property = "DiffusionCoefficient"
            if not self._has_settings_entry(settings, "StartTimeSlope"):
                settings.StartTimeSlope = self.start_time_fit_fs
        return super().get_input()

    def prerun(self):  # noqa F811
        """
        Constructs the final settings
        """
        self._parent_prerun()  # trajectory and atom_indices handled
        self._settings_to_list(self.settings.input, self._task)

        for settings in self.settings.input.MeanSquareDisplacement:
            if not self._has_settings_entry(settings, "MaxFrame"):
                max_dt_frames = self._get_max_dt_frames(self.max_correlation_time_fs)
                settings.MaxFrame = max_dt_frames

    def postrun(self):
        """
        Creates msd.txt, fit_msd.txt, and D.txt
        """
        time, msd = self.results.get_msd()
        with open(self.path + "/msd.txt", "w") as f:
            f.write("#Time(fs) MSD(ang^2)")
            for x, y in zip(time, msd):
                f.write(f"{x} {y}\n")

        _, fit_x, fit_y = self.results.get_linear_fit()
        with open(self.path + "/fit_msd.txt", "w") as f:
            f.write("#Time(fs), LinearFitToMSD(ang^2)")
            for x, y in zip(fit_x, fit_y):
                f.write(f"{x} {y}\n")

        D = self.results.get_diffusion_coefficient()
        with open(self.path + "/D.txt", "w") as f:
            f.write(f"{D}\n")


def viscosity_double_exponential(x, A, lam, tau1, tau2):
    return A * (lam * (1 - np.exp(-x / tau1)) + (1 - lam) * (1 - np.exp(-x / tau2)))


class AMSViscosityFromBinLogResults(AMSAnalysisResults):
    """Results class for AMSViscosityFromBinLogJob"""

    def get_viscosity_integral(self):
        """Extract the running viscosity integral.

        Returns a 2-tuple with 1D numpy arrays ``time, viscosity_integral``, with time in fs and viscosity_integral in Pa*s.
        """

        xy = self.get_xy("Integral")
        time = np.array(xy.x[0])  # fs
        y = np.array(xy.y)  # Pa s

        return time, y

    @requires_optional_package("scipy")
    def get_double_exponential_fit(self):
        """Perform a double exponential fit to the viscosity integral.

        The fitted function is of the form

        ``A * (lam * (1 - np.exp(-x / tau1)) + (1 - lam) * (1 - np.exp(-x / tau2)))``

        where ``A`` is the limiting value in the infinite time limit.

        Returns: a 3-tuple ``popt, time, prediction``, where

        - ``popt`` is a 4-tuple containing A, lam, tau1, and tau2,

        - ``time`` is in fs, and

        - ``prediction`` is the value of the above function vs time.
        """

        from scipy.optimize import curve_fit

        x, y = self.get_viscosity_integral()

        # The fit
        f, p0 = viscosity_double_exponential, (y[-1], 0.5, 0.1 * x[-1], 0.9 * x[-1])
        popt, _ = curve_fit(f, x, y, p0=p0)
        prediction = f(x, *popt)

        return popt, x, prediction


class AMSViscosityFromBinLogJob(AMSConvenientAnalysisJob):
    """A convenient class wrapping around the trajectory analysis ViscosityFromBinLog tool. Only runs with default input options (max correlation time 10% of the total trajectory)."""

    _result_type = AMSViscosityFromBinLogResults
    _parent_write_atoms = False
    _task = "AutoCorrelation"

    def __init__(
        self,
        previous_job=None,  # needs to be finished
        **kwargs,
    ):
        """
        previous_job: AMSJob
            An AMSJob with an MD trajectory. Note that the trajectory should have been equilibrated before it starts. It must have been run with the BinLog%PressureTensor option set.

        kwargs: dict
            Other options to AMSAnalysisJob

        """
        AMSConvenientAnalysisJob.__init__(self, previous_job=previous_job, **kwargs)

    def get_input(self):
        """
        Generate the input file
        """
        self._settings_to_list(self.settings.input, self._task)

        for settings in self.settings.input.AutoCorrelation:
            settings.Property = "ViscosityFromBinLog"
        return super().get_input()

    def prerun(self):  # noqa F811
        """
        Constructs the final settings
        """
        self._parent_prerun()  # trajectory and atom_indices handled

    def postrun(self):
        """
        Creates viscosity_integral, fit_viscosity_integral, and viscosity.txt
        """
        time, integral = self.results.get_viscosity_integral()
        with open(self.path + "/viscosity_integral.txt", "w") as f:
            f.write("#Time(fs) Integral(Pa*s)")
            for x, y in zip(time, integral):
                f.write(f"{x} {y}\n")

        popt, fit_x, fit_y = self.results.get_double_exponential_fit()
        with open(self.path + "/fit_viscosity_integral.txt", "w") as f:
            f.write("#Time(fs), DoubleExponentialFitToViscosityIntegral(Pa*s)")
            for x, y in zip(fit_x, fit_y):
                f.write(f"{x} {y}\n")

        eta = popt[0]
        with open(self.path + "/viscosity.txt", "w") as f:
            f.write(f"{eta}\n")


class AMSVACFResults(AMSAnalysisResults):
    """Results class for AMSVACFJob"""

    def get_vacf(self):
        """Extract the velocity autocorrelation function vs. time.

        Returns a 2-tuple ``time, y`` with ``time`` in fs.
        """
        xy = self.get_xy()
        time = np.array(xy.x[0])  # fs
        y = np.array(xy.y)

        return time, y

    def get_power_spectrum(self, max_freq=None):
        """Calculate the power spectrum as the Fourier transform of the velocity autocorrelation function.

        max_freq: float, optional
            The maximum frequency in cm^-1.

        Returns: A 2-tuple ``freq, y`` with ``freq`` in cm^-1 and ``y`` unitless.
        """

        max_freq = max_freq or self.job.max_freq or 5000
        xy = self.get_xy("Spectrum")
        freq = np.array(xy.x[0])
        y = np.array(xy.y)

        y = y[freq < max_freq]
        freq = freq[freq < max_freq]

        return freq, y

    def get_acf(self):
        return self.get_vacf()

    def get_spectrum(self, max_freq=None):
        return self.get_power_spectrum(max_freq=max_freq)


class AMSVACFJob(AMSConvenientAnalysisJob):
    """A class for calculating the velocity autocorrelation function and its power spectrum"""

    _result_type = AMSVACFResults
    _parent_write_atoms = True
    _task = "AutoCorrelation"

    def __init__(
        self,
        previous_job=None,  # needs to be finished
        max_correlation_time_fs=5000,  # fs
        max_freq=5000,  # cm^-1
        atom_indices=None,
        **kwargs,
    ):
        """
        previous_job: AMSJob
            An AMSJob with an MD trajectory. Note that the trajectory should have been equilibrated before it starts.

        max_correlation_time_fs: float
            Maximum correlation time in femtoseconds

        max_freq: float
            The maximum frequency for the power spectrum in cm^-1

        atom_indices: List[int]
            Atom indices (starting with 1). If None, use all atoms.

        """
        AMSConvenientAnalysisJob.__init__(self, previous_job=previous_job, atom_indices=atom_indices, **kwargs)

        self.max_correlation_time_fs = max_correlation_time_fs
        self.max_freq = max_freq

    def get_input(self):
        """
        Generate the input file
        """
        self.settings_to_list()

        for settings in self.settings.input.AutoCorrelation:
            settings.Property = "Velocities"
        return super().get_input()

    def prerun(self):  # noqa F811
        """
        Creates final settings
        """
        self._parent_prerun()  # trajectory and atom_indices handled
        self.settings_to_list()

        for settings in self.settings.input.AutoCorrelation:
            if not self._has_settings_entry(settings, "MaxFrame"):
                max_dt_frames = self._get_max_dt_frames(self.max_correlation_time_fs)
                settings.MaxFrame = max_dt_frames

    def postrun(self):
        """
        Creates vacf.txt and power_spectrum.txt
        """
        try:
            time, vacf = self.results.get_vacf()
            with open(self.path + "/vacf.txt", "w") as f:
                f.write("#Time(fs) VACF")
                for x, y in zip(time, vacf):
                    f.write(f"{x} {y}\n")

            freq, intens = self.results.get_power_spectrum()
            with open(self.path + "/power_spectrum.txt", "w") as f:
                f.write("#Frequency(cm^-1) Intensity(arb.units)")
                for x, y in zip(freq, intens):
                    f.write(f"{x} {y}\n")
        except:
            pass


class AMSDipoleDerivativeACFResults(AMSAnalysisResults):
    """Results class for AMSVACFJob"""

    def get_acf(self):
        xy = self.get_xy()
        time = np.array(xy.x[0])  # fs
        y = np.array(xy.y)

        return time, y

    def get_ir_spectrum(self, max_freq=None):
        max_freq = max_freq or self.job.max_freq or 5000
        xy = self.get_xy("Spectrum")
        freq = np.array(xy.x[0])
        y = np.array(xy.y)

        y = y[freq < max_freq]
        freq = freq[freq < max_freq]

        return freq, y

    def get_spectrum(self, max_freq=None):
        return self.get_ir_spectrum(max_freq=max_freq)


class AMSDipoleDerivativeACFJob(AMSConvenientAnalysisJob):
    """A class for calculating the velocity autocorrelation function and its power spectrum"""

    _result_type = AMSDipoleDerivativeACFResults
    _parent_write_atoms = True
    _task = "AutoCorrelation"

    def __init__(
        self,
        previous_job=None,  # needs to be finished
        max_correlation_time_fs=5000,  # fs
        max_freq=5000,  # cm^-1
        atom_indices=None,
        **kwargs,
    ):
        """
        previous_job: AMSJob
            An AMSJob with an MD trajectory. Note that the trajectory should have been equilibrated before it starts.

        max_correlation_time_fs: float
            Maximum correlation time in femtoseconds

        max_freq: float
            The maximum frequency for the power spectrum in cm^-1

        atom_indices: List[int]
            Atom indices (starting with 1). If None, use all atoms.

        """
        AMSConvenientAnalysisJob.__init__(self, previous_job=previous_job, atom_indices=atom_indices, **kwargs)

        self.max_correlation_time_fs = max_correlation_time_fs
        self.max_freq = max_freq

    def get_input(self):
        """
        Generate the input file
        """
        self.settings_to_list()

        for settings in self.settings.input.AutoCorrelation:
            settings.Property = "DipoleDerivativeFromCharges"
        return super().get_input()

    def prerun(self):  # noqa F811
        """
        Creates final settings
        """
        self._parent_prerun()  # trajectory and atom_indices handled
        self.settings_to_list()

        for settings in self.settings.input.AutoCorrelation:
            if not self._has_settings_entry(settings, "MaxFrame"):
                max_dt_frames = self._get_max_dt_frames(self.max_correlation_time_fs)
                settings.MaxFrame = max_dt_frames

    def postrun(self):
        """
        Creates dipolederivativeacf.txt and ir_spectrum.txt
        """
        try:
            time, acf = self.results.get_acf()
            with open(self.path + "/dipolederivativeacf.txt", "w") as f:
                f.write("#Time(fs) DipoleDerivativeACF")
                for x, y in zip(time, acf):
                    f.write(f"{x} {y}\n")

            freq, intens = self.results.get_ir_spectrum()
            with open(self.path + "/ir_spectrum.txt", "w") as f:
                f.write("#Frequency(cm^-1) Intensity(arb.units)")
                for x, y in zip(freq, intens):
                    f.write(f"{x} {y}\n")
        except:
            pass


class AMSRDFResults(AMSAnalysisResults):
    """Results class for AMSRDFJob"""

    def get_rdf(self):
        """
        Returns a 2-tuple r, rdf.

        r: numpy array of float (angstrom)
        rdf: numpy array of float
        """
        xy = self.get_xy()
        r = np.array(xy.x[0])
        y = np.array(xy.y)

        return r, y


class AMSRDFJob(AMSConvenientAnalysisJob):
    _result_type = AMSRDFResults
    _parent_write_atoms = False
    _task = "RadialDistribution"

    def __init__(
        self, previous_job=None, atom_indices=None, atom_indices_to=None, rmin=0.5, rmax=6.0, rstep=0.1, **kwargs
    ):
        """
        previous_job: AMSJob
            AMSJob with finished MD trajectory.

        atom_indices: List[int]
            Atom indices (starting with 1). If None, calculate RDF *from* all atoms.

        atom_indices_to: List[int]
            Atom indices (starting with 1). If None, calculate RDF *to* all atoms.

        rmin: float
            Minimum distance (angstrom)

        rmax: float
            Maximum distance (angstrom)

        rstep: float
            Bin width for the histogram (angstrom)
        """

        AMSConvenientAnalysisJob.__init__(self, previous_job=previous_job, atom_indices=atom_indices, **kwargs)
        self.atom_indices_to = atom_indices_to
        self.rmin = rmin
        self.rmax = rmax
        self.rstep = rstep

    def get_input(self):
        """
        Generate the input file
        """
        self._settings_to_list(self.settings.input, self._task)

        for settings in self.settings.input.RadialDistribution:
            if not self._has_settings_entry(settings, "AtomsFrom"):
                if self.atom_indices:
                    settings.AtomsFrom.Atom = self.atom_indices
            if not self._has_settings_entry(settings, "AtomsTo"):
                if self.atom_indices_to:
                    settings.Atom = self.atom_indices_to
            if not self._has_settings_entry(settings, "Range"):
                settings.Range = f"{self.rmin} {self.rmax} {self.rstep}"
        return super().get_input()

    def prerun(self):  # noqa F811
        """
        Creates the final settings. Do not call or override this method.
        """
        self._parent_prerun()
        self._settings_to_list(self.settings.input, self._task)

        prevjobs = self.previous_jobs
        if not isinstance(prevjobs, list):
            prevjobs = [prevjobs]
        main_mol = prevjobs[0].results.get_main_molecule()

        if not self.atom_indices:
            self.atom_indices = list(range(1, len(main_mol) + 1))
        if not self.atom_indices_to:
            self.atom_indices_to = list(range(1, len(main_mol) + 1))
        # Settings will be adjusted when get_input is called (again).

    def postrun(self):
        """
        Creates rdf.txt. Do not call or override this method.
        """
        try:
            r, gr = self.results.get_rdf()
            with open(self.path + "/rdf.txt", "w") as f:
                f.write("#r(angstrom) g(r)")
                for x, y in zip(r, gr):
                    f.write(f"{x} {y}\n")
        except:
            pass


class AMSConvenientAnalysisPerRegionResults(Results):
    def _getter(self, analysis_job_type, method, kwargs):
        assert (
            self.job.analysis_job_type is analysis_job_type
        ), f"{method}() can only be called for {analysis_job_type}, tried for type {self.job.analysis_job_type}"
        ret = {}
        for name, job in self.job.children.items():
            ret[name] = getattr(job.results, method)(**kwargs)
        return ret

    def get_diffusion_coefficient(self, **kwargs):
        return self._getter(AMSMSDJob, "get_diffusion_coefficient", kwargs)

    def get_msd(self, **kwargs):
        return self._getter(AMSMSDJob, "get_msd", kwargs)

    def get_linear_fit(self, **kwargs):
        return self._getter(AMSMSDJob, "get_linear_fit", kwargs)

    def get_vacf(self, **kwargs):
        return self._getter(AMSVACFJob, "get_vacf", kwargs)

    def get_power_spectrum(self, **kwargs):
        return self._getter(AMSVACFJob, "get_power_spectrum", kwargs)


class AMSConvenientAnalysisPerRegionJob(MultiJob):
    _result_type = AMSConvenientAnalysisPerRegionResults

    def __init__(self, previous_job, analysis_job_type, name=None, regions=None, per_element=False, **kwargs):
        MultiJob.__init__(self, children=OrderedDict(), name=name or "analysis_per_region")
        self.previous_job = previous_job
        self.analysis_job_type = analysis_job_type
        self.analysis_job_kwargs = kwargs
        self.regions_dict = regions
        self.per_element = per_element

    @staticmethod
    def get_regions_dict(molecule, per_element: bool = False):
        regions_dict = defaultdict(lambda: [])
        for i, at in enumerate(molecule, 1):
            regions = set([at.properties.region]) if isinstance(at.properties.region, str) else at.properties.region
            if len(regions) == 0:
                region_name = "NoRegion" if not per_element else f"NoRegion_{at.symbol}"
                regions_dict[region_name].append(i)
            for region in regions:
                region_name = region if not per_element else f"{region}_{at.symbol}"
                regions_dict[region_name].append(i)
            regions_dict["All"].append(i)
            if per_element:
                regions_dict[f"All_{at.symbol}"].append(i)

        return regions_dict

    def prerun(self):  # noqa F811
        regions_dict = self.regions_dict or self.get_regions_dict(
            self.previous_job.results.get_main_molecule(), per_element=self.per_element
        )

        for region in regions_dict:
            self.children[region] = self.analysis_job_type(
                previous_job=self.previous_job,
                name=region,
                atom_indices=regions_dict[region],
                **self.analysis_job_kwargs,
            )

    @staticmethod
    def get_mean_std_per_region(list_of_jobs, function_name, **kwargs):
        """
        list_of_jobs: List[AMSConvenientAnalysisPerRegionJob]
            List of jobs over which to average

        function_name: str
            e.g. 'get_msd', 'get_power_spectrum'

        """

        if isinstance(list_of_jobs, dict):
            list_of_jobs = [x for x in list_of_jobs.values()]

        all_x = defaultdict(lambda: [])
        all_y = defaultdict(lambda: [])
        for vacfjob in list_of_jobs:
            # let's calculate mean and std for each region
            results_dict = getattr(vacfjob.results, function_name)(**kwargs)
            for region_name, ret_value in results_dict.items():
                if np.isscalar(ret_value):
                    all_x[region_name].append(np.atleast_1d(ret_value))
                elif len(ret_value) == 2:
                    all_x[region_name].append(np.atleast_1d(ret_value[0]))
                    all_y[region_name].append(np.atleast_1d(ret_value[1]))

        mean_x = {}
        std_x = {}
        mean_y = {}
        std_y = {}
        for region_name, values in all_x.items():
            mean_x[region_name] = np.mean(values, axis=0)
            std_x[region_name] = np.std(values, axis=0)

            if len(all_y) > 0:
                mean_y[region_name] = np.mean(all_y[region_name], axis=0)
                std_y[region_name] = np.std(all_y[region_name], axis=0)

        if len(mean_y) > 0:
            return mean_x, std_x, mean_y, std_y
        else:
            return mean_x, std_x

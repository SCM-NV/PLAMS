"""Implementation of the AMSPipeCalculator class.

"""

from copy import deepcopy

import numpy as np

from scm.plams.core.functions import config
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.amsworker import AMSWorker
from scm.plams.interfaces.molecule.ase import fromASE, toASE

__all__ = ["AMSCalculator", "BasePropertyExtractor"]

try:
    from ase.calculators.calculator import Calculator, all_changes
    from ase.units import Bohr, Hartree
except ImportError:
    # empty interface if ase does not exist:
    __all__ = []

    class Calculator:
        def __init__(self, *args, **kwargs):
            raise NotImplementedError("AMSCalculator can not be used without ASE")

    all_changes = []
    Hartree = Bohr = 0


class BasePropertyExtractor:
    def __call__(self, ams_results, atoms):
        return self.extract(ams_results, atoms)

    def extract(self, ams_result, atoms):
        raise NotImplementedError

    def set_settings(self, settings):
        pass

    def check_settings(self, settings):
        return True


class EnergyExtractor(BasePropertyExtractor):
    name = "energy"

    def extract(self, ams_results, atoms):
        return ams_results.get_energy() * Hartree


def canonicalize_string(possible_string):
    try:
        return possible_string.lower().strip()
    except (AttributeError, TypeError):
        return possible_string


def is_ams_true(value):
    return canonicalize_string(value) in [True, "true", "yes"]


def is_ams_false(value):
    return not is_ams_true(value)


class ForceExtractor(BasePropertyExtractor):
    name = "forces"

    def extract(self, ams_results, atoms):
        return -ams_results.get_gradients() * Hartree / Bohr

    def set_settings(self, settings):
        settings.input.ams.Properties.Gradients = "Yes"
        return settings

    def check_settings(self, settings):
        return is_ams_true(settings.copy().input.ams.Properties.Gradients)


class StressExtractor(BasePropertyExtractor):
    name = "stress"

    def extract(self, ams_results, atoms):
        D = sum(atoms.get_pbc())
        if isinstance(atoms.get_pbc(), bool):
            D *= 3
        st = ams_results.get_stresstensor() * Hartree / Bohr**D
        xx, yy, zz, yz, xz, xy = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        if D >= 1:
            xx = st[0][0]
            if D >= 2:
                yy = st[1][1]
                xy = st[0][1]
                if D >= 3:
                    zz = st[2][2]
                    yz = st[1][2]
                    xz = st[0][2]
        return np.array([xx, yy, zz, yz, xz, xy])

    def set_settings(self, settings):
        settings.input.ams.Properties.StressTensor = "Yes"
        return settings

    def check_settings(self, settings):
        return is_ams_true(settings.copy().input.ams.Properties.StressTensor)


class AMSCalculator(Calculator):
    """
    ASE Calculator which can run any AMS engine (ADF, BAND, DFTB, ReaxFF, MLPotential, ForceField, ...).

    The settings are specified with a PLAMS ``Settings`` object in the same way as when running AMS through PLAMS.

    .. important::

        Before initializing the AMSCalculator you need to call ``plams.init()``:

        .. code-block:: python

            from scm.plams import *
            init()


    Parameters:

    settings  : Settings
                A Settings object representing the input for an AMSJob or AMSWorker.
                This also determines which `implemented_properties` are available:
                `settings.input.ams.properties.gradients`: `force`
                `settings.input.ams.properties.stresstensor`: `stress`
    name      : str, optional
                Name of the rundir of calculations done by this calculator. A counter
                is appended to the name for every calculation.
    amsworker : bool , optional
                If True, use the AMSWorker to set up an interactive session.
                The AMSWorker will spawn a seperate
                process (an amsdriver). In order to make sure this process is closed,
                either use AMSCalculator as a context manager or ensure that
                AMSCalculator.stop_worker() is called before python is finished:

                .. code-block:: python

                    with AMSCalculator(settings=settings, amsworker=True) as calc:
                        atoms.set_calculator(calc)
                        atoms.get_potential_energy()

                If False, use AMSJob to set up an io session (a normal AMS calculation storing all output on disk).
    restart   : bool , optional
                Allow the engine to restart based on previous calculations.
    molecule  : Molecule , optional
                A Molecule object for which the calculation has to be performed. If
                settings.input.ams.system is defined it overrides the molecule argument.
                If AMSCalculator.calculate(atoms = atoms) is called with an atoms argument
                it overrides any earlier definition of the system and remembers it.
    extractors: List[BasePropertyExtractor] , optional
                Define extractors for additional properties.


    Examples:
    """

    # counters are a dict as a class variable. This is to support deepcopying/multiple instances with the same name
    _counter = {}

    def __new__(cls, settings=None, name="", amsworker=False, restart=True, molecule=None, extractors=[]):
        """Dispatch object creation to AMSPipeCalculator or AMSJobCalculator depending on |amsworker|"""
        if cls == AMSCalculator:
            if amsworker:
                obj = object.__new__(AMSPipeCalculator)
            else:
                obj = object.__new__(AMSJobCalculator)
        else:
            obj = object.__new__(cls)
        return obj

    def __init__(self, settings=None, name="", amsworker=False, restart=True, molecule=None, extractors=[]):

        if not isinstance(settings, Settings):
            settings = Settings.from_dict(settings)
        else:
            settings = settings.copy()

        self.settings = settings.copy()
        self.amsworker = amsworker
        self._name = name
        self.restart = restart
        self.molecule = molecule
        extractors = settings.pop("Extractors", [])
        self.extractors = [EnergyExtractor(), ForceExtractor(), StressExtractor()]
        self.extractors += [e for e in extractors if not e in self.extractors]

        if "system" in self.settings.input.ams:
            mol_dict = AMSJob.settings_to_mol(settings)
            atoms = toASE(mol_dict[""]) if mol_dict else None
        elif molecule:
            atoms = toASE(molecule)
        else:
            atoms = None

        super().__init__()
        self.atoms = atoms

        self.prev_ams_results = None
        self.results = dict()
        self.properties_updated = False

    @property
    def counter(self):
        # this is needed for deepcopy/pickling etc
        if self.name not in self._counter:
            self.set_counter()
        self._counter[self.name] += 1
        return self._counter[self.name]

    def set_counter(self, value=0):
        self._counter[self.name] = value

    @property
    def implemented_properties(self):
        """Returns the list of properties that this calculator has implemented"""
        return [extractor.name for extractor in self.extractors if extractor.check_settings(self.settings)]

    @property
    def name(self):
        return self._name

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        """Calculate the requested properties. If atoms is not set, it will reuse the last known Atoms object."""
        if atoms is not None:
            # no need to redo the calculation, we already have everything.
            if (
                self.atoms == atoms
                and system_changes == []
                and all([p in self.results for p in properties])
                and not self.properties_updated
            ):
                return
            self.atoms = atoms.copy()
        self.properties_updated = False
        if self.atoms is None:
            raise ValueError("No atoms object was set.")

        if not config.init:
            raise RuntimeError("Before AMSCalculator can calculate results you need to call plams.init()")

        molecule = fromASE(self.atoms, set_charge=True)
        ams_results = self._get_ams_results(molecule, properties)
        if not ams_results.ok():
            self.results = dict()
            return

        self.results_from_ams_results(ams_results, self._get_job_settings(properties))
        self.prev_ams_results = ams_results

    def ensure_property(self, properties):
        """A list of ASE properties that the calculator will ensure are available from AMS or it gives an error."""
        if isinstance(properties, str):
            properties = [properties]
        for prop in properties:
            property_found = False
            for extractor in self.extractors:
                if prop == extractor.name:
                    self.settings = extractor.set_settings(self.settings.copy())
                    self.properties_updated = True
                property_found = True
            if not property_found:
                raise NotImplementedError(f"No extractor known for property {prop}")

    def results_from_ams_results(self, ams_results, job_settings):
        """Populates the self.results dictionary by having extractors act on an AMSResults object."""
        for extractor in self.extractors:
            if extractor.check_settings(job_settings):
                self.results[extractor.name] = extractor.extract(ams_results, self.atoms)

    def _get_ams_results(self, molecule, properties):
        raise NotImplementedError("Subclasses of AMSCalculator should implement this")

    def _get_job_settings(self, properties):
        """Returns a Settings object which ensures that an AMS calculation is run from which all requested
        properties can be extracted"""
        return self.settings.copy()

    def stop_worker(self):
        """Stops the amsworker if it exists"""
        if hasattr(self, "worker") and self.worker:
            stop = self.worker.stop()
            self.worker = None
            return stop
        else:
            # this is what AMSWorker.stop() would return if it was already stopped previously
            self.worker = None
            return (None, None)

    def clean_exit(self):
        """Function called by ASEPipeWorker to tell the Calculator to stop and clean up"""
        self.stop_worker()

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.stop_worker()

    @property
    def amsresults(self):
        if hasattr(self, "prev_ams_results"):
            return self.prev_ams_results


class AMSPipeCalculator(AMSCalculator):
    """This class should be instantiated through AMSCalculator with settings.Calculator.Pipe defined"""

    def __init__(self, settings=None, name="", amsworker=True, restart=True, molecule=None, extractors=[]):
        super().__init__(settings, name, amsworker, restart, molecule, extractors)

        self.worker_settings = self.settings.copy()
        if "Task" in self.worker_settings.input.ams:
            del self.worker_settings.input.ams.Task
        if "Properties" in self.worker_settings.input.ams:
            del self.worker_settings.input.ams.Properties
        self.worker = None

    def _get_ams_results(self, molecule, properties):
        job_settings = self._get_job_settings(properties)
        if self.worker is None:
            self.worker = AMSWorker(self.worker_settings, use_restart_cache=self.restart)
        # AMSWorker expects no engine definition at this point.
        s = Settings()
        s.input.ams = job_settings.input.ams
        if "amsworker" in job_settings:
            s.amsworker = job_settings.amsworker
        s.amsworker.prev_results = self.prev_ams_results
        job_settings = s
        # run the worker from _solve_from_settings
        return self.worker._solve_from_settings(
            name=self.name + str(self.counter), molecule=molecule, settings=job_settings
        )

    def __deepcopy__(self, memo):
        """The AMSWorker instance is not copied, but instead, all the copies use the same worker"""
        memo[id(self.worker)] = self.worker
        try:
            this_method = self.__deepcopy__
            self.__deepcopy__ = None
            copy = deepcopy(self, memo)
            self.__deepcopy__ = this_method
            return copy
        except Exception as e:
            self.__deepcopy__ = this_method
            raise e


class AMSJobCalculator(AMSCalculator):
    """This class should be instantiated through AMSCalculator with settings.Calculator.AMSJob defined"""

    def _get_ams_results(self, molecule, properties):
        job_settings = self._get_job_settings(properties)
        if self.restart and self.prev_ams_results:
            job_settings.input.ams.EngineRestart = self.prev_ams_results.rkfpath(file="engine")

        return AMSJob(name=self.name + str(self.counter), molecule=molecule, settings=job_settings).run()

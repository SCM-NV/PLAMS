from ..core.results import Results
from ..core.basejob import MultiJob
from ..core.errors import PlamsError
from ..interfaces.adfsuite.adf import ADFJob
from numpy import array as npa
from pickle import dump as pickleDump
from os.path import join as osPJ
from os.path import basename as osPB

__all__ = ['VibrationsResults' ,'VibrationsJob', 'IRJob']

try:
    from ase.vibrations import Vibrations as aseVib
    from ase.vibrations import Infrared as aseIR
    from ..interfaces.molecule.ase import *
except ImportError:
    __all__ = []


class VibrationsResults(Results):
    """
    A Class for handling numerical Vibrational analysis Results.

    Use ``get_ASEVib()`` to access the ``ase.vibrations.Vibrations`` object. Using e.g. ``self.get_ASEVib().summary()`` the results of the ASE frequency calculation are accessible.
    """
    def get_ASEVib(self):
        """Returns a ``ase.vibrations.Vibrations`` object if available"""
        if hasattr(self, '_vib'):
            return self._vib
        else:
            raise PlamsError('WARNING: No aseVib object present, probably Job %s failed?!' % self.job.name)


class VibrationsJob(MultiJob):
    """
    Class for calculating numerical Frequencies in parallel using arbitrary interfaces.
    This is achieved using the Vibrations class from ASE.

    *   ``name`` -- Name of the |MultiJob|
    *   ``molecule`` -- |Molecule| object (most propably in the chosen Methods optimized state)
    *   ``settings`` -- |Settings| instance for all Single-Point jobs to be run. Don't forget reference to a restart file if you want to save a lot of computer time!
    *   ``jobType`` -- |Job| Class you want to use.
    *   ``get_gradients`` -- Function name to retrieve gradients of the |Results| object of your chosen ``jobType.results``. Must take options ``energy_unit='eV'`` and ``dist_unit='Angstrom'``.
    *   ``aseVibOpt`` -- Options for ``ase.vibrations.Vibrations.__init__``

    The ``self.__init__()`` method calls on ase.vibrations.Vibrations.run to create all displacements, then adds them as children.
    After running, the ``self.postrun()`` method reads the gradients from all single-points and saves them according to the ASE scheme.
    Finally the ASE Vibrations object is made accessible in ``self.results.get_ASEVib`` for the retrieval of the ASE results.

    """
    _result_type = VibrationsResults

    def __init__(self, molecule, settings, jobType=ADFJob, get_gradients='get_gradients', aseVibOpt={}, name='plams.vib'):
        super().__init__(name=name)
        self.molecule = molecule
        self.settings = settings
        self.jobType = jobType
        self.aseVibOpt = aseVibOpt
        self.vibClass = aseVib

        self.get_grad = getattr(self.jobType._result_type, get_gradients)


    def new_children(self):
        #Create displaced molecules with ASE
        #Include our workingdir in the name, this results in the ASE .pckl files to be saved there
        #kinda hacky but works

        if len(self.children) == 0:
            add = []
            path = self.path
            self._vib = self.vibClass(toASE(self.molecule), name=osPJ(path,self.name), **self.aseVibOpt)
            for name, atoms in self._vib.iterdisplace():
                add.append(self.jobType(molecule=fromASE(atoms, properties=self.molecule.properties),
                            name=osPB(name), settings=self.settings))
            return add

        else:
            return None


    def postrun(self):
        #get gradients and save them to the pickle files for ASE
        path = self.path
        forces = {}
        for child in self.children:
            res = child.results
            name = child.name
            # don't rely on gradients beeing numpy arrays
            f = [ npa(vec) for vec in self.get_grad(res, energy_unit='eV', dist_unit='Angstrom') ]
            forces[name+'.pckl'] = -1.0 * npa(f)
        filename = osPJ(path, self.name+'.all.pckl')
        with open(filename, 'wb') as f:
            pickleDump(forces, f, protocol=2)

        self.results._vib = self._vib



class IRJob(VibrationsJob):
    """
    Subclass of |VibrationsJob| to calculate IR modes and intensities.

    Usage is the same as for the parent class, see |VibrationsJob|.

    Additional arguments:

    *   ``get_dipole_vector`` -- Function name to retrieve dipole vector of the |Results| object. Must take argument ``unit='au'``.
    """
    def __init__(self, molecule, settings, get_dipole_vector='get_dipole_vector', **kwargs):
        super().__init__(molecule, settings, **kwargs)
        self.get_dipole = getattr(self.jobType._result_type, get_dipole_vector)
        self.vibClass = aseIR



    def postrun(self):
        #get gradients and dipole and save them to the pickle files for ASE
        path = self.path
        save = {}
        for child in self.children:
            res = child.results
            name = child.name
            # don't rely on gradients beeing numpy arrays
            f = [ npa(vec) for vec in self.get_grad(res, energy_unit='eV', dist_unit='Angstrom') ]
            force = -1.0 * npa(f)
            dipole = npa(self.get_dipole(res, unit='au'))
            save[name+'.pckl'] = [force, dipole]
        filename = osPJ(path, self.name+'.all.pckl')
        with open(filename, 'wb') as f:
            pickleDump(save, f, protocol=2)

        self.results._vib = self._vib

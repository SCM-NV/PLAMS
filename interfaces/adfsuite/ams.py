import os
import numpy as np

from os.path import join as opj

from ...core.basejob import SingleJob
from ...core.errors import FileError, JobError, ResultsError, PTError
from ...core.functions import config, log, parse_heredoc
from ...core.private import sha256, UpdateSysPath
from ...core.results import Results
from ...core.settings import Settings, ig
from ...mol.molecule import Molecule
from ...mol.atom import Atom
from ...tools.kftools import KFFile
from ...tools.units import Units


__all__ = ['AMSJob', 'AMSResults']


class AMSResults(Results):
    """A specialized |Results| subclass for accessing the results of |AMSJob|."""

    def __init__(self, *args, **kwargs):
        Results.__init__(self, *args, **kwargs)
        self.rkfs = {}


    def collect(self):
        """Collect files present in the job folder. Use parent method from |Results|, then create an instance of |KFFile| for each ``.rkf`` file present in the job folder. Collect these files in ``rkfs`` dictionary, with keys being filenames without ``.rkf`` extension.

        The information about ``.rkf`` files generated by engines is taken from the main ``ams.rkf`` file.

        This method is called automatically during the final part of the job execution and there is no need to call it manually.
        """
        Results.collect(self)

        rkfname = 'ams.rkf'
        if rkfname in self.files:
            main = KFFile(opj(self.job.path, rkfname))
            n = main[('EngineResults','nEntries')]
            for i in range(1, n+1):
                title =  main[('EngineResults','Title({})'.format(i))]
                files =  main[('EngineResults','Files({})'.format(i))].splitlines()
                if files[0].endswith('.rkf'):
                    key = files[0][:-4]
                    self.rkfs[key] = KFFile(opj(self.job.path, files[0]))
            self.rkfs['ams'] = main

        else:
            log('WARNING: Main KF file {} not present in {}'.format(rkfname, self.job.path), 1)


    def refresh(self):
        """Refresh the contents of ``files`` list.

        Use the parent method from |Results|, then look at |KFFile| instances present in ``rkfs`` dictionary and check if they point to existing files. If not, try to reinstantiate them with current job path (that can happen while loading a pickled job after the entire job folder was moved).
        """
        Results.refresh(self)
        to_remove = []
        for key,val in self.rkfs.items():
            if not os.path.isfile(val.path):
                if os.path.dirname(val.path) != self.job.path:
                    guessnewpath = opj(self.job.path, os.path.basename(val.path))
                    if os.path.isfile(guessnewpath):
                        self.rkfs[key] = KFFile(guessnewpath)
                    else:
                        to_remove.append(key)
                else:
                    to_remove.append(key)
        for i in to_remove:
            del self.rkfs[i]


    def engine_names(self):
        """Return a list of all names of engine specific ``.rkf`` files. The identifier of the main result file (``'ams'``) is not present in the returned list, only engine specific names are listed.
        """
        self.refresh()
        ret = list(self.rkfs.keys())
        ret.remove('ams')
        return ret


    def rkfpath(self, file='ams'):
        """Return the absolute path of a chosen ``.rkf`` file.

        The *file* argument should be the identifier of the file to read. It defaults to ``'ams'``. To access a file called ``something.rkf`` you need to call this function with ``file='something'``. If there exists only one engine results ``.rkf`` file, you can call this function with ``file='engine'`` to access this file.
        """
        return self._access_rkf(lambda x: x.path, file)


    def readrkf(self, section, variable, file='ams'):
        """Read data from *section*/*variable* of a chosen ``.rkf`` file.

        The *file* argument should be the identifier of the file to read. It defaults to ``'ams'``. To access a file called ``something.rkf`` you need to call this function with ``file='something'``. If there exists only one engine results ``.rkf`` file, you can call this function with ``file='engine'`` to access this file.

        The type of the returned value depends on the type of *variable* defined inside KF file. It can be: single int, list of ints, single float, list of floats, single boolean, list of booleans or string.

        .. note::

            If arguments *section* or *variable* are incorrect (not present in the chosen file), the returned value is ``None``. Please mind the fact that KF files are case sensitive.

        """
        return self._access_rkf(lambda x: x.read(section, variable), file)


    def read_rkf_section(self, section, file='ams'):
        """Return a dictionary with all variables from a given *section* of a chosen ``.rkf`` file.

        The *file* argument should be the identifier of the file to read. It defaults to ``'ams'``. To access a file called ``something.rkf`` you need to call this function with ``file='something'``. If there exists only one engine results ``.rkf`` file, you can call this function with ``file='engine'`` to access this file.

        .. note::

            If *section* is not present in the chosen file, the returned value is an empty dictionary. Please mind the fact that KF files are case sensitive.

        """
        return self._access_rkf(lambda x: x.read_section(section), file)


    def get_rkf_skeleton(self, file='ams'):
        """Return a dictionary with the structure of a chosen ``.rkf`` file. Each key corresponds to a section name with the value being a set of variable names present in that section.

        The *file* argument should be the identifier of the file to read. It defaults to ``'ams'``. To access a file called ``something.rkf`` you need to call this function with ``file='something'``. If there exists only one engine results ``.rkf`` file, you can call this function with ``file='engine'`` to access this file.
        """
        return self._access_rkf(lambda x: x.get_skeleton(), file)


    def get_molecule(self, section, file='ams'):
        """Return a |Molecule| instance stored in a given *section* of a chosen ``.rkf`` file.

        The *file* argument should be the identifier of the file to read. It defaults to ``'ams'``. To access a file called ``something.rkf`` you need to call this function with ``file='something'``. If there exists only one engine results ``.rkf`` file, you can call this function with ``file='engine'`` to access this file.

        All data used by this method is taken from the chosen ``.rkf`` file. The ``molecule`` attribute of the corresponding job is ignored.
        """
        sectiondict = self.read_rkf_section(section, file)
        if sectiondict:
            return AMSResults._mol_from_rkf_section(sectiondict)


    def get_input_molecule(self):
        """Return a |Molecule| instance with the initial coordinates.

        All data used by this method is taken from ``ams.rkf`` file. The ``molecule`` attribute of the corresponding job is ignored.
        """
        return self.get_molecule('InputMolecule', 'ams')


    def get_main_molecule(self):
        """Return a |Molecule| instance with the final coordinates.

        All data used by this method is taken from ``ams.rkf`` file. The ``molecule`` attribute of the corresponding job is ignored.
        """
        return self.get_molecule('Molecule', 'ams')


    def get_history_molecule(self, step):
        """Return a |Molecule| instance with coordinates taken from a particular *step* in the ``History`` section of ``ams.rkf`` file.

        All data used by this method is taken from ``ams.rkf`` file. The ``molecule`` attribute of the corresponding job is ignored.
        """
        if 'ams' in self.rkfs:
            main = self.rkfs['ams']
            if 'History' not in main:
                raise KeyError("'History' section not present in {}".format(main.path))
            n = main.read('History', 'nEntries')
            if step > n:
                raise KeyError("Step {} not present in 'History' section of {}".format(step, main.path))
            coords = main.read('History', f'Coords({step})')
            coords = [coords[i:i+3] for i in range(0,len(coords),3)]
            if ('History', f'SystemVersion({step})') in main:
                version = main.read('History', f'SystemVersion({step})')
                if 'SystemVersionHistory' in main:
                    if ('SystemVersionHistory', 'blockSize') in main:
                        blockSize = main.read('SystemVersionHistory', 'blockSize')
                    else:
                        blockSize = 1
                    block = (version - 1) // blockSize + 1
                    offset = (version - 1) % blockSize
                    system = main.read('SystemVersionHistory', f'SectionNum({block})', return_as_list=True)[offset]
                else:
                    system = version
                mol = self.get_molecule(f'ChemicalSystem({system})')
                molsrc = f'ChemicalSystem({system})'
            else:
                mol = self.get_main_molecule()
                molsrc = 'Molecule'
            if len(mol) != len(coords):
                raise ResultsError(f'Coordinates taken from "History%Coords({step})" have incompatible length with molecule from {molsrc} section')
            for at, c in zip(mol, coords):
                at.move_to(c, unit='bohr')

            if ('History', f'LatticeVectors('+str(step)+')') in main:
                lattice = Units.convert(main.read('History', f'LatticeVectors('+str(step)+')'), 'bohr', 'angstrom')
                mol.lattice = [tuple(lattice[j:j+3]) for j in range(0,len(lattice),3)]

            if all(('History', i) in main for i in [f'Bonds.Index({step})', f'Bonds.Atoms({step})', f'Bonds.Orders({step})']):
                index = main.read('History', f'Bonds.Index({step})')
                atoms = main.read('History', f'Bonds.Atoms({step})')
                orders = main.read('History', f'Bonds.Orders({step})')
                for i in range(len(index)-1):
                    for j in range(index[i], index[i+1]):
                        mol.add_bond(mol[i+1], mol[atoms[j-1]], orders[j-1])
            return mol


    def get_engine_results(self, engine=None):
        """Return a dictionary with contents of ``AMSResults`` section from an engine results ``.rkf`` file.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        return self._process_engine_results(lambda x: x.read_section('AMSResults'), engine)


    def get_engine_properties(self, engine=None):
        """Return a dictionary with all the entries from ``Properties`` section from an engine results ``.rkf`` file.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        def properties(kf):
            n = kf.read('Properties', 'nEntries')
            ret = {}
            for i in range(1, n+1):
                tp = kf.read('Properties', 'Type({})'.format(i)).strip()
                stp = kf.read('Properties', 'Subtype({})'.format(i)).strip()
                val = kf.read('Properties', 'Value({})'.format(i))
                key = stp if stp.endswith(tp) else ('{} {}'.format(stp, tp) if stp else tp)
                ret[key] = val
            return ret
        return self._process_engine_results(properties, engine)


    def get_energy(self, unit='au', engine=None):
        """Return final energy, expressed in *unit*.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        return self._process_engine_results(lambda x: x.read('AMSResults', 'Energy'), engine) * Units.conversion_ratio('au', unit)


    def get_gradients(self, energy_unit='au', dist_unit='au', engine=None):
        """Return the gradients of the final energy, expressed in *energy_unit* / *dist_unit*.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        return np.asarray(self._process_engine_results(lambda x: x.read('AMSResults', 'Gradients'), engine)).reshape(-1,3) * Units.conversion_ratio('au', energy_unit) / Units.conversion_ratio('au', dist_unit)


    def get_stresstensor(self, engine=None):
        """Return the final stress tensor, expressed in atomic units.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        return np.asarray(self._process_engine_results(lambda x: x.read('AMSResults', 'StressTensor'), engine)).reshape(len(self.get_input_molecule().lattice),-1)


    def get_hessian(self, engine=None):
        """Return the Hessian matrix, i.e. the second derivative of the total energy with respect to the nuclear coordinates, expressed in atomic units.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        return np.asarray(self._process_engine_results(lambda x: x.read('AMSResults', 'Hessian'), engine)).reshape(3*len(self.get_input_molecule()),-1)


    def get_elastictensor(self, engine=None):
        """Return the elastic tensor, expressed in atomic units.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        et_flat = np.asarray(self._process_engine_results(lambda x: x.read('AMSResults', 'ElasticTensor'), engine))
        num_latvec = len(self.get_input_molecule().lattice)
        if num_latvec == 1:
            return et_flat.reshape(1,1)
        elif num_latvec == 2:
            return et_flat.reshape(3,3)
        else:
            return et_flat.reshape(6,6)


    def get_frequencies(self, unit='cm^-1', engine=None):
        """Return a numpy array of vibrational frequencies, expressed in *unit*.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        freqs = np.array(self._process_engine_results(lambda x: x.read('Vibrations', 'Frequencies[cm-1]'), engine))
        return freqs * Units.conversion_ratio('cm^-1', unit)


    def get_charges(self, engine=None):
        """Return the atomic charges, expressed in atomic units.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        return np.asarray(self._process_engine_results(lambda x: x.read('AMSResults', 'Charges'), engine))


    def get_dipolemoment(self, engine=None):
        """Return the electric dipole moment, expressed in atomic units.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        return np.asarray(self._process_engine_results(lambda x: x.read('AMSResults', 'DipoleMoment'), engine))


    def get_dipolegradients(self, engine=None):
        """Return the nuclear gradients of the electric dipole moment, expressed in atomic units. This is a (3*numAtoms x 3) matrix.
        """
        return np.asarray(self._process_engine_results(lambda x: x.read('AMSResults', 'DipoleGradients'), engine)).reshape(-1,3)


    def recreate_molecule(self):
        """Recreate the input molecule for the corresponding job based on files present in the job folder. This method is used by |load_external|.

        If ``ams.rkf`` is present in the job folder, extract data from the ``InputMolecule`` section.
        """
        if 'ams' in self.rkfs:
            return self.get_input_molecule()
        return None


    def recreate_settings(self):
        """Recreate the input |Settings| instance for the corresponding job based on files present in the job folder. This method is used by |load_external|.

        If ``ams.rkf`` is present in the job folder, extract user input and parse it back to a |Settings| instance using ``scm.input_parser`` module. Remove the ``system`` branch from that instance.
        """
        if 'ams' in self.rkfs:
            user_input = self.readrkf('General', 'user input')
            try:
                from scm.input_parser import input_to_settings
                inp = input_to_settings(user_input)
            except:
                log('Failed to recreate input settings from {}'.format(self.rkfs['ams'].path, 5))
                return None
            s = Settings()
            s.input = inp
            del s.input[ig('ams')][ig('system')]
            s.soft_update(config.job)
            return s
        return None


    def ok(self):
        """Check if the execution of the associated :attr:`job` was successful or not.
        See :meth:`Job.ok<scm.plams.core.basejob.Job.ok>` for more information."""
        return self.job.ok()

    @property
    def name(self):
        """Retrun the :attr:`job.name` of the job associated with this results instance."""
        return self.job.name



    #=========================================================================


    def _access_rkf(self, func, file='ams'):
        """A skeleton method for accessing any of the ``.rkf`` files produced by AMS.

        The *file* argument should be the identifier of the file to read. It defaults to ``'ams'``. To access a file called ``something.rkf`` you need to call this function with ``file='something'``. If there exists only one engine results ``.rkf`` file, you can call this function with ``file='engine'`` to access this file.

        The *func* argument has to be a function to call on a chosen ``.rkf`` file. It should take one argument, an instance of |KFFile|.
        """
        #Try unique engine:
        if file == 'engine':
            names = self.engine_names()
            if len(names) == 1:
                return func(self.rkfs[names[0]])
            else:
                raise ValueError("You cannot use 'engine' as 'file' argument if the engine results file is not unique. Please use the real name of the file you wish to read")

        #Try:
        if file in self.rkfs:
            return func(self.rkfs[file])

        #Try harder:
        filename = file + '.rkf'
        self.refresh()
        if filename in self.files:
            self.rkfs[file] = KFFile(opj(self.job.path, filename))
            return func(self.rkfs[file])

        #Surrender:
        raise FileError('File {} not present in {}'.format(filename, self.job.path))


    def _process_engine_results(self, func, engine=None):
        """A generic method skeleton for processing any engine results ``.rkf`` file. *func* should be a function that takes one argument (an instance of |KFFile|) and returns arbitrary data.

        The *engine* argument should be the identifier of the file you wish to read. To access a file called ``something.rkf`` you need to call this function with ``engine='something'``. The *engine* argument can be omitted if there's only one engine results file in the job folder.
        """
        names = self.engine_names()
        if engine is not None:
            if engine in names:
                return func(self.rkfs[engine])
            else:
                raise FileError('File {}.rkf not present in {}'.format(engine, self.job.path))
        else:
            if len(names) == 1:
                return func(self.rkfs[names[0]])
            else:
                raise ValueError("You need to specify the 'engine' argument when there are multiple engine result files present in the job folder")


    @staticmethod
    def _mol_from_rkf_section(sectiondict):
        """Return a |Molecule| instance constructed from the contents of the whole ``.rkf`` file section, supplied as a dictionary returned by :meth:`KFFile.read_section<scm.plams.tools.kftools.KFFile.read_section>`."""

        ret = Molecule()
        coords = [sectiondict['Coords'][i:i+3] for i in range(0,len(sectiondict['Coords']),3)]
        symbols = sectiondict['AtomSymbols'].split()
        for at, crd, sym in zip(sectiondict['AtomicNumbers'], coords, symbols):
            newatom = Atom(atnum=at, coords=crd, unit='bohr')
            if sym.startswith('Gh.'):
                sym = sym[3:]
                newatom.properties.ghost = True
            if '.' in sym:
                sym, name = sym.split('.', 1)
                newatom.properties.name = name
            ret.add_atom(newatom)
        if sectiondict['Charge'] != 0:
            ret.properties.charge = sectiondict['Charge']
        if 'nLatticeVectors' in sectiondict:
            ret.lattice = Units.convert([tuple(sectiondict['LatticeVectors'][i:i+3]) for i in range(0,len(sectiondict['LatticeVectors']),3)], 'bohr', 'angstrom')
        if 'EngineAtomicInfo' in sectiondict:
            suffixes = sectiondict['EngineAtomicInfo'].splitlines()
            for at, suffix in zip(ret, suffixes):
                at.properties.suffix = suffix
        return ret


#===========================================================================
#===========================================================================
#===========================================================================


class AMSJob(SingleJob):
    """A class representing a single computation with AMS driver. The corresponding results type is |AMSResults|.
    """
    _result_type = AMSResults


    def get_input(self):
        """Generate the input file. This method is just a wrapper around :meth:`_serialize_input`.

        Each instance of |AMSJob| or |AMSResults| present as a value in ``settings.input`` branch is replaced with an absolute path to ``ams.rkf`` file of that job.

        If you need to use a path to some engine specific ``.rkf`` file rather than the main ``ams.rkf`` file, you can to it by supplying a tuple ``(x, name)`` where ``x`` is an instance of |AMSJob| or |AMSResults| and ``name`` is a string with the name of the ``.rkf`` file you want. For example, ``(myjob, 'dftb')`` will transform to the absolute path to ``dftb.rkf`` file in ``myjob``'s folder, if such a file is present.

        Instances of |KFFile| are replaced with absolute paths to corresponding files.
        """

        special = {
            AMSJob: lambda x: x.results.rkfpath(),
            AMSResults: lambda x: x.rkfpath(),
            KFFile: lambda x: x.path,
            tuple: lambda x: AMSJob.tuple2rkf(x)
        }
        return self._serialize_input(special)


    def get_runscript(self):
        """Generate the runscript. Returned string is of the form::

            unset AMS_SWITCH_LOGFILE_AND_STDOUT
            AMS_JOBNAME=jobname AMS_RESULTSDIR=. $ADFBIN/ams [-n nproc] <jobname.in [>jobname.out]

        ``-n`` flag is added if ``settings.runscript.nproc`` exists. ``[>jobname.out]`` is used based on ``settings.runscript.stdout_redirect``.
        """
        ret  = 'unset AMS_SWITCH_LOGFILE_AND_STDOUT\n'
        ret += 'AMS_JOBNAME="{}" AMS_RESULTSDIR=. $ADFBIN/ams'.format(self.name)
        if 'nproc' in self.settings.runscript:
            ret += ' -n {}'.format(self.settings.runscript.nproc)
        ret += ' <"{}"'.format(self._filename('inp'))
        if self.settings.runscript.stdout_redirect:
            ret += ' >"{}"'.formatself._filename('out')
        ret += '\n\n'
        return ret


    def check(self):
        """Check if ``termination status`` variable from ``General`` section of main KF file equals ``NORMAL TERMINATION``."""
        try:
            status = self.results.readrkf('General', 'termination status')
        except:
            return False
        if 'NORMAL TERMINATION' in status:
            if 'errors' in status:
                log('Job {} reported errors. Please check the the output'.format(self.name), 1)
                return False
            if 'warnings' in status:
                log('Job {} reported warnings. Please check the the output'.format(self.name), 1)
            return True
        return False


    def hash_input(self):
        """Calculate the hash of the input file.

        All instances of |AMSJob| or |AMSResults| present as values in ``settings.input`` branch are replaced with hashes of corresponding job's inputs. Instances of |KFFile| are replaced with absolute paths to corresponding files.
        """
        special = {
            AMSJob: lambda x: x.hash_input(),
            AMSResults: lambda x: x.job.hash_input(),
            KFFile: lambda x: x.path,
            tuple: lambda x: AMSJob.tuple2rkf(x)
        }
        return sha256(self._serialize_input(special))


    #=========================================================================


    def _serialize_input(self, special):
        """Transform the contents of ``settings.input`` branch into string with blocks, keys and values.

        First, the contents of ``settings.input`` are extended with entries returned by :meth:`_serialize_molecule`. Then the contents of ``settings.input.ams`` are used to generate AMS text input. Finally, every other (than ``ams``) entry in ``settings.input`` is used to generate engine specific input.

        Special values can be indicated with *special* argument, which should be a dictionary having types of objects as keys and functions translating these types to strings as values.
        """

        def unspec(value):
            """Check if *value* is one of a special types and convert it to string if it is."""
            for spec_type in special:
                if isinstance(value, spec_type):
                    return special[spec_type](value)
            return value

        def serialize(key, value, indent, end='end'):
            """Given a *key* and its corresponding *value* from the |Settings| instance produce a snippet of the input file representing this pair.

            If the value is a nested |Settings| instance, use recursive calls to build the snippet for the entire block. Indent the result with *indent* spaces.
            """
            ret = ''
            if isinstance(value, Settings):
                ret += ' '*indent + key
                if '_h' in value:
                    ret += ' ' + unspec(value['_h'])
                ret += '\n'

                i = 1
                while ('_'+str(i)) in value:
                    ret += serialize('', value['_'+str(i)], indent+2)
                    i += 1

                for el in value:
                    if not el.startswith('_'):
                        if key.lower().startswith('engine') and el.lower() == 'input':
                            ret += serialize(el, value[el], indent+2, 'endinput')
                        else:
                            ret += serialize(el, value[el], indent+2)

                ret += ' '*indent + end+'\n'

            elif isinstance(value, list):
                for el in value:
                    ret += serialize(key, el, indent)
            elif value == '' or value is True:
                ret += ' '*indent + key + '\n'
            elif value is False or value is None:
                pass
            else:
                ret += ' '*indent + key + ' ' + str(unspec(value)) + '\n'
            return ret

        fullinput = self.settings.input.copy()

        #prepare contents of 'system' block(s)
        more_systems = self._serialize_molecule()
        if more_systems:
            if ig('system') in fullinput[ig('ams')]:
                #nonempty system block was already present in input.ams
                system = fullinput[ig('ams')][ig('system')]
                system_list = system if isinstance(system, list) else [system]

                system_list_set = Settings({(s._h if '_h' in s else ''):s   for s in system_list})
                more_systems_set = Settings({(s._h if '_h' in s else ''):s   for s in more_systems})

                system_list_set += more_systems_set
                system_list = list(system_list_set.values())
                system = system_list[0] if len(system_list) == 1 else system_list
                fullinput[ig('ams')][ig('system')] = system

            else:
                fullinput[ig('ams')][ig('system')] = more_systems[0] if len(more_systems) == 1 else more_systems

        txtinp = ''
        ams = fullinput.find_case('ams')

        #contents of the 'ams' block (AMS input) go first
        for item in fullinput[ams]:
            txtinp += serialize(item, fullinput[ams][item], 0) + '\n'

        #and then engines
        for engine in fullinput:
            if engine != ams:
                txtinp += serialize('engine '+engine, fullinput[engine], 0, end='endengine') + '\n'

        return txtinp


    def _serialize_molecule(self):
        """Return a list of |Settings| instances containing the information about one or more |Molecule| instances stored in the ``molecule`` attribute.

        Molecular charge is taken from ``molecule.properties.charge``, if present. Additional, atom-specific information to be put in ``atoms`` block after XYZ coordinates can be supplied with ``atom.properties.suffix``.

        If the ``molecule`` attribute is a dictionary, the returned list is of the same length as the size of the dictionary. Keys from the dictionary are used as headers of returned ``system`` blocks.
        """

        if self.molecule is None:
            return Settings()

        moldict = {}
        if isinstance(self.molecule, Molecule):
            moldict = {'':self.molecule}
        elif isinstance(self.molecule, dict):
            moldict = self.molecule
        else:
            raise JobError("Incorrect 'molecule' attribute of job {}. 'molecule' should be a Molecule, a dictionary or None".format(self.name))

        ret = []
        for name, molecule in moldict.items():
            newsystem = Settings()
            if name:
                newsystem._h = name

            if len(molecule.lattice) in [1,2] and molecule.align_lattice():
                log("The lattice of {} Molecule supplied for job {} did not follow the convention required by AMS. I rotated the whole system for you. You're welcome".format(name if name else 'main', self.name), 3)

            newsystem.atoms._1 = [atom.str(symbol=self._atom_symbol(atom), space=18, decimal=10,
                    suffix=(atom.properties.suffix if 'suffix' in atom.properties else '')) for atom in molecule]

            if molecule.lattice:
                newsystem.lattice._1 = ['{:16.10f} {:16.10f} {:16.10f}'.format(*vec) for vec in molecule.lattice]

            if ig('charge') in molecule.properties:
                newsystem.charge = molecule.properties[ig('charge')]

            ret.append(newsystem)

        return ret


    #=========================================================================


    @staticmethod
    def _atom_symbol(atom):
        """Return the atomic symbol of *atom*. Ensure proper formatting for AMSuite input taking into account ``ghost`` and ``name`` entries in ``properties`` of *atom*."""
        smb = atom.symbol if atom.atnum > 0 else ''  #Dummy atom should have '' instead of 'Xx'
        if 'ghost' in atom.properties and atom.properties.ghost:
            smb = ('Gh.'+smb).rstrip('.')
        if 'name' in atom.properties:
            smb = (smb+'.'+str(atom.properties.name)).lstrip('.')
        return smb


    @staticmethod
    def _tuple2rkf(arg):
        """Transform a pair ``(x, name)`` where ``x`` is an instance of |AMSJob| or |AMSResults| and ``name`` is a name of ``.rkf`` file (``ams`` or engine) to an absolute path to that ``.rkf`` file."""
        if len(arg) == 2 and isinstance(arg[1], str):
            if isinstance(arg[0], AMSJob):
                return arg[0].results.rkfpath(arg[1])
            if isinstance(arg[0], AMSResults):
                return arg[0].rkfpath(arg[1])
        return str(arg)


    @classmethod
    def from_inputfile(cls, filename: str, heredoc_delimit: str = 'eor', **kwargs) -> 'AMSJob':
        """Construct an :class:`AMSJob` instance from an ADF inputfile or runfile.

        If a runscript is provide than this method will attempt to extract the input file based
        on the heredoc delimiter (see *heredoc_delimit*).

        """
        try:
            from scm.input_parser.parse import input_to_settings
        except ImportError:  # Try to load the parser from $ADFHOME/scripting
            with UpdateSysPath():
                from scm.input_parser.parse import input_to_settings

        s = Settings()
        with open(filename, 'r') as f:
            inp_file = parse_heredoc(f.read(), heredoc_delimit)

        s.input = input_to_settings(inp_file, cls._command)
        if not s.input:
            raise JobError(f"from_inputfile: failed to parse '{filename}'")

        # Extract a molecule from the input settings
        mol = cls.settings_to_mol(s)

        # Create and return the Job instance
        if mol is not None:
            return cls(molecule=mol, settings=s, **kwargs)
        else:
            s.ignore_molecule = True
            return cls(settings=s, **kwargs)


    @staticmethod
    def settings_to_mol(s: Settings) -> dict:
        """Pop the `s.input.ams.system` block from a settings instance and convert it into a dictionary of molecules.

        The provided settings should be in the same style as the ones produced by the SCM input parser.
        Dictionary keys are taken from the header of each system block.
        The existing `s.input.ams.system` block is removed in the process, assuming it was present in the first place.

        """
        def read_mol(settings_block: Settings) -> Molecule:
            """Retrieve single molecule from a single `s.input.ams.system` block."""
            mol = Molecule()
            for atom in settings_block.atoms._1:
                # Extract arguments for Atom()
                symbol, x, y, z, *comment = atom.split(maxsplit=4)
                kwargs = {} if not comment else {'suffix': comment[0]}
                coords = float(x), float(y), float(z)

                try:
                    at = Atom(symbol=symbol, coords=coords, **kwargs)
                except PTError:  # It's either a ghost atom and/or an atom with a custom name
                    if symbol.startswith('Gh.'):  # Ghost atom
                        kwargs['ghost'], symbol = symbol.split('.', maxsplit=1)
                    if '.' in symbol:  # Atom with a custom name
                        symbol, kwargs['name'] = symbol.split('.', maxsplit=1)
                    at = Atom(symbol=symbol, coords=coords, **kwargs)

                mol.add_atom(at)

            # Add bonds
            for bond in settings_block.bondorders._1:
                _at1, _at2, _order = bond.split()
                at1, at2, order = mol[int(_at1)], mol[int(_at2)], float(_order)
                mol.add_bond(at1, at2, order)

            # Set the lattice vector if applicable
            if settings_block.lattice._1:
                mol.lattice = [tuple(float(j) for j in i.split()) for i in settings_block.lattice._1]

            # Set the molecular charge
            if settings_block.charge:
                mol.properties.charge = float(settings_block.charge)

            mol.properties.name = str(settings_block._h)
            return mol

        # Raises a KeyError if the `system` key is absent
        with s.supress_missing():
            try:
                settings_list = s.input.ams.pop('system')
            except KeyError:  # The block s.input.ams.system is absent
                return None

        # Create a new dictionary with system headers as keys and molecules as values
        moldict = {}
        for settings_block in settings_list:
            key = str(settings_block._h)
            if key in moldict:
                raise KeyError(f"Duplicate system headers found in s.input.ams.system: {repr(key)}")
            moldict[key] = read_mol(settings_block)

        return moldict

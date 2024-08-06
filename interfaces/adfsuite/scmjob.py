import os
from os.path import join as opj

import numpy as np

from scm.plams.core.basejob import SingleJob
from scm.plams.core.errors import FileError, JobError, PlamsError, ResultsError
from scm.plams.core.functions import log, parse_heredoc
from scm.plams.core.private import sha256
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.mol.atom import Atom
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.kftools import KFFile
from scm.plams.tools.units import Units


class SCMResults(Results):
    """Abstract class gathering common mechanisms for results of ADF Suite programs."""

    _kfext = ""

    def collect(self):
        """Collect files present in the job folder. Use parent method from |Results|, then create an instance of |KFFile| for the main KF file and store it as ``_kf`` attribute."""
        Results.collect(self)
        kfname = self.job.name + self.__class__._kfext
        if kfname in self.files:
            self._kf = KFFile(opj(self.job.path, kfname))
        else:
            log("WARNING: Main KF file {} not present in {}".format(kfname, self.job.path), 1)

    def refresh(self):
        """Refresh the contents of ``files`` list. Use parent method from |Results|, then look at all attributes that are instances of |KFFile| and check if they point to existing files. If not, try to reinstantiate them with current job path (that can happen while loading a pickled job after the entire job folder was moved)."""
        Results.refresh(self)
        to_remove = []
        for attr, val in self.__dict__.items():
            if isinstance(val, KFFile) and os.path.dirname(val.path) != self.job.path:
                guessnewpath = opj(self.job.path, os.path.basename(val.path))
                if os.path.isfile(guessnewpath):
                    self.__dict__[attr] = KFFile(guessnewpath)
                else:
                    to_remove.append(attr)
        for i in to_remove:
            del self.__dict__[i]

    def readkf(self, section, variable):
        """readkf(section, variable)
        Read data from *section*/*variable* of the main KF file.

        The type of the returned value depends on the type of *variable* defined inside KF file. It can be: single int, list of ints, single float, list of floats, single boolean, list of booleans or string.
        """
        if self._kfpresent():
            return self._kf.read(section, variable)
        raise FileError("File {} not present in {}".format(self.job.name + self.__class__._kfext, self.job.path))

    def newkf(self, filename):
        """newkf(filename)
        Create new |KFFile| instance using file *filename* in the job folder.

        Example usage::

            >>> res = someadfjob.run()
            >>> tape13 = res.newkf('$JN.t13')
            >>> print(tape13.read('Geometry', 'xyz'))

        """
        self.refresh()
        filename = filename.replace("$JN", self.job.name)
        if filename in self.files:
            return KFFile(opj(self.job.path, filename))
        else:
            raise FileError("File {} not present in {}".format(filename, self.job.path))

    def get_properties(self):
        """get_properties()
        Return a dictionary with all the entries from ``Properties`` section in the main KF file.
        """
        n = self.readkf("Properties", "nEntries")
        ret = {}
        for i in range(1, n + 1):
            tp = self.readkf("Properties", "Type({})".format(i)).strip()
            stp = self.readkf("Properties", "Subtype({})".format(i)).strip()
            val = self.readkf("Properties", "Value({})".format(i))
            key = stp if stp.endswith(tp) else ("{} {}".format(stp, tp) if stp else tp)
            ret[key] = val
        return ret

    def get_molecule(self, section, variable, unit="bohr", internal=False, n=1):
        """get_molecule(section, variable, unit='bohr', internal=False, n=1)
        Read molecule coordinates from *section*/*variable* of the main KF file.

        Returned |Molecule| instance is created by copying a molecule from associated |SCMJob| instance and updating atomic coordinates with values read from *section*/*variable*. The format in which coordinates are stored is not consistent for all programs or even for different sections of the same KF file. Sometimes coordinates are stored in bohr, sometimes in angstrom. The order of atoms can be either input order or internal order. These settings can be adjusted with *unit* and *internal* parameters. Some variables store more than one geometry, in those cases *n* can be used to choose the preferred one.
        """
        atnums = self._atomic_numbers_input_order()
        natoms = len(atnums)
        coords = self.readkf(section, variable)
        coords = [coords[i : i + 3] for i in range(0, len(coords), 3)]
        if len(coords) > natoms:
            if len(coords) < n * natoms:
                raise ResultsError(
                    "get_molecule() failed. Not enough data in {}%{} to extract geometry no {}".format(
                        section, variable, n
                    )
                )
            coords = coords[(n - 1) * natoms : n * natoms]
        if internal:
            mapping = self._int2inp()
            coords = [coords[mapping[i] - 1] for i in range(len(coords))]
        ret = Molecule()
        for z, crd in zip(atnums, coords):
            ret.add_atom(Atom(atnum=z, coords=crd, unit=unit))
        return ret

    def _get_single_value(self, section, variable, output_unit, native_unit="au"):
        """_get_single_value(section, variable, output_unit, native_unit='au')

        A small method template for all the single number "get_something()" methods extracting data from main KF file. Returned value is converted from *native_unit* to *output_unit*.
        """
        if self._kfpresent() and (section, variable) in self._kf:
            return Units.convert(self.readkf(section, variable), native_unit, output_unit)
        raise ResultsError("'{}%{}' not present in {}".format(section, variable, self._kfpath()))

    def _atomic_numbers_input_order(self):
        """_atomic_numbers_input_order()
        Return a list of atomic numbers, in the input order. Abstract method.
        """
        raise PlamsError("Trying to run an abstract method SCMResults._atomic_numbers_input_order()")

    def kfpath(self):
        """kfpath()
        Return the absolute path to the main KF file.
        """
        return self._kfpath()

    def _kfpath(self):
        """_kfpath()
        Return the absolute path to the main KF file.
        """
        return opj(self.job.path, self.job.name + self.__class__._kfext)

    def _kfpresent(self):
        """_kfpresent()
        Check if this instance has a valid ``_kf`` attribute.
        """
        return hasattr(self, "_kf") and isinstance(self._kf, KFFile)

    def _export_attribute(self, attr, other):
        """_export_attribute(attr, other)
        If *attr* is a KF file take care of a proper path. Otherwise use parent method. See :meth:`Results._copy_to<scm.plams.core.results.Results._copy_to>` for details.
        """
        if isinstance(attr, KFFile):
            oldname = os.path.basename(attr.path)
            newname = Results._replace_job_name(oldname, self.job.name, other.job.name)
            newpath = opj(other.job.path, newname)
            return KFFile(newpath) if os.path.isfile(newpath) else None
        else:
            return Results._export_attribute(self, attr, other)

    def _int2inp(self):
        """_int2inp()
        Obtain mapping from internal atom order to the input one. Abstract method.
        """
        raise PlamsError("Trying to run an abstract method SCMResults._int2inp()")

    def to_input_order(self, data):
        """to_input_order(self, data)
        Reorder any iterable *data* from the internal atom order to the input atom order. The length of *data* must be equal to the number of atoms, otherwise an exception is raised. Returned value is a container of the same type as *data*.
        """
        mapping = self._int2inp()
        if len(mapping) != len(data):
            raise PlamsError(
                "to_input_order() got an argument with incorrect length. Length must be equal to the number of atoms"
            )
        t = np.array if type(data) == np.ndarray else type(data)
        return t([data[mapping[i] - 1] for i in range(len(mapping))])

    def readarray(self, section: str, subsection: str, **kwargs) -> "np.ndarray":
        """Read data from *section*/*subsection* of the main KF file and return as NumPy array.

        All additional provided keyword arguments will be passed onto the numpy.array_ function.

        .. _numpy.array: https://docs.scipy.org/doc/numpy/reference/generated/numpy.array.html
        """
        return np.array(self.readkf(section, subsection), **kwargs)


class SCMJob(SingleJob):
    """Abstract class gathering common mechanisms for jobs with ADF Suite programs."""

    _result_type = SCMResults
    _top = ["title", "units", "define"]
    _command = ""
    _json_definitions = _command
    _subblock_end = "subend"
    _legacy = ["band", "dftb", "uff"]

    def __init__(self, **kwargs):
        SingleJob.__init__(self, **kwargs)
        if self.__class__._command in self.__class__._legacy:
            log(
                "LEGACY WARNING: Job {} uses executable '{}' which is not present in AMS2018. Please use AMSJob (unless you're running an older version of ADFSuite)".format(
                    self.name, self.__class__._command
                ),
                1,
            )

    def get_input(self):
        """Generate the input file. This method is just a wrapper around :meth:`_serialize_input`.

        Each instance of |SCMJob| or |SCMResults| present as a value in ``settings.input`` branch is replaced with an absolute path to the main KF file of that job.
        """
        special = {SCMJob: lambda x: x.results._kfpath(), SCMResults: lambda x: x._kfpath(), KFFile: lambda x: x.path}
        return self._serialize_input(special)

    def get_runscript(self):
        """Generate a runscript. Returned string is of the form::

            $AMSBIN/name [-n nproc] <jobname.in [>jobname.out]

        ``name`` is taken from the class attribute ``_command``. ``-n`` flag is added if ``settings.runscript.nproc`` exists. ``[>jobname.out]`` is used based on ``settings.runscript.stdout_redirect``.
        """
        s = self.settings.runscript
        ret = "$AMSBIN/" + self._command
        if "nproc" in s:
            ret += " -n " + str(s.nproc)
        ret += ' <"{}"'.format(self._filename("inp"))
        if s.stdout_redirect:
            ret += ' >"{}"'.format(self._filename("out"))
        ret += "\n\n"
        return ret

    def check(self):
        """Check if ``termination status`` variable from ``General`` section of main KF file equals ``NORMAL TERMINATION``."""
        try:
            status = self.results.readkf("General", "termination status")
        except:
            return False
        if "NORMAL TERMINATION" in status:
            if "errors" in status:
                log("Job {} reported errors. Please check the output".format(self.name), 1)
                return False
            if "warnings" in status:
                log("Job {} reported warnings. Please check the output".format(self.name), 1)
            return True
        return False

    def hash_input(self):
        """Calculate the hash of the input file.

        All instances of |SCMJob| or |SCMResults| present as values in ``settings.input`` branch are replaced with hashes of corresponding job's inputs.
        """
        special = {SCMJob: lambda x: x.hash_input(), SCMResults: lambda x: x.job.hash_input(), KFFile: lambda x: x.path}
        return sha256(self._serialize_input(special))

    def _serialize_input(self, special):
        """Transform all contents of ``setting.input`` branch into string with blocks, keys and values.

        On the highest level alphabetic order of iteration is modified: keys occuring in attribute ``_top`` are printed first. Special values can be indicated with *special* argument, which should be a dictionary having types of objects as keys and functions translating these types to strings as values.

        Automatic handling of ``molecule`` can be disabled with ``settings.ignore_molecule = True``.
        """

        def unspec(value):
            """Check if *value* is one of a special types and convert it to string if it is."""
            for spec_type in special:
                if isinstance(value, spec_type):
                    return special[spec_type](value)
            return value

        def serialize(key, value, indent):
            """Given a *key* and its corresponding *value* from the |Settings| instance produce a snippet of the input file representing this pair.

            If the value is a nested |Settings| instance, use recursive calls to build the snippet for the entire block. Indent the result with *indent* spaces.
            """
            ret = ""
            if isinstance(value, Settings):
                ret += " " * indent + key
                if "_h" in value:
                    ret += " " + unspec(value["_h"])
                ret += "\n"

                i = 1
                while ("_" + str(i)) in value:
                    ret += serialize("", value["_" + str(i)], indent + 2)
                    i += 1

                for el in value:
                    if not el.startswith("_"):
                        ret += serialize(el, value[el], indent + 2)

                if indent == 0:
                    ret += "end\n"
                else:
                    ret += " " * indent + self._subblock_end + "\n"
            elif isinstance(value, list):
                for el in value:
                    ret += serialize(key, el, indent)
            elif value == "" or value is True:
                ret += " " * indent + key + "\n"
            elif value is False:
                pass
            else:
                value = str(unspec(value))
                ret += " " * indent + key
                if key != "" and not value.startswith("="):
                    ret += " "
                ret += value + "\n"
            return ret

        use_molecule = ("ignore_molecule" not in self.settings) or (self.settings.ignore_molecule == False)
        if use_molecule:
            self._serialize_mol()

        inp = ""
        for item in self._top:
            item = self.settings.input.find_case(item)
            if item in self.settings.input:
                inp += serialize(item, self.settings.input[item], 0) + "\n"
        for item in self.settings.input:
            if item.lower() not in self._top:
                inp += serialize(item, self.settings.input[item], 0) + "\n"

        if use_molecule:
            self._remove_mol()
        return inp

    def _serialize_mol(self):
        """Process |Molecule| instance stored in ``molecule`` attribute and add it as relevant entries of ``settings.input`` branch. Abstract method."""
        raise PlamsError("Trying to run an abstract method SCMJob._serialize_mol()")

    def _remove_mol(self):
        """Remove from ``settings.input`` all entries added by :meth:`_serialize_mol`. Abstract method."""
        raise PlamsError("Trying to run an abstract method SCMJob._remove_mol()")

    @staticmethod
    def _atom_symbol(atom):
        """Return the atomic symbol of *atom*. Ensure proper formatting for ADFSuite input taking into account ``ghost`` and ``name`` entries in ``properties`` of *atom*."""
        smb = atom.symbol
        if "ghost" in atom.properties and atom.properties.ghost:
            smb = ("Gh." + smb).rstrip(".")
        if "name" in atom.properties:
            smb = (smb + "." + str(atom.properties.name)).lstrip(".")
        return smb

    @classmethod
    def from_inputfile(cls, filename: str, heredoc_delimit: str = "eor", **kwargs) -> "SCMJob":
        """Construct a :class:`SCMJob` instance from an ADF inputfile.

        If a runscript is provided then this method will attempt to extract the input file based
        on the heredoc delimiter (see *heredoc_delimit*).

        """
        from scm.plams.interfaces.adfsuite.inputparser import InputParserFacade

        s = Settings()
        with open(filename, "r") as f:
            inp_file = parse_heredoc(f.read(), heredoc_delimit)

        s.input = InputParserFacade().to_settings(cls._json_definitions or cls._command, inp_file)
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
    def settings_to_mol(s: Settings) -> None:
        """An abstract method for extracting molecules from input settings (see :meth:`SCMJob.from_inputfile`)."""
        return None

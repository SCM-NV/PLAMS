"""
Run DFTB+ with plams
contributed by Patrick Melix
based on code from Michal Handzlik

see documentation for an example
"""

from os.path import join as opj

from scm.plams.core.basejob import SingleJob
from scm.plams.core.errors import MoleculeError
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.units import Units

__all__ = ["DFTBPlusJob", "DFTBPlusResults"]


class DFTBPlusResults(Results):
    """A Class for handling DFTB+ Results."""

    _outfile = "detailed.out"
    _xyzout = "geo_end.xyz"
    _genout = "geo_end.gen"

    def get_molecule(self):
        """get_molecule()
        Return the molecule from the 'geo_end.gen' file. If there is an Error, try to read the 'geo_end.xyz' file.
        """
        try:
            # The .gen file contains the cell, ASE can read it
            mol = Molecule(self[self._genout], inputformat="ase")
        except MoleculeError:
            # Fallback if no ASE found
            mol = Molecule(filename=self[self._xyzout])
        except:
            mol = Molecule()
        return mol

    def get_energy(self, string="Total energy", unit="au"):
        """get_energy(string='Total energy', unit='au')
        Return the energy given in the output with the description *string*, expressed in *unit*. Defaults to ``Total energy`` and ``au``.
        """
        try:
            energy = float(self.grep_file(self._outfile, pattern=string + ":")[0].split()[2])
            energy = Units.convert(energy, "au", unit)
        except:
            energy = float("nan")
        return energy

    def get_atomic_charges(self):
        """get_atomic_charges()
        Returns a dictonary with atom numbers and their charges, ordering is the same as in the input.
        """
        try:
            atomic_charges = {}
            string = self.awk_file(
                self._outfile, script="/Net atomic charges/{do_print=1} NF==0 {do_print=0 }do_print==1 {print}"
            )
            for line in string:
                if string[0] == line or string[1] == line:
                    continue
                l = line.split()
                atomic_charges[l[0]] = float(l[1])
        except:
            atomic_charges = {}
        return atomic_charges


class DFTBPlusJob(SingleJob):
    """A class representing a single computational job with DFTB+.
    Only supports molecular coordinates, no support for lattice yet."""

    _result_type = DFTBPlusResults
    _filenames = {"inp": "dftb_in.hsd", "run": "$JN.run", "out": "$JN.out", "err": "$JN.err", "gen": "$JN.gen"}

    def get_input(self):
        """Transform all contents of ``setting.input`` branch into string with blocks, keys and values.

        Automatic handling of ``molecule`` can be disabled with ``settings.ignore_molecule = True``.
        """

        def parse(key, value, indent=""):
            if value is True:
                value = ""

            ret = indent + key
            if key != "" and value != "":
                ret += " ="

            if isinstance(value, Settings):
                if "_h" in value:
                    ret += " " + value["_h"]
                ret += " {\n"

                i = 1
                while ("_" + str(i)) in value:
                    ret += parse("", value["_" + str(i)], indent + "  ")
                    i += 1

                for el in value:
                    if not el.startswith("_"):
                        ret += parse(el, value[el], indent + "  ")
                ret += indent + "}"
            else:
                ret += " " + str(value)
            ret += "\n"
            return ret

        inp = ""
        use_molecule = ("ignore_molecule" not in self.settings) or (self.settings.ignore_molecule == False)
        if use_molecule:
            self._parsemol()

        for item in self.settings.input:
            inp += parse(item, self.settings.input[item])

        if use_molecule:
            self._removemol()
        return inp

    def _parsemol(self):
        # use ASE to write molecule if available
        if "ase" in Molecule._writeformat:
            filename = opj(self.path, self._filename("gen"))
            self.molecule.write(filename, outputformat="ase", format="gen")
            self.settings.input.geometry._h = "GenFormat"
            self.settings.input.geometry._1 = "<<< " + self._filename("gen")

        else:
            # Old way of handling gen-format ourselves, delete if ASE becomes obligatory
            atom_types = {}
            n = 1
            atoms_line = ""
            for atom in self.molecule:
                if atom.symbol not in atom_types:
                    atoms_line += atom.symbol + " "
                    atom_types[atom.symbol] = n
                    n += 1

            # check PBC
            lattice = []
            geomType = "C"
            for vec in self.molecule.lattice:
                if not all(isinstance(x, (int, float)) for x in vec):
                    raise ValueError("Non-Number in Lattice Vectors, not compatible with DFTBPlus")

                lattice.append(vec)
                geomType = "S"

            self.settings.input.geometry._h = "GenFormat"
            self.settings.input.geometry._1 = "%i %s" % (len(self.molecule), geomType)
            self.settings.input.geometry._2 = atoms_line
            self.settings.input.geometry._3 = ""

            for i, atom in enumerate(self.molecule):
                self.settings.input.geometry["_" + str(i + 4)] = ("%5i" % (i + 1)) + atom.str(
                    symbol=str(atom_types[atom.symbol])
                )

            if len(vec) > 0:
                j = i + 1
                # origin
                self.settings.input.geometry["_" + str(j + 4)] = "0.0 0.0 0.0"
                j += 1
                for i, vec in enumerate(lattice):
                    self.settings.input.geometry["_" + str(i + j + 4)] = "%f %f %f" % (vec)

    def _removemol(self):
        if "geometry" in self.settings.input:
            del self.settings.input.geometry

    def get_runscript(self):
        """dftb+ has to be in your $PATH!"""
        ret = "dftb+ "
        if self.settings.runscript.stdout_redirect:
            ret += " >" + self._filename("out")
        ret += "\n\n"
        return ret

    def check(self):
        """Returns true if 'ERROR!' is not found in the output."""
        s = self.results.grep_output("ERROR!")
        return len(s) == 0

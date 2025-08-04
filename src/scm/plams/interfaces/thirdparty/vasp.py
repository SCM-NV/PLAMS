from os.path import join as opj

from scm.plams.core.basejob import SingleJob
from scm.plams.core.errors import PlamsError
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.units import Units

__all__ = ["VASPJob", "VASPResults"]


class VASPResults(Results):
    """A class for VASP results."""

    def get_energy(self, index=-1, unit="a.u."):
        """Returns sigma->0 (!!!) energy without entropy."""
        s = self.grep_output("energy  without entropy=")[index]
        if not isinstance(index, slice):
            return Units.convert(float(s.split()[-1]), "eV", unit)
        else:
            return [Units.convert(float(x.split()[-1]), "eV", unit) for x in s]

    def get_dispersion_energy(self, index=-1, unit="a.u."):
        """Returns `Edisp (eV)` from the OUTCAR."""
        s = self.grep_output("Edisp (eV)")[index]
        if not isinstance(index, slice):
            return Units.convert(float(s.split()[-1]), "eV", unit)
        else:
            return [Units.convert(float(x.split()[-1]), "eV", unit) for x in s]


class VASPJob(SingleJob):
    """
    A class representing a single computational job with `VASP <https://www.vasp.at/>`

    * Set 'ignore_molecule' in ``self.settings`` to disable Molecule handling through ASE.
    * Set 'ignore_potcar' in ``self.settings`` to disable automatic `POTCAR` creation.
    * Set the path to the `POTCAR` files in  ``self.settings.input.potcar`` for automatic `POTCAR` creation.
    * If `POTCAR` files not matching the element symbol should be used, give a translation dict in ``self.settings.input.potcardict``.
        E.g. `{'Fe': 'Fe_pv'}`.
    * Settings branch ``input.incar`` is parsed into the `INCAR` file, ``input.xxx`` into the corresponding `XXX` file.
    * Use the PLAMS notation `_h`, `_1`, `_2`, ... to obtain keywords in specific order (e.g. for the KPOINTS file).
    """

    _command = "vasp_std"
    _filenames = {"inp": "INCAR", "run": "$JN.run", "out": "OUTCAR", "err": "$JN.err", "log": "$JN.log"}
    _result_type = VASPResults

    def get_input(self):
        """
        Transform all contents of ``input`` branch of ``settings`` into string.
        """

        def vaspstr(inp):
            # convert to VASP type string
            if isinstance(inp, Settings):
                raise PlamsError("Nested Settings Object not supported in VASP parser.")
            if isinstance(inp, bool):
                if inp:
                    return ".TRUE."
                else:
                    return ".FALSE."
            elif isinstance(inp, list):
                return " ".join([str(x) for x in inp])
            else:
                return str(inp).upper()

        def parse(key, value):
            ret = ""

            if isinstance(value, Settings):
                if "_h" in value:
                    ret += "{}\n".format(vaspstr(value["_h"]))
                i = 1
                while ("_" + str(i)) in value:
                    ret += "{}\n".format(vaspstr(value["_" + str(i)]))
                    i += 1
                for el in value:
                    if not el.startswith("_"):
                        ret += parse(el, value[el])

            elif isinstance(value, list):
                ret += "{}\n".format(vaspstr(value))
            else:
                ret += "{} = {}\n".format(key.upper(), vaspstr(value).upper())
            return ret

        use_molecule = ("ignore_molecule" not in self.settings) or (self.settings.ignore_molecule == False)
        use_potcar = ("ignore_potcar" not in self.settings) or (self.settings.ignore_potcar == False)
        if use_molecule:
            self._parsemol()
        if use_potcar:
            self._parsepotcar()

        # we need all first-level keys in uppercase
        tmp = Settings()
        for key in self.settings.input:
            tmp[key.upper()] = self.settings.input[key]

        for item in tmp:
            # POTCAR creation handled above
            if "POTCAR" in item:
                continue
            # INCAR
            elif item == "INCAR":
                inp = parse(item, tmp[item])
            else:
                with open(opj(self.path, item), "w") as f:
                    f.write(parse(item, tmp[item]))
        return inp

    def _parsemol(self):
        if "ase" in Molecule._writeformat:
            # ASE has a write function for VASP coordinate files, use that if possible
            filename = opj(self.path, "POSCAR")
            self.molecule.writease(filename)
        else:
            raise PlamsError("VASP Interface has no builtin Molecule support, install ASE. See Doc for details.")

    def _parsepotcar(self):
        tree = self.settings.input
        if "potcar" in tree:
            elements = [self.molecule.atoms[0].symbol]
            for atom in self.molecule.atoms[1:]:
                if not atom.symbol == elements[-1]:
                    elements.append(atom.symbol)
            if "potcardict" in tree:
                translate = dict(tree.potcardict)
            else:
                translate = {el: el for el in set(elements)}
            # open files
            files = [open(opj(tree.potcar, translate[el], "POTCAR"), "r") for el in elements]

            with open(opj(self.path, "POTCAR"), "w") as f:
                for potcar in files:
                    f.write(potcar.read())
                    potcar.close()
        else:
            raise PlamsError("VASP Interface needs the POTCAR path in self.settings.input.potcar.")

    def get_runscript(self):
        """
        Run VASP

        Overwrite ``self._command`` to change the default VASP Binary.
        """
        ret = self._command
        if self.settings.runscript.stdout_redirect:
            ret += " >" + self._filename("log")
        ret += "\n\n"
        return ret

    def check(self):
        """
        Look for the normal termination line in output. Note, that does not mean your calculation was successful!
        """
        termination = self.results.grep_output("General timing and")
        return bool(termination)

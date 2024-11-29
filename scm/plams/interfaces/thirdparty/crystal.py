from os.path import join as opj

import numpy as np

from scm.plams.core.basejob import SingleJob
from scm.plams.core.errors import PlamsError
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Molecule

__all__ = ["CrystalJob", "mol2CrystalConf"]


class CrystalJob(SingleJob):
    """
    A class representing a single computational job with `CRYSTAL <https://www.crystal.unito.it/>`
    Use Crystal version >= 14, lower versions have an even stricter and more complicated input.
    """

    _command = "crystal"
    _filenames = {"inp": "INPUT", "run": "$JN.run", "out": "$JN.out", "err": "$JN.err"}

    def get_input(self):
        """
        Transform all contents of ``input`` branch of ``settings`` into string.
        """

        # these are the geometry keys from the manual, first mandatory block
        _geom_keys = ["CRYSTAL", "SLAB", "POLYMER", "HELIX", "MOLECULE", "EXTERNAL", "DLVINPUT"]

        # second mandatory block: basis set
        # we use the option to end the geometry section with BASISSET instead of END,
        # because that is much more like PLAMS and way nicer to code.
        # Use CUSTOM keyword to give custom basis sets, see Crystal manual
        _basis_keys = ["BASISSET"]

        # third mandatory block: Hamiltonian and SCF control
        # this needs to be after geometry and basis and this is also the last block
        # so everything else goes in here
        _option_keys = ["OPTIONS"]

        def parse(key, value):
            ret = ""

            key = key.upper()
            if isinstance(value, Settings):
                if key in _geom_keys:
                    # geometry block should be available as a list of lines
                    if not isinstance(value, list):
                        raise PlamsError("Geometry block does not support subblocks")

                    # add geometry lines
                    if "_h" in value:
                        ret += "{}\n".format(value["_h"].upper())
                    i = 1
                    while ("_" + str(i)) in value:
                        ret += parse("", value["_" + str(i)])
                        i += 1

                    for el in value:
                        if not el.startswith("_"):
                            ret += parse(el, value[el])

                elif key in _basis_keys:
                    for subkey, item in value:
                        # we have a CUSTOM basis set
                        if not subkey.upper() == "CUSTOM":
                            raise PlamsError("BasisSet block only supports subblock CUSTOM")
                        if "_h" in item:
                            ret += "{}\n".format(item["_h"].upper())
                        i = 1
                        while ("_" + str(i)) in item:
                            ret += parse("", item["_" + str(i)])
                            i += 1

                        for el in item:
                            if not el.startswith("_"):
                                ret += parse("", item[el])
                        ret += "END\n"

                # option block has no start key
                elif key in _option_keys:
                    if "_h" in value:
                        ret += "{}\n".format(value["_h"].upper())
                    i = 1
                    while ("_" + str(i)) in value:
                        ret += parse("", value["_" + str(i)])
                        i += 1

                    for el in value:
                        if not el.startswith("_"):
                            ret += parse(el, value[el])
                else:
                    # one line per item
                    ret += "{}\n".format(key)
                    if "_h" in value:
                        ret += "{}\n".format(value["_h"].upper())
                    i = 1
                    while ("_" + str(i)) in value:
                        ret += parse("", value["_" + str(i)])
                        i += 1

                    for el in value:
                        if not el.startswith("_"):
                            ret += parse(el, value[el])
                    ret += "END\n"

            elif isinstance(value, str) and key in [*_geom_keys, *_basis_keys, *_option_keys]:
                ret += "{}\n".format(value.upper())

            elif isinstance(value, list):
                if key != "":
                    ret += "{}\n".format(key)
                for el in value:
                    ret += "{}\n".format(str(el).upper())

            elif key == "":
                ret += "{}\n".format(str(value).upper())

            elif value == "" or value is True:
                ret += "{}\n".format(key)

            elif value is False:
                pass

            else:
                ret += "{}\n{}\n".format(key, str(value).upper())

            return ret

        inp = ""
        use_molecule = ("ignore_molecule" not in self.settings) or (self.settings.ignore_molecule == False)
        if use_molecule:
            self._parsemol()

        # we need a certain ordering, so make a copy of the settings instance
        # and convert all first-level keys to uppercase
        tmp = Settings()
        for key in self.settings.input:
            tmp[key.upper()] = self.settings.input[key]

        # check for three blocks
        if not any([x in _geom_keys for x in tmp]):
            raise PlamsError("One geometry block is necessary for a Crystal Job")
        if not any([x in _basis_keys for x in tmp]):
            raise PlamsError("BasisSet block is necessary for a Crystal Job")

        # first title
        inp += "{}\n".format(self.name)

        # geometry block next
        for item in _geom_keys:
            if item in tmp:
                inp += parse(item, tmp[item])
                # do not close block, it is closed by the BASISSET keyword following
                del tmp[item]

        # basis set block next
        for item in _basis_keys:
            if item in tmp:
                inp += "BASISSET\n"
                inp += parse(item, tmp[item])
                del tmp[item]

        # everything else now
        for item in tmp:
            inp += parse(item, tmp[item])

        inp += "END"
        return inp

    def _parsemol(self):
        if "ase" in Molecule._writeformat:
            # ASE has a write function for Crystal coordinate files, use that if possible
            filename = opj(self.path, "fort.34")
            # The Crystal IO of ASE only supports filenames, not descirptors right now
            self.molecule.writease(filename)
            self.settings.input.external = True
        else:
            raise PlamsError(
                "Crystal Interface has no builtin Molecule support, install ASE or use function crystalMol2Conf() and set self.settings.ignore_molecule. See Doc for details."
            )

    def get_runscript(self):
        """
        Run Crystal.
        """
        ret = self._command
        ret += " < " + self._filename("inp")
        if self.settings.runscript.stdout_redirect:
            ret += " >" + self._filename("out")
        ret += "\n\n"
        return ret

    def check(self):
        """
        Look for the normal termination signal in output. Note, that does not mean your calculation was successful!
        """
        termination = self.results.grep_output("TERMINATION")
        return len(termination) == 1


def mol2CrystalConf(molecule):
    """
    Call this function to create a CRYSTAL-type input of your structure:

    Returns a given |Molecule| object as a geomkey and a list of strings that can be used to create a Settings instance for Crystal.

            >>> print(crystalMol2Conf(mol))
            'GEOMKEY', ['0 0 0', '1', 'lattice', 'nAtoms', 'ElementNumber1 X1 Y1 Z1','ElementNumber2 X2 Y2 Z2', ...]

    - IFLAG,IFHR and IFSO are returned as 0,0,0 by default with Symmetry group P1 (number 1). This should allow most calculations to run. The user needs to change them if he wants to take advantage of symmetry.
    - The geometry key is guessed from the number of lattice vectors. For special stuff change it by hand.
    - The number of lattice vectors in the given molecule should correspond to the dimensionality of the system. Do not fill them with zeros or unit vectors, this will result in a 3D-Periodic system with wrong fractional coordinates. So stick with the standard PLAMS way of doing things.
    """
    geomList = []

    # add lattice keyword
    if len(molecule.lattice) == 0:
        geomKey = "MOLECULE"
    elif len(molecule.lattice) == 1:
        geomKey = "POLYMER"
    elif len(molecule.lattice) == 2:
        geomKey = "SLAB"
    elif len(molecule.lattice) == 3:
        geomKey = "CRYSTAL"

    # add line for IFLAG,IFHR,IFSO: 0 0 0 always works with P1, only for CRYSTAL
    if geomKey == "CRYSTAL":
        geomList.append("0 0 0")
    # add a line for space group, assume P1 Symmetry because this always works
    geomList.append("1")

    # add lattice information: vector lengths and angles (if there is any)
    lengths = []
    angles = []
    lattice = molecule.lattice[:]
    for vec in lattice:
        lengths.append(np.linalg.norm(vec))
    if len(lattice) == 2:
        v1_u = lattice[0] / np.linalg.norm(lattice[0])
        v2_u = lattice[1] / np.linalg.norm(lattice[1])
        angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) / np.pi * 180.0
        angles.append(angle)
    elif len(lattice) == 3:
        for first in range(0, 3):
            second = first - 1
            third = first - 2
            v1_u = lattice[second] / np.linalg.norm(lattice[second])
            v2_u = lattice[third] / np.linalg.norm(lattice[third])
            angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) / np.pi * 180.0
            angles.append(angle)
    if len(lattice) > 0:
        geomList.append((("{} ") * (len(lengths) + len(angles))).format(*lengths, *angles))

    # add number of atoms
    geomList.append(str(len(molecule)))

    # now add atoms in fractional coordinates or real coords if no lattice
    # transpose lattice
    nDim = len(lattice)
    if nDim > 0:
        if nDim == 1:
            lattice.append(np.array([0.0, 1.0, 0.0]))
            lattice.append(np.array([0.0, 0.0, 1.0]))
        elif nDim == 2:
            lattice.append(np.array([0.0, 0.0, 1.0]))
        # transpose lattice, since PLAMS saves vectors as rows, we need it as columns
        lattice = np.transpose(lattice)
        latticeMatInv = np.linalg.inv(lattice)
        for atom in molecule:
            atomVec = np.dot(latticeMatInv, np.array(atom.coords)).tolist()
            # make the crystal coordinates go from 0 to 1
            for i in range(0, nDim):
                if atomVec[i] < 0.0:
                    while atomVec[i] < 0.0:
                        atomVec[i] += 1.0
                if atomVec[i] > 1.0:
                    while atomVec[i] > 1.0:
                        atomVec[i] -= 1.0
            geomList.append("{:<2}  {:>14f} {:>14f} {:>14f}".format(atom.atnum, *atomVec))
    else:
        for atom in molecule:
            geomList.append("{:<2}  {:}".format(atom.atnum, atom.str(symbol=False)))

    return geomKey, geomList

"""
This file contains a class that holds info about a ForceField patch file.
It is often necessary to create multiple patch files and combine them into a single one,
to be pased as input to an AMSJob.
This class can do that.
"""

__all__ = ["ForceFieldPatch", "forcefield_params_from_kf"]


class ForceFieldPatch:
    """
    Class representing an Amber format force field patch file, as created by AMS
    """

    def __init__(self, text=None):
        """
        Creates an instance of the class ForceFieldPatch
        """
        # Set all instance variables
        self.comment = ""

        self.types = []
        self.bondtypes = []
        self.angletypes = []
        self.dihedraltypes = []
        self.impropertypes = []
        self.ljtypes = []

        self.typelines = []
        self.bondlines = []
        self.anglelines = []
        self.dihedrallines = []
        self.improperlines = []
        self.ljlines = []

        # Read the sections, and set the parameters
        if text is not None:
            lines = [l + "\n" for l in text.split("\n")]
            self.comment = lines[0]
            self._set_types(lines)
            self._set_bonds(lines)
            self._set_angles(lines)
            self._set_dihedrals(lines)
            self._set_impropers(lines)
            self._set_ljparams(lines)

    def get_text(self):
        """
        Create the full patch text
        """
        block = ""
        if len(self.types) == 0:
            return block
        block += self.comment
        block += "".join(self.typelines)
        block += "\n\n"
        block += "".join(self.bondlines)
        block += "\n"
        block += "".join(self.anglelines)
        block += "\n"
        block += "".join(self.dihedrallines)
        block += "\n"
        block += "".join(self.improperlines)
        block += "\n\n\n"
        block += "MOD4      RE\n"
        block += "".join(self.ljlines)
        block += "\nEND\n"
        return block

    def copy(self):
        """
        Returns a copy of self
        """
        ret = self.__class__()
        ret.comment = self.comment

        ret.types = self.types.copy()
        ret.bondtypes = self.bondtypes.copy()
        ret.angletypes = self.angletypes.copy()
        ret.dihedraltypes = self.dihedraltypes.copy()
        ret.impropertypes = self.impropertypes.copy()
        ret.ljtypes = self.ljtypes.copy()

        ret.typelines = self.typelines.copy()
        ret.bondlines = self.bondlines.copy()
        ret.anglelines = self.anglelines.copy()
        ret.dihedrallines = self.dihedrallines.copy()
        ret.improperlines = self.improperlines.copy()
        ret.ljlines = self.ljlines.copy()

        return ret

    def clear(self):
        """
        Empty self
        """
        self.comment = ""
        self._set_types([])
        self._set_bonds([])
        self._set_angles([])
        self._set_dihedrals([])
        self._set_impropers([])
        self._set_ljparams([])

    def __len__(self):
        """
        Returns the size of the patch
        """
        return len(self.types)

    def __str__(self):
        """
        Returns the patch as text
        """
        return self.get_text()

    def __add__(self, other):
        """
        Combine two patch files
        """
        ret = self.copy()
        if len(ret.comment) == 0:
            ret.comment = other.comment

        # Combine the lines from each section
        ret.typelines += [other.typelines[i] for i, t in enumerate(other.types) if not t in self.types]
        ret.bondlines += [other.bondlines[i] for i, b in enumerate(other.bondtypes) if not b in self.bondtypes]
        ret.anglelines += [other.anglelines[i] for i, t in enumerate(other.angletypes) if not t in self.angletypes]
        ret.dihedrallines += [
            other.dihedrallines[i] for i, t in enumerate(other.dihedraltypes) if not t in self.dihedraltypes
        ]
        ret.improperlines += [
            other.improperlines[i] for i, t in enumerate(other.impropertypes) if not t in self.impropertypes
        ]
        ret.ljlines += [other.ljlines[i] for i, t in enumerate(other.ljtypes) if not t in self.types]

        # Now adjust all the parameter info (not the lines)
        ret.bondtypes += [b for b in other.bondtypes if not b in self.bondtypes]
        ret.angletypes += [a for a in other.angletypes if not a in self.angletypes]
        ret.dihedraltypes += [d for d in other.dihedraltypes if not d in self.dihedraltypes]
        ret.impropertypes += [imp for imp in other.impropertypes if not imp in self.impropertypes]
        ret.ljtypes += [t for t in other.ljtypes if not t in self.types]
        ret.types += [t for t in other.types if not t in self.types]

        return ret

    def read_from_kf(self, kf):
        """
        Read patch infor from kf
        """
        if len(self) > 0:
            self.clear()
        npatches = kf.read("AMSResults", "Config.nPatches")
        patch = ForceFieldPatch()
        if npatches > 0:
            patchtext = kf.read("AMSResults", "Config.FFPatch(1)")
            patch += ForceFieldPatch(patchtext)

        for key in vars(patch):
            if "type" in key or "lines" in key or "comment" in key:
                self.__dict__[key] = patch.__dict__[key]

    def write_to_kf(self, kf):
        """
        Write the patch info to KF
        """
        kf.write("AMSResults", "Config.nPatches", 1)
        kf.write("AMSResults", "Config.FFPatch(1)", str(self))

    def _set_types(self, lines):
        """
        Read the atom types
        """
        types = []
        typelines = []
        for line in lines[1:]:
            words = line.split()
            if len(words) == 0:
                break
            types.append(words[0])
            typelines.append(line)
        self.types = types
        self.typelines = typelines

    def _set_bonds(self, lines):
        """
        Set the bond parameters from the list of lines
        """
        b, blines = self._read_atoms(lines, nats=2)
        self.bondtypes = b
        self.bondlines = blines

    def _set_angles(self, lines):
        """
        Set the angle parameters from list of lines
        """
        a, alines = self._read_atoms(lines, nats=3)
        self.angletypes = a
        self.anglelines = alines

    def _set_dihedrals(self, lines):
        """
        Set the dihedral parameters from list of lines
        """
        d, dlines = self._read_atoms(lines, nats=4)
        self.dihedraltypes = d
        self.dihedrallines = dlines

    def _set_impropers(self, lines):
        """
        Set improper parameters from list of lines
        """
        imp, implines = self._read_atoms(lines, nats=4, improper=True)
        self.impropertypes = imp
        self.improperlines = implines

    def _set_ljparams(self, lines):
        """
        Set the Lennard-Jones paramters from list of lines
        """
        lj, ljlines = self._read_LJtypes(lines)
        self.ljtypes = lj
        self.ljlines = ljlines

    @staticmethod
    def _read_atoms(lines, nats=2, improper=False):
        """
        Read the atoms from the patch lines for bond, angle, dihedral parameters

        * ``nats`` - Integer: 2 for bond info, 3 for angl info, 4 for dihedral info
        """
        step = 3
        positions = [2 + (i * step) for i in range(nats)]
        intervals = [((i * step), (i * step) + step - 1) for i in range(nats)]

        atomlist = []
        line_list = []
        for line in lines[1:]:
            if len(line) - 1 < max(positions):
                continue
            symbols = [line[i] for i in positions]
            if not False in [v == "-" for v in symbols[:-1]] and symbols[-1] == " ":
                # Disregard impropers
                found_improper = nats == 4 and line[14] == " "
                if not improper:
                    if found_improper:
                        continue
                else:
                    if not found_improper:
                        continue
                # Extract the atoms, and sort them if they are not improper angle atoms
                atoms = [line[t[0] : t[1]] for t in intervals]
                if not improper:
                    if sorted([atoms[0], atoms[-1]]) != [atoms[0], atoms[-1]]:
                        atoms = atoms[::-1]
                atomlist.append(atoms)
                line_list.append(line)
        return atomlist, line_list

    @staticmethod
    def _read_LJtypes(lines):
        """
        Read the LJ info from the patch file lines
        """
        start = False
        ljtypes = []
        ljlines = []
        for line in lines[1:]:
            if "MOD4      RE" in line:
                start = True
                continue
            if "END" in line:
                start = False
            if start:
                words = line.split()
                if len(words) == 0:
                    continue
                ljtypes.append(words[0])
                ljlines.append(line)
        return ljtypes, ljlines


def forcefield_params_from_kf(kf):
    """
    Read the parameters from kf
    """
    charges = kf.read("AMSResults", "Charges")
    alltypes = kf.read("AMSResults", "AtomTyping.atomTypes").split("\x00")
    indices = kf.read("AMSResults", "AtomTyping.atomIndexToType")
    types = [alltypes[i - 1] for i in indices]

    # Read the force field patch
    patch = ForceFieldPatch()
    patch.read_from_kf(kf)
    if len(patch) == 0:
        return charges, types, None
    return charges, types, patch

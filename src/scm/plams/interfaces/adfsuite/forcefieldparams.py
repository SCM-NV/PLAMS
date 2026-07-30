"""
This file contains a class that holds info about a ForceField patch file.
It is often necessary to create multiple patch files and combine them into a single one,
to be pased as input to an AMSJob.
This class can do that.
"""

from typing import Optional, Tuple, Sequence, List

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from scm.plams.tools.kftools import KFFile

__all__ = ["ForceFieldPatch", "forcefield_params_from_kf"]


class ForceFieldPatchError(Exception):
    """
    Error in handling the ForceFieldPatch object
    """


class ForceFieldPatch:
    """
    Class representing an Amber format force field patch file, as created by AMS
    """

    def __init__(self, text: Optional[str] = None):
        """
        Creates an instance of the class ForceFieldPatch
        """
        # Set all instance variables
        self.comment = ""

        self.types: List[str] = []
        self.bondtypes: List[List[str]] = []
        self.angletypes: List[List[str]] = []
        self.dihedraltypes: List[List[str]] = []
        self.impropertypes: List[List[str]] = []
        self.ljtypes: List[str] = []

        self.typelines: List[str] = []
        self.bondlines: List[str] = []
        self.anglelines: List[str] = []
        self.dihedrallines: List[str] = []
        self.improperlines: List[str] = []
        self.ljlines: List[str] = []

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

    @classmethod
    def from_frcmod(cls, prm_string: str) -> "ForceFieldPatch":
        """
        From text of an frcmod file
        """
        patch_string = cls._convert_frcmod(prm_string)
        ret = ForceFieldPatch(patch_string)
        return ret

    def get_text(self) -> str:
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

    def copy(self) -> "ForceFieldPatch":
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

    def clear(self) -> None:
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

    def __len__(self) -> int:
        """
        Returns the size of the patch
        """
        return len(self.types)

    def __str__(self) -> str:
        """
        Returns the patch as text
        """
        return self.get_text()

    def __add__(self, other: "ForceFieldPatch") -> "ForceFieldPatch":
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

    def read_from_kf(self, kf: "KFFile", ipatch: int = 0) -> None:
        """
        Read patch infor from kf
        """
        if len(self) > 0:
            self.clear()
        npatches = kf.read_int("AMSResults", "Config.nPatches")
        patch = ForceFieldPatch()
        if npatches > 0:
            patchtext = kf.read_string("AMSResults", f"Config.FFPatch({ipatch + 1})")
            patch += ForceFieldPatch(patchtext)

        for key in vars(patch):
            if "type" in key or "lines" in key or "comment" in key:
                self.__dict__[key] = patch.__dict__[key]

    def get_fragment(self, types: Sequence[str]) -> "ForceFieldPatch":
        """
        Get the patch object for a subset of atom types

        Note: This does not handle starred types (C*)
        """

        def get_data(
            datatype: str, system_types: List[List[str]], system_typelines: List[str]
        ) -> Tuple[List[List[str]], List[str]]:
            """
            Get the relevant data
            """
            fragment_types: List[List[str]] = []
            fragment_typelines: List[str] = []
            for it, connected_types in enumerate(system_types):
                contypes = [t for t in connected_types if t in types]
                if len(contypes) == 0:
                    continue
                elif len(contypes) < len(connected_types):
                    continue
                    # raise ForceFieldPatchError("Cannot extract fragment: Dependency in {datatype} [{" ".join(contypes)}]")
                else:
                    fragment_types.append(connected_types)
                    fragment_typelines.append(system_typelines[it])
            return fragment_types, fragment_typelines

        ret = ForceFieldPatch()

        for atomtype in types:
            if atomtype not in self.types:
                raise ForceFieldPatchError(f"Requested type {atomtype} not present in parent patch")

        indices = [i for i, t in enumerate(self.types) if t in types]
        ret.types = [self.types[i] for i in indices]
        ret.typelines = [self.typelines[i] for i in indices]

        indices = [i for i, t in enumerate(self.ljtypes) if t in types]
        ret.ljtypes = [self.ljtypes[i] for i in indices]
        ret.ljlines = [self.ljlines[i] for i in indices]

        ret.bondtypes, ret.bondlines = get_data("bond", self.bondtypes, self.bondlines)
        ret.angletypes, ret.anglelines = get_data("angle", self.angletypes, self.anglelines)
        ret.dihedraltypes, ret.dihedrallines = get_data("dihedral", self.dihedraltypes, self.dihedrallines)
        ret.impropertypes, ret.improperlines = get_data("improper", self.impropertypes, self.improperlines)

        return ret

    def write_to_kf(self, kf: "KFFile") -> None:
        """
        Write the patch info to KF
        """
        kf.write("AMSResults", "Config.nPatches", 1)
        kf.write("AMSResults", "Config.FFPatch(1)", str(self))

    def _set_types(self, lines: Sequence[str]) -> None:
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

    def _set_bonds(self, lines: Sequence[str]) -> None:
        """
        Set the bond parameters from the list of lines
        """
        b, blines = self._read_atoms(lines, nats=2)
        self.bondtypes = b
        self.bondlines = blines

    def _set_angles(self, lines: Sequence[str]) -> None:
        """
        Set the angle parameters from list of lines
        """
        a, alines = self._read_atoms(lines, nats=3)
        self.angletypes = a
        self.anglelines = alines

    def _set_dihedrals(self, lines: Sequence[str]) -> None:
        """
        Set the dihedral parameters from list of lines
        """
        d, dlines = self._read_atoms(lines, nats=4)
        self.dihedraltypes = d
        self.dihedrallines = dlines

    def _set_impropers(self, lines: Sequence[str]) -> None:
        """
        Set improper parameters from list of lines
        """
        imp, implines = self._read_atoms(lines, nats=4, improper=True)
        self.impropertypes = imp
        self.improperlines = implines

    def _set_ljparams(self, lines: Sequence[str]) -> None:
        """
        Set the Lennard-Jones paramters from list of lines
        """
        lj, ljlines = self._read_LJtypes(lines)
        self.ljtypes = lj
        self.ljlines = ljlines

    @staticmethod
    def _read_atoms(lines: Sequence[str], nats: int = 2, improper: bool = False) -> Tuple[List[List[str]], List[str]]:
        """
        Read the atoms from the patch lines for bond, angle, dihedral parameters

        * ``nats`` - Integer: 2 for bond info, 3 for angl info, 4 for dihedral info
        """
        step = 3
        positions = [2 + (i * step) for i in range(nats)]
        intervals = [((i * step), (i * step) + step - 1) for i in range(nats)]

        atomlist: List[List[str]] = []
        line_list: List[str] = []
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
    def _read_LJtypes(lines: Sequence[str]) -> Tuple[List[str], List[str]]:
        """
        Read the LJ info from the patch file lines
        """
        start = False
        ljtypes: List[str] = []
        ljlines: List[str] = []
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

    @staticmethod
    def _convert_frcmod(prm_string: str) -> str:
        """
        Converts frcmod format to AMBER .dat format
        """
        sections = ForceFieldPatch._get_frcmod_sections(prm_string)

        # First find all atom types
        atomTypes = []
        for txt in sections[2:]:
            lines = txt.split("\n")
            for line in lines:
                single_term_atom_types = []
                if len(line) > 0:
                    single_term_atom_types = ForceFieldPatch._atom_types_from_frcmodline(line)
                    # Add these types to atomTypes, if they are not yet in there
                for atom_type in single_term_atom_types:
                    if not atom_type in atomTypes:
                        atomTypes.append(atom_type)

        # Then add them to the mass section
        # First: Is there already something in there?
        # If so, remove that from the atomTypes list
        remainingAtomTypes = ForceFieldPatch._non_present_atom_types(sections[0], atomTypes)
        # Now add the not yet present elements to the atomTypes list (with made up values)
        for atom_type in remainingAtomTypes:
            sections[0] += atom_type + " 1.008" + "         0.000"
            sections[0] += "               Added by AMS\n"

        # Do the same for the LJ section
        remainingAtomTypes = ForceFieldPatch._non_present_atom_types(sections[5], atomTypes)
        for atom_type in remainingAtomTypes:
            sections[5] += "  " + atom_type + "          0.0000  0.000"
            sections[5] += "             Added by AMS\n"

        # Now write all the sections in the appropriate format
        block = "Patch to GAFF force field produced by prmcheck\n"
        for i, section in enumerate(sections):
            if i == 5:
                block += "MOD4      RE\n"
            lines = section.split("\n")
            for line in lines:
                if len(line) > 0:
                    block += line + "\n"
            block += "\n"
            # There is a hydrophilic section before the bond-pars that can be skipped
            if i == 0:
                block += "\n"
            # Here the .dat reader calls a Dummy routine twice that reads an empty line (9 and 10)
            if i == 4:
                block += "\n\n"
        block += "END\n"

        return block

    @staticmethod
    def _get_frcmod_sections(txt: str) -> List[str]:
        """
        Separate the sections into a list
        """
        sections = 6 * [""]
        parts = txt.split("NONBON\n")
        sections[5] = parts[1]  # nonbonded text
        parts = parts[0].split("IMPROPER\n")
        sections[4] = parts[1]  # improper text
        parts = parts[0].split("DIHE\n")
        sections[3] = parts[1]  # dihedral text
        parts = parts[0].split("ANGLE\n")
        sections[2] = parts[1]  # angle text
        parts = parts[0].split("BOND\n")
        sections[1] = parts[1]  # bond text
        parts = parts[0].split("MASS\n")
        sections[0] = parts[1]  # mass text
        return sections

    @staticmethod
    def _atom_types_from_frcmodline(line: str) -> List[str]:
        """
        Read all the atomtypes from a frcmod line
        """
        # I have to look for a dash in specific places to find out the
        # number of types here and add them
        ntypes = 0
        if line[2:3] == "-":
            ntypes = 2
        if line[5:6] == "-":
            ntypes = 3
        if line[8:9] == "-":
            ntypes = 4
        if line[0:2] == "  " and line[2:4] != "  " and line[6:8] == "  ":
            ntypes = 1
        if ntypes == 1:
            # LJ section
            words = [line[2:4]]
        else:
            words = [line[k * 3 : k * 3 + 2] for k in range(ntypes)]
        return words

    @staticmethod
    def _non_present_atom_types(section: str, atomTypes: List[str]) -> List[str]:
        """
        First: Is there already something in there?
        If so, remove that from the atomTypes list
        """
        remainingAtomTypes = atomTypes[:]
        lines = section.split("\n")
        for line in lines:
            words = line.split()
            if len(words) == 0:
                continue
            atomType = words[0]
            remainingAtomTypes = [t for t in remainingAtomTypes if t != atomType]
        return remainingAtomTypes


def forcefield_params_from_kf(kf: "KFFile") -> Tuple[List[float], List[str], Optional[ForceFieldPatch]]:
    """
    Read the parameters from kf
    """
    charges = kf.read_reals("AMSResults", "Charges")
    alltypes = kf.read_string("AMSResults", "AtomTyping.atomTypes").split("\x00")
    indices = kf.read_ints("AMSResults", "AtomTyping.atomIndexToType")
    types = [alltypes[i - 1] for i in indices]

    # Read the force field patches
    npatches = kf.read_int("AMSResults", "Config.nPatches")
    if npatches == 0:
        return charges, types, None
    patch = ForceFieldPatch()
    for i in range(npatches):
        p = ForceFieldPatch()
        p.read_from_kf(kf, i)
        patch += p
    return charges, types, patch

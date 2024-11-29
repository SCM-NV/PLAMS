import numpy as np
from typing import Optional, Tuple, TYPE_CHECKING

from scm.plams.core.errors import MoleculeError
from scm.plams.core.settings import Settings

__all__ = ["Bond"]

if TYPE_CHECKING:
    from scm.plams.mol.atom import Atom


class Bond:
    """A class representing a bond between two atoms.

    An instance of this class has the following attributes:

    *   ``atom1`` and ``atom2`` -- two instances of |Atom| that form this bond
    *   ``order`` -- order of the bond. It is either an integer number or the floating point value stored in ``Bond.AR``, indicating an aromatic bond
    *   ``mol`` -- |Molecule| this bond belongs to
    *   ``properties`` -- |Settings| instance storing all other  information about this bond (initially it is populated with *\*\*other*)

    .. note::

        Newly created bond is **not** added to ``atom1.bonds`` or ``atom2.bonds``. Storing information about |Bond| in |Atom| is relevant only in the context of the whole |Molecule|, so this information is updated by :meth:`~Molecule.add_bond`.

    """

    AR = 1.5

    def __init__(self, atom1=None, atom2=None, order=1, mol=None, **other):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        self.mol = mol
        self.properties = Settings(other)

    def __str__(self):
        """Return a string representation of this bond."""
        return "({})--{:1.1f}--({})".format(str(self.atom1).strip(), self.order, str(self.atom2).strip())

    def __iter__(self):
        """Iterate over bonded atoms (``atom1`` first, then ``atom2``)."""
        yield self.atom1
        yield self.atom2

    def is_aromatic(self) -> bool:
        """Check if this bond is aromatic."""
        return self.order == Bond.AR

    def length(self, unit: str = "angstrom") -> float:
        """Return bond length, expressed in *unit*."""
        return self.atom1.distance_to(self.atom2, result_unit=unit)

    def as_vector(self, start: Optional["Atom"] = None, unit: str = "angstrom") -> Tuple[float, float, float]:
        """Return a vector between two atoms that form this bond. *start* can be used to indicate which atom should be the beginning of that vector. If not specified, ``self.atom1`` is used. Returned value if a tuple of length 3, expressed in *unit*."""
        if start:
            if start not in self:
                raise MoleculeError("Bond.as_vector: given atom is not a part of this bond")
            a, b = start, self.other_end(start)
        else:
            a, b = self.atom1, self.atom2
        return a.vector_to(b, result_unit=unit)

    def other_end(self, atom: "Atom") -> "Atom":
        """Return the atom on the other end of this bond with respect to *atom*. *atom* has to be one of the atoms forming this bond, otherwise an exception is raised."""
        if atom is self.atom1:
            return self.atom2
        elif atom is self.atom2:
            return self.atom1
        else:
            raise MoleculeError("Bond.other_end: invalid atom passed")

    def resize(self, moving_atom: "Atom", length: float, unit: str = "angstrom") -> None:
        """Change the length of this bond to *length* expressed in *unit* by moving *moving_atom*.

        *moving_atom* should be one of the atoms that form this bond. This atom is moved along the bond axis in such a way that new bond length equals *length*. If this bond is a part of a |Molecule| the whole part connected to *moving_atom* is moved.

        .. note::

            Calling this method on a bond that forms a ring within a molecule raises a |MoleculeError|.

        """

        if self.mol:
            self.mol.resize_bond(self, moving_atom, length, unit)
        else:
            bond_v = np.array(self.as_vector(start=moving_atom))
            trans_v = (1 - length / self.length(unit)) * bond_v
            moving_atom.translate(trans_v)

    def rotate(self, moving_atom: "Atom", angle: float, unit: str = "radian"):
        """Rotate part of the molecule containing *moving_atom* along axis defined by this bond by an *angle* expressed in *unit*.

        Calling this method makes sense only if this bond is a part of a |Molecule|. *moving_atom* should be one of the atoms that form this bond and it indicates which part of the molecule is rotated. A positive value of *angle* denotes counterclockwise rotation (when looking along the bond, from the stationary part of the molecule).

        .. note::

            Calling this method on a bond that forms a ring raises a |MoleculeError|.

        """
        if self.mol:
            self.mol.rotate_bond(self, moving_atom, angle, unit)

    def has_cell_shifts(self) -> bool:
        return (
            "suffix" in self.properties
            and isinstance(self.properties.suffix, str)
            and self.properties.suffix != ""
            and not self.properties.suffix.isspace()
        )

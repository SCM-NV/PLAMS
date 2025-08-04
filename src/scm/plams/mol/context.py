from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from scm.plams.mol.molecule import Molecule


class AsArrayContext:
    """A context manager for temporary inter-converting between PLAMS molecules and numpy arrays."""

    def __init__(self, mol: Molecule):
        self._atoms = mol.atoms
        self._from_array = mol.from_array

    def __call__(self, atom_subset=None):
        """Return cartesian coordinates of this molecule's atoms as a numpy array.

        *atom_subset* argument can be used to specify only a subset of atoms, it should be an iterable container with atoms belonging to this molecule.

        Returned value is a n*3 numpy array where n is the number of atoms in the whole molecule, or in *atom_subset*, if used.

        Alternatively, this function can be used in conjunction with the ``with`` statement,
        which automatically calls :meth:`Molecule.from_array` upon exiting the context manager.
        Note that the molecules' coordinates will be updated based on the array that was originally returned,
        so creating and operating on a copy thereof will not affect the original molecule.

        .. code-block:: python

            >>> from scm.plams import Molecule
            >>> mol = Molecule(...)
            >>> with mol.as_array as xyz_array:
            >>>     xyz_array += 5.0
            >>>     xyz_array[0] = [0, 0, 0]
            # Or equivalently
            >>> xyz_array = mol.as_array()
            >>> xyz_array += 5.0
            >>> xyz_array[0] = [0, 0, 0]
            >>> mol.from_array(xyz_array)
        """
        atom_subset = atom_subset or self._atoms

        try:
            at_len = len(atom_subset)
        except TypeError:  # atom_subset is an iterator
            count = -1
            shape = -1, 3
        else:
            count = at_len * 3
            shape = at_len, 3

        atom_iterator = itertools.chain.from_iterable(at.coords for at in atom_subset)
        xyz_array = np.fromiter(atom_iterator, count=count, dtype=float)
        xyz_array.shape = shape
        return xyz_array

    def __enter__(self, atom_subset=None):
        """Enter the context manager; return the Cartesian coordinate array."""
        self._atom_subset = atom_subset
        self._xyz_array = self.__call__(atom_subset=atom_subset)
        return self._xyz_array

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context manager; update the ``mol`` coordinates with those from the Cartesian coordinate array."""
        self._from_array(self._xyz_array, self._atom_subset)

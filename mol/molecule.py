import copy
import heapq
import io
import itertools
import math
import os
from collections import OrderedDict

import numpy as np

from scm.plams.core.errors import FileError, MoleculeError, PTError
from scm.plams.core.functions import log
from scm.plams.core.private import parse_action, smart_copy
from scm.plams.core.settings import Settings
from scm.plams.mol.atom import Atom
from scm.plams.mol.bond import Bond
from scm.plams.mol.context import AsArrayContext
from scm.plams.mol.pdbtools import PDBHandler, PDBRecord
from scm.plams.tools.geometry import axis_rotation_matrix, cell_angles, cell_lengths, distance_array, rotation_matrix
from scm.plams.tools.kftools import KFFile
from scm.plams.tools.periodic_table import PT
from scm.plams.tools.units import Units

input_parser_available = "AMSBIN" in os.environ
from typing import List, Optional, Tuple, overload

__all__ = ["Molecule"]


class Molecule:
    """A class representing the molecule object.

    An instance of this class has the following attributes:

    *   ``atoms`` -- list of |Atom| objects that belong to the molecule
    *   ``bonds`` -- list of |Bond| objects between atoms listed in ``atoms``
    *   ``lattice`` -- list of lattice vectors in case of periodic structures
    *   ``properties`` -- |Settings| instance storing all other information about the molecule

    .. note::

        Each |Atom| in ``atoms`` list and each |Bond| in ``bonds`` list has a reference to the parent molecule. Moreover, each atom stores the list of bonds it's a part of and each bond stores references to atoms it bonds. That creates a complex net of references between objects that are part of a molecule. Consistency of this data is crucial for proper functioning of many methods. Because of that it is advised not to modify contents of ``atoms`` and ``bonds`` by hand. When you need to alter your molecule, methods :meth:`add_atom`, :meth:`delete_atom`, :meth:`add_bond` and :meth:`delete_bond` can be used to ensure that all these references are updated properly.

    Creating a |Molecule| object for your calculation can be done in several ways. You can start with an empty molecule and manually add all atoms (and bonds, if needed)::

        mol = Molecule()
        mol.add_atom(Atom(atnum=1, coords=(0,0,0)))
        mol.add_atom(Atom(atnum=1, coords=(d,0,0)))

    This approach can be useful for building small molecules, especially if you wish to parametrize some of atomic coordinates (like in :ref:`simple_example`), but in general it's not very practical.
    If coordinates and atom numbers are available, instantiation can be done by passing a value to the `positions`, `numbers` and optionally the `lattice` arguments::

        xyz     = np.random.randn(10,3) # 10 atoms, 3 coordinates per atom
        numbers = 10*[6] # 10 carbon atoms. If left None, will initialize to dummy atoms
        lattice = [[1,2,3], [1,2,3]] # lattice should have a shape of {1,2,3}x3
        mol     = Molecule(positions=xyz, numbers=numbers, lattice=lattice)

    Alternatively, one can import atomic coordinates from some external file::

        mol = Molecule('xyz/Benzene.xyz')

    The constructor of a |Molecule| object accepts four arguments that can be used to supply this information from a file in your filesystem. *filename* should be a string with a path (absolute or relative) to such a file. *inputformat* describes the format of the file. Currently, the following formats are supported: ``xyz``, ``mol``, ``mol2`` and ``pdb``. If *inputformat* is ``ase`` the file reader engine of the ASE.io module is used, enabling you to read all input formats supported by :ref:`ASEInterface`. See :meth:`read` for further details. If the *inputformat* argument is not supplied, PLAMS will try to deduce it by examining the extension of the provided file, so in most of cases it is not needed to use *inputformat*, if only the file has the proper extension. Some formats (``xyz`` and ``pdb``) allow to store more than one geometry of a particular molecule within a single file. See the respective :meth:`read` function for details how to access them. All *other* keyword arguments will be passed to the appropriate read function for the selected or determined file format.

    If a |Molecule| is initialized from an external file, the path to this file (*filename* argument) is stored in ``properties.source``. The base name of the file (filename without the extension) is kept in ``properties.name``.

    It is also possible to write a molecule to a file in one of the formats mentioned above or using the ASE.io engine. See :meth:`write` for details.

    The ``lattice`` attribute is used to store information about lattice vectors in case of periodic structures. Some job types will automatically use that data while constructing input files. ``lattice`` should be a list of up to 3 vectors (for different types of periodicity: chain, slab or bulk), each of which needs to be a list or a tuple of 3 numbers.

    Lattice vectors can be directly read from and written to ``xyz`` files using the following convention (please mind the fact that this is an unofficial extension to the XYZ format):

    .. code-block:: none

        3

            H      0.000000      0.765440     -0.008360
            O      0.000000      0.000000      0.593720
            H      0.000000     -0.765440     -0.008360
        VEC1       3.000000      0.000000      0.000000
        VEC2       0.000000      3.000000      0.000000
        VEC3       0.000000      0.000000      3.000000

    For 1D (2D) periodicity please supply only ``VEC1`` (``VEC1`` and ``VEC2``). Writing lattice vectors to ``xyz`` files can be disabled by simply resetting the ``lattice`` attribute::

        mol.lattice = []


    The detailed description of all available methods is presented below. Many of these methods require arguments that are atoms belonging to the current molecule. It can by done by using a reference to an |Atom| object present it the ``atoms`` list, but not by passing a number of an atom (its position within ``atoms`` list). Unlike some other tools, PLAMS does not use integer numbers as primary identifiers of atoms. It is done to prevent problems when atoms within a molecule are reordered or some atoms are deleted. References to |Atom| or |Bond| objects can be obtained directly from ``atoms`` or ``bonds`` lists, or with dictionary-like bracket notation::

        >>> mol = Molecule('xyz/Ammonia.xyz')
        >>> mol.guess_bonds()
        >>> print(mol)
          Atoms:
            1         H      0.942179      0.000000     -0.017370
            2         H     -0.471089      0.815951     -0.017370
            3         N      0.000000      0.000000      0.383210
            4         H     -0.471089     -0.815951     -0.017370
          Bonds:
           (1)--1.0--(3)
           (2)--1.0--(3)
           (3)--1.0--(4)
        >>> at = mol[1]
        >>> print(at)
                 H      0.942179      0.000000     -0.017370
        >>> b = mol[(1,3)]
        >>> print(b)
        (         H      0.942179      0.000000     -0.017370 )--1.0--(         N      0.000000      0.000000      0.383210 )
        >>> b = mol[(1,4)]
        >>> print(b)
        None

    .. note::

        For the purpose of ``mol[i]`` notation, the numbering of atoms within a molecule starts with 1. Negative integers can be used to access atoms enumerated in the reversed order (``mol[-1]`` for the last atom etc.)

    However, if you feel more familiar with identifying atoms by natural numbers, you can use :meth:`set_atoms_id` to equip each atom of the molecule with ``id`` attribute equal to atom's position within ``atoms`` list. This method can also be helpful to track changes in your molecule during tasks that can reorder atoms.
    """

    def __init__(self, filename=None, inputformat=None, positions=None, numbers=None, lattice=None, **other):
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        self.lattice: List[List[float]] = []
        self.properties = Settings()

        if filename is not None:
            self.read(filename, inputformat, **other)
            self.properties.source = filename
            self.properties.name = os.path.splitext(os.path.basename(filename))[0]

        elif positions is not None:
            positions = np.array(positions)
            assert positions.ndim == 2, "`Positions` must be a 2d array"
            assert positions.shape[-1] == 3, "Inner dim of `positions` must be 3"
            if numbers is None:
                numbers = len(positions) * [0]
            assert len(numbers) == len(positions), "Length or `numbers` and `positions` does not match"
            for num, xyz in zip(numbers, positions):
                self.add_atom(Atom(atnum=num, coords=xyz))
            if lattice is not None:
                lattice = np.array(lattice)
                assert lattice.ndim == 2, "`Lattice` must be a 2d array"
                assert lattice.shape[0] <= 3, "`Lattice` should be a 3x3 vector at most"
                assert lattice.shape[-1] == 3, "Inner dim of `lattice` must be 3"
                self.lattice = lattice.tolist()
        # create as_array method as an object that supports both direct calling and use as context manager
        self._as_array = AsArrayContext(self)

    # ===========================================================================
    # ==== Atoms/bonds manipulation =============================================
    # ===========================================================================

    def copy(self, atoms=None):
        """Return a copy of the molecule. The copy has atoms, bonds and all other components distinct from the original molecule (it is so called "deep copy").

        By default the entire molecule is copied. It is also possible to copy only some part of the molecule, indicated by *atoms* argument. It should be a list of atoms that belong to the molecule. If used, only these atoms, together with any bonds between them, are copied and included in the returned molecule.
        """

        if atoms is None:
            atoms = self.atoms

        # _as_array is an object that contains a reference to the Molecule object and should thus be excluded from
        # the copy. It will be recreated on Molecule.__init__
        ret = smart_copy(self, owncopy=["properties"], without=["atoms", "bonds", "_as_array"])

        bro = {}  # mapping of original to copied atoms
        for at in atoms:
            at_copy = smart_copy(at, owncopy=["properties"], without=["mol", "bonds"])
            ret.add_atom(at_copy)
            bro[at] = at_copy

        for bo in self.bonds:
            if (bo.atom1 in bro) and (bo.atom2 in bro):
                bo_copy = smart_copy(bo, owncopy=["properties"], without=["atom1", "atom2", "mol"])
                bo_copy.atom1 = bro[bo.atom1]
                bo_copy.atom2 = bro[bo.atom2]
                ret.add_bond(bo_copy)

        return ret

    def add_molecule(self, other, copy=False, margin=-1):
        """Add some *other* molecule to this one::

            protein += water

        If *copy* is ``True``, *other* molecule is copied and the copy is added to this molecule. Otherwise, *other* molecule is directly merged with this one
        The ``properties`` of this molecule are :meth:`soft_updated<scm.plams.core.settings.Settings.soft_update>` with the  ``properties`` of the *other* molecules.

        margin: float
            If <0, keep the coordinates of ``other``. If >=0, all atoms in the ``other`` molecule will have *at least* this distance (in angstrom) to all atoms in ``self``.

        """
        if margin >= 0:
            if len(self) == 0:
                dx, dy, dz = 0, 0, 0
            else:
                dx, dy, dz = np.max(self.as_array(), axis=0) - np.min(other.as_array(), axis=0) + margin
            copy = True

        other = other.copy() if copy else other

        if margin >= 0:
            other.translate([dx, dy, dz])

        self.atoms += other.atoms
        self.bonds += other.bonds
        for atom in self.atoms:
            atom.mol = self
        for bond in self.bonds:
            bond.mol = self
        self.properties.soft_update(other.properties)

    def add_atom(self, atom, adjacent=None):
        """Add a new *atom* to the molecule.

        *atom* should be an |Atom| instance that does not belong to any molecule. Bonds between the new atom and other atoms of the molecule can be automatically added based on *adjacent* argument. It should be a list describing atoms of the molecule that the new atom is connected to. Each element of *adjacent* list can either be a pair ``(Atom, order)`` to indicate new bond's order (use ``Bond.AR`` for aromatic bonds) or an |Atom| instance (a single bond is created in this case).

        Example::

            mol = Molecule() #create an empty molecule
            h1 = Atom(symbol='H', coords=(1.0, 0.0, 0.0))
            h2 = Atom(symbol='H', coords=(-1.0, 0.0, 0.0))
            o = Atom(symbol='O', coords=(0.0, 1.0, 0.0))
            mol.add_atom(h1)
            mol.add_atom(h2)
            mol.add_atom(o)
            mol.add_atom(Atom(symbol='C', coords=(0.0, 0.0, 0.0)), adjacent=[h1, h2, (o,2)])

        """
        self.atoms.append(atom)
        atom.mol = self
        if adjacent is not None:
            for adj in adjacent:
                if isinstance(adj, Atom):
                    self.add_bond(atom, adj)
                else:
                    self.add_bond(atom, *adj)

    def delete_atom(self, atom):
        """Delete an *atom* from the molecule.

        *atom* should be an |Atom| instance that belongs to the molecule. All bonds containing this atom are removed too.

        Examples::

            #delete all hydrogens
            mol = Molecule('protein.pdb')
            hydrogens = [atom for atom in mol if atom.atnum == 1]
            for i in hydrogens: mol.delete_atom(i)

        ::

            #delete first two atoms
            mol = Molecule('geom.xyz')
            mol.delete_atom(mol[1])
            mol.delete_atom(mol[1]) #since the second atom of original molecule is now the first

        """
        if atom.mol != self:
            raise MoleculeError("delete_atom: passed atom should belong to the molecule")
        try:
            self.atoms.remove(atom)
        except:
            raise MoleculeError("delete_atom: invalid argument passed as atom")
        atom.mol = None
        for b in reversed(atom.bonds):
            self.delete_bond(b)

    def add_bond(self, arg1, arg2=None, order=1):
        """Add a new bond to the molecule.

        This method can be used in two different ways. You can call it with just one argument being a |Bond| instance (other arguments are then ignored)::

            >>> b = Bond(mol[2], mol[4], order=Bond.AR) #create aromatic bond between 2nd and 4th atom
            >>> mol.add_bond(b)

        The other way is to pass two atoms (and possibly bond order) and new |Bond| object will be created automatically::

            >>> mol.add_bond(mol[2], mol[4], order=Bond.AR)

        In both cases both atoms that are bonded have to belong to the molecule, otherwise an exception is raised.
        """
        if isinstance(arg1, Atom) and isinstance(arg2, Atom):
            newbond = Bond(arg1, arg2, order=order)
        elif isinstance(arg1, Bond):
            newbond = arg1
        else:
            raise MoleculeError("add_bond: invalid arguments passed")

        if newbond.atom1.mol == self and newbond.atom2.mol == self:
            newbond.mol = self
            self.bonds.append(newbond)
            newbond.atom1.bonds.append(newbond)
            newbond.atom2.bonds.append(newbond)
        else:
            raise MoleculeError("add_bond: bonded atoms have to belong to the molecule")

    def delete_bond(self, arg1, arg2=None):
        """Delete a bond from the molecule.

        Just like :meth:`add_bond`, this method accepts either a single argument that is a |Bond| instance, or two
        arguments being instances of |Atom|. In both cases objects used as arguments have to belong to the molecule.
        """
        if isinstance(arg1, Atom) and isinstance(arg2, Atom):
            delbond = self.find_bond(arg1, arg2)
        elif isinstance(arg1, Bond):
            delbond = arg1
        else:
            raise MoleculeError("delete_bond: invalid arguments passed")
        if delbond in self.bonds:
            delbond.mol = None
            self.bonds.remove(delbond)
            delbond.atom1.bonds.remove(delbond)
            delbond.atom2.bonds.remove(delbond)

    def delete_all_bonds(self):
        """Delete all bonds from the molecule."""
        self.bonds.clear()
        for atom in self:
            atom.bonds.clear()

    def find_bond(self, atom1, atom2):
        """Find and return a bond between *atom1* and *atom2*. Both atoms have to belong to the molecule. If no bond
        between chosen atoms exists, the returned value is ``None``."""
        if atom1.mol != self or atom2.mol != self:
            raise MoleculeError("find_bond: atoms passed as arguments have to belong to the molecule")
        for b in atom1.bonds:
            if atom2 is b.other_end(atom1):
                return b
        return None

    def set_atoms_id(self, start=1):
        """Equip each atom of the molecule with the ``id`` attribute equal to its position within ``atoms`` list.

        The starting value of the numbering can be set with *start* (starts at 1 by default).
        """
        for i, at in enumerate(self.atoms, start):
            at.id = i

    def unset_atoms_id(self):
        """Delete ``id`` attributes of all atoms."""
        for at in self.atoms:
            try:
                del at.id
            except AttributeError:
                pass

    def neighbors(self, atom):
        """Return a list of neighbors of *atom* within the molecule.

        *atom* has to belong to the molecule. Returned list follows the same order as the ``bonds`` attribute of *atom*.
        """
        if atom.mol != self:
            raise MoleculeError("neighbors: passed atom should belong to the molecule")
        return [b.other_end(atom) for b in atom.bonds]

    def bond_matrix(self):
        """Return a square numpy array with bond orders. The size of the array is equal to the number of atoms."""
        ret = np.zeros((len(self), len(self)))
        self.set_atoms_id(start=0)
        for b in self.bonds:
            i, j = b.atom1.id, b.atom2.id
            ret[i, j] = ret[j, i] = b.order
        self.unset_atoms_id()
        return ret

    def separate(self):
        """Separate the molecule into connected components.

        Returned is a list of new |Molecule| objects (all atoms and bonds are disjoint with the original molecule). Each element of this list is identical to one connected component of the base molecule. A connected component is a subset of atoms such that there exists a path (along one or more bonds) between any two atoms. Usually these connected components are molecules.

        Example::

            >>> mol = Molecule('xyz_dimers/NH3-H2O.xyz')
            >>> mol.guess_bonds()
            >>> print(mol)
              Atoms:
                1         N     -1.395591     -0.021564      0.000037
                2         H     -1.629811      0.961096     -0.106224
                3         H     -1.862767     -0.512544     -0.755974
                4         H     -1.833547     -0.330770      0.862307
                5         O      1.568501      0.105892      0.000005
                6         H      0.606736     -0.033962     -0.000628
                7         H      1.940519     -0.780005      0.000222
              Bonds:
               (5)--1.0--(7)
               (5)--1.0--(6)
               (1)--1.0--(3)
               (1)--1.0--(4)
               (1)--1.0--(2)
            >>> x = mol.separate()
            >>> for i in x: print(i)
              Atoms:
                1         N     -1.395591     -0.021564      0.000037
                2         H     -1.629811      0.961096     -0.106224
                3         H     -1.862767     -0.512544     -0.755974
                4         H     -1.833547     -0.330770      0.862307
              Bonds:
               (1)--1.0--(3)
               (1)--1.0--(4)
               (1)--1.0--(2)

              Atoms:
                1         O      1.568501      0.105892      0.000005
                2         H      0.606736     -0.033962     -0.000628
                3         H      1.940519     -0.780005      0.000222
              Bonds:
               (1)--1.0--(3)
               (1)--1.0--(2)

        """
        frags = []
        clone = self.copy()
        for at in clone:
            at._visited = False

        def dfs(start_v, mol):
            stack = [start_v]

            while stack:
                v = stack.pop()
                if not v._visited:
                    v._visited = True
                    v.mol = mol
                    for e in v.bonds:
                        e.mol = mol
                        u = e.other_end(v)
                        if not u._visited:
                            stack.append(u)

        for src in clone.atoms:
            if not src._visited:
                m = Molecule()
                dfs(src, m)
                frags.append(m)
                frags[-1].lattice = self.lattice

        for at in clone.atoms:
            del at._visited
            at.mol.atoms.append(at)
        for b in clone.bonds:
            b.mol.bonds.append(b)

        return frags

    def guess_bonds(self, atom_subset=None, dmax=1.28, metal_atoms: bool = True):
        """Try to guess bonds in the molecule based on types and positions of atoms.

        All previously existing bonds are removed. New bonds are generated based on interatomic distances and information about maximal number of bonds for each atom type (``connectors`` property, taken from |PeriodicTable|).

        The problem of finding molecular bonds for a given set of atoms in space does not have a general solution, especially considering the fact the chemical bond in itself is not a precisely defined concept. For every method, no matter how sophisticated, there will always be corner cases for which the method produces disputable results. Moreover, depending on the context (area of application) the desired solution for a particular geometry may vary. Please do not treat this method as an oracle always providing a proper solution. The algorithm used here gives very good results for geometries that are not very far from the optimal geometry, especially consisting of lighter atoms. All kinds of organic molecules, including aromatic ones, usually work very well. Problematic results can emerge for transition metal complexes, transition states, incomplete molecules etc.

        The algorithm used scales as *n log n* where *n* is the number of atoms.

        The *atom_subset* argument can be used to limit the bond guessing to a subset of atoms, it should be an iterable container with atoms belonging to this molecule.

        The *dmax* argument gives the maximum value for ratio of the bond length to the sum of atomic radii for the two atoms in the bond.

        metal_atoms : bool
            If True, bonds to metal atoms will be guessed. They are often useful for visualization.  The bond order for any bond to a metal atom will be set to 1.

        .. warning::

            This method works reliably only for geometries representing complete molecules. If some atoms are missing (for example, a protein without hydrogens) the resulting set of bonds would usually contain more bonds or bonds with higher order than expected.

        """

        class HeapElement:
            def __init__(self, order, ratio, atom1, atom2):
                eff_ord = order
                if order == 1.5:  # effective order for aromatic bonds
                    eff_ord = 1.15
                elif order == 1 and {atom1.symbol, atom2.symbol} == {"C", "N"}:
                    eff_ord = 1.11  # effective order for single C-N bond
                value = (eff_ord + 0.9) * ratio
                self.data = (value, order, ratio)
                self.atoms = (atom1, atom2)

            def unpack(self):
                val, o, r = self.data
                at1, at2 = self.atoms
                return val, o, r, at1, at2

            def __lt__(self, other):
                return self.data < other.data

            def __le__(self, other):
                return self.data <= other.data

            def __eq__(self, other):
                return self.data == other.data

            def __ne__(self, other):
                return self.data != other.data

            def __gt__(self, other):
                return self.data > other.data

            def __ge__(self, other):
                return self.data >= other.data

        def get_neighbors(atom_list, dmax):
            """adds attributes ._id, .free, and .cube to all atoms in atom_list"""
            cubesize = dmax * 2.1 * max([at.radius for at in atom_list])

            cubes = {}
            for i, at in enumerate(atom_list, 1):
                at._id = i
                at.free = at.connectors
                at.cube = tuple(map(lambda x: int(math.floor(x / cubesize)), at.coords))
                if at.cube in cubes:
                    cubes[at.cube].append(at)
                else:
                    cubes[at.cube] = [at]

            neighbors = {}
            for cube in cubes:
                neighbors[cube] = []
                for i in range(cube[0] - 1, cube[0] + 2):
                    for j in range(cube[1] - 1, cube[1] + 2):
                        for k in range(cube[2] - 1, cube[2] + 2):
                            if (i, j, k) in cubes:
                                neighbors[cube] += cubes[(i, j, k)]

            return neighbors

        def find_and_add_bonds(
            atom_list, neighbors, dmax, from_atoms_subset=None, to_atoms_subset=None, ignore_free=False
        ):
            if from_atoms_subset is None:
                from_atoms_subset = atom_list
            elif not all([x in atom_list for x in from_atoms_subset]):
                raise ValueError("from_atoms_subset must be a subset of atoms_subset")
            if to_atoms_subset is None:
                to_atoms_subset = atom_list
            elif not all([x in atom_list for x in to_atoms_subset]):
                raise ValueError("to_atoms_subset must be a subset of atoms_subset")

            heap = []
            for at1 in from_atoms_subset:
                if at1.free > 0 or ignore_free:
                    for at2 in neighbors[at1.cube]:
                        if not at2 in to_atoms_subset:
                            continue
                        if ignore_free:
                            if at2 in from_atoms_subset:
                                if at2._id <= at1._id:
                                    continue
                        else:
                            if at2.free <= 0 or at2._id <= at1._id:
                                continue
                        # the bond guessing is more accurate with smaller metallic radii
                        ratio = at1.distance_to(at2) / (
                            at1.radius * (1 - 0.1 * at1.is_metallic) + at2.radius * (1 - 0.1 * at2.is_metallic)
                        )
                        if ratio < dmax:
                            if ignore_free:
                                self.add_bond(at1, at2, 1)
                            else:
                                heap.append(HeapElement(0, ratio, at1, at2))
                                # I hate to do this, but I guess there's no other way :/ [MiHa]
                                if at1.atnum == 16 and at2.atnum == 8:
                                    at1.free = 6
                                elif at2.atnum == 16 and at1.atnum == 8:
                                    at2.free = 6
                                elif at1.atnum == 7:
                                    at1.free += 1
                                elif at2.atnum == 7:
                                    at2.free += 1
            if not ignore_free:
                heapq.heapify(heap)

                for at in atom_list:
                    if at.atnum == 7:
                        if at.free > 6:
                            at.free = 4
                        else:
                            at.free = 3

                while heap:
                    val, o, r, at1, at2 = heapq.heappop(heap).unpack()
                    step = 1 if o in [0, 2] else 0.5
                    if at1.free >= step and at2.free >= step:
                        o += step
                        at1.free -= step
                        at2.free -= step
                        if o < 3:
                            heapq.heappush(heap, HeapElement(o, r, at1, at2))
                        else:
                            self.add_bond(at1, at2, o)
                    elif o > 0:
                        if o == 1.5:
                            o = Bond.AR
                        self.add_bond(at1, at2, o)

                def dfs(atom, par):
                    atom.arom += 1000
                    for b in atom.bonds:
                        oe = b.other_end(atom)
                        if b.is_aromatic() and oe.arom < 1000:
                            if oe.arom > 2:
                                return False
                            if par and oe.arom == 1:
                                b.order = 2
                                return True
                            if dfs(oe, 1 - par):
                                b.order = 1 + par
                                return True

                for at in atom_list:
                    at.arom = len(list(filter(Bond.is_aromatic, at.bonds)))

                for at in atom_list:
                    if at.arom == 1:
                        dfs(at, 1)

        def cleanup_atom_list(atom_list):
            for at in atom_list:
                del at.cube, at.free, at._id
                if hasattr(at, "arom"):
                    del at.arom
                if hasattr(at, "_metalbondcounter"):
                    del at._metalbondcounter
                if hasattr(at, "_electronegativebondcounter"):
                    del at._electronegativebondcounter

        self.delete_all_bonds()
        atom_list = atom_subset or self.atoms

        neighbors = get_neighbors(atom_list, dmax)

        nonmetallic = [x for x in atom_list if not x.is_metallic]
        metallic = [x for x in atom_list if x.is_metallic]
        hydrogens = [x for x in atom_list if x.atnum == 1]
        potentially_ignore_metal_bonds = [x for x in atom_list if x.symbol in ["C", "N", "S", "P", "As"]]

        # first guess bonds for non-metals. This also captures bond orders.
        find_and_add_bonds(nonmetallic, neighbors, dmax=dmax)

        # add stray hydrogens
        stray_hydrogens = [x for x in hydrogens if len(x.bonds) == 0]
        find_and_add_bonds(
            nonmetallic,
            neighbors,
            from_atoms_subset=stray_hydrogens,
            to_atoms_subset=nonmetallic,
            ignore_free=True,
            dmax=dmax,
        )

        if metal_atoms:
            # bonds to metal atoms are very useful for visualization but
            # are not typically used in force fields. So provide an option to not add them

            # for obvious anions like carbonate, nitrate, sulfate, phosphate, and arsenate, do not allow metal atoms to bond to the central atom
            new_atom_list = []
            for at in atom_list:
                if at in potentially_ignore_metal_bonds:
                    if len([x for x in at.bonds if x.other_end(at).is_electronegative]) >= 3:
                        continue
                new_atom_list.append(at)
            find_and_add_bonds(
                atom_list,
                neighbors,
                from_atoms_subset=metallic,
                to_atoms_subset=new_atom_list,
                ignore_free=True,
                dmax=dmax,
            )

            # delete metal-metal bonds and metal-hydrogen bonds if the metal is bonded to enough electronegative atoms and not enough metal atoms
            # (this means that the metal is a cation, so bonds should almost never be drawn unless it's a dimetal complex or a hydride/H2 ligand, but that should be rare)
            for at in metallic:
                at._metalbondcounter = len([x for x in at.bonds if x.other_end(at).is_metallic])
                at._electronegativebondcounter = len([x for x in at.bonds if x.other_end(at).is_electronegative])
                if (
                    at._electronegativebondcounter >= 3
                    or (at._electronegativebondcounter >= 2 >= at._metalbondcounter)
                    or (at._electronegativebondcounter >= 1 and at._metalbondcounter <= 0)
                ):
                    bonds_to_delete = [b for b in at.bonds if b.other_end(at).is_metallic or b.other_end(at).atnum == 1]
                    for b in bonds_to_delete:
                        self.delete_bond(b)

        cleanup_atom_list(atom_list)

    def guess_system_charge(self):
        """
        Attempt to guess the charge of the full system based on connectivity
        """
        charge = sum(self.guess_atomic_charges(adjust_to_systemcharge=False))
        return charge

    def guess_atomic_charges(
        self, adjust_to_systemcharge=True, keep_hydrogen_charged=False, depth=1, electronegativities=None
    ):
        """
        Return a list of guessed charges, one for each atom, based on connectivity

        * ``depth`` -- The electronegativity of an atom is determined all its neighbors up to depth
        * ``electronegativities`` -- A dictionary containing electronegativity values for the electronegative elements

        Note: Fairly basic implementation that will not always yield reliable results
        """

        def get_electronegativity(atom, prevat, search_depth=None):
            """
            Get the electronegativity of atom by searching through the molecules
            """
            ens = get_electronegativities(atom, [self.index(prevat) - 1], search_depth)
            en = sum(ens) / len(ens) if len(ens) > 0 else 0.0
            return en

        def get_electronegativities(atom, prevats=[], search_depth=None):
            """
            Get the electronegativities of neighbors by searching through the molecules
            """
            en = [electronegativities[atom.symbol] if atom.symbol in electronegativities else None]
            en = [v for v in en if v is not None]
            if search_depth is not None:
                search_depth -= 1
                if search_depth <= 0:
                    return en
            neighbors = [at for at in self.neighbors(atom) if not self.index(at) - 1 in prevats]
            prevats = prevats + [self.index(atom) - 1]
            for other_at in neighbors:
                en += get_electronegativities(other_at, prevats, search_depth)
                prevats += [self.index(other_at) - 1]
            return en

        if electronegativities is None:
            # https://pubchem.ncbi.nlm.nih.gov/periodic-table/electronegativity/
            electronegativities = {
                "Te": 2.1,
                "P": 2.19,
                "At": 2.2,
                "C": 2.55,
                "Se": 2.55,
                "S": 2.58,
                "I": 2.66,
                "Br": 2.96,
                "N": 3.04,
                "Cl": 3.16,
                "O": 3.44,
                "F": 3.98,
            }
        charges = [0.0 for atom in self.atoms]

        # Negative charge to unsaturated electronegative atoms (C, N, O, P, S, etc)
        # I am giving them all as negative a charge as their number of max connectors, then compensated by the neighbors.
        # The exception is if the neighbor is also electro-negative. Then they cancel each other out.
        # electronegative = [PT.get_symbol(i) for i,values in enumerate(PT.data) if PT.get_electronegative(i) and (not PT.get_symbol(i)=='C')]
        electronegative = [PT.get_symbol(i) for i, values in enumerate(PT.data) if PT.get_electronegative(i)]
        # This will go very wrong for a lot of elements like P, Si, Al,
        # which appear to have a connectors value of 8
        electronegative = [s for s in electronegative if not PT.get_connectors(PT.get_atomic_number(s)) == 8]
        # First select the electronegative indices
        en_indices = []
        for i, atom in enumerate(self.atoms):
            if not atom.symbol in electronegative:
                continue
            # MiHa made a similar adjustment in PLAMS for the valence of S bonded to an O (which is 6 instead of 2).
            if "O" in [at.symbol for at in atom.neighbors()] and atom.symbol == "S":
                continue
            en_indices.append(i)
        # Then assign the charges to these atoms and their neighbors
        echarges = [0.0 for at in self.atoms]
        for i in en_indices:
            atom = self.atoms[i]
            echarges[i] -= atom.connectors
            for bond in atom.bonds:
                other_at = bond.other_end(atom)
                j = self.index(other_at) - 1
                # Here we make sure that C is only treated as electronegative towards non-electronegative neighbors
                # If there was any relative electronegativity data, this could be generalized
                if bond.order % 1 > 0.4 and bond.order % 1 < 0.6:
                    order = round(bond.order * 2) / 2  # Rounded to 0.5 only if it is very obviously a partial bond
                else:
                    order = round(bond.order)
                # Compare the electronegativity values to decide where the electrons go
                en_i = get_electronegativity(atom, other_at, search_depth=depth)
                en_j = get_electronegativity(other_at, atom, search_depth=depth)
                # en_i = electronegativities[atom.symbol]
                # en_j = electronegativities[other_at.symbol] if other_at.symbol in electronegativities else 0.
                if en_i <= en_j and j in en_indices:
                    # if atom.symbol == 'C' and j in en_indices :
                    echarges[i] += order
                else:
                    echarges[j] += order
        charges = [q + echarges[i] for i, q in enumerate(charges)]

        # Assign positive charge to H and Compensate charges to the rest (assumes that H-atoms have only single neighbor)
        q_hydrogens = [1.0 if atom.symbol == "H" and charges[i] < 1.0 else 0.0 for i, atom in enumerate(self.atoms)]
        for i, atom in enumerate(self.atoms):
            if atom.symbol == "H":
                continue
            if i in en_indices:
                continue
            neighbors = atom.neighbors()
            nhs = len([n for n in neighbors if n.symbol == "H"])
            q = -nhs
            q_hydrogens[i] += q
        charges = [q + q_hydrogens[i] for i, q in enumerate(charges)]

        # Formal charges are generally not on the H-atoms, so we will displace them all to their neighbors
        if not keep_hydrogen_charged:
            hydrogens = [i for i, at in enumerate(self.atoms) if at.symbol == "H"]
            for ih in hydrogens:
                q = charges[ih]
                neighbors = [at for at in self.neighbors(self.atoms[ih]) if not at.symbol == "H"]
                for at in neighbors:
                    charges[self.index(at) - 1] += q / (len(neighbors))
                if len(neighbors) > 0:
                    charges[ih] -= q

        # Check the assigned charges, and possibly adjust
        if adjust_to_systemcharge:
            if "charge" in self.properties.keys():
                molcharge = float(self.properties.charge)
                if sum(charges) != molcharge:
                    # If there is a problem with the charges, just take the highest charged atom, and change that
                    log(
                        "Warning: Guessed atomic charges (%i) using PLAMS do not match assigned charge molecule (%i)"
                        % (sum(charges), molcharge)
                    )
                    dq = sum(charges) - molcharge
                    if abs(dq) > 1:
                        log("Charges: %s" % (" ".join([str(q) for q in charges])))
                        raise MoleculeError(
                            "Guessed atomic charges (%i) using PLAMS do not match assigned charge molecule (%i)"
                            % (sum(charges), molcharge)
                        )
                    tmpcharges = charges.copy()
                    if dq < 0:
                        tmpcharges = [-q for q in charges]
                    ind = tmpcharges.index(max(tmpcharges))
                    tmpcharges[ind] -= abs(dq)
                    if dq < 0:
                        tmpcharges = [-q for q in tmpcharges]
                    charges = tmpcharges

        return charges

    def in_ring(self, arg):
        """Check if an atom or a bond belonging to this |Molecule| forms a ring. *arg* should be an instance of |Atom| or |Bond| belonging to this |Molecule|."""

        if (not isinstance(arg, (Atom, Bond))) or arg.mol != self:
            raise MoleculeError("in_ring: Argument should be a Bond or an Atom and it should be a part of the Molecule")

        def dfs(v, depth=0):
            v._visited = True
            for bond in v.bonds:
                if bond is not arg:
                    u = bond.other_end(v)
                    if u is arg and depth > 1:
                        u._visited = "cycle"
                    if not u._visited:
                        dfs(u, depth + 1)

        for at in self:
            at._visited = False

        if isinstance(arg, Atom):
            dfs(arg)
            ret = arg._visited == "cycle"
        else:
            dfs(arg.atom1)
            ret = arg.atom2._visited

        for at in self:
            del at._visited
        return ret

    def supercell(self, *args) -> "Molecule":
        """Return a new |Molecule| instance representing a supercell build by replicating this |Molecule| along its lattice vectors.

        One should provide in input an integer matrix :math:`T_{i,j}` representing the supercell transformation (:math:`\\vec{a}_i' = \sum_j T_{i,j}\\vec{a}_j`). The size of the matrix should match the number of lattice vectors, i.e. 3x3 for 3D periodic systems, 2x2 for 2D periodic systems and one number for 1D periodic systems. The matrix can be provided in input as either a nested list or as a numpy matrix.

        For a diagonal supercell expansion (i.e. :math:`T_{i \\neq j}=0`) one can provide in input n positive integers instead of a matrix, where n is number of lattice vectors in the molecule. e.g. This ``mol.supercell([[2,0],[0,2]])`` is equivalent to ``mol.supercell(2,2)``.

        The returned |Molecule| is fully distinct from the current one, in a sense that it contains a different set of |Atom| and |Bond| instances. However, each atom of the returned |Molecule| carries an additional information about its origin within the supercell. If ``atom`` is an |Atom| instance in the supercell, ``atom.properties.supercell.origin`` points to the |Atom| instance of the original molecule that was copied to create ``atom``, while ``atom.properties.supercell.index`` stores the tuple (with length equal to the number of lattice vectors) with cell index. For example, ``atom.properties.supercell.index == (2,1,0)`` means that ``atom`` is a copy of ``atom.properties.supercell.origin`` that was translated twice along the first lattice vector, once along the second vector, and not translated along the third vector.

        Example usage:

        .. code-block:: python

            >>> graphene = Molecule('graphene.xyz')
            >>> print(graphene)
              Atoms:
                1         C      0.000000      0.000000      0.000000
                2         C      1.230000      0.710000      0.000000
              Lattice:
                    2.4600000000     0.0000000000     0.0000000000
                    1.2300000000     2.1304224933     0.0000000000

            >>> graphene_supercell = graphene.supercell(2,2) # diagonal supercell expansion
            >>> print(graphene_supercell)
              Atoms:
                1         C      0.000000      0.000000      0.000000
                2         C      1.230000      0.710000      0.000000
                3         C      1.230000      2.130422      0.000000
                4         C      2.460000      2.840422      0.000000
                5         C      2.460000      0.000000      0.000000
                6         C      3.690000      0.710000      0.000000
                7         C      3.690000      2.130422      0.000000
                8         C      4.920000      2.840422      0.000000
              Lattice:
                    4.9200000000     0.0000000000     0.0000000000
                    2.4600000000     4.2608449866     0.0000000000

            >>> diamond = Molecule('diamond.xyz')
            >>> print(diamond)
              Atoms:
                1         C     -0.446100     -0.446200     -0.446300
                2         C      0.446400      0.446500      0.446600
              Lattice:
                    0.0000000000     1.7850000000     1.7850000000
                    1.7850000000     0.0000000000     1.7850000000
                    1.7850000000     1.7850000000     0.0000000000

            >>> diamond_supercell = diamond.supercell([[-1,1,1],[1,-1,1],[1,1,-1]])
            >>> print(diamond_supercell)
              Atoms:
                1         C     -0.446100     -0.446200     -0.446300
                2         C      0.446400      0.446500      0.446600
                3         C      1.338900      1.338800     -0.446300
                4         C      2.231400      2.231500      0.446600
                5         C      1.338900     -0.446200      1.338700
                6         C      2.231400      0.446500      2.231600
                7         C     -0.446100      1.338800      1.338700
                8         C      0.446400      2.231500      2.231600
              Lattice:
                    3.5700000000     0.0000000000     0.0000000000
                    0.0000000000     3.5700000000     0.0000000000
                    0.0000000000     0.0000000000     3.5700000000
        """

        def diagonal_supercell(*args):
            supercell_lattice = [tuple(n * np.array(vec)) for n, vec in zip(args, self.lattice)]
            cell_translations = [t for t in itertools.product(*[range(arg) for arg in args])]
            return supercell_lattice, cell_translations

        def general_supercell(S):
            determinant = int(round(np.linalg.det(S)))
            if determinant < 1:
                raise MoleculeError(
                    f"supercell: The determinant of the supercell transformation should be one or larger. Determinant: {determinant}."
                )

            supercell_lattice = [tuple(vec) for vec in S @ np.array(self.lattice)]

            max_supercell_index = np.max(abs(S))
            all_possible_translations = itertools.product(
                range(-max_supercell_index, max_supercell_index + 1), repeat=len(self.lattice)
            )
            S_inv = np.linalg.inv(S)

            tol = 1e-10
            cell_translations = []
            for index in all_possible_translations:
                fractional_coord = np.dot(index, S_inv)
                if all(fractional_coord > -tol) and all(fractional_coord < 1.0 - tol):
                    cell_translations.append(index)

            if len(cell_translations) != determinant:
                raise MoleculeError(
                    f"supercell: Failed to find the appropriate supercell translations. We expected to find {determinant} cells, but we found {len(cell_translations)}"
                )

            return supercell_lattice, cell_translations

        if len(args) == 0:
            raise MoleculeError("supercell: This function needs input arguments...")

        if all(isinstance(arg, int) for arg in args):
            # diagonal supercell expansion
            if len(args) != len(self.lattice):
                raise MoleculeError(
                    "supercell: The lattice has {} vectors, but {} arguments were given".format(
                        len(self.lattice), len(args)
                    )
                )
            supercell_lattice, cell_translations = diagonal_supercell(*args)

        elif len(args) == 1 and hasattr(args[0], "__len__"):
            # general_supercell
            try:
                S = np.array(args[0], dtype=int)
                assert S.shape == (len(self.lattice), len(self.lattice))
            except:
                n = len(self.lattice)
                raise MoleculeError(
                    f"supercell: For {n}D system the supercell method expects a {n}x{n} integer matrix (provided as a nested list or as numpy array) or {n} integers."
                )

            supercell_lattice, cell_translations = general_supercell(S)

        else:
            raise MoleculeError(f"supercell: invalid input {args}.")

        tmp = self.copy()
        for parent, son in zip(self, tmp):
            son.properties.supercell.origin = parent

        ret = Molecule()
        for index in cell_translations:
            newmol = tmp.copy()
            for atom in newmol:
                atom.properties.supercell.index = index
            newmol.translate(sum(i * np.array(vec) for i, vec in zip(index, self.lattice)))
            ret += newmol

        ret.lattice = supercell_lattice
        return ret

    def unit_cell_volume(self, unit="angstrom"):
        """Return the volume of the unit cell of a 3D system.

        *unit* is the unit of length, the cube of which will be used as the unit of volume.
        """
        if len(self.lattice) == 3:
            return (
                float(np.linalg.det(np.dstack([self.lattice[0], self.lattice[1], self.lattice[2]])))
                * Units.conversion_ratio("angstrom", unit) ** 3
            )
        elif len(self.lattice) == 2:
            return (
                float(np.linalg.norm(np.cross(self.lattice[0], self.lattice[1])))
                * Units.conversion_ratio("angstrom", unit) ** 2
            )
        elif len(self.lattice) == 1:
            return float(np.linalg.norm(self.lattice[0])) * Units.conversion_ratio("angstrom", unit)
        elif len(self.lattice) == 0:
            raise ValueError("Cannot calculate unit cell volume for a non-periodic molecule")
        else:
            raise ValueError("len(self.lattice) = {}, should be <=3.".format(len(self.lattice)))

    def cell_lengths(self, unit="angstrom"):
        """Return the lengths of the lattice vector. Returns a list with the same length as self.lattice"""

        return cell_lengths(self.lattice, unit=unit)

    def cell_angles(self, unit="degree"):
        """Return the angles between lattice vectors.

        unit : str
            output unit

        For 2D systems, returns a list [gamma]

        For 3D systems, returns a list [alpha, beta, gamma]
        """
        return cell_angles(self.lattice, unit=unit)

    def set_integer_bonds(self, action="warn", tolerance=10**-4):
        """Convert non-integer bond orders into integers.

        For example, bond orders of aromatic systems are no longer set to the non-integer
        value of ``1.5``, instead adopting bond orders of ``1`` and ``2``.

        The implemented function walks a set of graphs constructed from all non-integer bonds,
        converting the orders of aforementioned bonds to integers by alternating calls to
        :func:`math.ceil` and :func:`math.floor`.
        The implication herein is that both :math:`i` and :math:`i+1` are considered valid
        (integer) values for any bond order within the :math:`(i, i+1)` interval.
        Floats which can be represented exactly as an integer, *e.g.* :math:`1.0`,
        are herein treated as integers.

        Can be used for sanitizing any Molecules passed to the
        :mod:`rdkit<scm.plams.interfaces.molecule.rdkit>` module,
        as its functions are generally unable to handle Molecules with non-integer bond orders.

        By default this function will issue a warning if the total (summed) bond orders
        before and after are not equal to each other within a given *tolerance*.
        Accepted values are for *action* are ``"ignore"``, ``"warn"`` and ``"raise"``,
        which respectively ignore such cases, issue a warning or raise a :exc:`MoleculeError`.

        .. code-block:: python

            >>> from scm.plams import Molecule

            >>> benzene = Molecule(...)
            >>> print(benzene)
              Atoms:
                1         C      1.193860     -0.689276      0.000000
                2         C      1.193860      0.689276      0.000000
                3         C      0.000000      1.378551      0.000000
                4         C     -1.193860      0.689276      0.000000
                5         C     -1.193860     -0.689276      0.000000
                6         C     -0.000000     -1.378551      0.000000
                7         H      2.132911     -1.231437     -0.000000
                8         H      2.132911      1.231437     -0.000000
                9         H      0.000000      2.462874     -0.000000
               10         H     -2.132911      1.231437     -0.000000
               11         H     -2.132911     -1.231437     -0.000000
               12         H     -0.000000     -2.462874     -0.000000
              Bonds:
               (3)--1.5--(4)
               (5)--1.5--(6)
               (1)--1.5--(6)
               (2)--1.5--(3)
               (4)--1.5--(5)
               (1)--1.5--(2)
               (3)--1.0--(9)
               (6)--1.0--(12)
               (5)--1.0--(11)
               (4)--1.0--(10)
               (2)--1.0--(8)
               (1)--1.0--(7)

            >>> benzene.set_integer_bonds()
            >>> print(benzene)
              Atoms:
                1         C      1.193860     -0.689276      0.000000
                2         C      1.193860      0.689276      0.000000
                3         C      0.000000      1.378551      0.000000
                4         C     -1.193860      0.689276      0.000000
                5         C     -1.193860     -0.689276      0.000000
                6         C     -0.000000     -1.378551      0.000000
                7         H      2.132911     -1.231437     -0.000000
                8         H      2.132911      1.231437     -0.000000
                9         H      0.000000      2.462874     -0.000000
               10         H     -2.132911      1.231437     -0.000000
               11         H     -2.132911     -1.231437     -0.000000
               12         H     -0.000000     -2.462874     -0.000000
              Bonds:
               (3)--1.0--(4)
               (5)--1.0--(6)
               (1)--2.0--(6)
               (2)--2.0--(3)
               (4)--2.0--(5)
               (1)--1.0--(2)
               (3)--1.0--(9)
               (6)--1.0--(12)
               (5)--1.0--(11)
               (4)--1.0--(10)
               (2)--1.0--(8)
               (1)--1.0--(7)

        """
        # Ignore, raise or warn
        action_func = parse_action(action)

        ceil = math.ceil
        floor = math.floor
        func_invert = {ceil: floor, floor: ceil}

        def dfs(atom, func) -> None:
            """Depth-first search algorithm for integer-ifying the bond orders."""
            for b2 in atom.bonds:
                if b2._visited:
                    continue

                b2._visited = True
                b2.order = func(b2.order)  # func = ``math.ceil()`` or ``math.floor()``
                del bond_dict[b2]

                atom_new = b2.other_end(atom)
                dfs(atom_new, func=func_invert[func])

        def collect_and_mark_bonds(self):
            order_before = []
            order_before_append = order_before.append

            # Mark all non-integer bonds; floats which can be represented exactly
            # by an integer (e.g. 1.0 and 2.0) are herein treated as integers
            bond_dict = OrderedDict()  # An improvised OrderedSet (as it does not exist)
            for bond in self.bonds:
                order = bond.order
                order_before_append(order)
                if (
                    hasattr(bond.order, "is_integer") and not bond.order.is_integer()
                ):  # Checking for ``is_integer()`` catches both float and np.float
                    bond._visited = False
                    bond_dict[bond] = None
                else:
                    bond._visited = True
            return bond_dict, order_before

        bond_dict, order_before = collect_and_mark_bonds(self)

        while bond_dict:
            b1, _ = bond_dict.popitem()
            order = b1.order

            # Start with either ``math.ceil()`` if the ceiling is closer than the floor;
            # start with ``math.floor()`` otherwise
            delta_ceil, delta_floor = ceil(order) - order, floor(order) - order
            func = ceil if abs(delta_ceil) < abs(delta_floor) else floor

            b1.order = func(order)
            b1._visited = True
            dfs(b1.atom1, func=func_invert[func])
            dfs(b1.atom2, func=func_invert[func])

        # Remove the Bond._visited attribute
        order_after_sum = 0.0
        for bond in self.bonds:
            order_after_sum += bond.order
            del bond._visited

        # Check that the total (summed) bond order has not changed
        order_before_sum = sum(order_before)
        if abs(order_before_sum - order_after_sum) > tolerance:
            err = MoleculeError(
                f"Bond orders before and after not equal to tolerance {tolerance!r}:\n"
                f"before: sum(...) == {order_before_sum!r}\n"
                f"after: sum(...) == {order_after_sum!r}"
            )
            try:
                action_func(err)
            except MoleculeError as ex:  # Restore the initial bond orders
                for b, order in zip(self.bonds, reversed(order_before)):
                    b.order = order
                raise ex

    def index(self, value, start=1, stop=None):
        """Return the first index of the specified Atom or Bond.

        Providing an |Atom| will return its 1-based index, while a |Bond| returns a 2-tuple with the 1-based indices of its atoms.

        Raises a |MoleculeError| if the provided is not an Atom/Bond or if the Atom/bond is not part of the molecule.

        .. code:: python

            >>> from scm.plams import Molecule, Bond, Atom

            >>> mol = Molecule(...)
            >>> atom: Atom = Molecule[1]
            >>> bond: Bond = Molecule[1, 2]

            >>> print(mol.index(atom))
            1

            >>> print(mol.index(bond))
            (1, 2)

        """
        args = [start - 1 if start > 0 else start]
        if stop is not None:  # Correct for the 1-based indices used in Molecule
            args.append(stop - 1 if stop > 0 else stop)

        try:
            if isinstance(value, Atom):
                return 1 + self.atoms.index(value, *args)
            elif isinstance(value, Bond):
                return 1 + self.atoms.index(value.atom1, *args), 1 + self.atoms.index(value.atom2, *args)

        except ValueError as ex:  # Raised if the provided Atom/Bond is not in self
            raise MoleculeError(f"Provided {value.__class__.__name__} is not in Molecule").with_traceback(
                ex.__traceback__
            )
        else:  # Raised if value is neither an Atom nor Bond
            raise MoleculeError(f"'value' expected an Atom or Bond; observed type: '{value.__class__.__name__}'")

    def round_coords(self, decimals=0, inplace=True):
        """Round the Cartesian coordinates of this instance to *decimals*.

        By default, with ``inplace=True``, the coordinates of this instance are updated inplace.
        If ``inplace=False`` then a new copy of this Molecule is returned with its
        coordinates rounded.

        .. code:: python

            >>> from scm.plams import Molecule

            >>> mol = Molecule(...)
              Atoms:
                1         H      1.234567      0.000000      0.000000
                2         H      0.000000      0.000000      0.000000

            >>> mol_rounded = round_coords(mol)
            >>> print(mol_rounded)
              Atoms:
                1         H      1.000000      0.000000      0.000000
                2         H      0.000000      0.000000      0.000000

            >>> mol.round_coords(decimals=3)
            >>> print(mol)
              Atoms:
                1         H      1.234000      0.000000      0.000000
                2         H      0.000000      0.000000      0.000000

        """
        xyz = self.as_array()

        # Follow the convention used in ``ndarray.round()``: always return floats,
        # even if ndigits=None
        xyz_round = xyz.round(decimals=decimals)

        if inplace:
            self.from_array(xyz_round)
            return None
        else:
            mol_copy = self.copy()
            mol_copy.from_array(xyz_round)
            return mol_copy

    def get_connection_table(self):
        """
        Get a connection table with atom indices (starting at 0)
        """
        table = []
        for iat, at in enumerate(self.atoms):
            indices = sorted([self.index(n) - 1 for n in self.neighbors(at)])
            table.append(indices)
        return table

    def get_molecule_indices(self):
        """
        Use the bond information to identify submolecules

        Returns a list of lists of indices (e.g. for two methane molecules: [[0,1,2,3,4],[5,6,7,8,9]])
        """
        molecule_indices = []
        for at in self:
            at._visited = False

        def dfs(v, indices):
            """
            Depth first search of self starting at atom v, extending the list of connected atoms (indices)
            """
            v._visited = True
            for e in v.bonds:
                u = e.other_end(v)
                if not u._visited:
                    indices.append(self.index(u, start=indices[0] + 1) - 1)
                    dfs(u, indices)

        def bfs(v, indices):
            """
            Breadth first search of self starting at atom v, extending the list of connected atoms (indices)
            """
            # REB: Changed to bfs to avoid Pythons recursion error for big systems
            atoms = [v]
            while len(atoms) > 0:
                for v in atoms:
                    v._visited = True
                res = []
                [
                    res.append(e.other_end(v))
                    for v in atoms
                    for e in v.bonds
                    if not e.other_end(v)._visited and not e.other_end(v) in res
                ]
                atoms = res
                # atoms = [e.other_end(v) for v in atoms for e in v.bonds if not e.other_end(v)._visited]
                for u in atoms:
                    indices.append(self.index(u, start=indices[0] + 1) - 1)

        for iatom, src in enumerate(self.atoms):
            if not src._visited:
                indices = [iatom]
                # dfs(src, indices)
                bfs(src, indices)
                molecule_indices.append(sorted(indices))

        for at in self:
            del at._visited

        return molecule_indices

    def get_fragment(self, indices):
        """
        Return a submolecule from self
        """
        ret = Molecule()
        ret.lattice = self.lattice.copy()

        # First the atoms
        bro = {}
        for iat in indices:
            at = self.atoms[iat]
            at_copy = smart_copy(at, owncopy=["properties"], without=["mol", "bonds"])
            ret.add_atom(at_copy)
            bro[at] = at_copy

        # Then the bonds
        for bo in self.bonds:
            if (bo.atom1 in bro) and (bo.atom2 in bro):
                bo_copy = smart_copy(bo, owncopy=["properties"], without=["atom1", "atom2", "mol"])
                bo_copy.atom1 = bro[bo.atom1]
                bo_copy.atom2 = bro[bo.atom2]
                ret.add_bond(bo_copy)

        ret.properties.soft_update(self.properties)

        return ret

    def get_complete_molecules_within_threshold(self, atom_indices, threshold: float):
        """
        Returns a new molecule containing complete submolecules for any molecules
        that are closer than ``threshold`` to any of the atoms in ``atom_indices``.

        Note: This only works for non-periodic systems.

        atom_indices: list of int
            One-based indices of the atoms

        threshold : float
            Distance threshold for whether to include molecules
        """
        if self.lattice:
            raise ValueError("Cannot run get_complete_molecules_within_threshold() on a Molecule with a lattice")
        molecule_indices = self.get_molecule_indices()  # [[0,1,2],[3,4],[5,6],...]

        solvated_coords = self.as_array()
        zero_based_indices = [x - 1 for x in atom_indices]
        D = distance_array(solvated_coords, solvated_coords)[zero_based_indices]
        less_equal = np.less_equal(D, threshold)
        within_threshold = np.any(less_equal, axis=0)
        good_indices = [i for i, value in enumerate(within_threshold) if value]

        complete_indices = set()
        for indlist in molecule_indices:
            for ind in good_indices:
                if ind in indlist:
                    complete_indices = complete_indices.union(indlist)
                    break
        complete_indices = sorted(list(complete_indices))

        newmolecule = self.get_fragment([i for i in complete_indices])
        return newmolecule

    def locate_rings(self):
        """
        Find the rings in the structure
        """

        def bond_from_indices(ret, iat1, iat2):
            """
            Return a bond object from the atom indices
            """
            bond = None
            for bond in ret.bonds:
                indices = [i - 1 for i in ret.index(bond)]
                if iat1 in indices and iat2 in indices:
                    break
            if bond is None:
                raise Exception("This bond should exist (%i %i), but does not!" % (iat1, iat2))
            return bond

        # For each atom find the list of neighbors,
        # break the corresponding bond, and then find out if they are still
        # connected. If so, find the shortest bridge, using Dijkstra's algorithm
        conect = self.get_connection_table()
        atoms = [i for i, at in enumerate(self)]
        ret = self.copy()
        allrings = []
        oldneighbors = []
        for iat in atoms:
            neighbors = conect[iat]
            for iatn in neighbors:
                if iatn in oldneighbors:
                    continue
                bond = bond_from_indices(ret, iat, iatn)
                # Delete the bond
                ret.delete_bond(bond)
                mollist = ret.get_molecule_indices()
                for m in mollist:
                    connection = False
                    if iatn in m:
                        if iat in m:
                            connection = True
                        break

                if connection:
                    retconect = ret.get_connection_table()
                    rings = self.shortest_path_dijkstra(iat, iatn, conect=retconect)
                    for ring in rings:
                        ring.sort()
                        if not ring in allrings:
                            allrings.append(ring)
                    # Put the bond back
                    ret.add_bond(bond)
            oldneighbors.append(iat)

        allrings = [self.order_ring(ring) for ring in allrings]
        return allrings

    def locate_rings_acm(self):
        """
        Use the ACM algorithm to find rings
        """

        def find_cycle(tree, root, leaf):
            """
            Find the path from leaf to root in tree

            Note: The problem is that a bond is not added when a cycle has closed.
                  And even if we add that bond, I don't know how to find the path other than with Dijkstra
            """
            # First move from path to the top of the tree
            prev = leaf
            pathb = [leaf]
            while prev != root:
                prevlist = [b[0] for b in tree if b[1] == prev]
                if len(prevlist) == 0:
                    break
                prev = prevlist[0]
                pathb.append(prev)
            if prev == root:
                return pathb
            # Also move forward, from the top to the root
            prev = root
            pathf = [root]
            while prev not in pathb:
                prevlist = [b[0] for b in tree if b[1] == prev]
                if len(prevlist) == 0:
                    break
                prev = prevlist[0]
                pathf.append(prev)
            path = pathb[: pathb.index(prev)] + pathf[::-1]
            return path

        rings = []
        rings_sorted = []

        conect = self.get_connection_table()
        atoms_to_examine = [self.index(at) - 1 for at in self.atoms]
        iat = atoms_to_examine[0]
        atoms_in_tree = [iat]
        tree = [[] for i in range(len(self))]
        counter = 0
        while 1:
            # print ('%8i %8i %8i'%(counter,iat,len(rings)))
            counter += 1
            for jat in conect[iat]:
                if jat in atoms_in_tree:
                    # ring = find_cycle(tree,iat,jat)
                    for ring in self.shortest_path_dijkstra(iat, jat, conect=tree):
                        sorted_ring = sorted(ring)
                        if not sorted_ring in rings_sorted:
                            rings.append(ring)
                            rings_sorted.append(sorted_ring)
                else:
                    atoms_in_tree.append(jat)
                tree[iat].append(jat)
                tree[jat].append(iat)
                # tree.append([iat,jat])
                conect[iat] = [j for j in conect[iat] if not j == jat]
                conect[jat] = [j for j in conect[jat] if not j == iat]
            atoms_to_examine = [i for i in atoms_to_examine if not i == iat]
            if len(atoms_to_examine) == 0:
                break
            candidates = [i for i in atoms_in_tree if i in atoms_to_examine]
            iat = candidates[0] if len(candidates) > 0 else atoms_to_examine[0]
        return rings

    def order_ring(self, ring_indices):
        """
        Order the ring indices so that they are sequential along the ring
        """
        conect = self.get_connection_table()
        new_ring = [ring_indices[0]]
        for iat in new_ring:
            neighbors = [jat for jat in conect[iat] if jat in ring_indices]
            neighbors = [jat for jat in neighbors if not jat in new_ring]
            if len(neighbors) == 0:
                break
            else:
                new_ring.append(neighbors[0])
        return new_ring

    def locate_rings_networkx(self):
        """
        Obtain a list of ring indices using RDKit (same as locate_rings, but much faster)
        """
        import networkx

        matrix = self.bond_matrix()
        matrix = matrix.astype(np.int32)
        matrix[matrix > 0] = 1
        graph = networkx.from_numpy_array(matrix)
        rings = networkx.cycle_basis(graph)
        return rings

    def shortest_path_dijkstra(self, source, target, conect=None):
        """
        Find the shortest paths (can be more than 1)
        between a source atom and a target
        atom in a connection table

        * ``source`` -- Index of the source atom
        * ``target`` -- Index of the target atom
        """
        if conect is None:
            conect = self.get_connection_table()

        huge = 100000.0

        dist = {}
        previous = {}
        for v in range(len(conect)):
            dist[v] = huge
            previous[v] = []
        dist[source] = 0

        Q = [source]
        for iat in range(len(conect)):
            if iat != source:
                Q.append(iat)

        while len(Q) > 0:
            # vertex in Q with smallest distance in dist
            u = Q[0]
            if dist[u] == huge:
                return []
            u = Q.pop(0)
            if u == target:
                break

            # Select the neighbors of u, and loop over them
            neighbors = conect[u]
            for v in neighbors:
                if not v in Q:
                    continue
                alt = dist[u] + 1.0
                if alt == dist[v]:
                    previous[v].append(u)
                if alt < dist[v]:
                    previous[v] = [u]
                    dist[v] = alt
                    # Reorder Q
                    for i, vertex in enumerate(Q):
                        if vertex == v:
                            ind = i
                            break
                    Q.pop(ind)
                    for i, vertex in enumerate(Q):
                        if dist[v] < dist[vertex]:
                            ind = i
                            break
                    Q.insert(ind, v)

        bridgelist = [[u]]
        d = dist[u]
        for i in range(int(d)):
            paths = []
            for j, path in enumerate(bridgelist):
                prevats = previous[path[-1]]
                for at in prevats:
                    newpath = path + [at]
                    paths.append(newpath)
            bridgelist = paths

        return bridgelist

    # ===========================================================================
    # ==== Geometry operations ==================================================
    # ===========================================================================

    def translate(self, vector, unit="angstrom"):
        """Move the molecule in space by *vector*, expressed in *unit*.

        *vector* should be an iterable container of length 3 (usually tuple, list or numpy array). *unit* describes unit of values stored in *vector*.
        """
        xyz_array = self.as_array()
        ratio = Units.conversion_ratio(unit, "angstrom")
        xyz_array += np.array(vector) * ratio
        self.from_array(xyz_array)

    def rotate_lattice(self, matrix):
        """Rotate **only** lattice vectors of the molecule with given rotation *matrix*.

        *matrix* should be a container with 9 numerical values. It can be a list (tuple, numpy array etc.) listing matrix elements row-wise, either flat (``[1,2,3,4,5,6,7,8,9]``) or in two-level fashion (``[[1,2,3],[4,5,6],[7,8,9]]``).

        .. note::

            This method does not check if *matrix* is a proper rotation matrix.
        """
        matrix = np.array(matrix).reshape(3, 3)
        self.lattice = [tuple(np.dot(matrix, i)) for i in self.lattice]

    def rotate(self, matrix, lattice=False):
        """Rotate the molecule with given rotation *matrix*. If *lattice* is ``True``, rotate lattice vectors too.

        *matrix* should be a container with 9 numerical values. It can be a list (tuple, numpy array etc.) listing matrix elements row-wise, either flat (``[1,2,3,4,5,6,7,8,9]``) or in two-level fashion (``[[1,2,3],[4,5,6],[7,8,9]]``).

        .. note::

            This method does not check if *matrix* is a proper rotation matrix.
        """
        xyz_array = self.as_array()
        matrix = np.array(matrix).reshape(3, 3)
        xyz_array = xyz_array @ matrix.T
        self.from_array(xyz_array)
        if lattice:
            self.rotate_lattice(matrix)

    def align_lattice(self, convention="AMS", zero=1e-10):
        """Rotate the molecule in such a way that lattice vectors are aligned with the coordinate system.

        This method is meant to be used with periodic systems only. Using it on a |Molecule| instance with an empty ``lattice`` attribute has no effect.

        Possible values of the *convention* argument are:

        *   ``AMS`` (default) -- for 1D systems the lattice vector aligned with X axis. For 2D systems both lattice vectors aligned with XY plane. No constraints for 3D systems
        *   ``reax`` (convention used by `ReaxFF <https://www.scm.com/product/reaxff>`_) -- second lattice vector (if present) aligned with YZ plane. Third vector (if present) aligned with Z axis.

        *zero* argument can be used to specify the numerical tolerance for zero (used to determine if some vector is already aligned with a particular axis or plane).

        The returned boolean value indicates if any rotation happened.
        """
        dim = len(self.lattice)

        if dim == 0:
            log("NOTE: align_lattice called on a Molecule without any lattice", 5)
            return False

        rotated = False
        if convention == "AMS":
            if dim == 1 and (abs(self.lattice[0][1]) > zero or abs(self.lattice[0][2]) > zero):
                mat = rotation_matrix(self.lattice[0], [1.0, 0.0, 0.0])
                self.rotate(mat, lattice=True)
                rotated = True

            if dim == 2 and (abs(self.lattice[0][2]) > zero or abs(self.lattice[1][2]) > zero):
                mat = rotation_matrix(self.lattice[0], [1.0, 0.0, 0.0])
                self.rotate(mat, lattice=True)
                if abs(self.lattice[1][2]) > zero:
                    mat = rotation_matrix([0.0, self.lattice[1][1], self.lattice[1][2]], [0.0, 1.0, 0.0])
                    self.rotate(mat, lattice=True)
                rotated = True

        elif convention == "reax":
            if dim == 3 and (abs(self.lattice[2][0]) > zero or abs(self.lattice[2][1]) > zero):
                mat = rotation_matrix(self.lattice[2], [0.0, 0.0, 1.0])
                self.rotate(mat, lattice=True)
                rotated = True

            if dim >= 2 and abs(self.lattice[1][0]) > zero:
                mat = rotation_matrix([self.lattice[1][0], self.lattice[1][1], 0.0], [0.0, 1.0, 0.0])
                self.rotate(mat, lattice=True)
                rotated = True

        else:
            raise MoleculeError(
                "align_lattice: unknown convention: {}. Possible values are 'AMS' or 'reax'".format(convention)
            )
        return rotated

    def rotate_bond(self, bond, moving_atom, angle, unit="radian"):
        """Rotate part of this molecule containing *moving_atom* along axis defined by *bond* by an *angle* expressed in *unit*.

        *bond* should be chosen in such a way, that it divides the molecule into two parts (using a bond that forms a ring results in a |MoleculeError|). *moving_atom* has to belong to *bond* and is used to pick which part of the molecule is rotated. A positive angle denotes counterclockwise rotation (when looking along the bond, from the stationary part of the molecule).
        """
        if moving_atom not in bond:
            raise MoleculeError("rotate_bond: atom has to belong to the bond")

        atoms_to_rotate = {moving_atom}

        def dfs(v):
            for e in v.bonds:
                if e is not bond:
                    u = e.other_end(v)
                    if u not in atoms_to_rotate:
                        atoms_to_rotate.add(u)
                        dfs(u)

        dfs(moving_atom)

        if len(atoms_to_rotate) == len(self):
            raise MoleculeError("rotate_bond: chosen bond does not divide the molecule")

        other_end = bond.other_end(moving_atom)
        v = np.array(other_end.vector_to(moving_atom))
        rotmat = axis_rotation_matrix(v, angle, unit)
        trans = np.array(other_end.vector_to((0, 0, 0)))

        xyz_array = self.as_array(atom_subset=atoms_to_rotate)
        xyz_array += trans
        xyz_array = xyz_array @ rotmat.T
        xyz_array -= trans

        self.from_array(xyz_array, atom_subset=atoms_to_rotate)

    def resize_bond(self, bond, moving_atom, length, unit="angstrom"):
        """Change the length of *bond* to *length* expressed in *unit* by moving part of the molecule containing *moving_atom*

        *bond* should be chosen in such a way, that it divides the molecule into two parts (using a bond that forms a ring results in a |MoleculeError|). *moving_atom* has to belong to *bond* and is used to pick which part of the molecule is moved.
        """
        if moving_atom not in bond:
            raise MoleculeError("resize_bond: atom has to belong to the bond")

        atoms_to_move = {moving_atom}

        def dfs(v):
            for e in v.bonds:
                if e is not bond:
                    u = e.other_end(v)
                    if u not in atoms_to_move:
                        atoms_to_move.add(u)
                        dfs(u)

        dfs(moving_atom)

        if len(atoms_to_move) == len(self):
            raise MoleculeError("resize_bond: chosen bond does not divide molecule")

        bond_v = np.array(bond.as_vector(start=moving_atom))
        trans_v = (1 - length / bond.length(unit)) * bond_v

        xyz_array = self.as_array(atom_subset=atoms_to_move)
        xyz_array += trans_v
        self.from_array(xyz_array, atom_subset=atoms_to_move)

    def closest_atom(self, point, unit="angstrom") -> Atom:
        """Return the atom of the molecule that is the closest one to some *point* in space.

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*.
        """
        if isinstance(point, Atom):
            point = np.array(point.coords) * Units.conversion_ratio(unit, "angstrom")
        else:
            point = np.array(point) * Units.conversion_ratio(unit, "angstrom")

        xyz_array = self.as_array()
        dist_array = np.linalg.norm(point - xyz_array, axis=1)
        idx = dist_array.argmin()
        return self[idx + 1]

    def distance_to_point(self, point, unit="angstrom", result_unit="angstrom") -> float:
        """Calculate the distance between the molecule and some *point* in space (distance between *point* and :meth:`closest_atom`).

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*. Returned value is expressed in *result_unit*.
        """
        at = self.closest_atom(point, unit)
        return at.distance_to(point, unit, result_unit)

    def distance_to_mol(self, other, result_unit="angstrom", return_atoms=False):
        """Calculate the distance between the molecule and some *other* molecule.

        The distance is measured as the smallest distance between any atom of this molecule and any atom of *other* molecule. Returned distance is expressed in *result_unit*.

        If *return_atoms* is ``False``, only a single number is returned.  If *return_atoms* is ``True``, the method returns a tuple ``(distance, atom1, atom2)`` where ``atom1`` and ``atom2`` are atoms fulfilling the minimal distance, with atom1 belonging to this molecule and atom2 to *other*.
        """
        xyz_array1 = self.as_array()
        xyz_array2 = other.as_array()

        dist_array = distance_array(xyz_array1, xyz_array2)

        res = Units.convert(dist_array.min(), "angstrom", result_unit)
        if return_atoms:
            idx1, idx2 = np.unravel_index(dist_array.argmin(), dist_array.shape)
            atom1 = self[idx1 + 1]
            atom2 = other[idx2 + 1]
            return res, atom1, atom2
        return res

    def wrap(self, length, angle=2 * math.pi, length_unit="angstrom", angle_unit="radian"):
        """wrap(self, length, angle=2*pi, length_unit='angstrom', angle_unit='radian')

        Transform the molecule wrapping its x-axis around z-axis. This method is useful for building nanotubes or molecular wedding rings.

        Atomic coordinates are transformed in the following way:

        *   z coordinates remain untouched
        *   x axis gets wrapped around the circle centered in the origin of new coordinate system. Each segment of x axis of length *length* ends up as an arc of a circle subtended by an angle *angle*. The radius of this circle is R = *length*/*angle*.
        *   part of the plane between the x axis and the line y=R is transformed into the interior of the circle, with line y=R being squashed into a single point - the center of the circle.
        *   part of the plane above line y=R is dropped
        *   part of the plane below x axis is transformed into outside of the circle
        *   transformation is done in such a way that distances along y axis are preserved

        Before:

        .. image:: ../_static/wrap.*

        After:

        .. image:: ../_static/wrap2.*

        """
        length = Units.convert(length, length_unit, "angstrom")
        angle = Units.convert(angle, angle_unit, "radian")

        xy_array = self.as_array().T[0:2]
        if xy_array[0].ptp() > length:
            raise MoleculeError("wrap: x-extension of the molecule is larger than length")

        if angle < 0 or angle > 2 * math.pi:
            raise MoleculeError("wrap: angle must be between 0 and 2*pi")

        R = length / angle

        x_array = (R - xy_array[1]) * np.cos(xy_array[0] / R)
        y_array = (R - xy_array[1]) * np.sin(xy_array[0] / R)

        for at, x, y in zip(self.atoms, x_array, y_array):
            at.coords = (x, y, at.coords[-1])

    def get_center_of_mass(self, unit="angstrom"):
        """Return the center of mass of the molecule (as a tuple). Returned coordinates are expressed in *unit*."""
        mass_array = np.array([atom.mass for atom in self])
        xyz_array = self.as_array().T
        xyz_array *= mass_array
        center = xyz_array.sum(axis=1)
        center /= mass_array.sum()
        return tuple(center * Units.conversion_ratio("angstrom", unit))

    def get_masses(self, unit: Optional[str] = "amu"):
        """Return list of masses, by default in atomic mass units."""
        unit_conversion_coeff = Units.convert(1.0, "amu", unit)
        return [at.mass * unit_conversion_coeff for at in self.atoms]

    def get_mass(self, unit="amu"):
        """Return the mass of the molecule, by default in atomic mass units."""
        return sum([at.mass for at in self.atoms]) * Units.convert(1.0, "amu", unit)

    def get_density(self):
        """Return the density in kg/m^3"""
        vol = self.unit_cell_volume(unit="angstrom") * 1e-30  # in m^3
        mass = self.get_mass(unit="kg")
        return mass / vol

    def set_density(self, density: float):
        """
        Applies a uniform strain so that the density becomes ``density`` kg/m^3
        """
        assert len(self.lattice) == 3
        current_density = self.get_density()  # in kg/m^3
        strain = (current_density / density) ** (1 / 3.0)
        strain -= 1.0
        self.apply_strain([strain, strain, strain, 0, 0, 0], voigt_form=True)

    def get_formula(self, as_dict=False):
        """Calculate the molecular formula of the molecule according to the Hill system.

        Here molecular formula is a dictionary with keys being atomic symbols. The value for each key is the number of atoms of that type. If *as_dict* is ``True``, that dictionary is returned. Otherwise, it is converted into a string::

            >>> mol = Molecule('Ubiquitin.xyz')
            >>> print(m.get_formula(True))
            {'N': 105, 'C': 378, 'O': 118, 'S': 1, 'H': 629}
            >>> print(m.get_formula(False))
            C378H629N105O118S1

        """
        occ = {}
        for atom in self:
            if atom.symbol not in occ:
                occ[atom.symbol] = 0
            occ[atom.symbol] += 1
        if as_dict:
            return occ

        def string_for_sym(sym, occ):
            if sym not in occ:
                return ""
            else:
                n = occ.pop(sym)
                return sym if n == 1 else f"{sym}{n}"

        if "C" in occ:
            # for organic molecules C and H come first
            return (
                string_for_sym("C", occ)
                + string_for_sym("H", occ)
                + "".join(string_for_sym(sym, occ) for sym in sorted(occ))
            )
        else:
            # for inorganic systems the order is strictly alphabetic
            return "".join(string_for_sym(sym, occ) for sym in sorted(occ))

    def get_inertia_matrix(self, length_unit: Optional[str] = "angstrom", mass_unit: Optional[str] = "amu"):
        """Get the moments of inertia matrix.

        Args:
            length_unit (str, optional): unit for distance. Defaults to 'angstrom'.
            mass_unit (str, optional): unit for mass. Defaults to 'amu'.

        Returns:
            np.ndarray: 3x3 matrix with the inertia matrix
        """
        com = self.get_center_of_mass(unit=length_unit)
        positions = self.as_array()
        positions -= np.array(com)  # translate center of mass to origin
        masses = self.get_masses(unit=mass_unit)

        x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
        I11 = np.sum(masses * (y**2 + z**2))
        I22 = np.sum(masses * (x**2 + z**2))
        I33 = np.sum(masses * (x**2 + y**2))
        I12 = -np.sum(masses * x * y)
        I13 = -np.sum(masses * x * z)
        I23 = -np.sum(masses * y * z)

        inertia_matrix = np.array([[I11, I12, I13], [I12, I22, I23], [I13, I23, I33]])
        return inertia_matrix

    def get_moments_of_inertia(
        self,
        eigen_vectors: Optional[bool] = False,
        length_unit: Optional[str] = "angstrom",
        mass_unit: Optional[str] = "amu",
    ):
        """Get the moments of inertia along the principal axes (in amu*angstrom**2 by default). They are computed from the eigenvalues of the symmetric inertial tensor. Optionally the eigenvectors can be returned.

        Args:
            eigen_vectors (bool, optional): return also the eigen_vectors. Defaults to False.
            length_unit (str, optional): unit for distance. Defaults to 'angstrom'.
            mass_unit (str, optional): unit for mass. Defaults to 'amu'.

        Returns:
            np.ndarray or Tuple(np.ndarray, np.ndarray): moments of inertia [(3,)] and optionally the eigenvectors
        """
        I = self.get_inertia_matrix(length_unit=length_unit, mass_unit=mass_unit)

        evals, evecs = np.linalg.eigh(I)
        if eigen_vectors:
            return evals, evecs.transpose()
        else:
            return evals

    def get_gyration_radius(self, unit: Optional[str] = "angstrom") -> float:
        """Return the gyration radius of the molecule by default in angstrom. It gives information about the overall dimensions of the rotating molecule around its center of mass.

        Args:
            unit (str, optional): unit for distance. Defaults to 'angstrom'.

        Returns:
            float: gyration radius of the molecule in unit.
        """
        moment_of_inertia = self.get_moments_of_inertia(length_unit=unit)
        magnitude_momento_inertia = np.linalg.norm(moment_of_inertia)
        gyration_radius = np.sqrt(magnitude_momento_inertia / self.get_mass())
        return gyration_radius

    def apply_strain(self, strain, voigt_form=False):
        """Apply a strain deformation to a periodic system (i.e. with a non-empty ``lattice`` attribute).
        The atoms in the unit cell will be strained accordingly, keeping the fractional atomic coordinates constant.

        If ``voigt_form=False``, *strain* should be a container with n*n numerical values, where n is the number of ``lattice`` vectors. It can be a list (tuple, numpy array etc.) listing matrix elements row-wise, either flat (e.g. ``[e_xx, e_xy, e_xz, e_yx, e_yy, e_yz, e_zx, e_zy, e_zz]``) or in two-level fashion (e.g. ``[[e_xx, e_xy, e_xz],[e_yx, e_yy, e_yz],[e_zx, e_zy, e_zz]]``).
        If ``voigt_form=True``, *strain* should be passed in Voigt form (for 3D periodic systems: ``[e_xx, e_yy, e_zz, gamma_yz, gamma_xz, gamma_xy]``; for 2D periodic systems: ``[e_xx, e_yy, gamma_xy]``; for 1D periodic systems: ``[e_xx]``  with e_xy = gamma_xy/2,...). Example usage::

            >>> graphene = Molecule('graphene.xyz')
            >>> print(graphene)
              Atoms:
                1         C      0.000000      0.000000      0.000000
                2         C      1.230000      0.710141      0.000000
              Lattice:
                    2.4600000000     0.0000000000     0.0000000000
                    1.2300000000     2.1304224900     0.0000000000
            >>> graphene.apply_strain([0.1,0.2,0.0], voigt_form=True)])
              Atoms:
                1         C      0.000000      0.000000      0.000000
                2         C      1.353000      0.852169      0.000000
              Lattice:
                    2.7060000000     0.0000000000     0.0000000000
                    1.3530000000     2.5565069880     0.0000000000
        """

        n = len(self.lattice)

        if n == 0:
            raise MoleculeError("apply_strain: can only be used for periodic systems.")

        if n in [1, 2] and self.align_lattice(convention="AMS"):
            raise MoleculeError(
                "apply_strain: the lattice vectors should follow the convention of AMS (i.e. for 1D-periodic systems the lattice vector should be along the x-axis, while for 2D-periodic systems the two vectors should be on the XY plane. Consider using the align_lattice function."
            )

        def from_voigt_to_matrix(strain_voigt, n):
            if len(strain_voigt) != n * (n + 1) / 2:
                raise MoleculeError(
                    "apply_strain: strain for %i-dim periodic system needs %i-sized vector in Voigt format"
                    % (n, n * (n + 1) / 2)
                )

            strain_matrix = np.diag(strain_voigt[:n])
            if n == 2:
                strain_matrix[1, 0] = strain_voigt[2] / 2.0
                strain_matrix[0, 1] = strain_voigt[2] / 2.0
            elif n == 3:
                strain_matrix[1, 2] = strain_voigt[3] / 2.0
                strain_matrix[2, 1] = strain_voigt[3] / 2.0
                strain_matrix[0, 2] = strain_voigt[4] / 2.0
                strain_matrix[2, 0] = strain_voigt[4] / 2.0
                strain_matrix[0, 1] = strain_voigt[5] / 2.0
                strain_matrix[1, 0] = strain_voigt[5] / 2.0
            return strain_matrix

        if voigt_form:
            strain = from_voigt_to_matrix(strain, n)
        else:
            try:
                strain = np.array(strain).reshape(n, n)
            except Exception:
                raise MoleculeError("apply_strain: could not convert the strain to a (%i,%i) numpy array" % (n, n))

        if n == 1:
            lattice_mat = np.array([[self.lattice[0][0]]])
        else:
            lattice_mat = np.array(self.lattice)[:n, :n]

        strained_lattice = lattice_mat.dot(np.eye(n) + strain)
        coords = self.as_array()
        frac_coords_transf = np.linalg.inv(lattice_mat.T)
        fractional_coords = coords[:, :n] @ frac_coords_transf.T
        coords[:, :n] = (strained_lattice.T @ fractional_coords.T).T

        self.from_array(coords)
        self.lattice = [tuple(vec + [0.0] * (3 - len(vec))) for vec in strained_lattice.tolist()]

    def map_to_central_cell(self, around_origin=True):
        """Maps all atoms to the original cell. If *around_origin=True* the atoms will be mapped to the cell with fractional coordinates [-0.5,0.5], otherwise to the the cell in which all fractional coordinates are in the [0:1] interval."""

        n = len(self.lattice)
        if n == 0:
            raise MoleculeError("map_to_central_cell: can only be used for periodic systems.")
        elif n == 1:
            lattice_mat = np.array([[self.lattice[0][0]]])
        else:
            lattice_mat = np.array(self.lattice)[:n, :n]

        coords = self.as_array()
        frac_coords_transf = np.linalg.inv(lattice_mat.T)
        fractional_coords = coords[:, :n] @ frac_coords_transf.T
        if around_origin:
            shift = -np.rint(fractional_coords)
        else:
            shift = -np.floor(fractional_coords)
        fractional_coords_new = fractional_coords + shift
        coords[:, :n] = (lattice_mat.T @ fractional_coords_new.T).T
        self.from_array(coords)

        if any(b.has_cell_shifts() for b in self.bonds):
            # Fix cell shifts for bonds for atoms that were moved.
            for b in self.bonds:
                # Check if the mapping has moved the bonded atoms relative to each other.
                at1, at2 = self.index(b)
                at1 = at1 - 1
                at2 = at2 - 1  # -1 because np.array is indexed from 0
                relshift = (shift[at2, :n] - shift[at1, :n]).astype(int)
                if not np.all(relshift == 0):
                    # Relative position has changed: cell shifts need updating!
                    if b.has_cell_shifts():
                        # Grab the original cell shifts from the suffix. An empty
                        # suffix means "0 0 0" if at least one atom has cell shifts. If
                        # no atom has cell shifts and all suffixes are empty, the bonds
                        # are always taken to be to the closest image, but that case is
                        # covered by to top level if (above) already.
                        try:
                            cell_shifts = np.array([int(cs) for cs in b.properties.suffix.split()])
                        except Exception:
                            raise MoleculeError("Cell shifts in bond suffix are not all integers.")
                        if cell_shifts.size != n:
                            raise MoleculeError("Wrong number of cell shifts in bond suffix.")
                    else:
                        cell_shifts = np.array([0 for i in range(n)])
                    cell_shifts_new = cell_shifts - relshift
                    if np.all(cell_shifts_new == 0):
                        # All 0 cell shifts are not written out explicitly
                        if "suffix" in b.properties:
                            del b.properties.suffix
                    else:
                        b.properties.suffix = " ".join(str(cs) for cs in cell_shifts_new)

    def perturb_atoms(self, max_displacement=0.01, unit="angstrom", atoms=None):
        """Randomly perturb the coordinates of the atoms in the molecule.

        Each Cartesian coordinate is displaced by a random value picked out of a uniform distribution in the interval *[-max_displacement, +max_displacement]* (converted to requested *unit*).

        By default, all atoms are perturbed. It is also possible to perturb only part of the molecule, indicated by *atoms* argument. It should be a list of atoms belonging to the molecule.
        """
        s = Units.convert(max_displacement, "angstrom", unit)

        if atoms is None:
            atoms = self.atoms

        for atom in atoms:
            atom.translate(np.random.uniform(-s, s, 3))

    def perturb_lattice(self, max_displacement=0.01, unit="angstrom", ams_convention=True):
        """Randomly perturb the lattice vectors.

        The Cartesian components of the lattice vectors are changed by a random value picked out of a uniform distribution in the interval *[-max_displacement, +max_displacement]* (converted to requested *unit*).

        If *ams_convention=True* then for 1D-periodic systems only the x-component of the lattice vector is perturbed, and for 2D-periodic systems only the xy-components of the lattice vectors are perturbed.
        """
        s = Units.convert(max_displacement, "angstrom", unit)
        n = len(self.lattice)

        if n == 0:
            raise MoleculeError("perturb_lattice can only be applied to periodic systems")

        for i, vec in enumerate(self.lattice):
            if ams_convention:
                # For 1D systems we only want to perturb the first number. For 2D systems only the first 2 numbers of each vector.
                perturbed_vec = np.array(vec) + np.concatenate((np.random.uniform(-s, s, n), np.zeros(3 - n)))
            else:
                perturbed_vec = np.array(vec) + np.random.uniform(-s, s, 3)
            self.lattice[i] = tuple(perturbed_vec)

    def substitute(
        self, connector, ligand, ligand_connector, bond_length=None, steps=12, cost_func_mol=None, cost_func_array=None
    ):
        """Substitute a part of this molecule with *ligand*.

        *connector* should be a pair of atoms that belong to this molecule and form a bond. The first atom of *connector* is the atom to which the  ligand will be connected. The second atom of *connector* is removed from the molecule, together with all "further" atoms connected to it (that allows, for example, to substitute the whole functional group with another). Using *connector* that is a part or a ring triggers an exception.

        *ligand_connector* is a *connector* analogue, but for *ligand*. IT describes the bond in the *ligand* that will be connected with the bond in this molecule described by *connector*.

        If this molecule or *ligand* don't have any bonds, :meth:`guess_bonds` is used.

        After removing all unneeded atoms, the *ligand* is translated to a new position, rotated, and connected by bond with the core molecule. The new |Bond| is added between the first atom of *connector* and the first atom of *ligand_connector*. The length of that bond can be adjusted with *bond_length* argument, otherwise the default is the sum of atomic radii taken from |PeriodicTable|.

        Then the *ligand* is rotated along newly created bond to find the optimal position. The full 360 degrees angle is divided into *steps* equidistant rotations and each such rotation is evaluated using a cost function. The orientation with the minimal cost is chosen.

        The default cost function is:

        .. math::

            \sum_{i \in mol, j\in lig} e^{-R_{ij}}

        A different cost function can be also supplied by the user, using one of the two remaining arguments: *cost_func_mol* or *cost_func_array*. *cost_func_mol* should be a function that takes two |Molecule| instances: this molecule (after removing unneeded atoms) and ligand in a particular orientation (also without unneeded atoms) and returns a single number (the lower the number, the better the fit). *cost_func_array* is analogous, but instead of |Molecule| instances it takes two numpy arrays (with dimensions: number of atoms x 3) with coordinates of this molecule and the ligand. If both are supplied, *cost_func_mol* takes precedence over *cost_func_array*.

        """
        try:
            _is_atom = [isinstance(i, Atom) and i.mol is self for i in connector]
            assert all(_is_atom) and len(_is_atom) == 2
        except (TypeError, AssertionError) as ex:
            raise MoleculeError(
                "substitute: connector argument must be a pair of atoms that belong to the current molecule"
            ).with_traceback(ex.__traceback__)

        try:
            _is_atom = [isinstance(i, Atom) and i.mol is ligand for i in ligand_connector]
            assert all(_is_atom) and len(_is_atom) == 2
        except (TypeError, AssertionError) as ex:
            raise MoleculeError(
                "substitute: ligand_connector argument must be a pair of atoms that belong to ligand"
            ).with_traceback(ex.__traceback__)

        _ligand = ligand.copy()
        _ligand_connector = [_ligand[ligand.index(atom)] for atom in ligand_connector]

        if len(self.bonds) == 0:
            self.guess_bonds()
        if len(_ligand.bonds) == 0:
            _ligand.guess_bonds()

        def dfs(atom, stay, go, delete, msg):
            for N in atom.neighbors():
                if N is stay:
                    if atom is go:
                        continue
                    raise MoleculeError("substitute: {} is a part of a cycle".format(msg))
                if N not in delete:
                    delete.add(N)
                    dfs(N, stay, go, delete, msg)

        stay, go = connector
        stay_lig, go_lig = _ligand_connector

        # remove 'go' and all connected atoms from self
        atoms_to_delete = {go}
        dfs(go, stay, go, atoms_to_delete, "connector")
        for atom in atoms_to_delete:
            self.delete_atom(atom)

        # remove 'go_lig' and all connected atoms from _ligand
        atoms_to_delete = {go_lig}
        dfs(go_lig, stay_lig, go_lig, atoms_to_delete, "ligand_connector")
        for atom in atoms_to_delete:
            _ligand.delete_atom(atom)

        # move the _ligand such that 'go_lig' is in (0,0,0) and rotate it to its desired position
        vec = np.array(stay.vector_to(go))
        vec_lig = np.array(go_lig.vector_to(stay_lig))

        _ligand.translate(go_lig.vector_to((0, 0, 0)))
        _ligand.rotate(rotation_matrix(vec_lig, vec))

        # rotate the _ligand along the bond to create 'steps' copies
        angles = [i * (2 * np.pi / steps) for i in range(1, steps)]
        axis_matrices = [axis_rotation_matrix(vec, angle) for angle in angles]
        xyz_ligand = _ligand.as_array()
        xyz_ligands = np.array([xyz_ligand] + [xyz_ligand @ matrix for matrix in axis_matrices])

        # move all the _ligand copies to the right position
        if bond_length is None:
            bond_length = stay.radius + stay_lig.radius
        vec *= bond_length / np.linalg.norm(vec)
        position = np.array(stay.coords) + vec
        trans_vec = np.array(stay_lig.vector_to(position))
        xyz_ligands += trans_vec

        # find the best _ligand orientation
        if cost_func_mol:
            best_score = np.inf
            for lig in xyz_ligands:
                _ligand.from_array(lig)
                score = cost_func_mol(self, _ligand)
                if score < best_score:
                    best_score = score
                    best_lig = lig
        else:
            xyz_self = self.as_array()
            if cost_func_array:
                best = np.argmin([cost_func_array(xyz_self, i) for i in xyz_ligands])
            else:
                a, b, c = xyz_ligands.shape
                dist_matrix = distance_array(xyz_ligands.reshape(a * b, c), xyz_self)
                dist_matrix.shape = a, b, -1
                best = np.sum(np.exp(-dist_matrix), axis=(1, 2)).argmin()
            best_lig = xyz_ligands[best]

        # add the best _ligand to the molecule
        _ligand.from_array(best_lig)
        self.add_molecule(_ligand)
        self.add_bond(stay, stay_lig)

    def map_atoms_to_bonds(self):
        """
        Corrects for lattice displacements along bonds

        Uses a breadth-first search to map the atoms bond by bond (wrap)
        """

        def get_translation_vectors(iat, neighbors, molcoords):
            """
            Compute the translation of the neighbors
            """
            ones = np.ones((len(neighbors), 3))
            boxarray = box.reshape((1, 3)) * ones
            dcoords = np.rint((molcoords[neighbors] - molcoords[iat].reshape((1, 3)) * ones) / boxarray) * boxarray
            return dcoords, -np.rint(dcoords / boxarray)

        def get_translated_vectors(iat, neighbors, molcoords):
            """
            Compute the translation vectors for non-orthorhombic boxes
            """
            ones = np.ones((len(neighbors), 3))
            rfractional = s_inv @ molcoords[neighbors].transpose()
            diffvec = (molcoords[neighbors] - molcoords[iat].reshape((1, 3)) * ones).transpose()
            shift = np.rint(s_inv @ diffvec)
            rnew = (s @ (rfractional - shift)).transpose()
            shift = shift.transpose()
            return rnew, -shift

        # Only do something for a periodic system
        if len(self.lattice) == 0:
            return

        # Get the lattice vectors, and make sure there are 3 of them
        latticevecs = []
        for i in range(3):
            if i < len(self.lattice):
                latticevecs.append(self.lattice[i])
            else:
                v = [1.0e10 if j == i else 0.0 for j in range(3)]
                latticevecs.append(v)
        latticevecs = np.array(latticevecs)

        # Get the inversion matrix from fractional coordinates to cartesian and vice versa
        s = latticevecs.transpose()
        s_inv = np.linalg.inv(s)
        box = np.sqrt((latticevecs**2).sum(axis=1))

        mol_indices = [indices for indices in self.get_molecule_indices()]
        # connectivity = self.get_connection_table()
        coords = self.as_array()
        for imol, atoms in enumerate(mol_indices):
            periodic = False
            molcoords = coords[atoms].copy()
            mol_indices = {iat: i for i, iat in enumerate(atoms)}
            explored = []
            to_explore = [0]
            for iat in to_explore:
                if iat in explored:
                    continue
                explored.append(iat)
                # neighbors = [mol_indices[i] for i in connectivity[atoms[iat]]]
                neighbors = [mol_indices[self.index(at) - 1] for at in self.neighbors(self.atoms[atoms[iat]])]
                unique_indices = [i for i, jat in enumerate(neighbors) if not jat in to_explore]

                # Compute the translation for all the neighbors
                if (box - latticevecs.sum(axis=1)).any() < 1e-10:
                    # Orthorhombic (faster)
                    dcoords, shift = get_translation_vectors(iat, neighbors, molcoords)
                    newcoords = molcoords[neighbors] - dcoords
                else:
                    # Non-orthorhombic
                    newcoords, shift = get_translated_vectors(iat, neighbors, molcoords)
                    dcoords = newcoords - molcoords[neighbors]

                # Check for bonds across lattice (A neighbor that was already moved will require movement again)
                if True in [abs(dcoords[i]).sum() > 1e-14 for i, jat in enumerate(neighbors) if jat in explored]:
                    periodic = True
                    break

                # Now move only the unique neighbors
                neighbors = [neighbors[i] for i in unique_indices]
                molcoords[neighbors] = newcoords[unique_indices]
                # Update the cell shifts
                for j, jat in enumerate(neighbors):
                    bond = self.find_bond(self.atoms[atoms[iat]], self.atoms[atoms[jat]])
                    if bond.has_cell_shifts():
                        cell_shifts = np.array([int(cs) for cs in bond.properties.suffix.split()]) + shift[j]
                        if np.all(cell_shifts == 0):
                            # All 0 cell shifts are not written out explicitly
                            if "suffix" in bond.properties:
                                del bond.properties.suffix
                        else:
                            bond.properties.suffix = " ".join(str(int(cs)) for cs in cell_shifts)
                to_explore += neighbors
            if not periodic:
                coords[atoms] = molcoords
        self.from_array(coords)

    # ===========================================================================
    # ==== Magic methods ========================================================
    # ===========================================================================

    def __repr__(self):
        # get_formula(), but with C,H first and counts of 1 not present in the string
        syms = sorted([at.symbol for at in self])
        for at in "HC":
            syms = syms.count(at) * [at] + [i for i in syms if i != at]
        uniq = list(dict.fromkeys(syms))  # preserves order
        cnts = [syms.count(i) for i in uniq]
        cnts = [str(i) if i > 1 else "" for i in cnts]
        s = "".join(f"{at}{cnt}" for at, cnt in zip(uniq, cnts))
        return f"{self.__class__.__name__}('{s}' at {hex(id(self))})"

    def __len__(self):
        """The length of the molecule is the number of atoms."""
        return len(self.atoms)

    def __str__(self):
        return self.str()

    def str(self, decimal=6):
        """Return a string representation of the molecule.

        Information about atoms is printed in ``xyz`` format fashion -- each atom in a separate, enumerated line. Then, if the molecule contains any bonds, they are printed. Each bond is printed in a separate line, with information about both atoms and bond order. Example:

        .. code-block:: none

                  Atoms:
                    1         N       0.00000       0.00000       0.38321
                    2         H       0.94218       0.00000      -0.01737
                    3         H      -0.47109       0.81595      -0.01737
                    4         H      -0.47109      -0.81595      -0.01737
                  Bonds:
                    (1)----1----(2)
                    (1)----1----(3)
                    (1)----1----(4)
        """
        s = "  Atoms: \n"
        for i, atom in enumerate(self.atoms, 1):
            s += ("%5i" % (i)) + atom.str(decimal=decimal) + "\n"
        if len(self.bonds) > 0:
            for j, atom in enumerate(self.atoms, 1):
                atom._tmpid = j
            s += "  Bonds: \n"
            for bond in self.bonds:
                s += "   (%d)--%1.1f--(%d)\n" % (bond.atom1._tmpid, bond.order, bond.atom2._tmpid)
            for atom in self.atoms:
                del atom._tmpid
        if self.lattice:
            s += "  Lattice:\n"
            for vec in self.lattice:
                s += "    {:16.10f} {:16.10f} {:16.10f}\n".format(*vec)
        return s

    def __iter__(self):
        """Iterate over atoms."""
        return iter(self.atoms)

    @overload
    def __getitem__(self, key: int) -> Atom: ...
    @overload
    def __getitem__(self, key: Tuple[int, int]) -> Bond: ...
    def __getitem__(self, key):
        """The bracket notation can be used to access atoms or bonds directly.

        If *key* is a single int (``mymol[i]``), return i-th atom of the molecule. If *key* is a pair of ints (``mymol[(i,j)]``), return the bond between i-th and j-th atom (``None`` if such a bond does not exist). Negative integers can be used to access atoms enumerated in the reversed order.

        This notation is read only: things like ``mymol[3] = Atom(...)`` are forbidden.

        Numbering of atoms within a molecule starts with 1.
        """
        if hasattr(key, "__index__"):  # Available in all "int-like" objects; see PEP 357
            if key == 0:
                raise MoleculeError("Numbering of atoms starts with 1")
            if key < 0:
                return self.atoms[key]
            return self.atoms[key - 1]

        try:
            i, j = key
            return self.find_bond(self[i], self[j])
        except TypeError as ex:
            raise MoleculeError(f"Molecule: argument ({repr(key)}) of invalid type inside []").with_traceback(
                ex.__traceback__
            )
        except ValueError as ex:
            raise MoleculeError(f"Molecule: argument ({repr(key)}) of invalid size inside []").with_traceback(
                ex.__traceback__
            )

    def __add__(self, other):
        """Create a new molecule that is a sum of this molecule and some *other* molecule::

            newmol = mol1 + mol2

        The new molecule has atoms, bonds and all other elements distinct from both components. The ``properties`` of ``newmol`` are a copy of the ``properties`` of ``mol1`` :meth:`soft_updated<scm.plams.core.settings.Settings.soft_update>` with the ``properties`` of ``mol2``.
        """
        m = self.copy()
        m += other
        return m

    def __iadd__(self, other):
        """Copy *other* molecule and add the copy to this one."""
        self.add_molecule(other, copy=True)
        return self

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __round__(self, ndigits=None):
        """Magic method for rounding this instance's Cartesian coordinates; called by the builtin :func:`round` function."""
        ndigits = 0 if ndigits is None else ndigits
        return self.round_coords(ndigits, inplace=False)

    def __getstate__(self) -> dict:
        """Returns the object which is to-be pickled by, *e.g.*, :func:`pickle.dump`.
        As :class:`Molecule` instances are heavily nested objects,
        pickling them can raise a :exc:`RecursionError`.
        This issue is herein avoided relying on the :meth:`Molecule.as_dict()` method.
        See `Pickling Class Instances <https://docs.python.org/3/library/pickle.html#pickling-class-instances>`_
        for more details.
        """  # noqa
        return self.as_dict()

    def __setstate__(self, state: dict) -> None:
        """Counterpart of :meth:`Molecule.__getstate__`; used for unpickling molecules."""
        try:
            mol_new = self.from_dict(state)
            self.__dict__ = mol_new.__dict__

        # Raised if *state* is the result of a pickled Molecule created prior to the introduction
        # of Molecule.__getstate__()
        except TypeError:
            self.__dict__ = state
            return

        # Molecule.from_dict() always returns a new instance
        # Simply steal this instance's attributes and changed its Atoms/Bonds parent Molecule
        for at in self.atoms:
            at.mol = self
        for bond in self.bonds:
            bond.mol = self

    # ===========================================================================
    # ==== Converters ===========================================================
    # ===========================================================================

    def as_dict(self):
        """Store all information about the molecule in a dictionary.

        The returned dictionary is, in principle, identical to ``self.__dict__`` of the current instance, apart from the fact that all |Atom| and |Bond| instances in ``atoms`` and ``bonds`` lists are replaced with dictionaries storing corresponding information.

        This method is a counterpart of :meth:`from_dict`.
        """
        mol_dict = copy.copy(self.__dict__)
        atom_indices = {id(a): i for i, a in enumerate(mol_dict["atoms"])}
        bond_indices = {id(b): i for i, b in enumerate(mol_dict["bonds"])}
        atom_dicts = [copy.copy(a.__dict__) for a in mol_dict["atoms"]]
        bond_dicts = [copy.copy(b.__dict__) for b in mol_dict["bonds"]]
        for a_dict in atom_dicts:
            a_dict["bonds"] = [bond_indices[id(b)] for b in a_dict["bonds"] if id(b) in bond_indices]
            del a_dict["mol"]
        for b_dict in bond_dicts:
            b_dict["atom1"] = atom_indices[id(b_dict["atom1"])]
            b_dict["atom2"] = atom_indices[id(b_dict["atom2"])]
            del b_dict["mol"]
        mol_dict["atoms"] = atom_dicts
        mol_dict["bonds"] = bond_dicts
        return mol_dict

    @classmethod
    def from_dict(cls, dictionary):
        """Generate a new |Molecule| instance based on the information stored in a *dictionary*.

        This method is a counterpart of :meth:`as_dict`.
        """
        mol = cls()
        mol.__dict__ = copy.copy(dictionary)
        atom_dicts = mol.atoms
        bond_dicts = mol.bonds
        mol.atoms = []
        mol.bonds = []
        for a_dict in atom_dicts:
            a = Atom()
            a.__dict__ = a_dict
            a.mol = mol
            a.bonds = []
            mol.add_atom(a)
        for b_dict in bond_dicts:
            b = Bond(None, None)
            b_dict["atom1"] = mol.atoms[b_dict["atom1"]]
            b_dict["atom2"] = mol.atoms[b_dict["atom2"]]
            b.__dict__ = b_dict
            b.mol = mol
            mol.add_bond(b)
        return mol

    @classmethod
    def from_elements(cls, elements):
        """Generate a new |Molecule| instance based on a list of *elements*.

        By default it sets all coordinates to zero
        """
        mol = cls()
        for el in elements:
            at = Atom(symbol=el, coords=(0.0, 0.0, 0.0))
            mol.add_atom(at)
        return mol

    @property
    def as_array(self):
        """
        Property that can either be called directly as a method: ``mol.as_array()`` or as a context manager ``with mol.as_array``. Take care of when to add the parentheses.

        Return cartesian coordinates of this molecule's atoms as a numpy array.

        *atom_subset* argument can be used to specify only a subset of atoms, it should be an iterable container with atoms belonging to this molecule.

        Returned value is a n*3 numpy array where n is the number of atoms in the whole molecule, or in *atom_subset*, if used.

        Alternatively, this property can be used in conjunction with the ``with`` statement,
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
        return self._as_array

    def from_array(self, xyz_array, atom_subset=None):
        """Update the cartesian coordinates of this |Molecule|, containing n atoms, with coordinates provided by a (≤n)*3 numpy array *xyz_array*.

        *atom_subset* argument can be used to specify only a subset of atoms, it should be an iterable container with atoms belonging to this molecule. It should have the same length as the first dimension of *xyz_array*.
        """
        atom_subset = atom_subset or self.atoms
        for at, (x, y, z) in zip(atom_subset, xyz_array):
            at.coords = (x, y, z)

    def __array__(self, dtype=None):
        """A magic method for constructing numpy arrays.

        This method ensures that passing a |Molecule| instance to numpy.array_ produces an array of Cartesian coordinates (see :meth:`.Molecule.as_array`).
        The array `data type`_ can, optionally, be specified in *dtype*.

        .. _numpy.array: https://docs.scipy.org/doc/numpy/reference/generated/numpy.array.html
        .. _`data type`: https://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html
        """
        ret = self.as_array()
        return ret.astype(dtype, copy=False)

    # ===========================================================================
    # ==== File/format IO =======================================================
    # ===========================================================================

    def readxyz(self, f, geometry=1, **other):
        """XYZ Reader:

        The xyz format allows to store more than one geometry of a particular molecule within a single file.
        In such cases the *geometry* argument can be used to indicate which (in order of appearance in the file) geometry to import.
        Default is the first one (*geometry* = 1).
        """

        def newatom(line):
            lst = line.split()
            shift = 1 if (len(lst) > 4 and lst[0] == str(i)) else 0
            num = lst[0 + shift]
            if isinstance(num, str):
                new_atom = Atom(symbol=num, coords=(lst[1 + shift], lst[2 + shift], lst[3 + shift]))
            else:
                new_atom = Atom(atnum=num, coords=(lst[1 + shift], lst[2 + shift], lst[3 + shift]))
            if len(lst) > shift + 4:
                new_atom.properties.suffix = " ".join(line.split()[shift + 4 :])
            self.add_atom(new_atom)

        def newlatticevec(line):
            lst = line.split()
            self.lattice.append((float(lst[1]), float(lst[2]), float(lst[3])))

        if isinstance(f, list):
            f = io.StringIO("\n".join(f))
            log(
                "WARNING: DEPRECATED A list was passed as 'f' argument to the Molecule.readxyz method. 'f' should be a file, and not a list. Tip: consider using io.StringIO if you want to pass a string as input to this function",
                1,
            )

        fr = geometry
        begin, first, nohead = True, True, False

        while True:
            last_pos = f.tell()
            line = f.readline()
            if line == "":
                break

            if first:
                if line.strip() == "":
                    continue
                first = False
                try:
                    n = int(line.strip())
                    fr -= 1
                except ValueError:
                    nohead = True
                    newatom(line)
            elif nohead:
                if line.strip() == "":
                    break
                if "VEC" in line.upper():
                    newlatticevec(line)
                else:
                    newatom(line)
            elif fr != 0:
                try:
                    n = int(line.strip())
                    fr -= 1
                except ValueError:
                    continue
            else:
                if begin:
                    begin = False
                    i = 1
                    if line:
                        self.properties["comment"] = line.rstrip()
                else:
                    if i <= n:
                        newatom(line)
                        i += 1
                    elif "VEC" in line.upper():
                        newlatticevec(line)
                    else:
                        # If we get here, it means that we just processed a line behind the last atom or VEC line
                        # If this xyz file contains more than one molecule, this line might be the header of the next molecule.
                        # Let's move back one step, so that if we call this function again the file pointer will be at the right position
                        f.seek(last_pos)
                        break
        if not nohead and fr > 0:
            raise FileError(f"readxyz: cannot read frame from {f.name}")

    def writexyz(self, f, space=16, decimal=8):
        """
        f: file
            An open file handle.

        See also the write method: `molecule.write("my_molecule.xyz")`

        example:

        .. code-block:: python

            with open(path_to_xyz_molecule_file, 'w') as f:
                molecule.writexyz(f)
        """
        f.write(str(len(self)) + "\n")
        if "comment" in self.properties:
            comment = self.properties["comment"]
            if isinstance(comment, list):
                comment = comment[0]
            f.write(comment)
        f.write("\n")

        for at in self.atoms:
            f.write(at.str(space=space, decimal=decimal) + "\n")

        for i, vec in enumerate(self.lattice, 1):
            f.write(f"VEC{i} " + " ".join([f"{v:>{space}.{decimal}f}" for v in vec]) + "\n")

    def readmol(self, f, **other):

        comment = []
        for i in range(4):
            line = f.readline().rstrip()
            if line:
                spl = line.split()
                if spl[-1].lower() == "V2000".lower():
                    if len(line) == 39:
                        natom = int(line[0:3])
                        nbond = int(line[3:6])
                    else:
                        natom = int(spl[0])
                        nbond = int(spl[1])
                    for j in range(natom):
                        atomline = f.readline().rstrip()
                        if len(atomline) == 69:
                            crd = (float(atomline[:10]), float(atomline[10:20]), float(atomline[20:30]))
                            symb = atomline[31:34].strip()
                        else:
                            tmp = atomline.split()
                            crd = tuple(map(float, tmp[0:3]))
                            symb = tmp[3]
                        self.add_atom(Atom(symbol=symb, coords=crd))
                    for j in range(nbond):
                        bondline = f.readline().rstrip()
                        if len(bondline) == 21:
                            at1 = int(bondline[0:3])
                            at2 = int(bondline[3:6])
                            ordr = int(bondline[6:9])
                        else:
                            tmp = bondline.split()
                            at1 = int(tmp[0])
                            at2 = int(tmp[1])
                            ordr = int(tmp[2])
                        if ordr == 4:
                            ordr = Bond.AR
                        if at1 > natom or at2 > natom:
                            at1 = int(bondline[0:3])
                            at2 = int(bondline[3:6])
                            ordr = int(bondline[6:9])
                            if ordr == 4:
                                ordr = Bond.AR
                        self.add_bond(Bond(atom1=self[at1], atom2=self[at2], order=ordr))
                    break
                elif spl[-1] == "V3000":
                    raise FileError("readmol: Molfile V3000 not supported. Please convert")
                else:
                    comment.append(line)
        if comment:
            self.properties["comment"] = comment

    def writemol(self, f, **other):
        commentblock = ["\n"] * 3
        if "comment" in self.properties:
            comment = self.properties["comment"]
            if isinstance(comment, str):
                commentblock[0] = comment + "\n"
            elif isinstance(comment, list):
                comment = comment[0:3]
                while len(comment) < 3:
                    comment.append("")
                commentblock = [a + b for a, b in zip(comment, commentblock)]
        f.writelines(commentblock)

        self.set_atoms_id()

        f.write("%3i %2i  0  0  0  0  0  0  0  0999 V2000\n" % (len(self.atoms), len(self.bonds)))
        for at in self.atoms:
            f.write("%10.4f %9.4f %9.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n" % (at.x, at.y, at.z, at.symbol))
        for bo in self.bonds:
            order = bo.order
            if order == Bond.AR:
                order = 4
            f.write("%3i %2i %2i  0  0  0  0\n" % (bo.atom1.id, bo.atom2.id, round(order)))
        self.unset_atoms_id()
        f.write("M  END\n")

    def readmol2(self, f, **other):

        bondorders = {"1": 1, "2": 2, "3": 3, "am": 1, "ar": Bond.AR, "du": 0, "un": 1, "nc": 0}
        mode = ("", 0)
        for i, line in enumerate(f):
            line = line.rstrip()
            if not line:
                continue
            elif line[0] == "#":
                continue
            elif line[0] == "@":
                line = line.partition(">")[2]
                if not line:
                    raise FileError("readmol2: Error in %s line %i: invalid @ record" % (f.name, str(i + 1)))
                mode = (line, i)

            elif mode[0] == "MOLECULE":
                pos = i - mode[1]
                if pos == 1:
                    self.properties["name"] = line
                elif pos == 3:
                    self.properties["type"] = line
                elif pos == 4:
                    self.properties["charge_type"] = line
                elif pos == 5:
                    self.properties["flags"] = line
                elif pos == 6:
                    self.properties["comment"] = line

            elif mode[0] == "ATOM":
                spl = line.split()
                if len(spl) < 6:
                    raise FileError("readmol2: Error in %s line %i: not enough values in line" % (f.name, str(i + 1)))
                symb = spl[5].partition(".")[0]
                crd = tuple(map(float, spl[2:5]))
                newatom = Atom(symbol=symb, coords=crd, name=spl[1], type=spl[5])
                if len(spl) > 6:
                    newatom.properties["subst_id"] = spl[6]
                if len(spl) > 7:
                    newatom.properties["subst_name"] = spl[7]
                if len(spl) > 8:
                    newatom.properties["charge"] = float(spl[8])
                if len(spl) > 9:
                    newatom.properties["flags"] = spl[9]
                self.add_atom(newatom)

            elif mode[0] == "BOND":
                spl = line.split()
                if len(spl) < 4:
                    raise FileError("readmol2: Error in %s line %i: not enough values in line" % (f.name, str(i + 1)))
                try:
                    atom1 = self.atoms[int(spl[1]) - 1]
                    atom2 = self.atoms[int(spl[2]) - 1]
                except IndexError:
                    raise FileError("readmol2: Error in %s line %i: wrong atom ID" % (f.name, str(i + 1)))
                newbond = Bond(atom1, atom2, order=bondorders[spl[3]])
                if len(spl) > 4:
                    for flag in spl[4].split("|"):
                        newbond.properties[flag] = True
                self.add_bond(newbond)

    def writemol2(self, f, **other):

        def write_prop(name, obj, separator, space=0, replacement=None):
            form_str = "%-" + str(space) + "s"
            if name in obj.properties:
                f.write(form_str % str(obj.properties[name]))
            elif replacement is not None:
                f.write(form_str % str(replacement))
            f.write(separator)

        f.write("@<TRIPOS>MOLECULE\n")
        write_prop("name", self, "\n")
        f.write("%i %i\n" % (len(self.atoms), len(self.bonds)))
        write_prop("type", self, "\n")
        write_prop("charge_type", self, "\n")
        write_prop("flags", self, "\n")
        write_prop("comment", self, "\n")

        f.write("\n@<TRIPOS>ATOM\n")
        for i, at in enumerate(self.atoms, 1):
            f.write("%5i " % (i))
            write_prop("name", at, " ", 5, at.symbol + str(i + 1))
            f.write("%10.4f %10.4f %10.4f " % at.coords)
            write_prop("type", at, " ", 5, at.symbol)
            write_prop("subst_id", at, " ", 5)
            write_prop("subst_name", at, " ", 7)
            write_prop("charge", at, " ", 6)
            write_prop("flags", at, "\n")
            at.id = i

        f.write("\n@<TRIPOS>BOND\n")
        for i, bo in enumerate(self.bonds, 1):
            f.write("%5i %5i %5i %4s" % (i, bo.atom1.id, bo.atom2.id, "ar" if bo.is_aromatic() else bo.order))
            write_prop("flags", bo, "\n")

        self.unset_atoms_id()

    def readpdb(self, f, geometry=1, **other):
        """PDB Reader:

        The pdb format allows to store more than one geometry of a particular molecule within a single file.
        In such cases the *geometry* argument can be used to indicate which (in order of appearance in the file) geometry to import.
        The default is the first one (*geometry* = 1).
        """
        pdb = PDBHandler(f)
        models = pdb.get_models()
        if geometry > len(models):
            raise FileError("readpdb: There are only %i geometries in %s" % (len(models), f.name))

        symbol_columns = [70, 6, 7, 8]
        for i in models[geometry - 1]:
            if i.name in ["ATOM  ", "HETATM"]:
                x = float(i.value[0][24:32])
                y = float(i.value[0][32:40])
                z = float(i.value[0][40:48])
                for n in symbol_columns:
                    symbol = i.value[0][n : n + 2].strip()
                    try:
                        atnum = PT.get_atomic_number(symbol)
                        break
                    except PTError:
                        if n == symbol_columns[-1]:
                            raise FileError(
                                "readpdb: Unable to deduce the atomic symbol in the following line:\n%s"
                                % (i.name + i.value[0])
                            )
                self.add_atom(Atom(atnum=atnum, coords=(x, y, z)))

        return pdb

    def writepdb(self, f, **other):
        pdb = PDBHandler()
        pdb.add_record(PDBRecord("HEADER"))
        model = []
        for i, at in enumerate(self.atoms, 1):
            s = "ATOM  %5i                   %8.3f%8.3f%8.3f                      %2s  " % (
                i,
                at.x,
                at.y,
                at.z,
                at.symbol.upper(),
            )
            model.append(PDBRecord(s))
        pdb.add_model(model)
        pdb.add_record(pdb.calc_master())
        pdb.add_record(PDBRecord("END"))
        pdb.write(f)

    def hydrogen_to_deuterium(self):
        """
        Modifies the current molecule so that all hydrogen atoms get mass 2.014 by modifying the atom.properties.mass
        """

        for at in self:
            if at.atnum == 1:
                at.properties.mass = 2.014

    @staticmethod
    def _mol_from_rkf_section(sectiondict):
        """Return a |Molecule| instance constructed from the contents of the whole ``.rkf`` file section, supplied as a dictionary returned by :meth:`KFFile.read_section<scm.plams.tools.kftools.KFFile.read_section>`."""
        from scm.plams.interfaces.adfsuite.ams import AMSJob

        ret = Molecule()
        coords = [sectiondict["Coords"][i : i + 3] for i in range(0, len(sectiondict["Coords"]), 3)]
        symbols = sectiondict["AtomSymbols"]
        # If the dictionary was read from memory and not from file, this is already a list
        if isinstance(symbols, str):
            symbols = symbols.split()
        for crd, sym in zip(coords, symbols):
            if sym.startswith("Gh."):
                isghost = True
                _, sym = sym.split(".", 1)
            else:
                isghost = False
            if "." in sym:
                elsym, name = sym.split(".", 1)
                newatom = Atom(symbol=elsym, coords=crd, unit="bohr")
                newatom.properties.name = name
            else:
                newatom = Atom(symbol=sym, coords=crd, unit="bohr")
            if isghost:
                newatom.properties.ghost = True
            ret.add_atom(newatom)
        if "fromAtoms" in sectiondict and "toAtoms" in sectiondict and "bondOrders" in sectiondict:
            fromAtoms = (
                sectiondict["fromAtoms"] if isinstance(sectiondict["fromAtoms"], list) else [sectiondict["fromAtoms"]]
            )
            toAtoms = sectiondict["toAtoms"] if isinstance(sectiondict["toAtoms"], list) else [sectiondict["toAtoms"]]
            bondOrders = (
                sectiondict["bondOrders"]
                if isinstance(sectiondict["bondOrders"], list)
                else [sectiondict["bondOrders"]]
            )

            for iBond, (fromAt, toAt, bondOrder) in enumerate(zip(fromAtoms, toAtoms, bondOrders)):
                b = Bond(ret[fromAt], ret[toAt], bondOrder)
                if "latticeDisplacements" in sectiondict:
                    nLatVec = int(sectiondict["nLatticeVectors"])
                    cellShifts = sectiondict["latticeDisplacements"][iBond * nLatVec : iBond * nLatVec + nLatVec]
                    if not all(cs == 0 for cs in cellShifts):
                        b.properties.suffix = " ".join(str(cs) for cs in cellShifts)
                ret.add_bond(b)
        if sectiondict["Charge"] != 0:
            ret.properties.charge = sectiondict["Charge"]
        if "nLatticeVectors" in sectiondict:
            ret.lattice = Units.convert(
                [
                    tuple(sectiondict["LatticeVectors"][i : i + 3])
                    for i in range(0, len(sectiondict["LatticeVectors"]), 3)
                ],
                "bohr",
                "angstrom",
            )
        if "EngineAtomicInfo" in sectiondict:
            if len(ret) == 1:
                # Just one atom: Need to make the list of length 1 explicitly.
                suffixes = [sectiondict["EngineAtomicInfo"]]
            elif "\x00" in sectiondict["EngineAtomicInfo"]:
                # AMS>2020: Separated with C NULL characters.
                suffixes = sectiondict["EngineAtomicInfo"].split("\x00")
            else:
                # AMS<=2019: Separated with new line characters
                suffixes = sectiondict["EngineAtomicInfo"].splitlines()
            for at, suffix in zip(ret, suffixes):
                if suffix:
                    at.properties.soft_update(AMSJob._atom_suffix_to_settings(suffix))
        return ret

    def forcefield_params_from_rkf(self, filename):
        """
        Read all force field data from a forcefield.rkf file into self

        * ``filename`` -- Name of the RKF file that contains ForceField data
        """
        from scm.plams.interfaces.adfsuite.ams import AMSJob
        from scm.plams.interfaces.adfsuite.forcefieldparams import forcefield_params_from_kf

        # Read atom types and charges
        kf = KFFile(filename)
        if not "AMSResults" in kf.sections():
            raise FileError("filename has to be a forcefield.rkf file from a GAFF atomtyping calculation.")
        charges, types, patch = forcefield_params_from_kf(kf)

        # Place the data into self (any force field parameter data already there will be overwritten)
        suffixes = [at.properties.suffix if "suffix" in at.properties else "" for at in self]
        suffixes = [suf.lower() for suf in suffixes]
        for i, at in enumerate(self.atoms):
            suffix = suffixes[i] + "forcefield.type=%s forcefield.charge=%f" % (types[i], charges[i])
            at.properties.soft_update(AMSJob._atom_suffix_to_settings(suffix))
        if patch is not None:
            self.properties.forcefieldpatch = patch

    def readrkf(self, filename: str, section: str = "Molecule", **other):
        kf = KFFile(filename)
        sectiondict = kf.read_section(section)
        self.__dict__.update(Molecule._mol_from_rkf_section(sectiondict).__dict__)
        for at in self.atoms:
            at.mol = self
        for bo in self.bonds:
            bo.mol = self

    def readin(self, f, **other):
        """Read a file containing a System block used in AMS driver input files."""
        if not input_parser_available:
            raise NotImplementedError(
                "Reading from System blocks from AMS input files requires an AMS installation to be available."
            )
        from scm.plams.interfaces.adfsuite.inputparser import InputParserFacade
        from scm.plams.interfaces.adfsuite.ams import AMSJob

        sett = Settings()
        sett.input.AMS = Settings(InputParserFacade().to_dict("ams", f.read(), string_leafs=True))
        if "System" not in sett.input.AMS:
            raise ValueError("No System block found in file.")
        sysname = other.get("sysname", "")
        mols = AMSJob.settings_to_mol(sett)
        if sysname not in mols:
            raise KeyError(f'No System block with id "{sysname}" found in file.')
        self.__dict__.update(mols[sysname].__dict__)
        for at in self.atoms:
            at.mol = self
        for bo in self.bonds:
            bo.mol = self

    def writein(self, f, **other):
        """Write the Molecule instance to a file as a System block from the AMS driver input files."""
        from scm.plams.interfaces.adfsuite.ams import AMSJob

        f.write(AMSJob(molecule={other.get("sysname", ""): self}).get_input())

    def read(self, filename, inputformat=None, **other):
        """Read molecular coordinates from a file.

        *filename* should be a string with a path to a file. If *inputformat* is not ``None``, it should be one of supported formats or engines (keys occurring in the class attribute ``_readformat``). Otherwise, the format is deduced from the file extension. For files without an extension the `xyz` format is used.

        All *other* options are passed to the chosen format reader.
        """

        if inputformat is None:
            _, extension = os.path.splitext(filename)
            inputformat = extension.strip(".") if extension else "xyz"
        if inputformat in self.__class__._readformat:
            if inputformat == "rkf":
                return self.readrkf(filename, **other)
            else:
                with open(filename, "r") as f:
                    ret = self._readformat[inputformat](self, f, **other)
                return ret
        else:
            raise MoleculeError(f"read: Unsupported file format '{inputformat}'")

    def write(self, filename, outputformat=None, mode="w", **other):
        """Write molecular coordinates to a file.

        *filename* should be a string with a path to a file. If *outputformat* is not ``None``, it should be one of supported formats or engines (keys occurring in the class attribute ``_writeformat``). Otherwise, the format is deduced from the file extension. For files without an extension the `xyz` format is used.

        *mode* can be either 'w' (overwrites the file if the file exists) or 'a' (appends to the file if the file exists).

        All *other* options are passed to the chosen format writer.
        """

        if not mode in ["w", "a"]:
            raise ValueError(f"invalid mode {mode}")

        if outputformat is None:
            _, extension = os.path.splitext(filename)
            outputformat = extension.strip(".") if extension else "xyz"
        if outputformat in self.__class__._writeformat:
            with open(filename, mode) as f:
                self._writeformat[outputformat](self, f, **other)
        else:
            raise MoleculeError(f"write: Unsupported file format '{outputformat}'")

    # Support for the ASE engine is added if available by interfaces.molecules.ase
    _readformat = {"xyz": readxyz, "mol": readmol, "mol2": readmol2, "pdb": readpdb, "rkf": readrkf}
    _writeformat = {"xyz": writexyz, "mol": writemol, "mol2": writemol2, "pdb": writepdb}
    if input_parser_available:
        _readformat["in"] = readin
        _writeformat["in"] = writein

    def add_hatoms(self) -> "Molecule":
        """
        Adds missing hydrogen atoms to the current molecule.
        Returns a new Molecule instance.

        Example::

            >>> o = Molecule()
            >>> o.add_atom(Atom(atnum=8))
            >>> print(o)
              Atoms:
                1         O      0.000000       0.000000       0.000000
            >>> h2o = o.add_hatoms()
            >>> print(h2o)
              Atoms:
                1         O      0.000000       0.000000       0.000000
                2         H     -0.109259       0.893161       0.334553
                3         H      0.327778       0.033891      -0.901672

        """
        from subprocess import DEVNULL, Popen
        from tempfile import NamedTemporaryFile

        with NamedTemporaryFile(mode="w+", suffix=".xyz", delete=False) as f_in:
            self.writexyz(f_in)
            f_in.close()
            with NamedTemporaryFile(mode="w+", suffix=".xyz", delete=False) as f_out:
                f_out.close()
                amsprep = os.path.join(os.environ["AMSBIN"], "amsprep")
                p = Popen(
                    f"sh {amsprep} -t SP -m {f_in.name} -addhatoms -exportcoordinates {f_out.name}",
                    shell=True,
                    stdout=DEVNULL,
                )
                p.communicate()
                retmol = self.__class__(f_out.name)
                os.remove(f_out.name)
            os.remove(f_in.name)
        return retmol

    @staticmethod
    def rmsd(mol1, mol2, ignore_hydrogen=False, return_rotmat=False, check=True):
        """
        Uses the
        `Kabsch algorithm <https://en.wikipedia.org/wiki/Kabsch_algorithm>`_ to align and
        calculate the root-mean-square deviation of two systems' atomic positions.

        Assumes all elements and their order is the same in both systems, will check this if `check == True`.

        :Returns:

        rmsd : float
            Root-mean-square-deviation of atomic coordinates
        rotmat : ndarray
            If `return_rotmat` is `True`, will additionally return the rotation matrix
            that aligns `mol2` onto `mol1`.
        """

        def kabsch(x, y, rotmat=False):
            """
            Rotate a set of points `y` such that they are aligned with `x`
            using the Kabsch algorithm (thanks to Toon for the idea).
            Same as `scipy.spatial.transform.Rotation.align_vectors`.
            """
            x -= x.mean(0)
            y -= y.mean(0)
            covar = np.dot(x.T, y)
            U, S, Vt = np.linalg.svd(covar)
            det = np.ones(x.shape[-1])
            det[-1] = np.sign(np.linalg.det(Vt.dot(U)))
            R = np.einsum("ji,j,kj", Vt, det, U)
            y = y @ R
            rmsd = np.linalg.norm(x - y) * np.sqrt(1 / len(x))
            return (rmsd, R) if rotmat is True else rmsd

        assert len(mol1) == len(mol2), "Can only calculate the RMSD of same-sized molecules"
        if check:
            nums1 = np.array([i.atnum for i in mol1])
            nums2 = np.array([i.atnum for i in mol2])
            assert (
                nums1 == nums2
            ).all(), "\nAtoms are not the same (or not in the same order). Use `check==False` if you do not care about this.\n"
        if ignore_hydrogen is True:
            mol1 = [at.coords for at in mol1 if at.symbol != "H"]
            mol2 = [at.coords for at in mol2 if at.symbol != "H"]
        return kabsch(np.array(mol1), np.array(mol2), rotmat=return_rotmat)

    def align2mol(self, molecule_ref, ignore_hydrogen: bool = False, watch: bool = False):
        """
        align the molecule to a reference molecule, they should be same molecule type and same order of atoms
        it is an wrapper of the rmsd methods
        if watch = True show the molecules (before and after) in a Jupyter notebook
        """
        if watch:
            mol_initial = self.copy()
        rmsd_value, R = Molecule.rmsd(
            self, molecule_ref, ignore_hydrogen=ignore_hydrogen, return_rotmat=True, check=True
        )
        self.rotate(R, lattice=False)

        center_m_var = self.get_center_of_mass(unit="angstrom")
        center_m_ref = molecule_ref.get_center_of_mass(unit="angstrom")
        vector = np.array(center_m_ref) - np.array(center_m_var)
        self.translate(vector, unit="angstrom")

        if watch:
            import matplotlib.pyplot as plt

            from scm.plams.tools.plot import plot_molecule

            fig, ax = plt.subplots(1, 2)
            ax[0].set_title("before alignment")
            plot_molecule(molecule_ref, ax=ax[0], keep_axis=True)
            plot_molecule(mol_initial, ax=ax[0], keep_axis=True)
            ax[1].set_title("after alignment")
            plot_molecule(molecule_ref, ax=ax[1], keep_axis=True)
            plot_molecule(self, ax=ax[1], keep_axis=True)
            print(f"Root mean square deviation: {rmsd_value:0.3} Ang")

    @property
    def numbers(self) -> "np.ndarray":
        """Return an array of all atomic numbers in the Molecule. Can also be used to set all numbers at once."""
        return np.array([i.atnum for i in self])

    @numbers.setter
    def numbers(self, values):
        if len(values) != len(self):
            raise ValueError(
                f"Number of elements in array ({len(values)}) does not match the molecule size ({len(self)})."
            )
        [setattr(at, "atnum", value) for at, value in zip(self, values)]

    @property
    def symbols(self) -> "np.ndarray":
        """Return an array of all atomic symbols in the Molecule. Can also be used to set all symbols at once."""
        return np.array([i.symbol for i in self])

    @symbols.setter
    def symbols(self, values):
        if len(values) != len(self):
            raise ValueError(
                f"Number of elements in array ({len(values)}) does not match the molecule size ({len(self)})."
            )
        [setattr(at, "symbol", value) for at, value in zip(self, values)]

    def _get_bond_id(self, at1, at2, id_type):
        """
        at1: Atom in this molecule
        at2: Atom in this molecule
        id_type: str, 'IDname' or 'symbol'
        This function is called by get_unique_bonds()

        Returns: a 2-tuple, the key and a bool. The bool is True if the order was reversed.
        """
        at1key = getattr(at1, id_type)
        at2key = getattr(at2, id_type)
        if at1key < at2key:
            return at1key + "-" + at2key, False
        else:
            return at2key + "-" + at1key, True

    def get_unique_bonds(self, ignore_dict=None, id_type="symbol", index_start=1):
        """

        Returns a dictionary of all unique bonds in this molecule, where the
        key is the identifier and the value is a 2-tuple containing the 1-based
        indices of the atoms making up the bond (or 0-based indices if
        index_start == 0).

        ignore_dict : dict
            Bonds already existing in ignore_dict (as defined by the keys) will not be added to the returned dictionary

            Example: if id_type == 'symbol' and ignore_dict has a key 'C-C', then no C-C bond will be added to the return dictionary.

        id_type: str
            'symbol': The atomic symbols become the keys, e.g. 'C-H' (alphabetically sorted)

            'IDname': The IDname from molecule.set_local_labels() become the keys, e.g. 'an4va8478432bfl471baf74-knrq78jhkhq78fak111nf' (alphabetically sorted). Note: You must first call Molecule.set_local_labels()


        index_start : int
            If 1, indices are 1-based. If 0, indices are 0-based.

        """
        ret = {}
        ignore_dict = ignore_dict or {}
        for at in self:
            for b in at.bonds:
                other_atom = b.other_end(at)
                bondid, reverse = self._get_bond_id(at, other_atom, id_type)
                if bondid not in ignore_dict and bondid not in ret:
                    if reverse:
                        ret[bondid] = self.index(other_atom) - 1 + index_start, self.index(at) - 1 + index_start
                    else:
                        ret[bondid] = self.index(at) - 1 + index_start, self.index(other_atom) - 1 + index_start

        return ret

    def _get_angle_id(self, at1, at2, at3, id_type):
        at1key = getattr(at1, id_type)
        at2key = getattr(at2, id_type)
        at3key = getattr(at3, id_type)
        if at1key < at3key:
            return at1key + "-" + at2key + "-" + at3key, False
        else:
            return at3key + "-" + at2key + "-" + at1key, True

    def get_unique_angles(self, ignore_dict=None, id_type="symbol", index_start=1):
        """

        Returns a dictionary of all unique angles in this molecule, where the
        key is the identifier and the value is a 3-tuple containing the 1-based
        indices of the atoms making up the angle (or 0-based indices if
        index_start == 0). The central atom is the second atom.

        ignore_dict : dict
            Angles already existing in ignore_dict (as defined by the keys) will not be added to the returned dictionary

            Example: if id_type == 'symbol' and ignore_dict has a key 'C-C-C', then no C-C-C angle will be added to the return dictionary.

        id_type: str
            'symbol': The atomic symbols become the keys, e.g. 'C-C-C' (alphabetically sorted, the central atom in the middle)

            'IDname': The IDname from molecule.set_local_labels() become the keys, e.g. 'an4va8478432bfl471baf74-knrq78jhkhq78fak111nf-mf42918vslahf879bakfhk' (alphabetically sorted, the central atom in middle). Note: You must first call Molecule.set_local_labels()


        index_start : int
            If 1, indices are 1-based. If 0, indices are 0-based.

        """
        ret = {}
        ignore_dict = ignore_dict or {}
        for at in self:
            for b in at.bonds:
                at2 = b.other_end(at)
                for b2 in at2.bonds:
                    at3 = b2.other_end(at2)
                    if at == at3:
                        continue
                    angleid, reverse = self._get_angle_id(at, at2, at3, id_type)
                    if angleid not in ignore_dict and angleid not in ret:
                        if reverse:
                            ret[angleid] = (
                                self.index(at3) - 1 + index_start,
                                self.index(at2) - 1 + index_start,
                                self.index(at) - 1 + index_start,
                            )
                        else:
                            ret[angleid] = (
                                self.index(at) - 1 + index_start,
                                self.index(at2) - 1 + index_start,
                                self.index(at3) - 1 + index_start,
                            )

        return ret

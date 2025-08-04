#!/usr/bin/env python

from scm.plams.core.settings import Settings
from scm.plams.mol.atom import Atom
from scm.plams.mol.molecule import Bond, Molecule
from scm.plams.tools.periodic_table import PT
from scm.plams.trajectories.sdffile import SDFTrajectoryFile, get_molecule

__all__ = ["SDFHistoryFile"]


class SDFHistoryFile(SDFTrajectoryFile):
    """
    Class representing an SDF file containing a molecular simulation history with varying numbers of atoms

    An instance of this class has the following attributes:

    *   ``file_object`` -- A Python :py:class:`file` object, referring to the actual SDF file
    *   ``position``    -- The frame to which the cursor is currently pointing in the SDF file
    *   ``mode``        -- Designates whether the file is in read or write mode ('r' or 'w')
    *   ``elements``    -- The elements of the atoms in the system at the current frame

    An |SDFHistoryFile| object behaves very similar to a regular file object.
    It has read and write methods (:meth:`read_next` and :meth:`write_next`)
    that read and write from/to the position of the cursor in the ``file_object`` attribute.
    If the file is in read mode, an additional method :meth:`read_frame` can be used that moves
    the cursor to any frame in the file and reads from there.
    The amount of information stored in memory is kept to a minimum, as only information from the current frame
    is ever stored.

    Reading and writing to and from the files can be done as follows::

        >>> from scm.plams import SDFHistoryFile

        >>> sdf = SDFHistoryFile('old.sdf')
        >>> mol = sdf.get_plamsmol()

        >>> sdfout = SDFHistoryFile('new.sdf',mode='w')

        >>> for i in range(sdf.get_length()) :
        >>>     crd,cell = sdf.read_frame(i,molecule=mol)
        >>>     sdfout.write_next(molecule=mol)

    The above script reads information from the SDF file ``old.sdf`` into the |Molecule| object ``mol``
    in a step-by-step manner.
    The |Molecule| object is then passed to the :meth:`write_next` method of the new |SDFHistoryFile|
    object corresponding to the new sdf file ``new.sdf``.

    The exact same result can also be achieved by iterating over the instance as a callable

        >>> sdf = SDFHistoryFile('old.sdf')
        >>> mol = sdf.get_plamsmol()

        >>> sdfout = SDFHistoryFile('new.sdf',mode='w')

        >>> for crd,cell in sdf(mol) :
        >>>     sdfout.write_next(molecule=mol)

    This procedure requires all coordinate information to be passed to and from the |Molecule| object
    for each frame, which can be time-consuming.
    It is therefore also possible to bypass the |Molecule| object when reading through the frames::

        >>> sdf = SDFHistoryFile('old.sdf')

        >>> sdfout = SDFHistoryFile('new.sdf',mode='w')

        >>> for crd,cell in sdf :
        >>>     sdfout.write_next(coords=crd,elements=sdf.elements)

    By default the write mode will create a minimal version of the SDF file, containing only elements
    and coordinates.
    Additional information can be written to the file by supplying additional arguments
    to the :meth:`write_next` method.
    The additional keywords `step` and `energy` trigger the writing of a remark containing
    the molecule name, the step number, the energy, and the lattice vectors.

        >>> mol = Molecule('singleframe.sdf')

        >>> sdfout = SDFHistoryFile('new.sdf',mode='w')
        >>> sdfout.set_name('MyMol')

        >>> sdfout.write_next(molecule=mol, step=0, energy=5.)
    """

    def __init__(self, filename, mode="r", fileobject=None, ntap=None):
        """
        Initiates an SDFHistoryFile object

        * ``filename``   -- The path to the SDF file
        * ``mode``       -- The mode in which to open the SDF file ('r' or 'w')
        * ``fileobject`` -- Optionally, a file object can be passed instead (filename needs to be set to None)
        * ``ntap``       -- If the file is in write mode, the number of atoms can be passed here
        """
        SDFTrajectoryFile.__init__(self, filename, mode, fileobject, ntap)

        self.input_elements = self.elements[:]

    def _read_coordinates(self, molecule):
        """
        Read the coordinates from file, and place them in the molecule
        """
        # Read lines until the end
        lines = []
        while True:
            line = self.file_object.readline()
            if len(line) == 0:
                return None, None  # End of file is reached
            lines.append(line)
            if len(line) < 4:
                continue
            if line[:4] == "$$$$":
                break

        # Get the mol part
        mol, restlines = get_molecule(lines)
        elements = [at.symbol for at in mol.atoms]

        # Get the coordinates and cell
        cell = mol.lattice
        if len(cell) == 0:
            cell = None
        if len(mol.bonds) > 0:
            conect = {}
            for bond in mol.bonds:
                iat = min(mol.index(bond))
                jat = max(mol.index(bond))
                if not iat in conect:
                    conect[iat] = []
                conect[iat].append(jat)
            self.conect = conect

        # If the elements changed, update the molecule
        if elements != self.elements:
            self.elements = elements
            self.coords = mol.as_array()
            # Rebuild the molecule (bonds will disappear for now)
            if isinstance(molecule, Molecule):
                for at in reversed(molecule.atoms):
                    molecule.delete_atom(at)
                molecule.properties = Settings()
                for el in elements:
                    atom = Atom(PT.get_atomic_number(el))
                    molecule.add_atom(atom)
        else:
            self.coords[:, :] = mol.as_array()

        # Read the additional data
        if self.include_historydata:
            historydata = {}
            # First find all entries (entries can run over multiple lines)
            entries = [i for i, line in enumerate(restlines[:-1]) if line[:4] == ">  <"] + [len(restlines) - 1]
            for i, iline in enumerate(entries[:-1]):
                key = restlines[iline].split("<")[1].split(">")[0]
                value = "".join(restlines[iline + 1 : entries[i + 1] - 1])
                value = value.strip()
                # Try to turn this into a float or integer?
                if value.isdigit():
                    value = int(value)
                else:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                historydata[key] = value
            self.historydata = historydata

        if isinstance(molecule, Molecule):
            self._set_plamsmol(self.coords, cell, molecule)

        return self.coords, cell

    def write_next(self, coords=None, molecule=None, cell=[0.0, 0.0, 0.0], conect=None, historydata=None):
        """
        Write frame to next position in trajectory file

        * ``coords``   -- A list or numpy array of (``ntap``,3) containing the system coordinates
        * ``molecule`` -- A molecule object to read the molecular data from
        * ``cell``     -- A set of lattice vectors or cell diameters
        * ``conect``   -- A dictionary containing connectivity info (not used)
        * ``historydata`` -- A dictionary containing additional variables to be written to the comment line

        The ``historydata`` dictionary can contain for example:
        ('Step','Energy'), the frame number and the energy respectively

        .. note::

                Either ``coords`` or ``molecule`` are mandatory arguments
        """
        if isinstance(molecule, Molecule):
            coords, cell, elements = self._read_plamsmol(molecule)[:3]
            self.elements = elements
        cell = self._convert_cell(cell)

        if not isinstance(molecule, Molecule):
            # Create the molecule?
            molecule = Molecule()
            for el, crd in zip(self.elements, coords):
                atom = Atom(symbol=el, coords=crd)
                molecule.add_atom(atom)
            if cell is not None:
                molecule.lattice = cell.tolist()
            # Add the bonds
            bondlist = []
            if conect is not None:
                for iat, neighbors in conect.items():
                    for jat in neighbors:
                        indices = tuple(sorted([iat, jat]))
                        if not indices in bondlist:
                            bondlist.append(indices)
                            bond = Bond(molecule.atoms[iat], molecule.atoms[jat])
                            molecule.add_bond(bond)

        self._write_moldata(molecule, historydata)

        self.position += 1

#!/usr/bin/env python

from scm.plams.core.settings import Settings
import numpy
from scm.plams.mol.atom import Atom
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.periodic_table import PT
from scm.plams.trajectories.xyzfile import XYZTrajectoryFile, data_from_xyzcomment

__all__ = ["XYZHistoryFile"]


class XYZHistoryFile(XYZTrajectoryFile):
    """
    Class representing an XYZ file containing a molecular simulation history with varying numbers of atoms

    An instance of this class has the following attributes:

    *   ``file_object`` -- A Python :py:class:`file` object, referring to the actual XYZ file
    *   ``position``    -- The frame to which the cursor is currently pointing in the XYZ file
    *   ``mode``        -- Designates whether the file is in read or write mode ('r' or 'w')
    *   ``elements``    -- The elements of the atoms in the system at the current frame

    An |XYZHistoryFile| object behaves very similar to a regular file object.
    It has read and write methods (:meth:`read_next` and :meth:`write_next`)
    that read and write from/to the position of the cursor in the ``file_object`` attribute.
    If the file is in read mode, an additional method :meth:`read_frame` can be used that moves
    the cursor to any frame in the file and reads from there.
    The amount of information stored in memory is kept to a minimum, as only information from the current frame
    is ever stored.

    Reading and writing to and from the files can be done as follows::

        >>> from scm.plams import XYZHistoryFile

        >>> xyz = XYZHistoryFile('old.xyz')
        >>> mol = xyz.get_plamsmol()

        >>> xyzout = XYZHistoryFile('new.xyz',mode='w')

        >>> for i in range(xyz.get_length()) :
        >>>     crd,cell = xyz.read_frame(i,molecule=mol)
        >>>     xyzout.write_next(molecule=mol)

    The above script reads information from the XYZ file ``old.xyz`` into the |Molecule| object ``mol``
    in a step-by-step manner.
    The |Molecule| object is then passed to the :meth:`write_next` method of the new |XYZHistoryFile|
    object corresponding to the new xyz file ``new.xyz``.

    The exact same result can also be achieved by iterating over the instance as a callable

        >>> xyz = XYZHistoryFile('old.xyz')
        >>> mol = xyz.get_plamsmol()

        >>> xyzout = XYZHistoryFile('new.xyz',mode='w')

        >>> for crd,cell in xyz(mol) :
        >>>     xyzout.write_next(molecule=mol)

    This procedure requires all coordinate information to be passed to and from the |Molecule| object
    for each frame, which can be time-consuming.
    It is therefore also possible to bypass the |Molecule| object when reading through the frames::

        >>> xyz = XYZHistoryFile('old.xyz')

        >>> xyzout = XYZHistoryFile('new.xyz',mode='w')

        >>> for crd,cell in xyz :
        >>>     xyzout.write_next(coords=crd,elements=xyz.elements)

    By default the write mode will create a minimal version of the XYZ file, containing only elements
    and coordinates.
    Additional information can be written to the file by supplying additional arguments
    to the :meth:`write_next` method.
    The additional keywords `step` and `energy` trigger the writing of a remark containing
    the molecule name, the step number, the energy, and the lattice vectors.

        >>> mol = Molecule('singleframe.xyz')

        >>> xyzout = XYZHistoryFile('new.xyz',mode='w')
        >>> xyzout.set_name('MyMol')

        >>> xyzout.write_next(molecule=mol, step=0, energy=5.)
    """

    def __init__(self, filename, mode="r", fileobject=None, ntap=None):
        """
        Initiates an XYZHistoryFile object

        * ``filename``   -- The path to the XYZ file
        * ``mode``       -- The mode in which to open the XYZ file ('r' or 'w')
        * ``fileobject`` -- Optionally, a file object can be passed instead (filename needs to be set to None)
        * ``ntap``       -- If the file is in write mode, the number of atoms can be passed here
        """
        XYZTrajectoryFile.__init__(self, filename, mode, fileobject, ntap)

        self.input_elements = self.elements[:]

    def _is_endoffile(self):
        """
        If the end of file is reached, return coords and cell as None
        """
        end = False
        line = self.file_object.readline()
        if len(line) == 0:
            end = True
            return end
        nats = int(line.split()[0])
        for i in range(nats + 1):
            line = self.file_object.readline()
            if len(line) == 0:
                end = True
                break
        for i in range(self.nveclines):
            line = self.file_object.readline()
            if len(line) == 0:
                end = True
                break
        return end

    def _read_coordinates(self, molecule):
        """
        Read the coordinates at current step
        """
        # Find the number of atoms
        line = self.file_object.readline()
        if len(line) == 0:
            return None, None  # End of file is reached
        nats = int(line.split()[0])
        line = self.file_object.readline()

        # Handle the comment line
        cell = None
        historydata = data_from_xyzcomment(line)
        if "Lattice" in historydata:
            cell = historydata["Lattice"]
            del historydata["Lattice"]
        if self.include_historydata:
            self.historydata = historydata

        # Read coordinates and elements
        coords = []
        elements = []
        for i in range(nats):
            line = self.file_object.readline()
            words = line.split()
            coords.append([float(w) for w in words[1:4]])
            elements.append(words[0])

        # If the elements changed, update the molecule
        if elements != self.elements:
            self.elements = elements
            self.coords = numpy.array(coords)
            # Rebuild the molecule (bonds will disappear for now)
            if isinstance(molecule, Molecule):
                for at in reversed(molecule.atoms):
                    molecule.delete_atom(at)
                molecule.properties = Settings()
                for el in elements:
                    atom = Atom(PT.get_atomic_number(el))
                    molecule.add_atom(atom)
        else:
            self.coords[:] = coords

        # Possibly read lattice
        lattice = []
        for i in range(self.nveclines):
            line = self.file_object.readline()
            words = line.split()
            lattice.append([float(w) for w in words[1:]])
        if cell is None and len(lattice) > 0:
            cell = lattice

        # Assign the data to the molecule object
        if isinstance(molecule, Molecule):
            self._set_plamsmol(self.coords, cell, molecule)

        return coords, cell

    def write_next(
        self, coords=None, molecule=None, elements=None, cell=[0.0, 0.0, 0.0], conect=None, historydata=None
    ):
        """
        Write frame to next position in trajectory file

        * ``coords``   -- A list or numpy array of (``ntap``,3) containing the system coordinates
        * ``molecule`` -- A molecule object to read the molecular data from
        * ``elements`` -- The element symbols of the atoms in the system
        * ``cell``     -- A set of lattice vectors or cell diameters
        * ``energy``   -- An energy value to be written to the remark line
        * ``conect``   -- A dictionary containing connectivity info (not used)
        * ``historydata`` -- A dictionary containing additional variables to be written to the comment line

        The ``historydata`` dictionary can contain for example:
        ('Step','Energy'), the frame number and the energy respectively

        .. note::

                Either ``coords`` and ``elements`` or ``molecule`` are mandatory arguments
        """
        if isinstance(molecule, Molecule):
            coords, cell, elements = self._read_plamsmol(molecule)[:3]
        self.elements = elements
        cell = self._convert_cell(cell)

        self._write_moldata(coords, cell, historydata)

        self.position += 1

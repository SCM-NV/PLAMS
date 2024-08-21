#!/usr/bin/env python

import io

from scm.plams.core.errors import TrajectoryError
import numpy
from scm.plams.mol.atom import Atom
from scm.plams.mol.bond import Bond
from scm.plams.mol.molecule import Molecule
from scm.plams.trajectories.trajectoryfile import TrajectoryFile

__all__ = ["SDFTrajectoryFile", "create_sdf_string"]


class SDFTrajectoryFile(TrajectoryFile):
    """
    Class representing an SDF file containing a molecular trajectory

    An instance of this class has the following attributes:

    *   ``file_object`` -- A Python :py:class:`file` object, referring to the actual SDF file
    *   ``position``    -- The frame to which the cursor is currently pointing in the SDF file
    *   ``mode``        -- Designates whether the file is in read or write mode ('r' or 'w')
    *   ``ntap``        -- The number of atoms in the molecular system (needs to be constant throughout)
    *   ``elements``    -- The elements of the atoms in the system (needs to be constant throughout)

    An |SDFTrajectoryFile| object behaves very similar to a regular file object.
    It has read and write methods (:meth:`read_next` and :meth:`write_next`)
    that read and write from/to the position of the cursor in the ``file_object`` attribute.
    If the file is in read mode, an additional method :meth:`read_frame` can be used that moves
    the cursor to any frame in the file and reads from there.
    The amount of information stored in memory is kept to a minimum, as only information from the current frame
    is ever stored.

    Reading and writing to and from the files can be done as follows::

        >>> from scm.plams import SDFTrajectoryFile

        >>> sdf = SDFTrajectoryFile('old.sdf')
        >>> mol = sdf.get_plamsmol()

        >>> sdfout = SDFTrajectoryFile('new.sdf',mode='w')

        >>> for i in range(sdf.get_length()) :
        >>>     crd,cell = sdf.read_frame(i,molecule=mol)
        >>>     sdfout.write_next(molecule=mol)

    The above script reads information from the SDF file ``old.sdf`` into the |Molecule| object ``mol``
    in a step-by-step manner.
    The |Molecule| object is then passed to the :meth:`write_next` method of the new |SDFTrajectoryFile|
    object corresponding to the new sdf file ``new.sdf``.

    The exact same result can also be achieved by iterating over the instance as a callable

        >>> sdf =SDFZTrajectoryFile('old.sdf')
        >>> mol = sdf.get_plamsmol()

        >>> sdfout = SDFTrajectoryFile('new.sdf',mode='w')

        >>> for crd,cell in sdf(mol) :
        >>>     sdfout.write_next(molecule=mol)

    This procedure requires all coordinate information to be passed to and from the |Molecule| object
    for each frame, which can be time-consuming.
    It is therefore also possible to bypass the |Molecule| object when reading through the frames::

        >>> sdf = SDFTrajectoryFile('old.sdf')

        >>> sdfout = SDFTrajectoryFile('new.sdf',mode='w')
        >>> sdfout.set_elements(sdf.get_elements())

        >>> for crd,cell in sdf :
        >>>     sdfout.write_next(coords=crd)
        >>> sdfout.close()

    By default the write mode will create a minimal version of the SDF file, containing only elements
    and coordinates.
    Additional information can be written to the file by supplying additional arguments
    to the :meth:`write_next` method.
    The additional keywords `step` and `energy` trigger the writing of a remark containing
    the molecule name, the step number, the energy, and the lattice vectors.

        >>> mol = Molecule('singleframe.xyz')

        >>> sdfout = SDFTrajectoryFile('new.sdf',mode='w')
        >>> sdfout.set_name('MyMol')

        >>> sdfout.write_next(molecule=mol, step=0, energy=5.)
    """

    def __init__(self, filename, mode="r", fileobject=None, ntap=None):
        """
        Initiates an SDFTrajectoryFile object

        * ``filename``   -- The path to the SDF file
        * ``mode``       -- The mode in which to open the SDF file ('r' or 'w')
        * ``fileobject`` -- Optionally, a file object can be passed instead (filename needs to be set to None)
        * ``ntap``       -- If the file is in write mode, the number of atoms needs to be passed here
        """
        TrajectoryFile.__init__(self, filename, mode, fileobject, ntap)

        # SDF specific attributes
        self.name = "PlamsMol"
        self.conect = None

        # Specific SDF stuff
        self.include_historydata = False
        self.historydata = None

        # Required setup before frames can be read/written
        if self.mode == "r":
            self._read_header()
        elif self.mode == "a":
            self._move_cursor_to_append_pos()

    def store_historydata(self):
        """
        Additional data should be read from/written to file
        """
        self.include_historydata = True

    def set_name(self, name):
        """
        Sets the name of the system, in case an extensive write is requested

        *   ``name`` -- A string containing the name of the molecule
        """
        self.name = name

    def _read_header(self):
        """
        Set up info required for reading frames
        """
        # Read the first molecule
        lines = []
        while True:
            line = self.file_object.readline()
            if len(line) == 0:
                raise Exception("End not found")
            lines.append(line)
            if len(line) < 6:
                continue
            if line[:6] == "M  END":
                break
        mol, _ = get_molecule(lines)
        self.ntap = len(mol)
        self.elements = [at.symbol for at in mol]
        if self.coords.shape == (0, 3):
            self.coords = numpy.zeros((self.ntap, 3))

        # Rewind
        self.file_object.seek(0)

    def read_next(self, molecule=None, read=True):
        """
        Reads coordinates from the current position of the cursor and returns it

        * ``molecule`` -- |Molecule| object in which the new coordinates need to be stored
        * ``read``     -- If set to False the cursor will move to the next frame without reading
        """
        if not read and not self.firsttime:
            return self._move_cursor_without_reading()

        cell = None
        # Read the coordinates
        crd, cell = self._read_coordinates(molecule)
        if crd is None:
            return None, None  # End of file is reached

        if self.firsttime:
            self.firsttime = False

        self.position += 1

        return self.coords, cell

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
        if len(mol) != self.ntap:
            raise TrajectoryError("Number of atoms changes. Not implemented.")

        # Get the coordinates and cell
        cell = mol.lattice
        if len(cell) == 0:
            cell = None
        self.coords[:, :] = mol.as_array()
        if len(mol.bonds) > 0:
            conect = {}
            for bond in mol.bonds:
                iat = min(mol.index(bond))
                jat = max(mol.index(bond))
                if not iat in conect:
                    conect[iat] = []
                conect[iat].append((jat, bond.order))
            self.conect = conect

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

    def _is_endoffile(self):
        """
        If the end of file is reached, return coords and cell as None
        """
        end = False
        while True:
            line = self.file_object.readline()
            if len(line) == 0:
                end = True
                break
            if len(line) < 4:
                continue
            if line[:4] == "$$$$":
                break
        return end

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
            if self.position == 0:
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
                    for t in neighbors:
                        jat = t
                        if isinstance(t, tuple):
                            jat = t[0]
                        indices = tuple(sorted([iat, jat]))
                        if not indices in bondlist:
                            bondlist.append(indices)
                            bond = Bond(molecule[iat], molecule[jat])
                            molecule.add_bond(bond)

        self._write_moldata(molecule, historydata)

        self.position += 1

    def _write_moldata(self, molecule, historydata):
        """
        Write all molecular info to file
        """
        block = create_sdf_string(molecule, self.position, historydata)
        self.file_object.write(block)

    def _rewind_to_first_frame(self):
        """
        Rewind the file to the first frame
        """
        self.file_object.seek(0)
        self.firsttime = True
        self.position = 0

    def _rewind_n_frames(self, nframes):
        """
        Rewind the file by nframes frames
        """
        new_frame = self.position - nframes
        self._rewind_to_first_frame()
        for i in range(new_frame):
            self.read_next(read=False)


def create_sdf_string(molecule, step=None, historydata=None):
    """
    Write an SDF entry based on the elements and the coordinates of the atoms
    """
    energy = 0.0
    if historydata is not None:
        if "Energy" in historydata:
            energy = historydata["Energy"]
        # The conformer case
        elif "energies" in historydata:
            energy = historydata["energies"]
    else:
        historydata = {}

    if "Step" in historydata:
        step = historydata["Step"]

    block = "Energy = %.10f kcal/mol\n" % (energy)
    f = io.StringIO()
    molecule.writemol(f)
    f.seek(0)
    text = f.read()
    f.close()
    text = "\n".join(text.split("\n")[1:])
    block += text
    for key, item in historydata.items():
        text = ">  <%s>  (%i)\n" % (key, step)
        text += f"{item}\n\n"
        block += text
    block += "$$$$\n"
    return block


def get_molecule(lines):
    """
    Read a molecule object from the lines for a single frame
    """
    end = 0
    for i, line in enumerate(lines):
        if line[:6] == "M  END":
            end = i + 1
            break
    moltext = "".join(lines[:end])
    f = io.StringIO(moltext)
    mol = Molecule()
    mol.readmol(f)
    f.close()
    restlines = lines[end:]
    return mol, restlines

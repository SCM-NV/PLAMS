#!/usr/bin/env python

from scm.plams.core.errors import PlamsError
import numpy
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.kftools import KFFile
from scm.plams.tools.periodic_table import PeriodicTable
from scm.plams.tools.units import Units
from scm.plams.trajectories.trajectoryfile import TrajectoryFile

__all__ = ["RKFTrajectoryFile", "write_general_section", "write_molecule_section"]

bohr_to_angstrom = Units.conversion_ratio("bohr", "angstrom")


class RKFTrajectoryFile(TrajectoryFile):
    """
    Class representing an RKF file containing a molecular trajectory

    An instance of this class has the following attributes:

    *   ``file_object`` -- A PLAMS |KFFile| object, referring to the actual RKF file
    *   ``position``    -- The frame to which the cursor is currently pointing in the RKF file
    *   ``mode``        -- Designates whether the file is in read or write mode ('rb' or 'wb')
    *   ``ntap``        -- The number of atoms in the molecular system (needs to be constant throughout)
    *   ``elements``    -- The elements of the atoms in the system (needs to be constant throughout)
    *   ``conect``      -- The connectivity information of the current frame
    *   ``mddata``      -- Read mode only: A dictionary containing data from the MDHistory section in the RKF file
    *   ``read_lattice``-- Read mode only: Wether the lattice vectors will be read from the file
    *   ``read_bonds``  -- Wether the connectivity information will be read from the file
    *   ``saving_freq`` -- How often the 'wb' file is written (default: only when :meth:`close` is called)

    An |RKFTrajectoryFile| object behaves very similar to a regular file object.
    It has read and write methods (:meth:`read_next` and :meth:`write_next`)
    that read and write from/to the position of the cursor in the ``file_object`` attribute.
    If the file is in read mode, an additional method :meth:`read_frame` can be used that moves
    the cursor to any frame in the file and reads from there.
    The amount of information stored in memory is kept to a minimum, as only information from the latest frame
    is ever stored.

    Reading and writing to and from the files can be done as follows::

        >>> from scm.plams import RKFTrajectoryFile

        >>> rkf = RKFTrajectoryFile('ams.rkf')
        >>> mol = rkf.get_plamsmol()

        >>> rkfout = RKFTrajectoryFile('new.rkf',mode='wb')

        >>> for i in range(rkf.get_length()) :
        >>>     crd,cell = rkf.read_frame(i,molecule=mol)
        >>>     rkfout.write_next(molecule=mol)
        >>> rkfout.close()

    The above script reads information from the RKF file ``ams.rkf`` into the |Molecule| object ``mol``
    in a step-by-step manner.
    The |Molecule| object is then passed to the :meth:`write_next` method of the new |RKFTrajectoryFile|
    object corresponding to the new rkf file ``new.rkf``.

    The exact same result can also be achieved by iterating over the instance as a callable

        >>> rkf = RKFTrajectoryFile('ams.rkf')
        >>> mol = rkf.get_plamsmol()

        >>> rkfout = RKFTrajectoryFile('new.rkf',mode='wb')

        >>> for crd,cell in rkf(mol) :
        >>>     rkfout.write_next(molecule=mol)
        >>> rkfout.close()

    This procedure requires all coordinate information to be passed to and from the |Molecule| object
    for each frame, which can be time-consuming.
    Some time can be saved by bypassing the |Molecule| object::

        >>> rkf = RKFTrajectoryFile('ams.rkf')

        >>> rkfout = RKFTrajectoryFile('new.rkf',mode='wb')
        >>> rkfout.set_elements(rkf.get_elements())

        >>> for crd,cell in rkf :
        >>>     rkfout.write_next(coords=crd,cell=cell,conect=rkf.conect)
        >>> rkfout.close()

    The only mandatory argument to the :meth:`write_next` method is ``coords``.
    Further time can be saved by setting the ``read_lattice`` and ``read_bonds`` variables to False.

    By default the write mode will create a minimal version of the RKF file, containing only elements,
    coordinates, lattice, and connectivity information.
    This minimal file format can be read by AMSMovie.

    It is possible to store additional information, such as energies, velocities, and charges.
    To enable this, the method :meth:`store_mddata` needs to be called after creation,
    and a dictionary of mddata needs to be passed to the :meth:`write_next` method.
    When that is done, the AMS trajectory analysis tools can be used on the file.
    Restarting an MD run with such a file is however currently not possible::

        >>> rkf = RKFTrajectoryFile('ams.rkf')
        >>> rkf.store_mddata()
        >>> mol = rkf.get_plamsmol()

        >>> rkf_out = RKFTrajectoryFile('new.rkf',mode='wb')
        >>> rkf_out.store_mddata(rkf)

        >>> for i in range(len(rkf)) :
        >>>         crd,cell = rkf.read_frame(i,molecule=mol)
        >>>         rkf_out.write_next(molecule=mol,mddata=rkf.mddata)
        >>> rkf_out.close()
    """

    def __init__(self, filename, mode="rb", fileobject=None, ntap=None):
        """
        Initiates an RKFTrajectoryFile object

        * ``filename``   -- The path to the RKF file
        * ``mode``       -- The mode in which to open the RKF file ('rb' or 'wb')
        * ``fileobject`` -- Optionally, a file object can be passed instead (filename needs to be set to None)
        * ``ntap``       -- If the file is in write mode, the number of atoms needs to be passed here
        """
        # TODO: If the mddata option is set to True, then the file created here works with AMSMovie and the analysis tools.
        #      To also make is work for restarts, two things have to be added:
        #      1. The final velocities have to be converted from bohr/fs to bohr/au (1/41.341373336493)
        #         and stored in MDResuts%EndVelocities
        #      2. The final coordinates need to be copied to the Molecule section.

        self.position = 0
        if filename is not None:
            # fileobject = KFFile(filename,autosave=False,keep_file_open=True)
            fileobject = KFFile(filename, autosave=False)
            # fileobject = KFFile(filename,autosave=False,fastsave=True)
            # This fastsave option (no copying) was not worth it, so I removed it.
            if fileobject is None:
                raise PlamsError("KFFile %s not found." % (filename))
        self.file_object = fileobject
        self.mode = mode

        self.ntap = 0
        if ntap is not None:
            self.ntap = ntap
        self.firsttime = True
        self.coords = numpy.zeros((self.ntap, 3))  # Only for reading purposes,
        # to avoid creating the array each time
        # PLAMS molecule related settings
        self.elements = ["H"] * self.ntap
        self.current_molecule = None
        self.store_molecule = True  # Even if True, the molecule attribute is only stored during iteration

        # RKF specific attributes
        self.program = "trajectory"
        self.nvecs = 3
        self.latticevecs = numpy.zeros((3, 3))
        self.read_lattice = True  # Reading time can be saved by skipping the lattice info
        self.read_bonds = True
        self.cell = numpy.zeros((3, 3))
        self.conect = None
        self.timestep = None
        self.saving_freq = None  # By default the 'wb' file is only written upon closing
        # Saving more often is much slower.
        self.include_mddata = False
        self.mdhistory_name = "MDHistory"
        self.mddata = None
        self.mdunits = None
        self.mdblocksize = None
        self.mditems = None
        self.mdblockitems = None
        self._mdblock = {}
        self.include_historydata = False  # Any additional data along the history section will be stored
        self.historydata = None
        self.historyitems = None

        # Skip to the trajectory part of the file (only if in read mode, because coords are required in header)
        if self.mode == "rb":
            self._read_header()
        elif self.mode == "wb":
            new_sections = ["InputMolecule", "Molecule", "History", self.mdhistory_name, "MDResults"]
            new_sections += ["SystemVersionHistory"]
            sections = self.file_object.sections()
            if len(sections) > 0:
                for secname in sections:
                    if secname in new_sections or "ChemicalSystem" in secname:
                        self.file_object.delete_section(secname)
            #       raise PlamsError ('RKF file %s already exists'%(filename))
        elif self.mode == "ab":
            self._move_cursor_to_append_pos()
        else:
            raise PlamsError('Mode %s is invalid. Only "rb" and "wb" are allowed.' % (self.mode))

    def store_mddata(self, rkf=None):
        """
        Read/write an MDHistory section

        * ``rkf`` -- If in write mode an RKFTrajectoryFile object in read mode needs to be passed to extract unit info
        """
        self.include_mddata = True
        if "r" in self.mode or "a" in self.mode:
            self._set_mddata_units()
            self._set_mddata_items()
        elif "w" in self.mode:
            if rkf is not None:
                self.timestep = rkf.timestep
                self._set_mdunits(rkf.mdunits)

    def store_historydata(self):
        """
        Read/write non-standard entries in the History section
        """
        self.include_historydata = True

        if "r" in self.mode:
            # Set the History items
            sections = self.file_object.get_skeleton()
            section = "History"
            item_keys = [kn for kn in sections[section] if "ItemName" in kn]
            items = [self.file_object.read(section, kn) for kn in item_keys]
            standard_items = [
                "Coords",
                "nLatticeVectors",
                "LatticeVectors",
                "Bonds.Index",
                "Bonds.Atoms",
                "Bonds.Orders",
            ]
            self.historyitems = [item for item in items if not item in standard_items]

    def close(self, override_molecule_section_with_last_frame=True):
        """
        Execute all prior commands and cleanly close and garbage collect the RKF file
        """
        # Write the step info
        if self.timestep is not None and self.mode == "wb":
            self.file_object.write("MDResults", "StartStep", 0)
            self.file_object.write("MDResults", "StartTime[fs]", 0.0)
            nsteps = self.get_length()
            self.file_object.write("MDResults", "EndStep", nsteps - 1)
            self.file_object.write("MDResults", "EndTime[fs]", (nsteps - 1) * self.timestep)

        # Write to file
        if self.mode == "wb":
            if override_molecule_section_with_last_frame:
                # First write the last frame into the Molecule section
                self._rewrite_molecule()
            # Then write to file
            self.file_object.save()
        del self

    def _rewrite_molecule(self):
        """
        Overwrite the molecule section with the latest frame
        """
        molecule = self.get_plamsmol()
        crd, cell = self.read_last_frame(molecule=molecule)
        self._write_molecule_section(crd, cell, molecule=molecule)

    def _read_header(self, molecule_section="Molecule"):
        """
        Set up info required for reading frames
        """
        if not molecule_section in self.file_object:
            return
        self.elements = self.file_object.read(molecule_section, "AtomSymbols")
        # If read from memory and not from file (write mode), it is already a list
        if isinstance(self.elements, str):
            self.elements = self.elements.split()
        self.elements = [el.split(".")[0] for el in self.elements]
        if (self.mdhistory_name, "Time(1)") in self.file_object:
            times = self.file_object.read(self.mdhistory_name, "Time(1)")
            if isinstance(times, list):
                self.timestep = times[1]
        self.ntap = len(self.elements)
        self.coords = numpy.zeros((self.ntap, 3))

        # Set the lattice info
        if (molecule_section, "LatticeVectors") in self.file_object:
            self.latticevecs = numpy.array(self.file_object.read(molecule_section, "LatticeVectors"))
            self.nvecs = int(len(self.latticevecs) / 3)  # Why did I remove this line locally?!
            self.latticevecs = self.latticevecs.reshape((self.nvecs, 3))

    def _set_mddata_units(self):
        """
        Get the units for the mddata, if those are to be read
        """
        # Look for the items
        section = self.mdhistory_name
        sections = self.file_object.get_skeleton()
        if not self.mdhistory_name in sections:
            self.mdunits = {}
            return
        item_keys = [kn for kn in sections[section] if "ItemName" in kn]
        items = [self.file_object.read(section, kn) for kn in item_keys]

        # Get the data for each item
        unit_dic = {}
        for item in items:
            if "%s(units)" % (item) in self.file_object.get_skeleton()[section]:
                unit_dic[item] = self.file_object.read(section, "%s(units)" % (item))

        self.mdunits = unit_dic

    def _set_mddata_items(self):
        """
        Get all the items for the mddatam if those are to be read
        """
        sections = self.file_object.get_skeleton()
        section = self.mdhistory_name
        if not self.mdhistory_name in sections:
            self.mditems = []
            self.mdblockitems = []
            self.mdblocksize = 100
            return
        blocksize = self.file_object.read(section, "blockSize")
        item_keys = [kn for kn in sections[section] if "ItemName" in kn]
        items = [self.file_object.read(section, kn) for kn in item_keys]
        blockitems = []
        for item in items:
            dim = self.file_object.read(section, "%s(dim)" % (item))
            if dim == 1 and not self.file_object.read(section, "%s(perAtom)" % (item)):
                is_blockitem = True
                if (section, "%s(1)" % (item)) in self.file_object:
                    if isinstance(self.file_object.read(section, "%s(1)" % (item)), str):
                        is_blockitem = False
                if is_blockitem:
                    blockitems.append(item)
        items = [item for item in items if not item in blockitems]
        self.mdblocksize = blocksize
        self.mditems = items
        self.mdblockitems = blockitems

    def _move_cursor_to_append_pos(self):
        """
        Get the instance ready for appending
        """
        self._read_header()
        self.position = self.get_length()

    def _write_header(self, coords, cell, molecule=None):
        """
        Write Molecule info to file (elements, periodicity)
        """
        # First write the general section
        if "General" not in self.file_object:
            write_general_section(self.file_object, self.program)

        # Then write the input molecule
        self._update_celldata(cell)
        self._write_molecule_section(coords, cell, molecule=molecule)
        self._write_molecule_section(coords, cell, section="InputMolecule", molecule=molecule)
        if self.include_mddata:
            # Start setting up the MDHistory section as well
            self.mdblocksize = 100
            self.file_object.write(self.mdhistory_name, "blockSize", 100)

        # Now make sure that it is possible to read from the file as well
        self._read_header()

    def _update_celldata(self, cell):
        """
        Use the newly supplied cell to update the dimensionality of the system
        """
        shape = numpy.array(cell).shape
        if len(shape) == 2:
            self.nvecs = shape[0]
            self.latticevecs = numpy.zeros((self.nvecs, 3))  # Not really necessay

    def _write_molecule_section(self, coords, cell, section="Molecule", molecule=None):
        """
        Write the molecule section
        """
        write_molecule_section(self.file_object, coords, cell, self.elements, section, molecule)

    def _set_mdunits(self, mdunits):
        """
        Store the dictionary with MD Units
        """
        if self.include_mddata:
            self.mdunits = mdunits

    def get_plamsmol(self):
        """
        Extracts a PLAMS molecule object from the RKF file
        """
        section_dict = self.file_object.read_section("InputMolecule")
        if len(section_dict) == 0:
            section_dict = self.file_object.read_section("Molecule")
        plamsmol = Molecule._mol_from_rkf_section(section_dict)
        return plamsmol

    def read_frame(self, i, molecule=None):
        """
        Reads the relevant info from frame ``i`` and returns it, or stores it in ``molecule``

        * ``i``        -- The frame number to be read from the RKF file
        * ``molecule`` -- |Molecule| object in which the new coordinates need to be stored
        """
        # Read the cell data
        cell = None
        if self.read_lattice:
            try:
                cell = self._read_cell_data(i)
            except (KeyError, AttributeError):
                pass

        # Read the bond data
        conect = None
        if self.read_bonds:
            conect = self._read_bond_data(section="History", step=i)
        self.conect = conect

        # Read the coordinates, and possible pass them to molecule
        try:
            self._read_coordinates(i, molecule, cell)
            # This has changed self.coords behind the scenes
        except (KeyError, AttributeError):
            return None, None

        # Read and store any additional data in the history section
        if self.include_historydata:
            self._store_historydata_for_step(i)
        # Read and store all MDData for this frame
        try:
            if self.include_mddata:
                self._store_mddata_for_step(i)
        except AttributeError:  # this is triggered when self.file_object is None triggered via self.close()
            pass
        # Finalize
        if self.firsttime:
            self.firsttime = False

        self.position = i
        return self.coords, cell

    def _read_coordinates(self, i, molecule, cell):
        """
        Read the coordinates at step i, and possible pass them to molecule
        """
        if not self.coords.shape == (self.ntap, 3):
            raise PlamsError("coords attribute has been changed outside the class")
        coords = self.coords.reshape(self.ntap * 3)
        coords[:] = self.file_object.read("History", "Coords(%i)" % (i + 1))
        coords *= bohr_to_angstrom
        # This has changed self.coords behind the scenes

        # Create the molecule
        if isinstance(molecule, Molecule):
            cell_reduced = None
            if cell is not None:
                cell_reduced = cell[: self.nvecs]
            self._set_plamsmol(self.coords, cell_reduced, molecule)

    def _read_cell_data(self, i):
        """
        Read the cell data at step i
        """
        if not ("History", "LatticeVectors(%i)" % (i + 1)) in self.file_object:
            return None
        latticevecs = self.latticevecs.reshape(self.nvecs * 3)
        latticevecs[:] = self.file_object.read("History", "LatticeVectors(%i)" % (i + 1))  # * bohr_to_angstrom
        latticevecs *= bohr_to_angstrom
        # This changed self.latticevecs behind the scenes
        # self.cell[:self.nvecs] = latticevecs
        self.cell[: self.nvecs] = self.latticevecs
        cell = self.cell
        return cell

    def _read_bond_data(self, section, step=None):
        """
        Read the bond data from the rkf file
        """
        conect = None
        try:
            step_txt = ""
            if step is not None:
                step_txt = "(%i)" % (step + 1)
            if not ("History", "Bonds.Index%s" % (step_txt)) in self.file_object:
                return conect
            indices = self.file_object.read(section, "Bonds.Index%s" % (step_txt))
            connection_table = self.file_object.read(section, "Bonds.Atoms%s" % (step_txt))
            if isinstance(connection_table, int):
                connection_table = [connection_table]
            bond_orders = self.file_object.read(section, "Bonds.Orders%s" % (step_txt))
            if isinstance(bond_orders, float):
                bond_orders = [bond_orders]
            conect = {}
            for i, (start, end) in enumerate(zip(indices[:-1], indices[1:])):
                if end - start > 0:
                    # conect[i+1] = connection_table[start-1:end-1]
                    conect[i + 1] = [
                        (ia, o)
                        for ia, o in zip(connection_table[start - 1 : end - 1], bond_orders[start - 1 : end - 1])
                    ]
            # Now correct the connection table
            # conect = self.symmetrize_conect(conect)
        except (KeyError, AttributeError):
            pass
        return conect

    def get_regular_connection_table(self):
        """
        Get the connection table without the bond orders
        """
        conect = self.symmetrize_conect(self.conect)  # noqa, ToDo: review, return value or remove method?

    def _store_mddata_for_step(self, istep):
        """
        Store the data from the MDHistory section
        """
        if "w" in self.mode:
            return
        self.mddata = {}

        # Get the data for each item
        section = self.mdhistory_name
        for item in self.mditems:
            if (section, "%s(%i)" % (item, istep + 1)) in self.file_object:
                self.mddata[item] = self.file_object.read(section, "%s(%i)" % (item, istep + 1))
        for item in self.mdblockitems:
            block = int(istep / self.mdblocksize)
            pos = istep % self.mdblocksize
            # First check if this block was already stored in memory
            if item in self._mdblock:
                if block + 1 in self._mdblock[item]:
                    values = self._mdblock[item][block + 1]
                    self.mddata[item] = values[pos]
                    continue
            # If not, try to read the block
            if not (section, "%s(%i)" % (item, block + 1)) in self.file_object:
                continue
            values = self.file_object.read(section, "%s(%i)" % (item, block + 1))
            if isinstance(values, str):
                values = values.split()
            if not isinstance(values, list):
                values = [values]
            self._mdblock[item] = {block + 1: values}
            self.mddata[item] = values[pos]

    def _store_historydata_for_step(self, istep):
        """
        Store the extra data from the History section

        Note: Block format is not used in the History section
        """
        if "w" in self.mode:
            return
        if self.historydata is None:
            self.historydata = {}
        section = "History"
        for item in self.historyitems:
            if (section, "%s(%i)" % (item, istep + 1)) in self.file_object:
                self.historydata[item] = self.file_object.read(section, "%s(%i)" % (item, istep + 1))

    def _is_endoffile(self):
        """
        Reads and checks If the end of file is reached.
        """
        return ("History", "Coords(%i)" % (self.position + 1)) in self.file_object

    def read_next(self, molecule=None, read=True):
        """
        Reads coordinates and lattice vectors from the current position of the cursor and returns it

        * ``molecule`` -- |Molecule| object in which the new coordinates need to be stored
        * ``read``     -- If set to False the cursor will move to the next frame without reading
        """
        if not read and not self.firsttime:
            return self._move_cursor_without_reading()

        if self.firsttime:
            self.firsttime = False

        crd, vecs = self.read_frame(self.position, molecule)
        self.position += 1
        return crd, vecs

    def write_next(self, coords=None, molecule=None, cell=[0.0, 0.0, 0.0], conect=None, historydata=None, mddata=None):
        """
        Write frame to next position in trajectory file

        * ``coords``   -- A list or numpy array of (``ntap``,3) containing the system coordinates in angstrom
        * ``molecule`` -- A molecule object to read the molecular data from
        * ``cell``     -- A set of lattice vectors (or cell diameters for an orthorhombic system) in angstrom
        * ``conect``   -- A dictionary containing the connectivity info (e.g. {1:[2],2:[1]})
        * ``historydata`` -- A dictionary containing additional variables to be written to the History section
        * ``mddata``   -- A dictionary containing the variables to be written to the MDHistory section

        The ``mddata`` dictionary can contain the following keys:
        ('TotalEnergy', 'PotentialEnergy', 'Step', 'Velocities', 'KineticEnergy',
        'Charges', 'ConservedEnergy', 'Time', 'Temperature')

        The ``historydata`` dictionary can contain for example:
        ('Energy','Gradients','StressTensor')
        All values must be in atomic units
        Numpy arrays or lists of lists will be flattened before they are written to the file

        .. note::

                Either ``coords`` or ``molecule`` are mandatory arguments
        """
        # Check for common error in the arguments
        if coords is not None:
            if isinstance(coords, Molecule):
                raise PlamsError("The PLAMS molecule needs to be passed as the second argument (molecule)")

        if isinstance(molecule, Molecule):
            coords, cell, elements, conect, _ = self._read_plamsmol(molecule, read_props=False)
            if self.position == 0:
                self.elements = elements
        # Make sure that the cell consists of vectors
        cell = self._convert_cell(cell)
        if conect is not None:
            if len(conect) == 0:
                conect = None
        self.conect = conect

        # Include a check on the size of coords?
        if len(coords) != len(self.elements):
            raise PlamsError("The coordinates do not match the rest of the trajectory")

        # If this is the first step, write the header
        if self.position == 0:
            self._write_header(coords, cell, molecule)
            self.firsttime = False

        # Define some local variables
        step = self.position
        if mddata is not None:
            if "Step" in mddata:
                step = mddata["Step"]
        # Energy should be read from mddata first, otherwise from historydata, otherwise set to zero
        energy = self._set_energy(mddata, historydata)
        if not self.include_historydata or historydata is None:
            historydata = {}
        historydata["Energy"] = energy

        # Write the history section
        counter = 1
        counter = self._write_history_entry(step, coords, cell, conect, historydata, counter)

        if self.include_mddata and mddata is not None:
            self._write_mdhistory_entry(mddata)

        self.position += 1

        if self.saving_freq is not None:
            if self.position % self.saving_freq == 0:
                self.file_object.save()

    def _set_energy(self, mddata, historydata):
        """
        Looks if an energy is passed as input, and it not, sets to zero
        """
        energy = None
        if mddata is not None:
            if "PotentialEnergy" in mddata:
                energy = mddata["PotentialEnergy"]
        if energy is None:
            if historydata is not None:
                if "Energy" in historydata:
                    energy = historydata["Energy"]
        if energy is None:
            energy = 0.0
        return energy

    def _write_history_entry(self, step, coords, cell, conect, historydata=None, counter=1):
        """
        Write the full entry into the History section
        """
        self.file_object.write("History", "nEntries", self.position + 1)
        self.file_object.write("History", "currentEntryOpen", False)
        self._write_keydata_in_history("Step", counter, False, 1, self.position + 1, step)
        counter += 1
        crd = [float(c) / bohr_to_angstrom for coord in coords for c in coord]
        self._write_keydata_in_history("Coords", counter, True, 3, self.position + 1, crd)
        counter += 1
        # self._write_keydata_in_history('Energy', counter, False, 1, self.position+1, energy)
        # counter += 1
        if cell is not None:
            self._write_keydata_in_history("nLatticeVectors", counter, False, 1, self.position + 1, self.nvecs)
            counter += 1
            vecs = [float(v) / bohr_to_angstrom for vec in cell for v in vec]
            # I should probably rethink the dimension of the lattice vectors (generalize it)
            self._write_keydata_in_history("LatticeVectors", counter, False, [3, 3], self.position + 1, vecs)
            counter += 1

        if historydata is not None:
            counter = self._write_dictionary_to_history(historydata, "History", counter)
        # if gradients is not None :
        #        grd = [float(g) for grad in gradients for g in grad]
        #        self._write_keydata_in_history('Gradients', counter, True, 3, self.position+1, grd)

        # if stresstensor is not None :
        #        stre = [float(s) for stress in stresstensor for s in stress]
        #        self._write_keydata_in_history('StressTensor', counter, False, [3,3], self.position+1, stre)

        # Write the bond info
        if conect is not None:
            counter = self._write_bonds_in_history(conect, counter, len(coords))

        return counter

    def _write_bonds_in_history(self, conect, counter, nats):
        """
        Write the bond data into the history section
        """
        # Get the bond orders out of the connection table
        connections = {}
        orders = {}
        for k in conect.keys():
            connections[k] = [t[0] if isinstance(t, tuple) else t for t in conect[k]]
            orders[k] = [t[1] if isinstance(t, tuple) else 1.0 for t in conect[k]]

        numbonds = 0
        indices = [1]
        connection_table = []
        bond_orders = []
        for iat in range(1, nats + 1):
            connections = []
            if iat in conect:
                connections = conect[iat]
            neighbors = []
            bos = []
            for t in connections:
                jat = t
                bo = 1.0
                if isinstance(t, tuple):
                    jat = t[0]
                    bo = t[1]
                # Correct for double counting
                if jat <= iat:
                    continue
                neighbors.append(jat)
                bos.append(bo)
            numbonds += len(bos)
            indices.append(numbonds + 1)
            connection_table += neighbors
            bond_orders += bos

        self.file_object.write("History", "Bonds.Index(%i)" % (self.position + 1), indices)
        self.file_object.write("History", "ItemName(%i)" % (counter), "%s" % ("Bonds.Index"))
        counter += 1

        self.file_object.write("History", "Bonds.Atoms(%i)" % (self.position + 1), connection_table)
        self.file_object.write("History", "ItemName(%i)" % (counter), "%s" % ("Bonds.Atoms"))
        counter += 1

        self.file_object.write("History", "Bonds.Orders(%i)" % (self.position + 1), bond_orders)
        self.file_object.write("History", "ItemName(%i)" % (counter), "%s" % ("Bonds.Orders"))
        counter += 1

        return counter

    def _write_mdhistory_entry(self, mddata):
        """
        Write the entry in the MDHistory section
        """
        counter = 1
        counter = self._write_dictionary_to_history(mddata, self.mdhistory_name, counter)

    def _write_dictionary_to_history(self, data, section, counter=1):
        """
        Add the entries of a dictionary to a History section
        """
        self.file_object.write(section, "nEntries", self.position + 1)
        self.file_object.write(section, "currentEntryOpen", False)
        for key, var in data.items():
            # Make sure that the entry is either a scalar or a 1D list
            var = self._flatten_variable(var)
            peratom = False
            dim = 1
            if isinstance(var, list):
                if len(var) % len(self.elements) == 0:
                    dim = int(len(var) / len(self.elements))
                    peratom = True
                else:
                    dim = len(var)
            self._write_keydata_in_history(key, counter, peratom, dim, self.position + 1, var, section)
            counter += 1
        return counter

    def _flatten_variable(self, var):
        """
        Make sure that the variable is a Python 1D list (not numpy)
        """
        while True:
            if isinstance(var, list) or isinstance(var, numpy.ndarray):
                if len(var) == 0:
                    break
                if isinstance(var[0], list) or isinstance(var[0], numpy.ndarray):
                    var = [v for varitem in var for v in varitem]
                else:
                    if isinstance(var[0], numpy.int64):
                        var = [int(v) for v in var]
                    elif isinstance(var[0], numpy.float64):
                        var = [float(v) for v in var]
                    break
            else:
                if isinstance(var, numpy.int64):
                    var = int(var)
                elif isinstance(var, numpy.float64):
                    var = float(var)
                break
        return var

    def _write_keydata_in_history(self, key, i, perAtom, dim, step, values, section="History"):
        """
        Write all data about a key value in KFFile
        """
        # Some data only needs to be printed once
        printstartdata = False
        if step == 1:
            printstartdata = True

        # Block code: if the data is to be written as blocks, then step and values need to be replaced.
        if section == self.mdhistory_name:
            step, values = self._get_block_info(key, perAtom, dim, step, values, section)

        # The rest should be independent on format (block or individual)
        self.file_object.write(section, "%s(%i)" % (key, step), values)
        if printstartdata:
            self.file_object.write(section, "ItemName(%i)" % (i), "%s" % (key))
            self.file_object.write(section, "%s(perAtom)" % (key), perAtom)
            self.file_object.write(section, "%s(dim)" % (key), dim)
            if section == self.mdhistory_name and self.mdunits is not None:
                if key in self.mdunits:
                    self.file_object.write(section, "%s(units)" % (key), self.mdunits[key])

    def _get_block_info(self, key, perAtom, dim, step, values, section):
        """
        If the data is to be written as blocks, then step and values need to be replaced.
        """
        if dim == 1 and not perAtom and not (isinstance(values, list) or isinstance(values, str)):
            if not self.include_mddata:
                raise Exception("Set include_mddata to write the MD section.")
            iblock = int((step - 1) / self.mdblocksize) + 1
            if step % self.mdblocksize != 1:
                old_values = []
                if key in self._mdblock:
                    if iblock in self._mdblock[key]:
                        if len(self._mdblock[key][iblock]) == step % self.mdblocksize - 1:
                            old_values = self._mdblock[key][iblock]
                if len(old_values) == 0:
                    if (section, "%s(%i)" % (key, iblock)) in self.file_object:
                        old_values = self.file_object.read(section, "%s(%i)" % (key, iblock))
                        if not isinstance(old_values, list):
                            old_values = [old_values]
                values = old_values + [values]  # Values is a scalar
            else:
                self.file_object.write(section, "nBlocks", iblock)
            step = iblock
            self._mdblock[key] = {iblock: values}
            if not isinstance(values, list):
                self._mdblock[key] = {iblock: [values]}
        return step, values

    def rewind(self, nframes=None):
        """
        Rewind the file either by ``nframes`` or to the first frame

        *   ``nframes`` -- The number of frames to rewind
        """
        self.firsttime = True
        self.position = 0

    def get_length(self):
        """
        Get the number of frames in the file
        """
        nsteps = 0
        if "History" in self.file_object:
            nsteps = self.file_object.read("History", "nEntries")
        return nsteps

    def read_last_frame(self, molecule=None):
        """
        Reads the last frame from the file
        """
        nsteps = self.get_length()
        crd, cell = self.read_frame(nsteps - 1, molecule)
        return crd, cell


def write_general_section(rkf, program="plams"):
    """
    Write the General section of the RKF file
    """
    rkf.write("General", "file-ident", "RKF")
    rkf.write("General", "termination status", "NORMAL TERMINATION")
    rkf.write("General", "program", "%s" % (program))
    rkf.write("General", "user input", " ")


def write_molecule_section(rkf, coords=None, cell=None, elements=None, section="Molecule", molecule=None):
    """
    Write the molecule section

    Note: Currently does not write bonds
    """
    if molecule is not None:
        if coords is None:
            coords = molecule.as_array()
        if cell is None and len(molecule.lattice) > 0:
            cell = molecule.lattice
        if elements is None:
            elements = [at.symbol for at in molecule.atoms]

    # Then write the input molecule
    charge = 0.0
    if molecule is not None:
        if "charge" in molecule.properties:
            charge = float(molecule.properties.charge)
    element_numbers = [PeriodicTable.get_atomic_number(el) for el in elements]

    rkf.write(section, "nAtoms", len(elements))
    rkf.write(section, "AtomicNumbers", element_numbers)
    rkf.write(section, "AtomSymbols", elements)
    crd = [Units.convert(float(c), "angstrom", "bohr") for coord in coords for c in coord]
    rkf.write(section, "Coords", crd)
    rkf.write(section, "Charge", charge)
    if cell is not None:
        rkf.write(section, "nLatticeVectors", len(cell))
        vecs = [Units.convert(float(v), "angstrom", "bohr") for vec in cell for v in vec]
        rkf.write(section, "LatticeVectors", vecs)
    # Should it write bonds?
    # Write atom properties
    if molecule is not None:
        from scm.plams.interfaces.adfsuite.ams import AMSJob

        suffixes = [AMSJob._atom_suffix(at) for at in molecule]
        if any(s != "" for s in suffixes):
            rkf.write(section, "EngineAtomicInfo", "\x00".join(suffixes))
        # Add atomic charges
        charges = [at.properties.forcefield for at in molecule.atoms if "forcefield" in at.properties.keys()]
        charges = [float(s.charge) for s in charges if "charge" in s.keys()]
        if len(charges) == len(molecule):
            rkf.write(section, "Charges", charges)
        # Add region sections
        if "regions" in molecule.properties.keys():
            region_names = molecule.properties.regions.keys()
            rkf.write(section, "RegionNames", "\x00".join(region_names))
            rkf.write(
                section,
                "RegionProperties",
                "\x00".join(["\n".join(molecule.properties.regions[k]) for k in region_names]),
            )
        # Also add a bond section
        if len(molecule.bonds) > 0:
            bond_indices = [sorted([iat for iat in molecule.index(bond)]) for bond in molecule.bonds]
            atoms_from = [bond[0] for bond in bond_indices]
            atoms_to = [bond[1] for bond in bond_indices]
            orders = [float(bond.order) for bond in molecule.bonds]
            rkf.write(section, "fromAtoms", atoms_from)
            rkf.write(section, "toAtoms", atoms_to)
            rkf.write(section, "bondOrders", orders)

            # I also need to write the lattice displacements of the bonds
            # To do that, I need to compute them.
            # I could make a start with writing in zeros
            if cell is not None:
                lattice_displacements = compute_lattice_displacements(molecule)
                # lattice_displacements = [0 for i in range(len(cell)*len(molecule.bonds))]
                rkf.write(section, "latticeDisplacements", lattice_displacements)


def compute_lattice_displacements(molecule):
    """
    Determine which bonds are displaced along the periodic lattice, so that they are not at their closest distance
    """
    cell = numpy.array(molecule.lattice)
    nvecs = len(cell)

    # Get the difference vectors for the bonds
    nbonds = len(molecule.bonds)
    bond_indices = numpy.array([sorted([iat - 1 for iat in molecule.index(bond)]) for bond in molecule.bonds])
    coords = molecule.as_array()
    vectors = coords[bond_indices[:, 0]] - coords[bond_indices[:, 1]]

    # Project the vectors onto the lattice vectors
    celldiameters_sqrd = (cell**2).sum(axis=1)
    proj = (vectors.reshape((nbonds, 1, 3)) * cell.reshape(1, nvecs, 3)).sum(axis=2)
    proj = proj / celldiameters_sqrd.reshape((1, nvecs))

    # Now see what multiple they are of 0.5
    lattice_displacements = numpy.round(proj).astype(int)
    lattice_displacements = lattice_displacements.reshape((nvecs * nbonds)).tolist()
    return lattice_displacements

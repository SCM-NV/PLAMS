#!/usr/bin/env python

import os
from typing import List, Set, Union

from scm.plams.core.errors import PlamsError
from scm.plams.core.settings import Settings
import numpy
from scm.plams.mol.atom import Atom
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.periodic_table import PT
from scm.plams.trajectories.rkffile import RKFTrajectoryFile, bohr_to_angstrom, write_general_section

__all__ = ["RKFHistoryFile", "molecules_to_rkf", "rkf_filter_regions"]


class RKFHistoryFile(RKFTrajectoryFile):
    """
    Class representing an RKF file containing a molecular simulation history with varying numbers of atoms

    An instance of this class has the following attributes:

    *   ``file_object`` -- A PLAMS |KFFile| object, referring to the actual RKF file
    *   ``position``    -- The frame to which the cursor is currently pointing in the RKF file
    *   ``mode``        -- Designates whether the file is in read or write mode ('rb' or 'wb')
    *   ``elements``    -- The elements of the atoms in the system at the current frame
    *   ``conect``      -- The connectivity information of the current frame
    *   ``mddata``      -- Read mode only: A dictionary containing data from the MDHistory section in the RKF file
    *   ``read_lattice``-- Read mode only: Wether the lattice vectors will be read from the file
    *   ``read_bonds``  -- Wether the connectivity information will be read from the file

    An |RKFHistoryFile| object behaves very similar to a regular file object.
    It has read and write methods (:meth:`read_next` and :meth:`write_next`)
    that read and write from/to the position of the cursor in the ``file_object`` attribute.
    If the file is in read mode, an additional method :meth:`read_frame` can be used that moves
    the cursor to any frame in the file and reads from there.
    The amount of information stored in memory is kept to a minimum, as only information from the latest frame
    is ever stored.

    Reading and writing to and from the files can be done as follows::

        >>> from scm.plams import RKFHistoryFile

        >>> rkf = RKFHistoryFile('ams.rkf')
        >>> mol = rkf.get_plamsmol()

        >>> rkfout = RKFHistoryFile('new.rkf',mode='wb')

        >>> for i in range(rkf.get_length()) :
        >>>     crd,cell = rkf.read_frame(i,molecule=mol)
        >>>     rkfout.write_next(molecule=mol)
        >>> rkfout.close()

    The above script reads information from the RKF file ``ams.rkf`` into the |Molecule| object ``mol``
    in a step-by-step manner..
    The |Molecule| object is then passed to the :meth:`write_next` method of the new |RKFHistoryFile|
    object corresponding to the new rkf file ``new.rkf``.

    The exact same result can also be achieved by iterating over the instance as a callable

        >>> rkf = RKFHistoryFile('ams.rkf')
        >>> mol = rkf.get_plamsmol()

        >>> rkfout = RKFHistoryFile('new.rkf',mode='wb')

        >>> for crd,cell in rkf(mol) :
        >>>     rkfout.write_next(molecule=mol)
        >>> rkfout.close()

    This procedure requires all coordinate information to be passed to and from the |Molecule| object
    for each frame, which can be time-consuming.
    Some time can be saved by bypassing the |Molecule| object::

        >>> rkf = RKFHistoryFile('ams.rkf')

        >>> rkfout = RKFHistoryFile('new.rkf',mode='wb')

        >>> for crd,cell in rkf :
        >>>     rkfout.write_next(coords=crd,cell=cell,elements=rkf.elements,conect=rkf.conect)
        >>> rkfout.close()

    The only mandatory argument to the :meth:`write_next` method is ``coords``.
    Further time can be saved by setting the ``read_lattice`` and ``read_bonds`` variables to False.

    By default the write mode will create a minimal version of the RKF file, containing only elements,
    coordinates, lattice, and connectivity information.
    This minimal file format can be read by AMSMovie.

    If the original RKF file contains an MDHistory section (if it resulted from a MolecularGun simulation)
    it is possible to store the information from that section and write it to another file.
    To enable this, the method :meth:`store_mddata` needs to be called after creation,
    and a dictionary of mddata needs to be passed to the :meth:`write_next` method.
    When that is done, the AMS trajectory analysis tools can be used on the file.
    Restarting an MD run with such a file is however currently not possible::

        >>> rkf = RKFHistoryFile('ams.rkf')
        >>> rkf.store_mddata()
        >>> mol = rkf.get_plamsmol()

        >>> rkf_out = RKFHistoryFile('new.rkf',mode='wb')
        >>> rkf_out.store_mddata(rkf)

        >>> for i in range(rkf.get_length()) :
        >>>         crd,cell = rkf.read_frame(i,molecule=mol)
        >>>         rkf_out.write_next(molecule=mol,mddata=rkf.mddata)
        >>> rkf_out.close()
    """

    def __init__(self, filename, mode="rb", fileobject=None, ntap=None):
        """
        Initializes the RKFHistoryFile object

        * ``filename``   -- The path to the RKF file
        * ``mode``       -- The mode in which to open the RKF file ('rb' or 'wb')
        * ``fileobject`` -- Optionally, a file object can be passed instead (filename needs to be set to None)
        * ``ntap``       -- If the file is in write mode, the number of atoms can be passed here
        """
        self.added_atoms = None
        self.removed_atoms = None
        self.chemical_systems = None
        RKFTrajectoryFile.__init__(self, filename, mode, fileobject, ntap)

        self.input_elements = self.elements[:]
        self.versionhistory_length = 0
        self.system_version_elements = {}
        self.frame = -1
        self.props = None
        self.charge = None

        # This is for writing only
        self.system_version_props = {}
        self.system_version_charge = {}
        self.version_history_items = []

    def set_elements(self, elements):
        """
        Sets the elements attribute (needed in write mode).

        *   ``elements`` -- A list containing the element name of each atom
        """
        if self.position > 0:
            raise PlamsError("Elements should not be changed while reading/writing is already in progress")
        RKFTrajectoryFile.set_elements(self, elements)
        self.input_elements = elements

    def get_plamsmol(self):
        """
        Extracts a PLAMS molecule object from the RKF file
        """
        section_dict = self.file_object.read_section("ChemicalSystem(1)")
        if len(section_dict) == 0:
            section_dict = self.file_object.read_section("InputMolecule")
        if len(section_dict) == 0:
            section_dict = self.file_object.read_section("Molecule")
        plamsmol = Molecule._mol_from_rkf_section(section_dict)
        return plamsmol

    def _rewrite_molecule(self):
        """
        Overwrite the molecule section with the latest frame (called in close())
        """
        mol = self.get_plamsmol()
        crd, cell = self.read_last_frame(molecule=mol)
        self._write_molecule_section(crd, cell, molecule=mol)

    def _read_header(self):
        """
        Read the start molecule data from the InputMolecule section (not the Molecule section)
        """
        self.added_atoms = {}
        self.removed_atoms = {}
        self.chemical_systems = {0: 1}
        secname = "ChemicalSystem(1)"
        if not "SystemVersionHistory" in self.file_object.sections():
            if "InputMolecule" in self.file_object:
                secname = "InputMolecule"
            elif "Molecule" in self.file_object:
                secname = "Molecule"
        if not secname in self.file_object:
            return
        RKFTrajectoryFile._read_header(self, molecule_section=secname)

        # Now store the added and removed atoms along the trajectory
        # (This might be slow?)
        # I could also do it on the fly, but that may be messy when we move back and forth through the file
        if not "SystemVersionHistory" in self.file_object.sections():
            return
        # version = 0
        version = 1
        self._set_system_version_elements()
        for i in range(self.get_length()):
            if not ("History", "SystemVersion(%i)" % (i + 1)) in self.file_object:
                continue
            new_version = self.file_object.read("History", "SystemVersion(%i)" % (i + 1))
            if new_version == version:
                continue
            self.added_atoms[i] = {}
            self.removed_atoms[i] = {}
            # Now look for the added and removed atoms
            removed_atoms = []
            if ("SystemVersionHistory", "RemovedAtoms(%i)" % (new_version)) in self.file_object:
                removed_atoms = self.file_object.read("SystemVersionHistory", "RemovedAtoms(%i)" % (new_version))
            if not isinstance(removed_atoms, list):
                removed_atoms = [removed_atoms]
            added_atoms = []
            if ("SystemVersionHistory", "AddedAtoms(%i)" % (new_version)) in self.file_object:
                added_atoms = self.file_object.read("SystemVersionHistory", "AddedAtoms(%i)" % (new_version))
            if not isinstance(added_atoms, list):
                added_atoms = [added_atoms]
            # Now find the corresponding elements
            chemSysNum = self.file_object.read("SystemVersionHistory", "SectionNum(%i)" % (new_version))
            # Compare to the previous chemical system
            # prev_version = new_version - 1
            # prevChemSysNum = self.file_object.read('SystemVersionHistory','SectionNum(%i)'%(prev_version))
            prevChemSysNum = self.chemical_systems[max(self.chemical_systems.keys())]
            sectionname = "ChemicalSystem(%i)" % (prevChemSysNum)
            atnums = self.file_object.read(sectionname, "AtomicNumbers")
            if isinstance(atnums, int):
                atnums = [atnums]
            prev_elements = [PT.get_symbol(atnum) for atnum in atnums]
            sectionname = "ChemicalSystem(%i)" % (chemSysNum)
            atnums = self.file_object.read(sectionname, "AtomicNumbers")
            if isinstance(atnums, int):
                atnums = [atnums]
            elements = [PT.get_symbol(atnum) for atnum in atnums]
            ################
            # Read badly written files
            diff = (len(added_atoms) - len(removed_atoms)) - (len(elements) - len(prev_elements))
            if diff != 0:
                chemSysNum, elements = self._correct_chemical_system(
                    elements, prev_elements, added_atoms, removed_atoms
                )
            # print ('chemSysNum: ',i,chemSysNum, len(elements),added_atoms,removed_atoms)
            ####################
            removed_elements = [prev_elements[i - 1].split(".")[0] for i in removed_atoms]
            added_elements = [elements[i - 1].split(".")[0] for i in added_atoms]
            # Now store the elements
            for iat, el in zip(removed_atoms, removed_elements):
                self.removed_atoms[i][iat] = el
            for iat, el in zip(added_atoms, added_elements):
                self.added_atoms[i][iat] = el
            self.chemical_systems[i] = chemSysNum
            version = new_version

    def _correct_chemical_system(self, elements, prev_elements, added_atoms, removed_atoms):
        """
        Check if the referenced chemical system is correct, and if not, find one matching added/removed atoms
        """
        # print ('Searching for chemical system elsewhere')
        new_elements = [el for el in prev_elements]
        for ra in reversed(removed_atoms):
            new_elements.pop(ra - 1)
        for aa in added_atoms:
            new_elements.insert(aa - 1, "*")  # Add a wildcard
        chemSysNum = self._check_for_chemical_system(new_elements)
        if chemSysNum is None:
            print("Chemical system does not exist")
            chemSysNum = self._check_for_chemical_system(new_elements, compare_elements=False)
            if chemSysNum is None:
                raise Exception("Chemical system with this number of atoms does not exist")
        sectionname = "ChemicalSystem(%i)" % (chemSysNum)
        elements = [PT.get_symbol(atnum) for atnum in self.file_object.read(sectionname, "AtomicNumbers")]
        ###################
        # Check the regions
        section_dict = self.file_object.read_section(sectionname)
        plamsmol = Molecule._mol_from_rkf_section(section_dict)
        print("Checking the regions: ", 0, [k for k in plamsmol.atoms[0].properties])
        for aa in added_atoms:
            print(aa - 1, [k for k in plamsmol.atoms[aa - 1].properties])
        ###################
        return chemSysNum, elements

    # FIXME: The _write_header section writes the starting molecule to the Molecule section,
    #        not the final molecule (like the Fortran code does)
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
        # I think the InputMolecule is mandatory in this case
        self._write_molecule_section(coords, cell, section="InputMolecule", molecule=molecule)
        self.chemical_systems = {}
        if self.include_mddata:
            # Start setting up the MDHistory section as well
            self.mdblocksize = 100
            self.file_object.write(self.mdhistory_name, "blockSize", 100)

        self.added_atoms = {}
        self.removed_atoms = {}

    def _read_coordinates(self, i, molecule, cell):
        """
        Read the coordinates at step i
        """
        coords = self.coords.reshape(len(self.coords) * 3)

        update_molecule = False
        elements = self._read_elements_for_frame(i)
        if elements != self.elements:
            update_molecule = True
        if isinstance(molecule, Molecule):
            _, _, mol_elements, _, props = self._read_plamsmol(molecule)
            charge = None
            if "charge" in molecule.properties:
                charge = molecule.properties.charge
            if (
                mol_elements != self.elements or props != self.props or charge != self.charge
            ):  # This molecule has nothing to do with the previously read one
                update_molecule = True
            prev_frames = [iframe for iframe in self.chemical_systems.keys() if iframe <= i]
            ifr = 0
            if len(prev_frames) > 0:
                ifr = prev_frames[-1]
            next_frames = [iframe for iframe in self.chemical_systems.keys() if iframe > i]
            jfr = len(self)
            if len(next_frames) > 0:
                jfr = next_frames[0]
            # If the previous frame falls before the last chemical system or after the next one, update
            if self.frame < ifr or self.frame >= jfr:
                update_molecule = True

        if update_molecule:
            self.elements = elements
            self.frame = i
            self.coords = numpy.zeros((len(elements), 3))
            coords = self.coords.reshape((len(elements) * 3))
            # Rebuild the molecule (bonds will disappear for now)
            if isinstance(molecule, Molecule):
                # self.props = props
                secname = "ChemicalSystem(%i)" % (self.chemical_systems[ifr])
                if not "SystemVersionHistory" in self.file_object.sections():
                    secname = "InputMolecule"
                section_dict = self.file_object.read_section(secname)
                new_mol = Molecule._mol_from_rkf_section(section_dict)
                for at in reversed(molecule.atoms):
                    molecule.delete_atom(at)
                molecule.properties = Settings()
                for iel, el in enumerate(elements):
                    atom = new_mol.atoms[iel]
                    atom.bonds = []
                    # atom = Atom(PT.get_atomic_number(el))
                    molecule.add_atom(atom)
                _, _, _, _, self.props = self._read_plamsmol(molecule)
                if "charge" in new_mol.properties:
                    molecule.properties.charge = new_mol.properties.charge
                    self.charge = molecule.properties.charge
                # Now what if the elements found in this ChemicalSystem section
                # do not match the expected elements (self.elements)?
                new_elements = [at.symbol for at in molecule.atoms]
                if new_elements != self.elements:
                    print("Chemical section does not match", i)
                    atnums = [PT.get_atomic_number(el) for el in self.elements]
                    for iat, atnum in enumerate(atnums):
                        molecule.atoms[iat].atnum = atnum

        coords[:] = self.file_object.read("History", "Coords(%i)" % (i + 1))
        coords *= bohr_to_angstrom
        # This changes self.coords behind the scenes

        # Assign the data to the molecule object
        if isinstance(molecule, Molecule):
            cell_reduced = None
            if cell is not None:
                cell_reduced = cell[: self.nvecs]
            self._set_plamsmol(self.coords, cell_reduced, molecule)

    def _read_elements_for_frame(self, frame):
        """
        Use the added and removed atoms to read the elements at each frame
        """
        elements = self.input_elements[:]
        for i in range(frame + 1):
            if i in self.removed_atoms:
                # Insert them into the elements list, in the correct place
                elements = [el for iat, el in enumerate(elements) if not iat + 1 in self.removed_atoms[i]]
            if i in self.added_atoms:
                for iat, el in self.added_atoms[i].items():
                    elements.insert(iat - 1, el)
        return elements

    def write_next(
        self,
        coords=None,
        molecule=None,
        elements=None,
        cell=[0.0, 0.0, 0.0],
        conect=None,
        historydata=None,
        mddata=None,
    ):
        """
        Write frame to next position in trajectory file

        * ``coords``   -- A list or numpy array of (``ntap``,3) containing the system coordinates
        * ``molecule`` -- A molecule object to read the molecular data from
        * ``elements`` -- The element symbols of the atoms in the system
        * ``cell``     -- A set of lattice vectors (or cell diameters for an orthorhombic system)
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

                Either ``coords`` and ``elements`` or ``molecule`` are mandatory arguments
        """
        # Check for common error in the arguments
        if coords is not None:
            if isinstance(coords, Molecule):
                raise PlamsError("The PLAMS molecule needs to be passed as the second argument (molecule)")

        props = None
        charge = None
        if isinstance(molecule, Molecule):
            coords, cell, elements, conect, props = self._read_plamsmol(molecule)
            if "charge" in molecule.properties:
                charge = molecule.properties.charge
        # Make sure that the cell consists of three vectors
        cell = self._convert_cell(cell)
        if conect is not None:
            if len(conect) == 0:
                conect = None
        self.conect = conect

        # If this is the first step, write the header
        if self.position == 0:
            self.elements = elements
            self.props = props
            self.charge = charge
            self._write_header(coords, cell, molecule=molecule)
            # This is specific for the history-file object
            self.input_elements = elements[:]

        # Define some local variables
        counter = 1
        step = self.position
        if mddata is not None:
            if "Step" in mddata:
                step = mddata["Step"]
        # Energy should be read from mddata first, otherwise from historydata, otherwise set to zero
        energy = self._set_energy(mddata, historydata)
        if historydata is None or not self.include_historydata:
            historydata = {}
        historydata["Energy"] = energy

        # Write the history section
        counter = self._write_history_entry(step, coords, cell, conect, historydata)

        if self.include_mddata and mddata is not None:
            self._write_mdhistory_entry(mddata)

        # If a change took place, write it.
        if elements != self.elements or props != self.props or self.position == 0 or charge != self.charge:
            # version history has to be checked before the elements can be updated
            chemsysversion = self._check_for_chemical_system(elements, props=props, charge=charge)
            self._write_version_history(elements, coords, cell, props, charge)
            # Now update the elements and props
            self.elements = elements
            self.props = props
            # Write the molecule sections
            if chemsysversion is None:
                chemsysversion = len(self.system_version_elements)
                self._write_molecule_section(
                    coords, cell, section="ChemicalSystem(%i)" % (chemsysversion), molecule=molecule
                )
            self.chemical_systems[self.position] = chemsysversion
        counter = self._write_system_version_history_entry(counter)

        self.position += 1

        if self.saving_freq is not None:
            if self.position % self.saving_freq == 0:
                self.file_object.save()

    def _write_version_history(self, elements, coords, cell, props, charge):
        """
        Write the version history
        """
        # Enter the correct SystemVersionHistory entry into the History section
        self.versionhistory_length += 1
        version = self.versionhistory_length

        # Find the corresponding ChemicalSystem and write it
        chemsysversion = self._check_for_chemical_system(elements, props=props, charge=charge)
        if chemsysversion is None:
            # self.system_version_elements.append(elements[:])
            if len(self.system_version_elements) in self.system_version_elements.keys():
                raise Exception("self.system_version_elements is not consecutively numbered")
            self.system_version_props[len(self.system_version_elements)] = props
            self.system_version_charge[len(self.system_version_elements)] = charge
            self.system_version_elements[len(self.system_version_elements)] = elements[:]
            chemsysversion = len(self.system_version_elements)

        # Now add an entry to SystemVersionHistory
        added_atoms, removed_atoms = self._find_system_change(elements, props)
        self.file_object.write("SystemVersionHistory", "nEntries", version)
        self.file_object.write("SystemVersionHistory", "currentEntryOpen", False)
        if "SectionNum" not in self.version_history_items:
            self.version_history_items.append("SectionNum")
        svh_counter = self.version_history_items.index("SectionNum") + 1
        self._write_keydata_in_history(
            "SectionNum", svh_counter, False, 1, version, chemsysversion, "SystemVersionHistory"
        )
        # Avoid writing empty data into this section by storing the items.
        if len(added_atoms) > 0:
            if "AddedAtoms" not in self.version_history_items:
                self.version_history_items.append("AddedAtoms")
            svh_counter = self.version_history_items.index("AddedAtoms") + 1
            data = [iat + 1 for iat in added_atoms]
            self._write_keydata_in_history(
                "AddedAtoms", svh_counter, False, len(data), version, data, "SystemVersionHistory"
            )
        if len(removed_atoms) > 0:
            if "RemovedAtoms" not in self.version_history_items:
                self.version_history_items.append("RemovedAtoms")
            svh_counter = self.version_history_items.index("RemovedAtoms") + 1
            data = [iat + 1 for iat in removed_atoms]
            self._write_keydata_in_history(
                "RemovedAtoms", svh_counter, False, len(data), version, data, "SystemVersionHistory"
            )

    def _write_system_version_history_entry(self, counter):
        """
        Write the entry for the SystemVersionHistory into the History section
        """
        version = self.versionhistory_length
        self._write_keydata_in_history("SystemVersion", counter, False, 1, self.position + 1, version)
        counter += 1
        return counter

    def _find_system_change(self, elements, props):
        """
        Find out which atoms were added and/or deleted
        """
        # First find out which elements were removed
        removed_atoms = []
        position = 0
        for i, el in enumerate(self.elements):
            if position == len(elements):
                removed_atoms += [j for j in range(i, len(self.elements))]
                break
            # Get the props, if relevant
            if self.props is None:
                op = None
            else:
                op = self.props[i]
            if props is None:
                p = None
            else:
                p = props[position]
            # Find the changes
            if el == elements[position] and op == p:
                position += 1
            else:
                removed_atoms.append(i)

        # Then find out which elements were added
        added_atoms = [i for i in range(position, len(elements))]

        # Now store them (this is actually not really necessary)
        if len(removed_atoms) > 0:
            self.removed_atoms[self.position] = {}
        if len(added_atoms) > 0:
            self.added_atoms[self.position] = {}
        for iat in removed_atoms:
            self.removed_atoms[self.position][iat + 1] = self.elements[iat]
        for iat in added_atoms:
            self.added_atoms[self.position][iat + 1] = elements[iat]
        return added_atoms, removed_atoms

    def _check_for_chemical_system(self, elements, compare_elements=True, props=None, charge=None):
        """
        Check if the new chemical system was encountered before
        """
        iterator = self.system_version_elements.items()

        version = None
        for ind, (i, prev_elements) in enumerate(iterator):
            equal = True
            if len(prev_elements) != len(elements):
                continue
            if compare_elements:
                for el, pel in zip(elements, prev_elements):
                    if el != pel:
                        if el == "*":
                            continue
                        equal = False
                        break
            if props is not None and equal:
                # This only occurs when writing, otherwise system_version_props is empty
                if i in self.system_version_props.keys():
                    if props != self.system_version_props[i]:
                        equal = False
            if charge is not None and equal:
                if i in self.system_version_charge:
                    if charge != self.system_version_charge[i]:
                        equal = False
            if equal:
                version = i + 1
        return version

    def _set_system_version_elements(self):
        """
        Store all chemical systems from the file
        """
        self.system_version_elements = {}
        keys = [key for key in self.file_object.get_skeleton().keys() if "ChemicalSystem" in key]
        nums = [int(k.split("(")[1].split(")")[0]) - 1 for k in keys]
        for num, key in zip(nums, keys):
            atnums = self.file_object.read(key, "AtomicNumbers")
            if isinstance(atnums, int):
                atnums = [atnums]
            elements = [PT.get_symbol(atnum) for atnum in atnums]
            self.system_version_elements[num] = elements


def molecules_to_rkf(molecules: List[Molecule], fname, overwrite=False):
    """
    Convert a list of molecules to a .rkf file
    """

    if os.path.exists(fname) and overwrite:
        os.remove(fname)

    rkfout = RKFHistoryFile(fname, mode="wb")

    for i, mol in enumerate(molecules):
        rkfout.write_next(molecule=mol)
    rkfout.close()


def rkf_filter_regions(in_file, out_file, region: Union[str, Set[str]], mode="keep", overwrite=False):
    """
    Filter an ams.rkf trajectory file from region info. The new trajectory file only contains the atomic positions (not velocities).

    in_file: str
        path to an ams.rkf file

    out_file: str
        path to the output .rkf file

    region: str or set of str
        Names of regions

    mode: str
        'keep': keep all atoms in the named regions. 'remove': remove all atoms in the named regions.

    overwrite: bool
        Whether to overwrite out_file if it exists


    """
    rkf_in = RKFHistoryFile(in_file)
    mol = rkf_in.get_plamsmol()

    if os.path.exists(out_file) and overwrite:
        os.remove(out_file)

    rkf_out = RKFHistoryFile(out_file, mode="wb")
    rkf_out.store_mddata()

    if isinstance(region, str):
        region = set([region])

    for crd, cell in rkf_in(mol):
        new_mol = Molecule()
        new_mol.lattice = mol.lattice

        for at in mol:
            region_match = any(x in at.properties.region for x in region)
            if (mode == "keep" and region_match) or (mode == "remove" and not region_match):
                new_at = Atom(atnum=at.atnum, coords=at.coords, mddata={"PotentialEnergy": 0})
                new_at.properties = at.properties
                new_mol.add_atom(new_at)

        rkf_out.write_next(molecule=new_mol)

    rkf_out.close()

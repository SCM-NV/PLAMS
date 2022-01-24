#!/usr/bin/env python

import numpy

from ..mol.molecule import Molecule, Bond
from ..tools.periodic_table import PeriodicTable
from ..tools.kftools import KFFile
from ..tools.units import Units
#from ..interfaces.adfsuite.ams import AMSResults
from ..core.errors import PlamsError
from .trajectoryfile import TrajectoryFile

__all__ = ['RKFTrajectoryFile','write_general_section','write_molecule_section']

bohr_to_angstrom = Units.conversion_ratio('bohr','angstrom')

class RKFTrajectoryFile (TrajectoryFile) :
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
        def __init__ (self, filename, mode='rb', fileobject=None, ntap=None) :
                """
                Initiates an RKFTrajectoryFile object

                * ``filename``   -- The path to the RKF file
                * ``mode``       -- The mode in which to open the RKF file ('rb' or 'wb')
                * ``fileobject`` -- Optionally, a file object can be passed instead (filename needs to be set to None)
                * ``ntap``       -- If the file is in write mode, the number of atoms needs to be passed here
                """
                #TODO: If the mddata option is set to True, then the file created here works with AMSMovie and the analysis tools.
                #      To also make is work for restarts, two things have to be added:
                #      1. The final velocities have to be converted from bohr/fs to bohr/au (1/41.341373336493) 
                #         and stored in MDResuts%EndVelocities
                #      2. The final coordinates need to be copied to the Molecule section.

                self.position = 0
                if filename is not None :
                        #fileobject = KFFile(filename,autosave=False,keep_file_open=True)
                        fileobject = KFFile(filename,autosave=False)
                        #fileobject = KFFile(filename,autosave=False,fastsave=True)
                        # This fastsave option (no copying) was not worth it, so I removed it.
                        if fileobject is None :
                                raise PlamsError('KFFile %s not found.'%(rkfname))
                self.file_object = fileobject
                self.mode = mode

                self.ntap = 0
                if ntap is not None :
                        self.ntap = ntap
                self.firsttime = True
                self.coords = numpy.zeros((self.ntap,3))                # Only for reading purposes,
                                                                        # to avoid creating the array each time
                # PLAMS molecule related settings
                self.elements = ['H']*self.ntap
                self.current_molecule = None
                self.store_molecule = True # Even if True, the molecule attribute is only stored during iteration
        
                # RKF specific attributes
                self.nvecs = 3
                self.latticevecs = numpy.zeros((3,3))
                self.read_lattice = True               # Reading time can be saved by skipping the lattice info
                self.read_bonds = True
                self.cell = numpy.zeros((3,3))
                self.conect = None
                self.timestep = None
                self.include_mddata = False
                self.mddata = None
                self.mdunits = None
                self.include_historydata = False       # Any additional data along the history section will be stored
                self.historydata = None
                self.saving_freq = None                # By default the 'wb' file is only written upon closing
                                                       # Saving more often is much slower.
                # Skip to the trajectory part of the file (only if in read mode, because coords are required in header)
                if self.mode == 'rb' :
                        self._read_header()
                elif self.mode == 'wb' :
                        sections = self.file_object.sections()
                        if len(sections) > 0 : 
                                raise PlamsError ('RKF file %s already exists'%(filename))
                else :
                        raise PlamsError ('Mode %s is invalid. Only "rb" and "wb" are allowed.'%(self.mode))

        def store_mddata (self, rkf=None) :
                """
                Read/write an MDHistory section

                * ``rkf`` -- If in write mode an RKFTrajectoryFile object in read mode needs to be passed to extract unit info
                """
                self.include_mddata = True
                if 'r' in self.mode :
                        self._set_mddata_units()
                elif 'w' in self.mode :
                        if rkf is not None :
                                self.timestep = rkf.timestep
                                self._set_mdunits(rkf.mdunits)

        def store_historydata (self) :
                """
                Read/write non-standard entries in the History section
                """
                self.include_historydata = True

        def close (self) :
                """
                Execute all prior commands and cleanly close and garbage collect the RKF file
                """
                # Write the step info
                if self.timestep is not None and self.mode == 'wb' :
                        self.file_object.write('MDResults','StartStep',0)
                        self.file_object.write('MDResults','StartTime[fs]',0.)
                        nsteps = self.get_length()
                        self.file_object.write('MDResults','EndStep',nsteps-1)
                        self.file_object.write('MDResults','EndTime[fs]',(nsteps-1)*self.timestep)

                # Write to file
                if self.mode == 'wb' :
                        # First write the last frame into the Molecule section  
                        self._rewrite_molecule()
                        # Then write to file
                        self.file_object.save()
                del(self)

        def _rewrite_molecule (self) :
                """
                Overwrite the molecule section with the latest frame
                """
                crd,cell = self.read_last_frame()
                self._write_molecule_section(crd,cell)

        def _read_header (self, molecule_section='Molecule') :
                """
                Set up info required for reading frames
                """
                self.elements = self.file_object.read(molecule_section,'AtomSymbols')
                # If read from memory and not from file (write mode), it is already a list
                if isinstance(self.elements,str) :
                        self.elements = self.elements.split()
                self.elements = [el.split('.')[0] for el in self.elements]
                if self.file_object.reader is not None :
                        if ('MDHistory','Time(1)') in self.file_object :
                                times = self.file_object.read('MDHistory','Time(1)')
                                if isinstance(times,list) :
                                        self.timestep = times[1]
                self.ntap = len(self.elements)
                self.coords = numpy.zeros((self.ntap,3))
                try :
                        self.latticevecs = numpy.array(self.file_object.read(molecule_section,'LatticeVectors'))
                        #self.nvecs = int(len(self.cell)/3)
                        # New code 27-05-2020
                        self.nvecs = int(len(self.latticevecs)/3) # Why did I remove this line locally?!
                        self.latticevecs = self.latticevecs.reshape((self.nvecs,3))
                except (KeyError,AttributeError) :
                        pass

        def _set_mddata_units (self) :
                """
                Get the units for the mddata, if those are to be read
                """
                # Look for the items
                section = 'MDHistory'
                sections = self.file_object.get_skeleton()
                item_keys = [kn for kn in sections[section] if 'ItemName' in kn]
                items = [self.file_object.read(section,kn) for kn in item_keys]

                # Get the data for each item
                unit_dic = {}
                for item in items :
                        if '%s(units)'%(item) in self.file_object.reader._sections[section] :
                                unit_dic[item] = self.file_object.read(section,'%s(units)'%(item))

                self.mdunits = unit_dic

        def _write_header (self, coords, cell, molecule=None) :
                """
                Write Molecule info to file (elements, periodicity)
                """
                # First write the general section
                write_general_section(self.file_object)

                # Then write the input molecule
                self._update_celldata(cell)
                self._write_molecule_section(coords, cell, molecule=molecule)
                self._write_molecule_section(coords, cell, section='InputMolecule', molecule=molecule)
                if self.include_mddata :
                        # Start setting up the MDHistory section as well
                        self.file_object.write('MDHistory','blockSize',100)

                # Now make sure that it is possible to read from the file as well
                self._read_header()

        def _update_celldata (self, cell) :
                """
                Use the newly supplied cell to update the dimensionality of the system
                """
                shape = numpy.array(cell).shape
                if len(shape) == 2 :
                        self.nvecs = shape[0]
                        self.latticevecs = numpy.zeros((self.nvecs,3)) # Not really necessay

        def _write_molecule_section (self, coords, cell, section='Molecule', molecule=None) :
                """
                Write the molecule section
                """
                write_molecule_section(self.file_object,coords,cell,self.elements,section,molecule)

        def _set_mdunits (self, mdunits) :
                """
                Store the dictionary with MD Units
                """
                if self.include_mddata :
                        self.mdunits = mdunits

        def get_plamsmol (self) :
                """
                Extracts a PLAMS molecule object from the RKF file
                """
                section_dict = self.file_object.read_section('InputMolecule')
                if len(section_dict) == 0 :
                        section_dict = self.file_object.read_section('Molecule')
                plamsmol = Molecule._mol_from_rkf_section(section_dict)
                return plamsmol

        def read_frame (self, i, molecule=None) :
                """
                Reads the relevant info from frame ``i`` and returns it, or stores it in ``molecule``

                * ``i``        -- The frame number to be read from the RKF file
                * ``molecule`` -- |Molecule| object in which the new coordinates need to be stored
                """
                # Read the cell data
                cell = None
                if self.read_lattice :
                        try :
                                cell = self._read_cell_data (i)
                        except (KeyError,AttributeError) :
                                pass

                # Read the bond data
                conect = None
                bonds = None
                if self.read_bonds :
                        conect, bonds = self._read_bond_data(section='History', step=i)
                self.conect = conect

                # Read the coordinates, and possible pass them to molecule
                try :
                        self._read_coordinates(i, molecule, cell, bonds)
                        # This has changed self.coords behind the scenes
                except (KeyError,AttributeError) :
                        return None, None

                # Read and store any additional data in the history section
                if self.include_historydata :
                        self._store_historydata_for_step(i)
                # Read and store all MDData for this frame
                if self.include_mddata :
                        self._store_mddata_for_step(i)
                # Finalize
                if self.firsttime :
                        self.firsttime = False

                self.position = i
                return self.coords, cell

        def _read_coordinates (self, i, molecule, cell, bonds) :
                """
                Read the coordinates at step i, and possible pass them to molecule
                """
                if not self.coords.shape == (self.ntap,3) :
                        raise PlamsError('coords attribute has been changed outside the class')
                coords = self.coords.reshape(self.ntap*3)
                coords[:] = self.file_object.read('History', 'Coords(%i)'%(i+1))
                coords *= bohr_to_angstrom
                # This has changed self.coords behind the scenes

                # Create the molecule
                if isinstance(molecule,Molecule) :
                        cell_reduced = None
                        if cell is not None :
                                cell_reduced = cell[:self.nvecs]
                        self._set_plamsmol(self.coords,cell_reduced,molecule,bonds)

        def _read_cell_data (self, i) :
                """
                Read the cell data at step i
                """
                latticevecs = self.latticevecs.reshape(self.nvecs*3)
                latticevecs[:] = self.file_object.read('History','LatticeVectors(%i)'%(i+1)) #* bohr_to_angstrom
                latticevecs *= bohr_to_angstrom
                # This changed self.latticevecs behind the scenes
                #self.cell[:self.nvecs] = latticevecs
                self.cell[:self.nvecs] = self.latticevecs
                cell = self.cell
                return cell

        def _read_bond_data (self, section, step=None) :
                """
                Read the bond data from the rkf file
                """
                conect = None
                bonds = None
                try :
                        step_txt = ''
                        if step is not None :
                                step_txt = '(%i)'%(step+1)
                        #if ('History','Bonds.Index%s'%(step_txt)) in self.file_object :
                        indices = self.file_object.read(section,'Bonds.Index%s'%(step_txt))
                        connection_table = self.file_object.read(section,'Bonds.Atoms%s'%(step_txt))
                        if isinstance(connection_table,int) :
                                connection_table = [connection_table]
                        bond_orders = self.file_object.read(section,'Bonds.Orders%s'%(step_txt))
                        if isinstance(bond_orders,float) :
                                bond_orders = [bond_orders]
                        conect = {}
                        bonds = []
                        for i,(start,end) in enumerate(zip(indices[:-1],indices[1:])) :
                                if end-start > 0 :
                                        conect[i+1] = connection_table[start-1:end-1]
                                        for j in range(start-1,end-1) :
                                                bonds.append([i+1,connection_table[j],bond_orders[j]])
                        # Now correct the connection table
                        conect_sym = {}
                        for i, neighbors_i in conect.items() :
                                conect_sym[i] = neighbors_i
                                for j in neighbors_i :
                                        if not j in conect_sym.keys() :
                                                conect_sym[j] = []
                                        if j in conect.keys() : 
                                                conect_sym[j] = conect[j]
                                        if not i in conect_sym[j] :
                                                conect_sym[j].append(i)
                        conect = conect_sym     
                except (KeyError, AttributeError) :
                        pass
                return conect, bonds

        def _store_mddata_for_step (self, istep) :
                """
                Store the data from the MDHistory section
                """
                if self.mddata is None: self.mddata = {}
                section = 'MDHistory'

                # First get the block info
                blocksize = self.file_object.read(section, 'blockSize')
                nblocks = self.file_object.read(section, 'nBlocks')

                # Look for the items
                sections = self.file_object.get_skeleton()
                item_keys = [kn for kn in sections[section] if 'ItemName' in kn]
                items = [self.file_object.read(section,kn) for kn in item_keys]

                # Get the data for each item
                for item in items :
                        # First read the units (not needed?)
                        #units = ''
                        #if '%s(units)'%(item) in self.file_object.reader._sections[section] :
                        #        units = self.file_object.read(section,'%s(units)'%(item))

                        dim = self.file_object.read(section,'%s(dim)'%(item))
                        if dim == 1 and not self.file_object.read(section,'%s(perAtom)'%(item)):
                                # Stored in block format
                                block = int(istep/blocksize)
                                pos = istep%blocksize
                                values = self.file_object.read(section,'%s(%i)'%(item,block+1))
                                if not isinstance(values,list) : values = [values]
                                self.mddata[item] = values[pos]
                        else :  
                                self.mddata[item] = self.file_object.read(section,'%s(%i)'%(item,istep+1))

        def _store_historydata_for_step (self, istep) :
                """
                Store the extra data from the History section
                """
                if self.historydata is None: self.historydata = {}
                section = 'History'

                # Look for the items
                sections = self.file_object.get_skeleton()
                item_keys = [kn for kn in sections[section] if 'ItemName' in kn]
                items = [self.file_object.read(section,kn) for kn in item_keys]
                standard_items = ['Coords','nLatticeVectors','LatticeVectors','Bonds.Index','Bonds.Atoms','Bonds.Orders']
                items = [item for item in items if not item in standard_items]

                # Get the data for each item
                for item in items :
                        # Never stored in block format
                        self.historydata[item] = self.file_object.read(section,'%s(%i)'%(item,istep+1))

        def _is_endoffile (self) :
                """
                Reads and checks If the end of file is reached.
                """
                return ('History', 'Coords(%i)'%(self.position+1)) in self.file_object

        def read_next (self, molecule=None, read=True) :
                """
                Reads coordinates and lattice vectors from the current position of the cursor and returns it

                * ``molecule`` -- |Molecule| object in which the new coordinates need to be stored
                * ``read``     -- If set to False the cursor will move to the next frame without reading
                """
                if not read and not self.firsttime :
                        return self._move_cursor_without_reading()

                if self.firsttime :
                        self.firsttime = False

                crd, vecs = self.read_frame (self.position,molecule)
                self.position += 1
                return crd, vecs 

        def write_next (self, coords=None, molecule=None, cell=[0.,0.,0.], conect=None, historydata=None, mddata=None) :
                """
                Write frame to next position in trajectory file

                * ``coords``   -- A list or numpy array of (``ntap``,3) containing the system coordinates in angstrom
                * ``molecule`` -- A molecule object to read the molecular data from
                * ``cell``     -- A set of lattice vectors (or cell diameters for an orthorhombic system) in angstrom
                * ``conect``   -- A dictionary containing the connectivity info (e.g. {1:[2],2:[1]})
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
                if coords is not None :
                        if isinstance(coords,Molecule) :
                                raise PlamsError('The PLAMS molecule needs to be passed as the second argument (molecule)')

                if isinstance(molecule,Molecule) :
                        coords, cell, elements, conect, props = self._read_plamsmol(molecule)
                        if self.position == 0 : self.elements = elements
                # Make sure that the cell consists of vectors
                cell = self._convert_cell(cell)
                if conect is not None :
                        if len(conect) == 0 : conect = None
                self.conect = conect

                # Include a check on the size of coords?
                if len(coords) != len(self.elements) :
                        raise PlamsError('The coordinates do not match the rest of the trajectory')

                # If this is the first step, write the header
                if self.position == 0 :
                        self._write_header(coords,cell,molecule)

                # Define some local variables
                step = self.position
                if mddata is not None :
                        if 'Step' in mddata :
                                step = mddata['Step']
                # Energy should be read from mddata first, otherwise from historydata, otherwise set to zero
                energy = self._set_energy(mddata, historydata)
                if not self.include_historydata :
                        historydata = {}
                historydata['Energy'] = energy

                # Write the history section
                counter = 1
                counter = self._write_history_entry(step, coords, cell, conect, historydata, counter)

                if self.include_mddata and mddata is not None :
                        self._write_mdhistory_entry(mddata)

                self.position += 1

                if self.saving_freq is not None :
                        if self.position%self.saving_freq == 0 : self.file_object.save()

        def _set_energy (self, mddata, historydata) :
                """
                Looks if an energy is passed as input, and it not, sets to zero
                """
                energy = None
                if mddata is not None :
                        if 'PotentialEnergy' in mddata :
                                energy = mddata['PotentialEnergy']
                if energy is None :
                        if historydata is not None :
                                if 'Energy' in historydata :
                                        energy = historydata['Energy']
                if energy is None :     
                        energy = 0.
                return energy

        def _write_history_entry (self, step, coords, cell, conect, historydata=None, counter=1) :
                """
                Write the full entry into the History section
                """
                self.file_object.write('History','nEntries',self.position+1)
                self.file_object.write('History','currentEntryOpen',False)
                self._write_keydata_in_history('Step', counter, False, 1, self.position+1, step)
                counter += 1
                crd = [float(c)/bohr_to_angstrom for coord in coords for c in coord]
                self._write_keydata_in_history('Coords', counter, True, 3, self.position+1, crd)
                counter += 1
                #self._write_keydata_in_history('Energy', counter, False, 1, self.position+1, energy)
                #counter += 1
                if cell is not None :
                        self._write_keydata_in_history('nLatticeVectors', counter, False, 1, self.position+1, self.nvecs)
                        counter += 1
                        vecs = [float(v)/bohr_to_angstrom for vec in cell for v in vec]
                        # I should probably rethink the dimension of the lattice vectors (generalize it)
                        self._write_keydata_in_history('LatticeVectors', counter, False, [3,3], self.position+1, vecs)
                        counter += 1

                if historydata is not None :
                         counter = self._write_dictionary_to_history(historydata,'History',counter)
                #if gradients is not None :
                #        grd = [float(g) for grad in gradients for g in grad]
                #        self._write_keydata_in_history('Gradients', counter, True, 3, self.position+1, grd)

                #if stresstensor is not None :
                #        stre = [float(s) for stress in stresstensor for s in stress]
                #        self._write_keydata_in_history('StressTensor', counter, False, [3,3], self.position+1, stre)

                # Write the bond info
                if conect is not None :
                        counter = self._write_bonds_in_history(conect, counter, len(coords))

                return counter

        def _write_bonds_in_history (self, conect, counter, nats) :
                """
                Write the bond data into the history section
                """
                # Create the index list (correct for double counting)
                connection_table = [[at for at in conect[i+1] if at>i+1] if i+1 in conect else [] for i in range(nats)]
                numbonds = [len(neighbors) for neighbors in connection_table]

                indices = [sum(numbonds[:i])+1 for i in range(nats+1)]
                self.file_object.write('History','Bonds.Index(%i)'%(self.position+1),indices)
                self.file_object.write('History','ItemName(%i)'%(counter),'%s'%('Bonds.Index'))
                counter += 1

                # Flatten the connection table
                connection_table = [at for i in range(nats) for at in connection_table[i]]
                self.file_object.write('History','Bonds.Atoms(%i)'%(self.position+1),connection_table)
                self.file_object.write('History','ItemName(%i)'%(counter),'%s'%('Bonds.Atoms'))
                counter += 1

                bond_orders = [1.0 for bond in connection_table]
                self.file_object.write('History','Bonds.Orders(%i)'%(self.position+1),bond_orders)
                self.file_object.write('History','ItemName(%i)'%(counter),'%s'%('Bonds.Orders'))
                counter += 1

                return counter

        def _write_mdhistory_entry (self, mddata) :
                """
                Write the entry in the MDHistory section
                """
                counter = 1
                counter = self._write_dictionary_to_history(mddata,'MDHistory',counter)

        def _write_dictionary_to_history (self, data, section, counter=1) :
                """
                Add the entries of a dictionary to a History section 
                """
                self.file_object.write(section,'nEntries',self.position+1)
                self.file_object.write(section,'currentEntryOpen',False)
                for key, var in data.items() :
                        # Make sure that the entry is either a scalar or a 1D list
                        var = self._flatten_variable(var)
                        peratom = False
                        dim = 1
                        if isinstance(var,list) :
                                if len(var) % len(self.elements) == 0 :
                                        dim = int(len(var) / len(self.elements))
                                        peratom = True
                                else :
                                        dim = len(var)
                        self._write_keydata_in_history(key, counter, peratom, dim, self.position+1, var, section)
                        counter += 1
                return counter

        def _flatten_variable (self, var) :
                """
                Make sure that the variable is a Python 1D list (not numpy)
                """
                while(1) :
                        if isinstance(var,list) or isinstance(var,numpy.ndarray) :
                                if len(var) == 0 : break
                                if isinstance(var[0],list) or isinstance(var[0],numpy.ndarray) :
                                        var = [v for varitem in var for v in varitem]
                                else :
                                        if isinstance(var[0],numpy.int64) :
                                                var = [int(v) for v in var]
                                        elif isinstance(var[0],numpy.float64) :
                                                var = [float(v) for v in var]
                                        break
                        else :
                                if isinstance(var,numpy.int64) : 
                                        var = int(var)
                                elif isinstance(var,numpy.float64) :
                                        var = float(var)
                                break
                return var

        def _write_keydata_in_history(self, key, i, perAtom, dim, step, values, section='History') :
                """
                Write all data about a key value in KFFile
                """
                # Some data only needs to be printed once
                printstartdata = False
                if step == 1 :
                        printstartdata = True

                # Block code: if the data is to be written as blocks, then step and values need to be replaced.
                if section == 'MDHistory' :
                        step, values = self._get_block_info (key, perAtom, dim, step, values, section)
                                
                # The rest should be independent on format (block or individual)
                self.file_object.write(section,'%s(%i)'%(key,step),values)
                if printstartdata :
                        self.file_object.write(section,'ItemName(%i)'%(i),'%s'%(key))
                        self.file_object.write(section,'%s(perAtom)'%(key),perAtom)
                        self.file_object.write(section,'%s(dim)'%(key),dim)
                        if section == 'MDHistory' and self.mdunits is not None :
                                if key in self.mdunits :
                                        self.file_object.write(section,'%s(units)'%(key),self.mdunits[key])

        def _get_block_info (self, key, perAtom, dim, step, values, section) :
                """
                If the data is to be written as blocks, then step and values need to be replaced.
                """
                if dim==1 and not perAtom :
                        try :
                                blocksize = self.file_object.read(section,'blockSize')
                        except AttributeError :
                                raise Exception ('Set include_mddata to write the MD section.')
                        iblock = int((step-1)/blocksize) + 1
                        if step%blocksize != 1 :
                                try :
                                        old_values = self.file_object.read(section,'%s(%i)'%(key,iblock))
                                        if not isinstance(old_values,list) :
                                                old_values = [old_values]
                                except AttributeError :
                                        old_values = []
                                values = old_values + [values] # Values is a scalar
                        else :
                                self.file_object.write(section,'nBlocks',iblock)
                        step = iblock
                return step, values

        def rewind (self, nframes=None) :
                """
                Rewind the file either by ``nframes`` or to the first frame

                *   ``nframes`` -- The number of frames to rewind
                """
                self.firsttime = True
                self.position = 0

        def get_length (self) :
                """
                Get the number of frames in the file
                """
                nsteps = self.file_object.read('History', 'nEntries')
                return nsteps

        def read_last_frame (self, molecule=None) :
                """
                Reads the last frame from the file
                """
                nsteps = self.get_length()
                crd,cell = self.read_frame(nsteps-1, molecule)
                return crd, cell

def write_general_section (rkf) :
        """
        Write the General section of the RKF file
        """
        rkf.write('General','file-ident','RKF')
        rkf.write('General','termination status','NORMAL TERMINATION')
        rkf.write('General','program','ams')
        rkf.write('General','user input',' ')

def write_molecule_section (rkf, coords=None, cell=None, elements=None, section='Molecule', molecule=None) :
        """
        Write the molecule section

        Note: Currently does not write bonds
        """
        if molecule is not None :
                if coords is None :
                        coords = molecule.as_array()
                if cell is None and len(molecule.lattice)>0 :
                        cell = molecule.lattice
                if elements is None :
                        elements = [at.symbol for at in molecule.atoms]

        # Then write the input molecule
        charge = 0.
        if molecule is not None :
                if 'charge' in molecule.properties :
                        charge = molecule.properties.charge
        element_numbers = [PeriodicTable.get_atomic_number(el) for el in elements]

        rkf.write(section,'nAtoms',len(elements))
        rkf.write(section,'AtomicNumbers',element_numbers)
        rkf.write(section,'AtomSymbols',elements)
        crd = [Units.convert(float(c),'angstrom','bohr') for coord in coords for c in coord]
        rkf.write(section,'Coords',crd)
        rkf.write(section,'Charge',charge)
        if cell is not None :
                rkf.write(section,'nLatticeVectors',len(cell))
                vecs = [Units.convert(float(v),'angstrom','bohr') for vec in cell for v in vec]
                rkf.write(section,'LatticeVectors',vecs)
        # Should it write bonds?
        # Write atom properties
        if molecule is not None :
                suffixes = []
                present = False
                for at in molecule :
                        if 'suffix' in at.properties :
                                present = True
                                suffixes.append(at.properties.suffix)
                        else :
                                suffixes.append('')
                if present :
                        rkf.write(section,'EngineAtomicInfo','\x00'.join(suffixes))
                # Also add a bond section
                if len(molecule.bonds) > 0 :
                        bond_indices = [sorted([iat for iat in molecule.index(bond)]) for bond in molecule.bonds]
                        atoms_from = [bond[0] for bond in bond_indices]
                        atoms_to = [bond[1] for bond in bond_indices]
                        orders = [bond.order for bond in molecule.bonds]
                        rkf.write(section,'fromAtoms',atoms_from)
                        rkf.write(section,'toAtoms',atoms_to)
                        rkf.write(section,'bondOrders',orders)


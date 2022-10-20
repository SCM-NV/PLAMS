#!/usr/bin/env amspython
import os
import numpy as np
from typing import List, Union
from ...core.private import saferun
from ...mol.molecule import Molecule
from ...core.settings import Settings
from .ase import toASE, fromASE
import tempfile
import subprocess

__all__ = ['PackMolStructure', 'PackMol', 'packmol_liquid', 'packmol_solid_liquid', 'packmol_solid_liquid_mixture', 'packmol_mixture', 'packmol', 'packmol_on_slab', 'PackMolFailedException', 'packmol_microsolvation']

class PackMolFailedException(Exception):
    pass


class PackMolStructure:
    def __init__(self, molecule : Molecule, n_molecules:int=None, n_atoms:int=None, box_bounds:List[float]=None, density:float=None, fixed:bool=False, sphere:bool=False):
        """

        Class representing a packmol structure.

        molecule : Molecule
            The molecule

        n_molecules: int
            The number of molecules to insert

        n_atoms: int
            An approximate number of atoms to insert

        box_bounds: list of float
            [xmin, ymin, zmin, xmax, ymax, zmax] in angstrom. The min values should all be 0, i.e. [0., 0., 0., xmax, ymax, zmax]

        density: float
            Density in g/cm^3

        fixed: bool
            Whether the structure should be fixed at its original coordinates.

        sphere: bool
            Whether the molecules should be packed in a sphere. The radius is determined by getting the volume from the box bounds! Cannot be combined with ``fixed`` (``fixed`` takes precedence).

        """
        self.molecule = molecule
        if fixed:
            assert(n_molecules is None or n_molecules == 1)
            #assert(box_bounds is None)
            assert(density is None)
            self.n_molecules = 1
            if molecule.lattice and len(molecule.lattice) == 3:
                self.box_bounds = [0., 0., 0., molecule.lattice[0][0], molecule.lattice[1][1], molecule.lattice[2][2]]
            else:
                self.box_bounds = None
            self.fixed = True
            self.sphere = False
        else:
            if box_bounds and density:
                if n_molecules or n_atoms:
                    raise ValueError("Cannot set all n_molecules or n_atoms together with (box_bounds AND density)")
                n_molecules = self._get_n_molecules_from_density_and_box_bounds(self.molecule, box_bounds, density)
            assert(n_molecules or n_atoms)
            self.n_molecules = n_molecules or self._get_n_molecules(self.molecule, n_atoms)
            assert(box_bounds or density)
            self.box_bounds = box_bounds or self._get_box_bounds(self.molecule, self.n_molecules, density)
            self.fixed = False
            self.sphere = sphere

    def _get_n_molecules_from_density_and_box_bounds(self, molecule:Molecule, box_bounds:List[float], density:float):
        """ density in g/cm^3 """ 
        molecule_mass = molecule.get_mass(unit='g')
        volume_ang3 = self.get_volume(box_bounds)
        volume_cm3 = volume_ang3 * 1e-24
        n_molecules = int(density * volume_cm3 / molecule_mass)
        return n_molecules

    def get_volume(self, box_bounds=None):
        bb = box_bounds or self.box_bounds
        vol = (bb[3]-bb[0])*(bb[4]-bb[1])*(bb[5]-bb[2])
        return vol

    def _get_n_molecules(self, molecule:Molecule, n_atoms:int):
        return n_atoms // len(molecule)

    def _get_box_bounds(self, molecule:Molecule, n_molecules:int, density:float):
        mass = n_molecules * molecule.get_mass(unit='g')
        volume_cm3 = mass / density
        volume_ang3 = volume_cm3 * 1e24
        side_length = volume_ang3 ** (1/3.0)
        return [0., 0., 0., side_length, side_length, side_length]

    def get_input_block(self, fname, tolerance):
        if self.fixed:
            ret = f'''
            structure {fname}
            number 1
            fixed 0. 0. 0. 0. 0. 0.
            avoid_overlap yes
            end structure
            '''
        elif self.sphere:
            vol = self.get_volume()
            # vol = 4*pi*r^3 /3 
            # radius = (3*vol/(4*pi))**0.33333
            radius = (3*vol/(4*3.14159))**0.3333
            ret = f'''
            structure {fname}
              number {self.n_molecules}
              inside sphere 0. 0. 0. {radius}
            end structure
            '''
        else:
            box_string = f'{self.box_bounds[0]+tolerance/2} {self.box_bounds[1]+tolerance/2} {self.box_bounds[2]+tolerance/2} {self.box_bounds[3]-tolerance/2} {self.box_bounds[4]-tolerance/2} {self.box_bounds[5]-tolerance/2}'
            ret = f'''
            structure {fname}
              number {self.n_molecules}
              inside box {box_string}
            end structure

        '''
        return ret
        

class PackMol:

    def __init__(self, tolerance=2.0, structures:List[PackMolStructure]=None, executable=None):
        """
        Class for setting up and running packmol.

        tolerance: float
            The packmol tolerance (approximate minimum interatomic distance)

        structures: list of PackMolStructure
            Structures to insert

        executable: str
            Path to the packmol executable. If not specified, $AMSBIN/packmol.exe will be used.

        Note: users are not recommended to use this class directly, but
        instead use the ``packmol_mixture`` or ``packmol_solid_liquid_mixture``
        functions.

        """
        self.tolerance = tolerance
        self.filetype = 'xyz'
        self.output = 'packmol_output.xyz'
        self.structures = structures or []
        self.executable = executable or os.path.join(os.path.expandvars('$AMSBIN'), 'packmol.exe')
        assert(os.path.exists(self.executable))

    def add_structure(self, structure: PackMolStructure):
        self.structures.append(structure)

    def _get_complete_box_bounds(self):
        min_x = min(s.box_bounds[0] for s in self.structures if s.box_bounds is not None)
        min_y = min(s.box_bounds[1] for s in self.structures if s.box_bounds is not None)
        min_z = min(s.box_bounds[2] for s in self.structures if s.box_bounds is not None)
        max_x = min(s.box_bounds[3] for s in self.structures if s.box_bounds is not None)
        max_y = min(s.box_bounds[4] for s in self.structures if s.box_bounds is not None)
        max_z = min(s.box_bounds[5] for s in self.structures if s.box_bounds is not None)

        #return min_x, min_y, min_z, max_x+self.tolerance, max_y+self.tolerance, max_z+self.tolerance
        return min_x, min_y, min_z, max_x, max_y, max_z

    def _get_complete_lattice(self):
        """
            returns a 3x3 list using the smallest and largest x/y/z box_bounds for all structures
        """
        if any(s.sphere for s in self.structures):
            return []
        min_x, min_y, min_z, max_x, max_y, max_z = self._get_complete_box_bounds()
        return [[max_x-min_x, 0., 0.], [0., max_y-min_y, 0.], [0., 0., max_z-min_z]]

    def run(self):
        """
            returns: a Molecule with the packed structures
        """

        assert(os.path.exists(self.executable))

        output_molecule = Molecule()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_fname = os.path.join(tmpdir, 'output.xyz')
            input_fname = os.path.join(tmpdir, 'input.inp')
            with open(input_fname, 'w') as input_file:
                input_file.write(f'tolerance {self.tolerance}\n')
                input_file.write(f'filetype xyz\n')
                input_file.write(f'output {output_fname}\n')

                for i, structure in enumerate(self.structures):
                    structure_fname = os.path.join(tmpdir, f'structure{i}.xyz')
                    structure.molecule.write(structure_fname)
                    input_file.write(structure.get_input_block(structure_fname, tolerance=2.0))

            #with open(input_fname, 'r') as f:
            #    for line in f:
            #        print(line)

            # cannot feed stdin as a string into packmol for some reason
            # it seems to need a file
            my_input = open(input_fname, 'r')
            saferun(self.executable, stdin=my_input, stdout=subprocess.DEVNULL)
            my_input.close()

            if not os.path.exists(output_fname):
                raise PackMolFailedException("Packmol failed. It may work if you try a lower density.")
            output_molecule = Molecule(output_fname) # without periodic boundary conditions
            output_molecule.lattice = self._get_complete_lattice()

        return output_molecule

def packmol(molecules:Union[List[Molecule],Molecule], mole_fractions:List[float]=None, density:float=None, n_atoms:int=None, box_bounds:List[float]=None, n_molecules:Union[List[int],int]=None, sphere:bool=False, region_names:List[str]=None, return_details:bool=False, executable:str=None):
    """

        Create a fluid of the given ``molecules``. The function will use the
        given input parameters and try to obtain good values others. You *must*
        specify ``density`` and/or ``box_bounds``.

        molecules : Molecule or list of Molecule
            The molecules to pack

        mole_fractions : list of float
            The mole fractions (in the same order as ``molecules``). Cannot be combined with ``n_molecules``. If not given, an equal (molar) mixture of all components will be created.

        density: float
            The total density (in g/cm^3) of the fluid

        n_atoms: int
            The (approximate) number of atoms in the final mixture

        box_bounds: list of float (length 6)
            The box in which to pack the molecules. The box is orthorhombic and should be specified as [xmin, ymin, zmin, xmax, ymax, zmax]. The minimum values should all be set to 0, i.e. set box_bounds=[0., 0., 0., xmax, ymax, zmax]. If not specified, a cubic box of appropriate dimensions will be used.

        n_molecules : int or list of int
            The (exact) number of molecules for each component (in the same order as ``molecules``). Cannot be combined with ``mole_fractions``.

        region_names : str or list of str
            Populate the region information for each atom. Should have the same length and order as ``molecules``. My default the regions are named ``mol0``, ``mol1``, etc.

        return_details : bool
            Return a 2-tuple (Molecule, dict) where the dict has keys like 'n_molecules', 'mole_fractions', 'density', etc. They contain the actual details of the returned molecule, which may differ slightly from the requested quantities.

            Returned keys:

            * 'n_molecules': list of integer with actually added number of molecules
            * 'mole_fractions': list of float with actually added mole fractions
            * 'density': float, gives the density in g/cm^3
            * 'n_atoms': int, the number of atoms in the returned molecule
            * 'molecule_type_indices': list of int of length n_atoms. For each atom, give an integer index for which TYPE of molecule it belongs to.
            * 'molecule_indices': list of int of length n_atoms. For each atom, give an integer index for which molecule it belongs to
            * 'atom_indices_in_molecule': list of int of length n_atoms. For each atom, give an integer index for which position in the molecule it is.

        executable : str
            The path to the packmol executable. If not specified, ``$AMSBIN/packmol.exe`` will be used (which is the correct path for the Amsterdam Modeling Suite).

        Useful combinations:

        * ``mole_fractions``, ``density``, ``n_atoms``: Create a mixture with a given density and approximate number of atoms (a cubic box will be created)

        * ``mole_fractions``, ``density``, ``box_bounds``: Create a mixture with a given density inside a given box (the number of molecules will approximately match the density and mole fractions)

        * ``n_molecules``, ``density``: Create a mixture with the given number of molecules and density (a cubic box will be created)

        * ``n_molecules``, ``box_bounds``: Create a mixture with the given number of molecules inside the given box

        Example:

        .. code-block:: python

            packmol_mixture(molecules=[from_smiles('O'), from_smiles('C')], 
                            mole_fractions=[0.8, 0.2],
                            density=0.8, 
                            n_atoms=100)

        Returns: a Molecule or tuple (Molecule, dict)
            If return_details=False, return a Molecule. If return_details=True, return a tuple.


    """
    assert(not (n_atoms and n_molecules))
    assert(n_atoms or n_molecules or density)
    assert(density or box_bounds)
    assert(not (mole_fractions and n_molecules))

    def tolist(x):
        return x if isinstance(x, list) else [x]

    molecules = tolist(molecules)
    if mole_fractions is None:
        mole_fractions = [1.0/len(molecules)] * len(molecules)

    if n_molecules:
        n_molecules = tolist(n_molecules)

    
    xs = np.array(mole_fractions)
    atoms_per_mol = np.array([len(a) for a in molecules])
    masses = np.array([m.get_mass(unit='g') for m in molecules])

    coeffs = None

    if n_molecules:
        coeffs = np.int_(n_molecules)
    elif n_atoms:
        coeff_0 = n_atoms / np.dot(xs, atoms_per_mol)
        coeffs_floats = xs * coeff_0
        coeffs = np.int_(np.round(coeffs_floats))

    if (n_atoms or n_molecules) and density and not box_bounds:
        mass = np.dot(coeffs, masses)
        volume_cm3 = mass / density
        volume_ang3 = volume_cm3 * 1e24
        side_length = volume_ang3 ** (1/3.0)
        box_bounds =  [0., 0., 0., side_length, side_length, side_length]
    elif box_bounds and density and not n_molecules:
        volume_cm3 = (box_bounds[3]-box_bounds[0])*(box_bounds[4]-box_bounds[1])*(box_bounds[5]-box_bounds[2]) * 1e-24
        mass_g = volume_cm3 * density
        coeffs = mass_g / np.dot(xs, masses)
        coeffs = xs * coeffs
        coeffs = np.int_(np.round(coeffs))

    if coeffs is None:
        raise ValueError(f"Illegal combination of options: n_atoms={n_atoms}, n_molecules={n_molecules}, box_bounds={box_bounds}, density={density}")

    pm = PackMol(executable=executable)
    for i, (mol, n_mol) in enumerate(zip(molecules, coeffs)):
        if sphere:
            if i == 0:
                pm.add_structure(PackMolStructure(mol, n_molecules=n_mol, box_bounds=box_bounds, fixed=True, sphere=False))
            else:
                pm.add_structure(PackMolStructure(mol, n_molecules=n_mol, box_bounds=box_bounds, fixed=False, sphere=True))
        else:
            pm.add_structure(PackMolStructure(mol, n_molecules=n_mol, box_bounds=box_bounds))

    out = pm.run()

    # packmol returns the molecules sorted
    molecule_type_indices = [] # [0,0,0,...,1,1,1] # two different molecules with 3 and 5 atoms
    molecule_indices = [] # [0,0,0,1,1,1,2,2,2,....,58,58,58,58,58,59,59,59,59,59] # two different molecules with 3 and 5 atoms
    atom_indices_in_molecule = [] # [0,1,2,0,1,2,...,0,1,2,3,4,0,1,2,3,4] 
    current = 0
    for i, (mol, n_mol) in enumerate(zip(molecules, coeffs)):
        molecule_type_indices += [i]*n_mol*len(mol)
        atom_indices_in_molecule += list(range(len(mol)))*n_mol

        temp = list(range(current, current+n_mol))
        molecule_indices += list(np.repeat(temp, len(mol)))
        current += n_mol
    assert len(molecule_type_indices) == len(out)
    assert len(molecule_indices) == len(out)
    assert len(atom_indices_in_molecule) == len(out)

    details = {
        'n_molecules': coeffs.tolist(),
        'mole_fractions': (coeffs/np.sum(coeffs)).tolist() if np.sum(coeffs) > 0 else [0.0]*len(coeffs),
        'n_atoms': len(out),
        'molecule_type_indices': molecule_type_indices, # for each atom, indicate which type of molecule it belongs to by an integer index (starts with 0)
        'molecule_indices': molecule_indices, # for each atoms, indicate which molecule it belongs to by an integer index (starts with 0)
        'atom_indices_in_molecule': atom_indices_in_molecule,
    }
    try:
        details['density']  = out.get_density()*1e-3,
    except ValueError:
        details['density'] = None
        pass # if not periodict

    if region_names:
        region_names = tolist(region_names)
    else:
        region_names = [f'mol{i}' for i in range(len(molecules))]

    for at, molindex in zip(out, molecule_type_indices):
        at.properties.suffix = f'region={region_names[molindex]}'

    if return_details:
        return out, details

    return out

def packmol_mixture(molecules:List[Molecule], mole_fractions:List[float]=None, density:float=None, n_atoms:int=None, box_bounds:List[float]=None, n_molecules:List[int]=None, executable:str=None):
    """ Deprecated """
    return packmol(molecules=molecules, mole_fractions=mole_fractions, density=density, n_atoms=n_atoms, box_bounds=box_bounds, n_molecules=n_molecules, executable=executable)


def packmol_liquid(molecule:Molecule, density:float=None, n_atoms:int=None, box_bounds:List[float]=None, n_molecules:int=None, executable:str=None):
    """ Deprecated """
    return packmol(molecules=molecule, density=density, n_atoms=n_atoms, box_bounds=box_bounds, n_molecules=n_molecules, executable=executable)

def get_packmol_solid_liquid_box_bounds(slab:Molecule):
    slab_max_z = max(at.coords[2] for at in slab)
    slab_min_z = min(at.coords[2] for at in slab)
    liquid_min_z = slab_max_z
    liquid_max_z = liquid_min_z + slab.lattice[2][2] - (slab_max_z-slab_min_z)
    box_bounds = [0., 0., liquid_min_z+1.5, slab.lattice[0][0], slab.lattice[1][1], liquid_max_z-1.5]
    return box_bounds

def packmol_on_slab(slab:Molecule, molecules:Union[List[Molecule],Molecule], density:float, mole_fractions:List[float]=None, executable:str=None):
    """

    Creates a solid/liquid interface with an approximately correct density. The
    density is calculated for the volume not occupied by the slab (+ 1.5
    angstrom buffer at each side of the slab).

    Returns: a Molecule

    slab : Molecule
        The system must have a 3D lattice (including a vacuum gap along z) and be orthorhombic. The vacuum gap will be filled with the liquid.

    For the other arguments, see ``packmol_mixture``.

    Example:

    .. code-block:: python

        packmol_solid_liquid_mixture(slab=slab_3d_with_vacuum_gap, 
                                     molecules=[from_smiles('O'), from_smiles('C')], 
                                     mole_fractions=[0.8, 0.2], 
                                     density=0.8)

    """
    out = slab.copy()
    box_bounds = get_packmol_solid_liquid_box_bounds(out)
    liquid = packmol(molecules=molecules, mole_fractions=mole_fractions, density=density, box_bounds=box_bounds, executable=executable)
    out.add_molecule(liquid)

    for at in out:
        if at.coords[2] > out.lattice[2][2]:
            at.translate([0, 0, -out.lattice[2][2]])
    return out

def packmol_solid_liquid_mixture(slab:Molecule, molecules:List[Molecule], mole_fractions:List[float], density:float, executable:str=None):
    """ Deprecated """
    return packmol_on_slab(slab=slab, molecules=molecules, mole_fractions=mole_fractions, density=density, executable=executable)

def packmol_solid_liquid(slab:Molecule, molecule:Molecule, density:float, executable:str=None): 
    """ Deprecated """
    return packmol_on_slab(slab=slab, molecules=molecule, mole_fractions=[1.0], density=density, executable=executable)

def get_n_from_density_and_box_bounds(molecule, box_bounds, density):
    molecule_mass = molecule.get_mass(unit='g')
    volume_ang3 = (box_bounds[3]-box_bounds[0])*(box_bounds[4]-box_bounds[1])*(box_bounds[5]-box_bounds[2])
    volume_cm3 = volume_ang3 * 1e-24
    n_molecules = int(density * volume_cm3 / molecule_mass)
    return n_molecules

def packmol_microsolvation(solute:Molecule, solvent:Molecule, density:float=1.0, threshold:float=3.0, executable:str=None):
    ase_solute = toASE(solute)
    com = np.mean(ase_solute.get_positions(), axis=0)
    ase_solute.set_positions(ase_solute.get_positions() - com)
    box_bounds = [0, 0, 0] + list(np.max(ase_solute.get_positions(), axis=0) - np.min(ase_solute.get_positions(), axis=0) + 3*threshold)

    plams_solute = fromASE(ase_solute)
    n_solvent = get_n_from_density_and_box_bounds(solvent, box_bounds, density=density)

    plams_solvated = packmol([plams_solute, solvent], n_molecules=[1, n_solvent], box_bounds=box_bounds, sphere=True)

    plams_solvated.guess_bonds()
    molecule_indices = plams_solvated.get_molecule_indices() # [[0,1,2],[3,4],[5,6],...]

    ase_solvated = toASE(plams_solvated)
    D = ase_solvated.get_all_distances()[:len(solute)]
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

    newatoms = ase_solvated[complete_indices]
    newmolecule = fromASE(newatoms)

    return newmolecule






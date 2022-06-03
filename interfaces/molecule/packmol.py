#!/usr/bin/env amspython
import os
import numpy as np
from typing import List
from ...core.private import saferun
from ...mol.molecule import Molecule
from ...core.settings import Settings
import tempfile
import subprocess

__all__ = ['PackMolStructure', 'PackMol', 'packmol_liquid', 'packmol_solid_liquid', 'packmol_solid_liquid_mixture', 'packmol_mixture']

class PackMolStructure:
    def __init__(self, molecule : Molecule, n_molecules:int=None, n_atoms:int=None, box_bounds:List[float]=None, density:float=None, fixed:bool=False):
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

        """
        self.molecule = molecule
        if fixed:
            assert(n_molecules is None)
            assert(box_bounds is None)
            assert(density is None)
            self.n_molecules = 1
            if molecule.lattice and len(molecule.lattice) == 3:
                self.box_bounds = [0., 0., 0., molecule.lattice[0][0], molecule.lattice[1][1], molecule.lattice[2][2]]
            else:
                self.box_bounds = None
            self.fixed = True
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
        min_x = min(s.box_bounds[0] for s in self.structures)
        min_y = min(s.box_bounds[1] for s in self.structures)
        min_z = min(s.box_bounds[2] for s in self.structures)
        max_x = min(s.box_bounds[3] for s in self.structures)
        max_y = min(s.box_bounds[4] for s in self.structures)
        max_z = min(s.box_bounds[5] for s in self.structures)

        #return min_x, min_y, min_z, max_x+self.tolerance, max_y+self.tolerance, max_z+self.tolerance
        return min_x, min_y, min_z, max_x, max_y, max_z

    def _get_complete_lattice(self):
        """
            returns a 3x3 list using the smallest and largest x/y/z box_bounds for all structures
        """
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

            output_molecule = Molecule(output_fname) # without periodic boundary conditions
            output_molecule.lattice = self._get_complete_lattice()

        return output_molecule

def packmol_mixture(molecules:List[Molecule], mole_fractions:List[float]=None, density:float=None, n_atoms:int=None, box_bounds:List[float]=None, n_molecules:List[int]=None, executable:str=None):
    """

        Create a mixture of the given ``molecules``. The function will use the
        given input parameters and try to obtain good values others. You *must*
        specify ``density`` and/or ``box_bounds``.

        molecules : list of Molecule
            The molecules to pack

        mole_fractions : list of float
            The mole fractions (in the same order as ``molecules``). Cannot be combined with ``n_molecules``

        density: float
            The total density (in g/cm^3) of the mixture

        n_atoms: int
            The (approximate) number of atoms in the final mixture

        box_bounds: list of float (length 6)
            The box in which to pack the molecules. The box is orthorhombic and should be specified as [xmin, ymin, zmin, xmax, ymax, zmax]. The minimum values should all be set to 0, i.e. set box_bounds=[0., 0., 0., xmax, ymax, zmax]. If not specified, a cubic box of appropriate dimensions will be used.

        n_molecules : list of int
            The (exact) number of molecules for each component (in the same order as ``molecules``). Cannot be combined with ``mole_fractions``.

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


    """
    assert(not (n_atoms and n_molecules))
    assert(n_atoms or n_molecules or density)
    assert(density or box_bounds)
    assert(not (mole_fractions and n_molecules))

    
    xs = np.array(mole_fractions)
    atoms_per_mol = np.array([len(a) for a in molecules])
    masses = np.array([m.get_mass(unit='g') for m in molecules])

    coeffs = None

    if n_molecules:
        coeffs = np.int_(n_molecules)
    elif n_atoms:
        coeff_0 = n_atoms / np.dot(xs, atoms_per_mol)
        coeffs_floats = xs * coeff_0
        coeffs = np.int_(coeffs_floats)

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
        coeffs = np.int_(coeffs)

    if coeffs is None:
        raise ValueError(f"Illegal combination of options: n_atoms={n_atoms}, n_molecules={n_molecules}, box_bounds={box_bounds}, density={density}")

    pm = PackMol(executable=executable)
    for mol, n_mol in zip(molecules, coeffs):
        pm.add_structure(PackMolStructure(mol, n_molecules=n_mol, box_bounds=box_bounds))

    out = pm.run()
    return out


def packmol_liquid(molecule:Molecule, density:float=None, n_atoms:int=None, box_bounds:List[float]=None, n_molecules:int=None, executable:str=None):
    """

    Creates a liquid/gas of the provided ``molecule``. Returns: a Molecule.
    
    molecule : Molecule
        The molecule to pack

    The other arguments are described for ``packmol_mixture``.

    Examples:

    .. code-block:: python

        packmol_liquid(molecule=from_smiles('O'), density=1.0, n_atoms=100)
        packmol_liquid(molecule=from_smiles('O'), n_molecules=64, box_bounds=[0., 0., 0., 12.2, 12.2, 12.2])

    """

    if n_molecules:
        n_molecules = [n_molecules]
        mole_fractions = None
    else:
        mole_fractions = [1.0]
    return packmol_mixture([molecule], mole_fractions=mole_fractions, n_atoms=n_atoms, density=density, box_bounds=box_bounds, n_molecules=n_molecules, executable=executable)

def get_packmol_solid_liquid_box_bounds(slab:Molecule):
    slab_max_z = max(at.coords[2] for at in slab)
    slab_min_z = min(at.coords[2] for at in slab)
    liquid_min_z = slab_max_z
    liquid_max_z = liquid_min_z + slab.lattice[2][2] - (slab_max_z-slab_min_z)
    box_bounds = [0., 0., liquid_min_z+1.5, slab.lattice[0][0], slab.lattice[1][1], liquid_max_z-1.5]
    return box_bounds

def packmol_solid_liquid_mixture(slab:Molecule, molecules:List[Molecule], mole_fractions:List[float], density:float, executable:str=None):
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
    liquid = packmol_mixture(molecules=molecules, mole_fractions=mole_fractions, density=density, box_bounds=box_bounds, executable=executable)
    out.add_molecule(liquid)

    for at in out:
        if at.coords[2] > out.lattice[2][2]:
            at.translate([0, 0, -out.lattice[2][2]])
    return out

def packmol_solid_liquid(slab:Molecule, molecule:Molecule, density:float, executable:str=None): 
    """

    Creates a solid/liquid interface with an approximately correct density. The
    density is calculated for the volume not occupied by the slab (+ 1.5
    angstrom buffer at each side of the slab).

    Returns: a Molecule

    For details, see ``packmol_solid_liquid_mixture`` and ``packmol_mixture``

    Example:

    .. code-block:: python

        packmol_solid_liquid(slab=slab_3d_with_vacuum_gap, molecule=from_smiles('O'), density=1.0)

    """

    return packmol_solid_liquid_mixture(slab, [molecule], mole_fractions=[1.0], density=density, executable=executable)


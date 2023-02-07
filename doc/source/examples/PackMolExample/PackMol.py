#!/usr/bin/env amspython
# coding: utf-8

# ## Initial imports

from scm.plams import *
from ase.optimize import BFGS
from ase.build import molecule as ase_build_molecule
from ase.visualize.plot import plot_atoms
from ase.build import fcc111
import matplotlib.pyplot as plt


# ## Helper functions

def printsummary(mol, details=None):
    s = f'{len(mol)} atoms, density = {mol.get_density()*1e-3:.3f} g/cm^3, box = {mol.lattice[0][0]:.3f}, {mol.lattice[1][1]:.3f}, {mol.lattice[2][2]:.3f}, formula = {mol.get_formula()}'
    if details:
        s+= f'\n#added molecules per species: {details["n_molecules"]}, mole fractions: {details["mole_fractions"]}'
    print(s)
    
def show(mol, figsize=None, **kwargs):
    """ Show a molecule in a Jupyter notebook """
    plt.figure(figsize=figsize or (2,2))
    plt.axis('off')
    plot_atoms(toASE(mol), **kwargs)


# ## Liquid water (fluid with 1 component)
# First, create the gasphase molecule:

water = from_smiles('O')
show(water)


print('pure liquid from approximate number of atoms and exact density (in g/cm^3), cubic box with auto-determined size')
out = packmol(water, n_atoms=194, density=1.0)
printsummary(out)
out.write('water-1.xyz')
show(out)


print('pure liquid from approximate density (in g/cm^3) and an orthorhombic box')
out = packmol(water, density=1.0, box_bounds=[0., 0., 0., 8., 12., 14.])
printsummary(out)
out.write('water-2.xyz')
show(out)


print('pure liquid with explicit number of molecules and exact density')
out = packmol(water, n_molecules=64, density=1.0)
printsummary(out)
out.write('water-3.xyz')
show(out)


print('pure liquid with explicit number of molecules and box')
out = packmol(water, n_molecules=64, box_bounds=[0., 0., 0., 12., 13., 14.])
printsummary(out)
out.write('water-4.xyz')
show(out)


print('water-5.xyz: pure liquid in non-orthorhombic box (requires AMS2022 or later)')
# first place the molecules in a cuboid surrounding the desired lattice
# then gradually change into the desired lattice using refine_lattice()
# note that the molecules may become distorted by this procedure
lattice = [[10., 2., -1.], [-5., 8., 0.], [0., -2., 11.]]
temp_out = packmol(water, n_molecules=32, box_bounds=[
    0, 0, 0,
    max(lattice[i][0] for i in range(3))-min(lattice[i][0] for i in range(3)),
    max(lattice[i][1] for i in range(3))-min(lattice[i][1] for i in range(3)),
    max(lattice[i][2] for i in range(3))-min(lattice[i][2] for i in range(3))
])
out = refine_lattice(temp_out, lattice=lattice)
if out is not None:
    out.write('water-5.xyz')
    print('Top: system in surrounding orthorhombic box before calling refine_lattice(). Bottom: System in non-orthorhombic box after calling refine_lattice()')
    show(temp_out)
    show(out)


# ## Water-acetonitrile mixture (fluid with 2 or more components)
# Let's also create a single acetonitrile molecule:

acetonitrile = from_smiles('CC#N')
show(acetonitrile)


# Set the desired mole fractions and density. Here, the density is calculated as the weighted average of water (1.0 g/cm^3) and acetonitrile (0.76 g/cm^3) densities, but you could use any other density.

# MIXTURES
x_water = 0.666                # mole fraction
x_acetonitrile = 1-x_water     # mole fraction
density = (x_water*1.0 + x_acetonitrile*0.76) / (x_water + x_acetonitrile)  # weighted average of pure component densities

print(f'\nMIXTURES. x_water = {x_water:.3f}, x_acetonitrile = {x_acetonitrile:.3f}, target density = {density:.3f} g/cm^3\n')


# By setting ``return_details=True``, you can get information about the mole fractions of the returned system. They may not exactly match the mole fractions you put in.

print('2-1 water-acetonitrile from approximate number of atoms and exact density (in g/cm^3), cubic box with auto-determined size')
out, details = packmol(molecules=[water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], n_atoms=200, density=density, return_details=True)
printsummary(out, details)
out.write('water-acetonitrile-1.xyz')
show(out)


# The ``details`` is a dictionary as follows:

for k, v in details.items():
    print(f'{k}: {v}')


print('2-1 water-acetonitrile from approximate density (in g/cm^3) and box bounds')
out, details = packmol(molecules=[water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], box_bounds=[0, 0, 0, 13.2, 13.2, 13.2], density=density, return_details=True)
printsummary(out, details)
out.write('water-acetonitrile-2.xyz')
show(out)


print('2-1 water-acetonitrile from explicit number of molecules and density, cubic box with auto-determined size')
out, details = packmol(molecules=[water, acetonitrile], n_molecules=[32, 16], density=density, return_details=True)
printsummary(out, details)
out.write('water-acetonitrile-3.xyz')
show(out)


print('2-1 water-acetonitrile from explicit number of molecules and box')
out = packmol(molecules=[water, acetonitrile], n_molecules=[32, 16], box_bounds=[0, 0, 0, 13.2, 13.2, 13.2])
printsummary(out)
out.write('water-acetonitrile-4.xyz')
show(out)


# ## Solid-liquid or solid-gas interfaces
# First, create a slab using the ASE ``fcc111`` function

rotation = ('90x,0y,0z')  # sideview of slab
figsize = (3,3)
slab = fromASE(fcc111('Al', size=(4,6,3), vacuum=15.0, orthogonal=True, periodic=True))
show(slab, figsize=figsize, rotation=rotation)


print('water surrounding an Al slab, from an approximate density')
out = packmol_on_slab(slab, water, density=1.0)
printsummary(out)
out.write('al-water-pure.xyz')
show(out, figsize=figsize, rotation=rotation)


print('2-1 water-acetonitrile mixture surrounding an Al slab, from mole fractions and an approximate density')
out = packmol_on_slab(slab, [water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], density=density)
printsummary(out)
out.write('al-water-acetonitrile.xyz')
show(out, figsize=figsize, rotation=rotation)


# ## Microsolvation
# ``packmol_microsolvation`` can create a microsolvation sphere around a solute.

out = packmol_microsolvation(solute=acetonitrile, solvent=water, density=1.5, threshold=4.0)
# for microsolvation it's a good idea to have a higher density than normal to get enough solvent molecules
print(f"Microsolvated structure: {len(out)} atoms.")
out.write('acetonitrile-microsolvated.xyz')
show(out, figsize=figsize)


# ## Bonds, atom properties (force field types, regions, ...)
# 
# The ``packmol()`` function accepts the arguments ``keep_bonds`` and ``keep_atom_properties``. These options will keep the bonds defined for the constitutent molecules, as well as any atomic properties.
# 
# The bonds and atom properties are easiest to see by printing the System block for an AMS job:

water = from_smiles('O')
n2 = from_smiles('N#N')

# delete properties coming from from_smiles
for at in water:
    at.properties = Settings()
for at in n2:
    at.properties = Settings()
    
water[1].properties.region = "oxygen_atom"
water[2].properties.mass = 2.014   # deuterium
water.delete_bond(water[1, 2]) # delete bond between atoms 1 and 2 (O and H)


out = packmol([water, n2], n_molecules=[2, 1], density=0.5)
print(AMSJob(molecule=out).get_input())


# By default, the ``packmol()`` function assigns regions called ``mol0``, ``mol1``, etc. to the different added molecules. The ``region_names`` option lets you set custom names. 

out = packmol(
    [water, n2], n_molecules=[2, 1], density=0.5, 
    region_names=["water", "nitrogen_molecule"]
)
print(AMSJob(molecule=out).get_input())


# Below, we also set ``keep_atom_properties=False``, this will remove the previous regions (in this example "oxygen_atom") and mass. 

out = packmol(
    [water, n2], n_molecules=[2, 1], density=0.5, 
    keep_atom_properties=False
)
print(AMSJob(molecule=out).get_input())


# ``keep_bonds=False`` will additionally ignore any defined bonds:

out = packmol(
    [water, n2], n_molecules=[2, 1], density=0.5, 
    region_names=["water", "nitrogen_molecule"], 
    keep_bonds=False, keep_atom_properties=False
)
print(AMSJob(molecule=out).get_input())


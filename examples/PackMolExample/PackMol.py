#!/usr/bin/env amspython
# coding: utf-8

# ## Initial imports

from scm.plams import plot_molecule, from_smiles, Molecule
from scm.plams.interfaces.molecule.packmol import packmol
from ase.visualize.plot import plot_atoms
from ase.build import fcc111, bulk
import matplotlib.pyplot as plt
import importlib.metadata

AMS2025 = importlib.metadata.version("scm") >= "2024.201"
if AMS2025:
    from scm.plams import packmol_around


# ## Helper functions


def printsummary(mol, details=None):
    if details:
        density = details["density"]
    else:
        density = mol.get_density() * 1e-3
    s = f"{len(mol)} atoms, density = {density:.3f} g/cm^3"
    if mol.lattice:
        s += f", box = {mol.lattice[0][0]:.3f}, {mol.lattice[1][1]:.3f}, {mol.lattice[2][2]:.3f}"
    s += f", formula = {mol.get_formula()}"
    if details:
        s += f'\n#added molecules per species: {details["n_molecules"]}, mole fractions: {details["mole_fractions"]}'
    print(s)


# ## Liquid water (fluid with 1 component)
# First, create the gasphase molecule:

water = from_smiles("O")
plot_molecule(water)


print("pure liquid from approximate number of atoms and exact density (in g/cm^3), cubic box with auto-determined size")
out = packmol(water, n_atoms=194, density=1.0)
printsummary(out)
out.write("water-1.xyz")
plot_molecule(out)


print("pure liquid from approximate density (in g/cm^3) and an orthorhombic box")
out = packmol(water, density=1.0, box_bounds=[0.0, 0.0, 0.0, 8.0, 12.0, 14.0])
printsummary(out)
out.write("water-2.xyz")
plot_molecule(out)


print("pure liquid with explicit number of molecules and exact density")
out = packmol(water, n_molecules=64, density=1.0)
printsummary(out)
out.write("water-3.xyz")
plot_molecule(out)


print("pure liquid with explicit number of molecules and box")
out = packmol(water, n_molecules=64, box_bounds=[0.0, 0.0, 0.0, 12.0, 13.0, 14.0])
printsummary(out)
out.write("water-4.xyz")
plot_molecule(out)


if AMS2025:
    print("water-5.xyz: pure liquid in non-orthorhombic box (requires AMS2025 or later)")
    print("NOTE: Non-orthorhombic boxes may yield inaccurate results, always carefully check the output")
    # You can pack inside any lattice using the packmol_around function
    box = Molecule()
    box.lattice = [[10.0, 2.0, -1.0], [-5.0, 8.0, 0.0], [0.0, -2.0, 11.0]]
    out = packmol_around(box, molecules=[water], n_molecules=[32])
    out.write("water-5.xyz")
    plot_molecule(out)


if AMS2025:
    print("Experimental feature (AMS2025): guess density for pure liquid")
    print("Note: This density is meant to be equilibrated with NPT MD. It can be very inaccurate!")
    out = packmol(water, n_atoms=100)
    print(f"Guessed density: {out.get_density():.2f} kg/m^3")
    plot_molecule(out)


# ## Water-acetonitrile mixture (fluid with 2 or more components)
# Let's also create a single acetonitrile molecule:

acetonitrile = from_smiles("CC#N")
plot_molecule(acetonitrile)


# Set the desired mole fractions and density. Here, the density is calculated as the weighted average of water (1.0 g/cm^3) and acetonitrile (0.76 g/cm^3) densities, but you could use any other density.

# MIXTURES
x_water = 0.666  # mole fraction
x_acetonitrile = 1 - x_water  # mole fraction
# weighted average of pure component densities
density = (x_water * 1.0 + x_acetonitrile * 0.76) / (x_water + x_acetonitrile)

print("MIXTURES")
print(f"x_water = {x_water:.3f}")
print(f"x_acetonitrile = {x_acetonitrile:.3f}")
print(f"target density = {density:.3f} g/cm^3")


# By setting ``return_details=True``, you can get information about the mole fractions of the returned system. They may not exactly match the mole fractions you put in.

print(
    "2-1 water-acetonitrile from approximate number of atoms and exact density (in g/cm^3), "
    "cubic box with auto-determined size"
)
out, details = packmol(
    molecules=[water, acetonitrile],
    mole_fractions=[x_water, x_acetonitrile],
    n_atoms=200,
    density=density,
    return_details=True,
)
printsummary(out, details)
out.write("water-acetonitrile-1.xyz")
plot_molecule(out)


# The ``details`` is a dictionary as follows:

for k, v in details.items():
    print(f"{k}: {v}")


print("2-1 water-acetonitrile from approximate density (in g/cm^3) and box bounds")
out, details = packmol(
    molecules=[water, acetonitrile],
    mole_fractions=[x_water, x_acetonitrile],
    box_bounds=[0, 0, 0, 13.2, 13.2, 13.2],
    density=density,
    return_details=True,
)
printsummary(out, details)
out.write("water-acetonitrile-2.xyz")
plot_molecule(out)


print("2-1 water-acetonitrile from explicit number of molecules and density, cubic box with auto-determined size")
out, details = packmol(
    molecules=[water, acetonitrile],
    n_molecules=[32, 16],
    density=density,
    return_details=True,
)
printsummary(out, details)
out.write("water-acetonitrile-3.xyz")
plot_molecule(out)


print("2-1 water-acetonitrile from explicit number of molecules and box")
out = packmol(
    molecules=[water, acetonitrile],
    n_molecules=[32, 16],
    box_bounds=[0, 0, 0, 13.2, 13.2, 13.2],
)
printsummary(out)
out.write("water-acetonitrile-4.xyz")
plot_molecule(out)


if AMS2025:
    print("Experimental feature (AMS2025): guess density for mixture")
    print("Note: This density is meant to be equilibrated with NPT MD. It can be very inaccurate!")
    out = packmol([water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], n_atoms=100)
    print(f"Guessed density: {out.get_density():.2f} kg/m^3")
    plot_molecule(out)


# ## Pack inside sphere
#
# Set ``sphere=True`` to pack in a sphere (non-periodic) instead of in a periodic box. The sphere will be centered near the origin.

print("water in a sphere from exact density and number of molecules")
out, details = packmol(molecules=[water], n_molecules=[100], density=1.0, return_details=True, sphere=True)
printsummary(out, details)
print(f"Radius  of sphere: {details['radius']:.3f} ang.")
print(f"Center of mass xyz (ang): {out.get_center_of_mass()}")
out.write("water-sphere.xyz")
plot_molecule(out)


print(
    "2-1 water-acetonitrile in a sphere from exact density (in g/cm^3) and "
    "approximate number of atoms and mole fractions"
)
out, details = packmol(
    molecules=[water, acetonitrile],
    mole_fractions=[x_water, x_acetonitrile],
    n_atoms=500,
    density=density,
    return_details=True,
    sphere=True,
)
printsummary(out, details)
out.write("water-acetonitrile-sphere.xyz")
plot_molecule(out)


# ## Packing ions, total system charge
#
# The total system charge will be sum of the charges of the constituent molecules.
#
# In PLAMS, ``molecule.properties.charge`` specifies the charge:

ammonium = from_smiles("[NH4+]")  # ammonia.properties.charge == +1
chloride = from_smiles("[Cl-]")  # chloride.properties.charge == -1
print("3 water molecules, 3 ammonium, 1 chloride (non-periodic)")
print("Initial charges:")
print(f"Water: {water.properties.get('charge', 0)}")
print(f"Ammonium: {ammonium.properties.get('charge', 0)}")
print(f"Chloride: {chloride.properties.get('charge', 0)}")
out = packmol(molecules=[water, ammonium, chloride], n_molecules=[3, 3, 1], density=0.4, sphere=True)
tot_charge = out.properties.get("charge", 0)
print(f"Total charge of packmol-generated system: {tot_charge}")
out.write("water-ammonium-chloride.xyz")
plot_molecule(out)


# ## Microsolvation
# ``packmol_microsolvation`` can create a microsolvation sphere around a solute.

from scm.plams import packmol_microsolvation

out = packmol_microsolvation(solute=acetonitrile, solvent=water, density=1.5, threshold=4.0)
# for microsolvation it's a good idea to have a higher density than normal to get enough solvent molecules
print(f"Microsolvated structure: {len(out)} atoms.")
out.write("acetonitrile-microsolvated.xyz")

figsize = (3, 3)
plot_molecule(out, figsize=figsize)


# ## Solid-liquid or solid-gas interfaces
# First, create a slab using the ASE ``fcc111`` function

from scm.plams import plot_molecule, fromASE
from ase.build import fcc111

rotation = "90x,0y,0z"  # sideview of slab
slab = fromASE(fcc111("Al", size=(4, 6, 3), vacuum=15.0, orthogonal=True, periodic=True))
plot_molecule(slab, figsize=figsize, rotation=rotation)


print("water surrounding an Al slab, from an approximate density")
if AMS2025:
    out = packmol_around(slab, water, density=1.0)
    printsummary(out)
    out.write("al-water-pure.xyz")
    plot_molecule(out, figsize=figsize, rotation=rotation)


print("2-1 water-acetonitrile mixture surrounding an Al slab, from mole fractions and an approximate density")
if AMS2025:
    out = packmol_around(slab, [water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], density=density)
    printsummary(out)
    out.write("al-water-acetonitrile.xyz")
    plot_molecule(out, figsize=figsize, rotation=rotation)


from ase.build import surface

if AMS2025:
    print("water surrounding non-orthorhombic Au(211) slab, from an approximate number of molecules")
    print("NOTE: non-orthorhombic cell, results are approximate, requires AMS2025")
    slab = surface("Au", (2, 1, 1), 6)
    slab.center(vacuum=11.0, axis=2)
    slab.set_pbc(True)
    out = packmol_around(fromASE(slab), [water], n_molecules=[32], tolerance=1.8)
    out.write("Au211-water.xyz")
    plot_molecule(out, figsize=figsize, rotation=rotation)
    print(f"{out.lattice=}")


# ## Pack inside voids in crystals
#
# Use the ``packmol_around`` function. You can decrease ``tolerance`` if you need to pack very tightly. The default value for ``tolerance`` is 2.0.

from scm.plams import fromASE
from ase.build import bulk

bulk_Al = fromASE(bulk("Al", cubic=True).repeat((3, 3, 3)))
rotation = "-85x,5y,0z"
plot_molecule(bulk_Al, rotation=rotation, radii=0.4)


if AMS2025:
    out = packmol_around(
        current=bulk_Al,
        molecules=[from_smiles("[H]"), from_smiles("[He]")],
        n_molecules=[50, 20],
        tolerance=1.5,
    )
    plot_molecule(out, rotation=rotation, radii=0.4)
    printsummary(out)
    out.write("al-bulk-with-h-he.xyz")


# ## Bonds, atom properties (force field types, regions, ...)
#
# The ``packmol()`` function accepts the arguments ``keep_bonds`` and ``keep_atom_properties``. These options will keep the bonds defined for the constitutent molecules, as well as any atomic properties.
#
# The bonds and atom properties are easiest to see by printing the System block for an AMS job:

from scm.plams import Settings

water = from_smiles("O")
n2 = from_smiles("N#N")

# delete properties coming from from_smiles
for at in water:
    at.properties = Settings()
for at in n2:
    at.properties = Settings()

water[1].properties.region = "oxygen_atom"
water[2].properties.mass = 2.014  # deuterium
water.delete_bond(water[1, 2])  # delete bond between atoms 1 and 2 (O and H)


from scm.plams import AMSJob

out = packmol([water, n2], n_molecules=[2, 1], density=0.5)
print(AMSJob(molecule=out).get_input())


# By default, the ``packmol()`` function assigns regions called ``mol0``, ``mol1``, etc. to the different added molecules. The ``region_names`` option lets you set custom names.

out = packmol(
    [water, n2],
    n_molecules=[2, 1],
    density=0.5,
    region_names=["water", "nitrogen_molecule"],
)
print(AMSJob(molecule=out).get_input())


# Below, we also set ``keep_atom_properties=False``, this will remove the previous regions (in this example "oxygen_atom") and mass.

out = packmol([water, n2], n_molecules=[2, 1], density=0.5, keep_atom_properties=False)
print(AMSJob(molecule=out).get_input())


# ``keep_bonds=False`` will additionally ignore any defined bonds:

out = packmol(
    [water, n2],
    n_molecules=[2, 1],
    density=0.5,
    region_names=["water", "nitrogen_molecule"],
    keep_bonds=False,
    keep_atom_properties=False,
)
print(AMSJob(molecule=out).get_input())

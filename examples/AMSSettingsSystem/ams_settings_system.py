#!/usr/bin/env amspython
# coding: utf-8

# ## Initial imports

from scm.plams import *

# this line is not required in AMS2025+
init()


# ## Elements, coordinates, lattice vectors, and charge

# ### Manual molecule definition

molecule = Molecule()
molecule.add_atom(Atom(symbol="O", coords=(0, 0, 0)))
molecule.add_atom(Atom(symbol="H", coords=(1, 0, 0)))
molecule.add_atom(Atom(symbol="H", coords=(0, 1, 0)))


# To see the input that will be passed to AMS, create an AMSJob and print the input:

print(AMSJob(molecule=molecule).get_input())


# ### Lattice vectors: 1D-periodic
#
# For periodic systems in 1 dimension, the lattice vector must be along the x direction (with 0 components along y and z)

molecule.lattice = [[10, 0, 0]]
print(AMSJob(molecule=molecule).get_input())


# ### Lattice vectors: 2D-periodic
#
# For 2 dimensions, the two lattice vectors must lie in the xy plane (with 0 component along z).

molecule.lattice = [
    [10, 0, 0],
    [0, 11, 0],
]
print(AMSJob(molecule=molecule).get_input())


# ### Lattice vectors: 3D-periodic

molecule.lattice = [[10, 0, 0], [0, 11, 0], [-1, 0, 12]]
print(AMSJob(molecule=molecule).get_input())


# ### Delete lattice vectors

molecule.lattice = []
print(AMSJob(molecule=molecule).get_input())


# ### Charge

molecule.properties.charge = -1
print(AMSJob(molecule=molecule).get_input())


# To get the charge of a molecule, use ``molecule.properties.get("charge", 0)``. If the charge is not defined you will then get 0 as the charge.

my_charge = molecule.properties.get("charge", 0)
print(f"The charge is {my_charge}")


# Unset the charge:

if "charge" in molecule.properties:
    del molecule.properties.charge

my_charge = molecule.properties.get("charge", 0)
print(f"The charge is {my_charge}")


# ## Atomic properties: masses, regions, force field types ...
#
# In the AMS system block most atomic properties are given as a suffix at the end of the line.
#
# To access an individual atom, use for example ``molecule[1]``, which corresponds to the first atom. **Note that the indexing starts with 1**, unlike normal Python lists that start with 0!

# ### Isotopes (atomic masses)

molecule[2].properties.mass = 2.014
print(AMSJob(molecule=molecule).get_input())


# ### Regions
#
# Regions are used for example to
#
# * set special basis sets on a subset of atoms, or
# * apply a thermostat in molecular dynamics to only a subset of atoms,
# * visualize atoms easily in the AMS GUI,
# * and much more!
#
# Use Python sets to specify regions. In this way, one atom can belong to multiple regions.

molecule[1].properties.region = {"region1"}
molecule[2].properties.region = {"region1"}
molecule[3].properties.region = {"region1", "region2"}
print(AMSJob(molecule=molecule).get_input())


# ### Force field types
#
# Some force fields need to know the specific atom type and not just the chemical element. Use ``ForceField.Type`` for this when you use the ForceField engine:

molecule[1].properties.ForceField.Type = "OW"  # these types would depend on what type of force field you use!
molecule[2].properties.ForceField.Type = "HW"
molecule[3].properties.ForceField.Type = "HW"
print(AMSJob(molecule=molecule).get_input())


# ### Delete all atom-specific options
# Loop over the atoms and set ``atom.properties`` to an empty ``Settings()``:

for at in molecule:
    at.properties = Settings()

print(AMSJob(molecule=molecule).get_input())


# ## Bonds
#
# Most methods (DFT, DFTB, ML Potential, ReaxFF) ignore any specified bonds.
#
# When using force fields, you sometimes need to specify the bonds that connect atoms. Some force fields (UFF, GAFF) can automatically guess the correct types of bonds.
#
# So **most of the time you do not manually need to specify bonds**.
#
# If you **need** to specify bonds, it is easiest
#
# * to handle in the AMS GUI: use File -> Export Coordinates -> .in, and then load the file with ``molecule = Molecule("my_file.in")``
# * to use the ``from_smiles`` function to generate a molecule from SMILES code, for example ``molecule = from_smiles("O")`` for water.
#
# If you need to add bonds manually in PLAMS you can do it as follows:

molecule.add_bond(molecule[1], molecule[2], order=1.0)
molecule.add_bond(molecule[1], molecule[3], order=1.0)
print(AMSJob(molecule=molecule).get_input())


# ## Multiple systems
#
# Some tasks like NEB (nudged elastic band) require more than 1 system in the input file. This can be accomplished by using a Python dictionary.
#
# In AMS,
#
# * the "main system" has no name. It should have the key ``""`` (empty string) in the dictionary.
#
# * every additional system needs to have a name, that is used as the key in the dictionary.
#
# Let's first define two ``Molecule`` in the normal way:

molecule1 = Molecule()
molecule1.add_atom(Atom(symbol="O", coords=(0, 0, 0)))
molecule1.add_atom(Atom(symbol="H", coords=(1, 0, 0)))
molecule1.add_atom(Atom(symbol="H", coords=(0, 1, 0)))

molecule2 = Molecule()
molecule2.add_atom(Atom(symbol="O", coords=(0, 0, 0)))
molecule2.add_atom(Atom(symbol="H", coords=(3.33333, 0, 0)))
molecule2.add_atom(Atom(symbol="H", coords=(0, 5.555555, 0)))


# Then create the ``mol_dict`` dictionary:

mol_dict = {
    "": molecule1,  # main system, empty key (no name)
    "final": molecule2,  # for NEB, use "final" as the name for the other endpoint
}


# Pass the ``mol_dict`` as the ``molecule`` argument to ``AMSJob``:

print(AMSJob(molecule=mol_dict).get_input())


# Above we see that the main system is printed just as before. A second system block "system final" is also added with ``molecule2``.

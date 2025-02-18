#!/usr/bin/env amspython
# coding: utf-8

# This example shows how to perform a geometry optimization of a water molecule and compute
# the vibrational normal modes using GFN1-xTB.
#
# If you do not have
# a DFTB license, remove the line with DFTB settings and instead set
# ``settings.input.ForceField.Type = 'UFF'``

# ## Initial imports

from scm.plams import *

# this line is not required in AMS2025+
init()


# ## Initial structure

# You could also load the geometry from an xyz file:
# molecule = Molecule('path/my_molecule.xyz')
# or generate a molecule from SMILES:
# molecule = from_smiles('O')
molecule = Molecule()
molecule.add_atom(Atom(symbol="O", coords=(0, 0, 0)))
molecule.add_atom(Atom(symbol="H", coords=(1, 0, 0)))
molecule.add_atom(Atom(symbol="H", coords=(0, 1, 0)))


try:
    plot_molecule(molecule)  # plot molecule in a Jupyter Notebook in AMS2023+
except NameError:
    pass


# ## Calculation settings
#
# The calculation settings are stored in a ``Settings`` object, which is a type of nested dictionary.

settings = Settings()
settings.input.ams.Task = "GeometryOptimization"
settings.input.ams.Properties.NormalModes = "Yes"
settings.input.DFTB.Model = "GFN1-xTB"
# settings.input.ForceField.Type = 'UFF' # set this instead of DFTB if you do not have a DFTB license. You will then not be able to extract the HOMO and LUMO energies.


# ## Create an AMSJob

job = AMSJob(molecule=molecule, settings=settings, name="water_optimization")


# You can check the input to AMS by calling the ``get_input()`` method:

print("-- input to the job --")
print(job.get_input())
print("-- end of input --")


# ## Run the job

job.run()


# ## Main results files: ams.rkf and dftb.rkf
#
# The paths to the main binary results files ``ams.rkf`` and ``dftb.rkf`` can be retrieved as follows:

print(job.results.rkfpath(file="ams"))
print(job.results.rkfpath(file="engine"))


# ## Optimized coordinates

optimized_molecule = job.results.get_main_molecule()

print("Optimized coordinates")
print("---------------------")
print(optimized_molecule)
print("---------------------")


try:
    plot_molecule(optimized_molecule)  # plot molecule in a Jupyter Notebook in AMS2023+
except NameError:
    pass


# ## Optimized bond lengths and angle

# Unlike python lists, where the index of the first element is 0,
# the index of the first atom in the molecule object is 1.

bond_length = optimized_molecule[1].distance_to(optimized_molecule[2])
print("O-H bond length: {:.3f} angstrom".format(bond_length))


bond_angle = optimized_molecule[1].angle(optimized_molecule[2], optimized_molecule[3])
print("Bond angle  : {:.1f} degrees".format(Units.convert(bond_angle, "rad", "degree")))


# ## Calculation timing

timings = job.results.get_timings()

print("Timings")
print("-------")
for key, value in timings.items():
    print(f"{key:<20s}: {value:.3f} seconds")
print("-------")


# ## Energy

energy = job.results.get_energy(unit="kcal/mol")

print("Energy      : {:.3f} kcal/mol".format(energy))


# ## Vibrational frequencies

frequencies = job.results.get_frequencies(unit="cm^-1")

print("Frequencies")
print("-----------")
for freq in frequencies:
    print(f"{freq:.3f} cm^-1")
print("-----------")


# ##  Dipole moment

import numpy as np

try:
    dipole_moment = np.linalg.norm(np.array(job.results.get_dipolemoment()))
    dipole_moment *= Units.convert(1.0, "au", "debye")
    print("Dipole moment: {:.3f} debye".format(dipole_moment))
except KeyError:
    print("Couldn't extract the dipole moment")


# ## HOMO, LUMO, and HOMO-LUMO gap
#
# Note: The methods for extracting HOMO, LUMO, and HOMO-LUMO gap only exist in AMS2023 and later.

try:
    homo = job.results.get_homo_energies(unit="eV")[0]
    lumo = job.results.get_lumo_energies(unit="eV")[0]
    homo_lumo_gap = job.results.get_smallest_homo_lumo_gap(unit="eV")

    print("HOMO        : {:.3f} eV".format(homo))
    print("LUMO        : {:.3f} eV".format(lumo))
    print("HOMO-LUMO gap : {:.3f} eV".format(homo_lumo_gap))
except AttributeError:
    print("Methods to extract HOMO and LUMO require AMS2023 or later")
except KeyError:
    print("Couldn't extract the HOMO and LUMO.")


# ## Read results directly from binary .rkf files
#
# You can also read results directly from the binary .rkf files. Use the "expert mode" of the KFbrowser program that comes with AMS to find out which section and variable to read.
#
# Below, we show how to extract the ``AMSResults%Energy`` variable from the dftb.rkf file. This is the same number that was extracted previously using the ``job.results.get_energy()`` method.

energy = job.results.readrkf("AMSResults", "Energy", file="engine")
print(f"Energy from the engine .rkf file (in hartree): {energy}")

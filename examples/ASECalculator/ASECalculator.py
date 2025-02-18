#!/usr/bin/env amspython
# coding: utf-8

# ## Initial imports

from scm.plams import *
from scm.plams.interfaces.adfsuite.ase_calculator import AMSCalculator
from ase.optimize import BFGS
from ase.build import molecule as ase_build_molecule
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

# In this example AMS runs in AMSWorker mode, so we have no use for the PLAMS working directory
# Let's delete it after the calculations are done
config.erase_workdir = True

# this line is not required in AMS2025+
init()


# ## Construct an initial system
# Here, we use the ``molecule()`` from ``ase.build`` to construct an ASE Atoms object.
#
# You could also convert a PLAMS Molecule to the ASE format using ``toASE()``.

atoms = ase_build_molecule("CH3OH")
# alternatively:
# atoms = toASE(from_smiles('CO'))

atoms.set_pbc((True, True, True))  # 3D periodic
atoms.set_cell([4.0, 4.0, 4.0])  # cubic box

# plot the atoms
plt.figure(figsize=(2, 2))
plt.axis("off")
plot_atoms(atoms, scale=0.5)


# ## Set the AMS settings
#
# First, set the AMS settings as you normally would do:

s = Settings()
s.input.ams.Task = "SinglePoint"  # the geometry optimization is handled by ASE
s.input.ams.Properties.Gradients = "Yes"  # ensures the forces are returned
s.input.ams.Properties.StressTensor = "Yes"  # ensures the stress tensor is returned

# Engine definition, could also be used to set up ADF, ReaxFF, ...
s.input.ForceField.Type = "UFF"

# run in serial
s.runscript.nproc = 1


# ## Run the ASE optimizer

print("Initial coordinates:")
print(atoms.get_positions())

with AMSCalculator(settings=s, amsworker=True) as calc:
    atoms.calc = calc
    optimizer = BFGS(atoms)
    optimizer.run(fmax=0.27)  # optimize until forces are smaller than 0.27 eV/ang

print(f"Optimized energy (eV): {atoms.get_potential_energy()}")
print("Optimized coordinates:")
print(atoms.get_positions())

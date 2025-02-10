#!/usr/bin/env amspython
# coding: utf-8

# ## Initial Imports

import numpy as np
from scm.plams import Settings, Molecule, Atom, AMSJob, init

# this line is not required in AMS2025+
init()


# ## Setup Dimer
# Create Helium atoms and an array of interatomic distances at which to run calculation.

# type of atoms
atom1 = "He"
atom2 = "He"


# interatomic distance values
dmin = 2.2
dmax = 4.2
step = 0.2


# create a list with interatomic distances
distances = np.arange(dmin, dmax, step)
print(distances)


# ## Calculation Settings
#
# The calculation settins are stored in a `Settings` object.

# calculation parameters (single point, TZP/PBE+GrimmeD3)
sett = Settings()
sett.input.ams.task = "SinglePoint"
sett.input.adf.basis.type = "TZP"
sett.input.adf.xc.gga = "PBE"
sett.input.adf.xc.dispersion = "Grimme3"


# ## Create and Run Jobs
#
# For each interatomic distance, create a Helium dimer molecule with the required geometry then the single point energy calculation job. Run the job and extract the energy.

energies = []
for d in distances:
    mol = Molecule()
    mol.add_atom(Atom(symbol=atom1, coords=(0.0, 0.0, 0.0)))
    mol.add_atom(Atom(symbol=atom2, coords=(d, 0.0, 0.0)))
    job = AMSJob(molecule=mol, settings=sett, name=f"dist_{d:.2f}")
    job.run()
    energies.append(job.results.get_energy(unit="kcal/mol"))


# ## Results
#
# Print table of results of the distance against the calculated energy.

print("== Results ==")
print("d[A]    E[kcal/mol]")
for d, e in zip(distances, energies):
    print(f"{d:.2f}    {e:.3f}")


import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(3, 3))
ax.plot(distances, energies, ".-")
ax.set_xlabel("He-He distance (Å)")
ax.set_ylabel("Energy (kcal/mol)")

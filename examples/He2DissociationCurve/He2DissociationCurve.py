#!/usr/bin/env amspython
# coding: utf-8

# ## Initial Imports

import sys
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


# ## Calculation Settings
#
# The calculation settings are stored in a `Settings` object.

# calculation parameters (single point, TZP/PBE+GrimmeD3)
sett = Settings()
sett.input.ams.task = "SinglePoint"
sett.input.adf.basis.type = "TZP"
sett.input.adf.xc.gga = "PBE"
sett.input.adf.xc.dispersion = "Grimme3"


# ## Create and Run Jobs
#
# For each interatomic distance, create a Helium dimer molecule with the required geometry then the single point energy calculation job. Run the job and extract the energy.

jobs = []
for d in distances:
    mol = Molecule()
    mol.add_atom(Atom(symbol=atom1, coords=(0.0, 0.0, 0.0)))
    mol.add_atom(Atom(symbol=atom2, coords=(d, 0.0, 0.0)))
    job = AMSJob(molecule=mol, settings=sett, name=f"dist_{d:.2f}")
    jobs.append(job)
    job.run()


# ## Results
#
# Print table of results of the distance against the calculated energy.

print("== Results ==")
try:
    # For AMS2025+ can use JobAnalysis class to perform results analysis
    from scm.plams import JobAnalysis

    ja = (
        JobAnalysis(jobs=jobs, standard_fields=None)
        .add_field("Dist", lambda j: j.molecule[2].x, display_name="d[A]", fmt=".2f")
        .add_field("Energy", lambda j: j.results.get_energy(unit="kcal/mol"), display_name="E[kcal/mol]", fmt=".3f")
    )

    # Pretty-print if running in a notebook
    if "ipykernel" in sys.modules:
        ja.display_table()
    else:
        print(ja.to_table())

    energies = ja.Energy

except ImportError:

    energies = [j.results.get_energy(unit="kcal/mol") for j in jobs]

    print("d[A]    E[kcal/mol]")
    for d, e in zip(distances, energies):
        print(f"{d:.2f}    {e:.3f}")


import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(3, 3))
ax.plot(distances, energies, ".-")
ax.set_xlabel("He-He distance (Å)")
ax.set_ylabel("Energy (kcal/mol)")

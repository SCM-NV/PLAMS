#!/usr/bin/env amspython
# coding: utf-8

# Lets see how two alanine molecules orient themselves using CREST conformer generation.
# To do this we will constrain the system in a spherical region using the `SphericalWall` constraint.
# We start by setting up a system of two alanine molecules in a relatively small space.

# ## Initial imports

import scm.plams as plams
import sys
from scm.conformers import ConformersJob
from scm.conformers.plams.plot import plot_conformers
import numpy as np
import matplotlib.pyplot as plt
import os

# this line is not required in AMS2025+
plams.init()


# ## Single alanine molecule

smiles = "CC(N)C(=O)O"
alanine = plams.from_smiles(smiles)
plams.plot_molecule(alanine)


# ## Initial system: alanine dimer

# Pack two alanine molecules in a sphere with a density of 0.5 kg/L.

density = 0.5
mol = plams.packmol(alanine, n_molecules=2, density=density, sphere=True)


# Translate the molecule to be centered around the origin (needed for SphericalWall later):

mol.translate(-np.array(mol.get_center_of_mass()))


plams.plot_molecule(mol, rotation="0x,0y,90z")


# ## Calculation setup

# To determine the radius of the `SphericalWall` we measure the size of the initial dimer.

dists = plams.distance_array(mol, mol)
max_dist = np.max(dists)
diameter = 1.33 * max_dist
radius = diameter / 2
print(f"Largest distance between atoms: {max_dist:.3f} ang.")
print(f"Radius: {radius:.3f} ang.")


# Now we can set up the Crest conformer generation job, with the appropriate spherical wall constraining the molecules close together.

settings = plams.Settings()
settings.input.ams.EngineAddons.WallPotential.Enabled = "Yes"
settings.input.ams.EngineAddons.WallPotential.Radius = radius
settings.input.ams.Generator.Method = "CREST"
settings.input.ams.Output.KeepWorkDir = "Yes"
settings.input.ams.GeometryOptimization.MaxConvergenceTime = "High"
settings.input.ams.Generator.CREST.NCycles = 3  # at most 3 CREST cycles for this demo
settings.input.GFNFF = plams.Settings()


# ## Run the conformers job

# Now we can run the conformer generation job.

job = ConformersJob(molecule=mol, settings=settings)
job.run()
# ConformersJob.load_external("plams_workdir/conformers/conformers.rkf")  # load from disk instead of running the job


rkf = job.results.rkfpath()
print(f"Conformers stored in {rkf}")


# This job will run for approximately 15 minutes.

# ## Results
# Here we plot the three lowest-energy conformers.

plot_conformers(job)


# You can also open the conformers in AMSmovie to browse all conformers 1000+ conformers:


# Finally in AMS2025, you can also inspect the conformer data using the JobAnalysis tool.

try:
    from scm.plams import JobAnalysis

    ja = (
        JobAnalysis(standard_fields=None)
        .add_job(job)
        .add_field(
            "Id",
            lambda j: list(range(1, len(j.results.get_conformers()) + 1)),
            display_name="Conformer Id",
            expansion_depth=1,
        )
        .add_field(
            "Energies",
            lambda j: j.results.get_relative_energies("kcal/mol"),
            display_name="E",
            expansion_depth=1,
            fmt=".2f",
        )
        .add_field(
            "Populations",
            lambda j: j.results.get_boltzmann_distribution(298),
            display_name="P",
            expansion_depth=1,
            fmt=".3f",
        )
    )

    # Pretty-print if running in a notebook
    if "ipykernel" in sys.modules:
        ja.display_table(max_rows=20)
    else:
        print(ja.to_table())

except ImportError:
    pass

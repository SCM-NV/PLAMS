#!/usr/bin/env amspython

import sys

from scm.plams import AMSJob, Settings, init

# this line is not required in AMS2025+
init()

# Read from a previous result file
path = sys.argv[1]
job = AMSJob.load_external(path)
charges, types, patch = job.results.get_forcefield_params()
print(charges)
print(types)
print(patch)

# Add the atom types and charges to the molecule
assert job.molecule is not None
mol = job.molecule.copy()
for i, at in enumerate(mol.atoms):
    at.properties.ForceField.Charge = charges[i]
    at.properties.ForceField.Type = types[i]

# Write a patch file
with open("patch.dat", "w") as outfile:
    outfile.write(str(patch))
    outfile.close()

# Set up new job settings
settings = Settings()
settings.input.ForceField.Type = "GAFF"
settings.input.ForceField.ForceFieldPatchFile = "patch.dat"
settings.input.ams.Task = "SinglePoint"

# Create the new job, and write the input file
job = AMSJob(molecule=mol, settings=settings)
text = job.get_input()
with open("ams.in", "w") as outfile:
    outfile.write(text)
    outfile.close()

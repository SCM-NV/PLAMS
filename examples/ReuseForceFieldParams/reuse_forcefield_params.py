#!/usr/bin/env amspython
# coding: utf-8

# ## Initial Imports

import sys

from scm.plams import AMSJob, Settings, init, from_smiles

# this line is not required in AMS2025+
init()


# ## Run/Load Job with ForceField Information

# First run a reference calculation where charges are guessed (using DFTB by default):

ref_job = AMSJob.from_input("""
Task GeometryOptimization

GeometryOptimization
   Convergence Step=1.0e-3
End

System
   Atoms
      C 0.0 0.0 0.0
      O 1.13 0.0 0.0
      C 0.0 0.0 2.1
      O 1.13 0.0 1.9
   End
End

Engine ForceField
   Verbosity Verbose
   GuessCharges True
EndEngine
""")


ref_job.run()


# Alternatively, load a previously run calculation
# ref_job = AMSJob.load_external("./plams_workdir/plamsjob/ams.rkf")


# ## Reuse ForceField Parameters

# Extract the charges and types from the job results and add them as properties on the molecule:

charges, types, patch = ref_job.results.get_forcefield_params()


mol = ref_job.molecule[""].copy()

for i, at in enumerate(mol.atoms):
    at.properties.ForceField.Charge = charges[i]
    at.properties.ForceField.Type = types[i]


sett = Settings()
sett.input.AMS.Task = "SinglePoint"
sett.input.ForceField.Type = "UFF"


# Create a patch file if required
if patch:
    with open("patch.dat", "w") as outfile:
        outfile.write(str(patch))
        outfile.close()
    # For example with:
    # sett.input.ForceField.GAFF.ForceFieldPatchFile = "patch.dat"


job = AMSJob(molecule=mol, settings=sett)


print(job.get_input())


job.run()

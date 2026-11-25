#!/usr/bin/env amspython
# coding: utf-8

# ## Initialization

from scm.plams import Settings, Molecule, init, AMSJob
from scm.plams.recipes.adfnbo import ADFNBOJob

# this line is not required in AMS2025+
init()


# ## Define molecule


# mol = Molecule("methane.xyz")
def get_molecule(input_string):
    job = AMSJob.from_input(input_string)
    return job.molecule[""]


mol = get_molecule("""
System
    Atoms
         C      0.000000      0.000000      0.000000
         H      0.631600      0.631600      0.631600
         H      0.631600     -0.631600     -0.631600
         H     -0.631600      0.631600     -0.631600
         H     -0.631600     -0.631600      0.631600
    End
End
""")


# ## Create and run job

s = Settings()
s.input.AMS.Task = "SinglePoint"
s.input.ADF.basis.type = "DZP"
s.input.ADF.xc.lda = "SCF VWN"
s.input.ADF.relativity.level = "scalar"
s.adfnbo = ["write", "spherical", "fock"]

j = ADFNBOJob(molecule=mol, settings=s)
r = j.run()


# ## Print results

lines = r.get_output_chunk(begin="NATURAL BOND ORBITALS (Summary):", end="Charge unit", inc_begin=True, inc_end=True)
for line in lines:
    print(line)

#!/usr/bin/env amspython
# coding: utf-8

# ## Initialization

from scm.plams import *

# this line is not required in AMS2025+
init()


# ## Define molecule
# Normally you would read the molecule from an xyz file

# mol = Molecule('methanol.xyz')


mol = AMSJob.from_input("""
System
  Atoms
    C       0.000000000000       0.138569980000       0.355570700000
    O       0.000000000000       0.187935770000      -1.074466460000
    H       0.882876920000      -0.383123830000       0.697839450000
    H      -0.882876940000      -0.383123830000       0.697839450000
    H       0.000000000000       1.145042790000       0.750208830000
  End
End
""").molecule[""]


# ## Setup and run job
# Here we use the fast DFTB engine, for ADF it is recommended to disable symmetry

s = Settings()
s.input.ams.task = "SinglePoint"
s.input.ams.Properties.Gradients = "Yes"
# s.input.adf.basis.type = 'DZP'
# s.input.adf.symmetry = 'NOSYM'
# s.input.adf.xc.gga = 'PW91'
s.input.DFTB.Model = "GFN1-xTB"
s.runscript.nproc = 1

j = NumHessJob(name="test", molecule=mol, settings=s, jobtype=AMSJob, gradient=lambda x: x.get_gradients().reshape(-1))
r = j.run(jobrunner=JobRunner(parallel=True, maxjobs=8))


# ## Print results

print(r.get_hessian(mass_weighted=True))

#!/usr/bin/env amspython
import numpy as np
from scm.plams.recipes.numgrad import NumGradJob
from scm.plams import config, JobRunner, AMSJob, Settings, AMSResults, Units, init

# this line is not required in AMS2025+
init()

config.default_jobrunner = JobRunner(parallel=True, maxjobs=8)

mol = AMSJob.from_input(
    """
System
  Atoms
    C       0.000000000000       0.138569980000       0.355570700000
    O       0.000000000000       0.187935770000      -1.074466460000
    H       0.882876920000      -0.383123830000       0.697839450000
    H      -0.882876940000      -0.383123830000       0.697839450000
    H       0.000000000000       1.145042790000       0.750208830000
  End
End
"""
).molecule[""]


s = Settings()
s.input.AMS.Task = "SinglePoint"
s.input.ForceField.Type = "UFF"
s.runscript.nproc = 1

j = NumGradJob(name="numgrad", molecule=mol, settings=s, jobtype=AMSJob)
r = j.run()
r.wait()

s_analytical = s.copy()
s_analytical.input.AMS.Properties.Gradients = "Yes"
analytical_job = AMSJob(name="analytical", molecule=mol, settings=s_analytical)
analytical_job.run()

numerical_gradients = np.array(r.get_gradient(AMSResults.get_energy)).reshape(-1, 3)
analytical_gradients = np.array(analytical_job.results.get_gradients()).reshape(-1, 3) * Units.convert(
    1.0, "bohr^-1", "angstrom^-1"
)

print("Numerical Gradients (NumGradJob), hartree/angstrom:")
print(numerical_gradients)
print("Analytical Gradients, hartree/angstrom:")
print(analytical_gradients)

#!/usr/bin/env amspython
# coding: utf-8

# ## Initial Imports

import multiprocessing
from scm.plams import JobRunner, config, Settings, read_molecules, AMSJob, init

# this line is not required in AMS2025+
init()


# ## Configure Job Runner
#
# Set the default job runner to run in parallel. Run as many jobs simultaneously as there are cpu on the system. In addition, set the number of cores for each job to 1.

maxjobs = multiprocessing.cpu_count()
print("Running up to {} jobs in parallel simultaneously".format(maxjobs))


config.default_jobrunner = JobRunner(parallel=True, maxjobs=maxjobs)


config.job.runscript.nproc = 1


# ## Load Molecules
#
# Load set of molecules from directory containing xyz files.

molecules = read_molecules("molecules")


# ## Set Up and Run Jobs
#
# Configure the calculation settings in the `Settings` object. Run a geometry optimization job for each molecule in parallel.

settings = Settings()
settings.input.ams.Task = "GeometryOptimization"
settings.input.dftb.Model = "GFN1-xTB"


results = []
for name, molecule in sorted(molecules.items()):
    job = AMSJob(molecule=molecule, settings=settings, name=name)
    results.append(job.run())


# ## Results
#
# Print a table of results only for the successful calculations.

# Only print the results of the succesful caluclations:
for result in [r for r in results if r.ok()]:
    print("Energy for {:<12}: {:>10.3f} kcal/mol".format(result.name, result.get_energy(unit="kcal/mol")))

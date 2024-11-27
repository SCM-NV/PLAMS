import multiprocessing
from pathlib import Path

import scm.plams as plams


def folder_with_many_go_jobs():
    folder = Path("data_go_many-TEST")
    if folder.exists():
        return folder
    plams.init(folder=folder)
    maxjobs = 4 if multiprocessing.cpu_count() > 4 else 1
    plams.config.default_jobrunner = plams.JobRunner(parallel=True, maxjobs=maxjobs)
    plams.config.job.runscript.nproc = 1
    s = plams.Settings()
    s.input.ams.Task = "GeometryOptimization"
    s.input.ams.Properties.Gradients = True
    s.input.ams.GeometryOptimization.MaxIterations = 200
    s.input.forcefield
    smiles = ["CCO", "CCC", "O", "N"]
    atoms_dataset = {i: plams.from_smiles(i) for i in smiles}
    jobs = []
    for i, mol in atoms_dataset.items():
        job = plams.AMSJob(molecule=mol, settings=s, name=i)
        if i == "N":
            job.settings.input.ams.GeometryOptimization.MaxIterations = 100
        jobs.append(job)
        job.run()
    s.pop_nested("input.ams.GeometryOptimization")
    for i, mol in atoms_dataset.items():
        nsteps = 30
        if i == "N":
            nsteps = 20
        job = plams.AMSNVTJob(molecule=mol, samplingfreq=1, nsteps=nsteps, settings=s, name=i)
        jobs.append(job)
        job.run()

    plams.finish()
    return folder

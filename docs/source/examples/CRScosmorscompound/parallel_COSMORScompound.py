from scm.plams import from_smiles, JobRunner, config
from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob

config.default_jobrunner = JobRunner(parallel=True, maxjobs=8)  # Set the default jobrunner to be parallel
config.default_jobmanager.settings.hashing = None  # Disable rerun prevention
config.job.runscript.nproc = 1  # Number of cores for each job
config.log.stdout = 1  # Suppress plams output

rd_smiles = ["O", "CO"]
rd_names = ["H2O", "CO"]
molecules = {}
for name, smiles in zip(rd_names, rd_smiles):
    molecules[name] = from_smiles(smiles, nconfs=100, forcefield="uff")[0]  # lowest energy one in 100 conformers

results = []
for name, mol in molecules.items():
    job = ADFCOSMORSCompoundJob(
        molecule=mol,  # the initial structure
        coskf_name=name,  # a name to be used for coskf file
        coskf_dir="test_coskfs",  # a directory to put the .coskf files generated
        preoptimization="GFN1-xTB",  # perform preoptimize or not
        singlepoint=False,  # run a singlepoint in gasphase and solvation calculation without geometry optimization. Cannot be combined with ``preoptimization``
        name=name,
    )  # an optional name for the calculation directory
    results.append(job.run())

# preoptimization='GFN1-xTB', singlepoint=False => preoptimization(GFN1-xTB) -> optimization(ADF) -> COSMO
# preoptimization= None     , singlepoint=False => no preoptimization        -> optimization(ADF) -> COSMO
# preoptimization= None     , singlepoint=True  => no preoptimization        -> single point(ADF) -> COSMO

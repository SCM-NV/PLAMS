#!/usr/bin/env amspython
# coding: utf-8

# ## Generating coskf files from xyz

# The example will first load all the molecules in the folder ``compounds_xyz`` and then optimize the gas geometry using ADF, and perform the ADF COSMO calculation for each compound. When the calculations are finished, we will find all the .coskf files in the ``test_coskfs_xyz`` directory.

from scm.plams import from_smiles, read_molecules, init, JobRunner, config
from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob
import os

# this line is not required in AMS2025+
init()


# Enable the parallel calculation through `JobRunner`. Here, we'll assign one core to each job, and we can have up to eight jobs running all at once.

config.default_jobrunner = JobRunner(parallel=True, maxjobs=8)  # Set the default jobrunner to be parallel
config.default_jobmanager.settings.hashing = None  # Disable rerun prevention
config.job.runscript.nproc = 1  # Number of cores for each job
config.log.stdout = 1  # Suppress plams output


# Optional compound metadata keyed by molecule name. This metadata is passed to ADFCOSMORSCompoundJob via `mol_info`, written to the "Compound Data" section of the generated COSKF file, and used by AMSCrs GUI and `pyCRS.Database` to search compounds.

mol_info_by_name = {
    "CO": {"CAS": "630-08-0", "IUPAC": "Carbon monoxide", "Other Name": ""},
    "H2O": {"CAS": "7732-18-5", "IUPAC": "Water", "Other Name": ""},
}


molecules = read_molecules("./compounds_xyz")

results = []
for name, mol in molecules.items():
    mol_info = mol_info_by_name.get(name)
    job = ADFCOSMORSCompoundJob(molecule=mol, coskf_name=name, coskf_dir="test_coskfs_xyz", mol_info=mol_info)
    results.append(job.run())


for result in results:
    result.wait()


print(f"coskf files generated: {', '.join([f for f in os.listdir('./test_coskfs_xyz')])}")


# ## Generating .coskf files from smiles

# Now, we will specify the smiles and name of a set of compounds and generate the initial geometry of each compound using `from_smiles` function. With the setting, `nconfs=100` and `forcefield='uff'`, we will generate 100 conformers and find the one with the lowest energy using 'uff' forcefield. When the calculations are finished, we will find all the .coskf file in the ``test_coskfs_smiles`` directory.

rd_smiles = ["O", "CO"]
rd_names = ["H2O", "CO"]
molecules = {}
for name, smiles in zip(rd_names, rd_smiles):
    molecules[name] = from_smiles(smiles, nconfs=100, forcefield="uff")[0]  # lowest energy one in 100 conformers


# Lastly, we give this information to the `ADFCOSMORSCompoundJob` class, including the name of the coskf files as well as the directory in which we'll find them after the calculations complete.  Using the setting, `preoptimization='GFN1-xTB'` and `singlepoint=False`, it will utilize the DFTB for a quick pre-optimization. Subsequently, it will execute a gas phase optimization using ADF, followed by the solvation calculation.

results = []
for name, mol in molecules.items():
    mol_info = mol_info_by_name.get(name)
    job = ADFCOSMORSCompoundJob(
        molecule=mol,  # The initial structure
        coskf_name=name,  # a name to be used for coskf file
        coskf_dir="test_coskfs_smiles",  # a directory to put the .coskf files generated
        preoptimization="GFN1-xTB",  # perform preoptimize or not
        singlepoint=False,  # run a singlepoint in gasphase and solvation calculation without geometry optimization. Cannot be combined with `preoptimization`
        name=name,  # an optional name for the calculation directory
        mol_info=mol_info,  # compound information to be stored in the Compound Data section of the COSKF file
    )
    results.append(job.run())


for result in results:
    result.wait()


print(f"coskf files generated: {', '.join([f for f in os.listdir('./test_coskfs_smiles')])}")

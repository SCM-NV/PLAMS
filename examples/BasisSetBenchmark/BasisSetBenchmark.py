#!/usr/bin/env amspython
# coding: utf-8

# ## Initial Imports

import sys
import multiprocessing
from scm.plams import JobRunner, config, from_smiles, Settings, AMSJob, init
import numpy as np

# this line is not required in AMS2025+
init()


# ## Set Up Job Runner
# Set up job runner, running as many jobs as possible in parallel.

config.default_jobrunner = JobRunner(parallel=True, maxjobs=multiprocessing.cpu_count())
config.job.runscript.nproc = 1


# ## Set Up Molecules
# Create the molecules we want to use in our benchmark from SMILES.

# The molecules we want to use in our benchmark:
mol_smiles = {"Methane": "C", "Ethane": "C-C", "Ethylene": "C=C", "Acetylene": "C#C"}
molecules = {}
for name, smiles in mol_smiles.items():
    # Compute 10 conformers, optimize with UFF and pick the lowest in energy.
    molecules[name] = from_smiles(smiles, nconfs=10, forcefield="uff")[0]
    print(name, molecules[name])


# ## Initialize Calculation Settings
# Set up the settings which are common across jobs. The basis type is added later for each job.

common_settings = Settings()
common_settings.input.ams.Task = "SinglePoint"
common_settings.input.ams.System.Symmetrize = "Yes"
common_settings.input.adf.Basis.Core = "None"


basis = ["QZ4P", "TZ2P", "TZP", "DZP", "DZ", "SZ"]
reference_basis = "QZ4P"


# ## Run Calculations

results = {}
jobs = []
for bas in basis:
    for name, molecule in molecules.items():
        settings = common_settings.copy()
        settings.input.adf.Basis.Type = bas
        job = AMSJob(name=name + "_" + bas, molecule=molecule, settings=settings)
        jobs.append(job)
        results[(name, bas)] = job.run()


# ## Results
# Extract the energy from each calculation. Calculate the average absolute error in bond energy per atom for each basis set.

try:
    # For AMS2025+ can use JobAnalysis class to perform results analysis
    from scm.plams import JobAnalysis

    ja = (
        JobAnalysis(jobs=jobs, standard_fields=["Formula", "Smiles"])
        .add_settings_field(("Input", "ADF", "Basis", "Type"), display_name="Basis")
        .add_field("NAtoms", lambda j: len(j.molecule))
        .add_field(
            "Energy", lambda j: j.results.get_energy(unit="kcal/mol"), display_name="Energy [kcal/mol]", fmt=".2f"
        )
        .sort_jobs(["NAtoms", "Energy"])
    )

    ref_ja = ja.copy().filter_jobs(lambda data: data["InputAdfBasisType"] == "QZ4P")

    ref_energies = {f: e for f, e in zip(ref_ja.Formula, ref_ja.Energy)}

    def get_average_error(job):
        return abs(job.results.get_energy(unit="kcal/mol") - ref_energies[job.molecule.get_formula()]) / len(
            job.molecule
        )

    ja.add_field("AvErr", get_average_error, display_name="Average Error [kcal/mol]", fmt=".2f")

    # Pretty-print if running in a notebook
    if "ipykernel" in sys.modules:
        ja.display_table()
    else:
        print(ja.to_table())

except ImportError:

    average_errors = {}
    for bas in basis:
        if bas != reference_basis:
            errors = []
            for name, molecule in molecules.items():
                reference_energy = results[(name, reference_basis)].get_energy(unit="kcal/mol")
                energy = results[(name, bas)].get_energy(unit="kcal/mol")
                errors.append(abs(energy - reference_energy) / len(molecule))
                print("Energy for {} using {} basis set: {} [kcal/mol]".format(name, bas, energy))
            average_errors[bas] = sum(errors) / len(errors)


print("== Results ==")
print("Average absolute error in bond energy per atom")
for bas in basis:
    if bas != reference_basis:
        if ja:
            av = np.average(ja.copy().filter_jobs(lambda data: data["InputAdfBasisType"] == bas).AvErr)
        else:
            av = average_errors[bas]
        print("Error for basis set {:<4}: {:>10.3f} [kcal/mol]".format(bas, av))

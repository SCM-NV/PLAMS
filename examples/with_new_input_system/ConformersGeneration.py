#!/usr/bin/env amspython
from scm.conformers import ConformersJob, ConformersResults
from scm.inputs import DFTB, Conformers
from scm.plams import *

# This example shows how to use the AMS's Conformers tool via PLAMS


def print_results(results: ConformersResults):
    for i, (energy, mol) in enumerate(zip(results.get_energies("kcal/mol"), results.get_conformers())):
        print(f"Energy conformer {i} = {energy} [kcal/mol]")
        print(mol)


ethanol = from_smiles("OCC")

# Simple example: generate conformers for ethanol using default settings
# - default Task: Generate
# - default Generator Method: RDKit
# - default Engine: ForceField (with UFF)

job = ConformersJob(name="conformers_generation", molecule=ethanol)
job.run()

print("Conformers generated using RDKit and UFF:")
print_results(job.results)

# Re-optimize the conformers generated in the previous steps using the GFN1-xTB engine:

sett = Settings()

# In the 'ams' part of the settings input you can specify
# the input options for the Conformers tool, which are
# described in the Conformers tool user manual.

sett.input = Conformers()
sett.input.Task = "Optimize"
sett.input.InputConformersSet = [job.results.rkfpath()]

# You can specify the engine to be used (and the engine options) like you would
# for an AMSJob. See the AMSJob documentation for more details.

sett.input.Engine = DFTB()
sett.input.Engine.Model = "GFN1-xTB"

# Note: here we do not specify the input molecule because we are passing the results
# of a previous ConformersJob (ConformersResults) via the "InputConformersSet" input.

reoptimize_job = ConformersJob(name="optimize_conformers", settings=sett)
print(reoptimize_job.get_input())
reoptimize_job.run()

print("Conformers re-optimized using the more accurate GFN1-xTB method:")
print_results(reoptimize_job.results)

temperature = 273
print(f"Boltzmann distribution at {temperature} Kelvin")
for i, p in enumerate(reoptimize_job.results.get_boltzmann_distribution(temperature)):
    print(f"Conformer {i} probability: {p}")

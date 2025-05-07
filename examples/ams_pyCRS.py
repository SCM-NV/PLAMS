#!/usr/bin/env amspython
"""
- Solubility Calculation using pyCRS workflow with AMS2024+ and AMS2025+ Features
- Utilizes COSKFDatabase and CRSSystem from pyCRS (AMS2024+ Required)
- Use of hbc_from_MESP in ADFCOSMORSCompoundJob for COSMO-SAC DHB-MESP method (AMS2025+ Required)
"""

import os
import matplotlib.pyplot as plt
from pyCRS.Database import COSKFDatabase
from pyCRS.CRSManager import CRSSystem
from scm.plams import Settings, from_smiles, CRSJob, CRSResults
from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob


def solubility_pyCRS():

    db = COSKFDatabase("my_coskf_db.db")

    if db.get_compounds_id("water")[0] is None:
        db.add_compound("Water.coskf")

    if db.get_compounds_id("benzene")[0] is None:
        solute_smiles = "c1ccccc1"

        # Define molecular information for Benzene to be stored in the "Compound Data" section of the *.coskf* file
        mol_info = {}
        mol_info["IUPAC"] = "Benzene"
        mol_info["Other Name"] = None
        mol_info["CAS"] = "71-43-2"
        mol_info["SMILES"] = solute_smiles

        # Generate COSKF file with hbc_from_MESP for DHB-MESP in COSMO-SAC
        coskf_file = generate_coskf(solute_smiles, jobname="adf_benzene_densf", mol_info=mol_info, hbc_from_MESP=True)

        db.add_compound(coskf_file)
        db.add_physical_property("benzene", "meltingpoint", 278.7)
        db.add_physical_property("benzene", "hfusion", 9.91, unit="kJ/mol")

    plot_sigma_profile_pyCRS(name="benzene", db=db)

    crs = CRSSystem()

    mixture = {}
    mixture["water"] = 1.0
    mixture["benzene"] = 0.0
    temp = "273.15 283.15 10"

    solvent_density = 1.0
    additional_sett = Settings()
    additional_sett.input.property.DensitySolvent = solvent_density

    crs.add_Mixture(
        mixture=mixture,
        temperature=temp,
        database="my_coskf_db.db",
        problem_type="solubility",
        additional_sett=additional_sett,
    )
    crs.runCRSJob()

    plot_results(crs.outputs[0])


def generate_coskf(smiles, jobname=None, mol_info=None, hbc_from_MESP=False):
    molecule = from_smiles(smiles, nconfs=100, forcefield="uff")[0]
    job = ADFCOSMORSCompoundJob(name=jobname, molecule=molecule, mol_info=mol_info, hbc_from_MESP=hbc_from_MESP)
    job.run()
    return job.results.coskfpath()


def plot_results(job):
    if isinstance(job, CRSJob):
        res = job.results.get_results("SOLUBILITY")
    elif isinstance(job, CRSResults):
        res = job.get_results("SOLUBILITY")
    solubility_g_L = res["solubility g_per_L_solvent"][1]
    temperatures = res["temperature"]
    for temperature, sol_g_l in zip(temperatures, solubility_g_L):
        print(f"{temperature:.2f} {sol_g_l:.4f}")

    plt.plot(temperatures, solubility_g_L)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Solubility (g/L solvent)")
    plt.show()


def plot_sigma_profile_pyCRS(name, db):

    crs = CRSSystem()
    mixture = {name: 1.0}

    crs.add_Mixture(mixture=mixture, problem_type="PURESIGMAPROFILE", database=db)
    crs.runCRSJob()
    sigma = crs.outputs[0].get_sigma_profile()
    xlabel = "σ (e/A**2)"
    for profile in sigma:
        if profile == xlabel:
            continue
        plt.plot(sigma[xlabel], sigma[profile], label=profile.split(".")[0])

    plt.xlabel("σ (e/Å**2)")
    plt.ylabel("p(σ)")
    plt.legend()
    plt.show()


def main():
    solubility_pyCRS()


if __name__ == "__main__":
    main()

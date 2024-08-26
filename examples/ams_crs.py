#!/usr/bin/env amspython
import os

import matplotlib.pyplot as plt
from scm.plams import Settings, Units, from_smiles, CRSJob, init, finish
from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob


def solubility():
    # database can also be replaced with the output of "$AMSBIN/amspackages loc adfcrs" /ADFCRS-2018
    database = CRSJob.database()
    if not os.path.exists(database):
        raise OSError(f"The provided path does not exist. Exiting.")

    solute_smiles = "c1ccccc1"
    solute_coskf = generate_coskf(solute_smiles, "adf_benzene")  # generate files with ADF
    # solute_coskf = os.path.abspath('plams_workdir/adf_benzene/adf_benzene.coskf') # to not rerun the ADF calculation
    # solute_coskf = os.path.join(database, 'Benzene.coskf') # to load from database

    # You can also estimate the solute properties with the Property Prediction tool. See the Property Prediction example
    solute_properties = {"meltingpoint": 278.7, "hfusion": 9.91}  # experimental values for benzene, hfusion in kJ/mol

    solvent_coskf = os.path.join(database, "Water.coskf")
    solvent_density = 1.0

    s = Settings()
    s.input.property._h = "solubility"
    s.input.property.DensitySolvent = solvent_density
    s.input.temperature = "273.15 283.15 10"
    s.input.pressure = "1.01325 1.01325 10"

    s.input.compound = [Settings(), Settings()]

    s.input.compound[0]._h = solvent_coskf
    s.input.compound[0].frac1 = 1.0

    s.input.compound[1]._h = solute_coskf
    s.input.compound[1].nring = 6  # number of ring atoms benzene
    s.input.compound[1].meltingpoint = solute_properties["meltingpoint"]
    s.input.compound[1].hfusion = solute_properties["hfusion"] * Units.convert(
        1.0, "kJ/mol", "kcal/mol"
    )  # convert from kJ/mol to kcal/mol

    job = CRSJob(name="benzene_in_water", settings=s)
    job.run()

    plot_results(job)


def generate_coskf(smiles, jobname=None):
    molecule = from_smiles(smiles, nconfs=100, forcefield="uff")[0]
    job = ADFCOSMORSCompoundJob(name=jobname, molecule=molecule)
    job.run()
    plot_sigma_profile(job.results)
    return job.results.coskfpath()


def plot_results(job):
    res = job.results.get_results("SOLUBILITY")
    solubility_g_L = res["solubility g_per_L_solvent"][1]
    temperatures = res["temperature"]
    for temperature, sol_g_l in zip(temperatures, solubility_g_L):
        print(f"{temperature:.2f} {sol_g_l:.4f}")

    plt.plot(temperatures, solubility_g_L)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Solubility (g/L solvent)")
    plt.show()


def get_sigma_profile(coskf_file):
    s = Settings()
    s.input.property._h = "PURESIGMAPROFILE"
    s.input.compound._h = coskf_file
    job = CRSJob(name="sigma_profile", settings=s)
    res = job.run()
    return res.get_sigma_profile()


def plot_sigma_profile(results):
    coskf_path = results.coskfpath()
    sigma = get_sigma_profile(coskf_path)
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
    solubility()


if __name__ == "__main__":
    init()
    main()
    finish()

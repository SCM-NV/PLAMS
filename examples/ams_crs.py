#!/usr/bin/env amspython
from scm.plams import *
import numpy as np
import os
import matplotlib.pyplot as plt


def solubility():
    # replace with the output of "$AMSBIN/amspackages loc adfcrs" /ADFCRS-2018
    database = CRSJob.database()
    
    solute_smiles = 'c1ccccc1'
    solute_coskf = generate_coskf(solute_smiles, 'adf_benzene') # generate files with ADF 
    #solute_coskf = os.path.abspath('plams_workdir/adf_benzene/adf_benzene.coskf') # to not rerun the ADF calculation
    #solute_coskf = os.path.join(database, 'Benzene.coskf') # to load from database

    solute_properties = PropertyPrediction(solute_smiles).results # estimate with the property prediction tool
    #solute_properties = { 'meltingpoint': 278.7, 'hfusion': 9.91  } #experimental values for benzene, hfusion in kJ/mol

    solvent_coskf = os.path.join(database, 'Water.coskf')
    solvent_density = 1.0

    s = Settings()
    s.input.property._h = 'solubility'
    s.input.property.DensitySolvent = solvent_density
    s.input.temperature = "273.15 283.15 10"
    s.input.pressure = "1.01325 1.01325 10"

    s.input.compound = [Settings(), Settings()]

    s.input.compound[0]._h = solvent_coskf
    s.input.compound[0].frac1 = 1.0

    s.input.compound[1]._h = solute_coskf
    s.input.compound[1].meltingpoint = solute_properties['meltingpoint']
    s.input.compound[1].hfusion = solute_properties['hfusion'] * Units.convert(1.0, 'kJ/mol', 'kcal/mol') # convert from kJ/mol to kcal/mol

    job = CRSJob(name='benzene_in_water', settings=s)
    job.run()

    plot_results(job)

def generate_coskf(smiles, jobname=None):
    molecule = from_smiles(smiles, nconfs=100, forcefield='uff')[0]
    job = ADFCOSMORSCompoundJob(name=jobname, molecule=molecule)
    job.run()
    return job.results.coskfpath()

def plot_results(job):
    res = job.results.get_results('SOLUBILITY')
    solubility_g_L = res['solubility g_per_L_solvent'][1]
    temperatures = res['temperature']
    for temperature, sol_g_l in zip(temperatures, solubility_g_L):
        print(f'{temperature:.2f} {sol_g_l:.4f}')

    plt.plot(temperatures, solubility_g_L)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Solubility (g/L solvent)")
    plt.show()


def main():
    solubility()

if __name__ == '__main__':
    init()
    main()
    finish()


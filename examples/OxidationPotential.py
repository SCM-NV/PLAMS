from scm.plams import *
import os, sys

job_dir = './Screening'
if not os.path.exists(job_dir):
    os.makedirs(job_dir)

COSMORS_solvent_path = os.path.abspath('Dichloromethane.coskf')

results = {}
for mol_file in ["./molecules/NDI.xyz", "./molecules/NDI44.xyz", "./molecules/NDI55.xyz", "./molecules/NDI54.xyz", "./molecules/PDI.xyz"]:
    for method in ['screening']:
        job_name = None #if set to None, a name will be generated
    
        if job_name is None:
            job_name = os.path.basename(mol_file).split('.')[0] + '_' + method

        #calculation part
        init(path=job_dir, folder=job_name)

        workdir = config.default_jobmanager.workdir
        logfile = open(f'{workdir}/{job_name}_python.log', 'w')

        ox_potential_calc = OxidationPotentialCalculator(logfile=logfile)
        mol = Molecule(mol_file)
        oxpot = ox_potential_calc(mol, job_dir=workdir, method=method, COSMORS_solvent_path=COSMORS_solvent_path)
        results[job_name] = oxpot
        finish()

name_len = max(len(n) for n in results)
print(f'\n{"System".ljust(name_len)} | Oxidation potential')
for n, o in results.items():
    print(f'{n:{name_len}} | {o:.3f} eV')

print('\nCalculations coomplete!\a')
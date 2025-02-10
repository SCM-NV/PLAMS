#!/usr/bin/env amspython
# coding: utf-8

# ## Initial Imports

from scm.plams import Settings, ADFFragmentResults, Molecule, log, ADFFragmentJob, add_to_class, init

# this line is not required in AMS2025+
init()


# ## Helper Functions
# Add helper results extraction method.


@add_to_class(ADFFragmentResults)
def get_transfer_integrals(self):
    return self.job.full.results.read_rkf_section("TransferIntegrals", file="adf")


# ## Configure Settings
# Set up common settings for all 3 jobs, the specific settings for full system job. Note that in the interest of computational speed we use a minimal basis set. For more quantitatively meaningful results, you should use a larger basis.

common = Settings()
common.input.ams.Task = "SinglePoint"
common.input.adf.Basis.Type = "SZ"
common.input.adf.Basis.Core = "None"
common.input.adf.Symmetry = "NoSym"


full = Settings()
full.input.adf.transferintegrals = True


# ## Load Molecule
# Load benzene dimer from xyz file and separate it into 2 fragments. Alternatively one could simply load fragments from separate xyz files:
# ```
# mol1 = Molecule('fragment1.xyz')
# mol2 = Molecule('fragment2.xyz')
# ```

mol = Molecule("BenzeneDimer.xyz")
mol.guess_bonds()
fragments = mol.separate()
if len(fragments) != 2:
    log("ERROR: Molecule {} was split into {} fragments".format(mol.properties.name, len(fragments)))
    import sys

    sys.exit(1)
else:
    mol1, mol2 = fragments


# ## Run Job and Get Results

job = ADFFragmentJob(name="ADFTI", fragment1=mol1, fragment2=mol2, settings=common, full_settings=full)
results = job.run()


# TI is a dictionary with the whole TransferIntegrals section from adf.rkf
print("== Results ==")
TI = results.get_transfer_integrals()
for key, value in sorted(TI.items()):
    print("{:<28}: {:>12.6f}".format(key, value))

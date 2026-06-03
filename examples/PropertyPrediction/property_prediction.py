#!/usr/bin/env amspython
# coding: utf-8

# ## Initial imports

import pyCRS
import matplotlib.pyplot as plt
import math

import scm.plams as plams

try:
    from scm.plams import view  # view molecule using AMSview in a Jupyter Notebook in AMS2026+
except ImportError:
    from scm.plams import plot_molecule  # plot molecule in a Jupyter Notebook in AMS2023+

    def view(molecule, **kwargs):
        plot_molecule(molecule)


# ## Property prediction from SMILES (ethyl acetate)

from scm.plams import from_smiles

# smiles = 'CCO' # ethanol
smiles = "O=C(OCC)C"  # ethyl acetate
view(plams.from_smiles(smiles), width=300, height=300)


# ### Temperature-independent properties

print(f"SMILES: {smiles}\n")
mol = pyCRS.Input.read_smiles(smiles)

temperatures = [298.15, 308.15, 318.15, 328.15, 338.15]
pyCRS.PropPred.estimate(mol, temperatures=temperatures)

for prop, value in mol.properties.items():
    unit = pyCRS.PropPred.units[prop]
    print(f"{prop:<20s}: {value:.3f} {unit}")


# ### Temperature-dependent properties (vapor pressure)

prop = "vaporpressure"
unit = pyCRS.PropPred.units[prop]
temperatures_K, vaporpressures = mol.get_tdep_values(prop)
temperatures_C = [t - 273.15 for t in temperatures_K]  # convert to Celsius

plt.figure(figsize=(3, 3))
plt.plot(temperatures_C, vaporpressures)
plt.plot(temperatures_C, vaporpressures, ".")
plt.xlabel("Temperature (degree Celsius)")
plt.title(f"SMILES: {smiles}")
plt.ylabel(f"{prop} [{unit}]")


# ## Create .csv for multiple compounds
#
# Define a list of compounds by their SMILES strings. This example also shows how to only calculate a subset of all properties.
#
# Note: The SMILES string 'C' corresponds to methane which is too small to be used with the property prediction tool, so the results are given as 'nan' (not a number).

smiles_list = [
    "CCO",
    "CCOC",
    "OCCCN",
    "C",  # methane is too small to be used with property prediction and will return "nan"
    "C1=CC=C(C=C1)COCC2=CC=CC=C2",
]
temperatures = list(range(280, 340, 10))

mols = [pyCRS.Input.read_smiles(s) for s in smiles_list]

properties = ["boilingpoint", "criticaltemp", "hformstd"]

for mol in mols:
    pyCRS.PropPred.estimate(mol, properties, temperatures=temperatures)


def get_csv(mols, properties):
    header = "SMILES"
    for prop in properties:
        unit = pyCRS.PropPred.units[prop]
        if unit:
            unit = f" [{unit}]"
        else:
            unit = ""

        header += f",{prop}{unit}"
    ret = header + "\n"

    for mol in mols:
        s = f"{mol.smiles}"
        for prop in properties:
            value = mol.properties.get(prop, "")
            try:
                s += f",{value:.4f}"
            except TypeError:
                s += f",{value}"
        s += "\n"
        ret += s
    return ret


csv = get_csv(mols, properties)
print(csv)

# To write to a .csv file:
# with open('outputfile.csv', 'w') as f:
#    f.write(csv)


# ### Bar chart for multiple compounds
#
# Continuing from the previous example, you can also create e.g. a bar chart with the boiling points:

prop = "boilingpoint"
values = [mol.properties.get(prop, None) for mol in mols]
filtered = [(s, v) for s, v in zip(smiles_list, values) if not math.isnan(v)]
smiles_plot, values_plot = zip(*filtered)

plt.barh(smiles_plot, values_plot)
plt.title("Boiling point [K]")

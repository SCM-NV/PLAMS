#!/usr/bin/env amspython
# coding: utf-8

# ## Test 1. Balance simple reactions
# Here we present eleven different sets of reactants and products, not necessarily based on real reactions.
# The molecules can be passed as formulas, or as PLAMS Molecule objects.

from scm.plams import ReactionEquation
from scm.plams import from_smiles


# ### Aspirin in water

reactants = ["C9H8O4", from_smiles("O")]
products = ["C2H4O2", "C7H6O3"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# By default, a native method is used to compute the nullspace. Optionally, this can be done with sympy. The result should be the same

reaction.method = "sympy"


reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)


# ### Aspirin in water with additional reactants

reactants = ["C9H8O4", "H2O", "HO", "HO"]
products = ["C2H4O2", "C7H6O3", "H2O2"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# ### Aspirin in water with different products

reactants = ["C9H8O4", "H2O"]
products = ["CH2O", "C7H6O3"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# ### Oxydation of ethane

reactants = ["C2H6", "O2"]
products = ["CO2", "H2O"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# ### Exampe that is not easily balanced
# An example SMILES string for the reactant is O=CC(=O)c1cccc(c1)C(=O)O

reactants = ["C9O4H6", "OH"]
products = ["C2O2H3", "C7O3H2"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# ### Carbonmonoxide with carbondioxide and hydrogen

reactants = ["CO", "CO2", "H2"]
products = ["CH4", "H2O"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# ### Example of a reaction that cannot be balanced

reactants = ["FeS2", "HNO3"]
products = ["Fe2S3O12", "NO", "H2SO4"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)
print(reaction.message)


# ### Potassiumnitrate and methane

reactants = ["KNO3", "C"]
products = ["K2CO3", "CO", "N2"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# ### Pyrite and nitric acid

reactants = ["FeS2", "HNO3"]
products = ["Fe2S4O12", "N2H2"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# ### Reacting iron ligands

reactants = ["FeS2O6N5H3"]
products = ["Fe2S4O12N10H6"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# ### Example of a reaction with non-matching elements

reactants = ["KNO3", "C"]
products = ["Fe"]
reaction = ReactionEquation(reactants, products)
coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)
print(reaction.message)


# ## Test 2. Balance charged reactions

# ### Aspirin and hydroxide ion

reactants = ["C9H8O4", "OH"]
products = ["C2H3O2", "C7H6O3"]
reaction = ReactionEquation(reactants, products)
rcharges = [0, -1]
pcharges = [-1, 0]
reaction.set_charges(rcharges, pcharges)

coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# If the charges provided cannot result in a neutral reaction equation, no solution will be found

rcharges = [0, -1]
pcharges = [0, 0]
reaction.set_charges(rcharges, pcharges)

coeffs = reaction.balance()
print(reaction)
print("Coefficients: ", coeffs)


# ## Test 3. Reaction with many possible products
# Aspirin can react with water to form a very wide range of products, but most likely not all in the same reaction. Here, we supply many possible products at once, and then balance the equation towards each product in turn.

import numpy

reactant_smiles = ["O", "CC(=O)Oc1ccccc1C(=O)O"]
product_smiles = [
    "CC(=O)O",
    "O=C(O)c1ccccc1O",
    "OO",
    "O=C(O)O",
    "O=C=O",
    "CO",
    "C",
    "O=C(O)C1=CC(O)C=CC1=O",
    "O=C(O)C1=CC(O)C(O)C=C1O",
    "O=C(O)C(C=CO)=C(O)C=CO",
    "O=C(CO)Oc1ccccc1",
    "Oc1ccccc1",
    "O=C1C=CC=CC1",
    "OC1=CC(O)C=CC1",
    "C=CC=CC=C=O",
    "C1=C=CC=CC=1",
    "O=C=CO",
    "O=CO",
    "O=C=C(O)O",
    "O=C=C=O",
    "O=C1OC1=O",
    "O=COc1ccccc1",
    "O=C1CC=CC(O)C1",
    "C=CC(O)C(O)C=C=O",
    "OC1=CC=CC(O)C1=C(O)O",
    "O=C(O)C12C(=O)C1C=CC2O",
    "O=C(O)CO",
    "O=CC(=O)O",
    "O=C1C=CC(O)=C(C1)C(=O)O",
    "C=CC=CC(=O)O",
]


reactants = [from_smiles(smiles) for smiles in reactant_smiles]
products = [from_smiles(smiles) for smiles in product_smiles]

# Create the Reaction object with all the molecules
reaction = ReactionEquation(reactants, products)

print("Starting loop over products..")
nmols = len(reactants) + len(products)
nreactants = len(reactants)
for iprod, _ in enumerate(products):
    print(f"{iprod:8d} {product_smiles[iprod]:30s}: ", end="")
    min_coeffs = numpy.zeros(nmols)
    min_coeffs[nreactants + iprod] = 1
    coeffs = reaction.balance(min_coeffs)
    print(reaction)

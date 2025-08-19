import numpy
from scm.plams import from_smiles
from scm.plams import ReactionEquation


def main():
    """
    The main script
    """
    #######################################################################
    # Test 1
    # Ten different sets of reactants and products (rformulas and pformulas)
    # Can be passed as formulas or as PLAMS Molecule objects
    #######################################################################

    print("Test 1.\n Balance many small reactions")

    rformulas = [
        ["C9H8O4", from_smiles("O")],
        ["C9H8O4", "H2O", "HO", "HO"],
        ["C9H8O4", "H2O"],
        ["C2H6", "O2"],
        ["C9O4H6", "OH"],
        ["CO", "CO2", "H2"],
        ["FeS2", "HNO3"],
        ["KNO3", "C"],
        ["FeS2", "HNO3"],
        ["FeS2O6N5H3"],
    ]
    pformulas = [
        ["C2H4O2", "C7H6O3", "H2O2"],
        ["C2H4O2", "C7H6O3", "H2O2"],
        ["CH2O", "C7H6O3"],
        ["CO2", "H2O"],
        ["C2O2H3", "C7O3H2"],
        ["CH4", "H2O"],
        ["Fe2S3O12", "NO", "H2SO4"],
        ["K2CO3", "CO", "N2"],
        ["Fe2S4O12", "N2H2"],
        ["Fe2S4O12N10H6"],
    ]

    nreactions = len(rformulas)
    for i, (reactants, products) in enumerate(zip(rformulas, pformulas)):
        reaction = ReactionEquation(reactants, products)
        # reaction.method = 'sympy' # If sympy is installed, this can be used
        coeffs = reaction.balance()
        print("%8i %s" % (i, reaction))

    #######################################################################
    # Test 2
    # Now create a single set of reactants and a very large set of products.
    # Then loop over the products, and try to get a balanced reaction
    # using the reactants and the pool of other products
    #######################################################################

    print("Test 2.\n Balance a single large set of molecules towards each product in turn.")

    reactant_text = "O CC(=O)Oc1ccccc1C(=O)O"
    product_text = """CC(=O)O O=C(O)c1ccccc1O OO O=C(O)O O=C=O CO
    C O=C(O)C1=CC(O)C=CC1=O O=C(O)C1=CC(O)C(O)C=C1O O=C(O)C(C=CO)=C(O)C=CO O=C(CO)Oc1ccccc1
    Oc1ccccc1 O=C1C=CC=CC1 OC1=CC(O)C=CC1 C=CC=CC=C=O C1=C=CC=CC=1
    O=C=CO O=CO O=C=C(O)O O=C=C=O O=C1OC1=O
    O=COc1ccccc1 O=C1CC=CC(O)C1 C=CC(O)C(O)C=C=O OC1=CC=CC(O)C1=C(O)O O=C(O)C12C(=O)C1C=CC2O
    O=C(O)CO O=CC(=O)O O=C1C=CC(O)=C(C1)C(=O)O C=CC=CC(=O)O"""
    product_lines = product_text.split("\n")

    reactants = [from_smiles(smiles) for smiles in reactant_text.split()]
    psmiles = [smiles for line in product_lines for smiles in line.split()]
    products = [from_smiles(smiles) for smiles in psmiles]

    # Create the Reaction object with all the molecules
    reaction = ReactionEquation(reactants, products)

    print("Starting loop over products..")
    nmols = len(reactants) + len(products)
    nreactants = len(reactants)
    for iprod, product in enumerate(pformulas):
        print("%8i %20s: " % (iprod, psmiles[iprod]), end="")
        min_coeffs = numpy.zeros(nmols)
        min_coeffs[nreactants + iprod] = 1
        coeffs = reaction.balance(min_coeffs)
        print("%s" % (reaction))

    #######################################################################
    # Test 3
    # Test equation balancing with charged molecules
    #######################################################################

    print("Test 3.\n Balance an equation with charged molecules")

    for i in [0, 3]:
        reactants = rformulas[i]
        products = pformulas[i]
        reaction = ReactionEquation(reactants, products)
        rcharges = [1 if ir == 0 else 0 for ir, r in enumerate(reactants)]
        pcharges = [1 if ip == 0 else 0 for ip, p in enumerate(products)]
        reaction.set_charges(rcharges, pcharges)
        # reaction.method = 'sympy' # If sympy is installed, this can be used
        coeffs = reaction.balance()
        print("%8i %s" % (i, reaction))


if __name__ == "__main__":
    main()

import numpy
import pytest

from scm.plams.tools.reaction import ReactionEquation
from scm.plams.interfaces.molecule.rdkit import from_smiles


class TestEquationBalancing:

    @pytest.fixture(autouse=True)
    def skip_if_no_cbc_solver(self):
        """
        Check whether CBC solver is available for pyomo, and skip test with a warning if it is not available.
        """
        from pyomo.environ import SolverFactory

        if not SolverFactory("cbc").available(exception_flag=False):
            pytest.skip("Skipping test as cannot find CBC solver.")

    @pytest.fixture
    def reactions(self):
        """
        Reactions, defined as reactant and product molecules
        """
        return [
            {
                "rformulas": ["C9H8O4", from_smiles("O")],
                "pformulas": ["C2H4O2", "C7H6O3", "H2O2"],
                "coeffs": [1, 1, 1, 1, 0],
            },
            {
                "rformulas": ["C9H8O4", "H2O", "HO", "HO"],
                "pformulas": ["C2H4O2", "C7H6O3", "H2O2"],
                "coeffs": [1, 1, 0, 0, 1, 1, 0],
            },
            {
                "rformulas": ["C9H8O4", "H2O"],
                "pformulas": ["CH2O", "C7H6O3"],
                "coeffs": [1, 1, 2, 1],
            },
            {
                "rformulas": ["C2H6", "O2"],
                "pformulas": ["CO2", "H2O"],
                "coeffs": [2, 7, 4, 6],
            },
            {
                "rformulas": ["C9O4H6", "OH"],
                "pformulas": ["C2O2H3", "C7O3H2"],
                "coeffs": [9, 25, 23, 5],
            },
            {
                "rformulas": ["CO", "CO2", "H2"],
                "pformulas": ["CH4", "H2O"],
                "coeffs": [1, 0, 3, 1, 1],
            },
            {
                "rformulas": ["FeS2", "HNO3"],
                "pformulas": ["Fe2S3O12", "NO", "H2SO4"],
                "coeffs": None,
            },
            {
                "rformulas": ["KNO3", "C"],
                "pformulas": ["K2CO3", "CO", "N2"],
                "coeffs": [2, 4, 1, 3, 1],
            },
            {
                "rformulas": ["FeS2", "HNO3"],
                "pformulas": ["Fe2S4O12", "N2H2"],
                "coeffs": [2, 4, 1, 2],
            },
            {
                "rformulas": ["FeS2O6N5H3"],
                "pformulas": ["Fe2S4O12N10H6"],
                "coeffs": [2, 1],
            },
            {
                "rformulas": ["KNO3", "C"],
                "pformulas": ["Fe"],
                "coeffs": None,
            },
        ]

    @pytest.fixture
    def molecules(self):
        """
        A set of molecules
        """
        reactants = ["O", "CC(=O)Oc1ccccc1C(=O)O"]
        products = [
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
        coeffs = [
            [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [2, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [2, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [3, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [3, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [2, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [2, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [2, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0],
            [1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            [2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        ]
        return {"reactants": reactants, "products": products, "coeffs": coeffs}

    def test_many_small_reactions(self, reactions):
        """
        Ten different sets of reactants and products (rformulas and pformulas)
        Can be passed as formulas or as PLAMS Molecule objects
        """
        for i, reaction_mols in enumerate(reactions):
            reactants = reaction_mols["rformulas"]
            products = reaction_mols["pformulas"]
            reaction = ReactionEquation(reactants, products)
            # reaction.method = 'sympy' # If sympy is installed, this can be used
            coeffs = reaction.balance()
            print(f"{i:8d} {reaction}")
            assert (coeffs == numpy.array(reaction_mols["coeffs"])).all()

    def test_large_molecule_space(self, molecules):
        """
        Now create a single set of reactants and a very large set of products.
        Then loop over the products, and try to get a balanced reaction
        using the reactants and the pool of other products
        """
        reactants = [from_smiles(smiles) for smiles in molecules["reactants"]]
        products = [from_smiles(smiles) for smiles in molecules["products"]]

        # Create the Reaction object with all the molecules
        reaction = ReactionEquation(reactants, products)

        print("Starting loop over products..")
        nmols = len(reactants) + len(products)
        nreactants = len(reactants)
        for iprod, _ in enumerate(products):
            print(f"{iprod:8d} {molecules['products'][iprod]:20s}: ", end="")
            min_coeffs = numpy.zeros(nmols)
            min_coeffs[nreactants + iprod] = 1
            coeffs = reaction.balance(min_coeffs)
            print(f"{reaction}")
            assert (coeffs == numpy.array(molecules["coeffs"])[iprod]).all()

    def test_charged_molecules(self, reactions):
        """
        Test equation balancing with charged molecules
        """
        indices = [0, 3]
        results = [[1, 1, 1, 1, 0], None]
        for i, result in zip(indices, results):
            reactants = reactions[i]["rformulas"]
            products = reactions[i]["pformulas"]
            reaction = ReactionEquation(reactants, products)
            rcharges = [1 if ir == 0 else 0 for ir, r in enumerate(reactants)]
            pcharges = [1 if ip == 0 else 0 for ip, p in enumerate(products)]
            reaction.set_charges(rcharges, pcharges)
            # reaction.method = 'sympy' # If sympy is installed, this can be used
            coeffs = reaction.balance()
            print(f"{i:8d} {reaction}")

            assert (coeffs == numpy.array(result)).all()

            if reaction.reaction_charge != "Inconsistent system: No solution":
                assert reaction.reaction_charge == 0

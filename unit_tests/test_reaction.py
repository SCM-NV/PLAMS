import numpy as np
import pytest

from scm.plams.tools.reaction import ReactionEquation


def test_formula_inputs_do_not_require_structure_packages() -> None:
    reaction = ReactionEquation(["H2", "O2"], ["H2O"])
    assert reaction._rformulas == ["H2", "O2"]
    assert reaction._pformulas == ["H2O"]
    assert reaction._rcharges == [0, 0]
    assert reaction.rsmiles is None
    assert reaction.psmiles is None
    reaction._set_elements()
    reaction._set_matrix()
    np.testing.assert_array_equal(reaction.matrix @ np.array([2, 1, 2]), 0)


@pytest.mark.parametrize("representation", ["system", "molecule", "mixed"])
def test_structures_preserve_formula_charge_and_smiles(representation) -> None:
    pytest.importorskip("rdkit")
    if representation != "molecule":
        ChemicalSystem = pytest.importorskip("scm.base").ChemicalSystem
    from scm.plams import from_smiles

    smiles = ["[NH4+]", "[OH-]", "N", "O"]
    structures = [
        (
            ChemicalSystem.from_smiles(s, optimize_with_uff=False)
            if representation == "system" or (representation == "mixed" and i % 2 == 0)
            else from_smiles(s)
        )
        for i, s in enumerate(smiles)
    ]
    reaction = ReactionEquation(structures[:2], structures[2:])
    assert reaction._rcharges == [1, -1]
    assert reaction._pcharges == [0, 0]
    assert reaction.rsmiles is not None
    assert reaction.psmiles is not None
    reaction._set_elements()
    reaction._set_matrix()
    np.testing.assert_array_equal(reaction.matrix @ np.ones(4), 0)


def test_mixed_formulas_and_systems() -> None:
    pytest.importorskip("rdkit")
    ChemicalSystem = pytest.importorskip("scm.base").ChemicalSystem
    water = ChemicalSystem.from_smiles("O", optimize_with_uff=False)
    coordinates = np.array(water.coords, copy=True)
    reaction = ReactionEquation(["H2", "O2"], [water])
    assert reaction.rsmiles is None
    assert reaction.psmiles is not None
    assert reaction._pformulas == [water.formula()]
    np.testing.assert_array_equal(water.coords, coordinates)


def test_system_total_charge_is_not_inferred_from_smiles() -> None:
    pytest.importorskip("rdkit")
    ChemicalSystem = pytest.importorskip("scm.base").ChemicalSystem
    water = ChemicalSystem.from_smiles("O", optimize_with_uff=False)
    water.charge = 2
    reaction = ReactionEquation([water], ["H2O"])
    assert reaction._rcharges == [2]
    assert water.charge == 2


def test_balance_systems_with_solver() -> None:
    pytest.importorskip("rdkit")
    ChemicalSystem = pytest.importorskip("scm.base").ChemicalSystem
    pyomo = pytest.importorskip("pyomo.environ")
    if not pyomo.SolverFactory("cbc").available(exception_flag=False):
        pytest.skip("CBC solver is not available")
    reactants = [ChemicalSystem.from_smiles(s, optimize_with_uff=False) for s in ["C", "O=O"]]
    products = [ChemicalSystem.from_smiles(s, optimize_with_uff=False) for s in ["O=C=O", "O"]]
    reaction = ReactionEquation(reactants, products)
    np.testing.assert_array_equal(reaction.balance(), [1, 2, 1, 2])


def test_reject_unsupported_species() -> None:
    with pytest.raises(TypeError, match="Reaction species"):
        ReactionEquation([object()], ["H2O"])

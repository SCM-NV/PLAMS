import os
import pytest

import numpy as np

try:
    import dill as pickle
except ImportError:
    import pickle

from scm.plams import Atom, Molecule, MoleculeError, read_all_molecules_in_xyz_file
from scm.plams.interfaces.molecule.rdkit import from_smiles


@pytest.fixture
def BENZENE(xyz_folder):
    b = Molecule(xyz_folder / "benzene.xyz")
    b.guess_bonds()
    return b


def test_index(BENZENE):
    """Test :meth:`Molecule.index`."""
    atom = BENZENE[1]
    bond = BENZENE[1, 2]
    atom_test = Atom(coords=[0, 0, 0], symbol="H")

    assert BENZENE.index(atom) == 1
    assert BENZENE.index(bond) == (1, 2)

    try:
        BENZENE.index(None)  # None is of invalid type
    except MoleculeError:
        pass
    else:
        raise AssertionError("'BENZENE.index(None)' failed to raise a 'MoleculeError'")

    try:
        BENZENE.index(atom_test)  # atom_test is not in BENZENE
    except MoleculeError:
        pass
    else:
        raise AssertionError("'BENZENE.index(atom_test)' failed to raise a 'MoleculeError'")


def test_set_integer_bonds(BENZENE):
    """Test :meth:`Molecule.set_integer_bonds`."""
    ref1 = np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1, 1, 1], dtype=float)
    ref2 = np.array([1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1], dtype=float)

    benzene = BENZENE.copy()
    np.testing.assert_array_equal([b.order for b in benzene.bonds], ref1)

    benzene.set_integer_bonds()
    np.testing.assert_array_equal([b.order for b in benzene.bonds], ref2)


def test_round_coords(BENZENE):
    """Test :meth:`Molecule.round_coords`."""
    benzene = BENZENE.copy()
    ref1 = np.array(
        [
            [1.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
            [0.0, -1.0, 0.0],
            [2.0, -1.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 2.0, 0.0],
            [-2.0, 1.0, 0.0],
            [-2.0, -1.0, 0.0],
            [0.0, -2.0, 0.0],
        ]
    )
    ref2 = np.array(
        [
            [1.19, -0.69, 0.0],
            [1.19, 0.69, 0.0],
            [0.0, 1.38, 0.0],
            [-1.19, 0.69, 0.0],
            [-1.19, -0.69, 0.0],
            [-0.0, -1.38, 0.0],
            [2.13, -1.23, -0.0],
            [2.13, 1.23, -0.0],
            [0.0, 2.46, -0.0],
            [-2.13, 1.23, -0.0],
            [-2.13, -1.23, -0.0],
            [-0.0, -2.46, -0.0],
        ]
    )

    benzene2 = round(benzene)
    np.testing.assert_array_equal(benzene2, ref1)

    benzene.round_coords(decimals=2)
    np.testing.assert_allclose(benzene, ref2)


IMMUTABLE_TYPE = (int, float, tuple)
ATTR_EXCLUDE = frozenset({"mol", "bonds", "atoms", "atom1", "atom2", "_dummysymbol"})


def _compare_attrs(obj1, obj2, eval_eq=True):
    assert obj1 is not obj2

    for name, attr in vars(obj1).items():
        if name in ATTR_EXCLUDE:
            continue

        attr_ref = getattr(obj2, name)
        if eval_eq:
            assert attr == attr_ref
        if not isinstance(attr, IMMUTABLE_TYPE):
            assert attr is not attr_ref


def test_set_get_state(BENZENE, tmp_path):
    """Tests for :meth:`Molecule.__setstate__` and :meth:`Molecule.__getstate__`."""
    mol = BENZENE.copy()
    dill_new = tmp_path / "benzene_new.dill"

    with open(dill_new, "wb") as f:
        pickle.dump(mol, f)
    with open(dill_new, "rb") as f:
        mol_new = pickle.load(f)

    for mol in [mol_new, mol]:
        _compare_attrs(mol, BENZENE, eval_eq=False)

        for at, at_ref in zip(mol.atoms, BENZENE.atoms):
            _compare_attrs(at, at_ref)

        for bond, bond_ref in zip(mol.bonds, BENZENE.bonds):
            _compare_attrs(bond, bond_ref)


def test_read_multiple_molecules_from_xyz(xyz_folder):
    """Test for read_all_molecules_in_xyz_file"""

    filename = os.path.join(xyz_folder, "multiple_mols_in_xyz.xyz")

    mols = read_all_molecules_in_xyz_file(filename)

    assert len(mols) == 12

    for mol, n_atoms in zip(mols, [3, 5, 6, 3, 5, 6, 3, 5, 6]):
        assert len(mol.atoms) == n_atoms

    for mol, n_lattice_vec in zip(mols, [0, 0, 0, 0, 0, 0, 1, 2, 3, 1, 2, 3]):
        assert len(mol.lattice) == n_lattice_vec


def test_write_multiple_molecules_to_xyz(xyz_folder, tmp_path):
    """Test for append mode of Molecule.write and read_all_molecules_in_xyz_file"""

    new_xyz_file = tmp_path / "test_write_multiple_molecules.xyz"
    mols_ref = read_all_molecules_in_xyz_file(xyz_folder / "multiple_mols_in_xyz.xyz")

    assert len(mols_ref) > 1

    for mol in mols_ref:
        mol.write(new_xyz_file, mode="a")

    mols = read_all_molecules_in_xyz_file(new_xyz_file)
    assert len(mols_ref) == len(mols)

    for mol, mol_ref in zip(mols, mols_ref):
        for at, at_ref in zip(mol.atoms, mol_ref.atoms):
            assert at.symbol == at_ref.symbol
            np.testing.assert_allclose(at.coords, at_ref.coords, atol=1e-08)


def test_as_array_context():
    expected = np.array(
        [
            [-0.000816, 0.366378, -0.000000],
            [-0.812316, -0.183482, -0.000000],
            [0.813132, -0.182896, 0.000000],
        ]
    )
    mol = from_smiles("O")
    with mol.as_array as coord_array:
        np.testing.assert_allclose(coord_array, expected, rtol=1e-2)
        coord_array += 1
    np.testing.assert_allclose(mol.as_array(), expected + 1, rtol=1e-2)


def test_as_array_function():
    expected = np.array(
        [
            [-0.000816, 0.366378, -0.000000],
            [-0.812316, -0.183482, -0.000000],
            [0.813132, -0.182896, 0.000000],
        ]
    )
    mol = from_smiles("O")
    np.testing.assert_allclose(mol.as_array(), expected, rtol=1e-2)


def test_get_moments_of_inertia(BENZENE):
    benzene = BENZENE.copy()
    expected = np.array([86.81739308, 86.8173935, 173.63478658])
    np.testing.assert_allclose(benzene.get_moments_of_inertia(), expected, rtol=1e-2)


def test_get_gyration_radius(BENZENE):
    benzene = BENZENE.copy()
    expected = 1.6499992631225113
    np.testing.assert_allclose(benzene.get_gyration_radius(), expected, rtol=1e-2)


def test_separate():
    # previously this failed due to exceeding the maximum recursion depth
    # thus here we just make sure this doesn't throw an error now
    mol = Molecule(positions=[[float(i)] * 3 for i in range(1000)])
    for i in range(1, 1000):
        mol.add_bond(mol[i], mol[i + 1])
    mol.separate()

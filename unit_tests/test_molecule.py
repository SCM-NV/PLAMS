import os
import pytest
from abc import ABC, abstractmethod

import numpy as np

try:
    import dill as pickle
except ImportError:
    import pickle

from scm.plams.mol.atom import Atom
from scm.plams.mol.bond import Bond
from scm.plams.mol.molecule import Molecule, MoleculeError
from scm.plams.core.functions import read_all_molecules_in_xyz_file
from scm.plams.interfaces.molecule.rdkit import from_smiles


class MoleculeTestBase(ABC):
    """
    Base class which verifies that a loaded input molecule has the expected atoms, bonds, lattice and properties.
    """

    @abstractmethod
    @pytest.fixture
    def mol(self, folder): ...

    @property
    def expected_atoms(self):
        return []

    @property
    def expected_bonds(self):
        return []

    @property
    def expected_lattice(self):
        return []

    @property
    def expected_charge(self):
        return {}

    def test_init_from_file_has_atoms_bonds_lattice_and_properties_as_expected(self, mol: Molecule):
        assert [(at.symbol, *at.coords, at.properties.as_dict()) for at in mol.atoms] == self.expected_atoms
        assert [(mol.index(b), b.order, b.properties.as_dict()) for b in mol.bonds] == self.expected_bonds
        assert mol.lattice == self.expected_lattice
        assert mol.properties.charge == self.expected_charge


class TestWater(MoleculeTestBase):
    """
    Water system of just atoms
    """

    @pytest.fixture
    def mol(self, xyz_folder):
        return Molecule(xyz_folder / "water.xyz")

    @pytest.fixture
    def water(self, mol):
        return mol

    @pytest.fixture
    def hydroxide(self, xyz_folder):
        return Molecule(xyz_folder / "hydroxide.xyz")

    @property
    def expected_atoms(self):
        return [("O", 0.0, 0.0, 0.0, {}), ("H", 1.0, 0.0, 0.0, {}), ("H", 0.0, 1.0, 0.0, {})]

    def test_add_delete_atoms_and_bonds_happy(self, hydroxide, water):
        # Make peroxide from each molecule and assert the same

        # 1) Add H20 to OH-
        # 2) Delete H atom
        # 3) Add O-O bond
        peroxide1 = hydroxide.copy()
        peroxide1.add_molecule(water, copy=True, margin=1.0)
        h1 = peroxide1.atoms[-1]
        h2 = peroxide1.atoms[3]
        peroxide1.delete_atom(h1)
        peroxide1.add_bond(peroxide1[1], peroxide1[3])
        peroxide1.add_bond(peroxide1[1], peroxide1[2])
        peroxide1.add_bond(Bond(peroxide1[3], h2))
        oh1 = peroxide1[(3, 4)]
        assert h1.mol is None
        assert h2.mol == peroxide1
        assert oh1.mol == peroxide1

        # 1) Guess bonds
        # 2) Add O to H20
        # 3) Delete O-H bonds
        # 4) Add O-H and O-O bonds
        peroxide2 = water.copy()
        peroxide2.guess_bonds()
        o1 = Atom(symbol="O", coords=(2, 0, 0))
        peroxide2.add_atom(o1)
        oh1 = peroxide2[(1, 2)]
        peroxide2.delete_bond(peroxide2[1], peroxide2[2])
        peroxide2.delete_bond(peroxide2.bonds[0])
        peroxide2.add_bond(peroxide2[1], peroxide2[4])
        peroxide2.add_bond(peroxide2[1], peroxide2[2])
        peroxide2.add_bond(peroxide2[3], peroxide2[4])
        oh2 = peroxide2[(3, 4)]
        assert o1.mol == peroxide2
        assert oh1.mol is None
        assert oh2.mol == peroxide2

        # 1) Guess bonds
        # 2) Add H to H20
        # 3) Add O adjacent to O and H
        # 4) Delete H atom
        peroxide3 = water.copy()
        peroxide3.guess_bonds()
        peroxide3.add_atom(Atom(1, coords=(2, 1, 0)))
        peroxide3.add_atom(Atom(8, coords=(2, 0, 0)), adjacent=[peroxide3[1], (peroxide3[-1], 1)])
        peroxide3.delete_atom(peroxide3[3])

        # 1) Guess bonds
        # 2) Delete all bonds
        # 3) Copy molecule above
        peroxide4 = water.copy()
        peroxide4.guess_bonds()
        oh1 = peroxide4.bonds[0]
        peroxide4.delete_all_bonds()
        peroxide4 = peroxide3.copy()
        oh2 = peroxide4.bonds[-1]
        o1 = peroxide4.atoms[0]
        assert oh1.mol is None
        assert oh2.mol == peroxide4
        assert o1.mol == peroxide4

        # Assert the same
        assert peroxide1.label(3) == peroxide2.label(3) == peroxide3.label(3) == peroxide4.label(3)

    def test_add_delete_atoms_and_bonds_unhappy(self, water):
        water2 = water.copy()

        # Cannot add atom which is already part of this/another molecule
        with pytest.raises(MoleculeError):
            water.add_atom(water[1])
        with pytest.raises(MoleculeError):
            water.add_atom(water2[1])

        # Cannot delete atom which is not part of the molecule
        with pytest.raises(MoleculeError):
            water.delete_atom(water2[1])
        with pytest.raises(MoleculeError):
            water.delete_atom(Atom(1))

        # Cannot delete atom which is removed from the molecule's atom list
        h2 = water[-1]
        water.atoms.remove(h2)
        with pytest.raises(MoleculeError):
            water.delete_atom(h2)

        # Cannot add bond which has invalid arguments
        water2.guess_bonds()
        water = water2.copy()
        with pytest.raises(MoleculeError):
            water.add_bond("foo")

        # Cannot add bond which is a member of this/another molecule
        with pytest.raises(MoleculeError):
            water.add_bond(water.bonds[0])
        with pytest.raises(MoleculeError):
            water.add_bond(water2.bonds[0])

        # Cannot add bond which has atoms belonging to no/another molecule
        with pytest.raises(MoleculeError):
            water.add_bond(water2[1], water[2])
        with pytest.raises(MoleculeError):
            water.add_bond(water[2], Atom(1))

        # Cannot delete bonds which has atoms belonging to no/another molecule
        with pytest.raises(MoleculeError):
            water.delete_bond(water[1], water2[2])
        with pytest.raises(MoleculeError):
            water.delete_bond(Atom(1), water[1])
        with pytest.raises(MoleculeError):
            water.delete_bond(Bond(water[1], water[2]))

        # Cannot add bond which has invalid arguments
        with pytest.raises(MoleculeError):
            water.delete_bond("foo")

    def test_guess_bonds(self, mol):
        assert len(mol.bonds) == 0
        mol.guess_bonds()
        assert len(mol.bonds) == 2
        assert [(mol.index(b), b.order) for b in mol.bonds] == [((1, 3), 1), ((1, 2), 1)]

    def test_system_and_atomic_charge(self, mol):
        mol.guess_bonds()
        assert mol.guess_system_charge() == 0
        assert mol.guess_atomic_charges() == [0, 0, 0]
        mol.delete_atom(mol[3])
        assert mol.guess_system_charge() == -1
        assert mol.guess_atomic_charges() == [-1, 0]
        mol.add_atom(Atom(1, coords=(1, 0, 1)), adjacent=[mol[1]])
        mol.add_atom(Atom(1, coords=(1, 1, 0)), adjacent=[mol[1]])
        assert mol.guess_system_charge() == 1
        assert mol.guess_atomic_charges() == [1, 0, 0, 0]
        mol.properties.charge = 2
        assert mol.guess_atomic_charges() == [1, 1, 0, 0]
        mol.properties.charge = 3
        with pytest.raises(MoleculeError):
            assert mol.guess_atomic_charges() == [1, 1, 0, 0]


class TestNiO(MoleculeTestBase):
    """
    Periodic NiO system
    """

    @pytest.fixture
    def mol(self, xyz_folder):
        return Molecule(xyz_folder / "NiO.xyz")

    @property
    def expected_atoms(self):
        return [("Ni", 0.0, 0.0, 0.0, {}), ("O", 2.085, 2.085, 2.085, {})]

    @property
    def expected_lattice(self):
        return [[0.0, 2.085, 2.085], [2.085, 0.0, 2.085], [2.085, 2.085, 0.0]]

    def test_supercell(self, mol):
        supercell = mol.supercell(2, 3, 4)
        assert supercell.get_formula() == "Ni24O24"
        assert supercell.lattice == [(0.0, 4.17, 4.17), (6.255, 0.0, 6.255), (8.34, 8.34, 0.0)]
        with pytest.raises(MoleculeError):
            mol.supercell(2, 2)

    def test_unit_cell_volume(self, mol):
        assert mol.unit_cell_volume("bohr") == pytest.approx(122.33332352511559)

    def test_cell_lengths(self, mol):
        assert np.allclose(mol.cell_lengths("bohr"), [5.572113115975432, 5.572113115975432, 5.572113115975432])

    def test_cell_angles(self, mol):
        assert np.allclose(mol.cell_angles("radian"), [1.0471975511965976, 1.0471975511965976, 1.0471975511965976])

    class TestHydroxide(MoleculeTestBase):
        """
        Charged ion system
        """

        @pytest.fixture
        def mol(self, xyz_folder):
            mol = Molecule(xyz_folder / "hydroxide.xyz")
            mol.properties.charge = -1
            return mol

        @property
        def expected_atoms(self):
            return [("O", 1.0, 0.0, 0.0, {}), ("H", 0.0, 0.0, 0.0, {})]

        @property
        def expected_charge(self):
            return -1.0

    class TestBenzeneDimer(MoleculeTestBase):
        """
        System with atoms, bonds and properties
        """

        @pytest.fixture
        def mol(self, xyz_folder):
            mol = Molecule(xyz_folder / "benzene_dimer.xyz")
            mol.guess_bonds()
            for i, at in enumerate(mol):
                at.properties.adf.f = f"subsystem{(i // 12) + 1}"
            return mol

        @property
        def expected_atoms(self):
            return [
                ("C", -1.9000793, -0.01180491, -1.63051319, {"adf": {"f": "subsystem1"}}),
                ("C", 0.87349469, -0.01023939, -1.76821915, {"adf": {"f": "subsystem1"}}),
                ("C", -1.20564674, -1.21381033, -1.65931829, {"adf": {"f": "subsystem1"}}),
                ("C", -1.20758387, 1.19098295, -1.67103029, {"adf": {"f": "subsystem1"}}),
                ("C", 0.17915587, 1.19167022, -1.73987307, {"adf": {"f": "subsystem1"}}),
                ("C", 0.18110787, -1.21293268, -1.72799219, {"adf": {"f": "subsystem1"}}),
                ("H", -2.98312785, -0.01239662, -1.57479944, {"adf": {"f": "subsystem1"}}),
                ("H", -1.74970769, 2.12993064, -1.64644374, {"adf": {"f": "subsystem1"}}),
                ("H", 0.72018058, 2.1311289, -1.76828739, {"adf": {"f": "subsystem1"}}),
                ("H", 1.95680577, -0.00959486, -1.81797251, {"adf": {"f": "subsystem1"}}),
                ("H", 0.72362449, -2.15176759, -1.74701794, {"adf": {"f": "subsystem1"}}),
                ("H", -1.74626146, -2.15334648, -1.62565114, {"adf": {"f": "subsystem1"}}),
                ("C", -1.04276513, 0.00660438, 1.20777829, {"adf": {"f": "subsystem2"}}),
                ("C", 1.73079327, 0.00816985, 1.0700785, {"adf": {"f": "subsystem2"}}),
                ("C", -0.34843214, -1.19529353, 1.17943874, {"adf": {"f": "subsystem2"}}),
                ("C", -0.35038422, 1.20928613, 1.167558, {"adf": {"f": "subsystem2"}}),
                ("C", 1.036363, 1.21016618, 1.09888633, {"adf": {"f": "subsystem2"}}),
                ("C", 1.03830008, -1.19460875, 1.11059824, {"adf": {"f": "subsystem2"}}),
                ("H", -2.1260752, 0.00595978, 1.2575375, {"adf": {"f": "subsystem2"}}),
                ("H", -0.8929011, 2.14811978, 1.18659394, {"adf": {"f": "subsystem2"}}),
                ("H", 1.57697669, 2.14970259, 1.06522938, {"adf": {"f": "subsystem2"}}),
                ("H", 2.81384166, 0.0087615, 1.01437093, {"adf": {"f": "subsystem2"}}),
                ("H", 1.58042281, -2.13355664, 1.08602198, {"adf": {"f": "subsystem2"}}),
                ("H", -0.88945704, -2.13475088, 1.20786337, {"adf": {"f": "subsystem2"}}),
            ]

        @property
        def expected_bonds(self):
            return [
                ((13, 16), 1.5, {}),
                ((13, 15), 1.5, {}),
                ((2, 6), 1.5, {}),
                ((2, 5), 1.5, {}),
                ((15, 18), 1.5, {}),
                ((16, 17), 1.5, {}),
                ((4, 5), 1.5, {}),
                ((3, 6), 1.5, {}),
                ((14, 17), 1.5, {}),
                ((14, 18), 1.5, {}),
                ((1, 3), 1.5, {}),
                ((1, 4), 1.5, {}),
                ((13, 19), 1, {}),
                ((2, 10), 1, {}),
                ((16, 20), 1, {}),
                ((15, 24), 1, {}),
                ((6, 11), 1, {}),
                ((14, 22), 1, {}),
                ((5, 9), 1, {}),
                ((1, 7), 1, {}),
                ((18, 23), 1, {}),
                ((17, 21), 1, {}),
                ((4, 8), 1, {}),
                ((3, 12), 1, {}),
            ]

        def test_set_unset_atoms_id(self, mol):
            mol.set_atoms_id(10)
            expected = list(range(10, 34))
            assert [at.id for at in mol] == expected

            mol.delete_atom(mol[10])
            expected.remove(19)
            assert [at.id for at in mol] == expected

            mol.unset_atoms_id()
            mol.unset_atoms_id()

            with pytest.raises(AttributeError):
                mol[1].id

        def test_neighbours(self, mol):
            for at in mol:
                assert len(mol.neighbors(at)) == 3 if at.symbol == "C" else 1
            with pytest.raises(MoleculeError):
                mol.neighbors(Atom(1))

        def test_bond_matrix(self, mol):
            assert np.all(
                mol.bond_matrix()
                == np.array(
                    [
                        [0, 0, 1.5, 1.5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 1.5, 1.5, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [1.5, 0, 0, 0, 0, 1.5, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [1.5, 0, 0, 0, 1.5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1.5, 0, 1.5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1.5, 1.5, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 1.5, 0, 0, 1, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 1.5, 0, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 0, 0, 0, 0, 1.5, 0, 0, 0, 0, 0, 1],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 0, 0, 0, 1.5, 0, 0, 1, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 0, 1.5, 0, 0, 0, 0, 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 1.5, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    ]
                )
            )

        def test_separate_splits_dimer_into_two_molecules(self, mol):
            mol1, mol2 = mol.separate()
            assert mol1.label(3) == mol2.label(3)

        def test_in_ring(self, mol):
            assert mol.in_ring(mol[1])
            assert not mol.in_ring(mol[7])
            assert mol.in_ring(mol[(1, 3)])
            assert not mol.in_ring(mol[(1, 7)])

            with pytest.raises(MoleculeError):
                assert mol.in_ring(Atom(1))
            with pytest.raises(MoleculeError):
                assert mol.in_ring(Bond(Atom(1), Atom(1)))


class TestBenzene(MoleculeTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "benzene.xyz")
        mol.guess_bonds()
        return mol

    @property
    def expected_atoms(self):
        return [
            ("C", 1.1938602, -0.68927551, 0.0, {}),
            ("C", 1.1938602, 0.68927551, 0.0, {}),
            ("C", 0.0, 1.37855102, 0.0, {}),
            ("C", -1.1938602, 0.68927551, 0.0, {}),
            ("C", -1.1938602, -0.68927551, 0.0, {}),
            ("C", -0.0, -1.37855102, 0.0, {}),
            ("H", 2.13291126, -1.23143689, -0.0, {}),
            ("H", 2.13291126, 1.23143689, -0.0, {}),
            ("H", 0.0, 2.46287378, -0.0, {}),
            ("H", -2.13291126, 1.23143689, -0.0, {}),
            ("H", -2.13291126, -1.23143689, -0.0, {}),
            ("H", -0.0, -2.46287378, -0.0, {}),
        ]

    @property
    def expected_bonds(self):
        return [
            ((3, 4), 1.5, {}),
            ((5, 6), 1.5, {}),
            ((1, 6), 1.5, {}),
            ((2, 3), 1.5, {}),
            ((4, 5), 1.5, {}),
            ((1, 2), 1.5, {}),
            ((3, 9), 1, {}),
            ((6, 12), 1, {}),
            ((5, 11), 1, {}),
            ((4, 10), 1, {}),
            ((2, 8), 1, {}),
            ((1, 7), 1, {}),
        ]

    def test_index(self, mol):
        """Test :meth:`Molecule.index`."""
        atom = mol[1]
        bond = mol[1, 2]
        atom_test = Atom(coords=[0, 0, 0], symbol="H")

        assert mol.index(atom) == 1
        assert mol.index(bond) == (1, 2)

        try:
            mol.index(None)  # None is of invalid type
        except MoleculeError:
            pass
        else:
            raise AssertionError("'benzene.index(None)' failed to raise a 'MoleculeError'")

        try:
            mol.index(atom_test)  # atom_test is not in BENZENE
        except MoleculeError:
            pass
        else:
            raise AssertionError("'benzene.index(atom_test)' failed to raise a 'MoleculeError'")

    def test_set_integer_bonds(self, mol):
        """Test :meth:`Molecule.set_integer_bonds`."""
        ref1 = np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1, 1, 1], dtype=float)
        ref2 = np.array([1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1], dtype=float)

        np.testing.assert_array_equal([b.order for b in mol.bonds], ref1)

        mol.set_integer_bonds()
        np.testing.assert_array_equal([b.order for b in mol.bonds], ref2)

    def test_round_coords(self, mol):
        """Test :meth:`Molecule.round_coords`."""
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

        benzene2 = round(mol)
        np.testing.assert_array_equal(benzene2, ref1)

        mol.round_coords(decimals=2)
        np.testing.assert_allclose(mol, ref2)

    IMMUTABLE_TYPE = (int, float, tuple)
    ATTR_EXCLUDE = frozenset({"mol", "bonds", "atoms", "atom1", "atom2", "_dummysymbol"})

    def _compare_attrs(self, obj1, obj2, eval_eq=True):
        assert obj1 is not obj2

        for name, attr in vars(obj1).items():
            if name in self.ATTR_EXCLUDE:
                continue

            attr_ref = getattr(obj2, name)
            if eval_eq:
                assert attr == attr_ref
            if not isinstance(attr, self.IMMUTABLE_TYPE):
                assert attr is not attr_ref

    def test_set_get_state(self, mol, tmp_path):
        """Tests for :meth:`Molecule.__setstate__` and :meth:`Molecule.__getstate__`."""
        dill_new = tmp_path / "benzene_new.dill"
        mol2 = mol.copy()

        with open(dill_new, "wb") as f:
            pickle.dump(mol2, f)
        with open(dill_new, "rb") as f:
            mol3 = pickle.load(f)

        for m in [mol2, mol3]:
            self._compare_attrs(m, mol, eval_eq=False)

            for at, at_ref in zip(m.atoms, mol.atoms):
                self._compare_attrs(at, at_ref)

            for bond, bond_ref in zip(m.bonds, mol.bonds):
                self._compare_attrs(bond, bond_ref)

        assert mol.label(5) == mol2.label(5) == mol3.label(5)

    def test_get_moments_of_inertia(self, mol):
        expected = np.array([86.81739308, 86.8173935, 173.63478658])
        np.testing.assert_allclose(mol.get_moments_of_inertia(), expected, rtol=1e-2)

    def test_get_gyration_radius(self, mol):
        expected = 1.6499992631225113
        np.testing.assert_allclose(mol.get_gyration_radius(), expected, rtol=1e-2)


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


def test_read_multiple_molecules_from_pdb(pdb_folder):
    """
    Test for PDB reading
    """
    filenames = [os.path.join(pdb_folder, fn) for fn in os.listdir(pdb_folder)]

    mols = [(os.path.basename(fn).rstrip(".pdb"), Molecule(fn)) for fn in filenames]

    actual = {k: {"n_atoms": len(v.atoms), "n_lattice_vec": len(v.lattice)} for k, v in mols}

    expected = {
        "2kpq": {"n_atoms": 1531, "n_lattice_vec": 3},
        "1DYZ": {"n_atoms": 1024, "n_lattice_vec": 3},
        "MET": {"n_atoms": 4671, "n_lattice_vec": 3},
        "pentapeptide": {"n_atoms": 75, "n_lattice_vec": 0},
        "1BXU": {"n_atoms": 776, "n_lattice_vec": 3},
        "chymotrypsin": {"n_atoms": 69, "n_lattice_vec": 0},
    }

    assert actual == expected


def test_write_multiple_molecules_to_pdb(pdb_folder, tmp_path):
    """
    Test for PDB writing
    """
    new_pdb_file = tmp_path / "test_write_molecule.pdb"

    filenames = [os.path.join(pdb_folder, fn) for fn in os.listdir(pdb_folder)]
    mols_ref = [Molecule(fn) for fn in filenames]

    assert len(mols_ref) > 1

    mols = []
    for m in mols_ref:
        m.write(new_pdb_file)
        m.read(new_pdb_file)
        mols.append(m)

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


def test_separate():
    # previously this failed due to exceeding the maximum recursion depth
    # thus here we just make sure this doesn't throw an error now
    mol = Molecule(positions=[[float(i)] * 3 for i in range(1000)])
    for i in range(1, 1000):
        mol.add_bond(mol[i], mol[i + 1])
    mol.separate()

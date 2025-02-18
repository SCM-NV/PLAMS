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
from scm.plams.unit_tests.test_helpers import skip_if_no_ams_installation


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

        # 1) Add water to water
        # 2) Delete extra H atoms (with partial success)
        # 3) Add O-O bond and O-H bonds
        peroxide5 = water.copy()
        peroxide5.add_molecule(water, copy=True, margin=1.0)

        def atoms_to_delete():
            yield peroxide5[3]
            yield peroxide5[6]
            yield peroxide5[6]

        with pytest.raises(MoleculeError):
            peroxide5.delete_atoms(atoms_to_delete())

        peroxide5.add_bond(peroxide5[1], peroxide5[3])
        peroxide5.add_bond(peroxide5[1], peroxide5[2])
        peroxide5.add_bond(peroxide5[3], peroxide5[4])

        # Assert the same
        assert (
            peroxide1.label(3) == peroxide2.label(3) == peroxide3.label(3) == peroxide4.label(3) == peroxide5.label(3)
        )

    def test_add_delete_atoms_and_bonds_unhappy(self, water):
        water2 = water.copy()

        # Cannot add atom which is already part of this/another molecule
        with pytest.raises(MoleculeError):
            water.add_atom(water[1])
        with pytest.raises(MoleculeError):
            water.add_atom(water2[1])
        with pytest.raises(MoleculeError):
            water.add_atom(Atom(8, mol=1))

        # Cannot delete atom which is not part of the molecule
        with pytest.raises(MoleculeError):
            water.delete_atom(water2[1])
        with pytest.raises(MoleculeError):
            water.delete_atom(Atom(1))
        with pytest.raises(MoleculeError):
            water.delete_atom(Atom(8, mol=1))

        # Cannot delete atom which is removed from the molecule's atom list
        h2 = water[-1]
        water.atoms.remove(h2)
        with pytest.raises(MoleculeError):
            water.delete_atom(h2)

        # Cannot delete multiple atoms which do not form part of the molecule
        with pytest.raises(MoleculeError):
            water.delete_atoms([water2[1], water2[2]])

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
        with pytest.raises(MoleculeError):
            water.add_bond(Bond(mol=1))

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
        with pytest.raises(MoleculeError):
            water.delete_bond(Bond(mol=1))

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

    def test_add_hatoms(self, mol, hydroxide):
        skip_if_no_ams_installation()

        water = hydroxide.add_hatoms()
        assert mol.label(4) == water.label(4)

    def test_get_complete_molecules_within_threshold(self, mol):
        m0 = mol.get_complete_molecules_within_threshold([2], 0)
        m1 = mol.get_complete_molecules_within_threshold([2], 1)
        m2 = mol.get_complete_molecules_within_threshold([2], 2)

        assert m0.get_formula() == "H"
        assert m1.get_formula() == "HO"
        assert m2.get_formula() == "H2O"


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

    def test_locate_rings(self, mol):
        expected = [[0, 1, 2, 3, 4, 5]]

        assert_rings_equal(mol.locate_rings(), expected)
        assert_rings_equal(mol.locate_rings_acm(False), expected)
        assert_rings_equal(mol.locate_rings_acm(True), expected)
        assert_rings_equal(mol.locate_rings_networkx(False), expected)
        assert_rings_equal(mol.locate_rings_networkx(True), expected)

    def test_add_hatoms(self, mol):
        skip_if_no_ams_installation()

        mol.guess_bonds()
        mol2 = mol.copy()

        mol2.delete_atoms([at for at in mol2 if at.symbol == "H"])
        mol2 = mol2.add_hatoms()
        assert len(mol2.bonds) == len(mol.bonds) == 12
        assert mol.label(4) == mol2.label(4)


def assert_rings_equal(actual, expected):
    assert len(actual) == len(expected)

    actual_sorted = sorted([sorted(r) for r in actual], key=lambda x: tuple(x))
    expected_sorted = sorted([sorted(r) for r in expected], key=lambda x: tuple(x))

    for a, e in zip(actual_sorted, expected_sorted):
        assert a == e


class TestChlorophyll(MoleculeTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "chlorophyl1.xyz")
        mol.guess_bonds()
        return mol

    @property
    def expected_atoms(self):
        return [
            ("H", -10.766861, 8.365033, -6.045633, {}),
            ("H", -0.286744, -2.920043, 4.518224, {}),
            ("H", -7.415298, 0.608669, 1.334206, {}),
            ("H", -10.834144, 8.959881, -2.473094, {}),
            ("O", 0.058129, -6.185237, 4.067959, {}),
            ("H", -12.161073, 8.46221, -4.991137, {}),
            ("C", 2.283945, 0.730489, -2.240157, {}),
            ("C", 3.39891, -1.752488, 5.167145, {}),
            ("C", 3.402802, -2.296142, 3.670762, {}),
            ("C", 1.525009, 0.121913, -0.363318, {}),
            ("H", 9.806859, -1.467415, 3.441949, {}),
            ("H", 1.301458, -3.786127, 4.402401, {}),
            ("H", -4.195918, -6.015259, 0.693635, {}),
            ("C", 1.190414, -0.509741, 0.96756, {}),
            ("H", -4.48999, -4.508015, 4.190468, {}),
            ("H", 4.223896, -1.255521, 5.531767, {}),
            ("H", -6.3514, -1.079152, -1.014676, {}),
            ("C", 4.669974, 0.89791, -3.199088, {}),
            ("C", -7.345978, 2.064953, -2.23349, {}),
            ("H", 8.481015, 3.148244, -3.993377, {}),
            ("C", -8.493172, 3.219308, -2.194661, {}),
            ("C", -0.122313, 0.897138, -3.804585, {}),
            ("C", -2.529554, -5.567796, 3.768038, {}),
            ("C", -5.219917, -2.164816, 1.750015, {}),
            ("C", -7.608128, 0.958592, -1.200238, {}),
            ("H", -2.471951, -5.318337, 1.112253, {}),
            ("C", 9.74045, -0.431376, 0.845466, {}),
            ("N", 2.730995, 0.327746, -0.850865, {}),
            ("C", 3.304909, 0.888652, -3.310867, {}),
            ("H", -4.960662, -1.965396, 0.708109, {}),
            ("N", 6.225564, -0.949205, 0.620143, {}),
            ("H", -1.773843, 3.673051, 2.302649, {}),
            ("H", -8.230617, 1.366851, -0.488136, {}),
            ("H", -8.545374, 6.391313, -4.427422, {}),
            ("C", 0.852813, 0.563976, -2.533048, {}),
            ("C", -9.980858, 2.586231, -1.597982, {}),
            ("H", -7.562378, 4.628748, -3.404755, {}),
            ("H", 4.790817, 3.516315, -5.25694, {}),
            ("H", -10.547323, 3.488711, -1.615968, {}),
            ("H", -1.012162, 1.131489, -3.277309, {}),
            ("C", 6.624225, 0.412022, -2.312189, {}),
            ("C", 6.823042, 1.51135, -3.17925, {}),
            ("H", -0.84211, -1.555423, 1.075281, {}),
            ("H", 6.449132, 3.786982, -4.357694, {}),
            ("C", 7.925689, 2.628908, -3.160471, {}),
            ("H", -7.173342, -2.193418, 2.027235, {}),
            ("H", 7.318388, 3.486856, -2.658193, {}),
            ("H", -10.026741, 5.80902, -5.39413, {}),
            ("C", -6.606974, -1.318924, 1.922339, {}),
            ("C", -0.987339, 0.55172, 1.998737, {}),
            ("H", 2.667181, 1.201734, -4.183633, {}),
            ("C", 5.728988, 3.060378, -4.755714, {}),
            ("O", -1.754705, -0.942991, -1.117896, {}),
            ("C", 10.088243, 0.796957, 0.187929, {}),
            ("C", 6.425326, 2.228959, -5.969432, {}),
            ("H", -9.324222, 8.264518, -4.120834, {}),
            ("O", -1.30968, 1.748042, 1.470894, {}),
            ("H", -9.060533, -1.333128, -0.557864, {}),
            ("H", -10.936466, 9.680486, -4.838663, {}),
            ("C", -0.326369, -5.094058, 3.680961, {}),
            ("C", -10.204073, 7.629442, -3.985463, {}),
            ("C", -0.420381, -0.490072, 0.843535, {}),
            ("H", -10.745989, 5.132236, -2.829537, {}),
            ("H", 2.425485, -1.451884, 5.548839, {}),
            ("H", 0.055595, -3.018999, 2.078651, {}),
            ("H", -7.053729, 1.492742, -3.188656, {}),
            ("C", -10.953883, 7.882267, -2.707784, {}),
            ("H", 1.40338, -4.110744, 2.154364, {}),
            ("C", 7.468435, -0.46759, 0.041097, {}),
            ("C", -3.693547, -4.709673, 3.354146, {}),
            ("H", -0.246236, 0.089078, -4.436668, {}),
            ("C", 2.191908, -0.901031, 1.888952, {}),
            ("H", 7.102653, 1.449605, -5.912669, {}),
            ("N", 3.521576, -0.845287, 1.555414, {}),
            ("H", -10.241334, 1.679415, -2.289627, {}),
            ("H", -6.440342, 2.521538, -2.043528, {}),
            ("C", 8.865725, -2.075635, 3.346505, {}),
            ("H", 3.801483, -3.408068, 3.424338, {}),
            ("C", 7.601947, 0.100293, -1.188034, {}),
            ("H", -6.261545, -0.588856, 2.747129, {}),
            ("H", -5.40613, 0.403446, -0.428726, {}),
            ("C", -0.782589, -0.280752, -0.70216, {}),
            ("C", -5.312017, -3.652771, 2.037982, {}),
            ("C", -9.692578, 5.248955, -3.03941, {}),
            ("H", -3.512435, -4.186557, 0.049757, {}),
            ("H", -12.120777, 7.926221, -2.962663, {}),
            ("H", -8.718395, 3.645703, -4.260665, {}),
            ("H", -10.036967, 2.36433, -0.443992, {}),
            ("H", -6.027046, -3.971709, 1.350556, {}),
            ("O", -1.264415, 0.075426, 3.035442, {}),
            ("C", -4.069872, -4.488776, 2.120056, {}),
            ("H", -10.587002, 7.350663, -1.866547, {}),
            ("H", 10.637029, -1.017122, 1.16185, {}),
            ("H", -4.387733, -1.657275, 2.163674, {}),
            ("C", -6.529114, 0.008244, -0.604541, {}),
            ("C", 5.682598, 1.966172, -3.738103, {}),
            ("N", 5.307446, 0.046101, -2.292931, {}),
            ("H", -0.047618, 3.135137, 2.817219, {}),
            ("H", -5.987226, -3.643865, 2.889564, {}),
            ("H", -8.449786, 0.458982, -1.514008, {}),
            ("Mg", 4.595064, -0.269632, -0.224787, {}),
            ("H", -8.919284, -2.160391, 1.132701, {}),
            ("H", 8.109073, -2.027523, 4.044068, {}),
            ("C", -9.676433, 6.257077, -4.412732, {}),
            ("H", 5.71833, 1.493872, -6.445224, {}),
            ("C", 0.328322, -3.605438, 3.920669, {}),
            ("C", 1.910528, -1.915478, 3.032284, {}),
            ("O", -1.524799, -4.887849, 3.03888, {}),
            ("H", 11.016857, 1.181377, 0.072677, {}),
            ("H", 9.012089, -3.175713, 3.072764, {}),
            ("H", -2.306797, -5.585637, 4.89145, {}),
            ("H", 6.703209, 3.192731, -6.574592, {}),
            ("C", 5.545183, -1.567259, 2.815592, {}),
            ("C", 1.042307, -3.160158, 2.584764, {}),
            ("H", -8.057422, 3.842984, -1.414606, {}),
            ("C", -7.387919, -0.478765, 0.787443, {}),
            ("C", -3.519154, -5.004612, 0.852819, {}),
            ("C", -11.098484, 8.578374, -5.027057, {}),
            ("H", 5.961186, -2.242689, 3.741442, {}),
            ("H", -9.680979, -0.523255, 0.704715, {}),
            ("H", -1.170132, 2.472334, 3.552079, {}),
            ("C", 4.237676, -1.510186, 2.627258, {}),
            ("H", -9.223372, 5.916082, -2.261976, {}),
            ("H", 8.791811, 0.121368, -1.229643, {}),
            ("C", 0.402983, 0.230521, -1.220625, {}),
            ("C", -1.125977, 2.691073, 2.393361, {}),
            ("H", 0.195614, 1.707686, -4.385004, {}),
            ("H", 3.244867, -2.537885, 5.905112, {}),
            ("C", 7.982473, -1.395619, 2.093753, {}),
            ("C", -8.857802, -1.104676, 0.471232, {}),
            ("C", -8.693485, 4.079447, -3.344344, {}),
            ("H", 1.273952, -1.566403, 3.702863, {}),
            ("H", -2.621283, -6.743847, 3.361242, {}),
            ("C", 8.38621, -0.764389, 1.010433, {}),
            ("C", 6.536428, -1.334041, 1.872737, {}),
            ("H", 9.340226, 1.260172, -0.496771, {}),
            ("H", 8.904538, 2.401677, -2.825987, {}),
        ]

    @property
    def expected_bonds(self):
        return [
            ((27, 134), 1, {}),
            ((129, 134), 2, {}),
            ((107, 132), 1, {}),
            ((113, 122), 1, {}),
            ((35, 125), 1, {}),
            ((57, 126), 1, {}),
            ((53, 82), 1.5, {}),
            ((21, 131), 1, {}),
            ((50, 57), 1, {}),
            ((42, 96), 2, {}),
            ((54, 109), 1, {}),
            ((87, 131), 1, {}),
            ((69, 79), 2, {}),
            ((91, 117), 1, {}),
            ((69, 134), 1, {}),
            ((25, 100), 1, {}),
            ((60, 108), 1, {}),
            ((8, 16), 1, {}),
            ((25, 33), 1, {}),
            ((77, 103), 1, {}),
            ((19, 76), 1, {}),
            ((22, 71), 1, {}),
            ((7, 29), 1, {}),
            ((55, 73), 1, {}),
            ((10, 28), 1.5, {}),
            ((120, 130), 1, {}),
            ((31, 135), 1, {}),
            ((52, 96), 1, {}),
            ((113, 135), 2, {}),
            ((83, 91), 1, {}),
            ((83, 89), 1, {}),
            ((82, 125), 1.5, {}),
            ((61, 67), 1, {}),
            ((23, 70), 1, {}),
            ((22, 127), 1, {}),
            ((46, 49), 1, {}),
            ((10, 14), 1, {}),
            ((24, 83), 1, {}),
            ((41, 79), 1, {}),
            ((24, 94), 1, {}),
            ((45, 137), 1, {}),
            ((23, 108), 1, {}),
            ((41, 42), 1, {}),
            ((22, 40), 1, {}),
            ((67, 92), 1, {}),
            ((14, 72), 2, {}),
            ((10, 125), 1.5, {}),
            ((61, 104), 1, {}),
            ((36, 39), 1, {}),
            ((19, 25), 1, {}),
            ((6, 118), 1, {}),
            ((41, 97), 2, {}),
            ((58, 130), 1, {}),
            ((18, 97), 1, {}),
            ((9, 122), 1, {}),
            ((72, 74), 1, {}),
            ((72, 107), 1, {}),
            ((63, 84), 1, {}),
            ((25, 95), 1, {}),
            ((83, 99), 1, {}),
            ((8, 64), 1, {}),
            ((8, 128), 1, {}),
            ((18, 96), 1, {}),
            ((84, 131), 1, {}),
            ((21, 115), 1, {}),
            ((42, 45), 1, {}),
            ((24, 30), 1, {}),
            ((1, 118), 1, {}),
            ((56, 61), 1, {}),
            ((106, 114), 1, {}),
            ((129, 135), 1, {}),
            ((44, 52), 1, {}),
            ((2, 106), 1, {}),
            ((107, 114), 1, {}),
            ((12, 106), 1, {}),
            ((68, 114), 1, {}),
            ((8, 9), 1, {}),
            ((4, 67), 1, {}),
            ((62, 82), 1, {}),
            ((54, 136), 1, {}),
            ((27, 93), 1, {}),
            ((65, 114), 1, {}),
            ((31, 69), 1, {}),
            ((49, 116), 1, {}),
            ((14, 62), 1, {}),
            ((26, 117), 1, {}),
            ((11, 77), 1, {}),
            ((29, 51), 1, {}),
            ((55, 105), 1, {}),
            ((84, 123), 1, {}),
            ((20, 45), 1, {}),
            ((59, 118), 1, {}),
            ((52, 55), 1, {}),
            ((19, 21), 1, {}),
            ((116, 130), 1, {}),
            ((48, 104), 1, {}),
            ((24, 49), 1, {}),
            ((22, 35), 1, {}),
            ((34, 104), 1, {}),
            ((74, 122), 2, {}),
            ((60, 106), 1, {}),
            ((77, 110), 1, {}),
            ((23, 111), 1, {}),
            ((85, 117), 1, {}),
            ((50, 62), 1, {}),
            ((19, 66), 1, {}),
            ((49, 80), 1, {}),
            ((38, 52), 1, {}),
            ((9, 107), 1, {}),
            ((61, 118), 1, {}),
            ((45, 47), 1, {}),
            ((77, 129), 1, {}),
            ((7, 28), 1, {}),
            ((43, 62), 1, {}),
            ((36, 75), 1, {}),
            ((55, 112), 1, {}),
            ((15, 70), 1, {}),
            ((17, 95), 1, {}),
            ((36, 88), 1, {}),
            ((32, 126), 1, {}),
            ((121, 126), 1, {}),
            ((84, 104), 1, {}),
            ((95, 116), 1, {}),
            ((79, 124), 1, {}),
            ((67, 86), 1, {}),
            ((21, 36), 1, {}),
            ((81, 95), 1, {}),
            ((9, 78), 1, {}),
            ((3, 116), 1, {}),
            ((113, 119), 1, {}),
            ((13, 117), 1, {}),
            ((98, 126), 1, {}),
            ((102, 130), 1, {}),
            ((23, 133), 1, {}),
            ((37, 131), 1, {}),
            ((50, 90), 2.0, {}),
            ((70, 91), 2.0, {}),
            ((5, 60), 2.0, {}),
            ((18, 29), 2.0, {}),
            ((27, 54), 2.0, {}),
            ((7, 35), 2.0, {}),
            ((101, 74), 1, {}),
            ((101, 28), 1, {}),
            ((101, 31), 1, {}),
            ((101, 97), 1, {}),
        ]

    def test_locate_rings(self, mol):
        expected = [
            [6, 27, 9, 124, 34],
            [6, 27, 100, 96, 17, 28],
            [8, 106, 71, 73, 121],
            [9, 13, 61, 81, 124],
            [9, 13, 71, 73, 100, 27],
            [17, 95, 41, 40, 96],
            [30, 68, 133, 128, 134],
            [30, 68, 78, 40, 96, 100],
            [30, 100, 73, 121, 112, 134],
        ]

        assert_rings_equal(mol.locate_rings(), expected)
        assert_rings_equal(mol.locate_rings_acm(False), expected)
        assert_rings_equal(mol.locate_rings_acm(True), expected)
        assert_rings_equal(mol.locate_rings_networkx(True), expected)


class TestFragments(MoleculeTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "C7H8N2_fragments.xyz")
        return mol

    @property
    def expected_atoms(self):
        return [
            ("H", 37.97589166, 6.50362513, 6.10947932, {}),
            ("H", 22.94566879, 5.21944993, 6.11208918, {}),
            ("H", 24.24811668, 3.58937684, 3.1101599, {}),
            ("H", 26.46029574, 4.22221813, 1.9970365, {}),
            ("H", 28.11383033, 5.49145028, 3.3383532, {}),
            ("H", 33.53518443, 7.87366456, 5.79895021, {}),
            ("C", 36.81091549, 5.79838054, 2.99659108, {}),
            ("N", 35.79722124, 6.77372035, 4.90114805, {}),
            ("H", 27.40840017, 6.85618188, 6.76571824, {}),
            ("C", 35.70493885, 6.05732335, 2.17566355, {}),
            ("C", 34.68676663, 6.93124889, 4.1845686, {}),
            ("C", 34.61783813, 6.60177631, 2.80991515, {}),
            ("N", 25.60906993, 5.44707201, 5.57733345, {}),
            ("N", 32.40678698, 7.45865484, 4.11253188, {}),
            ("N", 28.93739263, 6.57981707, 5.45707573, {}),
            ("H", 37.67050745, 5.3014948, 2.59000904, {}),
            ("H", 35.78463338, 5.85496894, 1.08423825, {}),
            ("H", 33.70046648, 6.8213344, 2.30167641, {}),
            ("C", 37.95605417, 5.90776431, 5.20779435, {}),
            ("C", 23.41069028, 4.45252478, 5.49358156, {}),
            ("C", 24.74504443, 4.71605898, 4.8543161, {}),
            ("C", 25.04631067, 4.17412898, 3.6262066, {}),
            ("C", 26.2268722, 4.52245672, 2.98318138, {}),
            ("C", 26.77975995, 5.66497554, 5.02272776, {}),
            ("C", 27.10838378, 5.26861749, 3.73156573, {}),
            ("C", 33.47999143, 7.48497558, 4.78073603, {}),
            ("C", 36.82761129, 6.16730707, 4.30869146, {}),
            ("C", 27.73379633, 6.4261697, 5.80579523, {}),
            ("H", 29.31659198, 7.50946073, 5.34720982, {}),
            ("H", 32.44075585, 7.41892399, 3.10388547, {}),
        ]

    def test_add_hatoms(self, mol):
        skip_if_no_ams_installation()

        # Given two fragments with no bond information
        # When add H atoms
        mol2 = mol.add_hatoms()

        # Then adds hydrogens according to guessed bonds
        assert mol2.get_formula() == "C14H15N4"

        # But given bonding information
        mol3 = mol.copy()
        mol3.guess_bonds()
        mol3[19, 27].order = 1

        # When add H atoms
        mol4 = mol3.add_hatoms()

        # Then adds hydrogens which respect modified bond orders
        assert mol4.get_formula() == "C14H16N4"


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


def test_read_multiple_molecules_from_coskf(coskf_folder):
    """
    Test for COSKF reading
    """
    filenames = [os.path.join(coskf_folder, fn) for fn in os.listdir(coskf_folder)]

    mols = [(os.path.basename(fn).rstrip(".coskf"), Molecule(fn)) for fn in filenames]

    actual = {k: {"n_atoms": len(v.atoms), "n_bonds": len(v.bonds), "charge": v.properties.charge} for k, v in mols}

    expected = {
        "Water": {"n_atoms": 3, "n_bonds": 2, "charge": 0.0},
        "IL_cation_ethylpyridinium": {"n_atoms": 18, "n_bonds": 18, "charge": 1.0},
        "Benzene": {"n_atoms": 12, "n_bonds": 12, "charge": 0.0},
        "IL_anion_acetate": {"n_atoms": 7, "n_bonds": 6, "charge": -1.0},
    }

    assert actual == expected


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

import pytest

from scm.plams.mol.bond import Bond
from scm.plams.mol.atom import Atom
from scm.plams.core.errors import MoleculeError


class TestBond:
    """
    Test suite for the Bond class
    """

    @pytest.fixture
    def hydrogen_dimer_bond(self):
        """
        Bond in hydrogen dimer
        """
        return Bond(atom1=Atom(1), atom2=Atom(1, coords=(0.74, 0, 0)))

    @pytest.fixture
    def carbon_bonds(self):
        """
        Selection of carbon bond types:
            - carbon single bond
            - carbon aromatic bond
            - carbon double bond
            - carbon triple bond
        """
        return [
            Bond(atom1=Atom(6), atom2=Atom(6, coords=(1.54, 0, 0))),
            Bond(atom1=Atom(6), atom2=Atom(6, coords=(1.52, 0, 0)), order=1.5),
            Bond(atom1=Atom(6), atom2=Atom(6, coords=(1.48, 0, 0)), order=2),
            Bond(atom1=Atom(6), atom2=Atom(6, coords=(1.38, 0, 0)), order=3),
        ]

    def test_is_aromatic_as_expected(self, carbon_bonds):
        assert [c.is_aromatic() for c in carbon_bonds] == [False, True, False, False]

    def test_length_as_expected(self, hydrogen_dimer_bond):
        assert hydrogen_dimer_bond.length() == 0.74

    def test_as_vector_as_expected(self, hydrogen_dimer_bond):
        assert hydrogen_dimer_bond.as_vector(start=hydrogen_dimer_bond.atom1) == (0.74, 0, 0)
        assert hydrogen_dimer_bond.as_vector(start=hydrogen_dimer_bond.atom2) == (-0.74, 0, 0)
        with pytest.raises(MoleculeError):
            assert hydrogen_dimer_bond.as_vector(start=Atom(2))

    def test_other_end_as_expected(self, hydrogen_dimer_bond):
        assert hydrogen_dimer_bond.other_end(hydrogen_dimer_bond.atom1) == hydrogen_dimer_bond.atom2
        assert hydrogen_dimer_bond.other_end(hydrogen_dimer_bond.atom2) == hydrogen_dimer_bond.atom1
        with pytest.raises(MoleculeError):
            assert hydrogen_dimer_bond.other_end(Atom(2))

    def test_resize_without_molecule_translates_atom_as_expected(self, hydrogen_dimer_bond):
        hydrogen_dimer_bond.resize(hydrogen_dimer_bond.atom1, 0.75)
        hydrogen_dimer_bond.resize(hydrogen_dimer_bond.atom2, 0.76)
        assert hydrogen_dimer_bond.length() == 0.76
        assert hydrogen_dimer_bond.atom1.x == pytest.approx(-0.01)
        assert hydrogen_dimer_bond.atom2.x == pytest.approx(0.75)

    @pytest.mark.parametrize("suffix,expected", [(None, False), ("", False), (" ", False), ("foo", True)])
    def test_has_cell_shifts_as_expected(self, suffix, expected):
        bond = Bond(suffix=suffix)
        assert bond.has_cell_shifts() == expected

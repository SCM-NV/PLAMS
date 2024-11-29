import pytest

from scm.plams.core.errors import PTError
from scm.plams.mol.atom import Atom


class TestAtom:
    """
    Test suite for the Atom class
    """

    @pytest.fixture
    def helium_atom(self):
        """
        Free Helium atom
        """
        return Atom(2, "He", (10.0, 10.0, 10.0))

    def test_default_atom_has_dummy_values(self):
        atom = Atom()

        assert atom.symbol == "Xx"
        assert atom.atnum == 0
        assert atom.coords == (0, 0, 0)
        assert atom.mass == 0
        assert atom.radius == 0
        assert atom.connectors == 0
        assert not atom.is_metallic
        assert not atom.is_electronegative

    def test_unhappy_atom_raises_error(self):
        with pytest.raises(PTError):
            Atom(symbol=42)
        with pytest.raises(TypeError):
            Atom(atnum="foo")
        with pytest.raises(TypeError):
            Atom(symbol="C", coords="foobar")

    def test_get_set_xyz_coords_as_expected(self):
        atom = Atom()
        atom.x = 10
        atom.y = 11
        atom.z = 12

        assert atom.x == 10
        assert atom.y == 11
        assert atom.z == 12

    @pytest.mark.parametrize("symbol,expected", [("H", "H"), ("Zn", "Zn"), ("zn", "Zn")])
    def test_get_set_symbol_as_expected(self, symbol, expected):
        atom = Atom()
        atom.symbol = symbol

        assert atom.symbol == expected

    def test_get_properties_as_expected(self, helium_atom):
        assert helium_atom.mass == 4.0026
        assert helium_atom.radius == 0.99
        assert helium_atom.connectors == 0
        assert not helium_atom.is_metallic
        assert not helium_atom.is_electronegative

    @pytest.mark.parametrize("unit,expected", [["angstrom", (11, 12, 13)], ["nm", (20, 30, 40)]])
    def test_translate_as_expected(self, helium_atom, unit, expected):
        # When translate by vector
        helium_atom.translate(vector=(1.0, 2.0, 3.0), unit=unit)

        # Then vector applied to coordinates
        assert helium_atom.coords == expected

    @pytest.mark.parametrize("unit,expected", [["angstrom", (1, 2, 3)], ["nm", (10, 20, 30)]])
    def test_move_to_as_expected(self, helium_atom, unit, expected):
        # When move to point
        helium_atom.move_to(point=(1.0, 2.0, 3.0), unit=unit)

        # Then coordinates set to point
        assert helium_atom.coords == expected

    @pytest.mark.parametrize(
        "unit,result_unit,expected",
        [
            ["angstrom", "angstrom", 7.0710678118654755],
            ["nm", "angstrom", 225.61028345356956],
            ["nm", "nm", 22.56102834535696],
            ["angstrom", "nm", 0.7071067811865476],
        ],
    )
    def test_distance_to_as_expected(self, helium_atom, unit, result_unit, expected):
        # When calculate distance to point
        distance = helium_atom.distance_to(point=(13.0, 14.0, 15.0), unit=unit, result_unit=result_unit)

        # Then distance correct
        assert distance == pytest.approx(expected)

    @pytest.mark.parametrize(
        "unit,result_unit,expected",
        [
            ["angstrom", "angstrom", (3.0, 4.0, 5.0)],
            ["nm", "angstrom", (120, 130, 140)],
            ["nm", "nm", (12, 13, 14)],
            ["angstrom", "nm", (0.3, 0.4, 0.5)],
        ],
    )
    def test_vector_to_as_expected(self, helium_atom, unit, result_unit, expected):
        # When calculate vector to point
        distance = helium_atom.vector_to(point=(13.0, 14.0, 15.0), unit=unit, result_unit=result_unit)

        # Then vector correct
        assert distance == pytest.approx(expected)

    @pytest.mark.parametrize(
        "unit1,unit2,result_unit,expected",
        [
            ["angstrom", "angstrom", "radian", 0.908184792647238],
            ["nm", "angstrom", "radian", 1.9627880696609752],
            ["nm", "nm", "radian", 1.7851442256767078],
            ["angstrom", "nm", "radian", 1.8669850379272077],
            ["angstrom", "nm", "degree", 106.97036308730092],
        ],
    )
    def test_angle_as_expected(self, helium_atom, unit1, unit2, result_unit, expected):
        # When calculate angle to points
        angle = helium_atom.angle(
            point1=(-3, 4, 5), point2=(6, -7, 8), point1unit=unit1, point2unit=unit2, result_unit=result_unit
        )

        # Then angle correct
        assert angle == pytest.approx(expected)

    @pytest.mark.parametrize(
        "matrix,expected",
        [
            [[1, 2, 3, 4, 5, 6, 7, 8, 9], (60.0, 150.0, 240.0)],
            [[[1, 2, 3], [4, 5, 6], [7, 8, 9]], (60.0, 150.0, 240.0)],
        ],
    )
    def test_rotate_as_expected(self, helium_atom, matrix, expected):
        # When rotate with given matrix
        helium_atom.rotate(matrix)

        # Then coords correct
        assert helium_atom.coords == pytest.approx(expected)

    def test_transformation_round_trip_leaves_crds_unchanged(self, helium_atom):
        # When perform operations and their inverses
        helium_atom.move_to((100, 100, 100))
        helium_atom.translate((50, 50, 50))
        helium_atom.rotate([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        helium_atom.rotate([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        helium_atom.translate((-50, -50, -50))

        # Then atom position unchanged
        assert helium_atom.coords == pytest.approx((100, 100, 100))

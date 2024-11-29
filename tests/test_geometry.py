import pytest
import numpy as np

from scm.plams.tools.geometry import (
    rotation_matrix,
    axis_rotation_matrix,
    distance_array,
    angle,
    dihedral,
    cell_shape,
    cell_lengths,
    cell_angles,
    cellvectors_from_shape,
)

HALF_PI = np.pi / 2


class TestGeometryTools:
    """
    Test suite for geometry tools
    """

    @pytest.mark.parametrize(
        "vec1,vec2,expected",
        [
            [
                [1, 0, 0],
                [0, 1, 0],
                [
                    [0, -1, 0],
                    [1, 0, 0],
                    [0, 0, 1],
                ],
            ],
            [
                [2, 0, 0],
                [0, 0, -1],
                [
                    [0, 0, 1],
                    [0, 1, 0],
                    [-1, 0, 0],
                ],
            ],
            [
                [2, 0, 0],
                [6, 0, 0],
                [
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                ],
            ],
            [
                [2, 0, 0],
                [-6, 0, 0],
                [
                    [-1, 0, 0],
                    [0, -1, 0],
                    [0, 0, -1],
                ],
            ],
            [
                [1.1, -2.2, 3.3],
                [-4.4, 5.5, 6.6],
                [
                    [0.72464129, 0.49324661, -0.48125127],
                    [0.31050314, 0.38974556, 0.86699838],
                    [0.61520956, -0.77769286, 0.12927111],
                ],
            ],
        ],
        ids=["xy_unit", "xz", "xx", "xx_negative", "arbitrary"],
    )
    def test_rotation_matrix_as_expected(self, vec1, vec2, expected):
        actual = rotation_matrix(vec1, vec2)
        assert np.allclose(actual, expected)

    @pytest.mark.parametrize(
        "vec,angle,expected",
        [
            [
                [1, 0, 0],
                90,
                [
                    [1, 0, 0],
                    [0, 0, -1],
                    [0, 1, 0],
                ],
            ],
            [
                [2, 2, 0],
                90,
                [
                    [0.5, 0.5, 0.707106781],
                    [0.5, 0.5, -0.707106781],
                    [-0.707106781, 0.707106781, 0],
                ],
            ],
            [
                [0, 0, 1],
                -180,
                [
                    [-1, 0, 0],
                    [0, -1, 0],
                    [0, 0, 1],
                ],
            ],
        ],
        ids=["x_unit_90", "xy_90", "z_-180"],
    )
    def test_axis_rotation_matrix_as_expected(self, vec, angle, expected):
        actual = axis_rotation_matrix(vec, angle, "degree")
        assert np.allclose(actual, expected)

    @pytest.mark.parametrize(
        "ar,expected",
        [
            [
                [[0, 0, 0]],
                [[0, 1, 5, 13]],
            ],
            [
                [[6, 8, 0], [-9, -8, -24]],
                [
                    [10, 9.43398113, 5, 19.20937271],
                    [26.85144316, 27.20294102, 29.39387691, 14],
                ],
            ],
        ],
        ids=["origin", "arbitrary"],
    )
    def test_distance_array_as_expected(self, ar, expected):
        actual = distance_array(ar, [[0, 0, 0], [1, 0, 0], [3, 4, 0], [-3, -4, -12]])
        assert np.allclose(actual, expected)

    @pytest.mark.parametrize(
        "vec1,vec2,expected",
        [
            [[1, 0, 0], [0, 1, 0], 90],
            [[2, 0, 0], [0, 0, -1], 90],
            [[2, 0, 0], [6, 0, 0], 0],
            [[2, 0, 0], [-6, 0, 0], 180],
            [[1.1, -2.2, 3.3], [-4.4, 5.5, 6.6], 83.00232996984745],
        ],
        ids=["xy_unit", "xz", "xx", "xx_negative", "arbitrary"],
    )
    def test_angle_as_expected(self, vec1, vec2, expected):
        actual = angle(vec1, vec2, "degree")
        assert actual == pytest.approx(expected)

    @pytest.mark.parametrize(
        "p1,p2,p3,p4,expected",
        [
            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 1, 0], 180],
            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1], 90],
            [[0, 0, 0], [1, 2, 3], [-42, 64, 12], [100, -41, 43], 193.81249587213918],
        ],
        ids=["colinear", "orthogonal", "arbitrary"],
    )
    def test_dihedral_as_expected(self, p1, p2, p3, p4, expected):
        actual = dihedral(p1, p2, p3, p4, "degree")
        assert actual == pytest.approx(expected)

    @pytest.fixture(
        params=[
            {
                "lattice": [
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                ],
                "lengths": [1, 1, 1],
                "angles": [HALF_PI, HALF_PI, HALF_PI],
            },
            {
                "lattice": [
                    [-4, 0, 0],
                    [0, -4, 0],
                    [0, 0, -3],
                ],
                "lengths": [4, 4, 3],
                "angles": [HALF_PI, HALF_PI, HALF_PI],
            },
            {
                "lattice": [
                    [1, 0, 0],
                    [3, 4, 0],
                    [0, 0, 2],
                ],
                "lengths": [1, 5, 2],
                "angles": [HALF_PI, HALF_PI, 0.9272952180016123],
            },
        ],
        ids=["simple_cubic", "tetragonal", "monoclinic"],
    )
    def cell(self, request):
        return request.param

    def test_cell_shape_as_expected(self, cell):
        actual = cell_shape(cell["lattice"])
        assert np.allclose(actual, cell["lengths"] + cell["angles"])

    def test_cell_lengths_as_expected(self, cell):
        actual = cell_lengths(cell["lattice"])
        assert np.allclose(actual, cell["lengths"])

    def test_cell_angles_as_expected(self, cell):
        actual = cell_angles(cell["lattice"], "radian")
        assert np.allclose(actual, cell["angles"])

    def test_cellvectors_from_shape_with_lengths_and_angles_recreates_cell_as_expected(self, cell):
        actual = cellvectors_from_shape(cell["lengths"] + cell["angles"])
        assert np.allclose(actual, np.abs(cell["lattice"]))

    def test_cellvectors_from_shape_with_lengths_only_recreates_cell_assuming_orthogonal_vectors(self, cell):
        actual = cellvectors_from_shape(cell["lengths"])

        if np.allclose(cell["angles"], HALF_PI):
            assert np.allclose(actual, np.abs(cell["lattice"]))
        else:
            x, y, z = actual
            assert np.linalg.norm(x) == pytest.approx(cell["lengths"][0])
            assert np.linalg.norm(y) == pytest.approx(cell["lengths"][1])
            assert np.linalg.norm(z) == pytest.approx(cell["lengths"][2])
            assert np.dot(x, y) == pytest.approx(0)
            assert np.dot(y, z) == pytest.approx(0)
            assert np.dot(z, x) == pytest.approx(0)

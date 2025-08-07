import PIL
import pytest
from unittest.mock import patch
import numpy as np

from scm.plams.tools.view import view, ViewConfig, _get_view_plane
from scm.plams.mol.molecule import Molecule

try:
    from scm.libbase import UnifiedChemicalSystem as ChemicalSystem, UnifiedLattice as Lattice

    _has_scm_chemsys = True
except ImportError:
    _has_scm_chemsys = False


class TestView:

    @pytest.fixture
    def water(self, xyz_folder):
        water = Molecule(xyz_folder / "water.xyz")
        water.guess_bonds()
        return water

    @pytest.mark.parametrize(
        "view_config, expected",
        [
            (
                None,
                "-transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 -1.000000 -fixedatomsize -hideregions -showunitcell thickness 0.05 -batch",
            ),
            (
                ViewConfig(width=100, height=100, normal=(1.0, 0.0, 0.0)),
                "-transparent -scmgeometry 100x100 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 1.000000 0.000000 0.000000 -fixedatomsize -hideregions -showunitcell thickness 0.05 -batch",
            ),
            (
                ViewConfig(
                    fixed_atom_size=False,
                    show_atom_labels=True,
                    atom_label_color="#FFFFFF",
                    atom_label_size=2,
                    show_regions=True,
                ),
                "-transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 -1.000000 -atomlabel AtomType -labelcolor #FFFFFF -labelsize 2 -showunitcell thickness 0.05 -batch",
            ),
            (
                ViewConfig(show_unit_cell_edges=True, unit_cell_edge_thickness=0.2),
                "-transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 -1.000000 -fixedatomsize -hideregions -showunitcell thickness 0.2 -batch",
            ),
            (
                ViewConfig(show_unit_cell_faces=True, show_lattice_vectors=True),
                "-transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 1 -viewplane 0.000000 0.000000 -1.000000 -fixedatomsize -hideregions -showunitcell faces -batch",
            ),
            (
                ViewConfig(dpi=600),
                "-transparent -scmgeometry 800x400 -dpi 600 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 -1.000000 -fixedatomsize -hideregions -showunitcell thickness 0.05 -batch",
            ),
        ],
        ids=[
            "default",
            "width_height_normal",
            "atom_labels_show_regions",
            "unit_cell_edges",
            "unit_cell_faces_lattice_vectors",
            "pic_dpi",
        ],
    )
    def test_view_command_passed_to_amsview(self, view_config, expected, water):
        with patch("scm.plams.tools.view.run_with_timeout") as mock_run_with_timeout:
            # This call will fail to generate the image due to the mock
            # but we are just testing that the command to AMSview is generated properly
            try:
                view(
                    water,
                    config=view_config or ViewConfig(),
                )
            except PIL.UnidentifiedImageError:
                pass

            mock_run_with_timeout.assert_called()

            called_args, _ = mock_run_with_timeout.call_args
            command = str.join(" ", called_args[0][4:])
            assert command == expected

    @pytest.mark.parametrize(
        "lattice_as_basis, normal, lattice, expected",
        [
            (False, (1, 0, 0), [[5, 0, 0], [0, 5, 0], [0, 0, 5]], "1.000000 0.000000 0.000000"),
            (False, (0, 2, 0), [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 1.000000 0.000000"),
            (False, (0, 0, 1.5), [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.000000 0.000000 1.000000"),
            (True, (1, 0, 0), [[5, 0, 0], [0, 5, 0], [0, 0, 5]], "1.000000 0.000000 0.000000"),
            (True, (0, 2, 0), [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.707107 -0.707107 0.000000"),
            (True, (0, 0, 1.5), [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.000000 -0.178885 0.983870"),
            (False, (2, 0, 0), [[1, 0, 0]], "1.000000 0.000000 0.000000"),
            (False, (3, 0, 0), [[5, 5, 0], [5, -5, 0]], "1.000000 0.000000 0.000000"),
            (False, (0, 0, 4), [[5, 5, 0], [5, -5, 0]], "0.000000 0.000000 1.000000"),
            (True, (2, 0, 0), [[1, 0, 0]], "1.000000 0.000000 0.000000"),
            (True, (3, 0, 0), [[5, 5, 0], [5, -5, 0]], "0.707107 0.707107 0.000000"),
            (True, (0, 0, 4), [[5, 5, 0], [5, -5, 0]], "0.000000 0.000000 1.000000"),
        ],
    )
    def test_get_view_plane_with_normal(self, lattice_as_basis, normal, lattice, expected):
        mols = []
        box = Molecule()
        box.lattice = lattice
        mols.append(box)
        if _has_scm_chemsys:
            box_cs = ChemicalSystem()
            box_cs.lattice = Lattice(np.array(lattice))
            mols.append(box_cs)

        for mol in mols:
            actual = _get_view_plane(mol, ViewConfig(lattice_as_basis=lattice_as_basis, normal=normal))
            assert actual == expected

    @pytest.mark.parametrize(
        "lattice_as_basis, direction, lattice, expected",
        [
            (False, "along_x", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "-1.000000 0.000000 0.000000"),
            (False, "along_y", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 -1.000000 0.000000"),
            (False, "along_z", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 0.000000 -1.000000"),
            (True, "along_x", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "-0.707107 -0.707107 0.000000"),
            (True, "along_y", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "-0.707107 0.707107 0.000000"),
            (True, "along_z", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 0.000000 -1.000000"),
            (False, "tilt_x", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-0.990148 0.099015 0.099015"),
            (False, "small_tilt_y", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.049875 -0.997509 0.049875"),
            (False, "large_tilt_z", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.192450 0.192450 -0.962250"),
            (True, "tilt_x", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-0.975055 -0.121556 0.185721"),
            (True, "small_tilt_y", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.563589 -0.824928 0.043150"),
            (True, "large_tilt_z", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.082627 0.359045 -0.929656"),
            (False, "corner_x", [[1, 0, 0]], "-0.577350 0.577350 0.577350"),
            (False, "corner_y", [[5, 5, 0], [5, -5, 0]], "0.577350 -0.577350 0.577350"),
            (False, "corner_z", [[5, 5, 0], [5, -5, 0]], "0.577350 0.577350 -0.577350"),
            (True, "corner_x", [[1, 0, 0]], "-0.577350 0.577350 0.577350"),
            (True, "corner_y", [[5, 5, 0], [5, -5, 0]], "0.000000 0.816497 0.577350"),
            (True, "corner_z", [[5, 5, 0], [5, -5, 0]], "0.816497 0.000000 -0.577350"),
        ],
    )
    def test_get_view_plane_with_happy_direction(self, lattice_as_basis, direction, lattice, expected):
        mols = []
        box = Molecule()
        box.lattice = lattice
        mols.append(box)
        if _has_scm_chemsys:
            box_cs = ChemicalSystem()
            box_cs.lattice = Lattice(np.array(lattice))
            mols.append(box_cs)

        for mol in mols:
            actual = _get_view_plane(mol, ViewConfig(lattice_as_basis=lattice_as_basis, direction=direction))
            assert actual == expected

    @pytest.mark.parametrize(
        "direction", ["p", "xx", "_x", "along_--x", "along__x", "view_x", "view_+x", "small_along_x", "med_tilt_x"]
    )
    def test_get_view_plane_with_unhappy_direction(self, direction, water):
        with pytest.raises(ValueError):
            _get_view_plane(water, ViewConfig(direction=direction))

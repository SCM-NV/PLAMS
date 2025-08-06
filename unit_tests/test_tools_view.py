import PIL
import pytest
from unittest.mock import patch

from scm.plams.tools.view import view, ViewConfig, RepresentationConfig, PeriodicConfig, PictureConfig
from scm.plams.mol.molecule import Molecule


class TestView:

    @pytest.fixture
    def water(self, xyz_folder):
        water = Molecule(xyz_folder / "water.xyz")
        water.guess_bonds()
        return water

    @pytest.mark.parametrize(
        "view_config, repr_config, periodic_config, picture_config, expected",
        [
            (
                None,
                None,
                None,
                None,
                "-transparent -scmgeometry 800x400 -dpi 300 -padding 0.0 -showlatticevectors 0 -viewplane 0.0 0.0 1.0 -fixedatomsize -hideregions -showunitcell thickness 0.05 -batch",
            ),
            (
                ViewConfig(width=100, height=100, normal=(1.0, 0.0, 0.0)),
                None,
                None,
                None,
                "-transparent -scmgeometry 100x100 -dpi 300 -padding 0.0 -showlatticevectors 0 -viewplane 1.0 0.0 0.0 -fixedatomsize -hideregions -showunitcell thickness 0.05 -batch",
            ),
            (
                None,
                RepresentationConfig(
                    fixed_atom_size=False,
                    show_atom_labels=True,
                    atom_label_color="#FFFFFF",
                    atom_label_size=2,
                    show_regions=True,
                ),
                None,
                None,
                "-transparent -scmgeometry 800x400 -dpi 300 -padding 0.0 -showlatticevectors 0 -viewplane 0.0 0.0 1.0 -atomlabel AtomType -labelcolor #FFFFFF -labelsize 2 -showunitcell thickness 0.05 -batch",
            ),
            (
                None,
                None,
                PeriodicConfig(show_edges=True, edge_thickness=0.2),
                None,
                "-transparent -scmgeometry 800x400 -dpi 300 -padding 0.0 -showlatticevectors 0 -viewplane 0.0 0.0 1.0 -fixedatomsize -hideregions -showunitcell thickness 0.2 -batch",
            ),
            (
                None,
                None,
                PeriodicConfig(show_faces=True, show_lattice_vectors=True),
                None,
                "-transparent -scmgeometry 800x400 -dpi 300 -padding 0.0 -showlatticevectors 1 -viewplane 0.0 0.0 1.0 -fixedatomsize -hideregions -showunitcell faces -batch",
            ),
            (
                    None,
                    None,
                    None,
                    PictureConfig(dpi=600),
                    "-transparent -scmgeometry 800x400 -dpi 600 -padding 0.0 -showlatticevectors 0 -viewplane 0.0 0.0 1.0 -fixedatomsize -hideregions -showunitcell thickness 0.05 -batch",
            ),
        ],
        ids=[
            "default",
            "view_width_height_normal",
            "repr_atom_labels_show_regions",
            "per_edges",
            "per_faces_lattice_vectors",
            "pic_dpi"
        ],
    )
    def test_view_command_passed_to_amsview(
        self, view_config, repr_config, periodic_config, picture_config, expected, water
    ):
        with patch("scm.plams.tools.view.run_with_timeout") as mock_run_with_timeout:
            # This call will fail to generate the image due to the mock
            # but we are just testing that the command to AMSview is generated properly
            try:
                view(
                    water,
                    view=view_config or ViewConfig(),
                    representation=repr_config or RepresentationConfig(),
                    periodic=periodic_config or PeriodicConfig(),
                    picture=picture_config,
                )
            except PIL.UnidentifiedImageError:
                pass

            mock_run_with_timeout.assert_called()

            called_args, _ = mock_run_with_timeout.call_args
            command = str.join(" ", called_args[0][4:])
            assert command == expected

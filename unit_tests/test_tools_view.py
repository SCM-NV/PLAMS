import os
import pytest
from unittest.mock import patch, MagicMock
import numpy as np

from scm.plams.interfaces.adfsuite.errors import AMSExecutionError
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.tools.view import view, ViewConfig, _AmsViewBackend, _AmsViewXvfbBackend, _XvfbManager, _AsePlotBackend
from scm.plams.mol.molecule import Molecule
from test_helpers import skip_if_windows

try:
    from scm.libbase import ChemicalSystem, Lattice

    _has_scm_chemsys = True
except ImportError:
    _has_scm_chemsys = False


@pytest.fixture
def water(xyz_folder):
    water = from_smiles("O")
    water.guess_bonds()
    return water


class TestView:

    def test_backends_cache(self, water):
        # Given backend cache with no successful backends
        import scm.plams.tools.view as viewer

        viewer._view_backends_cache = {
            "amsview": (_AmsViewBackend(), False, RuntimeError("something went wrong")),
            "amsview_xvfb": (_AmsViewXvfbBackend(), False, RuntimeError("something also went wrong")),
            "ase_plot": (_AsePlotBackend(), False, RuntimeError("something else went wrong")),
        }

        # When view
        # Then raises error
        with pytest.raises(RuntimeError):
            view(water, backend="auto")
        with pytest.raises(RuntimeError):
            view(water, backend="amsview")

        # Given backend cache with successful backend
        viewer._view_backends_cache["ase_plot"] = (_AsePlotBackend(), True, None)

        # When view
        # Then succeeds
        view(water, backend="auto")

        viewer._view_backends_cache = None


class TestAmsViewBackend:

    backend = _AmsViewBackend

    def test_check_available(self, monkeypatch):
        monkeypatch.setenv("AMSBIN", "test/ams")
        with patch("subprocess.run") as mock_run:
            response = MagicMock()
            response.stderr = "something went wrong"
            mock_run.return_value = response
            with pytest.raises(AMSExecutionError):
                self.backend.check_available()

    @pytest.mark.parametrize(
        "view_config, expected",
        [
            (
                ViewConfig(),
                "foo.in -transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 1.000000 -fixedatomsize -hideregions -showunitcell 0.05 -save bar.png -batch",
            ),
            (
                ViewConfig(width=100, height=100, normal=(1.0, 0.0, 0.0)),
                "foo.in -transparent -scmgeometry 100x100 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 1.000000 0.000000 0.000000 -fixedatomsize -hideregions -showunitcell 0.05 -save bar.png -batch",
            ),
            (
                ViewConfig(
                    fixed_atom_size=False,
                    show_atom_labels=True,
                    atom_label_color="#FFFFFF",
                    atom_label_size=2,
                    show_regions=True,
                ),
                "foo.in -transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 1.000000 -atomlabel Element -labelcolor #FFFFFF -labelsize 2 -showunitcell 0.05 -save bar.png -batch",
            ),
            (
                ViewConfig(show_unit_cell_edges=True, unit_cell_edge_thickness=0.2),
                "foo.in -transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 1.000000 -fixedatomsize -hideregions -showunitcell 0.2 -save bar.png -batch",
            ),
            (
                ViewConfig(show_unit_cell_faces=True, show_lattice_vectors=True),
                "foo.in -transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 1 -viewplane 0.000000 0.000000 1.000000 -fixedatomsize -hideregions -showunitcell faces -save bar.png -batch",
            ),
            (
                ViewConfig(dpi=600),
                "foo.in -transparent -scmgeometry 800x400 -dpi 600 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 1.000000 -fixedatomsize -hideregions -showunitcell 0.05 -save bar.png -batch",
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
    def test_get_command(self, view_config, expected, water):
        command = self.backend.get_command(water, view_config, input_path="foo.in", img_path="bar.png")

        command = str.join(" ", command[1:])
        assert command == expected

    @pytest.mark.parametrize(
        "normal_basis, normal, lattice, expected",
        [
            ("xyz", (1, 0, 0), [[5, 0, 0], [0, 5, 0], [0, 0, 5]], "1.000000 0.000000 0.000000"),
            ("xyz", (0, 2, 0), [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 1.000000 0.000000"),
            ("xyz", (0, 0, 1.5), [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.000000 0.000000 1.000000"),
            ("abc", (1, 0, 0), [[5, 0, 0], [0, 5, 0], [0, 0, 5]], "1.000000 0.000000 0.000000"),
            ("abc", (0, 2, 0), [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.707107 -0.707107 0.000000"),
            ("abc", (0, 0, 1.5), [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.000000 -0.178885 0.983870"),
            ("xyz", (2, 0, 0), [[1, 0, 0]], "1.000000 0.000000 0.000000"),
            ("xyz", (3, 0, 0), [[5, 5, 0], [5, -5, 0]], "1.000000 0.000000 0.000000"),
            ("xyz", (0, 0, 4), [[5, 5, 0], [5, -5, 0]], "0.000000 0.000000 1.000000"),
            ("abc", (2, 0, 0), [[1, 0, 0]], "1.000000 0.000000 0.000000"),
            ("abc", (3, 0, 0), [[5, 5, 0], [5, -5, 0]], "0.707107 0.707107 0.000000"),
            ("abc", (0, 0, 4), [[5, 5, 0], [5, -5, 0]], "0.000000 0.000000 1.000000"),
        ],
    )
    def test_get_view_plane_with_normal(self, normal_basis, normal, lattice, expected):
        mols = []
        box = Molecule()
        box.lattice = lattice
        mols.append(box)
        if _has_scm_chemsys:
            box_cs = ChemicalSystem()
            box_cs.lattice = Lattice(np.array(lattice))
            mols.append(box_cs)

        for mol in mols:
            actual = self.backend.get_view_plane(mol, ViewConfig(normal_basis=normal_basis, normal=normal))
            assert [float(v) for v in actual.split()] == [float(v) for v in expected.split()]

    @pytest.mark.parametrize(
        "direction, lattice, expected",
        [
            ("along_x", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "1.000000 0.000000 0.000000"),
            ("along_y", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 1.000000 0.000000"),
            ("along_z", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 0.000000 1.000000"),
            ("along_a", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.707107 0.707107 0.000000"),
            ("along_b", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.707107 -0.707107 0.000000"),
            ("along_c", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 0.000000 1.000000"),
            ("along_pca1", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "1.000000 0.000000 0.000000"),
            ("along_pca2", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 0.707107 0.707107"),
            ("along_pca3", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 -0.707107 0.707107"),
            ("tilt_x", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-0.990148 0.099015 0.099015"),
            ("small_tilt_y", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.049875 -0.997509 0.049875"),
            ("large_tilt_z", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.192450 0.192450 -0.962250"),
            ("tilt_a", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-0.975055 -0.121556 0.185721"),
            ("small_tilt_b", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.563589 -0.824928 0.043150"),
            ("large_tilt_c", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.082627 0.359045 -0.929656"),
            ("tilt_pca1", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-0.990148 0.000000 0.140028"),
            ("small_tilt_pca2", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.049875 -0.740613 -0.670078"),
            ("large_tilt_pca3", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.192450 0.816497 -0.544331"),
            ("corner_x", [[1, 0, 0]], "-0.577350 0.577350 0.577350"),
            ("corner_y", [[5, 5, 0], [5, -5, 0]], "0.577350 -0.577350 0.577350"),
            ("corner_z", [[5, 5, 0], [5, -5, 0]], "0.577350 0.577350 -0.577350"),
            ("corner_a", [[1, 0, 0]], "-0.577350 0.577350 0.577350"),
            ("corner_b", [[5, 5, 0], [5, -5, 0]], "0.000000 0.816497 0.577350"),
            ("corner_c", [[5, 5, 0], [5, -5, 0]], "0.816497 0.000000 -0.577350"),
            ("corner_pca1", [[1, 0, 0]], "-0.577350 0.000000 0.816497"),
            ("corner_pca2", [[5, 5, 0], [5, -5, 0]], "0.577350 -0.816497 0.000000"),
            ("corner_pca3", [[5, 5, 0], [5, -5, 0]], "0.577350 0.816497 0.000000"),
        ],
    )
    def test_get_view_plane_with_happy_direction(self, direction, lattice, expected):
        mols = []
        # Set up (unphysical) atoms to test PCA
        #   C
        # C-C-C-C
        #   C
        box = Molecule(
            positions=[[-100, 0, 10], [-98, 0, 10], [-96, 0, 10], [-94, 0, 10], [-98, 1, 11], [-98, -1, 9]],
            numbers=[6, 6, 6, 6, 6, 6],
        )
        box.lattice = lattice
        mols.append(box)
        if _has_scm_chemsys:
            box_cs = ChemicalSystem(AMSJob(box).get_input())
            mols.append(box_cs)

        for mol in mols:
            actual = self.backend.get_view_plane(mol, ViewConfig(direction=direction))
            assert [float(v) for v in actual.split()] == [float(v) for v in expected.split()]

    @pytest.mark.parametrize(
        "direction", ["p", "xx", "_x", "along_--x", "along__x", "view_x", "view_+x", "small_along_x", "med_tilt_x"]
    )
    def test_get_view_plane_with_unhappy_direction(self, direction, water):
        with pytest.raises(ValueError):
            self.backend.get_view_plane(water, ViewConfig(direction=direction))


class TestAmsViewXvfbBackend(TestAmsViewBackend):

    backend = _AmsViewXvfbBackend

    def test_check_available(self, monkeypatch):
        with patch("shutil.which") as mock_which, patch("subprocess.run") as mock_run:
            mock_which.return_value = "foo/Xvfb"

            def raise_filenotfounderror(*args, **kwargs):
                raise FileNotFoundError("Cannot find xvfb")

            mock_run.side_effect = raise_filenotfounderror
            with pytest.raises(RuntimeError):
                self.backend.check_available()


class TestXvfbManager:

    @pytest.fixture(autouse=True)
    def reset_singleton(self):
        skip_if_windows()
        _XvfbManager._instance = None

    def test_get_command(self, monkeypatch):
        # Given parameters
        manager = _XvfbManager(size=(800, 400), color_depth=16)
        manager._write_file_descriptor = 42

        # When get command
        # Then values passed to xvfb command
        assert " ".join(manager._command) == "Xvfb -nolisten tcp -screen 0 800x400x16 -displayfd 42"

    def test_alive(self):
        # Given manager with no process
        manager = _XvfbManager()
        manager._proc = None
        # When check alive
        # Then fails
        assert not manager.alive

        # Given manager with ongoing process
        # When check alive
        # Then passes
        manager._proc = MagicMock()
        manager._proc.poll.return_value = None
        assert manager.alive

        # Given manager with finished process
        # When check alive
        # Then fails and stdout/stderr captured
        manager._proc.poll.return_value = 0
        manager._proc.communicate.return_value = ("foo", "bar")
        assert not manager.alive
        assert manager._stdout == "foo"
        assert manager._stderr == "bar"

    def test_await_display_number(self):
        # Given manager with failed process
        manager = _XvfbManager(startup_timeout=0.2)
        manager._read_file_descriptor, manager._write_file_descriptor = os.pipe()
        manager._proc = MagicMock()
        manager._proc.poll.return_value = 1
        manager._proc.communicate.return_value = ("foo", "bar")
        assert manager.display is None

        # When wait for display number
        # Then errors
        with pytest.raises(RuntimeError):
            manager._await_display_number()

        # Given manager with long ongoing process
        manager._proc.poll.return_value = None

        # When wait for display number
        # Then errors
        with pytest.raises(TimeoutError):
            manager._await_display_number()

        # Given manager with process writing to file descriptor
        os.write(manager._write_file_descriptor, b"99\n")

        # When wait for display number
        # Then set
        manager._await_display_number()
        assert manager.display_number == 99
        assert manager.display == ":99"

        os.close(manager._read_file_descriptor)
        os.close(manager._write_file_descriptor)

    def test_check_xvfb(self):
        # Given manager
        manager = _XvfbManager()

        # When Xvfb not on path
        # Then check fails
        with patch("shutil.which") as mock_which:
            mock_which.return_value = None
            with pytest.raises(RuntimeError):
                manager.check_xvfb()

        with patch("shutil.which") as mock_which:
            mock_which.return_value = "foo/Xvfb"

            # When Xvfb cannot be run
            # Then check fails
            with patch("subprocess.run") as mock_run:

                def raise_filenotfounderror(*args, **kwargs):
                    raise FileNotFoundError("Cannot find xvfb")

                mock_run.side_effect = raise_filenotfounderror
                with pytest.raises(RuntimeError):
                    assert manager.check_xvfb()

            with patch("subprocess.run") as mock_run:
                mock_ret = MagicMock()
                mock_run.return_value = mock_ret
                # When does not have displayfd
                # Then check fails
                with pytest.raises(RuntimeError):
                    manager.check_xvfb()

                mock_ret.stderr = """\
                -nocursor              disable the cursor
                -core                  generate core dump on fatal error
                -displayfd fd          file descriptor to write display number to when ready to connect
                """
                # Otherwise passes
                manager.check_xvfb()

    def test_start(self):
        # Given manager failing xvfb command
        manager = _XvfbManager(startup_timeout=0.2, startup_retries=3)

        # Given manager with xvfb command which initially fails then succeeds
        with patch("subprocess.Popen") as mock_popen, patch("shutil.which") as mock_which, patch(
            "subprocess.run"
        ) as mock_run:
            mock_which.return_value = "foo/Xvfb"
            mock_run_ret = MagicMock()
            mock_run_ret.stderr = "-displayfd fd"
            mock_run.return_value = mock_run_ret
            calls = [0]

            def fails_once_then_passes(*args, **kwargs):
                calls[0] += 1
                mock_ret = MagicMock()
                if calls[0] == 1:
                    mock_ret.communicate.return_value = ("foo", "bar")
                    return mock_ret
                elif calls[0] == 2:
                    mock_ret.poll.return_value = None
                    os.write(manager._write_file_descriptor, b"99\n")
                    return mock_ret
                return

            mock_popen.side_effect = fails_once_then_passes

            # When start
            manager.start()

            # Then display set
            assert manager.display_number == 99
            assert manager.display == ":99"

            manager.start()
            manager.start()

    def test_session(self):
        # Given manager
        manager = _XvfbManager()

        # When not started
        with pytest.raises(RuntimeError):

            # Then session raises error
            with manager.session():
                pass

        # When started
        with patch("subprocess.Popen") as mock_popen, patch("shutil.which") as mock_which, patch(
            "subprocess.run"
        ) as mock_run:
            mock_which.return_value = "foo/Xvfb"
            mock_run_ret = MagicMock()
            mock_run_ret.stderr = "-displayfd fd"
            mock_run.return_value = mock_run_ret

            def setup(*args, **kwargs):
                mock_ret = MagicMock()
                mock_ret.poll.return_value = None
                os.write(manager._write_file_descriptor, b"99\n")
                return mock_ret

            mock_popen.side_effect = setup

            manager.start()

        # Then session sets display env var
        user_env = os.environ.copy()
        user_env["TEST"] = "foo"
        with manager.session(env=user_env) as env:
            assert env["TEST"] == "foo"
            assert env["DISPLAY"] == ":99"
        with manager.session() as env:
            assert env["DISPLAY"] == ":99"

    def test_stop(self):
        # Given manager
        manager = _XvfbManager()

        # When stopped without being started
        # Then call is a no-op
        manager.stop()
        assert manager.display is None

        # When started
        with patch("subprocess.Popen") as mock_popen, patch("shutil.which") as mock_which, patch(
            "subprocess.run"
        ) as mock_run:
            mock_which.return_value = "foo/Xvfb"
            mock_run_ret = MagicMock()
            mock_run_ret.stderr = "-displayfd fd"
            mock_run.return_value = mock_run_ret

            def setup(*args, **kwargs):
                mock_ret = MagicMock()
                mock_ret.poll.return_value = None
                os.write(manager._write_file_descriptor, b"99\n")
                return mock_ret

            mock_popen.side_effect = setup

            manager.start()
            assert manager.display_number == 99

        # Then stop clears display variable
        manager.stop()
        assert manager.display_number is None

    def test_singleton(self):
        # Given two managers
        manager1 = _XvfbManager()
        manager2 = _XvfbManager()

        # When start one manager
        with patch("subprocess.Popen") as mock_popen, patch("shutil.which") as mock_which, patch(
            "subprocess.run"
        ) as mock_run:
            mock_which.return_value = "foo/Xvfb"
            mock_run_ret = MagicMock()
            mock_run_ret.stderr = "-displayfd fd"
            mock_run.return_value = mock_run_ret

            def setup(*args, **kwargs):
                mock_ret = MagicMock()
                mock_ret.poll.return_value = None
                os.write(manager1._write_file_descriptor, b"99\n")
                return mock_ret

            mock_popen.side_effect = setup

            manager1.start()

        # Then second manager is also on the same display
        assert manager2.display_number == manager1.display_number == 99


class TestAsePlotBackend:

    backend = _AsePlotBackend

    @pytest.mark.parametrize(
        "direction, lattice, expected",
        [
            ("along_x", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "-90y"),
            ("along_y", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "90x"),
            ("along_z", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], ""),
            ("along_a", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "90.00x,-45.00y,-90.00z"),
            ("along_b", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "-90.00x,-45.00y,90.00z"),
            ("along_c", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], ""),
            ("tilt_x", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "45.00x,81.95y,90.00z"),
            ("small_tilt_y", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-87.14x,-2.86y,90.00z"),
            ("large_tilt_z", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "168.69x,-11.10y,-7.22z"),
            ("tilt_a", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-33.21x,77.18y,-90.00z"),
            ("small_tilt_b", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-87.01x,-34.30y,90.00z"),
            ("large_tilt_c", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "158.88x,-4.74y,-90.00z"),
            ("corner_x", [[1, 0, 0]], "45.00x,35.26y,90.00z"),
            ("corner_y", [[5, 5, 0], [5, -5, 0]], "-45.00x,-35.26y,90.00z"),
            ("corner_z", [[5, 5, 0], [5, -5, 0]], "135.00x,-35.26y,-90.00z"),
            ("corner_a", [[1, 0, 0]], "45.00x,35.26y,90.00z"),
            ("corner_b", [[5, 5, 0], [5, -5, 0]], "54.74x,0.00y,0.00z"),
            ("corner_c", [[5, 5, 0], [5, -5, 0]], "180.00x,-54.74y,90.00z"),
        ],
    )
    def test_get_view_rotation_with_happy_direction(self, direction, lattice, expected):
        mols = []
        box = Molecule()
        box.lattice = lattice
        mols.append(box)
        if _has_scm_chemsys:
            box_cs = ChemicalSystem()
            box_cs.lattice = Lattice(np.array(lattice))
            mols.append(box_cs)

        for mol in mols:
            actual = self.backend.get_view_rotation(mol, ViewConfig(direction=direction))

            # assert just on x and y - the scipy version seems to lead to variations in the final z rotation
            def decompose(s: str):
                if s == "":
                    return 0, 0, 0
                vals = {"x": 0, "y": 0, "z": 0}
                for part in s.split(","):
                    ax = part[-1]
                    val = float(part[:-1])
                    vals[ax] = val
                return tuple(vals.values())

            actual_x, actual_y, _ = decompose(actual)
            exp_x, exp_y, _ = decompose(expected)
            assert actual_x == exp_x
            assert actual_y == exp_y

    def test_view(self):
        # Given a series of molecules, chemical systems and options, check that the view method generates an image (not the contents)
        single_water = from_smiles("O")

        single_water_in_1d_box = single_water.copy()
        single_water_in_1d_box.lattice = [[10.0, 0.0, 0.0]]
        single_water_in_1d_box.atoms[0].properties.region = {"O"}

        single_water_in_2d_box = single_water_in_1d_box.copy()
        single_water_in_2d_box.lattice = [[10.0, 0.0, 0.0], [0.0, 12.0, 0.0]]
        single_water_in_2d_box.atoms[1].properties.region = {"H"}

        single_water_in_3d_box = single_water_in_2d_box.copy()
        single_water_in_3d_box.lattice = [[10.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 15.0]]
        single_water_in_3d_box.atoms[2].properties.region = {"H"}

        single_water_in_3d_non_orthorhombic_box = single_water_in_3d_box.copy()
        single_water_in_3d_non_orthorhombic_box.lattice = [[10.0, 2.0, -1.0], [-5.0, 8.0, 0.0], [0.0, -2.0, 11.0]]
        single_water_in_3d_non_orthorhombic_box.atoms[0].properties.region = {"water"}
        single_water_in_3d_non_orthorhombic_box.atoms[1].properties.region = {"water"}
        single_water_in_3d_non_orthorhombic_box.atoms[2].properties.region = {"water"}

        # These are not necessarily pretty, just functional to test the options!
        configs = []
        configs.append(ViewConfig())
        configs.append(
            ViewConfig(
                padding=1,
                direction="tilt_z",
                dpi=300,
                fixed_atom_size=True,
                show_atom_labels=True,
                atom_label_type="Element",
                atom_label_color="#FFFFFF",
                atom_label_size=2,
                show_regions=True,
                show_unit_cell_edges=True,
                unit_cell_edge_thickness=2,
                show_lattice_vectors=True,
            )
        )
        configs.append(
            ViewConfig(
                normal=(0.3, 0.3, 0.3),
                normal_basis="abc",
                fixed_atom_size=False,
                show_atom_labels=True,
                atom_label_type="Name",
                show_unit_cell_faces=True,
            )
        )

        molecules = [
            single_water,
            single_water_in_1d_box,
            single_water_in_2d_box,
            single_water_in_3d_box,
            single_water_in_3d_non_orthorhombic_box,
        ]
        if _has_scm_chemsys:
            chem_systems = []
            for m in molecules:
                chem_systems.append(ChemicalSystem(AMSJob(molecule=m).get_input()))
            molecules += chem_systems

        for m in molecules:
            for c in configs:
                view(m, c, backend="ase_plot")

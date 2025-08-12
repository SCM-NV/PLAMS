import os
import pytest
from unittest.mock import patch, MagicMock
import numpy as np

from scm.plams.interfaces.adfsuite.errors import AMSExecutionError
from scm.plams.tools.view import view, ViewConfig, _AmsViewBackend, _AmsViewXvfbBackend, _XvfbManager, _AsePlotViewBackend
from scm.plams.mol.molecule import Molecule
from test_helpers import skip_if_no_ams_installation

try:
    from scm.libbase import UnifiedChemicalSystem as ChemicalSystem, UnifiedLattice as Lattice

    _has_scm_chemsys = True
except ImportError:
    _has_scm_chemsys = False


@pytest.fixture
def water(xyz_folder):
    water = Molecule(xyz_folder / "water.xyz")
    water.guess_bonds()
    return water


class TestView:

    def test_backends_cache(self, water):
        # Given backend cache with no successful backends
        view._backends = {
            "amsview": (_AmsViewBackend(), False, RuntimeError("something went wrong")),
            "amsview_xvfb": (_AmsViewXvfbBackend(), False, RuntimeError("something also went wrong")),
            "ase_plot": (_AsePlotViewBackend(), False, RuntimeError("something else went wrong"))
        }

        # When view
        # Then raises error
        with pytest.raises(RuntimeError):
            view(water, backend="auto")
        with pytest.raises(RuntimeError):
            view(water, backend="amsview")

        # Given backend cache with successful backend
        view._backends["ase_plot"] = (_AsePlotViewBackend(), True, None)

        # When view
        # Then succeeds
        view(water, backend="auto")


class TestAmsViewBackend:

    backend = _AmsViewBackend

    def test_check_available(self, monkeypatch):
        with patch("scm.plams.tools.view.run_with_timeout") as mock_run_with_timeout, patch(
            "subprocess.run"
        ) as mock_run:
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
                "foo.in -save bar.png -transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 -1.000000 -fixedatomsize -hideregions -showunitcell thickness 0.05 -batch",
            ),
            (
                ViewConfig(width=100, height=100, normal=(1.0, 0.0, 0.0)),
                "foo.in -save bar.png -transparent -scmgeometry 100x100 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 1.000000 0.000000 0.000000 -fixedatomsize -hideregions -showunitcell thickness 0.05 -batch",
            ),
            (
                ViewConfig(
                    fixed_atom_size=False,
                    show_atom_labels=True,
                    atom_label_color="#FFFFFF",
                    atom_label_size=2,
                    show_regions=True,
                ),
                "foo.in -save bar.png -transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 -1.000000 -atomlabel AtomType -labelcolor #FFFFFF -labelsize 2 -showunitcell thickness 0.05 -batch",
            ),
            (
                ViewConfig(show_unit_cell_edges=True, unit_cell_edge_thickness=0.2),
                "foo.in -save bar.png -transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 -1.000000 -fixedatomsize -hideregions -showunitcell thickness 0.2 -batch",
            ),
            (
                ViewConfig(show_unit_cell_faces=True, show_lattice_vectors=True),
                "foo.in -save bar.png -transparent -scmgeometry 800x400 -dpi 300 -padding 0.000000 -showlatticevectors 1 -viewplane 0.000000 0.000000 -1.000000 -fixedatomsize -hideregions -showunitcell faces -batch",
            ),
            (
                ViewConfig(dpi=600),
                "foo.in -save bar.png -transparent -scmgeometry 800x400 -dpi 600 -padding 0.000000 -showlatticevectors 0 -viewplane 0.000000 0.000000 -1.000000 -fixedatomsize -hideregions -showunitcell thickness 0.05 -batch",
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
            assert actual == expected

    @pytest.mark.parametrize(
        "direction, lattice, expected",
        [
            ("along_x", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "-1.000000 0.000000 0.000000"),
            ("along_y", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 -1.000000 0.000000"),
            ("along_z", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 0.000000 -1.000000"),
            ("along_a", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "-0.707107 -0.707107 0.000000"),
            ("along_b", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "-0.707107 0.707107 0.000000"),
            ("along_c", [[5, 5, 0], [5, -5, 0], [0, 0, 10]], "0.000000 0.000000 -1.000000"),
            ("tilt_x", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-0.990148 0.099015 0.099015"),
            ("small_tilt_y", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.049875 -0.997509 0.049875"),
            ("large_tilt_z", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.192450 0.192450 -0.962250"),
            ("tilt_a", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "-0.975055 -0.121556 0.185721"),
            ("small_tilt_b", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.563589 -0.824928 0.043150"),
            ("large_tilt_c", [[10, 2, -1.0], [-5, 8, 0], [0, -2, 11]], "0.082627 0.359045 -0.929656"),
            ("corner_x", [[1, 0, 0]], "-0.577350 0.577350 0.577350"),
            ("corner_y", [[5, 5, 0], [5, -5, 0]], "0.577350 -0.577350 0.577350"),
            ("corner_z", [[5, 5, 0], [5, -5, 0]], "0.577350 0.577350 -0.577350"),
            ("corner_a", [[1, 0, 0]], "-0.577350 0.577350 0.577350"),
            ("corner_b", [[5, 5, 0], [5, -5, 0]], "0.000000 0.816497 0.577350"),
            ("corner_c", [[5, 5, 0], [5, -5, 0]], "0.816497 0.000000 -0.577350"),
        ],
    )
    def test_get_view_plane_with_happy_direction(self, direction, lattice, expected):
        mols = []
        box = Molecule()
        box.lattice = lattice
        mols.append(box)
        if _has_scm_chemsys:
            box_cs = ChemicalSystem()
            box_cs.lattice = Lattice(np.array(lattice))
            mols.append(box_cs)

        for mol in mols:
            actual = self.backend.get_view_plane(mol, ViewConfig(direction=direction))
            assert actual == expected

    @pytest.mark.parametrize(
        "direction", ["p", "xx", "_x", "along_--x", "along__x", "view_x", "view_+x", "small_along_x", "med_tilt_x"]
    )
    def test_get_view_plane_with_unhappy_direction(self, direction, water):
        with pytest.raises(ValueError):
            self.backend.get_view_plane(water, ViewConfig(direction=direction))

    def test_generate_image(self, water):
        skip_if_no_ams_installation()

        img = self.backend.generate_image(water, ViewConfig())

        assert img is not None


class TestAmsViewXvfbBackend(TestAmsViewBackend):

    backend = _AmsViewXvfbBackend

    def test_generate_image(self, water):
        pytest.skip("Skipping as Xvfb not installed")


class TestXvfbManager:

    @pytest.fixture(autouse=True)
    def reset_singleton(self):
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
        with pytest.raises(RuntimeError):
            manager.check_xvfb()

        with patch("shutil.which") as mock_which:
            mock_which.return_value = "foo/Xvfb"

            # When Xvfb cannot be run
            # Then check fails
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

        # When start
        # Then all attempts fail
        with pytest.raises(RuntimeError):
            manager.start()

        # Given manager with succeeding xvfb command
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

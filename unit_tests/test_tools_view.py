import os
import PIL
import pytest
from unittest.mock import patch, MagicMock
import numpy as np

from scm.plams.tools.view import view, ViewConfig, _get_view_plane, _XvfbManager
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
    def test_view_command_passed_to_amsview(self, view_config, expected, water, monkeypatch):

        with patch("scm.plams.tools.view.run_with_timeout") as mock_run_with_timeout, patch(
            "subprocess.run"
        ) as mock_run:
            # This call will fail to generate the image due to the mock
            # but we are just testing that the command to AMSview is generated properly
            response = MagicMock()
            response.stderr = None
            response.stdout = "release=2025.204"
            mock_run.return_value = response
            monkeypatch.setenv("AMSBIN", "dummy_value")
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

                mock_ret.stderr = b"""\
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
            mock_run_ret.stderr = b"-displayfd fd"
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
            mock_run_ret.stderr = b"-displayfd fd"
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
            mock_run_ret.stderr = b"-displayfd fd"
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

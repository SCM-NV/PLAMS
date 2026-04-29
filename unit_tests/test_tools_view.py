import os
import subprocess
import pytest
import threading
import time
from unittest.mock import patch, MagicMock
import numpy as np
from PIL import Image as PilImage

from scm.plams.interfaces.adfsuite.errors import AMSExecutionError
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.tools.view import (
    view,
    ViewConfig,
    _AMSViewManager,
    _AmsViewBackend,
    _AmsViewXvfbBackend,
    _XvfbManager,
    _AsePlotBackend,
    _LazyViewBackend,
)
from scm.plams.mol.molecule import Molecule
from test_helpers import skip_if_windows

try:
    from scm.base import ChemicalSystem, Lattice

    _has_scm_chemsys = True
except ImportError:
    _has_scm_chemsys = False


@pytest.fixture
def water(xyz_folder):
    water = from_smiles("O")
    water.guess_bonds()
    return water


class TestLazyViewBacked:

    @pytest.fixture
    def backend(self):
        mock = MagicMock()
        mock.check_available.return_value = True
        return _LazyViewBackend(mock)

    def test_new_lazy_view_backend_uninitialized(self, backend):
        assert backend._available is None
        assert backend._error is None

    def test_lazy_backend_checked_available_only_once(self, backend):
        assert backend.is_available()
        assert backend.is_available()
        assert backend.backend.check_available.call_count == 1


class TestView:

    @staticmethod
    def get_new_view_backends_cache():
        return {
            "amsview": _LazyViewBackend(_AmsViewBackend()),
            "amsview_xvfb": _LazyViewBackend(_AmsViewXvfbBackend()),
            "ase_plot": _LazyViewBackend(_AsePlotBackend()),
        }

    def test_backends_cache2(self, water):
        # Given backend cache
        import scm.plams.tools.view as viewer

        viewer._view_backends_cache = self.get_new_view_backends_cache()

        def check_available_errors():
            raise RuntimeError("Something went wrong.")

        # When view with auto and first backend available
        with patch.object(_AmsViewBackend, "check_available", return_value=True):
            with patch.object(_AmsViewBackend, "generate_image", return_value=MagicMock()):
                view(water, backend="auto")

            # Then only first backend checked
            assert viewer._view_backends_cache["amsview"]._available
            assert viewer._view_backends_cache["amsview_xvfb"]._available is None
            assert viewer._view_backends_cache["ase_plot"]._available is None

        # When view with a specific later backend
        with patch.object(_AsePlotBackend, "check_available", return_value=True):
            with patch.object(_AsePlotBackend, "generate_image", return_value=MagicMock()):
                view(water, backend="ase_plot")

            # Then only specific backend checked
            assert viewer._view_backends_cache["amsview"]._available
            assert viewer._view_backends_cache["amsview_xvfb"]._available is None
            assert viewer._view_backends_cache["ase_plot"]._available

        # When view with auto and first backend unavailable
        viewer._view_backends_cache = self.get_new_view_backends_cache()

        with patch.object(_AmsViewBackend, "check_available", side_effect=check_available_errors):
            with patch.object(_AmsViewXvfbBackend, "check_available", return_value=True):
                with patch.object(_AmsViewXvfbBackend, "generate_image", return_value=MagicMock()):
                    view(water, backend="auto")

            # Then backends checked until one is available
            assert not viewer._view_backends_cache["amsview"]._available
            assert viewer._view_backends_cache["amsview_xvfb"]._available
            assert viewer._view_backends_cache["ase_plot"]._available is None

            # Then unavailable backend errors
            with pytest.raises(RuntimeError):
                view(water, backend="amsview")

        # When view with auto and no backends available
        viewer._view_backends_cache = self.get_new_view_backends_cache()
        with patch.object(_AmsViewBackend, "check_available", side_effect=check_available_errors):
            with patch.object(_AmsViewXvfbBackend, "check_available", side_effect=check_available_errors):
                with patch.object(_AsePlotBackend, "check_available", side_effect=check_available_errors):
                    # Then auto backend errors
                    with pytest.raises(RuntimeError):
                        view(water, backend="auto")

                    # Then specific backend errors
                    with pytest.raises(RuntimeError):
                        view(water, backend="amsview")

                    # Then all backend checkeds
                    assert not viewer._view_backends_cache["amsview"]._available
                    assert not viewer._view_backends_cache["amsview_xvfb"]._available
                    assert not viewer._view_backends_cache["ase_plot"]._available

        viewer._view_backends_cache = self.get_new_view_backends_cache()


class TestAmsViewBackend:

    backend = _AmsViewBackend

    def test_check_available(self, monkeypatch):
        monkeypatch.setenv("AMSBIN", "test/ams")
        with patch.object(_AMSViewManager, "_generate_image", side_effect=RuntimeError("something went wrong")):
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

    def test_open_window_saves_with_manager_then_opens_window(self, water, tmp_path):
        if self.backend is not _AmsViewBackend:
            pytest.skip("open_window is disabled for this backend")

        input_path = tmp_path / "water.in"
        input_path.write_text("")
        image = MagicMock()

        with patch.object(_AMSViewManager, "_generate_image", return_value=image) as mock_generate, patch.object(
            self.backend, "write_system_input", return_value=str(input_path)
        ), patch.object(self.backend, "run_command") as mock_run_command:
            actual = self.backend.generate_image(water, ViewConfig(open_window=True))

        assert actual is image
        save_config = mock_generate.call_args.args[1]
        assert not save_config.open_window
        assert save_config.timeout == 10
        command, open_config = mock_run_command.call_args.args
        assert open_config.open_window
        assert "-save" not in command

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


class _FakeStdin:

    def __init__(self, fail_on_write=False):
        self.commands = []
        self.closed = False
        self.fail_on_write = fail_on_write

    def write(self, data):
        if self.fail_on_write:
            raise BrokenPipeError("broken pipe")
        self.commands.append(data)

    def flush(self):
        pass

    def close(self):
        self.closed = True


class _FakeProcess:

    def __init__(self, returncode=None, fail_on_write=False, wait_raises=False):
        self.returncode = returncode
        self.stdin = _FakeStdin(fail_on_write=fail_on_write)
        self.terminated = False
        self.killed = False
        self.wait_raises = wait_raises

    def poll(self):
        return self.returncode

    def terminate(self):
        self.terminated = True

    def kill(self):
        self.killed = True
        self.returncode = -9

    def wait(self, timeout=None):
        if self.wait_raises:
            raise subprocess.TimeoutExpired("amsview", timeout)
        self.returncode = 0


class TestAMSViewManager:

    @pytest.fixture(autouse=True)
    def close_singleton(self):
        _AMSViewManager.close_instance()
        yield
        _AMSViewManager.close_instance()

    def test_singleton(self):
        assert _AMSViewManager() is _AMSViewManager()

    def test_reuses_process(self, water):
        manager = _AMSViewManager()
        proc = _FakeProcess()
        image = MagicMock()

        with patch("subprocess.Popen", return_value=proc) as mock_popen, patch.object(
            _AMSViewManager, "_wait_for_image"
        ), patch.object(_AmsViewBackend, "load_and_resize_image", return_value=image):
            assert manager._generate_image(water, ViewConfig(timeout=1)) is image
            assert manager._generate_image(water, ViewConfig(timeout=1)) is image

        assert mock_popen.call_count == 1
        assert proc.stdin.commands[0].startswith('"')
        assert "-save" in proc.stdin.commands[0]
        assert len(proc.stdin.commands) == 2

    def test_restarts_dead_process_on_next_call(self, water):
        manager = _AMSViewManager()
        dead_proc = _FakeProcess(returncode=1)
        live_proc = _FakeProcess()

        with patch("subprocess.Popen", side_effect=[dead_proc, live_proc]) as mock_popen, patch.object(
            _AMSViewManager, "_wait_for_image"
        ), patch.object(_AmsViewBackend, "load_and_resize_image", return_value=MagicMock()):
            manager._ensure_started()
            manager._generate_image(water, ViewConfig(timeout=1))

        assert mock_popen.call_count == 2
        assert len(live_proc.stdin.commands) == 1

    def test_broken_stdin_closes_and_raises(self, water):
        manager = _AMSViewManager()
        broken_proc = _FakeProcess(fail_on_write=True)

        with patch("subprocess.Popen", return_value=broken_proc) as mock_popen, patch.object(
            _AMSViewManager, "_wait_for_image"
        ):
            with pytest.raises(AMSExecutionError):
                manager._generate_image(water, ViewConfig(timeout=1))

        assert mock_popen.call_count == 1
        assert broken_proc.terminated
        assert manager._proc is None

    def test_timeout_closes_and_raises(self, water):
        manager = _AMSViewManager()
        proc = _FakeProcess()

        with patch("subprocess.Popen", return_value=proc) as mock_popen, patch.object(
            _AMSViewManager, "_wait_for_image", side_effect=TimeoutError("timeout")
        ):
            with pytest.raises(AMSExecutionError):
                manager._generate_image(water, ViewConfig(timeout=1))

        assert mock_popen.call_count == 1
        assert proc.terminated
        assert manager._proc is None

    def test_close_terminates_then_kills_if_needed(self):
        manager = _AMSViewManager()
        proc = _FakeProcess(wait_raises=True)
        manager._proc = proc

        manager.close(timeout=0.01)
        manager.close(timeout=0.01)

        assert proc.stdin.closed
        assert proc.terminated
        assert proc.killed
        assert manager._proc is None

    def test_close_allows_later_restart(self, water):
        manager = _AMSViewManager()
        first_proc = _FakeProcess()
        second_proc = _FakeProcess()

        with patch("subprocess.Popen", side_effect=[first_proc, second_proc]) as mock_popen, patch.object(
            _AMSViewManager, "_wait_for_image"
        ), patch.object(_AmsViewBackend, "load_and_resize_image", return_value=MagicMock()):
            manager._generate_image(water, ViewConfig(timeout=1))
            manager.close()
            manager._generate_image(water, ViewConfig(timeout=1))

        assert mock_popen.call_count == 2
        assert first_proc.terminated
        assert len(second_proc.stdin.commands) == 1

    def test_context_manager_closes(self):
        proc = _FakeProcess()
        with patch("subprocess.Popen", return_value=proc):
            with _AMSViewManager() as manager:
                manager._ensure_started()

        assert proc.terminated

    def test_thread_safe_generate_image(self, tmp_path):
        manager = _AMSViewManager()
        active = 0
        max_active = 0
        temp_paths = [tmp_path / f"{idx}.in" for idx in range(3)]
        for path in temp_paths:
            path.write_text("")

        def send_command(command):
            nonlocal active, max_active
            active += 1
            max_active = max(max_active, active)
            time.sleep(0.05)
            active -= 1

        with patch.object(
            _AmsViewBackend, "write_system_input", side_effect=[str(path) for path in temp_paths]
        ), patch.object(
            _AmsViewBackend, "get_image_path", return_value=(str(tmp_path / "image.png"), False)
        ), patch.object(
            manager, "_send_command", side_effect=send_command
        ), patch.object(
            manager, "_wait_for_image"
        ), patch.object(
            _AmsViewBackend, "load_and_resize_image", return_value=MagicMock()
        ):
            threads = [
                threading.Thread(
                    target=manager._generate_image,
                    args=(MagicMock(), ViewConfig(timeout=1, normal=(0.0, 0.0, 1.0))),
                )
                for _ in range(3)
            ]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()

        assert max_active == 1

    def test_send_command_quotes_arguments(self):
        manager = _AMSViewManager()
        proc = _FakeProcess()
        manager._proc = proc

        manager._send_command(
            ["/ams/bin/amsview", "path with spaces/file.in", "-viewplane", "0 0 1", "-labelcolor", "#00AAFF"]
        )

        assert proc.stdin.commands[0] == '"path with spaces/file.in" "-viewplane" "0 0 1" "-labelcolor" "#00AAFF"\n'

    def test_wait_for_image_accepts_fast_completed_render(self, tmp_path):
        manager = _AMSViewManager()
        manager._proc = _FakeProcess()
        img_path = tmp_path / "image.png"
        PilImage.new("RGB", (1, 1)).save(img_path)

        manager._wait_for_image(str(img_path), ViewConfig(timeout=1), None)

    def test_wait_for_image_rejects_unchanged_existing_picture_path(self, tmp_path):
        manager = _AMSViewManager()
        manager._proc = _FakeProcess()
        img_path = tmp_path / "image.png"
        PilImage.new("RGB", (1, 1)).save(img_path)
        previous_stat = os.stat(img_path)

        with pytest.raises(TimeoutError):
            manager._wait_for_image(str(img_path), ViewConfig(timeout=0), previous_stat)

    def test_wait_for_image_waits_until_image_is_readable(self, tmp_path):
        manager = _AMSViewManager()
        manager._proc = _FakeProcess()
        img_path = tmp_path / "image.png"
        img_path.write_bytes(b"not yet a png")
        calls = [0]

        def make_readable(path):
            calls[0] += 1
            if calls[0] == 1:
                return False
            PilImage.new("RGB", (1, 1)).save(path)
            return True

        with patch.object(manager, "_is_readable_image", side_effect=make_readable):
            manager._wait_for_image(str(img_path), ViewConfig(timeout=1), None)

    def test_close_instance_closes_singleton(self):
        manager = _AMSViewManager()

        with patch.object(manager, "close") as mock_close:
            _AMSViewManager.close_instance()

        mock_close.assert_called_once()


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

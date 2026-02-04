import os
import re
import subprocess
from typing import (
    Optional,
    Tuple,
    Union,
    TYPE_CHECKING,
    Literal,
    Sequence,
    List,
    Dict,
    Generator,
    Any,
    TypeVar,
    Type,
    ClassVar,
)
from typing_extensions import Self

import numpy as np
from dataclasses import dataclass, replace
from threading import Lock
import select
import time
from contextlib import contextmanager
import atexit
import shutil
from abc import ABC, abstractmethod
from tempfile import NamedTemporaryFile

from scm.plams.core.functions import requires_optional_package, log
from scm.plams.interfaces.adfsuite.errors import AMSExecutionError
from scm.plams.interfaces.adfsuite.utils import requires_ams
from scm.plams.mol.molecule import Molecule
from scm.plams.core.private import run_with_timeout
from scm.plams.tools.units import Units

try:
    from scm.base import ChemicalSystem

    _has_scm_chemsys = True
except ImportError:
    _has_scm_chemsys = False

if TYPE_CHECKING:
    from PIL import Image as PilImage

TBackend = TypeVar("TBackend", bound="_ViewBackend")

__all__ = ["ViewConfig", "view"]

ViewDirections = Literal[
    "along_x",
    "along_y",
    "along_z",
    "along_a",
    "along_b",
    "along_c",
    "along_pca1",
    "along_pca2",
    "along_pca3",
    "tilt_x",
    "tilt_y",
    "tilt_z",
    "tilt_a",
    "tilt_b",
    "tilt_c",
    "tilt_pca1",
    "tilt_pca2",
    "tilt_pca3",
    "small_tilt_x",
    "small_tilt_y",
    "small_tilt_z",
    "small_tilt_a",
    "small_tilt_b",
    "small_tilt_c",
    "small_tilt_pca1",
    "small_tilt_pca2",
    "small_tilt_pca3",
    "large_tilt_x",
    "large_tilt_y",
    "large_tilt_z",
    "large_tilt_a",
    "large_tilt_b",
    "large_tilt_c",
    "large_tilt_pca1",
    "large_tilt_pca2",
    "large_tilt_pca3",
    "corner_x",
    "corner_y",
    "corner_z",
    "corner_a",
    "corner_b",
    "corner_c",
    "corner_pca1",
    "corner_pca2",
    "corner_pca3",
]

Backends = Literal["amsview", "amsview_xvfb", "ase_plot", "auto"]


@dataclass
class ViewConfig:
    """
    Configuration for view settings

    :param width: width of the image in pixels, defaults to ``800``
    :param height: height of the image in pixels, defaults to ``400``
    :param padding: padding around system in Angstrom, defaults to ``0.0`` (can be negative)
    :param direction: direction to view system along, selected from a series of preset values, defaults to ``along_z``
    :param normal: orientation of the normal to the view plane, takes precedence over direction when specified, defaults to ``None``
    :param normal_basis: whether to use cartesian axes, ``xyz``, lattice vectors (where applicable), ``abc``, or principal component analysis vectors ``pca``, as the basis for the normal to the view plane, defaults to ``xyz``
    :param dpi: resolution of any saved image in dots per inch, defaults to ``300``
    :param picture_path: optional path for the location to save the generated image file, defaults to ``None``
    :param fixed_atom_size: use the same radius for all elements (except Hydrogen), defaults to ``True``
    :param show_atom_labels: display text label on each atom, defaults to ``False``
    :param atom_label_type: property used for atom labels, defaults to ``Element``
    :param atom_label_color: hexadecimal color code for atom labels, defaults to ``#000000`` i.e. black
    :param atom_label_size: scale atom labels by the given factor, to make them larger or smaller, defaults to ``1.0``
    :param guess_bonds: guess bonds before viewing, defaults to ``False``
    :param show_regions: display translucent spheres on atoms according to their regions, defaults to ``False``
    :param show_unit_cell_edges: display unit cell for periodic systems using semi-transparent edges, defaults to ``True``
    :param unit_cell_edge_thickness: specify thickness of the displayed unit cell boundary, defaults to ``0.05``
    :param show_unit_cell_faces: display unit cell for periodic systems using semi-transparent faces, defaults to ``False``
    :param show_lattice_vectors: display the lattice vectors for periodic systems, defaults to ``False``
    :param backend: program to use as a backend to generate images, defaults to ``auto`` i.e. any available program
    :param timeout: kill visualization process after given time in seconds, defaults to ``10`` if window is not opened, otherwise no limit
    :param open_window: open AMSview in a dedicated window if ``True``, otherwise render image offscreen, defaults to ``False``
    """

    # Image/viewpoint
    width: int = 800
    height: int = 400
    padding: float = 0.0
    direction: Optional[ViewDirections] = "along_z"
    normal: Optional[Tuple[float, float, float]] = None
    normal_basis: Literal["xyz", "abc", "pca"] = "xyz"

    # Picture
    dpi: int = 300
    picture_path: Optional[Union[str, os.PathLike]] = None

    # Atom/molecule/bond etc. representation
    fixed_atom_size: bool = True
    show_atom_labels: bool = False
    atom_label_type: Literal["Element", "AtomType", "Name"] = "Element"
    atom_label_color: str = "#000000"
    atom_label_size: float = 1.0
    guess_bonds: bool = False
    show_regions: bool = False

    # Periodic
    show_unit_cell_edges: bool = True
    unit_cell_edge_thickness: float = 0.05
    show_unit_cell_faces: bool = False
    show_lattice_vectors: bool = False

    # Program
    backend: Literal[Backends] = "auto"
    timeout: Optional[int] = None
    open_window: bool = False

    def __post_init__(self) -> None:
        if self.timeout is None:
            self.timeout = 10 if not self.open_window else None

    def validate(self) -> None:
        """
        Check if config values are valid for AMSview

        :raises ValueError: if any value in the config is invalid
        """
        if not isinstance(self.width, int) or self.width < 0:
            raise ValueError(f"width must be a positive integer, but was '{self.width}'")
        if not isinstance(self.height, int) or self.height < 0:
            raise ValueError(f"height must be a positive integer, but was '{self.height}'")
        if not isinstance(self.padding, (int, float)):
            raise ValueError(f"padding must be a numeric value, but was '{self.padding}'")
        if self.direction and (
            not isinstance(self.direction, str) or self.direction not in ViewDirections.__args__  # type: ignore[attr-defined]
        ):
            raise ValueError(
                f"direction must be one of: '{', '.join(ViewDirections.__args__)}'; but was '{self.direction}'"  # type: ignore[attr-defined]
            )
        if self.normal and (
            not isinstance(self.normal, Sequence)
            or len(self.normal) != 3
            or not all(isinstance(v, (int, float)) for v in self.normal)
        ):
            raise ValueError(f"normal must be a sequence of three numeric values, but was '{self.normal}'")
        if not self.direction and not self.normal:
            raise ValueError("direction or normal must be specified")
        if not isinstance(self.normal_basis, str) or self.normal_basis not in ["xyz", "abc", "pca"]:
            raise ValueError(f"normal_basis must be one of: 'xyz', 'ase', 'pca', but was '{self.normal_basis}'")

        if not isinstance(self.dpi, int) or self.dpi < 0:
            raise ValueError(f"dpi must be a positive integer, but was '{self.dpi}'")
        if self.picture_path and not isinstance(self.picture_path, (str, os.PathLike)):
            raise ValueError(f"picture_path must be string or pathlike, but was '{self.picture_path}'")

        if not isinstance(self.fixed_atom_size, bool):
            raise ValueError(f"fixed_atom_size must be a boolean value, but was '{self.fixed_atom_size}'")
        if not isinstance(self.show_atom_labels, bool):
            raise ValueError(f"show_atom_labels must be a boolean value, but was '{self.show_atom_labels}'")
        if not isinstance(self.atom_label_type, str) or self.atom_label_type not in [
            "Element",
            "AtomType",
            "Name",
        ]:
            raise ValueError(
                f"atom_label_type must be one of: 'Element', 'AtomType', 'Name', but was '{self.atom_label_type}'"
            )
        if not isinstance(self.atom_label_color, str) or not bool(
            re.fullmatch(r"#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})", self.atom_label_color)
        ):
            raise ValueError(
                f"atom_label_color must be a color hex code (starting with #), but was '{self.atom_label_color}'"
            )
        if not isinstance(self.atom_label_size, (int, float)):
            raise ValueError(f"atom_label_size must be a numeric value, but was '{self.atom_label_size}'")
        if not isinstance(self.show_regions, bool):
            raise ValueError(f"show_regions must be a boolean value, but was '{self.show_regions}'")

        if not isinstance(self.show_unit_cell_edges, bool):
            raise ValueError(f"show_unit_cell_edges must be a boolean value, but was '{self.show_unit_cell_edges}'")
        if not isinstance(self.unit_cell_edge_thickness, (int, float)) or self.unit_cell_edge_thickness < 0:
            raise ValueError(
                f"unit_cell_edge_thickness must be a positive numeric value, but was '{self.unit_cell_edge_thickness}'"
            )
        if not isinstance(self.show_unit_cell_faces, bool):
            raise ValueError(f"show_unit_cell_faces must be a boolean value, but was '{self.show_unit_cell_faces}'")
        if not isinstance(self.show_lattice_vectors, bool):
            raise ValueError(f"show_lattice_vectors must be a boolean value, but was '{self.show_lattice_vectors}'")

        if not isinstance(self.backend, str) or self.backend not in Backends.__args__:  # type: ignore[attr-defined]
            raise ValueError(
                f"backend must be one of: '{', '.join(Backends.__args__)}'; but was '{self.backend}'"  # type: ignore[attr-defined]
            )
        if self.timeout and (not isinstance(self.timeout, int) or self.timeout < 0):
            raise ValueError(f"timeout must be a positive integer, but was '{self.timeout}'")
        if not isinstance(self.open_window, bool):
            raise ValueError(f"open_window must be a boolean value, but was '{self.open_window}'")


_view_backends_cache: Optional[Dict[str, Tuple["_ViewBackend", bool, Optional[Exception]]]] = None


@requires_optional_package("PIL")
def view(
    system: Union[Molecule, "ChemicalSystem"],
    config: Optional[ViewConfig] = None,
    *,
    width: Optional[int] = None,
    height: Optional[int] = None,
    padding: Optional[float] = None,
    direction: Optional[ViewDirections] = None,
    fixed_atom_size: Optional[bool] = None,
    show_atom_labels: Optional[bool] = None,
    atom_label_type: Optional[Literal["Element", "AtomType", "Name"]] = None,
    guess_bonds: Optional[bool] = None,
    show_regions: Optional[bool] = None,
    show_unit_cell_edges: Optional[bool] = None,
    show_lattice_vectors: Optional[bool] = None,
    picture_path: Optional[Union[str, os.PathLike]] = None,
    backend: Optional[Backends] = None,
    open_window: Optional[bool] = None,
) -> "PilImage.Image":
    """
    View a chemical system or molecule in a Jupyter notebook by generating an image using AMSview/ASE

    :param system: molecule or chemical system to visualize
    :param config: configuration for view
    :param width: override for width of the image in pixels
    :param height: override for height of the image in pixels
    :param padding: override for padding around system in Angstrom
    :param direction: override for direction to view system along
    :param fixed_atom_size: override to use the same radius for all elements (except Hydrogen)
    :param show_atom_labels: override to display text label on each atom
    :param atom_label_type: override for property used for atom labels
    :param guess_bonds: override for guessing bonds before viewing
    :param show_regions: override to display translucent spheres on atoms according to their regions
    :param show_unit_cell_edges: override to display unit cell for periodic systems using semi-transparent edges
    :param show_lattice_vectors: override to display the lattice vectors for periodic systems
    :param picture_path: override for path for the location to save the generated image file
    :param backend: override for program to use as a backend to generate images
    :param open_window: override to open AMSview in a dedicated window
    :return: image of the molecule generated using AMSView
    """
    global _view_backends_cache
    # Set up config objects, applying any config overrides from the keyword args
    config = config or ViewConfig()
    if width is not None:
        config.width = width
    if height is not None:
        config.height = height
    if padding is not None:
        config.padding = padding
    if direction is not None:
        config.direction = direction

    if fixed_atom_size is not None:
        config.fixed_atom_size = fixed_atom_size
    if show_atom_labels is not None:
        config.show_atom_labels = show_atom_labels
    if atom_label_type is not None:
        config.atom_label_type = atom_label_type
    if guess_bonds is not None:
        config.guess_bonds = guess_bonds
    if show_regions is not None:
        config.show_regions = show_regions

    if show_unit_cell_edges is not None:
        config.show_unit_cell_edges = show_unit_cell_edges
    if show_lattice_vectors is not None:
        config.show_lattice_vectors = show_lattice_vectors

    if picture_path is not None:
        config.picture_path = picture_path

    if backend is not None:
        config.backend = backend
    if open_window is not None:
        config.open_window = open_window
        config.timeout = 10 if not config.open_window else None

    # On first call check which backends are available
    if _view_backends_cache is None:

        def check_backend_available(b: TBackend) -> Tuple[TBackend, bool, Optional[Exception]]:
            try:
                b.check_available()
                return b, True, None
            except Exception as ex:
                return b, False, ex

        backends = {
            "amsview": check_backend_available(_AmsViewBackend()),
            "amsview_xvfb": check_backend_available(_AmsViewXvfbBackend()),
            "ase_plot": check_backend_available(_AsePlotBackend()),
        }
        _view_backends_cache = backends
    else:
        backends = _view_backends_cache

    # On subsequent calls get the available backend
    if config.backend != "auto" and config.backend not in backends:
        raise ValueError(f"View backend '{config.backend}' not recognised")

    if config.backend == "auto":
        available_backends = [v for v in backends.values() if v[1]]
        if not any(available_backends):
            errors = "\n\t".join([f"{k}: {err}" for k, (_, __, err) in backends.items()])
            raise RuntimeError(f"No backends available for view.\nErrors were:\n\t{errors}")
        else:
            selected_backend, _, __ = available_backends[0]
    else:
        selected_backend, available, error = backends[config.backend]
        if not available:
            raise RuntimeError(f"Backend '{config.backend}' not available for view.\nError was: {error}")

    # Validation to help prevent crashing due to bad options
    config.validate()

    if config.guess_bonds:
        system = system.copy()
        system.guess_bonds()

    # Render image with backend
    img = selected_backend.generate_image(system, config)

    return img


class _ViewBackend(ABC):
    """
    Abstract base class for a viewer for a molecule/chemical system
    """

    @classmethod
    def check_available(cls) -> None:
        """
        Check whether this backend is available on the current system, otherwise raise an error
        """
        try:
            he_atom = Molecule(positions=[[0, 0, 0]], numbers=[2])
            cls.generate_image(he_atom, ViewConfig())
        except Exception as ex:
            raise AMSExecutionError("Could not generate test image", ex)

    @classmethod
    @abstractmethod
    def generate_image(cls, system: Union[Molecule, "ChemicalSystem"], config: ViewConfig) -> "PilImage.Image":
        """
        Generate image file for the given system

        :param system: molecule or chemical system to visualize
        :param config: configuration for view
        """

    @staticmethod
    def get_view_plane(system: Union[Molecule, "ChemicalSystem"], config: ViewConfig) -> str:
        """
        Get view plane for system from config as a normal vector

        :param system: molecule or chemical system to visualize
        :param config: configuration for view
        """
        # get preset or explicitly specified normal
        if config.normal:
            normal = np.array(config.normal)
            use_lattice_basis = config.normal_basis == "abc"
            use_pca_basis = config.normal_basis == "pca"
        elif config.direction:
            # N.B. currently AMSview only accepts the normal to the view plane as input
            # this restricts slightly what we can support
            # e.g. 'along_-z' (0, 0, 1) and 'along_z' (0, 0, -1) render the same, hence only z is supported for clarity

            # parse direction
            parts = config.direction.split("_")

            # extract main axis
            if len(parts) == 1 or not all(parts):
                raise ValueError(f"direction '{config.direction}' not recognized")
            else:
                main_view_axis = parts[-1]
                if main_view_axis not in ["x", "y", "z", "a", "b", "c", "pca1", "pca2", "pca3"]:
                    raise ValueError(
                        f"direction '{config.direction}' not recognized: '{main_view_axis}' cannot be parsed"
                    )

            # determine whether we are using cartesian axes or lattice vectors as basis
            use_lattice_basis = main_view_axis in ["a", "b", "c"]
            use_pca_basis = main_view_axis in ["pca1", "pca2", "pca3"]

            # orientate based on keyword and main axis
            keywords = parts[:-1]
            if len(keywords) == 2:
                modifier = keywords[0]
                main_keyword = keywords[1]
            else:
                modifier = None
                main_keyword = keywords[0]

            if (
                len(keywords) > 2
                or (modifier and main_keyword != "tilt")
                or (modifier and modifier not in ["small", "large"])
            ):
                raise ValueError(
                    f"direction '{config.direction}' not recognized: '{'_'.join(keywords)}' cannot be parsed"
                )

            if main_keyword == "along":
                main_value = 1.0
                other_value = 0.0
            elif main_keyword == "tilt":
                main_value = -1.0
                other_value = 0.1 if modifier is None else (0.05 if modifier == "small" else 0.2)
            elif main_keyword == "corner":
                main_value = -1.0
                other_value = 1.0
            else:
                raise ValueError(f"direction '{config.direction}' not recognized: '{main_keyword}' cannot be parsed")

            if main_view_axis == "x" or main_view_axis == "a" or main_view_axis == "pca1":
                normal = np.array([main_value, other_value, other_value])
            elif main_view_axis == "y" or main_view_axis == "b" or main_view_axis == "pca2":
                normal = np.array([other_value, main_value, other_value])
            elif main_view_axis == "z" or main_view_axis == "c" or main_view_axis == "pca3":
                normal = np.array([other_value, other_value, main_value])
            else:
                raise ValueError(f"direction '{config.direction}' not recognized: '{main_view_axis}' cannot be parsed")
        else:
            raise ValueError("direction or normal must be specified")

        # use cartesian basis, lattice vector basis or principal component analysis basis as required
        basis = np.identity(3)
        if use_lattice_basis:
            if isinstance(system, Molecule) and system.lattice:
                for i, vec in enumerate(system.lattice):
                    basis[:, i] = np.array(vec)
            elif _has_scm_chemsys and isinstance(system, ChemicalSystem):
                for i, vec in enumerate(system.lattice.vectors):
                    basis[:, i] = np.array(vec)
        if use_pca_basis:
            skip = False
            if isinstance(system, Molecule):
                coords = system.as_array()
                masses = np.array([1e-3 if m == 0 else m for m in system.get_masses()])  # small correction for 0 masses
            elif _has_scm_chemsys and isinstance(system, ChemicalSystem):
                coords = system.coords
                masses = np.array([1e-3 if at.mass == 0 else at.mass for at in system.atoms])
            else:
                # fall-back to cartesian
                skip = True

            # centre coords on origin
            if not skip and len(coords) > 0 and len(masses) > 0:
                total_mass = masses.sum()
                com = (masses[:, None] * coords).sum(axis=0) / total_mass
                recentred_coords = coords - com

                # calculate covariance matrix and compute eigen vectors and values
                covariance = (recentred_coords * masses[:, None]).T @ recentred_coords / total_mass
                eigen_vals, eigen_vecs = np.linalg.eigh(covariance)

                # sort in order of increasing variance
                order = np.argsort(eigen_vals)[::-1]
                basis = eigen_vecs[:, order]
                if np.linalg.det(basis) < 0:
                    basis[:, 2] *= -1
        basis = basis / np.linalg.norm(basis, axis=0, keepdims=True)

        # convert to cartesian basis and normalize
        normal_cartesian_basis = basis @ normal
        normal_cartesian_basis /= np.linalg.norm(normal_cartesian_basis)

        return " ".join([f"{v:.6f}" for v in normal_cartesian_basis])


class _AmsViewBackend(_ViewBackend):
    """
    Viewer for molecule/chemical system using AMSview backend.
    """

    @classmethod
    @requires_ams(minimum_version="2025.204")
    def check_available(cls) -> None:
        super().check_available()

    @classmethod
    def get_command(
        cls, system: Union[Molecule, "ChemicalSystem"], config: ViewConfig, input_path: str, img_path: str
    ) -> List[str]:
        """
        Generate command for AMSview

        :param system: molecule or chemical system to visualize
        :param config: configuration for view
        :param input_path: path to .in file for system
        :param img_path: path to output image file
        """
        command = [
            os.path.expandvars("$AMSBIN/amsview"),
            input_path,
            "-transparent",
            "-scmgeometry",
            f"{config.width}x{config.height}",
            "-dpi",
            str(config.dpi),
            "-padding",
            f"{Units.convert(config.padding, 'angstrom', 'bohr'):.6f}",
            "-showlatticevectors",
            str(int(config.show_lattice_vectors)),
            "-viewplane",
            cls.get_view_plane(system, config),
        ]
        if config.fixed_atom_size:
            command += ["-fixedatomsize"]
        if not config.show_regions:
            command += ["-hideregions"]
        if config.show_atom_labels:
            command += [
                "-atomlabel",
                config.atom_label_type,
                "-labelcolor",
                config.atom_label_color,
                "-labelsize",
                str(config.atom_label_size),
            ]
        if config.show_unit_cell_faces:
            command += ["-showunitcell", "faces"]
        elif config.show_unit_cell_edges:
            command += ["-showunitcell", f"{config.unit_cell_edge_thickness}"]
        else:
            command += ["-showunitcell", "hide"]

        if not config.open_window:
            command += ["-save", img_path, "-batch"]

        return command

    @classmethod
    def run_command(cls, command: List[str], config: ViewConfig) -> None:
        """
        Execute command
        """
        run_with_timeout(command, timeout=config.timeout)

    @classmethod
    def generate_image(cls, system: Union[Molecule, "ChemicalSystem"], config: ViewConfig) -> "PilImage.Image":
        from PIL import Image as PilImage

        # Write temporary input file
        with NamedTemporaryFile(mode="w", suffix=".in", delete=False) as input_file:
            input_path = input_file.name
            if isinstance(system, Molecule):
                system.writein(input_file)
            elif _has_scm_chemsys and isinstance(system, ChemicalSystem):
                input_file.write(str(system))
            else:
                raise ValueError(
                    f"System must be a PLAMS Molecule or a ChemicalSystem, but was {type(system).__name__}"
                )

        if config.picture_path:
            img_path = str(config.picture_path)
        else:
            with NamedTemporaryFile(mode="wb", suffix=".png", delete=False) as img_file:
                img_path = img_file.name

        # Build and execute command
        command = cls.get_command(system, config, input_path, img_path)
        try:
            # For open-window, we run command twice, once to generate the image and the second to view
            # as both cannot be combined without the window auto-closing
            if config.open_window:
                save_config = replace(config, open_window=False, timeout=10)
                save_command = cls.get_command(system, save_config, input_path, img_path)
                cls.run_command(save_command, save_config)
            cls.run_command(command, config)

            # Open image file and resize, making sure to maintain aspect ratio as AMSView may not generate with precise dimensions
            img = PilImage.open(img_path)
            img_width, img_height = img.size
            aspect_ratio = img_width / img_height
            resized_img = img.resize(
                (config.width, int(np.ceil(config.width / aspect_ratio))),
                resample=PilImage.Resampling.LANCZOS,
                reducing_gap=3.0,
            )
        except subprocess.CalledProcessError as ex:
            raise AMSExecutionError(" ".join(command), ex.stderr)
        finally:
            os.remove(input_path)
            if not config.picture_path:
                os.remove(img_path)

        return resized_img


class _AmsViewXvfbBackend(_AmsViewBackend):

    @classmethod
    def check_available(cls) -> None:
        _XvfbManager.check_xvfb()
        super().check_available()

    @classmethod
    def run_command(cls, command: List[str], config: ViewConfig) -> None:
        env = os.environ.copy()
        env["SCM_OPENGL_SOFTWARE"] = "1"

        manager = _XvfbManager()
        manager.start()

        with manager.session(env=env):
            run_with_timeout(command, timeout=config.timeout, env=env)

    @classmethod
    def generate_image(cls, system: Union[Molecule, "ChemicalSystem"], config: ViewConfig) -> "PilImage.Image":
        # do not open the AMSview window with xvfb, otherwise it will hang
        if config.open_window:
            config = replace(config, open_window=False, timeout=10)
        return super().generate_image(system, config)


class _XvfbManager:
    """
    Singleton class to manage xvfb virtual display server.
    Allows graphical programs (like AMSview) to run on headless linux server.

    See: https://github.com/ponty/PyVirtualDisplay for inspiration
    """

    _instance: ClassVar[Optional[Self]] = None
    _initialized: bool = False
    _lock: ClassVar[Lock] = Lock()
    xvfb: ClassVar[str] = "Xvfb"

    def __new__(cls: Type[Self], *args: Any, **kwargs: Any) -> Self:
        with cls._lock:
            if cls._instance is None:
                cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(
        self,
        size: Tuple[int, int] = (1024, 768),
        color_depth: int = 24,
        startup_timeout: int = 30,
        startup_retries: int = 3,
    ):
        if self._initialized:
            return
        self._screen = 0
        self._size = size
        self._color_depth = color_depth
        self._startup_timeout = startup_timeout
        self._startup_retries = startup_retries
        self._started = False
        self._read_file_descriptor: Optional[int] = None
        self._write_file_descriptor: Optional[int] = None
        self._proc: Optional[subprocess.Popen] = None
        self._stdout = None
        self._stderr = None
        self.display_number: Optional[int] = None
        self._initialized = True

    @classmethod
    def check_xvfb(cls) -> None:
        """
        Check if Xvfb is installed and runnable
        """
        if not shutil.which(cls.xvfb):
            raise RuntimeError("Could not find Xvfb, please install it or add it to the PATH")
        try:
            ret = subprocess.run([cls.xvfb, "-help"], capture_output=True, check=True, text=True)
            helptext = ret.stderr or ret.stdout or ""
            if "-displayfd" not in helptext:
                raise RuntimeError(
                    "Found version of Xvfb does not have '-displayfd' support, please update and try again"
                )
        except (subprocess.CalledProcessError, FileNotFoundError) as ex:
            raise RuntimeError(f"Could not successfully run 'Xvfb -help'. Error was: {ex}")

    def start(self) -> None:
        """
        Starts display. If already started this is a no-op.
        """
        with self._lock:
            if self._started:
                return

            self.check_xvfb()

            log("Starting Xvfb...", 3)

            retry_count = 0
            while True:
                self._read_file_descriptor, self._write_file_descriptor = os.pipe()
                try:
                    self._proc = subprocess.Popen(
                        self._command,
                        pass_fds=[self._write_file_descriptor],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                    )
                    self._await_display_number()
                    break
                except (RuntimeError, TimeoutError, FileNotFoundError) as ex:
                    log(f"Xvfb failed to start. Error was {ex}", 5)
                    time.sleep(0.05)
                    retry_count += 1
                    self._kill()
                    if retry_count >= self._startup_retries:
                        log(
                            f"Xvfb failed to start after {retry_count} retries with a {self._startup_retries}s timeout",
                            3,
                        )
                        raise RuntimeError(
                            f"Xvfb failed to start after {retry_count} retries with a {self._startup_timeout}s timeout. Command was: {self._command}. Last error was: {self._stderr}."
                        )
                finally:
                    os.close(self._read_file_descriptor)
                    os.close(self._write_file_descriptor)

            atexit.register(self._stop)
            self._started = True
            log("Xvfb started", 3)

    def _await_display_number(self) -> None:
        """
        Wait for display number to be written to the file descriptor.
        """
        buffer = b""
        start_time = time.monotonic()

        while True:
            # Poll file descriptor every 0.1s
            ready, _, _ = select.select([self._read_file_descriptor], [], [], 0.1)

            # Check that process did not die before sending anything
            if not self.alive:
                raise RuntimeError(f"Xvfb closed. Command was: {self._command}. Error was: {self._stderr}")

            if self._read_file_descriptor is not None and self._read_file_descriptor in ready:
                chunk = os.read(self._read_file_descriptor, 1024)
                if not chunk:  # EOF
                    break
                buffer += chunk

                # Stop if newline found
                if b"\n" in chunk:
                    buffer = buffer.split(b"\n", 1)[0]
                    break

            # Final timeout to prevent hanging
            if time.monotonic() - start_time >= self._startup_timeout:
                raise TimeoutError(f"Xvfb timed out. Command was: {self._command}. Error was: {self._stderr}")

        self.display_number = int(buffer.decode("ascii", errors="replace"))

    def _kill(self) -> None:
        """
        Kill Xvfb subprocess if it is running
        """
        if self._proc and self.alive:
            try:
                self._proc.kill()
            except (ProcessLookupError, OSError):
                pass

            self._proc.wait()

    @property
    def _command(self) -> List[str]:
        """
        Command to start Xvfb
        """
        return [
            self.xvfb,
            "-nolisten",
            "tcp",
            "-screen",
            f"{self._screen}",
            f"{self._size[0]}x{self._size[1]}x{self._color_depth}",
            "-displayfd",
            str(self._write_file_descriptor),
        ]

    @property
    def alive(self) -> bool:
        """
        Check whether Xvfb process is alive and running
        """
        if not self._proc:
            return False

        rc = self._proc.poll()
        if rc is not None:
            self._stdout, self._stderr = self._proc.communicate()
        return rc is None

    @property
    def display(self) -> Optional[str]:
        """
        Get current display value, if available
        """
        return f":{self.display_number}" if self.display_number else None

    @contextmanager
    def session(self, env: Optional[Dict[str, str]] = None) -> Generator[Dict[str, str], None, None]:
        """
        Exclusively use a running Xvfb virtual server display for a process.
        Takes a lock so the display can only be used in this session context, and sets the ``DISPLAY`` environment variable.
        Uses the provided user environment, or takes a copy of current environment.
        The set of environment variables is then returned, and can be used in a subprocess call.

        :param env: optional user environment, defaults to ``None``
        """
        with self._lock:
            env = env if env else os.environ.copy()

            if not self.display:
                raise RuntimeError(
                    f"Xvfb display is not set, server has{' not' if not self._started else ''} been started"
                )

            env["DISPLAY"] = self.display
            yield env

    def stop(self) -> None:
        """
        Stops display. If not started this is a no-op.
        This is automatically registered to fire ``atexit``, when display is started.
        """
        if not self._started:
            return

        log("Stopping Xvfb...", 3)
        self._stop()
        log("Xvfb stopped.", 3)

    def _stop(self) -> None:
        """
        Stops display without logging. If not started this is a no-op.
        This is automatically registered to fire ``atexit``, when display is started.
        """
        if not self._started:
            return

        self._kill()
        self.display_number = None
        self._started = False


class _AsePlotBackend(_ViewBackend):

    @classmethod
    @requires_optional_package("ase")
    @requires_optional_package("matplotlib")
    @requires_optional_package("scipy")
    def check_available(cls) -> None:
        super().check_available()

    @classmethod
    def generate_image(cls, system: Union[Molecule, "ChemicalSystem"], config: ViewConfig) -> "PilImage.Image":
        from PIL import Image as PilImage
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        import matplotlib.colors as colors
        from ase.visualize.plot import Matplotlib
        from scm.plams.interfaces.molecule.ase import toASE

        # Convert system to ASE atoms
        if isinstance(system, Molecule):
            ase_atoms = toASE(system)
        elif _has_scm_chemsys and isinstance(system, ChemicalSystem):
            ase_atoms = system.to_ase_atoms()
        else:
            raise ValueError(f"System must be a PLAMS Molecule or a ChemicalSystem, but was {type(system).__name__}")

        # Get image path for (temporary) file
        if config.picture_path:
            img_path = str(config.picture_path)
        else:
            with NamedTemporaryFile(mode="wb", suffix=".png", delete=False) as img_file:
                img_path = img_file.name

        # Create axis
        fig, ax = plt.subplots(figsize=(config.width / config.dpi, config.height / config.dpi), dpi=config.dpi)

        # Set options
        rotation = cls.get_view_rotation(system, config)
        radii = [0.31 if at.symbol == "H" else 0.5 for at in ase_atoms] if config.fixed_atom_size else None
        show_unit_cell = 2 if config.show_unit_cell_edges or config.show_lattice_vectors else 0

        try:
            # Equivalent of plot_atoms from ASE, but gives us more flexibility
            import matplotlib.pyplot as plt

            if ax is None:
                ax = plt.gca()
            plotter = Matplotlib(ase_atoms, ax, rotation=rotation, radii=radii, show_unit_cell=show_unit_cell)
            plotter.write()

            # Reduce atom circle outer line width and unit cell thickness
            for patch in ax.patches:
                if isinstance(patch, plt.Circle):
                    patch.set_linewidth(0.2)
                elif isinstance(patch, patches.PathPatch):
                    patch.set_linewidth(config.unit_cell_edge_thickness)

            # Set unit cell line width
            for line in ax.lines:
                line.set_linewidth(config.unit_cell_edge_thickness)

            # Draw lattice vectors
            if config.show_lattice_vectors:
                # Find the corner to place the lattice vectors on, choose the bottom leftmost-corner (after possible rotation)
                vertices = plotter.cell_vertices
                bottom_mask = vertices[:, 1] <= np.median(vertices[:, 1])
                bottom_left_idx = np.where(bottom_mask)[0][np.argmin(vertices[bottom_mask, 0])]
                bottom_left = plotter.cell_vertices[bottom_left_idx, :]

                # Find the lattice displacements for this corner
                # This relies on the fact that we know the following order for the vertices:
                # 0 0 0 o
                # 0 0 1 c
                # 0 1 0 b
                # 0 1 1 b + c
                # 1 0 0 a
                # 1 0 1 a + c
                # 1 1 0 a + b
                # 1 1 1 a + b + c
                bottom_left_displacements = [int(d) for d in format(bottom_left_idx, "03b")]

                # Draw vectors to the three nearest vertices using these displacements
                for i, color in zip(range(3), ["r", "g", "b"]):
                    vertex_displacements = [j for j in bottom_left_displacements]
                    vertex_displacements[i] = 0 if bottom_left_displacements[i] == 1 else 1
                    vertex_idx = int("".join(str(d) for d in vertex_displacements), 2)
                    vertex = plotter.cell_vertices[vertex_idx, :]
                    if not np.allclose(vertex - bottom_left, 0):
                        ax.arrow(
                            bottom_left[0],
                            bottom_left[1],
                            vertex[0] - bottom_left[0],
                            vertex[1] - bottom_left[1],
                            head_width=0,
                            head_length=0,
                            color=color,
                        )

                # Remove unit cell if not required
                if not config.show_unit_cell_edges:
                    to_remove = [p for p in ax.patches if isinstance(p, patches.PathPatch)]
                    for patch in to_remove:
                        patch.remove()

            # Draw unit cell faces
            if config.show_unit_cell_faces:
                verts = plotter.cell_vertices
                faces = [
                    [verts[j, :2] for j in [0, 1, 3, 2]],
                    [verts[j, :2] for j in [4, 5, 7, 6]],
                    [verts[j, :2] for j in [0, 1, 5, 4]],
                    [verts[j, :2] for j in [2, 3, 7, 6]],
                    [verts[j, :2] for j in [1, 3, 7, 5]],
                    [verts[j, :2] for j in [0, 2, 6, 4]],
                ]
                for face in faces:
                    poly = patches.Polygon(face, fill=True, color="grey", alpha=0.05)
                    ax.add_patch(poly)

            # Draw regions
            if config.show_regions:
                if isinstance(system, Molecule):
                    regions = [
                        [at.properties.region] if isinstance(at.properties.region, str) else list(at.properties.region)
                        for at in system
                    ]
                elif _has_scm_chemsys and isinstance(system, ChemicalSystem):
                    regions = [system.get_regions_of_atom(at) for at in system]

                try:
                    import matplotlib.colormaps as colormaps
                except ImportError:
                    import matplotlib.cm as colormaps

                cmap = colormaps.get_cmap("tab10")
                color_counter = 0
                region_cmap: Dict[str, colors.Colormap] = {}
                for patch in ax.patches:
                    if isinstance(patch, plt.Circle):
                        # Patches are not necessarily in the same order as the atoms
                        # so using the centre of the patch we find the correct atom
                        # not super efficient, but probably not using regions on very large systems...
                        x, y = patch.center
                        d_sq = (plotter.positions[:, 0] - x) ** 2 + (plotter.positions[:, 1] - y) ** 2
                        atom_idx = np.argmin(d_sq)
                        if d_sq[atom_idx] < 1e-6:  # arbitrary tolerance
                            atom_regions = regions[atom_idx]
                            for region in atom_regions:
                                if region in region_cmap:
                                    region_color = region_cmap[region]
                                else:
                                    region_color = cmap.colors[color_counter]
                                    region_cmap[region] = region_color
                                    color_counter += 1
                                region_patch = patches.Circle(
                                    (x, y), patch.radius * 1.5, alpha=0.2, color=region_color, linewidth=0
                                )
                                ax.add_patch(region_patch)

            # Add labels
            atom_type_counts: Dict[str, int] = {}
            if config.show_atom_labels:
                for i, at in enumerate(ase_atoms):
                    x, y = plotter.positions[i, :2]
                    l = at.symbol
                    if config.atom_label_type == "Name":
                        if l in atom_type_counts:
                            atom_type_counts[l] += 1
                        else:
                            atom_type_counts[l] = 1
                        l += f"({atom_type_counts[l]})"
                    ax.text(
                        x,
                        y,
                        l,
                        ha="center",
                        va="center",
                        fontsize=4 * config.atom_label_size,
                        color=config.atom_label_color,
                    )

            # Add padding to image
            if config.padding:
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                ax.set_xlim(xlim[0] - config.padding, xlim[1] + config.padding)
                ax.set_ylim(ylim[0] - config.padding, ylim[1] + config.padding)

            ax.axis("off")

            fig.savefig(img_path, dpi=config.dpi, bbox_inches="tight")

            img = PilImage.open(img_path)
            img_width, img_height = img.size
            aspect_ratio = img_width / img_height
            resized_img = img.resize(
                (config.width, int(np.ceil(config.width / aspect_ratio))),
                resample=PilImage.Resampling.LANCZOS,
                reducing_gap=3.0,
            )
        finally:
            plt.close(fig)
            if not config.picture_path:
                os.remove(img_path)

        return resized_img

    @classmethod
    def get_view_rotation(cls, system: Union[Molecule, "ChemicalSystem"], config: ViewConfig) -> str:
        """
        Get view rotation for system from config as Euler angles

        :param system: molecule or chemical system to visualize
        :param config: configuration for view
        """
        import warnings
        from scipy.spatial.transform import Rotation

        normal = np.array([float(v) for v in cls.get_view_plane(system, config).split(" ")])

        # Calculate rotation to z-axis (ASE default), then convert to Euler angles
        xyz = np.eye(3)
        if np.allclose(normal, xyz[0, :]):
            return "-90y"
        elif np.allclose(normal, xyz[0, :]):
            return "90y"
        if np.allclose(normal, xyz[1, :]):
            return "90x"
        elif np.allclose(normal, xyz[1, :]):
            return "-90x"
        elif np.allclose(normal, xyz[2, :]):
            return ""
        elif np.allclose(normal, -xyz[2, :]):
            return "180x,180z"

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="Optimal rotation is not uniquely or poorly defined for the given sets of vectors."
            )
            rotation, _ = Rotation.align_vectors([xyz[2, :]], [normal])
            angles = Rotation.from_matrix(rotation.as_matrix()).as_euler("xyz", degrees=True)

        return ",".join(f"{ang:.2f}{ax}" for ang, ax in zip(angles, "xyz"))

import os
import re
from typing import Optional, Tuple, Union, TYPE_CHECKING, Literal, Sequence
import numpy as np
from dataclasses import dataclass

from scm.plams.core.functions import requires_optional_package
from scm.plams.interfaces.adfsuite.utils import requires_ams
from scm.plams.mol.molecule import Molecule
from scm.plams.core.private import run_with_timeout
from scm.plams.tools.units import Units

try:
    from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

    _has_scm_chemsys = True
except ImportError:
    _has_scm_chemsys = False

if TYPE_CHECKING:
    from PIL import Image as PilImage

__all__ = ["ViewConfig", "view"]


@dataclass
class ViewConfig:
    """
    Configuration for view settings for AMSview

    :ivar width: width of the image in pixels, defaults to ``800``
    :ivar height: height of the image in pixels, defaults to ``400``
    :ivar padding: padding around system in Angstrom, defaults to ``0`` (can be negative)
    :ivar direction: direction to view system along, selected from a series of preset values, defaults to ``along_z``
    :ivar normal: orientation of the normal to the view plane, takes precedence over direction when specified, defaults to ``None``
    :ivar lattice_as_basis: whether to use lattice vectors (where available) as the view basis, otherwise uses cartesian axes, defaults to ``False``
    :ivar dpi: resolution of any saved image in dots per inch, defaults to ``300``
    :ivar picture_path: optional path for the location to save the generated image file, defaults to ``None``
    :ivar fixed_atom_size: use the same radius for all elements (except Hydrogen), defaults to ``True``
    :ivar show_atom_labels: display text label on each atom, defaults to ``False``
    :ivar atom_label_type: property used for atom labels, defaults to ``AtomType``
    :ivar atom_label_color: hexadecimal color code for atom labels, defaults to ``#000000`` i.e. black
    :ivar atom_label_size: scale atom labels by the given factor, to make them larger or smaller, defaults to ``1.0``
    :ivar show_regions: display translucent spheres on atoms according to their regions, defaults to ``False``
    :ivar show_unit_cell_edges: display unit cell for periodic systems using semi-transparent edges, defaults to ``True``
    :ivar unit_cell_edge_thickness: specify thickness of the displayed unit cell boundary, defaults to ``0.05``
    :ivar show_unit_cell_faces: display unit cell for periodic systems using semi-transparent faces, defaults to ``False``
    :ivar show_lattice_vectors: display the lattice vectors for periodic systems, defaults to ``False``
    :ivar timeout: kill AMSView process after given time in seconds, defaults to ``10`` if window is not opened, otherwise no limit
    :ivar open_window: open AMSview in a dedicated window if ``True``, otherwise render image offscreen, defaults to ``False``
    """

    # Image/viewpoint
    width: int = 800
    height: int = 400
    padding: int = 0
    direction: Optional[
        Literal[
            "along_x",
            "along_y",
            "along_z",
            "tilt_x",
            "tilt_y",
            "tilt_z",
            "small_tilt_x",
            "small_tilt_y",
            "small_tilt_z",
            "large_tilt_x",
            "large_tilt_y",
            "large_tilt_z",
            "corner_x",
            "corner_y",
            "corner_z",
        ]
    ] = "along_z"
    normal: Optional[Tuple[float, float, float]] = None
    lattice_as_basis: bool = False

    # Picture
    dpi: int = 300
    picture_path: Optional[Union[str, os.PathLike]] = None

    # Atom/molecule/bond etc. representation
    fixed_atom_size: bool = True
    show_atom_labels: bool = False
    atom_label_type: Literal["AtomType", "Element", "Name", "SurfaceRadius"] = "AtomType"
    atom_label_color: str = "#000000"
    atom_label_size: float = 1.0
    show_regions: bool = False

    # Periodic
    show_unit_cell_edges: bool = True
    unit_cell_edge_thickness: float = 0.05
    show_unit_cell_faces: bool = False
    show_lattice_vectors: bool = False

    # Program
    timeout: Optional[int] = None
    open_window: bool = False

    def __post_init__(self):
        if self.timeout is None:
            self.timeout = 10 if not self.open_window else None

    def validate(self):
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
        # further validated in _get_view_plane
        if self.direction and (not isinstance(self.direction, str)):
            raise ValueError(f"direction must be a string value, but was '{self.direction}'")
        if self.normal and (
            not isinstance(self.normal, Sequence)
            or len(self.normal) != 3
            or not all(isinstance(v, (int, float)) for v in self.normal)
        ):
            raise ValueError(f"normal must be a sequence of three numeric values, but was '{self.normal}'")
        if not self.direction and not self.normal:
            raise ValueError(f"direction or normal must be specified")
        if not isinstance(self.lattice_as_basis, bool):
            raise ValueError(f"lattice_as_basis must be a boolean value, but was '{self.lattice_as_basis}'")

        if not isinstance(self.dpi, int) or self.dpi < 0:
            raise ValueError(f"dpi must be a positive integer, but was '{self.dpi}'")
        if self.picture_path and not isinstance(self.picture_path, (str, os.PathLike)):
            raise ValueError(f"picture_path must be string or pathlike, but was '{self.picture_path}'")

        if not isinstance(self.fixed_atom_size, bool):
            raise ValueError(f"fixed_atom_size must be a boolean value, but was '{self.fixed_atom_size}'")
        if not isinstance(self.show_atom_labels, bool):
            raise ValueError(f"show_atom_labels must be a boolean value, but was '{self.show_atom_labels}'")
        if not isinstance(self.atom_label_type, str) or self.atom_label_type not in [
            "AtomType",
            "Element",
            "Name",
            "SurfaceRadius",
        ]:
            raise ValueError(
                f"atom_label_type must be one of: 'AtomType', 'Element', 'Name', 'SurfaceRadius', but was '{self.atom_label_type}'"
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

        if self.timeout and (not isinstance(self.timeout, int) or self.timeout < 0):
            raise ValueError(f"timeout must be a positive integer, but was '{self.timeout}'")
        if not isinstance(self.open_window, bool):
            raise ValueError(f"open_window must be a boolean value, but was '{self.open_window}'")


@requires_optional_package("PIL")
@requires_ams(minimum_version="2025.204")
def view(
    system: Union[Molecule, "ChemicalSystem"],
    config: Optional[ViewConfig] = None,
    *,
    width: Optional[int] = None,
    height: Optional[int] = None,
    padding: Optional[float] = None,
    direction: Optional[Literal["foo"]] = None,
    lattice_as_basis: Optional[bool] = None,
    fixed_atom_size: Optional[bool] = None,
    show_atom_labels: Optional[bool] = None,
    show_regions: Optional[bool] = None,
    show_unit_cell_edges: Optional[bool] = None,
    show_lattice_vectors: Optional[bool] = None,
    picture_path: Optional[Union[str, os.PathLike]] = None,
    open_window: Optional[bool] = None,
) -> "PilImage.Image":
    """
    View a chemical system or molecule in a Jupyter notebook by generating an image using AMSview.

    :param system: molecule or chemical system to visualize
    :param config: configuration for AMSview
    :param width: override for width of the image in pixels
    :param height: override for height of the image in pixels
    :param padding: override for padding around system in Angstrom
    :param direction: override for direction to view system along
    :param lattice_as_basis: override for whether to use lattice vectors (where available) as the view basis
    :param fixed_atom_size: override to use the same radius for all elements (except Hydrogen)
    :param show_atom_labels: override to display text label on each atom
    :param show_regions: override to display translucent spheres on atoms according to their regions
    :param show_unit_cell_edges: override to display unit cell for periodic systems using semi-transparent edges
    :param show_lattice_vectors: override to display the lattice vectors for periodic systems
    :param picture_path: override for path for the location to save the generated image file
    :param open_window: override to open AMSview in a dedicated window
    :return: image of the molecule generated using AMSView
    """
    from tempfile import NamedTemporaryFile
    from PIL import Image as PilImage

    # Set up config objects, applying any config overrides from the keyword args
    config = config or ViewConfig()
    if width:
        config.width = width
    if height:
        config.height = height
    if padding:
        config.padding = padding
    if lattice_as_basis:
        config.lattice_as_basis = lattice_as_basis
    if direction:
        config.direction = direction

    if fixed_atom_size:
        config.fixed_atom_size = fixed_atom_size
    if show_atom_labels:
        config.show_atom_labels = show_atom_labels
    if show_regions:
        config.show_regions = show_regions

    if show_unit_cell_edges:
        config.show_edges = show_unit_cell_edges
    if show_lattice_vectors:
        config.show_lattice_vectors = show_lattice_vectors

    if picture_path:
        config.path = picture_path

    if open_window:
        config.open_window = open_window
        config.timeout = 10 if not config.open_window else None

    # Validation to help prevent AMSView crashing due to bad options
    config.validate()

    # Write temporary input file
    with NamedTemporaryFile(mode="w", suffix=".in", delete=False) as input_file:
        input_path = input_file.name
        if isinstance(system, Molecule):
            system.writein(input_file)
        elif _has_scm_chemsys and isinstance(system, ChemicalSystem):
            input_file.write(str(system))
        else:
            raise ValueError(f"System must be a PLAMS Molecule or a ChemicalSystem, but was {type(system).__name__}")

    if config.picture_path:
        img_path = config.picture_path
    else:
        with NamedTemporaryFile(mode="wb", suffix=".png", delete=False) as img_file:
            img_path = img_file.name

    # Build and execute command
    try:
        command = [
            os.path.expandvars("$AMSBIN/amsview"),
            input_path,
            "-save",
            img_path,
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
            _get_view_plane(system, config),
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
            command += ["-showunitcell", f"thickness {config.unit_cell_edge_thickness}"]
        else:
            command += ["-showunitcell", "hide"]

        if not config.open_window:
            command += ["-batch"]

        env = os.environ.copy()
        env["SCM_OPENGL_SOFTWARE"] = "1"
        run_with_timeout(command, timeout=config.timeout, env=env)

        # Open image file and resize, making sure to maintain aspect ratio as AMSView may not generate with precise dimensions
        img = PilImage.open(img_path)
        img_width, img_height = img.size
        aspect_ratio = img_width / img_height
        img = img.resize(
            (config.width, int(np.ceil(config.width / aspect_ratio))),
            resample=PilImage.Resampling.LANCZOS,
            reducing_gap=3.0,
        )
    finally:
        os.remove(input_path)
        if not config.picture_path:
            os.remove(img_path)

    return img


def _get_view_plane(system: Union[Molecule, "ChemicalSystem"], config: ViewConfig) -> str:
    # use cartesian basis or lattice vector basis as required
    basis = np.identity(3)
    if config.lattice_as_basis:
        if isinstance(system, Molecule) and system.lattice:
            for i, vec in enumerate(system.lattice):
                basis[:, i] = np.array(vec)
        elif _has_scm_chemsys and isinstance(system, ChemicalSystem):
            for i, vec in enumerate(system.lattice.vectors):
                basis[:, i] = np.array(vec)
    basis = basis / np.linalg.norm(basis, axis=0, keepdims=True)

    # get preset or explicitly specified normal
    if config.normal:
        normal = np.array(config.normal)
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
            if main_view_axis not in ["x", "y", "z"]:
                raise ValueError(f"direction '{config.direction}' not recognized: '{main_view_axis}' cannot be parsed")

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
            raise ValueError(f"direction '{config.direction}' not recognized: '{'_'.join(keywords)}' cannot be parsed")

        if main_keyword == "along":
            main_value = -1.0
            other_value = 0.0
        elif main_keyword == "tilt":
            main_value = -1.0
            other_value = 0.1 if modifier is None else (0.05 if modifier == "small" else 0.2)
        elif main_keyword == "corner":
            main_value = -1.0
            other_value = 1.0
        else:
            raise ValueError(f"direction '{config.direction}' not recognized: '{main_keyword}' cannot be parsed")

        if main_view_axis == "x":
            normal = np.array([main_value, other_value, other_value])
        elif main_view_axis == "y":
            normal = np.array([other_value, main_value, other_value])
        elif main_view_axis == "z":
            normal = np.array([other_value, other_value, main_value])
        else:
            raise ValueError(f"direction '{config.direction}' not recognized: '{main_view_axis}' cannot be parsed")
    else:
        raise ValueError(f"direction or normal must be specified")

    # convert to cartesian basis and normalize
    normal_cartesian_basis = basis @ normal
    normal_cartesian_basis /= np.linalg.norm(normal_cartesian_basis)

    return " ".join([f"{v:.6f}" for v in normal_cartesian_basis])

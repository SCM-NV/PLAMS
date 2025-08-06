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

__all__ = ["ViewConfig", "PictureConfig", "RepresentationConfig", "PeriodicConfig", "ProgramConfig", "view"]


@dataclass
class ViewConfig:
    """
    Configuration for view settings for AMSview

    :ivar width: width of the image in pixels, defaults to ``800``
    :ivar height: height of the image in pixels, defaults to ``400``
    :ivar padding: padding around system in Angstrom, defaults to ``0`` (can be negative)
    :ivar direction: direction to view system along, selected from a series of preset values
    :ivar normal: orientation of the normal to the view plane, defaults ``(0, 0, 1)`` i.e. in the x-y plane
    :ivar rotation: rotation angles in degrees around the x, y and z axes, defaults to ``(0, 0, 0)``
    """

    width: int = 800
    height: int = 400
    padding: int = 0
    direction: Optional[Literal["foo"]] = None
    normal: Tuple[float, float, float] = (0.0, 0.0, 1.0)
    rotation: Tuple[float, float, float] = (0.0, 0.0, 0.0)

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
        if self.direction and (not isinstance(self.direction, str) or self.direction not in ["foo"]):
            raise ValueError(f"direction must be one of: 'foo', but was '{self.direction}'")
        if (
            not isinstance(self.normal, Sequence)
            or len(self.normal) != 3
            or not all(isinstance(v, (int, float)) for v in self.normal)
        ):
            raise ValueError(f"normal must be a sequence of three numeric values, but was '{self.normal}'")
        if (
            not isinstance(self.rotation, Sequence)
            or len(self.rotation) != 3
            or not all(isinstance(v, (int, float)) for v in self.rotation)
        ):
            raise ValueError(f"rotation must be a sequence of three numeric values, but was '{self.rotation}'")


@dataclass
class PictureConfig:
    """
    Configuration for the saved picture settings for AMSview

    :ivar dpi: resolution of any saved image in dots per inch, defaults to ``300``
    :ivar path: optional path for the location to save the generated image file, defaults to ``None``
    """

    dpi: int = 300
    path: Optional[Union[str, os.PathLike]] = None

    def validate(self):
        """
        Check if config values are valid for AMSview

        :raises ValueError: if any value in the config is invalid
        """
        if not isinstance(self.dpi, int) or self.dpi < 0:
            raise ValueError(f"dpi must be a positive integer, but was '{self.dpi}'")
        if self.path and not isinstance(self.path, (str, os.PathLike)):
            raise ValueError(f"path must be string or pathlike, but was '{self.path}'")


@dataclass
class RepresentationConfig:
    """
    Configuration for atoms/bonds/molecule representation for AMSview

    :ivar fixed_atom_size: use the same radius for all elements (except Hydrogen), defaults to ``True``
    :ivar show_atom_labels: display text label on each atom, defaults to ``False``
    :ivar atom_label_type: property used for atom labels, defaults to ``AtomType``
    :ivar atom_label_color: hexadecimal color code for atom labels, defaults to ``#000000`` i.e. black
    :ivar atom_label_size: scale atom labels by the given factor, to make them larger or smaller, defaults to ``1.0``
    :ivar show_regions: display translucent spheres on atoms according to their regions, defaults to ``False``
    """

    fixed_atom_size: bool = True
    show_atom_labels: bool = False
    atom_label_type: Literal["AtomType", "Element", "Name", "SurfaceRadius"] = "AtomType"
    atom_label_color: str = "#000000"
    atom_label_size: float = 1.0
    show_regions: bool = False

    def validate(self):
        """
        Check if config values are valid for AMSview

        :raises ValueError: if any value in the config is invalid
        """
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


@dataclass
class PeriodicConfig:
    """
    Configuration for periodic settings for AMSview

    :ivar show_edges: display unit cell for periodic systems using semi-transparent edges, defaults to ``True``
    :ivar edge_thickness: specify thickness of the displayed unit cell boundary, defaults to ``0.05``
    :ivar show_faces: display unit cell for periodic systems using semi-transparent faces, defaults to ``False``
    :ivar show_lattice_vectors: display the lattice vectors for periodic systems, defaults to ``False``
    """

    show_edges: bool = True
    edge_thickness: float = 0.05
    show_faces: bool = False
    show_lattice_vectors: bool = False

    def validate(self):
        """
        Check if config values are valid for AMSview

        :raises ValueError: if any value in the config is invalid
        """
        if not isinstance(self.show_edges, bool):
            raise ValueError(f"show_edges must be a boolean value, but was '{self.show_edges}'")
        if not isinstance(self.edge_thickness, (int, float)) or self.edge_thickness < 0:
            raise ValueError(f"edge_thickness must be a positive numeric value, but was '{self.edge_thickness}'")
        if not isinstance(self.show_faces, bool):
            raise ValueError(f"show_faces must be a boolean value, but was '{self.show_faces}'")
        if not isinstance(self.show_lattice_vectors, bool):
            raise ValueError(f"show_lattice_vectors must be a boolean value, but was '{self.show_lattice_vectors}'")


@dataclass
class ProgramConfig:
    """
    Configuration for AMSview program

    :ivar timeout: kill AMSView process after given time in seconds, defaults to ``10`` if window is not opened, otherwise no limit
    :ivar open_window: open AMSview in an dedicated window if ``True``, otherwise render image offscreen, defaults to ``False``
    """

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
        if self.timeout and (not isinstance(self.timeout, int) or self.timeout < 0):
            raise ValueError(f"timeout must be a positive integer, but was '{self.timeout}'")
        if not isinstance(self.open_window, bool):
            raise ValueError(f"open_window must be a boolean value, but was '{self.open_window}'")


@requires_optional_package("PIL")
@requires_ams(minimum_version="2025.204")
def view(
    system: Union[Molecule, "ChemicalSystem"],
    view: Optional[ViewConfig] = None,
    representation: Optional[RepresentationConfig] = None,
    periodic: Optional[PeriodicConfig] = None,
    picture: Optional[PictureConfig] = None,
    program: Optional[ProgramConfig] = None,
    *,
    width: Optional[int] = None,
    height: Optional[int] = None,
    direction: Optional[Literal["foo"]] = None,
    fixed_atom_size: Optional[bool] = None,
    show_atom_labels: Optional[bool] = None,
    show_regions: Optional[bool] = None,
    show_unit_cell_edges: Optional[bool] = None,
    show_lattice_vectors: Optional[bool] = None,
    picture_path: Optional[Union[str, os.PathLike]] = None,
    open_window: Optional[bool] = None,
) -> "PilImage.Image":
    """
    View a chemical system or molecule in a Jupyter notebook by generating an image using AMSView.

    :param system: molecule or chemical system to visualize
    :param view: configuration for display and viewpoint settings
    :param representation: configuration for how system is represented e.g. atoms, bonds, and molecules are drawn
    :param periodic: configuration for displaying the unit cell and lattice vectors
    :param picture: configuration related to image output properties
    :param program: configuration related to AMSview program settings
    :param width: override for width of the image in pixels
    :param height: override for height of the image in pixels
    :param direction: override for direction to view system along
    :param fixed_atom_size: override to use the same radius for all elements (except Hydrogen)
    :param show_atom_labels: override to display text label on each atom
    :param show_regions: override to display translucent spheres on atoms according to their regions
    :param show_unit_cell_edges: override to display unit cell for periodic systems using semi-transparent edges
    :param show_lattice_vectors: override to display the lattice vectors for periodic systems
    :param picture_path: override for path for the location to save the generated image file
    :return: image of the molecule generated using AMSView
    """
    from tempfile import NamedTemporaryFile
    from PIL import Image as PilImage

    # Set up config objects, applying any config overrides from the keyword args
    view_config = view or ViewConfig()
    if width:
        view_config.width = width
    if height:
        view_config.height = height
    if direction:
        view_config.direction = direction

    repr_config = representation or RepresentationConfig()
    if fixed_atom_size:
        repr_config.fixed_atom_size = fixed_atom_size
    if show_atom_labels:
        repr_config.show_atom_labels = show_atom_labels
    if show_regions:
        repr_config.show_regions = show_regions

    periodic_config = periodic or PeriodicConfig()
    if show_unit_cell_edges:
        periodic_config.show_edges = show_unit_cell_edges
    if show_lattice_vectors:
        periodic_config.show_lattice_vectors = show_lattice_vectors

    picture_config = picture or PictureConfig()
    if picture_path:
        picture_config.path = picture_path

    program_config = program or ProgramConfig()
    if open_window:
        program_config.open_window = open_window
        program_config.timeout = 10 if not program_config.open_window else None

    # Validation to help prevent AMSView crashing due to bad options
    view_config.validate()
    repr_config.validate()
    periodic_config.validate()
    picture_config.validate()
    program_config.validate()

    # Write temporary input file
    with NamedTemporaryFile(mode="w", suffix=".in", delete=False) as input_file:
        input_path = input_file.name
        if isinstance(system, Molecule):
            system.writein(input_file)
        elif _has_scm_chemsys and isinstance(system, ChemicalSystem):
            input_file.write(str(system))
        else:
            raise ValueError(f"System must be a PLAMS Molecule or a ChemicalSystem, but was {type(system).__name__}")

    if picture_config.path:
        img_path = picture_config.path
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
            f"{view_config.width}x{view_config.height}",
            "-dpi",
            str(picture_config.dpi),
            "-padding",
            str(Units.convert(view_config.padding, "angstrom", "bohr")),
            "-showlatticevectors",
            str(int(periodic_config.show_lattice_vectors)),
            "-viewplane",
            " ".join([str(v) for v in view_config.normal]),
        ]
        if repr_config.fixed_atom_size:
            command += ["-fixedatomsize"]
        if not repr_config.show_regions:
            command += ["-hideregions"]
        if repr_config.show_atom_labels:
            command += [
                "-atomlabel",
                repr_config.atom_label_type,
                "-labelcolor",
                repr_config.atom_label_color,
                "-labelsize",
                str(repr_config.atom_label_size),
            ]
        if periodic_config.show_faces:
            command += ["-showunitcell", "faces"]
        elif periodic_config.show_edges:
            command += ["-showunitcell", f"thickness {periodic_config.edge_thickness}"]
        else:
            command += ["-showunitcell", "hide"]

        if not program_config.open_window:
            command += ["-batch"]

        env = os.environ.copy()
        env["SCM_OPENGL_SOFTWARE"] = "1"
        run_with_timeout(command, timeout=program_config.timeout, env=env)

        # Open image file and resize, making sure to maintain aspect ratio as AMSView may not generate with precise dimensions
        img = PilImage.open(img_path)
        img_width, img_height = img.size
        aspect_ratio = img_width / img_height
        img = img.resize(
            (view_config.width, int(np.ceil(view_config.width / aspect_ratio))),
            resample=PilImage.Resampling.LANCZOS,
            reducing_gap=3.0,
        )
    finally:
        os.remove(input_path)
        if not picture_config.path:
            os.remove(img_path)

    return img

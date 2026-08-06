from typing import (
    List,
    Optional,
    Tuple,
    Union,
    TYPE_CHECKING,
    Dict,
    Any,
    Literal,
    cast,
    Sequence,
)
import numpy as np

from scm.plams.core.errors import MissingOptionalPackageError
from scm.plams.core.functions import requires_optional_package
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.units import Units

try:
    from scm.base import ChemicalSystem

    _has_scm_chemsys = True
except ImportError:
    _has_scm_chemsys = False

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import ase
    from os import PathLike
    from PIL import Image as PilImage
    from scm.plams.recipes.md.trajectoryanalysis import AMSMSDJob

__all__ = [
    "linear_fit_extrapolate_to_0",
    "plot_band_structure",
    "plot_phonons_band_structure",
    "plot_phonons_dos",
    "plot_phonons_thermodynamic_properties",
    "plot_energy_landscape",
    "plot_molecule",
    "plot_image_grid",
    "plot_correlation",
    "plot_msd",
    "plot_work_function",
    "plot_grid_molecules",
]


@requires_optional_package("scipy")
def linear_fit_extrapolate_to_0(x: Sequence[float], y: Sequence[float]) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Perform a linear regression on ``x`` and ``y`` and return the fit extended to ``x = 0``.

    x: sequence of float
        X values for the linear regression.

    y: sequence of float
        Y values for the linear regression.

    Returns: tuple
        ``fit_x``, ``fit_y``, ``slope``, ``intercept``.

    If ``0`` is already present in ``x``, it is not appended a second time.
    """
    from scipy.stats import linregress

    result = linregress(x, y)
    fit_x_values = list(x)
    if 0 not in fit_x_values:
        fit_x_values.append(0.0)
    fit_x = np.array(fit_x_values, dtype=float)
    fit_y = result.slope * fit_x + result.intercept

    return fit_x, fit_y, result.slope, result.intercept


@requires_optional_package("matplotlib")
def plot_band_structure(
    x: List[float],
    y_spin_up: np.ndarray,
    y_spin_down: Optional[np.ndarray] = None,
    labels: Optional[List[str]] = None,
    fermi_energy: Optional[float] = None,
    zero: Optional[Union[Literal["fermi", "vbm", "vbmax", "cbm", "cbmin"], float]] = None,
    ax: Optional["plt.Axes"] = None,
) -> "plt.Axes":
    """
    Plots an electronic band structure from DFTB, BAND, or QuantumEspresso engines with matplotlib.

    To control the appearance of the plot you need to call ``plt.ylim(bottom, top)``, ``plt.title(title)``, etc.
    manually outside this function.

    x: list of float
        Returned by AMSResults.get_band_structure()

    y_spin_up: 2D numpy array of float
        Returned by AMSResults.get_band_structure()

    y_spin_down: 2D numpy array of float. If None, the spin down bands are not plotted.
        Returned by AMSResults.get_band_structure()

    labels: list of str
        Returned by AMSResults.get_band_structure()

    fermi_energy: float
        Returned by AMSResults.get_band_structure(). Should have the same unit as ``y``.

    zero: None or float or one of 'fermi', 'vbm', 'vbmax', 'cbm', 'cbmin'
        Shift the curves so that y=0 is at the specified value. If None, no shift is performed. 'fermi', 'vbmax', and 'cbmin' require that the ``fermi_energy`` is not None. Note: 'vbmax' and 'cbmin' calculate the zero as the highest (lowest) eigenvalue smaller (greater) than or equal to ``fermi_energy``. This is NOT necessarily equal to the valence band maximum or conduction band minimum as calculated by the compute engine.

    Additional parameters:

    ``ax``: matplotlib axis
        The axis. If None, one will be created
    """
    import matplotlib.pyplot as plt

    if zero is None:
        zero_value = 0.0
    elif zero == "fermi":
        assert fermi_energy is not None
        zero_value = fermi_energy
    elif zero in ["vbm", "vbmax"]:
        assert fermi_energy is not None
        zero_value = y_spin_up[y_spin_up <= fermi_energy].max()
        if y_spin_down is not None:
            zero_value = max(zero_value, y_spin_down[y_spin_down <= fermi_energy].max())
    elif zero in ["cbm", "cbmin"]:
        assert fermi_energy is not None
        zero_value = y_spin_up[y_spin_up >= fermi_energy].min()
        if y_spin_down is not None:
            zero_value = min(zero_value, y_spin_down[y_spin_down <= fermi_energy].min())
    else:
        raise ValueError(
            f"When specified, zero must be a float or one of: 'fermi', 'vbm', 'vbmax', 'cbm', 'cbmin'; but was '{zero}'"
        )

    labels = labels or []

    for i, label in enumerate(labels):
        if label:
            label = (
                label.replace("GAMMA", "\\Gamma")
                .replace("DELTA", "\\Delta")
                .replace("LAMBDA", "\\Lambda")
                .replace("SIGMA", "\\Sigma")
            )
            labels[i] = f"${label}$"

    if ax is None:
        _, ax = plt.subplots()

    ax.plot(x, y_spin_up - zero_value, "-")
    if y_spin_down is not None:
        ax.plot(x, y_spin_down - zero_value, "--")

    tick_x: List[float] = []
    tick_labels: List[str] = []
    for xx, ll in zip(x, labels):
        if ll:
            if len(tick_x) == 0:
                tick_x.append(xx)
                tick_labels.append(ll)
                continue
            if np.isclose(xx, tick_x[-1]):
                if ll != tick_labels[-1]:
                    tick_labels[-1] += f",{ll}"
            else:
                tick_x.append(xx)
                tick_labels.append(ll)

    for xx in tick_x:
        ax.axvline(xx)

    if fermi_energy is not None:
        ax.axhline(fermi_energy - zero_value, linestyle="--")

    ax.set_xticks(ticks=tick_x, labels=tick_labels)

    return ax


@requires_optional_package("matplotlib")
def plot_phonons_band_structure(
    x: List[float],
    y: np.ndarray,
    labels: Optional[List[str]] = None,
    zero: Optional[float] = None,
    ax: Optional["plt.Axes"] = None,
) -> "plt.Axes":
    """
    Plots a phonons band structure from DFTB, BAND or QuantumEspresso engines with matplotlib.

    To control the appearance of the plot you need to call ``plt.ylim(bottom, top)``, ``plt.title(title)``, etc.
    manually outside this function.

    x: list of float
        Returned by AMSResults.get_phonons_band_structure()

    y: 2D numpy array of float
        Returned by AMSResults.get_phonons_band_structure()

    labels: list of str
        Returned by AMSResults.get_phonons_band_structure()

    zero: None or float
        Shift the curves so that y=0 is at the specified value. If None, no shift is performed.

    Additional parameters:

    ``ax``: matplotlib axis
        The axis. If None, one will be created
    """
    import matplotlib.pyplot as plt

    if zero is None:
        zero = 0

    labels = labels or []

    for i, label in enumerate(labels):
        if label:
            label = (
                label.replace("GAMMA", "\\Gamma")
                .replace("DELTA", "\\Delta")
                .replace("LAMBDA", "\\Lambda")
                .replace("SIGMA", "\\Sigma")
            )
            labels[i] = f"${label}$"

    if ax is None:
        _, ax = plt.subplots()

    ax.plot(x, y - zero, "-")

    tick_x: List[float] = []
    tick_labels: List[str] = []
    for xx, ll in zip(x, labels):
        if ll:
            if len(tick_x) == 0:
                tick_x.append(xx)
                tick_labels.append(ll)
                continue
            if np.isclose(xx, tick_x[-1]):
                if ll != tick_labels[-1]:
                    tick_labels[-1] += f",{ll}"
            else:
                tick_x.append(xx)
                tick_labels.append(ll)

    for xx in tick_x:
        ax.axvline(xx, dashes=[2, 2], color="gray")

    ax.set_xticks(ticks=tick_x, labels=tick_labels)

    return ax


@requires_optional_package("matplotlib")
def plot_phonons_dos(
    energy: List[float],
    total_dos: List[float],
    dos_per_species: Dict[str, List[float]],
    dos_per_atom: Dict[str, List[float]],
    dos_type: str = "total",
    ax: Optional["plt.Axes"] = None,
) -> "plt.Axes":
    """
    Plots the phonons DOS from DFTB, BAND or QuantumEspresso engines with matplotlib.

    To control the appearance of the plot you need to call ``plt.ylim(bottom, top)``, ``plt.title(title)``, etc.
    manually outside this function.

    energy: list of float
        Returned by AMSResults.get_phonons_thermodynamic_properties()

    total_dos: list of float
        Returned by AMSResults.get_phonons_thermodynamic_properties()

    dos_per_species: dictionary of list of float
        Returned by AMSResults.get_phonons_thermodynamic_properties()

    dos_per_atom: dictionary of list of float
        Returned by AMSResults.get_phonons_thermodynamic_properties()

    dos_type: str
        Specifies the kind of plot to show. Possible options:
            - "total": Total DOS.
            - "species": DOS decomposed by species.
            - "atom": DOS decomposed by atom.

    Additional parameters:

    ``ax``: matplotlib axis
        The axis. If None, one will be created
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    if dos_type == "total":
        ax.plot(energy, total_dos, color="black", label="Total DOS", linestyle="-", zorder=1)

    elif dos_type == "species":
        ax.plot(
            energy,
            total_dos,
            color="black",
            label="Total DOS",
            linestyle="-",
            zorder=-1,
        )
        for i, (l, v) in enumerate(dos_per_species.items()):
            ax.plot(energy, v, label=f"pDOS {l}", dashes=[3, i + 1, 2], zorder=i)

    elif dos_type == "atoms":
        ax.plot(
            energy,
            total_dos,
            color="black",
            label="Total DOS",
            linestyle="-",
            zorder=-1,
        )
        for i, (l, v) in enumerate(dos_per_atom.items()):
            ax.plot(energy, v, label=f"pDOS {l}", dashes=[3, i + 1, 2], zorder=i)

    else:
        raise ValueError("Invalid dos_type. Must be 'total', 'species', or 'atom'.")

    plt.legend()

    return ax


@requires_optional_package("matplotlib")
def plot_phonons_thermodynamic_properties(
    temperature: List[float],
    properties: Dict[str, List[float]],
    units: Dict[str, str],
    ax: Optional["plt.Axes"] = None,
) -> "plt.Axes":
    """
    Plots the phonons thermodynamic properties from DFTB, BAND or QuantumEspresso engines with matplotlib.

    To control the appearance of the plot you need to call ``plt.ylim(bottom, top)``, ``plt.title(title)``, etc.
    manually outside this function.

    temperature: list of float
        Returned by AMSResults.get_phonons_thermodynamic_properties()

    properties: dictionary of list of float
        Returned by AMSResults.get_phonons_thermodynamic_properties()

    units: dictionary of str
        Returned by AMSResults.get_phonons_thermodynamic_properties()

    Additional parameters:

    ``ax``: matplotlib axis
        The axis. If None, one will be created
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    for i, (label, prop) in enumerate(properties.items()):
        ax.plot(
            temperature,
            prop,
            label=label + " (" + units[label] + ")",
            linestyle="-",
            lw=2,
            zorder=1,
        )

    plt.legend()

    return ax


@requires_optional_package("matplotlib")
@requires_optional_package("ase")
def plot_molecule(
    molecule: Union["ase.Atoms", Molecule],
    figsize: Optional[Union[int, Tuple[int]]] = None,
    ax: Optional["plt.Axes"] = None,
    keep_axis: bool = False,
    **kwargs: Any,
) -> "plt.Axes":
    """Show a molecule in a Jupyter notebook"""
    import matplotlib.pyplot as plt
    from ase.visualize.plot import plot_atoms
    from scm.plams.interfaces.molecule.ase import toASE

    if isinstance(molecule, Molecule):
        molecule = toASE(molecule)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize or (2, 2))

    plot_atoms(molecule, ax=ax, **kwargs)

    if keep_axis:
        ax.axis("on")
    else:
        ax.axis("off")

    return ax


@requires_optional_package("rdkit")
def plot_grid_molecules(
    molecules: List[Union[Molecule, "ChemicalSystem"]],
    legends: Optional[List[str]] = None,
    molsPerRow: int = 2,
    subImgSize: Tuple[int, int] = (200, 200),
    ax: Optional["plt.Axes"] = None,
    save_svg_path: Optional[Union[str, "PathLike"]] = None,
    **kwargs: Any,
) -> Union["PilImage.Image", "plt.Axes", str]:
    """Plot series of molecules in a grid using RDKit

    :param ax: if provided molecules are plotted in these axes, note that the quality of the image might reduce, defaults to None
    :param save_svg_path: pathlike of the file, with formats .svg to save it as image, it returns the svg string, defaults to None
    :return: an image of the molecules
    :rtype: pil.Image or plt.Axes or string
    """
    from rdkit.Chem import Draw, rdchem
    from rdkit.Chem.Draw import IPythonConsole  # type: ignore[attr-defined]
    from scm.plams.interfaces.molecule.rdkit import _rdmol_for_image

    # guess bonds, the bonds will be included in the RDKit molecule
    for m in molecules:
        if len(m.bonds) == 0:
            m.guess_bonds()
    molecules = [_rdmol_for_image(m, remove_hydrogens=False) for m in molecules]

    if ax is not None or save_svg_path is not None:
        if hasattr(rdchem.Mol, "_repr_svg_"):
            IPythonConsole.UninstallIPythonRenderer()
    else:
        if not hasattr(rdchem.Mol, "_repr_svg_"):
            IPythonConsole.InstallIPythonRenderer()
        IPythonConsole.ipython_useSVG = True
    if ax is not None:
        IPythonConsole.ipython_useSVG = False
        kwargs["useSVG"] = False
    if save_svg_path is not None:
        kwargs["useSVG"] = True

    img = Draw.MolsToGridImage(
        mols=molecules,
        molsPerRow=molsPerRow,  # Number of molecules per row
        subImgSize=subImgSize,  # Size of each individual image
        legends=legends,
        **kwargs,
    )

    if save_svg_path is not None:
        if not isinstance(img, str):
            raise TypeError(
                f"{type(img)=} but expected str, most likely it is due to previously using ipy_useSVG=True in a notebook"
            )
        with open(save_svg_path, "w") as f:
            f.write(img)
        return img

    if ax is not None and save_svg_path is None:
        image_data = np.array(img, dtype=np.int32)
        ax.imshow(image_data)
        return ax
    return img


@requires_optional_package("matplotlib")
def plot_image_grid(
    images: Dict[str, "PilImage.Image"],
    rows: Optional[int] = None,
    cols: Optional[int] = None,
    figsize: Optional[Tuple[float, float]] = None,
    show_labels: bool = True,
    save_path: Optional[Union[str, "PathLike"]] = None,
) -> np.ndarray:
    """Plot a dictionary of images in a matplotlib grid.

    :param images: dictionary with labels as keys and images as values; iteration order determines image order in the grid
    :param rows: number of rows in the grid; if ``None``, infer from ``cols`` and number of images
    :param cols: number of columns in the grid; if ``None``, infer from ``rows`` and number of images
    :param figsize: matplotlib figure size; if ``None``, uses a grid-proportional default
    :param show_labels: whether to show labels above images; labels are taken from dictionary keys
    :param save_path: optional path to save the plotted grid image using matplotlib ``savefig``
    :return: 2D numpy array of matplotlib axes with shape ``(rows, cols)``
    :rtype: np.ndarray
    """
    import matplotlib.pyplot as plt

    items = list(images.items())
    n_images = len(items)

    if n_images == 0:
        raise ValueError("images must contain at least one image")

    if rows is not None and rows <= 0:
        raise ValueError(f"rows must be a positive integer when provided, but got {rows}")
    if cols is not None and cols <= 0:
        raise ValueError(f"cols must be a positive integer when provided, but got {cols}")

    if rows is None and cols is None:
        cols = int(np.ceil(np.sqrt(n_images)))
        rows = int(np.ceil(n_images / cols))
    elif rows is None:
        rows = int(np.ceil(n_images / cols))  # type: ignore[operator]
    elif cols is None:
        cols = int(np.ceil(n_images / rows))

    grid_size = rows * cols  # type: ignore[operator]
    if n_images > grid_size:
        raise ValueError(f"Grid of shape ({rows}, {cols}) can hold at most {grid_size} images, but got {n_images}")

    if figsize is None:
        figsize = ((4.0 * cols), (4.0 * rows))  # type: ignore[operator]
    fig, axes = plt.subplots(rows, cols, figsize=figsize)  # type: ignore[arg-type]
    axes = np.array(axes, dtype=object).reshape(rows, cols)  # type: ignore[arg-type]

    for ax in axes.flat:
        ax.axis("off")

    for i, (key, image) in enumerate(items):
        row, col = divmod(i, cols)  # type: ignore[operator]
        ax = cast(Any, axes[row, col])
        ax.imshow(image)
        if show_labels:
            ax.set_title(key)

    if save_path is not None:
        fig.savefig(save_path)

    return axes


def get_correlation_xy(
    job1: Union[AMSJob, List[AMSJob]],
    job2: Union[AMSJob, List[AMSJob]],
    section: str,
    variable: str,
    alt_section: Optional[str] = None,
    alt_variable: Optional[str] = None,
    file: str = "ams",
    multiplier: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    def tolist(x: Any) -> List:
        if isinstance(x, list):
            return x
        return [x]

    job1 = tolist(job1)
    job2 = tolist(job2)

    alt_section = alt_section or section
    alt_variable = alt_variable or variable

    data1 = []
    data2 = []
    for j1, j2 in zip(job1, job2):
        try:
            d1 = cast(
                Union[List[float], float],
                j1.results.readrkf(section, variable, file=file),
            )
        except KeyError:
            d1 = cast(
                Union[List[float], float],
                j1.results.get_history_property(variable, history_section=section),
            )
        d1a = np.ravel(d1) * multiplier

        try:
            d2 = cast(
                Union[List[float], float],
                j2.results.readrkf(alt_section, alt_variable, file=file),
            )
        except KeyError:
            d2 = cast(
                Union[List[float], float],
                j2.results.get_history_property(alt_variable, history_section=alt_section),
            )
        d2a = np.ravel(d2) * multiplier

        data1.extend(list(d1a))
        data2.extend(list(d2a))

    return np.array(data1), np.array(data2)


@requires_optional_package("matplotlib")
def plot_correlation(
    job1: Union[AMSJob, List[AMSJob]],
    job2: Union[AMSJob, List[AMSJob]],
    section: str,
    variable: str,
    alt_section: Optional[str] = None,
    alt_variable: Optional[str] = None,
    file: str = "ams",
    multiplier: float = 1.0,
    unit: Optional[str] = None,
    save_txt: Optional[str] = None,
    ax: Optional["plt.Axes"] = None,
    show_xy: bool = True,
    show_linear_fit: bool = True,
    show_mad: bool = True,
    show_rmsd: bool = True,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
) -> "plt.Axes":
    """

    Plot a correlation plot from AMS .rkf files

    job1: AMSJob or List[AMSJob]
        Job(s) plotted on x-axis

    job2: AMSJob or List[AMSJob]
        job2: Job(s) plotted on y-axis

    section: str
        section: section to read on .rkf files

    variable: str
        variable: variable to read

    alt_section: str
        Section to read on .rkf files for job2. If not specified it will be the same as ``section``

    alt_variable : str
        Variable to read for job2. If not specified it will be the same as ``variable``.

    file: str, optional
        file: "ams" or "engine", defaults to "ams"

    multiplier: float, optional
        multiplier: Numbers will be multiplied by this number, defaults to 1.0

    unit: str, optional
        unit: unit will be shown in the plot, defaults to None

    save_txt: str, optional
        save_txt: If not None, save the xy data to this text file, defaults to None

    ax: matplotlib axis, optional
        ax: matplotlib axis, defaults to None

    show_xy: bool, optional
        show_xy: Whether to show y=x line, defaults to True

    show_linear_fit: bool, optional
        show_linear_fit: Whether to perform and show a linear fit, defaults to True

    show_mad: bool, optional
        show_mad: Whether to show mean absolute deviation, defaults to True

    show_rmsd: bool, optional
        show_rmsd: Whether to show root-mean-square deviation, defaults to True

    xlabel: str, optional
        xlabel: The x-label. If not given will be a list of job names, defaults to None

    ylabel: str, optional
        ylabel: THe y-label. If not given will be al ist of job names, defaults to None

    Returns: A matplotlib axis

    """

    import matplotlib.pyplot as plt

    def tolist(x: Any) -> List:
        if isinstance(x, list):
            return x
        return [x]

    job1 = tolist(job1)
    job2 = tolist(job2)

    alt_section = alt_section or section
    alt_variable = alt_variable or variable

    data1, data2 = get_correlation_xy(job1, job2, section, variable, alt_section, alt_variable, file, multiplier)

    def add_unit(s: str) -> str:
        if unit is not None:
            return f"{s} ({unit})"

        return s

    if ax is None:
        fig, ax = plt.subplots()

    complete_data = np.stack((data1, data2), axis=1)

    min_data = np.min(complete_data)
    max_data = np.max(complete_data)
    min_max = np.array([min_data, max_data])

    legend = []
    title = [f"{section}%{variable}"]
    if show_xy:
        ax.plot(min_max, min_max, "-")
        legend.append("y=x")

    stats_title = ""

    if show_mad:
        mad = np.mean(np.abs(data2 - data1))
        stats_title += add_unit(f" MAD: {mad:.5f}")

    if show_rmsd:
        rmsd = np.sqrt(np.mean((data2 - data1) ** 2))
        stats_title += add_unit(f" RMSD: {rmsd:.5f}")

    linear_fit_title = None
    if show_linear_fit:
        try:
            from scipy.stats import linregress
        except ImportError:
            raise MissingOptionalPackageError("scipy")

        result = linregress(data1, data2)
        min_max_linear_fit = result.slope * min_max + result.intercept
        r2 = result.rvalue**2
        ax.plot(min_max, min_max_linear_fit, "-")
        legend.append("Fit")
        stats_title += f" R^2: {r2:.3f}"
        linear_fit_title = f"Linear fit slope={result.slope:.3f} intercept={result.intercept:.3f}"

    if stats_title:
        title.append(stats_title)

    if linear_fit_title:
        title.append(linear_fit_title)

    ax.plot(data1, data2, ".")
    legend.append("data")

    if xlabel is None:
        xlabel = ", ".join(x.name for x in job1)
        if len(xlabel) > 40:
            xlabel = xlabel[:35] + "..."
    if ylabel is None:
        ylabel = ", ".join(x.name for x in job2)
        if len(ylabel) > 40:
            ylabel = ylabel[:35] + "..."

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.set_title("\n".join(title))
    ax.legend(legend)

    ax.set_box_aspect(1)
    ax.set_xlim(*min_max)
    ax.set_ylim(*min_max)

    if save_txt is not None:
        np.savetxt(save_txt, complete_data, header=f"{xlabel} {ylabel}")

    return ax


@requires_optional_package("matplotlib")
def plot_msd(
    job: "AMSMSDJob",
    start_time_fit_fs: Optional[float] = None,
    ax: Optional["plt.Axes"] = None,
) -> "plt.Axes":
    """
    job: AMSMSDJob
        The job for which to plot the results

    start_time_fit_fs: float
        The start time (in fs) for which to perform the linear fit

    ax: matplotlib axis
        The axis. If None, one will be created

    Returns: matplotlib axis
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    time, msd = job.results.get_msd()
    fit_result, fit_x, fit_y = job.results.get_linear_fit(start_time_fit_fs=start_time_fit_fs)
    # the diffusion coefficient can also be calculated as fit_result.slope/6 (ang^2/fs)
    diffusion_coefficient = job.results.get_diffusion_coefficient(start_time_fit_fs=start_time_fit_fs)  # m^2/s
    ax.plot(time, msd, label="MSD")
    ax.plot(fit_x, fit_y, label=f"Linear fit slope={fit_result.slope:.5f} ang^2/fs")
    ax.legend()
    ax.set_xlabel("Correlation time (fs)")
    ax.set_ylabel("Mean square displacement (ang^2)")
    ax.set_title(f"MSD: Diffusion coefficient = {diffusion_coefficient:.2e} m^2/s")

    return ax


@requires_optional_package("matplotlib")
def plot_work_function(
    coordinate: np.ndarray,
    planarAverage: np.ndarray,
    macroscopicAverage: np.ndarray,
    Efermi: float,
    Vbulk: float,
    Vvacuum: Tuple[float, float],
    WF: Tuple[float, float],
    ax: Optional["plt.Axes"] = None,
) -> "plt.Axes":
    """
    Plots an Electrostatic Potential Profile from AMS-QE with matplotlib.

    To control the appearance of the plot you need to call ``plt.ylim(bottom, top)``, ``plt.title(title)``, etc.
    manually outside this function.

    ``coordinate``: 1D array of float.
        Returned by AMSResults.get_work_function_results().

    ``planarAverage``: 1D array of float.
        Returned by AMSResults.get_work_function_results().

    ``macroscopicAverage``: 1D array of float.
        Returned by AMSResults.get_work_function_results(). Should have the same unit as ``planarAverage``.

    ``Efermi``: float.
        Returned by AMSResults.get_work_function_results(). Should have the same unit as ``planarAverage``.

    ``Vbulk``: float.
        Returned by AMSResults.get_work_function_results(). Should have the same unit as ``planarAverage``.

    ``Vvacuum``: Tuple[float,float].
        Returned by AMSResults.get_work_function_results(). Should have the same unit as ``planarAverage``.

    ``WF``: Tuple[float,float].
        Returned by AMSResults.get_work_function_results(). Should have the same unit as ``planarAverage``.

    Additional parameters:

    ``ax``: matplotlib axis
        The axis. If None, one will be created

    Returns: matplotlib axis
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    ax.set_xlabel("Length", fontsize=13)
    ax.set_ylabel("Energy", fontsize=13)

    x0 = min(coordinate)
    y0 = min(planarAverage)
    x1 = max(coordinate)

    ax.plot(coordinate, planarAverage, color="red", linestyle="-.", lw=2, zorder=1)
    ax.plot(coordinate, macroscopicAverage, color="blue", linestyle="-", lw=2, zorder=2)
    ax.text(x0 + 0.8 * (x1 - x0), y0, "Planar\nAverage", fontsize=11, color="red")
    ax.text(x0 + 0.0 * (x1 - x0), y0, "Macroscopic\nAverage", fontsize=11, color="blue")
    ax.axhline(y=Efermi, color="black", linestyle="dashed", linewidth=1)
    ax.text(x0 + 0.05 * (x1 - x0), Efermi + 0.1, "E. Fermi", fontsize=11, color="black")
    ax.axhline(y=Vbulk, color="black", linestyle="dashed", linewidth=1)
    ax.text(x0 + 0.05 * (x1 - x0), Vbulk + 0.1, "Pot. bulk", fontsize=11, color="black")

    # If the material is symmetric:
    if abs(Vvacuum[0] - Vvacuum[1]) < 1e-3 or abs(Vvacuum[0] - Vvacuum[1]) < 1e-3:
        ax.axhline(y=Vvacuum[0], color="black", linestyle="dashed", linewidth=1)
        ax.text(min(coordinate), Vvacuum[0] + 0.1, "Pot. vacuum", fontsize=11, color="black")

        head_length = 0.4
        ax.arrow(
            x0 + 1.0 * (x1 - x0),
            Efermi,
            0.0,
            Vvacuum[1] - Efermi - head_length,
            head_width=0.3,
            head_length=head_length,
            fc="black",
            ec="black",
        )
        ax.text(
            x0 + 0.98 * (x1 - x0),
            (Vvacuum[1] + Efermi) / 2,
            f"WF={WF[1]:.1f} eV",
            fontsize=11,
            color="black",
            horizontalalignment="right",
        )

    # Otherwise:
    else:
        ax.plot(
            [x0, x0 + 0.3 * (x1 - x0)],
            [Vvacuum[0], Vvacuum[0]],
            color="black",
            linestyle="dashed",
            linewidth=1,
        )
        ax.text(x0, Vvacuum[0] + 0.1, "Pot. vacuum", fontsize=11, color="black")

        ax.plot(
            [x1, x1 - 0.3 * (x1 - x0)],
            [Vvacuum[1], Vvacuum[1]],
            color="black",
            linestyle="dashed",
            linewidth=1,
        )
        ax.text(
            x1 - 0.3 * (x1 - x0),
            Vvacuum[1] + 0.1,
            "Pot. vacuum",
            fontsize=11,
            color="black",
        )

        head_length = 0.4
        ax.arrow(
            x0 + 0.0 * (x1 - x0),
            Efermi,
            0.0,
            Vvacuum[0] - Efermi - head_length,
            head_width=0.3,
            head_length=head_length,
            fc="black",
            ec="black",
        )
        ax.text(
            x0 + 0.02 * (x1 - x0),
            (Vvacuum[0] + Efermi) / 2,
            f"WF={WF[0]:.1f} eV",
            fontsize=11,
            color="black",
        )
        ax.arrow(
            x0 + 1.0 * (x1 - x0),
            Efermi,
            0.0,
            Vvacuum[1] - Efermi - head_length,
            head_width=0.3,
            head_length=head_length,
            fc="black",
            ec="black",
        )
        ax.text(
            x0 + 0.98 * (x1 - x0),
            (Vvacuum[1] + Efermi) / 2,
            f"WF={WF[1]:.1f} eV",
            fontsize=11,
            color="black",
            horizontalalignment="right",
        )

    return ax


@requires_optional_package("matplotlib")
def plot_energy_landscape(
    energy_landscape,
    ax: Optional["plt.Axes"] = None,
    unit: str = "eV",
    landscape_width: float = 0.45,
    spacing: float = 1.0,
    ts_color: str = "red",
    min_color: str = "black",
    connector_color: str = "black",
    connector_linestyle: Any = (0, (4, 4)),
    label_states: bool = True,
    layout: str = "auto",
    force_iterations: int = 200,
) -> "plt.Axes":
    """Plot an energy landscape returned by ``AMSResults.get_energy_landscape()``.

    State energies are shown relative to the global minimum in the requested
    unit. Local minima are plotted in ``min_color`` and transition states in
    ``ts_color``. Dashed connectors are drawn between connected states.

    Parameters
    ----------
    energy_landscape
        Energy landscape object returned by ``AMSResults.get_energy_landscape()``.
    ax
        Matplotlib axis to draw on. If ``None``, a new figure and axis are created.
    unit
        Energy unit used for the y-axis and for converting relative energies from Hartree.
    landscape_width
        Width of the horizontal line segment drawn for each state.
    spacing
        Horizontal spacing between consecutive plotted states within a component.
    ts_color
        Color used for transition-state level segments and their labels.
    min_color
        Color used for local-minimum level segments and their labels.
    connector_color
        Color used for the dashed lines connecting related states.
    connector_linestyle
        Matplotlib linestyle specification for the state-connection lines.
    label_states
        If ``True``, annotate each state with its integer state ID.
    layout
        Layout strategy used to order the states along the horizontal axis.
        Supported values are ``"auto"``, ``"dfs"``, ``"bfs"``,
        ``"longest_path"``, ``"force"``, and ``"crossings"``.
    force_iterations
        Number of relaxation iterations used by the ``"force"`` layout.

    Returns
    -------
    matplotlib.axes.Axes
        The axis containing the plotted energy landscape.
    """
    import matplotlib.pyplot as plt
    from collections import deque
    import itertools

    states = list(energy_landscape)

    if ax is None:
        _, ax = plt.subplots()

    if len(states) == 0:
        ax.set_ylabel(f"Relative Energy ({unit})")
        ax.set_xticks([])
        return ax

    state_map = {state.id: state for state in states}
    # Build a graph representation of the landscape from the TS reactant/product links.
    adjacency = {state.id: set() for state in states}
    edge_pairs = set()
    for state in states:
        if state.reactants is not None:
            adjacency[state.id].add(state.reactants.id)
            adjacency[state.reactants.id].add(state.id)
            edge_pairs.add(tuple(sorted((state.id, state.reactants.id))))
        if state.products is not None:
            adjacency[state.id].add(state.products.id)
            adjacency[state.products.id].add(state.id)
            edge_pairs.add(tuple(sorted((state.id, state.products.id))))

    def _state_sort_key(state_id):
        state = state_map[state_id]
        return (state.isTS, state.energy, state.id)

    def _component(start_id, remaining_ids):
        component = set()
        stack = [start_id]
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            stack.extend(adjacency[node] & remaining_ids)
        return component

    def _component_start(component):
        endpoints = [node for node in component if len(adjacency[node] & component) <= 1]
        candidates = endpoints if endpoints else list(component)
        return min(candidates, key=_state_sort_key)

    def _connected_components():
        # Layout disconnected reaction networks independently before placing them side by side.
        remaining_ids = {state.id for state in states}
        components = []
        while remaining_ids:
            start_id = min(remaining_ids, key=_state_sort_key)
            component = _component(start_id, remaining_ids)
            components.append(component)
            remaining_ids -= component
        components.sort(key=lambda component: _state_sort_key(_component_start(component)))
        return components

    def _ordered_components(orderer):
        ordered_ids = []
        for component in _connected_components():
            ordered_ids.extend(orderer(component))
        return ordered_ids

    def _dfs_order(component):
        # Follow one branch deeply before backtracking, which often matches a reaction path view.
        start_id = _component_start(component)
        visited = set()
        ordered_ids = []

        def _visit(node):
            visited.add(node)
            ordered_ids.append(node)
            neighbors = sorted(adjacency[node] & component, key=_state_sort_key)
            for neighbor in neighbors:
                if neighbor not in visited:
                    _visit(neighbor)

        _visit(start_id)
        return ordered_ids

    def _bfs_order(component):
        # Expand level by level from one endpoint to keep nearby states grouped together.
        start_id = _component_start(component)
        queue = deque([start_id])
        visited = {start_id}
        ordered_ids = []
        while queue:
            node = queue.popleft()
            ordered_ids.append(node)
            neighbors = sorted(adjacency[node] & component, key=_state_sort_key)
            for neighbor in neighbors:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        return ordered_ids

    def _farthest(start_id, component):
        distances = {start_id: 0}
        parents = {start_id: None}
        queue = deque([start_id])
        while queue:
            node = queue.popleft()
            for neighbor in sorted(adjacency[node] & component, key=_state_sort_key):
                if neighbor not in distances:
                    distances[neighbor] = distances[node] + 1
                    parents[neighbor] = node
                    queue.append(neighbor)
        farthest_id = max(
            distances, key=lambda node: (distances[node], -int(state_map[node].isTS), -state_map[node].energy, -node)
        )
        return farthest_id, distances, parents

    def _longest_path_order(component):
        # Use the graph backbone first, then attach side branches around that main path.
        start_id = _component_start(component)
        end_a, _, _ = _farthest(start_id, component)
        end_b, _, parents = _farthest(end_a, component)

        backbone = []
        node = end_b
        while node is not None:
            backbone.append(node)
            node = parents[node]
        backbone.reverse()

        ordered_ids = []
        placed = set()

        def _add_side_branch(root_id, blocked):
            neighbors = sorted((adjacency[root_id] & component) - blocked - placed, key=_state_sort_key)
            for neighbor in neighbors:
                ordered_ids.append(neighbor)
                placed.add(neighbor)
                _add_side_branch(neighbor, blocked | {root_id})

        for node in backbone:
            if node not in placed:
                ordered_ids.append(node)
                placed.add(node)
            _add_side_branch(node, set(backbone))

        leftovers = sorted(component - placed, key=_state_sort_key)
        ordered_ids.extend(leftovers)
        return ordered_ids

    def _crossings_for_order(order):
        # Score a 1D ordering by counting connector crossings and total connector length.
        positions = {state_id: idx for idx, state_id in enumerate(order)}
        edge_list = sorted(edge_pairs)
        crossings = 0
        edge_length = 0
        for idx, (a1, b1) in enumerate(edge_list):
            x1, x2 = sorted((positions[a1], positions[b1]))
            edge_length += x2 - x1
            for a2, b2 in edge_list[idx + 1 :]:
                if len({a1, b1, a2, b2}) < 4:
                    continue
                y1, y2 = sorted((positions[a2], positions[b2]))
                if (x1 < y1 < x2 < y2) or (y1 < x1 < y2 < x2):
                    crossings += 1
        return crossings, edge_length, tuple(order)

    def _crossings_order(component):
        # For small graphs try all permutations; otherwise improve a good initial guess locally.
        component_list = sorted(component, key=_state_sort_key)
        if len(component_list) <= 8:
            return list(min(itertools.permutations(component_list), key=_crossings_for_order))

        best = min(
            (_dfs_order(component), _bfs_order(component), _longest_path_order(component)),
            key=_crossings_for_order,
        )
        improved = list(best)
        improved_flag = True
        while improved_flag:
            improved_flag = False
            best_score = _crossings_for_order(improved)
            for i in range(len(improved) - 1):
                candidate = improved.copy()
                candidate[i], candidate[i + 1] = candidate[i + 1], candidate[i]
                candidate_score = _crossings_for_order(candidate)
                if candidate_score < best_score:
                    improved = candidate
                    best_score = candidate_score
                    improved_flag = True
        return improved

    def _force_order(component):
        # Relax 1D positions with attractive edges and repulsive nodes, then sort by the relaxed positions.
        component_ids = sorted(component, key=_state_sort_key)
        initial = _longest_path_order(component)
        x_pos = {state_id: float(index) for index, state_id in enumerate(initial)}
        ideal_gap = max(spacing, 1.0)
        for _ in range(max(force_iterations, 1)):
            delta = {state_id: 0.0 for state_id in component_ids}
            for i, left in enumerate(component_ids):
                for right in component_ids[i + 1 :]:
                    distance = x_pos[right] - x_pos[left]
                    if abs(distance) < 1e-8:
                        distance = 1e-8
                    repulsion = 0.02 / abs(distance)
                    delta[left] -= repulsion
                    delta[right] += repulsion
            for left, right in edge_pairs:
                if left in component and right in component:
                    distance = x_pos[right] - x_pos[left]
                    attraction = 0.08 * (distance - ideal_gap)
                    delta[left] += attraction
                    delta[right] -= attraction
            for state_id in component_ids:
                state = state_map[state_id]
                delta[state_id] += 0.03 * Units.convert(state.energy - min(s.energy for s in states), "hartree", unit)
            for state_id in component_ids:
                x_pos[state_id] += delta[state_id]
            centered = np.mean([x_pos[state_id] for state_id in component_ids])
            for state_id in component_ids:
                x_pos[state_id] -= centered
        return sorted(component_ids, key=lambda state_id: (x_pos[state_id],) + _state_sort_key(state_id))

    layout_orderers = {
        "dfs": lambda: _ordered_components(_dfs_order),
        "bfs": lambda: _ordered_components(_bfs_order),
        "longest_path": lambda: _ordered_components(_longest_path_order),
        "force": lambda: _ordered_components(_force_order),
        "crossings": lambda: _ordered_components(_crossings_order),
    }

    if layout == "auto":
        # Compare all available strategies and keep the one with the cleanest connector pattern.
        candidate_orders = {name: orderer() for name, orderer in layout_orderers.items()}
        ordered_ids = min(
            candidate_orders.values(),
            key=lambda order: _crossings_for_order(order)
            + (sum(abs(i - order.index(state.id)) for i, state in enumerate(states)),),
        )
    else:
        if layout not in layout_orderers:
            raise ValueError(
                "Unsupported layout '{}'. Choose from 'auto', 'dfs', 'bfs', 'longest_path', 'force', or 'crossings'.".format(
                    layout
                )
            )
        ordered_ids = layout_orderers[layout]()

    ordered_states = [state_map[state_id] for state_id in ordered_ids]
    # Convert the chosen ordering into actual x-coordinates, leaving a gap between components.
    component_gap = max(1.5 * spacing, spacing + landscape_width)
    x_map = {}
    current_x = 0.0
    for component in _connected_components():
        component_order = [state_id for state_id in ordered_ids if state_id in component]
        for offset, state_id in enumerate(component_order):
            x_map[state_id] = current_x + offset * spacing
        current_x = x_map[component_order[-1]] + component_gap

    reference_energy = min(state.energy for state in states)
    relative_energies = {state.id: Units.convert(state.energy - reference_energy, "hartree", unit) for state in states}
    half_width = landscape_width / 2.0

    for state in ordered_states:
        x_pos = x_map[state.id]
        y_pos = relative_energies[state.id]
        color = ts_color if state.isTS else min_color
        ax.hlines(y=y_pos, xmin=x_pos - half_width, xmax=x_pos + half_width, colors=color, linewidth=1.6, zorder=3)
        if label_states:
            ax.text(x_pos, y_pos, str(state.id), ha="center", va="bottom", color=color, fontsize=11, zorder=4)

    plotted_pairs = set()
    # Draw each connection only once, even though TS links are visible from both endpoints.
    for state in ordered_states:
        x_pos = x_map[state.id]
        y_pos = relative_energies[state.id]
        for other in (state.reactants, state.products):
            if other is None:
                continue
            pair = tuple(sorted((state.id, other.id)))
            if pair in plotted_pairs:
                continue
            other_x = x_map[other.id]
            other_y = relative_energies[other.id]
            if other_x < x_pos:
                x1, x2 = other_x + half_width, x_pos - half_width
                y1, y2 = other_y, y_pos
            else:
                x1, x2 = x_pos + half_width, other_x - half_width
                y1, y2 = y_pos, other_y
            ax.plot([x1, x2], [y1, y2], color=connector_color, linestyle=connector_linestyle, linewidth=1.2, zorder=2)
            plotted_pairs.add(pair)

    ax.set_ylabel(f"Relative Energy ({unit})")
    ax.set_xticks([])
    x_values = [x_map[state.id] for state in ordered_states]
    ax.set_xlim(min(x_values) - half_width - 0.2, max(x_values) + half_width + 0.2)

    y_min = min(relative_energies.values())
    y_max = max(relative_energies.values())
    padding = 0.05 * max(y_max - y_min, 1.0)
    ax.set_ylim(y_min - padding, y_max + padding)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)

    return ax

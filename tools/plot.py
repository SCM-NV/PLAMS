from typing import List, Optional, Tuple, Union, TYPE_CHECKING

import numpy as np
from scm.plams.core.errors import MissingOptionalPackageError
from scm.plams.core.functions import requires_optional_package
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.mol.molecule import Molecule

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    from os import PathLike
    from PIL import Image as PilImage

__all__ = [
    "plot_band_structure",
    "plot_phonons_band_structure",
    "plot_phonons_dos",
    "plot_phonons_thermodynamic_properties",
    "plot_molecule",
    "plot_correlation",
    "plot_msd",
    "plot_work_function",
    "plot_grid_molecules",
]


@requires_optional_package("matplotlib")
def plot_band_structure(x, y_spin_up, y_spin_down=None, labels=None, fermi_energy=None, zero=None, ax=None):
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

    zero: None or float or one of 'fermi', 'vbmax', 'cbmin'
        Shift the curves so that y=0 is at the specified value. If None, no shift is performed. 'fermi', 'vbmax', and 'cbmin' require that the ``fermi_energy`` is not None. Note: 'vbmax' and 'cbmin' calculate the zero as the highest (lowest) eigenvalue smaller (greater) than or equal to ``fermi_energy``. This is NOT necessarily equal to the valence band maximum or conduction band minimum as calculated by the compute engine.

    Additional parameters:

    ``ax``: matplotlib axis
        The axis. If None, one will be created
    """
    import matplotlib.pyplot as plt

    if zero is None:
        zero = 0
    elif zero == "fermi":
        assert fermi_energy is not None
        zero = fermi_energy
    elif zero in ["vbm", "vbmax"]:
        assert fermi_energy is not None
        zero = y_spin_up[y_spin_up <= fermi_energy].max()
        if y_spin_down is not None:
            zero = max(zero, y_spin_down[y_spin_down <= fermi_energy].max())
    elif zero in ["cbm", "cbmax"]:
        assert fermi_energy is not None
        zero = y_spin_up[y_spin_up >= fermi_energy].min()
        if y_spin_down is not None:
            zero = min(zero, y_spin_down[y_spin_down <= fermi_energy].min())

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

    ax.plot(x, y_spin_up - zero, "-")
    if y_spin_down is not None:
        ax.plot(x, y_spin_down - zero, "--")

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
        ax.axhline(fermi_energy - zero, linestyle="--")

    ax.set_xticks(ticks=tick_x, labels=tick_labels)

    return ax


@requires_optional_package("matplotlib")
def plot_phonons_band_structure(x, y, labels=None, zero=None, ax=None):
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
def plot_phonons_dos(energy, total_dos, dos_per_species, dos_per_atom, dos_type="total", ax=None):
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
        ax.plot(energy, total_dos, color="black", label="Total DOS", linestyle="-", zorder=-1)
        for i, (l, v) in enumerate(dos_per_species.items()):
            ax.plot(energy, v, label=f"pDOS {l}", dashes=[3, i + 1, 2], zorder=i)

    elif dos_type == "atoms":
        ax.plot(energy, total_dos, color="black", label="Total DOS", linestyle="-", zorder=-1)
        for i, (l, v) in enumerate(dos_per_atom.items()):
            ax.plot(energy, v, label=f"pDOS {l}", dashes=[3, i + 1, 2], zorder=i)

    else:
        raise ValueError("Invalid dos_type. Must be 'total', 'species', or 'atom'.")

    plt.legend()

    return ax


@requires_optional_package("matplotlib")
def plot_phonons_thermodynamic_properties(temperature, properties, units, ax=None):
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
        ax.plot(temperature, prop, label=label + " (" + units[label] + ")", linestyle="-", lw=2, zorder=1)

    plt.legend()

    return ax


@requires_optional_package("matplotlib")
@requires_optional_package("ase")
def plot_molecule(molecule, figsize=None, ax=None, keep_axis: bool = False, **kwargs):
    """Show a molecule in a Jupyter notebook"""
    import matplotlib.pyplot as plt
    from ase.visualize.plot import plot_atoms
    from scm.plams.interfaces.molecule.ase import toASE

    if isinstance(molecule, Molecule):
        molecule = toASE(molecule)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize or (2, 2))

    plot_atoms(molecule, ax=ax, **kwargs)

    if not keep_axis:
        ax.axis("off")

    return ax


@requires_optional_package("rdkit")
def plot_grid_molecules(
    molecules: List[Molecule],
    legends: Optional[List[str]] = None,
    molsPerRow: int = 2,
    subImgSize: Tuple[int, int] = (200, 200),
    ax: Optional["plt.Axes"] = None,
    save_svg_path: Optional[Union[str, "PathLike"]] = None,
    **kwargs,
) -> Union["PilImage", "plt.Axes", str]:
    """Plot series of molecules in a grid using RDKit

    :param ax: if provided molecules are plotted in these axes, note that the quality of the image might reduce, defaults to None
    :param save_svg_path: pathlike of the file, with formats .svg to save it as image, it returns the svg string, defaults to None
    :return: an image of the molecules
    :rtype: pil.Image or plt.Axes or string
    """
    from rdkit.Chem import Draw, rdchem
    from rdkit.Chem.Draw import IPythonConsole
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


def get_correlation_xy(
    job1: Union[AMSJob, List[AMSJob]],
    job2: Union[AMSJob, List[AMSJob]],
    section: str,
    variable: str,
    alt_section: Optional[str] = None,
    alt_variable: Optional[str] = None,
    file: str = "ams",
    multiplier: float = 1.0,
) -> Tuple:
    def tolist(x):
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
            d1 = j1.results.readrkf(section, variable, file=file)
        except KeyError:
            d1 = j1.results.get_history_property(variable, history_section=section)
        d1 = np.ravel(d1) * multiplier

        try:
            d2 = j2.results.readrkf(alt_section, alt_variable, file=file)
        except KeyError:
            d2 = j2.results.get_history_property(alt_variable, history_section=alt_section)
        d2 = np.ravel(d2) * multiplier

        data1.extend(list(d1))
        data2.extend(list(d2))

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
    ax=None,
    show_xy: bool = True,
    show_linear_fit: bool = True,
    show_mad: bool = True,
    show_rmsd: bool = True,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
):
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

    def tolist(x):
        if isinstance(x, list):
            return x
        return [x]

    job1 = tolist(job1)
    job2 = tolist(job2)

    alt_section = alt_section or section
    alt_variable = alt_variable or variable

    data1, data2 = get_correlation_xy(job1, job2, section, variable, alt_section, alt_variable, file, multiplier)

    def add_unit(s: str):
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
def plot_msd(job, start_time_fit_fs=None, ax=None):
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
    ax.plot(fit_x, fit_y, label="Linear fit slope={:.5f} ang^2/fs".format(fit_result.slope))
    ax.legend()
    ax.set_xlabel("Correlation time (fs)")
    ax.set_ylabel("Mean square displacement (ang^2)")
    ax.set_title("MSD: Diffusion coefficient = {:.2e} m^2/s".format(diffusion_coefficient))

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
    ax=None,
):
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
            "WF=" + "%.1f" % WF[1] + " eV",
            fontsize=11,
            color="black",
            horizontalalignment="right",
        )

    # Otherwise:
    else:
        ax.plot([x0, x0 + 0.3 * (x1 - x0)], [Vvacuum[0], Vvacuum[0]], color="black", linestyle="dashed", linewidth=1)
        ax.text(x0, Vvacuum[0] + 0.1, "Pot. vacuum", fontsize=11, color="black")

        ax.plot([x1, x1 - 0.3 * (x1 - x0)], [Vvacuum[1], Vvacuum[1]], color="black", linestyle="dashed", linewidth=1)
        ax.text(x1 - 0.3 * (x1 - x0), Vvacuum[1] + 0.1, "Pot. vacuum", fontsize=11, color="black")

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
            x0 + 0.02 * (x1 - x0), (Vvacuum[0] + Efermi) / 2, "WF=" + "%.1f" % WF[0] + " eV", fontsize=11, color="black"
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
            "WF=" + "%.1f" % WF[1] + " eV",
            fontsize=11,
            color="black",
            horizontalalignment="right",
        )

    return ax

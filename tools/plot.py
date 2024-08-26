from scm.plams.mol.molecule import Molecule
from scm.plams.interfaces.adfsuite.ams import AMSJob
from typing import Tuple, Union, List

__all__ = ["plot_band_structure", "plot_molecule", "plot_correlation"]


def plot_band_structure(x, y_spin_up, y_spin_down=None, labels=None, fermi_energy=None, zero=None, show=False):
    """
    Plots an electronic band structure from DFTB or BAND with matplotlib.

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

    show: bool
        If True, call plt.show() at the end
    """
    import matplotlib.pyplot as plt
    import numpy as np

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

    fig, ax = plt.subplots()

    plt.plot(x, y_spin_up - zero, "-")
    if y_spin_down is not None:
        plt.plot(x, y_spin_down - zero, "--")

    tick_x = []
    tick_labels = []
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
        plt.axvline(xx)

    if fermi_energy is not None:
        plt.axhline(fermi_energy - zero, linestyle="--")

    plt.xticks(ticks=tick_x, labels=tick_labels)

    if show:
        plt.show()


def plot_molecule(molecule, figsize=None, ax=None, keep_axis: bool = False, **kwargs):
    """Show a molecule in a Jupyter notebook"""
    import matplotlib.pyplot as plt
    from ase.visualize.plot import plot_atoms
    from scm.plams.interfaces.molecule.ase import toASE

    if isinstance(molecule, Molecule):
        molecule = toASE(molecule)

    if not ax:
        plt.figure(figsize=figsize or (2, 2))

    plot_atoms(molecule, ax=ax, **kwargs)

    if not keep_axis:
        if ax:
            ax.axis("off")
        else:
            plt.axis("off")


def get_correlation_xy(
    job1: Union[AMSJob, List[AMSJob]],
    job2: Union[AMSJob, List[AMSJob]],
    section: str,
    variable: str,
    alt_section: str = None,
    alt_variable: str = None,
    file: str = "ams",
    multiplier: float = 1.0,
) -> Tuple:
    import numpy as np

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

    data1 = np.array(data1)
    data2 = np.array(data2)

    return data1, data2


def plot_correlation(
    job1: Union[AMSJob, List[AMSJob]],
    job2: Union[AMSJob, List[AMSJob]],
    section: str,
    variable: str,
    alt_section: str = None,
    alt_variable: str = None,
    file: str = "ams",
    multiplier: float = 1.0,
    unit: str = None,
    save_txt: bool = None,
    ax=None,
    show_xy: bool = True,
    show_linear_fit: bool = True,
    show_mad: bool = True,
    show_rmsd: bool = True,
    xlabel: str = None,
    ylabel: str = None,
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
    import numpy as np

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
        from scipy.stats import linregress

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

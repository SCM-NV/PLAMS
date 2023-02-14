from ..mol.molecule import Molecule

__all__ = ['plot_band_structure', 'plot_molecule']

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
    elif zero == 'fermi':
        assert fermi_energy is not None
        zero = fermi_energy 
    elif zero in ['vbm', 'vbmax']:
        assert fermi_energy is not None
        zero = y_spin_up[y_spin_up <= fermi_energy].max()
        if y_spin_down is not None:
            zero = max(zero, y_spin_down[y_spin_down <= fermi_energy].max())
    elif zero in ['cbm', 'cbmax']:
        assert fermi_energy is not None
        zero = y_spin_up[y_spin_up >= fermi_energy].min()
        if y_spin_down is not None:
            zero = min(zero, y_spin_down[y_spin_down <= fermi_energy].min())

    labels = labels or []

    fig, ax = plt.subplots()

    plt.plot(x, y_spin_up-zero, '-')
    if y_spin_down is not None:
        plt.plot(x, y_spin_down-zero, '--')

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
                    tick_labels[-1] += f',{ll}'
            else:
                tick_x.append(xx)
                tick_labels.append(ll)

    for xx in tick_x:
        plt.axvline(xx)

    if fermi_energy is not None:
        plt.axhline(fermi_energy-zero, linestyle='--')

    plt.xticks(ticks=tick_x, labels=tick_labels)

    if show:
        plt.show()



def plot_molecule(molecule, figsize=None, ax=None, **kwargs):
    """ Show a molecule in a Jupyter notebook """
    from ase.visualize.plot import plot_atoms
    import matplotlib.pyplot as plt
    from ..interfaces.molecule.ase import toASE

    if isinstance(molecule, Molecule):
        molecule = toASE(molecule)

    if not ax:
        plt.figure(figsize=figsize or (2,2))

    plot_atoms(molecule, ax=ax, **kwargs)

    if ax:
        ax.axis('off')
    else:
        plt.axis('off')



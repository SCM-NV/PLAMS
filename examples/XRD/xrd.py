#!/usr/bin/env amspython
# coding: utf-8

# ## Initial imports

from scm.plams import *

try:
    from ase import Atoms
    from pymatgen.core.structure import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
except ImportError as e:
    print(
        "You need ASE and pymatgen installed in the AMS python environment to run this example. Install the package for m3gnet to do this."
    )
    print(e)
    exit(1)


# ## Create ASE atoms object for BaTiO3

at = Atoms(
    symbols=[
        "Ba",
        "Ti",
        "O",
        "O",
        "O",
    ],
    scaled_positions=[
        [
            0.0,
            0.0,
            0.0,
        ],
        [0.5, 0.5, 0.5],
        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0],
    ],
    cell=[4.01, 4.01, 4.01],
    pbc=(True, True, True),
)
plot_molecule(at, rotation="-5x,5y,0z")
# show in Jupyter notebook


# ## Save ASE Atoms to .cif format

fname = "batio3.cif"
at.write(fname)


# ## Load .cif in pymatgen and calculate XRD

# Available radiation sources:

print(f"Available radiation sources: {XRDCalculator.AVAILABLE_RADIATION}")


# Let's choose Cu K-alpha (default):

structure = Structure.from_file(fname)
xrd_calc = XRDCalculator(wavelength="CuKa")
xrd_calc.show_plot(structure)


pattern = xrd_calc.get_pattern(structure)
print("2*Theta Intensity hkl d_hkl(angstrom)")
for two_theta, intensity, hkls, d_hkl in zip(pattern.x, pattern.y, pattern.hkls, pattern.d_hkls):
    hkl_tuples = [hkl["hkl"] for hkl in hkls]
    for hkl in hkl_tuples:
        label = ", ".join(map(str, hkl))
        print(f"{two_theta:.2f} {intensity:.2f} {hkl} {d_hkl:.3f}")

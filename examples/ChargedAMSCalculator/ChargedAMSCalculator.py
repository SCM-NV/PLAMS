#!/usr/bin/env amspython
# coding: utf-8

# ## Initial imports

from scm.plams import *
from scm.plams.interfaces.adfsuite.ase_calculator import AMSCalculator
from ase import Atoms
from ase.visualize.plot import plot_atoms

# this line is not required in AMS2025+
init()


# ## Example 1: Total system charge
#
# ### Create the charged molecule (ion)
# Create a charged ion using using `ase.Atoms` and setting the `info` dictionary.

atoms = Atoms("OH", positions=[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
# define a total charge
atoms.info["charge"] = -1

plot_atoms(atoms)


# ### Set the AMS settings
#
# First, set the AMS settings as you normally would do in PLAMS:

settings = Settings()
settings.input.ADF  # Use ADF with the default settings
settings.input.ams.Task = "SinglePoint"


# ### Run AMS through the ASE Calculator
#
# Below, the ``amsworker=False`` (default) will cause AMS to run in standalone mode. This means that all input and output files will be stored on disk.

atoms.calc = AMSCalculator(settings=settings, name="total_charge", amsworker=False)

energy = atoms.get_potential_energy()  # calculate the energy of a charged ion
print(f"Energy: {energy:.3f} eV")  # ASE uses eV as energy unit


# ### Access the input file
#
# ``atoms.calc.amsresults`` contains the corresponding PLAMS AMSResults object.
#
# ``atoms.calc.amsresults.job`` contains the corresponding PLAMS AMSJob object. This object has, for example, the ``get_input()`` method to access the input to AMS.
#
# **Note**: These are actually properties of the Calculator, not the Atoms! So if you run more calculations with the same calculator you will **overwrite** the AMSResults in ``atoms.calc.amsresults``!
#
# AMS used the following input:

print(atoms.calc.amsresults.job.get_input())


# ### Access the binary .rkf results files and use PLAMS AMSResults methods
#
# Access the paths to the binary results files:

ams_rkf = atoms.calc.amsresults.rkfpath(file="ams")
print(ams_rkf)


# If you prefer, you can use the PLAMS methods to access results like the energy:

energy2 = atoms.calc.amsresults.get_energy(unit="eV")
print(f"Energy: {energy2:.3f} eV")


# ## Example 2: Define atomic charges
#
# ### Construct a charged ion with atomic charges

atoms = Atoms("OH", positions=[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]], charges=[-1, 0])

plot_atoms(atoms)


# ### Run AMS

calc = AMSCalculator(settings=settings, name="atomic_charges")
atoms.calc = calc

atoms.get_potential_energy()  # calculate the energy of a charged ion


# AMS only considers the total charge of the system and not the individual atomic charges. PLAMS thus reuses the results of the previous calculation since the calculation is for the same chemical system. Both input options are allowed. If both input options are used, the total charge is the sum of both.

print(calc.amsresults.job.get_input())


# ## Example 3: Set the charge in the AMS System block
#
# ### Set the charge in the AMS System block
# A charge can be set for the calculator in the settings object.

atoms = Atoms("OH", positions=[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

settings = Settings()
settings.input.ADF  # Use ADF with the default settings
settings.input.ams.Task = "SinglePoint"
settings.input.ams.System.Charge = -1

calc = AMSCalculator(settings=settings, name="default_charge")
atoms.calc = calc
atoms.get_potential_energy()  # calculate the energy of a charged ion
print(calc.amsresults.job.get_input())


# In this case, the charge of the `Atoms` object is no longer used.

atoms = Atoms(
    "OH",
    positions=[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
)
atoms.info["charge"] = 100

settings = Settings()
settings.input.ADF  # Use ADF with the default settings
settings.input.ams.Task = "SinglePoint"
settings.input.ams.System.Charge = -1

calc = AMSCalculator(settings=settings, name="default_charge_overridden")
atoms.calc = calc
atoms.get_potential_energy()  # calculate the energy of a charged ion
print(calc.amsresults.job.get_input())

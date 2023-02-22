#!/usr/bin/env amspython
# coding: utf-8

# ## Initial imports

from scm.plams import *
from ase import Atoms
from ase.visualize.plot import plot_atoms

# Before running AMS jobs, you need to call init()
init()


# ## Construct a charged ion
# Create a charged ion using using `ase.Atoms` and setting the `info` dictionairy.

atoms = Atoms('OH',
              positions = [[1.0,0.0,0.0],[0.0,0.0,0.0]]
             )
#define a total charge
atoms.info['charge'] = -1

plot_atoms(atoms, scale=0.5);


# ## Set the AMS settings
# 
# First, set the AMS settings as you normally would do:

settings = Settings()
settings.input.ADF #Use ADF with the default settings
settings.input.ams.Task = "SinglePoint"


# ## Run AMS

calc = AMSCalculator(settings = settings, name='total_charge')
atoms.calc = calc

atoms.get_potential_energy() #calculate the energy of a charged ion


# AMS used the following input file:

print(calc.amsresults.job.get_input())


# ## Construct a charged ion with atomic charges

atoms = Atoms('OH',
              positions = [[1.0,0.0,0.0],[0.0,0.0,0.0]],
              charges = [-1, 0]
             )

plot_atoms(atoms, scale=0.5);


# ## Run AMS 

calc = AMSCalculator(settings = settings, name='atomic_charges')
atoms.calc = calc

atoms.get_potential_energy() #calculate the energy of a charged ion


# AMS only considers the total charge of the system and not the individual atomic charges. PLAMS thus reuses the results of the previous calculation since the calculation is for the same chemical system. Both input options are allowed. If both input options are used, the total charge is the sum of both.

print(calc.amsresults.job.get_input())


# ## Setting the charge as a calculator property
# A charge can be set for the calculator in the settings object. 

atoms = Atoms('OH',
              positions = [[1.0,0.0,0.0],[0.0,0.0,0.0]]
             )

settings = Settings()
settings.input.ADF #Use ADF with the default settings
settings.input.ams.Task = "SinglePoint"
settings.input.ams.System.Charge = -1

calc = AMSCalculator(settings = settings, name='default_charge')
atoms.calc = calc
atoms.get_potential_energy() #calculate the energy of a charged ion
print(calc.amsresults.job.get_input())


# In this case, the charge of the `Atoms` object is no longer used.

atoms = Atoms('OH',
              positions = [[1.0,0.0,0.0],[0.0,0.0,0.0]],
             )
atoms.info['charge'] = 100

settings = Settings()
settings.input.ADF #Use ADF with the default settings
settings.input.ams.Task = "SinglePoint"
settings.input.ams.System.Charge = -1

calc = AMSCalculator(settings = settings, name='default_charge_overridden')
atoms.calc = calc
atoms.get_potential_energy() #calculate the energy of a charged ion
print(calc.amsresults.job.get_input())


# ## Finish PLAMS

finish()


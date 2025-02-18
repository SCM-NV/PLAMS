#!/usr/bin/env amspython
# coding: utf-8

# ## Initial imports

import scm.plams as plams
from scm.plams.tools.plot import plot_correlation, get_correlation_xy
import matplotlib.pyplot as plt

# this line is not required in AMS2025+
plams.init()


# ## Define two engines to compare
#
# Here we choose GFNFF and GFN1-xTB

e1 = plams.Settings()
e1.input.GFNFF

e2 = plams.Settings()
e2.input.DFTB.Model = "GFN1-xTB"


# Let's use a glycine molecule generated from SMILES:

glycine = plams.from_smiles("C(C(=O)O)N")
plams.plot_molecule(glycine)


# Run a single-point calculation storing the Gradients (negative forces):

sp = plams.Settings()
sp.input.ams.Task = "SinglePoint"
sp.input.ams.Properties.Gradients = "Yes"
sp.runscript.nproc = 1  # run in serial


job1 = plams.AMSJob(settings=sp + e1, name="glycine-engine1", molecule=glycine)
job2 = plams.AMSJob(settings=sp + e2, name="glycine-engine2", molecule=glycine)


job1.run()
job2.run()


plot_correlation(
    job1,
    job2,
    section="AMSResults",
    variable="Gradients",
    file="engine",
)


# To get the actual numbers, use ``get_correlation_xy``:

x, y = plams.tools.plot.get_correlation_xy(job1, job2, section="AMSResults", variable="Gradients", file="engine")
print("x")
print(x)
print("y")
print(y)


# ## Compare multiple jobs

smiles_list = ["CC=C", "CCCO", "C(C(=O)O)N"]
names = ["propene", "propanol", "glycine"]
molecules = [plams.from_smiles(x) for x in smiles_list]
for mol in molecules:
    plams.plot_molecule(mol)


jobs1 = [plams.AMSJob(settings=sp + e1, name="e1" + name, molecule=mol) for name, mol in zip(names, molecules)]
jobs2 = [plams.AMSJob(settings=sp + e2, name="e2" + name, molecule=mol) for name, mol in zip(names, molecules)]


for job in jobs1 + jobs2:
    job.run()


# The correlation plot can be plotted as before. You can also add a unit conversion to get your preferred units, and add custom xlabel and ylabel:

unit = "eV/angstrom"
multiplier = plams.Units.convert(1.0, "hartree/bohr", unit)

plot_correlation(
    jobs1,
    jobs2,
    section="AMSResults",
    variable="Gradients",
    file="engine",
    xlabel="Engine 1",
    ylabel="Engine 2",
    unit=unit,
    multiplier=multiplier,
)


plot_correlation(
    jobs1,
    jobs2,
    section="AMSResults",
    variable="Charges",
    file="engine",
    xlabel="Engine 1",
    ylabel="Engine 2",
)


# ## Use Task Replay to compare multiple frames from a trajectory

# The forces from an MD job can be stored with ``writeenginegradients=True``

md = plams.AMSNVEJob(
    settings=e1,
    name="nve-md-e1",
    molecule=glycine,
    velocities=400,
    nsteps=100,
    samplingfreq=10,
    writeenginegradients=True,
)
md.run()


#

# When using the Replay task, set ``Properties.Gradients`` to get the forces:

replay_s = plams.Settings()
replay_s.input.ams.Task = "Replay"
replay_s.input.ams.Properties.Gradients = "Yes"
replay_s.input.ams.Replay.File = md.results.rkfpath()
replay = plams.AMSJob(settings=e2 + replay_s, name="replay-e2")
replay.run()


# For the MD job the gradients (negative forces) are stored in ``History%EngineGradients``, but for the Replay job they are stored in ``History%Gradients``. Use the ``alt_variable`` to specify the variable for the second job:

plot_correlation(md, replay, section="History", variable="EngineGradients", alt_variable="Gradients", file="ams")

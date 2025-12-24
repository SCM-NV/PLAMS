Worked Example
--------------

Initial imports
~~~~~~~~~~~~~~~

.. code:: ipython3

   import scm.plams as plams
   import os
   import numpy as np

   try:
       from scm.plams import view  # view molecule using AMSview in a Jupyter Notebook in AMS2026+
   except ImportError:
       from scm.plams import plot_molecule  # plot molecule in a Jupyter Notebook in AMS2023+

       def view(molecule, **kwargs):
           plot_molecule(molecule)

Create initial benzene system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   benzene_box = plams.packmol(plams.from_smiles("c1ccccc1", forcefield="uff"), density=0.875, n_molecules=40)
   view(benzene_box, width=300, height=300, direction="tilt_z", padding=-2)

.. figure:: ViscosityGreenKubo_files/ViscosityGreenKubo_3_0.png

Equilibration MD simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Short equilibration MD simulation at 300 K for 5 ps with the GAFF force field and automatic atom typing.

.. code:: ipython3

   s = plams.Settings()
   s.input.ForceField.Type = "GAFF"
   s.input.ForceField.AnteChamberIntegration = "Yes"
   eq_job = plams.AMSNVTJob(
       settings=s,
       molecule=benzene_box,
       name="equilibration",
       nsteps=5000,
       timestep=1.0,
       temperature=300,
       thermostat="Berendsen",
   )
   eq_job.run();

::

   [24.12|16:10:09] JOB equilibration STARTED
   [24.12|16:10:09] JOB equilibration RUNNING
   [24.12|16:10:51] JOB equilibration FINISHED
   [24.12|16:10:51] JOB equilibration SUCCESSFUL

Production MD simulation
~~~~~~~~~~~~~~~~~~~~~~~~

Enable the “Bin log” for Pressure tensor and Time. This stores these quantities at every time step of the simulation, even if the structure is not stored. It is required to be able to calculate viscosity from the Green-Kubo relation.

Here, we run a very short MD simulation of 200 ps. Typically, to converge the viscosity you would need longer simulations.

The time step is quite large at 1 fs. You may want to decrease it (and increase the number of MD steps correspondingly).

.. code:: ipython3

   prod_job = plams.AMSNVTJob.restart_from(
       eq_job,
       name="production",
       nsteps=200000,
       timestep=1,
       temperature=300,
       thermostat="NHC",
       tau=100,
       binlog_pressuretensor=True,
       binlog_time=True,
   )

.. code:: ipython3

   prod_job.run()

::

   [24.12|16:10:51] JOB production STARTED
   [24.12|16:10:51] JOB production RUNNING
   [24.12|16:37:59] JOB production FINISHED
   [24.12|16:38:00] JOB production SUCCESSFUL





   <scm.plams.interfaces.adfsuite.ams.AMSResults at 0x328d62730>

Calculate the viscosity autocorrelation integral
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   from scm.plams.recipes.md.trajectoryanalysis import AMSViscosityFromBinLogJob

   job = AMSViscosityFromBinLogJob(previous_job=prod_job)
   job.run();

::

   [24.12|16:38:00] JOB plamsjob STARTED
   [24.12|16:38:00] JOB plamsjob RUNNING
   [24.12|16:38:24] JOB plamsjob FINISHED


   /path/plams/recipes/md/trajectoryanalysis.py:275: RuntimeWarning: overflow encountered in exp
     return A * (lam * (1 - np.exp(-x / tau1)) + (1 - lam) * (1 - np.exp(-x / tau2)))


   [24.12|16:38:25] JOB plamsjob SUCCESSFUL

Plot the results
~~~~~~~~~~~~~~~~

The viscosity is the limiting value of the viscosity integral.

The viscosity integral **must** converge to a near-constant value. If it is not close to converging, then you should run a longer MD simulation.

To better estimate the viscosity, a double-exponential function is fitted to the viscosity integral curve. This makes it easier to get a numerical value for the viscosity.

:math:`\eta(t) = A(\lambda (1 - \exp(-t/\tau _1) + (1-\lambda) (1- \exp(-t/\tau _2))`

The parameters :math:`A`, :math:`\lambda`, :math:`\tau _1`, and :math:`\tau _2` are fitted and available as ``popt[0]``, ``popt[1]``, ``popt[2]``, and ``popt[3]``.

.. code:: ipython3

   x, viscosity_integral = job.results.get_viscosity_integral()
   popt, x_fit, fit_viscosity_integral = job.results.get_double_exponential_fit()

::

   /path/plams/recipes/md/trajectoryanalysis.py:275: RuntimeWarning: overflow encountered in exp
     return A * (lam * (1 - np.exp(-x / tau1)) + (1 - lam) * (1 - np.exp(-x / tau2)))

.. code:: ipython3

   import matplotlib.pyplot as plt

   fig, ax = plt.subplots()
   ax.set_title(f"Limiting value: {popt[0]:.5f} Pa s")
   ax.plot(x, viscosity_integral, color="tab:blue")
   ax.plot(x_fit, fit_viscosity_integral, color="tab:red")
   ax.set_xlabel("Correlation Time (fs)")
   ax.set_ylabel("Viscosity (Pa s)");

.. figure:: ViscosityGreenKubo_files/ViscosityGreenKubo_13_0.png

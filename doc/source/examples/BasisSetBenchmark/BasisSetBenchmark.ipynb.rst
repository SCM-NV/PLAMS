Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

   import multiprocessing
   from scm.plams import JobRunner, config, from_smiles, Settings, AMSJob, init

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/examples/BasisSetBenchmark/plams_workdir

Set Up Job Runner
~~~~~~~~~~~~~~~~~

Set up job runner, running as many jobs as possible in parallel.

.. code:: ipython3

   config.default_jobrunner = JobRunner(parallel=True, maxjobs=multiprocessing.cpu_count())
   config.job.runscript.nproc = 1

Set Up Molecules
~~~~~~~~~~~~~~~~

Create the molecules we want to use in our benchmark from SMILES.

.. code:: ipython3

   # The molecules we want to use in our benchmark:
   mol_smiles = {"Methane": "C", "Ethane": "C-C", "Ethylene": "C=C", "Acetylene": "C#C"}
   molecules = {}
   for name, smiles in mol_smiles.items():
       # Compute 10 conformers, optimize with UFF and pick the lowest in energy.
       molecules[name] = from_smiles(smiles, nconfs=10, forcefield="uff")[0]
       print(name, molecules[name])

::

   Methane   Atoms: 
       1         C       0.000000       0.000000       0.000000
       2         H       0.538912       0.762358      -0.599295
       3         H       0.731244      -0.596616       0.583182
       4         H      -0.567129      -0.670302      -0.678108
       5         H      -0.703028       0.504560       0.694220
     Bonds: 
      (1)--1.0--(2)
      (1)--1.0--(3)
      (1)--1.0--(4)
      (1)--1.0--(5)

   Ethane   Atoms: 
       1         C      -0.757196      -0.040522       0.044605
       2         C       0.757196       0.040522      -0.044605
       3         H      -1.205222       0.185290      -0.945970
       4         H      -1.130281       0.694397       0.788688
       5         H      -1.061719      -1.061491       0.357407
       6         H       1.205222      -0.185290       0.945971
       7         H       1.130281      -0.694396      -0.788689
       8         H       1.061719       1.061491      -0.357406
     Bonds: 
      (1)--1.0--(2)
      (1)--1.0--(3)
      (1)--1.0--(4)
      (1)--1.0--(5)
      (2)--1.0--(6)
      (2)--1.0--(7)
      (2)--1.0--(8)

   Ethylene   Atoms: 
       1         C       0.664485       0.027988      -0.023685
       2         C      -0.664485      -0.027988       0.023685
       3         H       1.253433      -0.878614       0.070299
       4         H       1.167038       0.980564      -0.156575
       5         H      -1.253433       0.878614      -0.070299
       6         H      -1.167038      -0.980564       0.156575
     Bonds: 
      (1)--2.0--(2)
      (1)--1.0--(3)
      (1)--1.0--(4)
      (2)--1.0--(5)
      (2)--1.0--(6)

   Acetylene   Atoms: 
       1         C      -0.587409       0.175060      -0.002211
       2         C       0.587409      -0.094463       0.002211
       3         H      -1.618985       0.411721      -0.006095
       4         H       1.618985      -0.331124       0.006094
     Bonds: 
      (1)--3.0--(2)
      (1)--1.0--(3)
      (2)--1.0--(4)

Initialize Calculation Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set up the settings which are common across jobs. The basis type is added later for each job.

.. code:: ipython3

   common_settings = Settings()
   common_settings.input.ams.Task = "SinglePoint"
   common_settings.input.ams.System.Symmetrize = "Yes"
   common_settings.input.adf.Basis.Core = "None"

.. code:: ipython3

   basis = ["QZ4P", "TZ2P", "TZP", "DZP", "DZ", "SZ"]
   reference_basis = "QZ4P"

Run Calculations
~~~~~~~~~~~~~~~~

.. code:: ipython3

   results = {}
   for bas in basis:
       for name, molecule in molecules.items():
           settings = common_settings.copy()
           settings.input.adf.Basis.Type = bas
           job = AMSJob(name=name + "_" + bas, molecule=molecule, settings=settings)
           results[(name, bas)] = job.run()

::

   [10.02|15:01:11] JOB Methane_QZ4P STARTED
   [10.02|15:01:11] JOB Ethane_QZ4P STARTED
   [10.02|15:01:11] JOB Ethylene_QZ4P STARTED
   [10.02|15:01:11] JOB Methane_QZ4P RUNNING
   [10.02|15:01:11] JOB Acetylene_QZ4P STARTED
   [10.02|15:01:11] JOB Methane_TZ2P STARTED
   [10.02|15:01:11] JOB Ethane_TZ2P STARTED
   [10.02|15:01:11] JOB Ethylene_TZ2P STARTED
   [10.02|15:01:11] JOB Acetylene_TZ2P STARTED
   [10.02|15:01:11] JOB Methane_TZP STARTED
   ... (PLAMS log lines truncated) ...

Results
~~~~~~~

Extract the energy from each calculation. Calculate the average absolute error in bond energy per atom for each basis set.

.. code:: ipython3

   average_errors = {}
   for bas in basis:
       if bas != reference_basis:
           errors = []
           for name, molecule in molecules.items():
               reference_energy = results[(name, reference_basis)].get_energy(unit="kcal/mol")
               energy = results[(name, bas)].get_energy(unit="kcal/mol")
               errors.append(abs(energy - reference_energy) / len(molecule))
               print("Energy for {} using {} basis set: {} [kcal/mol]".format(name, bas, energy))
           average_errors[bas] = sum(errors) / len(errors)

::

   [10.02|15:01:11] JOB Acetylene_TZP RUNNING
   [10.02|15:01:11] Waiting for job Methane_QZ4P to finish
   [10.02|15:01:11] JOB Methane_DZ RUNNING
   [10.02|15:01:11] JOB Ethane_DZP RUNNING
   [10.02|15:01:11] JOB Ethylene_DZP RUNNING
   [10.02|15:01:11] JOB Acetylene_DZP RUNNING
   [10.02|15:01:11] JOB Methane_SZ RUNNING
   [10.02|15:01:11] JOB Ethane_DZ RUNNING
   [10.02|15:01:11] JOB Ethane_SZ RUNNING
   [10.02|15:01:11] JOB Ethylene_DZ RUNNING
   [10.02|15:01:11] JOB Acetylene_DZ RUNNING
   ... (PLAMS log lines truncated) ...
   Energy for Methane using TZ2P basis set: -572.110159165262 [kcal/mol]
   [10.02|15:01:15] Waiting for job Ethane_QZ4P to finish
   Energy for Ethane using TZ2P basis set: -971.8820186844459 [kcal/mol]
   Energy for Ethylene using TZ2P basis set: -769.4329031250381 [kcal/mol]
   Energy for Acetylene using TZ2P basis set: -555.6672902509043 [kcal/mol]
   Energy for Methane using TZP basis set: -571.0448969099549 [kcal/mol]
   Energy for Ethane using TZP basis set: -970.0758887573307 [kcal/mol]
   Energy for Ethylene using TZP basis set: -767.3275176578105 [kcal/mol]
   [10.02|15:01:18] Waiting for job Acetylene_TZP to finish
   Energy for Acetylene using TZP basis set: -552.9562856742521 [kcal/mol]
   Energy for Methane using DZP basis set: -569.1190156251034 [kcal/mol]
   [10.02|15:01:19] Waiting for job Ethane_DZP to finish
   Energy for Ethane using DZP basis set: -966.0916443143674 [kcal/mol]
   Energy for Ethylene using DZP basis set: -764.4132984010868 [kcal/mol]
   Energy for Acetylene using DZP basis set: -550.6461805496328 [kcal/mol]
   Energy for Methane using DZ basis set: -560.9344313073021 [kcal/mol]
   [10.02|15:01:19] Waiting for job Ethane_DZ to finish
   Energy for Ethane using DZ basis set: -951.1666971758781 [kcal/mol]
   Energy for Ethylene using DZ basis set: -750.1745108422972 [kcal/mol]
   Energy for Acetylene using DZ basis set: -537.1008020388887 [kcal/mol]
   Energy for Methane using SZ basis set: -723.550123154895 [kcal/mol]
   Energy for Ethane using SZ basis set: -1216.914233427825 [kcal/mol]
   Energy for Ethylene using SZ basis set: -934.6558200110123 [kcal/mol]
   Energy for Acetylene using SZ basis set: -647.5029836817757 [kcal/mol]

.. code:: ipython3

   print("== Results ==")
   print("Average absolute error in bond energy per atom")
   for bas in basis:
       if bas != reference_basis:
           print("Error for basis set {:<4}: {:>10.3f} [kcal/mol]".format(bas, average_errors[bas]))

::

   == Results ==
   Average absolute error in bond energy per atom
   Error for basis set TZ2P:      0.170 [kcal/mol]
   Error for basis set TZP :      0.537 [kcal/mol]
   Error for basis set DZP :      1.024 [kcal/mol]
   Error for basis set DZ  :      3.339 [kcal/mol]
   Error for basis set SZ  :     27.683 [kcal/mol]

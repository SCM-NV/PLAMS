Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

   import sys
   import multiprocessing
   from scm.plams import JobRunner, config, from_smiles, Settings, AMSJob, init
   import numpy as np

   # this line is not required in AMS2025+
   init();

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
   jobs = []
   for bas in basis:
       for name, molecule in molecules.items():
           settings = common_settings.copy()
           settings.input.adf.Basis.Type = bas
           job = AMSJob(name=name + "_" + bas, molecule=molecule, settings=settings)
           jobs.append(job)
           results[(name, bas)] = job.run()

::

   [25.02|15:23:48] JOB Methane_QZ4P STARTED
   [25.02|15:23:48] JOB Ethane_QZ4P STARTED
   [25.02|15:23:48] JOB Ethylene_QZ4P STARTED
   [25.02|15:23:48] JOB Acetylene_QZ4P STARTED
   [25.02|15:23:48] JOB Methane_TZ2P STARTED
   [25.02|15:23:48] JOB Ethane_TZ2P STARTED
   [25.02|15:23:48] JOB Ethylene_TZ2P STARTED
   [25.02|15:23:48] JOB Acetylene_TZ2P STARTED
   [25.02|15:23:48] JOB Methane_TZP STARTED
   [25.02|15:23:48] JOB Ethane_TZP STARTED
   ... (PLAMS log lines truncated) ...

Results
~~~~~~~

Extract the energy from each calculation. Calculate the average absolute error in bond energy per atom for each basis set.

.. code:: ipython3

   try:
       # For AMS2025+ can use JobAnalysis class to perform results analysis
       from scm.plams import JobAnalysis

       ja = (
           JobAnalysis(jobs=jobs, standard_fields=["Formula", "Smiles"])
           .add_settings_field(("Input", "ADF", "Basis", "Type"), display_name="Basis")
           .add_field("NAtoms", lambda j: len(j.molecule))
           .add_field(
               "Energy", lambda j: j.results.get_energy(unit="kcal/mol"), display_name="Energy [kcal/mol]", fmt=".2f"
           )
           .sort_jobs(["NAtoms", "Energy"])
       )

       ref_ja = ja.copy().filter_jobs(lambda data: data["InputAdfBasisType"] == "QZ4P")

       ref_energies = {f: e for f, e in zip(ref_ja.Formula, ref_ja.Energy)}

       def get_average_error(job):
           return abs(job.results.get_energy(unit="kcal/mol") - ref_energies[job.molecule.get_formula()]) / len(
               job.molecule
           )

       ja.add_field("AvErr", get_average_error, display_name="Average Error [kcal/mol]", fmt=".2f")

       # Pretty-print if running in a notebook
       if "ipykernel" in sys.modules:
           ja.display_table()
       else:
           print(ja.to_table())

   except ImportError:

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

   [25.02|15:23:48] Waiting for job Methane_QZ4P to finish
   [25.02|15:23:48] JOB Ethylene_TZ2P RUNNING
   [25.02|15:23:48] JOB Acetylene_TZ2P RUNNING
   [25.02|15:23:48] JOB Ethane_TZ2P RUNNING
   [25.02|15:23:48] JOB Acetylene_TZP RUNNING
   [25.02|15:23:48] JOB Acetylene_DZ RUNNING
   [25.02|15:23:48] JOB Ethane_DZP RUNNING
   [25.02|15:23:48] JOB Ethylene_DZP RUNNING
   [25.02|15:23:48] JOB Methane_DZP RUNNING
   [25.02|15:23:48] JOB Methane_SZ RUNNING
   [25.02|15:23:48] JOB Ethane_TZP RUNNING
   ... (PLAMS log lines truncated) ...
   [25.02|15:23:50] Waiting for job Ethane_QZ4P to finish
   [25.02|15:23:52] Waiting for job Methane_TZP to finish
   [25.02|15:23:54] Waiting for job Ethane_DZP to finish

======= ====== ===== ====== ================= ========================
Formula Smiles Basis NAtoms Energy [kcal/mol] Average Error [kcal/mol]
======= ====== ===== ====== ================= ========================
C2H2    C#C    DZ    4      -537.10           4.91
C2H2    C#C    DZP   4      -550.65           1.53
C2H2    C#C    TZP   4      -552.96           0.95
C2H2    C#C    TZ2P  4      -555.67           0.27
C2H2    C#C    QZ4P  4      -556.76           0.00
C2H2    C#C    SZ    4      -647.50           22.69
CH4     C      DZ    5      -560.93           2.34
CH4     C      DZP   5      -569.12           0.70
CH4     C      TZP   5      -571.04           0.32
CH4     C      TZ2P  5      -572.11           0.10
CH4     C      QZ4P  5      -572.63           0.00
CH4     C      SZ    5      -723.55           30.18
C2H4    C=C    DZ    6      -750.17           3.37
C2H4    C=C    DZP   6      -764.41           1.00
C2H4    C=C    TZP   6      -767.33           0.51
C2H4    C=C    TZ2P  6      -769.43           0.16
C2H4    C=C    QZ4P  6      -770.41           0.00
C2H4    C=C    SZ    6      -934.66           27.37
C2H6    CC     SZ    8      -1216.91          30.49
C2H6    CC     DZ    8      -951.17           2.73
C2H6    CC     DZP   8      -966.09           0.87
C2H6    CC     TZP   8      -970.08           0.37
C2H6    CC     TZ2P  8      -971.88           0.14
C2H6    CC     QZ4P  8      -973.02           0.00
======= ====== ===== ====== ================= ========================

.. code:: ipython3

   print("== Results ==")
   print("Average absolute error in bond energy per atom")
   for bas in basis:
       if bas != reference_basis:
           if ja:
               av = np.average(ja.copy().filter_jobs(lambda data: data["InputAdfBasisType"] == bas).AvErr)
           else:
               av = average_errors[bas]
           print("Error for basis set {:<4}: {:>10.3f} [kcal/mol]".format(bas, av))

::

   == Results ==
   Average absolute error in bond energy per atom
   Error for basis set TZ2P:      0.170 [kcal/mol]
   Error for basis set TZP :      0.537 [kcal/mol]
   Error for basis set DZP :      1.024 [kcal/mol]
   Error for basis set DZ  :      3.339 [kcal/mol]
   Error for basis set SZ  :     27.683 [kcal/mol]

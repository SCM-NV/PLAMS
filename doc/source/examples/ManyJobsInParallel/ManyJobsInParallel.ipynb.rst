Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

   import multiprocessing
   from scm.plams import JobRunner, config, Settings, read_molecules, AMSJob, init

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/examples/ManyJobsInParallel/plams_workdir

Configure Job Runner
~~~~~~~~~~~~~~~~~~~~

Set the default job runner to run in parallel. Run as many jobs simultaneously as there are cpu on the system. In addition, set the number of cores for each job to 1.

.. code:: ipython3

   maxjobs = multiprocessing.cpu_count()
   print("Running up to {} jobs in parallel simultaneously".format(maxjobs))

::

   Running up to 12 jobs in parallel simultaneously

.. code:: ipython3

   config.default_jobrunner = JobRunner(parallel=True, maxjobs=maxjobs)

.. code:: ipython3

   config.job.runscript.nproc = 1

Load Molecules
~~~~~~~~~~~~~~

Load set of molecules from directory containing xyz files.

.. code:: ipython3

   molecules = read_molecules("molecules")

Set Up and Run Jobs
~~~~~~~~~~~~~~~~~~~

Configure the calculation settings in the ``Settings`` object. Run a geometry optimization job for each molecule in parallel.

.. code:: ipython3

   settings = Settings()
   settings.input.ams.Task = "GeometryOptimization"
   settings.input.dftb.Model = "GFN1-xTB"

.. code:: ipython3

   results = []
   for name, molecule in sorted(molecules.items()):
       job = AMSJob(molecule=molecule, settings=settings, name=name)
       results.append(job.run())

::

   [10.02|15:50:35] JOB Acetic_acid STARTED
   [10.02|15:50:35] JOB Benzene STARTED
   [10.02|15:50:35] JOB Butane STARTED
   [10.02|15:50:35] JOB Ethane STARTED
   [10.02|15:50:35] JOB Ethanol STARTED
   [10.02|15:50:35] JOB Butane RUNNING
   [10.02|15:50:35] JOB Formic_acid STARTED
   [10.02|15:50:35] JOB Acetic_acid RUNNING
   [10.02|15:50:35] JOB Methanol STARTED
   [10.02|15:50:35] JOB Water STARTED

Results
~~~~~~~

Print a table of results only for the successful calculations.

.. code:: ipython3

   # Only print the results of the succesful caluclations:
   for result in [r for r in results if r.ok()]:
       print("Energy for {:<12}: {:>10.3f} kcal/mol".format(result.name, result.get_energy(unit="kcal/mol")))

::

   [10.02|15:50:35] Waiting for job Acetic_acid to finish
   [10.02|15:50:35] JOB Formic_acid RUNNING
   [10.02|15:50:35] JOB Benzene RUNNING
   [10.02|15:50:35] JOB Water RUNNING
   [10.02|15:50:35] JOB Ethane RUNNING
   [10.02|15:50:35] JOB Ethanol RUNNING
   [10.02|15:50:35] JOB Methanol RUNNING
   [10.02|15:50:36] JOB Water FINISHED
   [10.02|15:50:36] JOB Formic_acid FINISHED
   [10.02|15:50:36] JOB Butane FINISHED
   [10.02|15:50:36] JOB Acetic_acid FINISHED
   ... (PLAMS log lines truncated) ...
   [10.02|15:50:36] Waiting for job Benzene to finish
   [10.02|15:50:37] Waiting for job Ethanol to finish
   Energy for Acetic_acid :  -9913.297 kcal/mol
   Energy for Benzene     : -12039.482 kcal/mol
   Energy for Butane      :  -8699.182 kcal/mol
   Energy for Ethane      :  -4686.354 kcal/mol
   Energy for Ethanol     :  -7629.287 kcal/mol
   Energy for Formic_acid :  -7890.662 kcal/mol
   Energy for Methanol    :  -5621.724 kcal/mol
   Energy for Water       :  -3618.401 kcal/mol

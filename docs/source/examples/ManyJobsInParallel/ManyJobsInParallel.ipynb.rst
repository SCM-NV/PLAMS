Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

    import multiprocessing
    from scm.plams import JobRunner, config, Settings, read_molecules, AMSJob

Configure Job Runner
~~~~~~~~~~~~~~~~~~~~

Set the default job runner to run in parallel. Run as many jobs
simultaneously as there are cpu on the system. In addition, set the
number of cores for each job to 1.

.. code:: ipython3

    maxjobs = multiprocessing.cpu_count()
    print("Running up to {} jobs in parallel simultaneously".format(maxjobs))


.. parsed-literal::

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

Configure the calculation settings in the ``Settings`` object. Run a
geometry optimization job for each molecule in parallel.

.. code:: ipython3

    settings = Settings()
    settings.input.ams.Task = "GeometryOptimization"
    settings.input.dftb.Model = "GFN1-xTB"

.. code:: ipython3

    results = []
    for name, molecule in sorted(molecules.items()):
        job = AMSJob(molecule=molecule, settings=settings, name=name)
        results.append(job.run())


.. parsed-literal::

    [13.08|11:29:05] JOB Acetic_acid STARTED
    [13.08|11:29:05] JOB Benzene STARTED
    [13.08|11:29:05] JOB Butane STARTED
    [13.08|11:29:05] JOB Ethane STARTED
    [13.08|11:29:05] JOB Ethanol STARTED
    [13.08|11:29:05] JOB Acetic_acid RUNNING
    [13.08|11:29:05] JOB Formic_acid STARTED
    [13.08|11:29:05] JOB Benzene RUNNING
    [13.08|11:29:05] JOB Methanol STARTED
    [13.08|11:29:05] JOB Butane RUNNING
    [13.08|11:29:05] JOB Water STARTED
    [13.08|11:29:05] JOB Ethane RUNNING
    [13.08|11:29:05] JOB Ethanol RUNNING


Results
~~~~~~~

Print a table of results only for the successful calculations.

.. code:: ipython3

    # Only print the results of the succesful caluclations:
    for result in [r for r in results if r.ok()]:
        print("Energy for {:<12}: {:>10.3f} kcal/mol".format(result.name, result.get_energy(unit="kcal/mol")))


.. parsed-literal::

    [13.08|11:29:05] JOB Formic_acid RUNNING
    [13.08|11:29:05] Waiting for job Acetic_acid to finish
    [13.08|11:29:05] JOB Methanol RUNNING
    [13.08|11:29:05] JOB Water RUNNING
    [13.08|11:29:05] JOB Acetic_acid FINISHED
    [13.08|11:29:05] JOB Acetic_acid SUCCESSFUL
    [13.08|11:29:05] Waiting for job Benzene to finish
    [13.08|11:29:05] JOB Benzene FINISHED
    [13.08|11:29:05] JOB Benzene SUCCESSFUL
    [13.08|11:29:05] Waiting for job Butane to finish
    [13.08|11:29:05] JOB Butane FINISHED
    [13.08|11:29:05] JOB Butane SUCCESSFUL
    [13.08|11:29:05] Waiting for job Ethane to finish
    [13.08|11:29:05] JOB Ethane FINISHED
    [13.08|11:29:05] JOB Ethane SUCCESSFUL
    [13.08|11:29:05] Waiting for job Ethanol to finish
    [13.08|11:29:06] JOB Ethanol FINISHED
    [13.08|11:29:06] JOB Ethanol SUCCESSFUL
    [13.08|11:29:06] Waiting for job Formic_acid to finish
    [13.08|11:29:06] JOB Formic_acid FINISHED
    [13.08|11:29:06] JOB Formic_acid SUCCESSFUL
    [13.08|11:29:06] Waiting for job Methanol to finish
    [13.08|11:29:06] JOB Methanol FINISHED
    [13.08|11:29:06] JOB Methanol SUCCESSFUL
    [13.08|11:29:06] Waiting for job Water to finish
    [13.08|11:29:06] JOB Water FINISHED
    [13.08|11:29:06] JOB Water SUCCESSFUL
    Energy for Acetic_acid :  -9913.297 kcal/mol
    Energy for Benzene     : -12039.482 kcal/mol
    Energy for Butane      :  -8699.182 kcal/mol
    Energy for Ethane      :  -4686.354 kcal/mol
    Energy for Ethanol     :  -7629.287 kcal/mol
    Energy for Formic_acid :  -7890.662 kcal/mol
    Energy for Methanol    :  -5621.724 kcal/mol
    Energy for Water       :  -3618.401 kcal/mol



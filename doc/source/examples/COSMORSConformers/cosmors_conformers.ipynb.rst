Worked Example
--------------

Initial imports
~~~~~~~~~~~~~~~

.. code:: ipython3

   from scm.plams.recipes.adfcosmorsconformers import ADFCOSMORSConfJob, ADFCOSMORSConfFilter
   from scm.plams import Molecule, from_smiles, Settings, init, config, JobRunner
   from scm.plams import ConformersJob

   # this line is not required in AMS2025+
   init();

::

   PLAMS working folder: /path/plams/examples/COSMORSConformers/plams_workdir.002

.. code:: ipython3

   config.default_jobrunner = JobRunner(parallel=True, maxjobs=8)  # Set the default jobrunner to be parallel
   config.default_jobmanager.settings.hashing = None  # Disable rerun prevention
   config.job.runscript.nproc = 1  # Number of cores for each job
   config.log.stdout = 1  # Suppress plams output

Set up conformer generator
~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we input the acetic acid molecule with the ``from_smiles`` function

.. code:: ipython3

   mol = from_smiles("CC(=O)O")

Now, we’ll specify a conformer generator (identical to the default) that generates only 50 initial structures:

.. code:: ipython3

   conf_sett = Settings()
   conf_sett.input.AMS.Generator.RDKit
   conf_sett.input.AMS.Generator.RDKit.InitialNConformers = 50
   conf_job = ConformersJob(name="conformers_uff", molecule=mol, settings=conf_sett)

Let’s also specify an additional step to add to the default workflow. Here, we’ll add a DFTB geometry optimization.

.. code:: ipython3

   dftb_sett = Settings()
   dftb_sett.input.DFTB
   dftb_sett.input.AMS.Task = "Optimize"

The final thing we need to specify are filters. Let’s make three filters, the first to take a maximum of 20 conformers with a maximum energy range of 22 kcal/mol, the second with 10 conformers and 12 kcal/mol and the third with 5 conformers and 7 kcal/mol.

.. code:: ipython3

   # ADFCOSMORSConfFilter(max number of conformers, max energy range)
   fil1 = ADFCOSMORSConfFilter(20, 22)  # applied to UFF
   fil2 = ADFCOSMORSConfFilter(10, 12)  # applied to DFTB
   fil3 = ADFCOSMORSConfFilter(5, 7)  # applied to ADF gas phase

Run COSMO-RS conformers job
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, we give this information to the ``ADFCOSMORSConfJob`` class. We also specify the name of the coskf files as well as the directory in which we’ll find them after the calculations complete.

.. code:: ipython3

   job = ADFCOSMORSConfJob(
       mol,
       conf_gen=conf_job,
       first_filter=fil1,
       additional=[(dftb_sett, fil2)],
       final_filter=fil3,
       coskf_name="acetic_acid",
       coskf_dir="test_coskfs",
   )
   job.run()
   job.results.wait()

::

   [12.09|16:17:48] JOB plamsjob STARTED
   [12.09|16:17:48] Waiting for job plamsjob to finish
   [12.09|16:17:48] JOB plamsjob/conformers_uff STARTED
   [12.09|16:17:48] JOB plamsjob/additional_1 STARTED
   [12.09|16:17:48] JOB plamsjob/adf_conformers STARTED
   [12.09|16:17:48] JOB plamsjob/adf_filter STARTED
   [12.09|16:17:48] Waiting for job adf_filter to finish
   [12.09|16:17:48] Waiting for job conformers_uff to finish
   [12.09|16:17:48] Waiting for job additional_1 to finish
   [12.09|16:17:48] Waiting for job adf_conformers to finish
   [12.09|16:17:52] JOB plamsjob/conformers_uff SUCCESSFUL
   [12.09|16:17:54] JOB plamsjob/additional_1 SUCCESSFUL
   [12.09|16:24:45] JOB plamsjob/adf_conformers SUCCESSFUL
   [12.09|16:24:45] JOB plamsjob/adf_filter SUCCESSFUL
   [12.09|16:24:45] JOB plamsjob/replay STARTED
   ... (PLAMS log lines truncated) ...

.. code:: ipython3

   !amsview test_coskfs/acetic_acid_0.coskf

.. code:: ipython3

   !amsview test_coskfs/acetic_acid_1.coskf

.. _ADFCOSMORSCompound:

ADF: Task COSMO-RS Compound
===========================

The |ADFCOSMORSCompoundJob| class generates results identical to the "Task COSMO-RS Compound" in the AMS ADF graphical user interface.  This python interface allows users to easily generate the .coskf files for one or many structures.  A possible usage is given in :ref:`ams_crs_workflow`.

Example: generating .coskf files for a set of compounds from xyz
------------------------------------------------------------------

:download:`Download compounds_xyz file: compounds_xyz.zip </examples/CRScosmorscompound/compounds_xyz.zip>`

The example will load all the molecules in the folder ``compounds_xyz`` and then optimize the gas geometry using ADF, and perform the ADF COSMO calculation for each compound. When the calculation is finished, we will find all the .coskf file in the ``test_coskfs`` directory.

.. code-block:: python

    from scm.plams import from_smiles, read_molecules
    from scm.plams import from_smiles, read_molecules
    from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob

    molecules = read_molecules('./compounds_xyz')

    for name, mol in molecules.items():
        job = ADFCOSMORSCompoundJob(molecule = mol, coskf_name = name, coskf_dir = 'test_coskfs')
        job.run()


Example: generating .coskf files for a set of compounds from smiles with parallel calculation
----------------------------------------------------------------------------------------------

First, we'll import the necessary classes and enable the parallel calculation through |JobRunner|. Here, we'll assign one core to each job, and we can have up to eight jobs running all at once.

.. code-block:: python

    from scm.plams import from_smiles, JobRunner, config
    from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob

    config.default_jobrunner = JobRunner(parallel=True, maxjobs=8)        # Set the default jobrunner to be parallel
    config.default_jobmanager.settings.hashing = None                     # Disable rerun prevention
    config.job.runscript.nproc = 1                                        # Number of cores for each job
    config.log.stdout = 1                                                 # Suppress plams output

Now, we will specify the smiles and name of a set of compounds and generate the initial geometry of each compound using |from_smiles| function. With the setting, ``nconfs=100`` and ``forcefield='uff'``, we will generate 100 conformers and find the one with the loweset energy using 'uff' forcefield. It's worth to notice that we can also generate a set of mutiple conformers through the |ADFCOSMORSConfJob| class as illustrate in the next `example <./ADFCOSMORSConformers.html>`__.

.. code-block:: python

    rd_smiles = ['O'  ,'CO']
    rd_names  = ['H2O','CO']
    molecules={}
    for name, smiles in zip(rd_names, rd_smiles):
        molecules[name] = from_smiles(smiles, nconfs=100, forcefield='uff')[0] #lowest energy one in 100 conformers

Lastly, we give this information to the |ADFCOSMORSCompoundJob| class, including the name of the coskf files as well as the directory in which we'll find them after the calculations complete.  Using the setting, ``preoptimization='GFN1-xTB'`` and ``singlepoint=False``, it will utilize the DFTB for a quick pre-optimization. Subsequently, it will execute a gas phase optimization using ADF, followed by the solvation calculation.

.. code-block:: python

    results = []
    for name, mol in molecules.items():
        job = ADFCOSMORSCompoundJob(
            molecule        = mol,                # The initial structure
            coskf_name      = name,               # a name to be used for coskf file
            coskf_dir       = 'test_coskfs',      # a directory to put the .coskf files generated
            preoptimization = 'GFN1-xTB',         # perform preoptimize or not 
            singlepoint     = False,              # run a singlepoint in gasphase and solvation calculation without geometry optimization. Cannot be combined with ``preoptimization``
            name            = name              ) # an optional name for the calculation directory
        results.append(job.run())

    finish()

In the ``test_coskfs`` directory, we will find the ``H2O.coskf`` and ``CO.coskf`` files.


The entire script can be seen/copied when expanded below.

.. raw:: html

   <details>
   <summary style="color:#008000;cursor:pointer">[show/hide code]</summary>

.. literalinclude:: /examples/CRScosmorscompound/parallel_COSMORScompound.py
    :language: python
    
.. raw:: html

     </details>



Source code for ``ADFCOSMORSCompound``
--------------------------------------

.. raw:: html

   <details>
   <summary style="color:#008000;cursor:pointer">[show/hide code]</summary>

.. literalinclude:: ../../../recipes/adfcosmorscompound.py

.. raw:: html

     </details>

Brief API Documentation
-----------------------

.. automodule:: scm.plams.recipes.adfcosmorscompound
    :no-special-members:
    :exclude-members: _result_type, __init__, new_children, postrun, _get_radii, adf_settings

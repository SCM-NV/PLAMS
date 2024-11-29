.. _General:

Getting Started
===============

This section contains introductory information about installing and running PLAMS.

For quick-start guides on a wider range of topics within PLAMS, see the :ref:`examples`.

Overview
--------

PLAMS (Python Library for Automating Molecular Simulation) is a flexible and extensible toolkit for streamlining molecular simulation workflows.

It simplifies and automates the process of configuring, running and analyzing computational chemistry calculations.
The key features of PLAMS are:

- **Amsterdam Modeling Suite (AMS) Integration**: Full support for interacting with AMS programs
- **Parallel Processing**: Run jobs in parallel without any need for separate parallelization scripts
- **Scheduler Integration**: Integration with job schedulers such as SLURM making large-scale computations easier to manage
- **Automatic File and Folder Organization**: PLAMS automatically handles file organization, preventing overwrites and ensuring clean data flows
- **Controllable Re-runs and Restarts**: Efficiently manage job executions by preventing redundant runs and easily restarting from crash points if needed
- **Output processing**: Extract, post-process, and analyze results, ensuring that only relevant data is used for further calculations or workflows
- **Compatibility with Chemistry Tools**: Includes built-in interfaces for popular programs and packages such as ASE, RDKit, Dirac, ORCA, CP2K, DFTB+ and Crystal and more

Quick Start
-----------

PLAMS is available to all users of AMS "out of the box" as part of the `AMS Python Stack <../Scripting/Python_Stack/Python_Stack.html>`__, which can be accessed with the ``$AMSBIN/amspython`` command.

For most use-cases, no specific installation outside of AMS is required. For usage outside of ``amspython``, please see the :ref:`installation guide <installation>` below.

To get started with PLAMS, import ``scm.plams`` into your python script or jupyter notebook.
Then, follow one of the :ref:`examples <examples>` to help create your script.

For example, the following is based upon :ref:`Geometry Optimization of Water <WaterOptimizationExample>`,
and also makes use of `PISA <../pisa/index.html>`__, also included in ``amspython``.

.. code:: ipython3

    # water_opt.py
    from scm.plams import from_smiles, AMSJob
    from scm.input_classes import drivers, engines

    water = from_smiles("O")

    driver = drivers.AMS()
    driver.Task = "GeometryOptimization"
    driver.Properties.NormalModes = "Yes"
    driver.Engine = engines.ForceField()
    driver.Engine.Type = "UFF"

    job = AMSJob(molecule=water, settings=driver, name="water_opt")
    results = job.run()

    print("Optimized geometry:")
    print(results.get_main_molecule())

Running the command ``$AMSBIN/amspython water_opt.py`` produces the successful output:

.. parsed-literal::

    JOB water_opt RUNNING
    JOB water_opt FINISHED
    JOB water_opt SUCCESSFUL
    Optimized geometry:
      Atoms:
        1         O      -0.000360       0.403461       0.000000
        2         H      -0.783821      -0.202431       0.000000
        3         H       0.784180      -0.201030       0.000000
      Bonds:
       (1)--1.0--(2)
       (1)--1.0--(3)

For more advanced workflows including usage of other AMS engines, see the other :ref:`examples <examples>`.


.. _installation:

Installation Guide
------------------

PLAMS and all its required and optional dependencies are included as part of the `AMS Python Stack <../Scripting/Python_Stack/Python_Stack.html>`__.
This is the easiest way to use PLAMS, as it requires no additional installation process.

However, if you want to use PLAMS outside of ``amspython``, since ``AMS2024.103`` PLAMS is available on `PyPI <https://pypi.org/project/plams>`__
and so can be installed via the ``pip`` python package installer.

To install the latest version of PLAMS into your python environment, simply run ``pip install plams``.
To install a specific version of PLAMS (e.g. ``2025.101``), run ``pip install plams==2025.101``.

By default, PLAMS only installs a minimal set of required packages on installation using pip.
For additional functionality, further optional packages are required.
Since ``AMS2025``, these are available for installation through extra dependency groups with pip.

The available groups are:

- **chem**: for chemistry packages such as ``RDKit``, ``ase``
- **analysis**: for packages used to analyse and plot results of calculations e.g. ``scipy``, ``matploblib``, ``networkx``
- **ams**: for technical packages for use with the AMS interface

One or more of these can be installed using the command ``pip install 'plams[chem,analysis,ams]'``.

Users of the AMS will also have to install the ``scm.amspipe`` package using the command ``pip install $AMSHOME/scripting/scm/amspipe``.


What's new in PLAMS for AMS2025?
--------------------------------------

Added
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Support for AMS ``ChemicalSystem`` within |AMSJob| and |AMSResults|. |AMSJob| can accept a ``ChemicalSystem`` as an input system, and the methods :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_system`, :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_input_system` and :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_main_system` on |AMSResults| return a ``ChemicalSystem``. These provide the option to use a ``ChemicalSystem`` in place of a PLAMS ``Molecule``.
* Support for work functions through :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_work_function_results` and :func:`~scm.plams.tools.plot.plot_work_function`

Changed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Calling |init| and |finish| functions in a script is now optional
* Functions for optional packages (e.g. RDKit, ASE) are available even when these packages are not installed, but will raise an |MissingOptionalPackageError| when called
* :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_main_ase_atoms` also includes atomic charges
* Global ``config`` is initialized with |ConfigSettings| instead of loading from the standard ``plams_defaults`` file (see |global-settings|)

* :attr:`~scm.plams.core.basejob.Job.status` is a ``JobStatus`` string enum
* Supercell and RDKit properties are no longer serialized to AMS input

Fixed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* ``charge`` property on a |Molecule| is a numeric instead of string type when loading molecule from a file
* :meth:`~scm.plams.mol.molecule.Molecule.delete_all_bonds` removes the reference molecule from the removed bond instances
* :meth:`~scm.plams.core.basejob.SingleJob.load` returns the correctly loaded job
* :meth:`~scm.plams.interfaces.adfsuite.ams.AMSJob.check` handles a ``NoneType`` status, returning ``False``

Deprecated
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* PLAMS launch script is deprecated in favour of simply running with ``amspython``

Removed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Legacy ``BANDJob``, ``DFTBJob``, ``UFFJob``, ``MOPACJob``, ``ReaxFFJob``, ``CSHessianADFJob`` and ``ADFJob`` have been removed. These were deprecated since AMS2019, and replaced by |AMSJob|.
* Exception classes ``AMSPipeDecodeError``, ``AMSPipeError``, ``AMSPipeInvalidArgumentError``, ``AMSPipeLogicError``, ``AMSPipeRuntimeError``, ``AMSPipeUnknownArgumentError``, ``AMSPipeUnknownMethodError``, ``AMSPipeUnknownVersionError``, were moved from ``scm.plams`` to ``scm.amspipe``.

What's new in PLAMS for AMS2024?
--------------------------------------

* :ref:`Packmol interface <PackmolInterface>` has been extended to pack in crystal voids and to get the total system charge from the sum of the constituent molecules

* Additions to |AMSResults|: get_normal_modes(), get_polarizability(), get_ir_spectrum(), get_ir_spectrum_md(), get_frequency_spectrum(), get_force_constants()

* Additions to |Molecule|: get_moments_of_inertia(), get_gyration_radius(), align2mol()

What's new in PLAMS for AMS2023?
--------------------------------------

* The :ref:`AMSCalculator` class for running any AMS engine with ASE (see: :ref:`ASECalculatorExample`)

* Classes for calculating :ref:`reduction and oxidation potentials  <RedoxExample>` with ADF and optionally COSMO-RS

* The :ref:`ADFCOSMORSCompoundJob <ADFCOSMORSCompound>` class for running jobs equivalent to "Task COSMO-RS Compound" in the AMS GUI. Such a job generates a .coskf file for use with COSMO-RS.

* The calculation of the :ref:`vibronic density of states<fcf_dos>` has been added to PLAMS.

* Classes for running and restarting :ref:`molecular dynamics (MD) jobs with AMS <AMSMDJob>`

* A class for generating and analyzing :ref:`conformers <conformers_interface>`

* :ref:`Quick jobs <Quickjobs>`, like for example the ``preoptimize()`` function let you quickly optimize a Molecule

* :ref:`Packmol interface <PackmolInterface>` for generating liquid and gas mixtures, solid-liquid interfaces, and microsolvation spheres

* :ref:`FileFormatConversionTools` for converting VASP, Gaussian, or Quantum ESPRESSO output to ams.rkf and engine.rkf files that can be opened with the AMS GUI

* :ref:`PlottingTools` for plotting a molecule or ASE Atoms inside a Jupyter notebook

* :ref:`PlottingTools` for plotting the :ref:`electronic band structure <BandStructureExample>`

* Additions to |AMSResults|: get_homo_energies(), get_lumo_energies, get_smallest_homo_lumo_gap()

* Additions to |Molecule|: guess_atomic_charges(), set_density(), get_unique_bonds(), get_unique_angles()

* Many new :ref:`examples`

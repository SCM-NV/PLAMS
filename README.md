PLAMS - Python Library for Automating Molecular Simulation
==========================================================

Overview
--------

PLAMS is a flexible and extensible toolkit for streamlining molecular simulation workflows.

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

PLAMS is available to all users of AMS "out of the box" as part of the [AMS Python Stack](https://www.scm.com/doc/Scripting/Python_Stack/Python_Stack.html), which can be accessed with the ``$AMSBIN/amspython`` command.

For most use-cases, no specific installation outside of AMS is required. For usage outside of `amspython`, please see the installation guide below.

To get started with PLAMS, import `scm.plams` into your python script or jupyter notebook.
Then, follow one of the examples to help create your script e.g.

```python
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
```

Running the command `$AMSBIN/amspython water_opt.py` produces the successful output:

```
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
```

For more advanced workflows including usage of other AMS engines, see the other examples.


Installation Guide
------------------

As mentioned, PLAMS and all its required and optional dependencies are included as part of the AMS python stack.
This is the easiest way to use PLAMS, as it requires no additional installation process.

However, if you want to use PLAMS outside of `amspython`, since `AMS2024.103` PLAMS is available on [PyPI](https://pypi.org/project/plams).
and so can be installed via the `pip` python package installer.

To install the latest version of PLAMS into your python environment, simply run `pip install plams`.
To install a specific version of PLAMS (e.g. `2025.101`), run `pip install plams==2025.101`.

By default, PLAMS only installs a minimal set of required packages on installation using pip.
For additional functionality, further optional packages are required.
Since `AMS2025`, these are available for installation through extra dependency groups with pip.

The available groups are:

- **chem**: for chemistry packages such as `RDKit`, `ase`
- **analysis**: for packages used to analyse and plot results of calculations e.g. `scipy`, `matploblib`, `networkx`
- **ams**: for technical packages for use with the AMS interface

One or more of these can be installed using the command `pip install 'plams[chem,analysis,ams]'`.

Users of the AMS will also have to install the `scm.amspipe` package using the command `pip install $AMSHOME/scripting/scm/amspipe`.

A final option is to download PLAMS directly from the [GitHub page](https://github.com/SCM-NV/PLAMS).
[Released versions](https://github.com/SCM-NV/PLAMS/releases) are available since `AMS2024.103`.
The latest (unreleased) development version can be downloaded from the [trunk branch](https://github.com/SCM-NV/PLAMS/archive/refs/heads/trunk.zip).
Once the downloaded zip file has been extracted, navigate to its location and run `pip install .` to install into your python environment.

License
-------

See [LICENSE.md](https://github.com/SCM-NV/PLAMS/blob/trunk/LICENSE.md) for details.

Contributing
------------

Contributions are welcome. See [CONTRIBUTING.md](https://github.com/SCM-NV/PLAMS/blob/trunk/CONTRIBUTING.md) for details.
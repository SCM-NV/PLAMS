RASPA
-------------------------

(*contributed by* `Patrick Melix <https://www.researchgate.net/profile/Patrick_Melix>`_\)

.. currentmodule:: scm.plams.interfaces.thirdparty.raspa

The runscript written by PLAMS expects the environment variable *$RASPA_DIR* to be set. It should point to the 
location RASPA_ is installed in (*$RASPA_DIR/bin/simulate* being the compiled executable).

Input
~~~~~~~~~~~~~~

The RASPA_ input is not easily wrapped into a |Settings| class because it does not know section endings and
sections don't have a common open/close statement.
Like the other interfaces, the RASPA_ input file is generated using the ``input`` branch of the job settings.
An example is given here::

    settings = Settings()
    settings.input.simulationtype = 'MonteCarlo'
    settings.input.numberofcycles = 100
    settings.input.printevery = 10
    settings.input.forcefield = 'GenericMOFs'
    settings.input.usechargesfromciffile = True
    settings.input.framework._h = 0
    settings.input.framework.frameworkname = 'someMOF'
    settings.input.framework.unitcells = '4 4 4'
    settings.input.externaltemperature = 298.0
    settings.input.component['0'] = Settings()
    settings.input.component['0'].moleculename = 'helium'
    settings.input.component['0'].moleculedefinition = 'TraPPE'
    settings.input.component['0'].widomprobability = 1.0
    settings.input.component['0'].createnumberofmolecules = 0
    settings.input.component._1 = Settings()
    settings.input.component._1.moleculename = 'methane'
    settings.input.component._1.moleculedefinition = 'TraPPE'
    settings.input.component._1.widomprobability = 1.0
    settings.input.component._1.createnumberofmolecules = 0


The input generated during the execution of the RASPA_ job is similar to: ::

    numberofcycles 100
    printevery 10
    simulationtype MonteCarlo
    usechargesfromciffile yes

    componentexternaltemperature 298.0
    forcefield GenericMOFs
    framework 0
      frameworkname someMOF
      unitcells 4 4 4

    component 0 MoleculeName helium
      createnumberofmolecules 0
      moleculedefinition TraPPE
      widomprobability 1.0

    component 1 MoleculeName methane
      createnumberofmolecules 0
      moleculedefinition TraPPE
      widomprobability 1.0


The *Component* sections differ from the standard behaviour and are printed by the input parser as shown in this example.
They can be inserted using the dictionary-like or dot-like notation as shown above. The ``MoleculeName`` key is obligatory.

As in other interfaces, the ``_h`` key results in the value being printed along the section title.
No other special sections are defined.

If you need some files to be copied to the actual execution directory, pass them to the constructor using the ``copy=`` option.
Alternatively you can symlink them using the ``symlink=`` option. See the API below.


Molecule parsing
~~~~~~~~~~~~~~~~~~~~

Molecule handling is not supported yet. Prepare the needed files and pass them to the constructor using the ``copy`` or ``symlink`` option.

For a more detailed description of the *RASPA* input see the documentation in the RASPA_ GitHub repository.

.. _RASPA: https://github.com/iRASPA/RASPA2

Loading jobs
~~~~~~~~~~~~~~~~~~~

Calculations done without PLAMS can be loaded using the |load_external| functionality. The |RaspaResults| class **does not** support reading input files into |Settings| objects yet.

Just do ``RaspaJob.load_external(path)`` to get use the job inside PLAMS.


API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: RaspaJob

.. autoclass:: RaspaResults

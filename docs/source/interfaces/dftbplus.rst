DFTB+
-------------------------

(*contributed by* `Patrick Melix <https://www.researchgate.net/profile/Patrick_Melix>`_\)

.. currentmodule:: scm.plams.interfaces.thirdparty.dftbplus

DFTB+ is a density-functional tight-binding implementation.
More information about DFTB can be found on its `official website <http://www.dftb-plus.info>`_.

PLAMS offers a simple and incomplete DFTB+ interface.
It is so far capable of handling molecular calculations.
The relevant classes are |DFTBPlusJob| and |DFTBPlusResults|.


Preparing a calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Preparing an instance of |DFTBPlusJob| follows the general principles for |SingleJob|.
Information adjusting the input file is stored in the ``myjob.settings.input`` branch.
The geometry of your system can be supplied via the class |Molecule|.
If the Atomic Simulation Environment (`ASE <https://wiki.fysik.dtu.dk/ase/index.html>`_) is installed, the molecule is written using its engine.
Otherwise this interface provies a basic routine.
Note that the molecule is in this case transformed into the ``GenFormat`` with the ``C`` (cluster) or ``S`` (supercell) option, meaning the internal routine can handle clusters and supercell systems.
The option ``F`` for fractional coordinates is not available! See `the manual <http://www.dftb-plus.info/documentation/>`_ for further information on the different geometry-input types.




.. _dftbplus-input:

Input
+++++

Input files for DFTB+ are either in HSD (*human-friendly structured data*) or XML format.
This interface will produce HSD format input files.
See `the manual <http://www.dftb-plus.info/documentation/>`_ for further information on keywords and structure.
The input file must be named *dftb_in.hsd* and is therefore created using this name.
Note that many values have standard settings, those will all be printed to *dftb_pin.hsd* when you start a calculation.
Check both files for errors when having problems.

HSD input files are organized using different properties.
Each property is represented by a key.
The key can be of type *logical*, *integer*, *real*, *string*, *property list* or *method type*.
The last two begin and end with curly brackets.
Since *property lists* and *method types* are very similar in notation we require a way of representing that style in the tree-like structure of |Settings|.
This is done by automatically creating either *property lists* or *method types* whenever a key has a subkey.

*   If a key has a subkey and a key ``_h``, the string assigned to ``_h`` will be used as a method name.
*   If a key has one or more subkeys it will be created as a *property list*::

        #sets the hamiltonian property to be of a method type named DFTB
        myjob.setting.input.hamiltonian._h = 'DFTB'
        #sets the key parserversion of the property list parserversion to 4
        myjob.setting.input.parseroptions.parserversion = '4'

Empty *method types* can be created by not giving any subkeys to a key except ``_h``.



.. _dftbplus-runscript:

Runscript
+++++++++

The runscript will call the binary ``dftb+``, so make sure it is in your ``$PATH``.
No options are supported.
The standard output is redirected to ``$JN.out``, errors to ``$JN.err``.


Results extraction
~~~~~~~~~~~~~~~~~~~~~~~~~
DFTB+ creates multiple outputfiles, none of them are renamed.
See ``detailed.out`` for the results of your calculation.
Resulting geometries are saved in ``.xyz`` and ``.gen`` format by DFTB+.
Other files might be created depending on your calculation type.

General text processing methods from |Results| can be used to obtain data from results files.
At the moment only three functions for result extraction are defined:

*   Read the total energy :meth:`~DFTBPlusResults.get_energy`
*   Get the molecule from the ``.xyz`` file :meth:`~DFTBPlusResults.get_molecule`
*   Read the atomic charges :meth:`~DFTBPlusResults.get_atomic_charges`

Example
~~~~~~~~~~~~~~~~~~~~~~~~~
::

    common = Settings()

    common.input.driver._h = 'ConjugateGradient'
    common.input.hamiltonian._h = 'DFTB'
    common.input.hamiltonian.scc = 'yes'
    common.input.hamiltonian.mixer._h = 'Broyden'
    common.input.hamiltonian.mixer.mixingparameter = '0.2'
    common.input.hamiltonian.slaterkosterfiles._h = 'Type2Filenames'
    common.input.hamiltonian.slaterkosterfiles.prefix = '"~/SLAKO/mio-1-1/"'
    common.input.hamiltonian.slaterkosterfiles.separator = '"-"'
    common.input.hamiltonian.slaterkosterfiles.suffix = ".skf"
    common.input.hamiltonian.slaterkosterfiles.lowercasetypename = 'No'
    common.input.hamiltonian.maxangularmomentum.c = '"p"'
    common.input.hamiltonian.maxangularmomentum.h = '"s"'
    common.input.parseroptions.parserversion = '4'

    mol = Molecule(filename='mol.xyz') # read Molecule from mol.xyz

    job = DFTBPlusJob(name='plamstest', molecule=mol, settings=common)
    jobres = job.run()

    energy = jobres.get_energy(string='Fermi energy', unit='ev')
    print(energy)

    mol = jobres.get_molecule()
    print(mol)

    atomic_charges = jobres.get_atomic_charges()
    print(atomic_charges)


API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: DFTBPlusJob(molecule=None, name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type
.. autoclass:: DFTBPlusResults

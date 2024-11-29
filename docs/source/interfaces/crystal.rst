Crystal
-------------------------

(*contributed by* `Patrick Melix <https://www.researchgate.net/profile/Patrick_Melix>`_\)

.. currentmodule:: scm.plams.interfaces.thirdparty.crystal

More information about CRYSTAL can be found on its `official website <http://www.crystal.unito.it>`_.

PLAMS offers a simple CRYSTAL interface which does not offer access to all possible input types of CRYSTAL just most.
CRYSTAL14 was used by the developer, but as far as the developer can tell the new input features from CRYSTAL17 can be achieved with this interface.
Older CRYSTAL versions are more restrictive with the input, so they have not been tested.
The relevant classes are |CrystalJob| and |Results|.


Preparing a calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Preparing an instance of |CRYSTALJob| follows the general principles for |SingleJob|.
Information adjusting the input file is stored in the ``myjob.settings.input`` branch.
The geometry of your system can be supplied via the class |Molecule| and will be parsed into a ``fort.34`` file using the *EXTERNAL* keyword (only if you have `ASE <https://wiki.fysik.dtu.dk/ase/index.html>`_ installed).
It can also be supplied to the ``myjob.settings.input`` branch by using the function |mol2CrystalConf| to create a CRYSTAL-type input of your structure inside the input file (set ``myjob.settings.ignore_molecule = True``).
Consult `the manual <http://http://www.crystal.unito.it/documentation.php>`_ for further information on the different input options of CRYSTAL.



.. crystal-input:

Input
+++++

Settings must contain at least (case insensitive):

- one geometry key ('CRYSTAL','SLAB','POLYMER','HELIX','MOLECULE','EXTERNAL','DLVINPUT'), if the |Molecule| parser of |CRYSTALJob| is used the 'EXTERNAL' keyword is added automatically. For using |mol2CrystalConf|, see below).
- one basis key ('BASISSET')
- one option_key ('options' and anything else)

The ordering inside the geometry-, basis- and option-block can be controlled::

    settings.input.basisset._h = 'top line'
    settings.input.basisset._1 = 'first line'
    settings.input.basisset._2 = 'second line'

and so on.
Note that you can also pass lists to the ordered version.
Every item will end up in one line.

To make input nicer, the 'options' key will never be printed, since the input does not allow an opening statement for this block.
This way you can use::

    settings.input.options.bla = 'bla'
    settings.input.options.test = ''
    settings.input.options._h = 'FIRST'

without the 'options' being printed, but the section will still be closed with an 'END'.


.. crystal-runscript:


Runscript
+++++++++

The command ``crystal`` should point to the crystal binary or a runscript (so make sure it is in your ``$PATH``), that the input can be piped to.
Modify ``CrystalJob._command`` if necessary.
PLAMS will not clean up the mess of files that crystal produces.
If you want that, your runscript should do it for you.
Standard output is written to ``$JN.out``.

Molecule parsing
~~~~~~~~~~~~~~~~~~~~~~~~~
The |Molecule| parser of |CRYSTALJob| makes use of the ``EXTERNAL`` keyword by writing the |Molecule| instance to the file *fort.34*.
This is done using the ``write_crystal`` method from ASE.

If you do not have ASE installed you cannot use this feature.
In this case use the function |mol2CrystalConf| to create a list to pass to ``myjob.settings.input.<geomKey>``.
Remember to set ``myjob.settings.ignore_molecule = True`` if you are doing this.


Results extraction
~~~~~~~~~~~~~~~~~~~~~~~~~
There is no special results extraction yet, use the standard methods from the |Results| class.

Example
~~~~~~~~~~~~~~~~~~~~~~~~~
::

     mol = Molecule('test.xyz')
     common = Settings()
     geom = ['0 0 0',
     '194',
     '2.456 6.696',
     '2',
     '6 0.0     0.0     0.25',
     '6 0.333333  0.666667  0.75',
     'SLAB',
     '0 0 1',
     '1 1',
     'OPTGEOM',
     'FULLOPTG',
     'ENDGEOM']

     common.input.basisset = 'STO-3G'
     common.input.options.shrink = '9 18'
     common.input.options.scfdir = ''
     common.input.options._h = 'RHF'
     common.input.options.dft.exchange = 'PBE'
     common.input.options.dft.correlat = 'PBE'
     common.input.options.maxcycle = 250
     common.input.options.fmixing = 90
     #common.input.options.test = True

     common.input.crystal = geom
     common.ignore_molecule = True

     #alternative 1
     #common.input.crystal = mol2CrystalConf(mol)
     #common.ignore_molecule = True

     #alternative 2
     #just pass the molecule to CrystalJob if you have ASE

     job = CrystalJob(name='crystaltest', settings=common, molecule=mol)
     jobres = job.run()


API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: CrystalJob(name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type
.. autofunction:: mol2CrystalConf(molecule)

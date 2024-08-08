.. _band-fragment-recipe:

BAND fragment job
--------------------

.. currentmodule:: scm.plams.recipes.bandfragment

In this module a dedicated job type for Energy Decomposition Analysis in BAND is defined.
Such an analysis is performed on a periodic system divided into 2 fragments and consists of a minimum of 3 separate Band runs: one for each fragment and one for full system. See 

We define a new job type |BANDFragmentJob|, by subclassing |ADFFragmentJob|, which in turn is a subclass of |MultiJob|.
The constructor (``__init__``) of this new job takes 2 more arguments (``fragment1`` and ``fragment2``) and one optional argument ``full_settings`` for additional input keywords that are used **only** in the full system calculation. Furthermore, you can specify
the optimized fragment geometries using ``fragment1_opt`` and ``fragment2_opt`` for single-point calculations to also obtain the preparation energies.

In the |prerun| method two fragment jobs and the full system job are created with the proper settings and molecules.
They are then added to the ``children`` list.

The dedicated |Results| subclass for |BANDFragmentJob| does not provide too much additional functionality.
It simply redirects the usual |AMSResults| methods to the results of the full system calculation. The pEDA results can be obtained using the ``get_energy_decomposition`` method. It will return a dictionary with the available energy decomposition terms.

The source code of the whole module with both abovementioned classes:

API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: BANDFragmentJob()
.. autoclass:: BANDFragmentResults()


Example
~~~~~~~~~~~~~~~~~~~~~~~~~

A basic example using `ASE to build a surface slab <https://wiki.fysik.dtu.dk/ase/ase/build/surface.html>`_ and perform a BAND fragment calculation:

.. code-block:: python
    
    #build the surface
    from ase.build import fcc111, add_adsorbate
    mol = fcc111('Au', size=(2,2,3))
    add_adsorbate(mol, 'H', 1.5, 'ontop')
    mol.center(vacuum=10.0, axis=2)

    # separate the fragments
    surface = mol.copy()
    symbols = surface.get_chemical_symbols()
    del surface[[i for i in range(len(symbols)) if symbols[i] != 'Au']]
    adsorbate = mol.copy()
    del adsorbate[[i for i in range(len(symbols)) if symbols[i] == 'Au']]

    # optional: load the optimized molecules
    from ase import io
    surface_opt = io.read('surface_opt.xyz')
    adsorbate_opt = io.read('adsorbate_opt.xyz')
    assert len(surface_opt) == len(surface)
    assert len(adsorbate_opt) == len(adsorbate)

    # settings for job
    base_settings = Settings()
    base_settings.input.ams.task = 'SinglePoint'
    base_settings.input.band.basis.type = 'TZP'
    base_settings.input.band.basis.core = 'Medium'
    base_settings.input.band.dos.calcdos = 'No'
    base_settings.input.band.kspace.regular.numberofpoints = '5 5 1'
    base_settings.input.band.beckegrid.quality = 'Good'
    base_settings.input.band.zlmfit.quality = 'Good'
    base_settings.input.band.usesymmetry = False
    base_settings.input.band.xc.gga = 'PBE'
    base_settings.input.band.xc.dispersion = 'Grimme4'

    eda_settings = plams.Settings()
    eda_settings.input.band.peda = ''

    eda_job = BANDFragmentJob(fragment1=fromASE(surface), fragment2=fromASE(adsorbate),
                            settings=base_settings, full_settings=eda_settings,
                            fragment1_opt=fromASE(surface_opt), fragment2_opt=fromASE(adsorbate_opt))

    results = eda_job.run()
    eda_res = job.results.get_energy_decomposition()
    print("{:<20} {:>10}".format('Term', 'Energy [kJ/mol]'))
    for key, value in eda_res.items():
        print("{:<20} {:>10.4f}".format(key, value))

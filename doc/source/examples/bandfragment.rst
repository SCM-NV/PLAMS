.. _band-fragment-recipe:

BAND fragment job
--------------------

.. currentmodule:: scm.plams.recipes.bandfragment

In this module a dedicated job type for Energy Decomposition Analysis in BAND is defined.
Such an analysis is performed on a periodic system divided into 2 fragments and consists of a minimum of 3 separate BAND runs: one for each fragment and one for full system. See also |ADFFragmentJob|.

We define a new job type |BANDFragmentJob|, by subclassing |ADFFragmentJob|, which in turn is a subclass of |MultiJob|.
The constructor (``__init__``) of this new job takes 2 more arguments (``fragment1`` and ``fragment2``) and one optional argument ``full_settings`` for additional input keywords that are used **only** in the full system calculation. Furthermore, you can specify
the optimized fragment geometries using ``fragment1_opt`` and ``fragment2_opt`` for single-point calculations to also obtain the preparation energies.

In the |prerun| method two fragment jobs and the full system job are created with the proper settings and molecules.
They are then added to the ``children`` list.

The dedicated |Results| subclass for |BANDFragmentJob| does not provide too much additional functionality.
It simply redirects the usual |AMSResults| methods to the results of the full system calculation. The pEDA results can be obtained using the ``get_energy_decomposition`` method. It will return a dictionary with the available energy decomposition terms.

A derived subclass |NOCVBandFragmentJob| is also provided. It can be usefull for generating NOCV plots after the PEDA-NOCV calculation.

The source code of the whole module with both abovementioned classes:

API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: BANDFragmentJob()
.. autoclass:: BANDFragmentResults()
.. autoclass:: NOCVBandFragmentJob()


Example
~~~~~~~~~~~~~~~~~~~~~~~~~

A basic example using `ASE to build a surface slab <https://wiki.fysik.dtu.dk/ase/ase/build/surface.html>`_ and perform a BAND fragment calculation: :download:`bandfrag_test.py <../../../examples/scripting/bandfrag_test.py>`

.. literalinclude:: ../../../examples/scripting/bandfrag_test.py
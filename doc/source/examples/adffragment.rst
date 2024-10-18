.. _adf-fragment-recipe:

ADF fragment job
--------------------
.. currentmodule:: scm.plams.recipes.adffragment

In this module a dedicated job type for ADF fragment analysis is defined.
Such an analysis is performed on a molecular system divided into 2 fragments and consists of 3 separate ADF runs: one for each fragment and one for full system.

We define a new job type ``ADFFragmentJob`` by extending |MultiJob|.
The constructor (``__init__``) of this new job takes 2 more arguments (``fragment1`` and ``fragment2``) and one optional argument ``full_settings`` for additional input keywords that are used **only** in the full system calculation.

In the |prerun| method two fragment jobs and the full system job are created with the proper settings and molecules.
They are then added to the ``children`` list.

The dedicated |Results| subclass for ``ADFFragmentJob`` does not provide too much additional functionality.
It simply redirects the usual |AMSResults| methods to the results of the full system calculation.


API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: ADFFragmentJob()
.. autoclass:: ADFFragmentResults()

An example usage: (:download:`ethene.xyz <../../../examples/ADFFrag/ethene.xyz>`,
:download:`butadiene.xyz <../../../examples/ADFFrag/butadiene.xyz>`,
:download:`adffrag_test.py <../../../examples/ADFFrag/adffrag.py>`)

.. literalinclude:: ../../../examples/ADFFrag/ethene.xyz
.. literalinclude:: ../../../examples/ADFFrag/butadiene.xyz
.. literalinclude:: ../../../examples/ADFFrag/adffrag.py

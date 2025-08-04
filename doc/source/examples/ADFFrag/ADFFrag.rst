.. _ADFFragExample:

.. _adf-fragment-recipe:

ADF fragment job
================
.. currentmodule:: scm.plams.recipes.adffragment

In this module a dedicated job type for ADF fragment analysis is defined.
Such an analysis is performed on a molecular system divided into 2 fragments and consists of 3 separate ADF runs: one for each fragment and one for full system.

We define a new job type ``ADFFragmentJob`` by extending |MultiJob|.
The constructor (``__init__``) of this new job takes 2 more arguments (``fragment1`` and ``fragment2``) and one optional argument ``full_settings`` for additional input keywords that are used **only** in the full system calculation.

In the |prerun| method two fragment jobs and the full system job are created with the proper settings and molecules.
They are then added to the ``children`` list.

The dedicated |Results| subclass for ``ADFFragmentJob`` does not provide too much additional functionality.
It simply redirects the usual |AMSResults| methods to the results of the full system calculation.



.. literalinclude:: ../../../../recipes/adffragment.py

.. include:: ADFFrag.common_header.rst
.. include:: adffrag.ipynb.rst
.. include:: ADFFrag.common_footer.rst


API
-----------------------

.. autoclass:: ADFFragmentJob()
.. autoclass:: ADFFragmentResults()



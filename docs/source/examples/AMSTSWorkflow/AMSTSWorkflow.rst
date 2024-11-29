.. _AMSTSWorkflowExample:

AMS transition state workflow
===============================

Example illustrating a Diels-Alder addition between cyclopentadiene and acrylonitrile.

In this example, a preliminary biased UFF MD simulation is run to orient the two reactant molecules near the transition state.

Then a transition state search is performed.

The PESExploration LandscapeRefinement is used to find the two minima on either side of the transition state.

.. seealso::

    The :ref:`AMSPlumedMDExample` example is similar, but in that example the reaction happens during the MD simulation. That type of simulation can be performed when the potential is reactive.

    In the current example, the MD simulation is performed with UFF which is not reactive.

.. seealso::

    Example: :ref:`MoleculeSubstitutionExample`

.. include:: AMSTSWorkflow.common_header.rst
.. include:: diels_alder_addition.ipynb.rst
.. include:: AMSTSWorkflow.common_footer.rst
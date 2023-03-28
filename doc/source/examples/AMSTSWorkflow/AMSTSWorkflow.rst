.. _AMSTSWorkflow:

AMS transition state workflow
===============================

Example illustrating a Diels-Alder addition between cyclopentadiene and acrylonitrile.

In this example, a preliminary biased UFF MD simulation is run to orient the two reactant molecules near the transition state.

Then a transition state search is performed.

The PESExploration LandscapeRefinement is used to find the two minima on either side of the transition state.

.. seealso::

    The :ref:`AMSPlumedMD` example is similar, but in that example the reaction happens during the MD simulation. That type of simulation can be performed when the potential is reactive.

    In the current example, the MD simulation is performed with UFF which is not reactive.

.. seealso::

    Example: :ref:`MoleculeSubstitution`

To follow along, either

* Download :download:`diels_alder_addition.py` (run as ``$AMSBIN/amspython diels_alder_addition.py``).
* Download :download:`diels_alder_addition.ipynb` (see also: how to install `Jupyterlab <../../../Scripting/Python_Stack/Python_Stack.html#install-and-run-jupyter-lab-jupyter-notebooks>`__ in AMS)

.. include:: diels_alder_addition.rst.include

Complete Python code
----------------------------

.. literalinclude:: diels_alder_addition.py
    :language: python

.. _ASECalculatorExample:
.. _AMSCalculatorExample:

AMSCalculator: ASE geometry optimizer with AMS forces
===============================================================

**Note**: This example requires AMS2023 or later.

.. note::

    Follow this example only if you need to use ASE. If you just want to run a
    normal AMS single-point or geometry optimization, do not go through ASE but
    instead see the :ref:`GeoOptWaterExample` example.

Example illustrating how to use the :ref:`AMSCalculator`. The ``BFGS`` geometry optimizer from ASE is used together with AMS-calculated forces (negative gradients).

In this example, the AMS **driver** is replaced by ASE tools, that use the AMS
**engines** (ADF, BAND, DFTB, ForceField, ...) to calculate energies and forces.

.. seealso::

    In the :ref:`CustomASECalculatorExample` example, the AMS *engines* are
    instead replaced with external ASE calculators, that can be coupled to AMS
    *driver* tasks (GeometryOptimization, TransitionStateSearch,
    MolecularDynamics, …).

.. note::

    This example launches AMS in "AMSworker" mode. This means that AMS is only started at the beginning of the calculation.

To follow along, either

* Download :download:`AMSCalculatorWorkerMode.py` (run as ``$AMSBIN/amspython AMSCalculatorWorkerMode.py``).
* Download :download:`AMSCalculatorWorkerMode.ipynb` (see also: how to install `Jupyterlab <../../../Scripting/Python_Stack/Python_Stack.html#install-and-run-jupyter-lab-jupyter-notebooks>`__ in AMS)


.. include:: AMSCalculatorWorkerMode.rst.include


Complete Python code
----------------------------

.. literalinclude:: AMSCalculatorWorkerMode.py
	:language: python






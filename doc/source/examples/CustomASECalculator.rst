.. _CustomASECalculatorExample:

Engine ASE: AMS geometry optimizer with forces from any ASE calculator 
==========================================================================

**Note**: This example requires AMS2023 or later.

This example illustrates how set up an AMS calculation using any external ASE calculator as the engine.

In this example, the EMT calculator included with ASE is used.

.. note::

    In this example, the AMS *engine* is replaced by an ASE calculator.

    In the :ref:`AMSCalculatorExample`, :ref:`iPIExample`, and
    :ref:`SellaExample` examples, the AMS *driver* is replaced by external
    codes that use the AMS *engine* energies and forces.

**Example usage:** (:download:`Download CustomASECalculator.py <../../../examples/CustomASECalculator.py>`)

.. literalinclude:: ../../../examples/CustomASECalculator.py
	:language: python


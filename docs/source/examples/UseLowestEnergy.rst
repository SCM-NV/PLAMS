.. _UseLowestEnergy:

Hybrid engine: Use lowest energy
=====================================

**Note**: This example requires AMS2023 or later.

If you are unsure of the lowest energy spin state for a structure, the safest
way is to try multiple different possibilities. Starting from AMS2023, this can
be done with the `AMS Hybrid engine <../../Hybrid/EngineOptions.html>`__
together with the option DynamicFactors=UseLowestEnergy.

In this example, the hybrid engine is used to replay a PES Scan of a dissociating hydrogen peroxide molecule.
For each bond length, both singlet and triplet states are evaluated.
The resulting energy-vs-bond-length curve will only give the lowest energy for each bond length. This can
be useful if you want to import such a bond scan into for example ParAMS for parametrization.

You may also try to optimize the spin using the OptimizeSpinRound feature of
the ADF engine. This can be done in combination with an ElectronicTemperature.
The energies are then not as accurate.

**Example usage:** (:download:`Download UseLowestEnergy.py <../../../examples/UseLowestEnergy.py>`)

.. image:: /_static/pesplot.png

.. literalinclude:: ../../../examples/UseLowestEnergy.py
	:language: python


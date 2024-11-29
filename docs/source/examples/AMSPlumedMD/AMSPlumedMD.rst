.. _AMSPlumedMDExample:

AMS biased MD / PLUMED
=======================

Example illustrating how to run AMS MD with

* AMS Restraints

* AMS EngineAddon WallPotential

* Plumed restraints

During this simulation, the reaction H₂CO₃(g) → H₂O(g) + CO₂(g) is enforced by
a PLUMED MovingRestraint between one of the hydrogen and oxygen atoms.

The WallPotential keeps the two product molecules near each other.

.. include:: AMSPlumedMD.common_header.rst
.. include:: ams_plumed.ipynb.rst
.. include:: AMSPlumedMD.common_footer.rst
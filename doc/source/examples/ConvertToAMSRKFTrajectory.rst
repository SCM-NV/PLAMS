.. _ConvertToAMSRKFTrajectory:

Convert to ams.rkf trajectory with bond guessing
==========================================================

**Note**: This example requires AMS2023 or later.

To convert a trajectory from a different format to ams.rkf, you can for example
use the below script.

This may be useful if you for example have a trajectory in VASP or .xyz format
and want to use the AMS GUI to visualize the trajectory or run e.g.
ChemTraYzer2 reaction analysis.

The bond guessing is done by running a "replay" job with engine LennardJones.


**Example usage:** (:download:`Download ConvertToAMSRKFTrajectory.py <../../../examples/ConvertToAMSRKFTrajectory.py>`)

.. literalinclude:: ../../../examples/ConvertToAMSRKFTrajectory.py
	:language: python


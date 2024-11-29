.. _IRSpectrumFromMDH2ODimer:

IR spectrum H2O dimer from MD
==========================================

Run an MD simulation with GFN1-xTB, store the trajectory at every timestep, and compute the IR spectrum using the AMS trajectory analysis tools.
Results can be read in with the GUI, but are also written in human readable format. Executing the python script below takes a while (say 15 minutes - 1 hour) on a desktop computer.

**Example usage:** (:download:`Download get_irspectrum.py <../../../examples/get_irspectrum.py>`)

.. literalinclude:: ../../../examples/get_irspectrum.py
	:language: python


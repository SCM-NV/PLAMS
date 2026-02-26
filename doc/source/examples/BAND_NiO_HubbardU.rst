.. _Band_NiO_HubbardU:

BAND: NiO with DFT+U
=====================================

The hubbard method is a way to calculate band gaps for metal oxides for which a normal GGA fails to predict a band gap.

**Example usage:**

* Download :download:`BAND_NiO_HubbardU.py <../../../examples/BAND_NiO_HubbardU.py>` (run as ``$AMSBIN/amspython BAND_NiO_HubbardU.py``).

.. literalinclude:: ../../../examples/BAND_NiO_HubbardU.py
	:language: python

Results
-------

The following results appear upon execution.

.. parsed-literal::

   Top of valence band:       -6.67 eV
   Bottom of conduction band: -5.41 eV
   Band gap:                  1.25 eV

A folder called ``NiO`` is also created containing other files.

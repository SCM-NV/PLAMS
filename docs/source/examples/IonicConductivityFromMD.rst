.. _IonicConductivityFromMD:

Ionic conductivity from ams.rkf trajectory
==========================================

First, an ams.rkf trajectory with ions in the system needs to be calculated, as is done, for example, in the
`ionic conductivity Tutorial <../../Tutorials/MolecularDynamicsAndMonteCarlo/IonicConductivity.html>`__.
Next, the ions in the system and their respective charges need to be identified.

**Example usage:** (:download:`Download get_charged_ions.py <../../../examples/get_charged_ions.py>`)

.. code-block:: none

   $AMSBIN/amspython get_charged_ions.py /path/to/ams.rkf


.. literalinclude:: ../../../examples/get_charged_ions.py
	:language: python

Running the above script results in a small input file (charges.in) containing the ion information. This input file can be used in the next script, to compute the ionic conductivity.

**Example usage:** (:download:`Download get_ionic_conductivity.py <../../../examples/get_ionic_conductivity.py>`)

.. code-block:: none

   $AMSBIN/amspython get_ionic_conductivity.py /path/to/ams.rkf <charges.in

.. literalinclude:: ../../../examples/get_ionic_conductivity.py
	:language: python


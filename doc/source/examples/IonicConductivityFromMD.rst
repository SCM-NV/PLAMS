.. _IonicConductivityFromMD:

Ionic conductivity from ams.rkf trajectory
==========================================

First, an ams.rkf trajectory with ions in the system needs to be calculated, as is done, for example, in the
`Ionic Conductivity Tutorial <../../Tutorials/MolecularDynamicsAndMonteCarlo/IonicConductivity.html>`__.
Next, the ionic conductivity for the ions in the solution can be computed. The script below first identifies the ions in the molecular system and guesses their charges (for metals, it uses a simple charge equilibration scheme). It then runs the AMS trajectory analysis tool to compute the ionic conductivity for these ions.

**Example usage:** (:download:`Download get_ionic_conductivity.py <../../../examples/IonicConductivity/get_ionic_conductivity.py>`)

.. code-block:: bash

   $AMSBIN/amspython get_ionic_conductivity.py /path/to/ams.rkf

To create an example ``ams.rkf`` file, you can use this PLAMS script (:download:`Download NaClwater.py <../../../examples/IonicConductivity/NaClwater.py>`). This script is designed for efficiency. For a more reliable setup, follow the `Ionic Conductivity Tutorial <../../Tutorials/MolecularDynamicsAndMonteCarlo/IonicConductivity.html>`__.

.. literalinclude:: ../../../examples/IonicConductivity/get_ionic_conductivity.py
   :language: python

**Results**

If you have used the ``ams.rkf`` generated from running ``$AMSBIN/amspython NaClwater.py``,
the output from ``get_ionic_conductivity.py`` should look something like:

.. parsed-literal::
     Average temperature 298.166533 K
     Ion        N     Charge
     Na        5    1.00000
     Cl        5   -1.00000
     [20.02|15:06:43] JOB plamsjob STARTED
     [20.02|15:06:43] JOB plamsjob RUNNING
     [20.02|15:06:43] JOB plamsjob FINISHED
     [20.02|15:06:43] JOB plamsjob SUCCESSFUL
     Na 2.4739268145933062e-09 m2/s
     [20.02|15:06:43] JOB plamsjob STARTED
     [20.02|15:06:43] Renaming job plamsjob to plamsjob.002
     [20.02|15:06:43] JOB plamsjob.002 RUNNING
     [20.02|15:06:43] JOB plamsjob.002 FINISHED
     [20.02|15:06:43] JOB plamsjob.002 SUCCESSFUL
     Cl 3.2624546292849828e-09 m2/s
     Ionic conductivity:     2.6075090545e+01 Siemens/m

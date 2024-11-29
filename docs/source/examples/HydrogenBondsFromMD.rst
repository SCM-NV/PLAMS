.. _HydrogenBondsFromMD:

Hydrogen bonds from MD
======================

SCM Python tools define hydrogen bonds as X1-H---X2 patterns with an maximum X1-X2 distance of 3.2 Angstrom, a maximum H-X distance of 2.5 Angstrom, and maximum angle of 36 degrees. Here, we present a script that extracts the number of hydrogen bonds for a selected set of atoms along a trajectory, and creates a histogram. The atoms are selected via an input file named ``indices.txt``.

**Example usage:** 
(:download:`Download get_hbonds.py <../../../examples/get_hbonds.py>`)

The script ``get_hbonds.py`` accepts the name of an MD trajectory RKF file, and the name of a text file with the indices (and/or element names) of the atoms for which the number of hydrogen bonds are to be computed. The script can be called from the command line as follows.

.. code-block:: bash

   amspython get_hbonds.py path/to/ams.rkf path/to/indices.txt

In its simplest form, the input file indices.txt will contain only one or more element names.

.. code-block:: bash

   O C

A more complicated file will contain specific atom indices. At the bottom of this page an example script (``get_water_indices.py``) is presented that creates such a file for all oxygen atoms that belong to a water molecule. 

The script ``get_hbonds.py``.

.. literalinclude:: ../../../examples/get_hbonds.py
	:language: python

(:download:`Download get_water_indices.py <../../../examples/get_water_indices.py>`)

The script ``get_water_indices`` can create the file ``indices.txt`` that is required by the ``get_hbonds.py`` script, for the specific case where the indices of all water oxygen atoms are required. The script accepts the name of the MD trajectory RKF file.

.. literalinclude:: ../../../examples/get_water_indices.py
	:language: python


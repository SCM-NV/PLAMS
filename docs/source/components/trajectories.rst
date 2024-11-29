Trajectories
-------------------------

The main output of a molecular simulation is often a series of conformations of a molecular system,
as successively produced by the software over time.

Example applications that produce such a trajectory are molecular dynamics (MD) simulations, 
and geometry optimizations.

The conformations consist of coordinates for all the atoms in the system, and - depending on the format -
information on its periodicity, atomic elements and the location of covalent bonds.

Some file formats include

* XYZ format (``xyz``), which is a human readable format that stores elements and coordinates for every conformation.

* Compact binary DCD format (``dcd``), which stores only coordinates and periodic lattice vectors,

* Binary RKF format (``rkf``) which is the native format of the Amsterdam Modeling Suite, and is able to store any desired information on the system.

* Binary .traj (``traj``) format, which is the native format used by the Atomic Simulation Environment

For **trajectory conversion tools**, see :ref:`FileFormatConversionTools`

PLAMS provides the possibility to analyze and alter single conformations of a molecular system through the :ref:`Molecule` object .

The ``trajectories`` module provides the possibility to extract such a conformation from a trajectory file 
and load it into a :ref:`Molecule` object,
as well as write (altered) conformations to a new trajectory file of one in the three formats described above.



.. toctree::

    xyz
    rkf
    dcd
    trajectory

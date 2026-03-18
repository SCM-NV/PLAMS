.. _ASEInterface:

ASE interface
~~~~~~~~~~~~~

The Atomic Simulation Environment (`ASE <https://wiki.fysik.dtu.dk/ase>`_) is a set of tools and Python modules for setting up, manipulating, running, visualizing, and analyzing atomistic simulations.
The PLAMS interface to ASE is limited to handling |Molecule| objects.
It features access to the ``ase.io`` module for reading/writing |Molecule| objects and two functions that translate PLAMS |Molecule| objects into ASE ``Atoms`` objects and vice versa (see :meth:`~scm.plams.mol.molecule.Molecule.readase` and :meth:`~scm.plams.mol.molecule.Molecule.writease`).

Please refer to the ASE documentation to see what can be done with ASE ``Atoms`` and its I/O module.

.. automodule:: scm.plams.interfaces.molecule.ase

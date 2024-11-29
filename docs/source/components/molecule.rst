.. _Molecule:

Molecule
-------------------------

In this chapter the PLAMS module responsible for handling molecular geometries is presented.
Information about atomic coordinates can be read from (or written to) files of various types: ``xyz``, ``pdb``, ``mol`` or ``mol2``.
PLAMS not only extracts the relevant data from those files, but also tries to "understand" the structure of the underlying molecule in terms of atoms and bonds between them, allowing a variety of simple operations like moving or rotating some parts of the molecule, splitting it into multiple parts, merging two molecules etc.
For importing coordinates in other formats or performing more advanced operations you can use :ref:`ASEInterface` and :ref:`RDKitmol`.




.. toctree::

    mol_api
    atombond
    mol_rdkit
    mol_ase
    mol_packmol

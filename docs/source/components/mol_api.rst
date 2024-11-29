Molecule
~~~~~~~~~~~~~~~~~~~~~~~~~


.. autoclass :: scm.plams.mol.molecule.Molecule
    :exclude-members: __weakref__, __copy__, reorder


Atom labeling
+++++++++++++

This subsection describes API of ``identify`` module, which is used to assign unique names to atoms in a molecule.
Unique atom names are used in |Molecule| labeling (:meth:`~scm.plams.mol.molecule.Molecule.label`).
.. and in the method restoring the order of atoms (:meth:`~scm.plams.mol.molecule.Molecule.reorder`).
All the functions below, apart from :func:`~scm.plams.mol.identify.label_atoms`, are for internal use and they are not visible in the main PLAMS namespace.


.. currentmodule:: scm.plams.mol.identify

.. autofunction:: label_atoms

..
    .. autofunction:: find_permutation
.. autofunction:: molecule_name
.. autofunction:: initialize
.. autofunction:: iterate
.. autofunction:: clear
.. autofunction:: new_name
.. autofunction:: knock
.. autofunction:: twist
.. autofunction:: bend
.. autofunction:: unique_atoms

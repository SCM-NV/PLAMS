.. _PackmolInterface:

Packmol interface
~~~~~~~~~~~~~~~~~~

Packmol (`Packmol website <https://m3g.github.io/packmol/download.shtml>`__) is a program for creating liquid or gas mixtures. The PLAMS interface only supports

* uniform mixtures
* solid/liquid interfaces
* packing inside voids in a crystal
* microsolvation

The following functions eixst:

* ``packmol`` (for fluids with 1 or more components)
* ``packmol_on_slab`` (for solid/liquid or solid/gas interfaces with 1 or more components in the fluid)
* ``packmol_in_void`` (for packing molecules inside crystal voids)
* ``packmol_microsolvation`` (for microsolvation of a solute with a solvent)

See the :ref:`Packmol example <PackMolExample>` for all the ways these functions can be used.

The above functions accept an ``executable`` argument, which should
contain the path to the packmol program. If it is not specified, the path to
the packmol program included with the Amsterdam Modeling Suite will be used.

.. currentmodule:: scm.plams.interfaces.molecule.packmol

.. autofunction:: packmol

.. autofunction:: packmol_on_slab

.. autofunction:: packmol_in_void

.. autofunction:: packmol_microsolvation

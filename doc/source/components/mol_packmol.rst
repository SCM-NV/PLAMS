.. _PackmolInterface:

Packmol interface
~~~~~~~~~~~~~~~~~~

Packmol (`Packmol website <https://m3g.github.io/packmol/download.shtml>`__) is a program for creating liquid or gas mixtures. The PLAMS interface only supports

* uniform mixtures
* solid/liquid interfaces
* packing inside voids in a crystal
* microsolvation

The following functions exist:

* ``packmol`` (for fluids with 1 or more components)
* ``packmol_around`` (for fluids with 1 or more components, used to pack around an existing system in AMS2025+)
* ``packmol_microsolvation`` (for microsolvation of a solute with a solvent)

In AMS2025, we recommend using ``packmol_around`` over ``packmol_on_slab`` and ``packmol_in_void``:

* ``packmol_on_slab`` (deprecated, for solid/liquid or solid/gas interfaces with 1 or more components in the fluid)
* ``packmol_in_void`` (deprecated, for packing molecules inside crystal voids)

See the `Packmol example <../../PythonExamples/pack-mol-example/index.html>`__ for all the ways these functions can be used.

The above functions accept an ``executable`` argument, which should
contain the path to the packmol program. If it is not specified, the path to
the packmol program included with the Amsterdam Modeling Suite will be used.

.. currentmodule:: scm.plams.interfaces.molecule.packmol

.. autofunction:: packmol

.. autofunction:: packmol_around

.. autofunction:: packmol_on_slab

.. autofunction:: packmol_in_void

.. autofunction:: packmol_microsolvation

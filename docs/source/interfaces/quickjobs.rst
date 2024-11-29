.. _Quickjobs:

Quick jobs
----------------------

.. currentmodule:: scm.plams.interfaces.adfsuite.quickjobs

Quick jobs that take in a structure and return a modified structure.
For example, the ``preoptimize()`` function lets you quickly optimize a
molecule without storing any results on disk. This mimics the preoptimize
function of the AMS GUI.

You can specify the ``model`` string as a shorthand for engine settings:

* 'UFF': the UFF force field
* 'GFNFF': the GFNFF force field
* 'ANI-2x': The ANI-2X machine learning potential

Alternatively, you can specify ``settings`` which will override the model.

.. autofunction:: preoptimize

.. autofunction:: refine_density

.. autofunction:: refine_lattice

.. autofunction:: shakemd

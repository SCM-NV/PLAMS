.. _NumHessExample:

Numerical Hessian
=================

This module implements a simple scheme for calculating a numerical Hessian matrix.
We define a new job type ``NumHessJob`` by extending |MultiJob|.
The constructor (``__init__``) of this new job accepts several new arguments and simply stores them.
These new arguments define the initial |Molecule|, the type of job used for single point calculations (``jobtype``), the size and unit of displacement step and the way of extracting gradients from single point results.

Then the |prerun| method takes the given |Molecule| and creates multiple copies of it, each one with one atom displaced along one axis.
For each of these molecules an instance of single point job is created and stored in the ``children`` dictionary.
Settings of ``NumHessJob`` are directly passed to children jobs, so creating a ``NumHessJob`` strongly resembles creating a regular single point job.

The dedicated |Results| subclass for ``NumHessJob`` takes care of extracting the Hessian from results of all single point jobs.
The returned Hessian can optionally be mass weighted.

The source code of the whole module with both aforementioned classes:

.. literalinclude:: ../../../../src/scm/plams/recipes/numhess.py
   :language: python

.. include:: NumHess.common_header.rst
.. include:: NumHess.ipynb.rst
.. include:: NumHess.common_footer.rst

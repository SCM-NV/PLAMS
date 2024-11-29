Numerical gradients
--------------------

This module implements a simple numerical differentiation scheme with respect to molecular coordinates.
We define a new job type ``NumGradJob`` by extending |MultiJob|.
The constructor (``__init__``) of this new job accepts several new arguments and simply stores them.
These new arguments define the initial |Molecule|, the type of job used for single point calculations (``jobtype``), the accuracy of the numerical differentiation scheme (``npoints``) and size and unit of displacement step.

Then the |prerun| method takes the given |Molecule| and creates multiple copies of it, each one with one atom displaced along one axis.
For each of these molecules an instance of single point job is created and stored in the ``children`` dictionary.
Settings of ``NumGradJob`` are directly passed to children jobs, so creating a ``NumGradJob`` strongly resembles creating a regular single point job.

The dedicated |Results| subclass for ``NumGradJob`` takes care of extracting the gradient from the results of all single point jobs.
Any function that takes results of a single point job and returns a single number can be used with the ``get_gradient`` method defined in  ``NumGradResults``

The source code of the whole module with both abovementioned classes:

.. literalinclude:: ../../../recipes/numgrad.py

An example usage:

.. literalinclude:: ../../../examples/numgrad.py


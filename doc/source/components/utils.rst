Utilities
-------------------------

Presented here is a small set of useful utility tools that can come in handy in various contexts in your scripts.
They are simple, standalone objects always present in the main namespace.

.. contents:: :local:

What is characteristic for the |PeriodicTable| and |Units| classes described below is that they are meant to be used in a bit different way than all other PLAMS classes.
Usually one takes a class (like |DiracJob|), creates an instance of it (``myjob = DiracJob(...)``) and executes some of its methods (``r = myjob.run()``).
In contrast, utility classes are designed in a way similar to the so-called singleton design pattern.
That means it is not possible to create instances of these classes.
The class itself serves for "one and only instance" and all methods should be called using the class as the calling object::

    >>> x = PeriodicTable()
    PTError: Instances of PeriodicTable cannot be created
    >>> s = PeriodicTable.get_symbol(20)
    >>> print(s)
    Ca

Periodic Table
~~~~~~~~~~~~~~~~~~~~~~~~~

Import path::

    scm.plams.tools.periodic_table

.. autoclass:: scm.plams.tools.periodic_table.PeriodicTable
    :exclude-members: __weakref__

Units
~~~~~~~~~~~~~~~~~~~~~~~~~

Import path::

    scm.plams.tools.units

.. autoclass:: scm.plams.tools.units.Units
    :exclude-members: __weakref__


Geometry tools
~~~~~~~~~~~~~~~~~~~~~~~~~

A small module with simple functions related to 3D geometry operations.

Import path::

    scm.plams.tools.geometry

.. automodule:: scm.plams.tools.geometry

.. _FileFormatConversionTools:

File format conversion tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A small module for converting VASP output to AMS-like output, and for converting ASE .traj trajectory files to the .rkf format.

Import path::

    scm.plams.tools.converters

.. automodule:: scm.plams.tools.converters

.. _ReactionEnergies:

Reaction energies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*New in AMS2026*: The ``balance`` function is new in AMS2026. For usage, see
the :ref:`BalanceReactionEquationsExample` example.

Import paths::

    scm.plams.tools.reaction
    scm.plams.tools.reaction_energies

.. autofunction:: scm.plams.tools.reaction.balance

.. autoclass:: scm.plams.tools.reaction.ReactionEquation


Older functions:

.. automodule:: scm.plams.tools.reaction_energies

.. _PlottingTools:

Plotting tools
~~~~~~~~~~~~~~~~~~

.. seealso::

    * :ref:`BandStructureExample`

Tools for creating plots with matplotlib.

Import path::

    scm.plams.tools.plot

The :mod:`scm.plams.tools.plot` module also contains small reusable helpers for
common analysis tasks. For example, ``linear_fit_extrapolate_to_0`` performs a
linear regression and returns the fitted line extended to ``x = 0``.

Example::

    from scm.plams.tools.plot import linear_fit_extrapolate_to_0

    fit_x, fit_y, slope, intercept = linear_fit_extrapolate_to_0(
        [1.0, 2.0, 3.0],
        [3.0, 5.0, 7.0],
    )

.. automodule:: scm.plams.tools.plot

.. _PostprocessResults:

Postprocess results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tools for postprocessing the results.

Import path::

    scm.plams.tools.postprocess_results

The :mod:`scm.plams.tools.postprocess_results` module contains helpers such as
``moving_average`` for smoothing paired ``x``/``y`` data and ``broaden_results``
for constructing broadened spectra.

Example::

    from scm.plams.tools.postprocess_results import moving_average

    avg_x, avg_y = moving_average([1.0, 2.0, 3.0], [3.0, 5.0, 7.0], window=2)

.. automodule:: scm.plams.tools.postprocess_results

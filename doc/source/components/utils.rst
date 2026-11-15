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
the `Balance reaction equations example <../../PythonExamples/balance-reaction-equations/index.html>`__.

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

    * `Band structure with BAND <../../PythonExamples/bandstructure-band/index.html>`__

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

Reaction networks
^^^^^^^^^^^^^^^^^

:func:`~scm.plams.tools.plot.plot_reaction_network` plots minima as circles and
transition states as squares, with directed ``reactant -> TS -> product``
connections. Each transition state remains a separate node, including when
several transition states connect the same minima. Arrowheads are placed at
connector midpoints. The stored orientation is preserved; it does not imply
irreversible kinetics.

Small networks use connectivity-based placement and state IDs. Larger networks
use a deterministic formula-grouped layout. Marker areas adapt to network size;
visible labels set an additional minimum marker size. When the function creates
a small-network figure, it also expands the canvas to maintain the requested
gap between markers. Supplying an existing axes leaves its figure size intact.

For example, read an existing energy landscape and save an energy-colored network::

    from scm.plams import AMSJob, plot_reaction_network

    job = AMSJob.load_external("example.results")
    landscape = job.results.get_energy_landscape()
    ax = plot_reaction_network(landscape, color_by="energy", cmap="jet")
    ax.figure.savefig("reaction_network.png", dpi=180)

The default ``color_by="formula"`` groups nonsurface fragments by formula and
marks adsorbed fragments with ``*``. Isomers with identical fragment formulas
share a color; their nodes are not collapsed. Use ``surface_only=True`` with a
named ``surface_region`` to retain only connected, fully surface-bound states.
Bond guessing does not modify the original molecules. It uses ``scm.base`` when
available and native PLAMS for nonperiodic molecules otherwise; periodic
molecules require ``scm.base``.

Energy coloring uses a shared scale for minima and transition states, relative
to the lowest plotted state. Input energies and ``energy_reference`` are in
hartree; ``unit`` controls the colorbar and ``energy_limits`` units. To compare
colors across networks, pass the same ``energy_limits``. Changing color mode
does not change the layout; explicit ``positions`` can also be reused.

The function needs the optional Matplotlib and NetworkX packages. Computing the
tree layout also needs SciPy. Saving is controlled by the caller, and the
returned axes may be further customized. This diagram shows connectivity;
:func:`~scm.plams.tools.plot.plot_energy_landscape` provides an energy-level view.

Landscapes can also be read directly from exported files without running AMS::

    from scm.plams import AMSResults, plot_reaction_network

    medium = AMSResults.EnergyLandscape.from_ams_inp("EnergyLandscape_HCNO.inp")
    ax = plot_reaction_network(medium, color_by="energy", cmap="jet")

    large = AMSResults.EnergyLandscape.from_xyz("EnergyLandscape_MetOH_Cu111.xyz")
    ax = plot_reaction_network(large, surface_only=True, color_by="energy")

The AMS INP reader accepts optional ``energy``, ``count``,
``preFactReact`` and ``preFactProd`` System-header metadata.
Energies are in hartree and directional prefactors in inverse seconds. Older files
may omit energies; the reader retains these as ``NaN``. An optional ``energies``
mapping, keyed by System label, supplies or overrides energies in hartree.
The compact XYZ reader preserves energies, counts, atom regions, lattice vectors,
directional prefactors when present, and named transition-state connections.
An omitted ``ts`` flag denotes a minimum. Both readers resolve links after all
states are loaded, so forward references and arbitrary state labels are supported.

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

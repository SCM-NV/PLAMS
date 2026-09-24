COSMO-RS
--------

.. currentmodule:: scm.plams.interfaces.adfsuite.crs

COSMO-RS can be run from PLAMS using the |CRSJob| class and the corresponding |CRSResults|,
both respectively being subclasses of |SCMJob| and |SCMResults|.

There are three ways to define a CRS job:

.. code-block:: text

    Input representation
          |
          +-- Settings
          |     |
          |     +-- job = CRSJob(settings=settings)
          |
          +-- typed CRS model
          |     |
          |     +-- job = CRSJob(settings=crs)
          |
          +-- CRS input builder
                |
                +-- job = builder.to_job()
                |       [equivalent to CRSJob(settings=builder.to_settings())]
                |
                +-- settings = builder.to_settings()
                |       +-- job = CRSJob(settings=settings)
                |
                +-- crs = builder.to_inputs()
                        +-- job = CRSJob(settings=crs)

For most workflows, use :meth:`CRSJob.input_builder` followed by
``builder.to_job()``. Use ``to_settings()`` or ``to_inputs()`` when you need
to inspect or modify the generated input before creating the job.

.. note:: There is also `a tutorial showing full code examples
   <../../COSMO-RS/PLAMS_COSMO-RS_scripting.html>`__ available in the COSMO-RS
   documentation. There are several templates available that can easily be
   customized for other problem types, workflows, etc.

Input builders
~~~~~~~~~~~~~~

Use :meth:`CRSJob.input_builder` to create a property-specific builder. The
builder provides the input keys, compound roles, and calculation modes
supported by the selected property type.

Basic example
^^^^^^^^^^^^^

The following example creates and runs a pure sigma-profile calculation
with :meth:`CRSJob.input_builder`:

.. code-block:: python

    from scm.plams import CRSJob

    builder = CRSJob.input_builder("PURESIGMAPROFILE", method="COSMO-RS")
    builder.nprofile = 50
    builder.sigmamax = 0.025
    builder.add_compound_from_adfcrs_database("Water.coskf")

    job = builder.to_job()
    results = job.run()

This builder generates the following COSMO-RS input:

.. code-block:: none

    method COSMO-RS

    compound /path/to/file.coskf
    end

    property puresigmaprofile
        nprofile 50
        sigmamax 0.025
    end

Configuring the input builder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The input builder is configured through the property type, method, mode, and
compounds.

.. code-block:: python

    builder = CRSJob.input_builder("SOLUBILITY")
    builder.method = "COSMO-RS"
    builder.mode = "solid"
    builder.temperature = 250
    builder.add_solvent_from_adfcrs_database("Water.coskf", frac1=1.0)
    builder.add_solute_from_adfcrs_database("Benzene.coskf", meltingpoint=278.7, hfusion=2.37)


Inspecting and setting builder input
=====================================

Use ``builder.describe()`` to inspect the available input keys, their types,
units, and descriptions, as well as supported modes and required keys:

.. code-block:: python

    for line in builder.describe(include_values=True):
        print(line)

The output is similar to:

.. code-block:: text

    SOLUBILITY: Solubility of solutes in a solvent mixture or under gas-pressure conditions.
    method [top_level]: COSMO-RS
    densitysolvent [property] float [kg/L]: Density of the solvent. value: None
    massfraction [top_level] bool: Use mass fractions; by default, fractions are interpreted as molar fractions. value: None
    pressure [top_level] float_list [bar]: Pressure value: None
    temperature [top_level] float_list [Kelvin]: Temperature value: 250
    mode: solid, liquid, gas
    compound_keys [compound]: frac1, density, meltingpoint, hfusion, cpfusion, pvap, tvap, vp_equation, vp_params
    required_keys [top_level]: temperature
    required_keys [compound: solvent]: frac1
    required_keys [compound: solute]: meltingpoint, hfusion
    hint [densitysolvent, density]: Used for molar-volume estimates, volume-based solubility, and Henry's-law results; falls back to COSMO volume estimates.
    hint [meltingpoint, hfusion, cpfusion]: Used for solid-solute fusion corrections; the heat capacity is optional.
    hint [pvap, tvap, vp_equation, vp_params]: Used for gas-phase pseudochemical potential corrections in VLE, flash, and Henry's-law calculations.

The labels indicate where inputs belong and which values are required:

* ``[property]`` identifies a key in the ``PROPERTY`` block of the CRS input.
* ``[top_level]`` identifies a key at the top level of the CRS input.
* ``mode`` lists the calculation modes supported by the selected property.
* ``compound_keys [compound]`` lists the compound keys that affect the result for the selected property type.
* ``required_keys`` lists mandatory keys for the current property, mode, and
  compound role. Other listed keys are optional.
* ``hint`` gives additional guidance for selected inputs, such as when they are used or how to choose their values.

Use ``include_values=True`` to also show the current builder values.
A value of ``None`` indicates that the input has not been set on the builder.

Use ``builder.describe(include_compound_details=True)`` to also show the
types, units, and descriptions of compound keys.

Set values directly on the builder. For gas-phase solubility, set the
temperature, specify the solute partial pressure in bar, optionally provide the
solvent density in kg/L, and select gas mode:

.. code-block:: python

    builder.mode = "gas"
    builder.temperature = 298.15
    builder.pressure = 1.01325
    builder.densitysolvent = 1.0

Changing the mode may change the required keys. Call ``builder.describe()``
again to inspect the updated requirements.


Calculation modes
==================

The following table summarizes the supported modes and their defaults.
Properties not listed here do not support a calculation mode.

.. list-table::
    :header-rows: 1
    :widths: 35 65

    * - Property type
      - Available modes
    * - ``SOLUBILITY``
      - ``solid`` (default), ``liquid``, ``gas``
    * - ``PURESOLUBILITY``
      - ``solid`` (default), ``liquid``, ``gas``
    * - ``BINMIXCOEF``
      - ``isotherm`` (default), ``isobar``, ``flashpoint``
    * - ``TERNARYMIX``
      - ``isotherm`` (default), ``isobar``, ``flashpoint``
    * - ``COMPOSITIONLINE``
      - ``isotherm`` (default), ``isobar``, ``flashpoint``


Adding compounds
================

The available compound roles and the number of compounds required depend on
the selected property type.

The builder checks compound roles and count limits when compounds are added,
and validates required counts and data when generating input.

Pass ``include_compounds=False`` to ``to_settings()`` or ``to_inputs()``
to omit compounds and skip compound validation during conversion.

For mixture properties, use ``add_compound()`` with a COSKF file path or
``add_compound_from_adfcrs_database()`` with an ADFCRS database filename:

.. code-block:: python

    from scm.plams import CRSJob

    builder = CRSJob.input_builder("TERNARYMIX", temperature=298.15)

    builder.add_compound_from_adfcrs_database("Water.coskf")
    builder.add_compound_from_adfcrs_database("Ethanol.coskf")
    builder.add_compound_from_adfcrs_database("Benzene.coskf")

    settings = builder.to_settings()

For solvent and solute roles, use ``add_solvent()`` and ``add_solute()`` or
``add_solvent_from_adfcrs_database`` and ``add_solute_from_adfcrs_database``:

.. code-block:: python

    builder = CRSJob.input_builder("ACTIVITYCOEF", temperature=298.15)
    builder.add_solvent_from_adfcrs_database("Water.coskf", frac1=1.0)
    builder.add_solute_from_adfcrs_database("Benzene.coskf")
    crs = builder.to_inputs()


Supported methods
==================

The builder uses ``COSMO-RS`` by default. Use ``builder.method`` to select
another method. To list the supported methods:

.. code-block:: python

    print(CRSJob.methods())

Each method uses its default parameter set. To use a different parameter set,
apply a preset to the generated input as described in
:ref:`crs_method_parameter_presets`.


Creating and running a job
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use ``to_job()`` when no additional settings are required before running the
job:

.. code-block:: python

    job = builder.to_job()
    results = job.run()


Converting to Settings
^^^^^^^^^^^^^^^^^^^^^^^

Use ``to_settings()`` when you want to inspect or modify the generated PLAMS
|Settings| object before creating the job:

.. code-block:: python

    settings = builder.to_settings()
    job = CRSJob(settings=settings)
    results = job.run()


Converting to typed CRS inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :meth:`CRSInputBuilder.to_inputs` for advanced workflows. It converts the
builder state to a typed :class:`scm.inputs.CRS` model, which gives you access
to CRS input options that are not exposed directly by the property-specific
builder.

For example, the following code changes an LLE convergence tolerance and enables
debug output:

.. code-block:: python

    from scm.plams import CRSJob

    builder = CRSJob.input_builder("LLE", temperature=298.15)
    builder.add_compound_from_adfcrs_database("Water.coskf", frac1=0.33)
    builder.add_compound_from_adfcrs_database("Ethanol.coskf", frac1=0.33)
    builder.add_compound_from_adfcrs_database("Benzene.coskf", frac1=0.34)

    crs = builder.to_inputs()
    crs.TECHNICAL.LLE.eps_g = 1.0e-5
    crs.TECHNICAL.LLE.debug = True

    job = CRSJob(settings=crs)
    results = job.run()

The object returned by ``to_inputs()`` is an independent
:class:`scm.inputs.CRS` model. Changes to this object do not update the builder.
Create the job from the modified model, as shown above, instead of calling
``builder.to_job()``.

For a general introduction to typed input models, see the
"AMS input models" Python example in the AMS documentation.


.. _crs_method_parameter_presets:

Applying method parameter presets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each CRS method uses a default parameter set. For example, ``COSMO-RS`` uses
``adf-combi2005``, and ``COSMOSAC2013`` uses ``2013-adf-xiong``.
To use a different parameter set, first inspect the available presets:

.. code-block:: python

    from scm.plams import CRSJob

    for method, parameter_sets in CRSJob.get_parameter_set_options().items():
        print(method, parameter_sets)

Example output:

.. code-block:: text

    COSMO-RS ('adf-combi2005', 'adf-combi1998', 'adf-lei-2018', 'klamt', 'mopac-pm6')
    COSMOSAC2013 ('2013-adf-xiong', '2013-adf-pure-xiong')
    COSMOSACDHB ('dhb-adf-chen',)
    COSMOSACDHB-MESP ('dhb-adf-mesp',)
    COSMOSAC2016 ('2016-adf-chen',)
    COSMOSAC2010 ('2010-hsieh',)
    COSMOSAC2007 ('2007-wang',)


Select a ``parameter_set`` from the returned tuple and apply it to the CRS
input:

.. code-block:: python

    crs = builder.to_inputs()
    CRSJob.apply_parameter_set_to_inputs(crs, "adf-lei-2018")

    # Optionally override individual parameters.
    crs.CRSPARAMETERS.chb = 9000.0

    job = CRSJob(settings=crs)

You can apply the same parameter set to PLAMS settings:

.. code-block:: python

    settings = builder.to_settings()
    CRSJob.apply_parameter_set_to_settings(settings, "adf-lei-2018")

    settings.input.crsparameters.chb = 9000.0

    job = CRSJob(settings=settings)

Both functions modify and return the supplied object. They set the method
and replace its parameter blocks, removing parameter blocks that are absent
from the selected preset. Other settings, including compounds and the
selected property, are preserved.


ADF and CRSJob
~~~~~~~~~~~~~~

A workflow is presented in the `ADF and COSMO-RS workflow Python example <../../PythonExamples/ams-crs-workflow/index.html>`__.
In this workflow, we follow the usual procedure of generating the inputs required to run COSMO-RS and COSMO-SAC calculations.

.. _parameters: ../../COSMO-RS/COSMO-RS_and_COSMO-SAC_parameters.html


Data analysis
~~~~~~~~~~~~~

Use :meth:`CRSResults.get_results` to read the results section of a ``.crskf``
file as a dictionary. If no section is supplied, PLAMS uses the property section
from the most recent calculation.

.. code-block:: python

    from scm.plams import CRSJob

    builder = CRSJob.input_builder("PURESIGMAPROFILE")
    builder.nprofile = 50
    builder.sigmamax = 0.025
    builder.add_compound_from_adfcrs_database("Water.coskf")

    results = builder.to_job().run()
    data = results.get_results()

    print(data["section"])
    print(data["filename"])
    print(data["chdval"])
    print(data["profil"])

A complete overview of all available sections and keys can be printed from the
KF file skeleton:

.. code-block:: python

    print(results._kf.get_skeleton())


API
~~~

.. autoclass:: CRSJob
    :no-private-members:
    :no-special-members:
    :exclude-members: input_builder

    .. py:staticmethod:: CRSJob.input_builder(property_type, *, method="COSMO-RS", mode=None, **kwargs)

        Return a property-specific CRS input builder.

        The ``property_type`` argument selects the builder class. The returned
        builder validates the keys, modes, and compound roles supported by that
        property type.

.. autoclass:: CRSResults
    :members:
        get_multispecies_dist,
        get_structure_energy,
        get_activity_coefficient,
        get_energy,
        get_sigma_profile,
        get_sigma_potential,
        get_results,
        get_prop_names
    :no-private-members:
    :no-special-members:

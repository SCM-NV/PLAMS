COSMO-RS
--------

.. currentmodule:: scm.plams.interfaces.adfsuite.crs

COSMO-RS can be run from PLAMS using the |CRSJob| class and the corresponding |CRSResults|,
both respectively being subclasses of |SCMJob| and |SCMResults|.

.. note:: There is also `a tutorial showing full code examples <../../COSMO-RS/PLAMS_COSMO-RS_scripting.html>`__ available in the COSMO-RS documentation. There are several templates available that can easily be customized for other problem types, workflows, etc.

Input builders
~~~~~~~~~~~~~~

The recommended way to prepare COSMO-RS input from PLAMS is
:meth:`CRSJob.input_builder`. This method returns a property-specific input
builder. The builder exposes the input keys, compound roles, and calculation
modes supported by the selected property type.

For most workflows, create a builder and call ``to_job()`` to create a
:class:`CRSJob` directly. You can also convert the builder to a |Settings|
object or to a typed :class:`scm.inputs.CRS` model when you need to inspect or
modify the generated input.

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
        frac1 1.0
    end

    property puresigmaprofile
        nprofile 50
        sigmamax 0.025
    end


Supported methods
^^^^^^^^^^^^^^^^^

The supported COSMO-RS methods can be inspected with:

.. code-block:: python

    print(CRSJob.methods())


Calculation modes
^^^^^^^^^^^^^^^^^

Some property types support calculation modes. For example, solubility can be
set up for a gas-phase solute using ``mode="gas"``:

.. code-block:: python

    builder = CRSJob.input_builder(
        "SOLUBILITY",
        mode="gas",
        temperature=298.15,
        pressure=1.01325
    )
    builder.add_solvent_from_adfcrs_database("Water.coskf", frac1=1.0)
    builder.add_solute_from_adfcrs_database("Benzene.coskf")


Inspecting builder input
^^^^^^^^^^^^^^^^^^^^^^^^

Use ``builder.describe()`` to inspect the input keys accepted by a builder:

.. code-block:: python

    for line in builder.describe(include_values=True):
        print(line)


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


Working with multiple compounds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Many COSMO-RS property types require more than one compound. With
:meth:`CRSJob.input_builder`, compounds are added with the role-specific methods
supported by the selected property type. The builder checks that compounds are
added with roles and counts supported by that property type.

For mixture properties, use ``add_compound()``:

.. code-block:: python

    from scm.plams import CRSJob

    builder = CRSJob.input_builder("TERNARYMIX", temperature=298.15)

    builder.add_compound_from_adfcrs_database("Water.coskf", frac1=0.33)
    builder.add_compound_from_adfcrs_database("Ethanol.coskf", frac1=0.33)
    builder.add_compound_from_adfcrs_database("Benzene.coskf", frac1=0.34)

    job = builder.to_job()
    results = job.run()

For solvent and solute roles, use ``add_solvent()`` and ``add_solute()``:

.. code-block:: python

    builder = CRSJob.input_builder("SOLUBILITY", temperature=298.15)
    builder.add_solvent_from_adfcrs_database("Water.coskf", frac1=1.0)
    builder.add_solute_from_adfcrs_database(
        "Benzene.coskf",
        meltingpoint=278.7,
        hfusion=2.37,
    )

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

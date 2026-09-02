COSMO-RS
--------

.. currentmodule:: scm.plams.interfaces.adfsuite.crs

COSMO-RS can be run from PLAMS using the |CRSJob| class and the corresponding |CRSResults|,
both respectively being subclasses of |SCMJob| and |SCMResults|.

.. note:: There is also `a tutorial showing full code examples <../../COSMO-RS/PLAMS_COSMO-RS_scripting.html>`__ available in the COSMO-RS documentation. There are several templates available that can easily be customized for other problem types, workflows, etc.

Settings
~~~~~~~~

A COSMO-RS job is configured through the PLAMS |Settings| object. For most
workflows, use :meth:`CRSJob.input_builder` to create these settings. The input
builder provides a property-specific interface and converts the builder state to
|Settings| with ``to_settings()`` or directly to a :class:`CRSJob` with
``to_job()``.

For advanced workflows, you can also construct or modify |Settings| manually.


Input builders
^^^^^^^^^^^^^^

The recommended way to prepare COSMO-RS input from PLAMS is
:meth:`CRSJob.input_builder`. This method returns a property-specific input
builder. The builder exposes the input keys, compound roles, and calculation
modes supported by the selected property type.

For example, the following COSMO-RS input performs a pure sigma-profile
calculation:

.. code-block:: none

    compound /path/to/file.coskf
        frac1 1.0
    end

    property puresigmaprofile
        nprofile 50
        sigmamax 0.025
    end

The same input can be prepared with :meth:`CRSJob.input_builder`:

.. code-block:: python

    from scm.plams import CRSJob

    builder = CRSJob.input_builder("PURESIGMAPROFILE")
    builder.nprofile = 50
    builder.sigmamax = 0.025
    builder.add_compound_from_adfcrs_database("Water.coskf")

    job = builder.to_job()
    results = job.run()

Use ``builder.describe()`` to inspect the input keys accepted by a builder:

.. code-block:: python

    for line in builder.describe(include_values=True):
        print(line)

Some property types support calculation modes. For example, solubility can be
set up for a gas-phase solute using ``mode="gas"``:

.. code-block:: python

    builder = CRSJob.input_builder("SOLUBILITY", mode="gas", temperature=298.15)
    builder.add_solvent_from_adfcrs_database("Water.coskf", frac1=1.0)
    builder.add_solute_from_adfcrs_database("Benzene.coskf")

    settings = builder.to_settings()

The supported methods can be inspected from Python:

.. code-block:: python

    CRSJob.methods()


Settings with multiple compounds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    settings = builder.to_settings()


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

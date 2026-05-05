COSMO-RS
--------

(*contributed by* `Bas van Beek <https://www.researchgate.net/profile/Bas_van_Beek>`_\)

.. currentmodule:: scm.plams.interfaces.adfsuite.crs

COSMO-RS can be run from PLAMS using the |CRSJob| class and the corresponding |CRSResults|,
both respectively being subclasses of |SCMJob| and |SCMResults|.

.. note:: There is also `a tutorial showing full code examples <../../COSMO-RS/PLAMS_COSMO-RS_scripting.html>`__ available in the COSMO-RS documentation.  There are several templates available that can easily be customized for other problem types, workflows, etc.

Settings
~~~~~~~~

For example, considering the following input file for a COSMO-RS
sigma-profile calculation [`1 <../../COSMO-RS/Analysis.html#sigma-profile>`_]:

.. code::

    compound /path/to/file.coskf
        frac1 1.0
    end

    property puresigmaprofile
        nprofile 50
        sigmamax 0.025
    end

    temperature 298.15

The input file displayed above corresponds to the following settings:

.. code:: python

    from scm.plams import CRSJob

    s = CRSJob.property_block("PURESIGMAPROFILE", nprofile=50, sigmamax=0.025)
    s.input.compound = [CRSJob.compound_block("/path/to/file.coskf", frac1=1.0)]
    s.input.temperature = 298.15

    my_job = CRSJob(settings=s)
    my_results = my_job.run()

Alternatively one can create a :class:`CRSJob` instance from a runscript created by,
for example, the ADF GUI (``File -> Save as``).

.. code:: python

    from scm.plams import CRSJob

    filename = "path/to/my/crs/inputfile.run"
    my_job = CRSJob.from_inputfile(filename)
    my_results = my_job.run()


Settings with multiple compounds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

More often than not one is interested in the properties of
multi-component mixtures (*e.g.* a dissolved solute).
In such cases one has to pass multiple ``compound`` blocks to
the input file, which is somewhat problematic as Python dictionaries
(including |Settings|) can only contain a set of unique keys.

This problem can be resolved by setting ``compound`` to a list of blocks created
with :meth:`CRSJob.compound_block`. Each item within this list is expanded into
its own ``compound`` block once :meth:`CRSJob.run` creates the actual input file.

Example |Settings| with three compounds:

.. code:: python

    from scm.plams import CRSJob

    compound1 = CRSJob.compound_block("/path/to/compound1.coskf", frac1=0.33)
    compound2 = CRSJob.compound_block("/path/to/compound2.coskf", frac1=0.33)
    compound3 = CRSJob.compound_block("/path/to/compound3.coskf", frac1=0.33)

    s = CRSJob.property_block("ACTIVITYCOEF")
    s.input.compound = [compound1, compound2, compound3]

    my_job = CRSJob(settings=s)
    my_results = my_job.run()

Which yields the following input:

.. code::

    compound /path/to/compound1.coskf
        frac1 0.33
    end

    compound /path/to/compound2.coskf
        frac1 0.33
    end

    compound /path/to/compound3.coskf
        frac1 0.33
    end


ADF and CRSJob
~~~~~~~~~~~~~~

A workflow is presented in the `ADF and COSMO-RS workflow Python example <../../PythonExamples/ams-crs-workflow/index.html>`__.
In this workflow, we follow the usual procedure of generating the inputs required to run COSMO-RS and COSMO-SAC calculations.


Property and method blocks
~~~~~~~~~~~~~~~~~~~~~~~~~~

A large number of configurable parameters_ is available for COSMO-RS and COSMO-SAC.
Use :meth:`CRSJob.property_block` to create the property-specific defaults, and use
:meth:`CRSJob.method_block` to select a COSMO-RS or COSMO-SAC method. These helpers
validate the requested property type, method name, and method-parameter names.

For example, a solubility calculation with the COSMO-SAC defaults can be created as follows:

.. code:: python

    from scm.plams import CRSJob

    s = CRSJob.property_block("SOLUBILITY")
    s += CRSJob.method_block("COSMOSAC", include_defaults=True)
    s.input.temperature = "273.15 283.15 10"
    s.input.compound = [
        CRSJob.compound_block("/path/to/Water.coskf", frac1=1.0),
        CRSJob.compound_block(
            "/path/to/Benzene.coskf",
            frac1=0.0,
            meltingpoint=278.7,
            hfusion=2.37,
        ),
    ]

The supported property types and methods can be inspected from Python:

.. code:: python

    CRSJob.property_types()
    CRSJob.property_type_metadata("SOLUBILITY", as_summary=True)
    CRSJob.methods()

.. _parameters: ../../COSMO-RS/COSMO-RS_and_COSMO-SAC_parameters.html


Data analyses and plotting
~~~~~~~~~~~~~~~~~~~~~~~~~~

As COSMO-RS can produce a large variety of data series,
a number of specialized methods are available in the :class:`CRSResults` for their extraction and analysis.
The resulting data is stored in either a dictionary of NumPy arrays or (optionally) a `Pandas DataFrame`_.

The sigma-profile and sigma-potential data can be further customized by altering
the ``subsection`` argument.
:meth:`CRSResults.get_results` returns the calculated properties for the current CRS result section.
For activity coefficient jobs, :meth:`CRSResults.get_activity_coefficient` and
:meth:`CRSResults.get_energy` provide direct access to commonly used scalar results.

A complete overview of all available sections and subsections can be acquired
by calling the :meth:`.KFFile.get_skeleton` method of the KF binary file
(*e.g.* :code:`print(my_results._kf.get_skeleton())`).

============================ =======================================
Quantity                     Method for data extraction
============================ =======================================
`Sigma profile`_             :meth:`CRSResults.get_sigma_profile`
`Sigma potential`_           :meth:`CRSResults.get_sigma_potential`
Activity coefficient         :meth:`CRSResults.get_activity_coefficient`
Solvation energy             :meth:`CRSResults.get_energy`
Other calculated properties  :meth:`CRSResults.get_results`
============================ =======================================

If the `Matplotlib <https://matplotlib.org/>`_ package is installed, the resulting data can be plotted by passing
it to the :meth:`CRSResults.plot` method (*e.g.* :code:`CRSResults.plot(my_sigma_profile)`):

.. code:: python

    from scm.plams import CRSJob
    import numpy as np

    s = CRSJob.property_block("PURESIGMAPROFILE", nprofile=50, sigmamax=0.025)
    s.input.compound = [CRSJob.compound_block("/path/to/Water.coskf", frac1=1.0)]
    s.input.temperature = 298.15

    my_job = CRSJob(settings=s)
    my_results = my_job.run()

    my_sigma_profile = my_results.get_sigma_profile()
    with np.printoptions(threshold=0, edgeitems=5):
        print(my_sigma_profile)

    my_results.plot(my_sigma_profile)

.. image:: ../_static/sigma_profile.png


.. _`Sigma profile`: ../../COSMO-RS/Analysis.html#sigma-profile
.. _`Sigma potential`: ../../COSMO-RS/Analysis.html#sigma-potential

.. _`Pandas DataFrame`: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html


API
~~~

.. autoclass:: CRSJob
    :no-private-members:
    :no-special-members:

.. autoclass:: CRSResults
    :no-private-members:
    :no-special-members:

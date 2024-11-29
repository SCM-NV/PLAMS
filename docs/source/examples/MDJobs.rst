.. _AMSMDJob:

AMS Molecular Dynamics PLAMS jobs
-------------------------------------

.. seealso::

    **Example**: `Molecular dynamics with Python tutorial <../../Tutorials/MolecularDynamicsAndMonteCarlo/MDintroPython/intro.html>`__

.. currentmodule:: scm.plams.recipes.md.amsmdjob

In AMS2023, the following special Jobs exist to simplify running MD simulations (and continuing/restarting MD simulations):

* ``AMSMDJob``. A general class for which you can set up most common MD options using arguments to the constructor.

* ``AMSNVEJob`` for NVE simulations

* ``AMSNVTJob`` for NVT simulations

* ``AMSNPTJob`` for NPT simulations

* ``AMSMDScanDensityJob`` for running MD deformation simulations while isotropically scaling the density

* ``AMSNVESpawnerJob`` is a special MultiJob that runs several NVE simulations with initial velocities taken from evenly spaced frames in a previous job.

Some default values are different from AMS. For example, the checkpoint
frequency is set to a higher number, and the thermostat constant ``tau`` is
automatically set to 400 times the timestep, by default.

Always check the input by calling ``job.get_input()`` before running a job.

In general, the jobs and results classes inherit from |AMSJob| and |AMSResults|.

The NVE, NVT, and NPT classes simply remove the below options if they are set:

.. csv-table::
    
    allowed input block, AMSMDJob, NVE, NVT, NPT
    Thermostat, yes, no, yes, yes
    Barostat, yes, no, no, yes
    Nanoreactor, yes, no, no, no
    Deformation, yes, no, no, no


The following jobs help with the postanalysis:

* ``AMSMSDJob`` for calculating mean square displacement (MSD)

* ``AMSRDFJob`` for calculating radial distribution functions (RDF)

* ``AMSVACFJob`` for calculating velocity autocorrelation functions (VACF)


AMSMDJob API
~~~~~~~~~~~~~

.. autoclass:: AMSMDJob
    :exclude-members: _result_type
    :no-private-members:

AMSNVEJob API
~~~~~~~~~~~~~~

This class uses the same arguments as ``AMSMDJob``.

.. autoclass:: AMSNVEJob
    :no-private-members:

AMSNVTJob API
~~~~~~~~~~~~~~

This class uses the same arguments as ``AMSMDJob``.

.. autoclass:: AMSNVTJob
    :no-private-members:


AMSNPTJob API
~~~~~~~~~~~~~~

This class uses the same arguments as ``AMSMDJob``.

``AMSNPTResults`` inherits from ``AMSResults``.

.. autoclass:: AMSNPTJob
    :no-private-members:


.. autoclass:: AMSNPTResults
    :no-private-members:

.. currentmodule:: scm.plams.recipes.md.nvespawner

AMSNVESpawnerJob API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: AMSNVESpawnerJob
    :no-private-members:

.. currentmodule:: scm.plams.recipes.md.scandensity

AMSMDScanDensityJob API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: AMSMDScanDensityJob
    :no-private-members:


.. currentmodule:: scm.plams.recipes.md.trajectoryanalysis

AMSRDFJob API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: AMSRDFJob
    :no-private-members:

AMSMSDJob API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: AMSMSDJob
    :no-private-members:

AMSVACFJob API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: AMSVACFJob
    :no-private-members:


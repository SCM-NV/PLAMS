.. _AMSCalculator:

ASE Calculator for AMS
-----------------------

.. currentmodule:: scm.plams.interfaces.adfsuite.ase_calculator

Introduction
~~~~~~~~~~~~~

The |AMSCalculator| class is an ASE Calculator that lets you use any AMS engine (ADF, BAND, DFTB, ReaxFF, MLPotential, ForceField, ...)
together with ASE.

.. seealso::

    **Example**: :ref:`ASECalculatorExample` 

.. important::

    Before using an AMSCalculator, you must call the ``init()`` function from PLAMS. This is done automatically if you use the ``plams`` binary.

AMS settings
~~~~~~~~~~~~~~~~

A |Settings| object is used to set up the AMS settings in the same way as for normal PLAMS jobs.

.. seealso::

    :ref:`Preparing input for AMS jobs in PLAMS <AMS_preparing_input>`

ASE results
...........

Currently only the energy, forces and stress tensor are provided through the ASE interface.
All other results can be accessed through ``AMSCalculator.ams_results``, which is an |AMSResults| object.
Use ``AMSCalculator.ensure_property('forces')`` and ``AMSCalculator.ensure_property('stress')`` to ensure
these ASE properties are computed regardless of whether it was originally requested in the |Settings| object.

Charge
......

There is currently no universal interface in ASE for the total charge of a system and is instead considered to be Calculator specific.
The easiest way to set the charge a calculation with the |AMSCalculator| is to define ``Atoms.info['charge']``.
Additionally, when the charge needs to be treated extensively w.r.t. manipulations of the ``Atoms`` object in ASE, the initial charge of each atom can also be set.
The total charge is thus obtained as ``sum(Atoms.get_initial_charges())+Atoms.info['charge']``.
See the ASE documentation for details on initial charges and info.

.. seealso::
   **Example**: :ref:`ChargedAMSCalculatorExample`

AMS standalone and worker mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AMS can run in two modes: standalone and worker.

In **standalone mode** (``amsworker=False``), AMS is started for every new structure, and stores the normal AMS output. Use this mode if:

* you need to access all results, *or*
* you run a calculation not supported by the AMSworker (e.g., a PESScan)

In **AMSworker mode** (``amsworker=True``), AMS is only started once. This significantly reduces overhead, especially for fast engines like ReaxFF. The downside is that you cannot access all results. Use this mode if:

* you only need access to energy, forces, and stress tensor, *and*
* you only need to run single point calculations, for example to use the ``ase.optimize.BFGS`` geometry optimizer.

.. technical::


    In **AMSworker mode**,  make sure to use |AMSCalculator| as a context manager, or explicitly
    call ``AMSCalculator.stop_worker()`` when you want to terminate the AMS process. 


Technical
~~~~~~~~~~~~~~~~~~~~~~~~~

.. technical::

    When creating a deepcopy of |AMSCalculator|
    the |AMSWorker| is not copied and instead every copy of |AMSCalculator| has a reference to the
    same |AMSWorker|.

AMSCalculator API
~~~~~~~~~~~~~~~~~

.. autoclass:: AMSCalculator
    :exclude-members: __init__, __weakref__, results_from_ams_results
    :no-private-members:


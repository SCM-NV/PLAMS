.. _job_runners:

Job runners
-------------

.. currentmodule:: scm.plams.core.jobrunner

Job runners have already been mentioned in previous chapters about jobs and results.
Here we sum up all that information and introduce a basic |JobRunner| object together with its subclass |GridRunner| which is meant for interacting with queueing systems that manage resources on computer clusters.

Job runners in PLAMS are very simple objects, both from user's perspective and in terms of internal architecture.
They have no methods that are meant to be called in your scripts, apart from constructors.
Job runners are supposed to be created (with some parameters adjusting their behavior) and passed to the |run| method as parameters (or placed as ``config.default_jobrunner``).

Local job runner
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: JobRunner
    :exclude-members: __weakref__, __metaclass__

.. technical::

    Similarly to the |Results| class, the proper behavior of |JobRunner| and its subclasses (also the ones defined by the user) is ensured using a metaclass.
    For the sake of completeness we present here a brief specification of all involved elements:

    .. autoclass:: _MetaRunner
    .. autofunction:: _limit
    .. autofunction:: _in_thread

Remote job runner
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: GridRunner
    :exclude-members: __weakref__, __metaclass__
.. _job_manager:

Job manager
-------------

.. currentmodule:: scm.plams.core.jobmanager

Job manager is the "commander" of PLAMS environment.
It creates the structure of the working folder, manages its contents, and keeps track of all jobs you run.

Every instance of |JobManager| is tied to a working folder.
This folder is created when |JobManager| instance is initialized and all the jobs managed by that instance have their job folders inside the working folder.
By default, the working folder is located in the directory where a script is executed and is called ``plams_workdir`` (``plams_workdir.002`` if ``plams_workdir`` already existed).
You should not change job manager's working folder after it has been created.

Within a PLAMS environment, there is an instance of |JobManager| stored in ``config.default_jobmanager``.
This can be left to default, or explicitly initialized, for example with:

.. code-block:: python

    # Either of the following are equivalent
    config.default_jobmanager = JobManager(config.jobmanager, "my/path", "my_folder")
    init(path="my/path", folder="my_folder")

Aside from this, it is rare to interact explicitly with a |JobManager| instance (create it manually, call any of its methods, explore its data etc.).
All interactions are handled automatically from |run| or other methods.

.. technical::

   Usually there is no need to use any other job manager than the default one.
   Splitting work between multiple instances of |JobManager| may lead to some problems (different instances don't communicate, so |RPM| does not work properly).

   However, it is possible to manually create another instance of |JobManager| (with a different working folder) and use it for part of your jobs (by passing it as ``jobmanager`` keyword argument to |run|).
   If you decide to do so, make sure to pass all instances of |JobManager| you manually created to |finish| (as a list).

   An example application for that could be running jobs within your script on many different machines (for example via SSH) and having a separate |JobManager| on each of them.



.. _rerun-prevention:

Rerun prevention
~~~~~~~~~~~~~~~~~~~~~~~~~

In some situations, for example when running many automatically generated small jobs, it may happen that two or more jobs are identical -- they have the same input files.
PLAMS has a built-in mechanism to detect such situations and avoid unnecessary work.

During |run|, just before the actual job execution, a unique identifier (called *hash*) of a job is calculated.
Job manager stores all hashes of previously started jobs and checks if the hash of the job you are just running has already occurred.
If such a situation is detected, the execution of the current job is skipped and results of the previous job are used.
Results from previous job's folder can be either copied or linked to the current job's folder, based on ``link_files`` key in **previous** job's ``settings``.

.. note::

    Linking is done using hard links.
    Windows does not support hard links so if you are running PLAMS under Windows results are always copied.

The crucial part of the whole rerun prevention logic is a properly working :meth:`~scm.plams.core.basejob.Job.hash` function.
It is a function that takes the whole job instance and produces its hash.
The hashing function needs to produce different hashes for different jobs and exactly the same hashes for jobs that do exactly the same work.
It is far from trivial to come up with the scheme that works well for all kind of external binaries, since the technical details about job preparation can differ a lot.
Currently implemented method works based on calculating SHA256 hash of input and/or runscript contents.
The value of ``hashing`` key in job manager's ``settings`` can be one of the following: ``'input'``, ``'runscript'``, ``'input+runscript'`` (or ``None`` to disable the rerun prevention).

If you decide to implement your own hashing method, it can be done by overriding :meth:`~scm.plams.core.basejob.SingleJob.hash_input` and/or meth:`~scm.plams.core.basejob.SingleJob.hash_runscript`.

.. warning::

    It may happen that two jobs with the same input and runscript files correspond to different jobs (for example, if they rely on some external file that is supplied using relative path).
    Sometimes it's even a desired behavior to run multiple different copies of the same job (for example, multiple MD trajectories with the same starting point and random initial velocities).
    If you are experiencing problems (PLAMS refuses to run a job, because it was already run in the past), you can disable the rerun prevention with ``config.default_jobmanager.settings.hashing = None``.

Hashing is disabled for |MultiJob| instances since they don't have inputs and runscripts.
Of course single jobs that are children of multijobs are hashed in the normal way, so trying to run exactly the same multijob twice will not trigger rerun prevention on the multijob level, but rather for every children job separately, effectively preventing any doubled work.



.. _pickling:

Pickling
~~~~~~~~~~~~~~~~~~~~~~~~~

The lifetime of the whole PLAMS environment is limited to a single script.
That means every PLAMS script you run uses its own independent job manager, working folder and ``config`` settings.
These objects are initialized at the beginning of the script and they cease to exist when the script ends.
Also all the settings adjustments (apart from those done by editing :ref:`global-settings`) are local and they affect only the current script.

As a consequence of that, the |JobManager| of the current script is not aware of any jobs that had been run in past scripts.
But often it would be very useful to import a previously run job to the current script and use its results or build some new jobs based on it.
For that purpose PLAMS offers data preservation for job objects.
Every time an execution of a job successfully finishes (see :ref:`job-life-cycle`), the whole job object is saved to a ``.dill`` file using Python serialization called :mod:`pickling<pickle>`.
Such a ``.dill`` file can be loaded ("unpickled") in future scripts using |load| function::

    oldjob = load('/home/user/science/plams_workdir/myjob/myjob.dill')

This operation brings back the old |Job| instance in (almost) the same state it was just after its execution finished.

.. note::

    The default Python pickling package :mod:`pickle<pickle>` is not powerful enough to handle some of common PLAMS objects.
    Fortunately, the `dill <https://pypi.python.org/pypi/dill>`_ package provides an excellent replacement for ``pickle``, following the same interface and being able to save and load almost everything.
    It is strongly recommended to use ``dill`` to ensure proper work of PLAMS data preservation logic.
    However, if ``dill`` is not installed for the Python interpreter you're using to run PLAMS, the regular ``pickle`` package will be used instead (which can work if your |Job| objects are not too fancy, but in most cases it will probably fail).
    Please use ``dill``, it's free, easy to get and awesome.


The pickling mechanism follows references in pickled object.
That means if an object you are trying to pickle contains a reference to another object (just like a |Job| instance has a reference to a |Results| instance), that other object is saved too.
Thanks to that there are no "empty" references in your objects after unpickling.
However, every |Job| instance in PLAMS has a reference to a job manager, which in turns has references to all other jobs, so pickling one job would effectively mean pickling the whole environment.
To avoid that, every |Job| instance needs to be prepared for pickling by removing references to "global" objects, as well as some local attributes (path to the job folder for example).
During loading, all the removed data is replaced with "proper" values (current job manager, current path to the job folder etc.).

.. note::

    There is a way of expanding the mechanism explained above.
    If your |Job| object has an attribute with reference to an object you don't want to save together with the job, you may add this object's name to job's ``_dont_pickle`` list::

        myjob.something = some_big_and_clumsy_object_you_dont_want_to_pickle
        myjob._dont_pickle.append('something')

    That way big clumsy object will not be stored in the ``.dill`` file.
    After loading such a ``.dill`` file the value of ``myjob.something`` will simply be ``None``.

    The ``_dont_pickle`` list is an attribute of every |Job| instance, initially an empty list.
    It does not contain names of attributes that are always removed (like ``jobmanager``), it's meant only for additional ones defined by the user (see :meth:`Job.__getstate__<scm.plams.core.basejob.Job.__getstate__>`)


As mentioned above, pickling a job happens at the very end of |run|.
The decision if a job should be pickled is based on the ``pickle`` key in job's ``settings``, so it can be adjusted for each job separately.
If you wish not to pickle a particular job just set ``myjob.settings.pickle = False``.
Of course the global default ``config.job.pickle`` can also be used.

If you modify a job or its corresponding |Results| instance after it has been pickler, these changes are not going to be reflected in the ``.dill`` file, since it was created before the changes happened.
To update the state of the ``.dill`` file to include such changes you need to repickle the job manually by calling ``myjob.pickle()`` after doing your changes.

.. note::

    Not all Python objects can be properly pickled, so you need to be careful what other objects your |Job| (or its |Results|) store references to.

The |Results| instance associated with the job is saved together with it.
However, these results don't contain all files produced by the job execution, but only relative paths to them.
For that reason the ``.dill`` file is not enough to fully restore the job object if you want to extract or process the results.
All other files present in the job folder are needed so that |Results| instance can see them.
So if you want to copy a previously run job to another location make sure to copy *the whole* job folder (including subdirectories).

A job loaded with |load| is **not** registered in the current job manager.
That means it does not get its own subfolder in the current working folder, it never gets renamed and no |cleaning| is done on |finish|.
However, it is added to the hash registry, so it is visible to |RPM|.

In case of a |MultiJob| all the information about children jobs is stored in parent's ``.dill`` file so loading a |MultiJob| results in loading all its children.
Each child job can have its own ``.dill`` file containing information about that particular job only.
When pickling, the ``parent`` attribute of a job is erased, so loading a child job does not result in loading its parent (and all other children).



.. _restarting:

Restarting scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~

Pickling and rerun prevention combine together into a handy restart mechanism.
When your script tries to do something "illegal", an exception is raised and the script gets terminated by the Python interpreter.
Usually it is caused by a mistake in the script (a typo, using wrong variable, accessing wrong element of a list etc.).
In such a case one would like to correct the script and run it again.
But some jobs in the terminated script may had already been run and successfully finished before the exception occurred.
It would be a waste of time to run those jobs again in the corrected script if they are meant to produce exactly the same results as previously.
The solution is to load all successful jobs from the old script at the beginning of the new one and let |RPM| do the rest.
But having to go to the old script's working folder and manually get paths to all ``.dill`` files present there would be cumbersome.
Fortunately, one can use |load_all| function which takes a path to the main working folder of some finished PLAMS run and loads all ``.dill`` files present there.
So when you edit your crashed script to remove mistakes you can add just one |load_all| call at the beginning.
Then you run your corrected script and no unnecessary work is done: all the finished jobs are loaded from the previous run, the current run tries to run the same jobs again, but |RPM| detects that and copies/links old jobs' folders into the current main working folder.

.. note::
    Please remember that rerun prevention checks the hash of the job after the |prerun| method is executed.
    So when you attempt to run a job identical to the one previously run (in the same script, or imported from a previous run), its |prerun| method is executed anyway, even if the rest of :ref:`job-life-cycle` is skipped.


API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: JobManager
    :exclude-members: __weakref__

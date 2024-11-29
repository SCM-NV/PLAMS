Results
-------

.. currentmodule:: scm.plams.core.results

The goal of |Results| object is to take care of the job folder after the execution is finished: gather information about produced files, help to manage them and extract data of interest from them.
Every |Job| instance has an associated |Results| instance created automatically on job creation and stored in its ``results`` attribute.

From the technical standpoint, |Results| class is the part of PLAMS environment responsible for thread safety and proper synchronization in parallel job execution.


Files in the job folder
~~~~~~~~~~~~~~~~~~~~~~~

Directly after the execution of a job is finished (see :ref:`job-life-cycle`), the job folder gets scanned by :meth:`~Results.collect` method.
All files present in the job folder, including files in subfolders, are gathered in a list stored in ``files`` attribute of the |Results| instance.
Entries in this list correspond to paths to files relative to the job folder, so files on the top level are stored by their names and files in subfolders by something like ``childjob/childjob.out``.

.. note::

    Files produced by |pickling| are excluded from the ``files`` list.
    Every file with ``.dill`` extension is simply ignored by |Results|.

If you need an absolute path to some file, the bracket notation known from dictionaries is defined for |Results| objects.
When supplied with an entry from ``files`` list, it returns the absolute path to that file.
The bracket notation is read-only::

    >>> r = j.run()
    >>> print(r.files)
    ['plamsjob.err', 'plamsjob.in', 'plamsjob.out', 'plamsjob.run']
    >>> print(r['plamsjob.out'])
    /home/user/plams.12345/plamsjob/plamsjob.in
    >>> r['newfile.txt'] = '/home/user/abc.txt'
    TypeError: 'Results' object does not support item assignment

In the bracket notation, and in every other context regarding |Results|, whenever you need to pass a string with a filename, a shortcut ``$JN`` can be used for the job name::

    >>> r.rename('$JN.out', 'outputfile')
    >>> r.grep_file('$JN.err', 'NORMAL TERMINATION')
    >>> print(r['$JN.run'])
    /home/user/plams.12345/plamsjob/plamsjob.run

Some produce produce fixed name files during execution.
If one wants to automatically rename those files it can be done with ``_rename_map`` class attribute -- a dictionary defining which files should be renamed and how.
Renaming is done during :meth:`~Results.collect`.

In the generic |Results| class ``_rename_map`` is an empty dictionary.


.. _parallel:

Synchronization of parallel job executions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the main advantages of PLAMS is the ability to run jobs in parallel.
The whole job execution logic is designed in a way that does not require a special parallel script for a parallel workflow execution.
Exactly the same scripts can be used for both serial and parallel execution.

However, it is important to have a basic understanding of how parallelism in PLAMS works to avoid potential deadlocks and maximize the performance of your scripts.

To run your job in parallel you need to use a parallel job runner::

    pjr = JobRunner(parallel=True)
    myresults = myjob.run(jobrunner=pjr)

Parallelism is not something that is "on" or "off" for the entire script: within one script you can use multiple job runners, some of them may be parallel and some may be serial.
However, if you wish to always use the same |JobRunner| instance, it is convenient to set is as default at the beginning of your script::

    config.default_jobrunner = JobRunner(parallel=True)

All |run| calls without explicit ``jobrunner`` argument will now use that instance.

When you run a job using a serial job runner, all steps of |run| (see :ref:`job-life-cycle`) are done in the main thread and |Results| instance is returned at the end.
On the other hand, when a parallel job runner is used, a new thread is spawned at the beginning of |run| and all further work is done in this thread.
Meanwhile the main thread proceeds with the next part of the script.
The important thing is that the |run| method called in the main thread returns a |Results| instance and allows the whole script to proceed even though the job is still running in a separate thread.
This |Results| instance acts as a "guardian" protecting the job from being accessed while it is still running.
Every time you call a method of a |Results| instance, the guardian checks the status of the job and, if the job is not yet finished, forces the thread from which the call was done to wait.
Thanks to that there is no need to explicitly put synchronization points in the script -- results requests serve for that purpose.

.. warning::

    You should **NEVER** access results in any other way than by a **method** of some |Results| instance.

The |Results| class is designed in such a way, that each of its methods automatically gets wrapped with the access guardian when a |Results| instance is created.
That behavior holds for any |Results| subclasses and new methods defined by user, so no need to worry about guardian when extending |Results| functionality or subclassing it.
Also |binding_decorators| recognize when you try to use them with |Results| and act accordingly.
Methods whose names end with two underscores, as well as :meth:`~Results.refresh`, :meth:`~Results.collect`, :meth:`~Results._clean` are not wrapped with the guardian.
The guardian gives special privileges (earlier access) to |postrun| and :meth:`~scm.plams.core.basejob.Job.check` (see :ref:`prerun-postrun`).

If you never request any results of your job and just want to run it, |finish| method works as a global synchronization point.
It waits for all spawned threads to end before cleaning the environment and exiting your script.

Examples
++++++++

This section provides a handful of examples together with an explanation of common pitfalls and good practices one should keep in mind when writing parallel PLAMS scripts.

Let us start with a simple parallel script that takes all ``.xyz`` files in a given folder and for each one calculates the dipole moment magnitude using a single point ADF calculation:

.. code-block:: python
    :linenos:

    config.default_jobrunner = JobRunner(parallel=True)
    config.log.stdout = 1
    folder = '/home/user/xyz'
    molecules = read_molecules(folder)

    s = Settings()
    s.input.ams.task = 'singlepoint'
    s.input.adf.basis.type = 'DZP'
    s.input.adf.xc.gga = 'PBE'

    jobs = [AMSJob(molecule=molecules[name], name=name, settings=s) for name in sorted(molecules)]
    results = [job.run() for job in jobs]

    for r in results:
        dipole_vec = r.readrkf('AMSResults', 'DipoleMoment', file='engine')
        dipole_magn = sum([a*a for a in dipole_vec])**0.5
        print('{}\t\t{}'.format(r.job.name, dipole_magn))

For an explanation purpose let us assume that ``/home/user/xyz`` contains three files: ``ammonia.xyz``, ``ethanol.xyz``, ``water.xyz``.
When you run this script the standard output will look something like:

.. code-block:: none

    [17:31:52] PLAMS working folder: /home/user/plams_workdir
    [17:31:52] JOB 'ammonia' STARTED
    [17:31:52] JOB 'ethanol' STARTED
    [17:31:52] JOB 'water' STARTED
    [17:31:52] Waiting for job ammonia to finish
    [17:31:56] JOB 'water' SUCCESSFUL
    [17:31:56] JOB 'ammonia' SUCCESSFUL
    ammonia     0.5949493793257181
    [17:31:56] Waiting for job ethanol to finish
    [17:32:01] JOB 'ethanol' SUCCESSFUL
    ethanol     0.5946259677193089
    water       0.7082267171673067

As you can see, print statements from line 17 are mixed with automatic logging messages.
Let us examine in more detail what causes such a behavior.
To do so we will follow the main thread.
In line 11 an alphabetically sorted list of jobs is created, so the job named ``'ethanol'`` will come after ``'ammonia'`` and before ``'water'``.
Line 12 is a for loop that goes along the list of jobs, runs each of them and collects their |Results| instances in a new list called ``results``.
If we were using a serial job runner, all the computational work would happen in line 12: the ``'ethanol'`` job would start only when ``'ammonia'`` was finished, ``'water'`` would wait for ``'ethanol'`` and the main thread would proceed to the next line only when ``'water'`` is done.

In our case, however, we are using a parallel job runner.
The first job (``'ammonia'``) is started and quickly moves to a separate thread, allowing the main thread to proceed to another instruction, which in this case is the |run| method of the ``'ethanol'`` job.
Thanks to that all three jobs are started almost immediately one after another, corresponding |Results| are gathered and the main thread proceeds to line 14, while the three jobs are running "in the background", handled by separate threads.
Now the main thread goes along the ``results`` list (which follows the same order as ``jobs``) and tries to obtain a dipole vector from each job.
It uses ``readkf`` method of |Results| instance associated with the ``'ammonia'`` job and since this job is still running, the main thread hangs and waits for the job to finish (*"Waiting for job ammonia to finish"*).
Meanwhile we can see that the ``'water'`` job ends and this fact is logged.
Quickly after that also the ``'ammonia'`` job finishes and the main thread obtains ``dipole_vec``, calculates ``dipole_magn`` and prints it.
Now the ``for`` loop in line 14 continues, this time for the ``'ethanol'`` job.
This job seems to be a bit longer than ``'ammonia'``, so it is still running and the main thread again hangs on the ``readkf`` method (*"Waiting for job ethanol to finish"*).
After finally obtaining the dipole vector of ethanol, calculating the magnitude and printing it, the ``for`` loop goes on with its last iteration, the ``'water'`` job.
This time there is no need to wait since the job is already finished - the result is calculated and printed immediately.

Knowing that, let us wonder what would happen if the order of jobs was different.
If ``'ethanol'`` was the first job on the list, by the time its results would be obtained and printed, both other jobs would have finished, so no further waiting would be needed.
On the other hand, if the order was ``'water'``--``'ammonia'``--``'ethanol'``, the main thread would have to wait every time when executing line 15.

The most important lesson from the above is: the order in which you start jobs does not matter (too much), it is the order of results requests that makes the difference.
Of course in our very simple example it influences only the way in which results are mixed with log messages, but in more complicated workflows it can directly affect the total runtime of your script.

By the way, to avoid print statements being mixed with logging messages one could first store the data and print it only when all the results are ready::

    to_print = []
    for r in results:
        dipole_vec = r.readkf('Properties', 'Dipole')
        dipole_magn = sum([a*a for a in dipole_vec])**0.5
        to_print += [(r.job.name, dipole_magn)]
    for nam, dip in to_print:
        print('{}\t\t{}'.format(nam, dip))

Another way could be disabling logging to the standard output by putting ``config.log.stdout = 0`` at the beginning of the script (see |log|).

Coming back to the main topic of our considerations, as we have seen above, parallelism in PLAMS is driven by results request.
Not only the order of requests is important, but also (probably even more important) the place from which they are made.
To picture this matter we will use the following script that performs geometry optimization followed by frequencies calculation of the optimized geometry:

.. code-block:: python
    :linenos:

    config.default_jobrunner = JobRunner(parallel=True)

    go = AMSJob(name='GeomOpt', molecule=Molecule('geom.xyz'))
    go.settings.input.ams.task = 'GeometryOptimization'
    ... #other settings adjustments for geometry optimisation
    go_results = go.run()

    opt_geo = go_results.get_main_molecule()

    freq = AMSJob(name='Freq', molecule=opt_geo)
    freq.settings.input.ams.properties.NormalModes = 'Yes'
    ... #other settings adjustments for frequency run
    freq_results = freq.run()

    do_other_work() # further part of the script, independent of GeomOpt and Freq

Again let us follow the main thread.
In line 8 we can see a results request for the optimized geometry from "GeomOpt" job.
The main thread will wait for that job to finish before preparing the "Freq" job and running it.
That means ``do_other_work()``, whatever it is, will not start before "GeomOpt" is done, even though it could, since it is independent of GeomOpt and Freq results.
This is bad.
The main thread is wasting time that could be used for ``do_other_work()`` on idle waiting.
We need to fix the script:

.. code-block:: python
    :linenos:

    config.default_jobrunner = JobRunner(parallel=True)

    go = AMSJob(name='GeomOpt', molecule=Molecule('geom.xyz'))
    go.settings.input.ams.task = 'GeometryOptimization'
    ... #other settings adjustments for geometry optimisation
    go_results = go.run()

    freq = AMSJob(name='Freq')
    freq.settings.input.ams.properties.NormalModes = 'Yes'
    ... #other settings adjustments for frequency run

    @add_to_instance(freq)
    def prerun(self):
        self.molecule = go_results.get_main_molecule()

    freq_results = freq.run()

    do_other_work() # further part of the script, independent of GeomOpt and Freq

The results request (``go_results.get_main_molecule()``) have been moved from the main script to the |prerun| method of the "Freq" job.
The |prerun| method is executed in job's thread rather than in the main thread.
That means the main thread starts the "Freq" job immediately after starting the "GeomOpt" job and then directly proceeds to ``do_other_work()``.
Meanwhile in the thread spawned for the "Freq" job the result request for molecule is made and only that thread waits for "GeomOpt" to finish.

As seen in the above example, it is extremely important to properly configure jobs that are dependent (setup of one depends on results of another).
Resolving all such dependencies in job's thread rather than in the main thread guarantees that waiting for results is done only by the code that really needs them.

.. note::

    In some cases dependencies between job are not easily expressed via methods of |Results| (for example, one job sets up some environment that is later used by another job).
    In such cases one can use job's ``depend`` attribute to explicitly tell the job about other jobs it has to wait for.
    Adding ``job2`` to ``job1.depend`` is roughly equivalent to putting ``job2.results.wait()`` in ``job1`` |prerun|.


To sum up all the above considerations, here is the rule of thumb on how to write properly working parallel PLAMS scripts:

1.  Request results as late as possible, preferably just before using them.
2.  If possible, avoid requesting results in the main thread.
3.  Place the result request in the thread in which that data is later used.



.. _cleaning:

Cleaning job folder
~~~~~~~~~~~~~~~~~~~~~~~~~

The |Results| instance associated with a job is responsible for cleaning the job folder (removing files that are no longer needed).
Cleaning is done automatically, twice for each job.

First cleaning is done during |run|, just after :meth:`~scm.plams.core.basejob.Job.check` and before |postrun|.
The value adjusting this first cleaning is taken from ``myjob.settings.keep`` and should be either a string or a list (see below).

This cleaning is intended for situations when your jobs produce large files that you don't need for further processing.
Running many of such jobs could deplete the disk space and cause the whole script to crash.
If you wish to immediately get rid of some files produced by your jobs (without having a chance to do anything with them), use the first cleaning.

In the majority of cases it is sufficient to use the second cleaning, which is performed at the end of your script, or when |finish| method is explicitly called.
It is adjusted by ``myjob.settings.save``.
You can use the second cleaning to remove files that you no longer need after you extracted relevant data earlier in your script.

The argument passed to :meth:`~Results._clean` (in other words the value that is supposed to be kept in ``myjob.settings.keep`` and ``myjob.settings.save``) can be one of the following:

*   ``'all'`` -- nothing is removed, cleaning is skipped.
*   ``'none'`` or ``[]`` or ``None`` -- everything is removed from the job folder.
*   list of strings -- list of filenames to be kept.
    Shortcut ``$JN`` can be used here, as well as \*-wildcards.
    For example ``['geo.*', '$JN.out', 'logfile']`` will keep ``[jobname].out``, ``logfile`` and all files whose names start with ``geo.`` and remove everything else from the job folder.
*   list of strings with the first element ``'-'`` -- reversed behavior to the above, listed files will be removed.
    For example ``['-', 't21.*', '$JN.err']`` will remove ``[jobname].err`` and all files whose names start with ``t21.``



Cleaning for multijobs
+++++++++++++++++++++++++

Cleaning happens for every job run with PLAMS, either single job or multijob.
That means, for example, that a single job that is a child of some multijob will have its job folder cleaned by two different |Results| instances: it's own |Results| and its parent's |Results|.
Those two cleanings can interfere with each other.
Hence it is a good practice to set cleaning only on one level (either in a parent job or in children jobs) and disable cleaning on the other level, by using ``'all'``.

Another shortcut can be used for cleaning in multijobs: ``$CH`` is expanded with every possible child name.
For example, if you have a multijob ``mj`` with 5 single job children (``child1``, ``child2`` and so on) and you wish to keep only input and output files of children jobs you can set::

    mj.settings.save = ['$CH/$CH.in', '$CH/$CH.out']

It is equivalent to::

    mj.settings.save = ['child1/child1.in', 'child2/child2.in', ... , 'child1/child1.out', 'child2/child2.out', ...]

As you can see above, while cleaning a multijob folder you have to keep in mind that files in subfolders are kept as relative paths.

API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    These functions work for all Results type. See the |AMSResults| class for functions extracting or manipulating results from AMS jobs.

.. autoclass:: Results
    :exclude-members: __weakref__, __metaclass__

.. technical::

    Other parts of ``results`` module described below are responsible for giving |Results| class its unique behavior described in |parallel|.
    They are presented here for the sake of completeness, from a user's perspective this information is not too relevant.

    .. autoclass:: ApplyRestrict
        :special-members:

    .. autofunction::  _caller_name_and_arg
    .. autofunction:: _privileged_access



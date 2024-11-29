AMS worker
----------

.. currentmodule:: scm.plams.interfaces.adfsuite.amsworker

The PLAMS interface to the AMS driver and its engines through the |AMSJob| and |AMSResults| questions is technically quite similar to how one would run calculations from the command line or from the GUI: the AMS input files are written to disk, the AMS driver starts up, reads its input, performs the calculations and writes the results to disk in form of human readable text files as well as machine readable binary files, usually in :ref:`KF format<kf-files>`. This setup has the advantage that any calculation that can be performed with AMS can be setup from PLAMS as an |AMSJob|, and that any result from any calculation can be accessed from PLAMS through the corresponding |AMSResults| instance. Furthermore, the resulting files on disk can often visualized using the AMS GUI, as if the job had been set up and run through the graphical user interface. As such, this way of running AMS offers maximum flexibility and convenience to users.

However, for simple and fast jobs where we only care about some basic results, this flexibility comes at a cost: input files need to be created on disk, a process is launched, possibly reading all kinds of configuration and parameter files. The process writes more files to disk, which we later need to open again to extract (in the worst case) just a single number. The overhead might be irrelevant for sufficiently slow engines, but for a very fast force field this overhead can easily become the performance bottleneck.

Starting with the AMS2019.3 release, the AMS driver implements a special task, in which the running process listens for calculation requests on a named pipe (FIFO) and communicates the results of the calculations back on another pipe. This avoids the overhead of starting processes and eliminates all file based I/O. You can find more information about the pipe interface in the AMS driver in the `corresponding part of the documentation <../../AMS/Input_Output.html#pipe-interface>`_. In PLAMS the |AMSWorker| class is used to represent this running AMS driver process. The |AMSWorker| class handles all communication with the process and hides the technical details of underlying `communication protocol <../../AMS/Pipe_protocol.html>`_.

Consider the following short PLAMS script, that calculates and prints the total GFN1-xTB energy for all molecules found in a folder full of xyz-files. Using the regular |AMSJob|, this can be written as:

.. code-block:: python

      molecules = read_molecules('folder/with/xyz/files')

      sett = Settings()
      sett.input.ams.Task = 'SinglePoint'
      sett.input.dftb.Model = 'GFN1-xTB'

      for name, mol in molecules.items():
          results = AMSJob(name=name, molecule=mol, settings=sett).run()
          print('Energy of {} = {}'.format(name, results.get_energy()))

In order to switch this script over to using the |AMSWorker|, we need to make only a couple of changes:

.. code-block:: python

      molecules = read_molecules('folder/with/xyz/files')

      sett = Settings()
      sett.input.dftb.Model = 'GFN1-xTB'

      with AMSWorker(sett) as worker:
          for name, mol in molecules.items():
              results = worker.SinglePoint(name, mol)
              print('Energy of {} = {}'.format(name, results.get_energy()))

With the first |AMSJob| based version, both the |Settings| instance and the |Molecule| instance were passed into the constructor of the |AMSJob|, while the |AMSWorker| constructor only accepts the |Settings| instance. The |Molecule| instance is only later passed into the :meth:`SinglePoint<AMSWorker.SinglePoint>` method. This shows the basic usage of the |AMSWorker| class: create it once, supplying the desired |Settings|, and use these fixed settings for calculations on multiple molecules. It is *not* possible to change the |Settings| on an already running |AMSWorker| instance. If you have to switch |Settings|, you need to create a new |AMSWorker| with the new settings. It therefore only makes sense to use the |AMSWorker| if one has to do calculations on many molecules *using the same settings*.

Note that when using the |AMSWorker| the type of ``results`` in the above example is not actually |AMSResults| anymore: the call to :meth:`SinglePoint<AMSWorker.SinglePoint>` returns an instance of |AMSWorkerResults|, which only implements a small subset of the methods available in the full |AMSResults|. This is the concession we have to make for using |AMSWorker| instead of |AMSJob|: after all the |AMSJob| class has many methods to extract arbitrary data from the result files of an AMS calculation. Since none of these files exist when directly communicating with the AMS process over a pipe, the |AMSWorkerResults| class supports none of these methods.

Given these restrictions we recommend that users first try the traditional route of running the AMS driver via the |AMSJob| class, and only switch to the |AMSWorker| alternative if they observe a significant slowdown due to the startup and I/O cost. The overhead is likely only relevant for simple tasks (single points, geometry optimizations) using rather fast engines such as semi-empirical methods and force fields.

In case the worker process fails to start up or terminates unexpectedly, an |AMSWorkerError| exception will be raised. The standard output and standard error output from the failed worker process is stored in the ``stdout`` and ``stderr`` attributes in |AMSWorkerError|.
If an |AMSWorkerError| or |AMSPipeRuntimeError| exception occurs during :meth:`SinglePoint<AMSWorker.SinglePoint>`, it will be internally caught and stored in the ``error`` attribute of the returned |AMSWorkerResults| object for further inspection.
These two types of exceptions are typically related to the calculation being performed (the combination of the |Molecule| and |Settings|), so they are not allowed to propagate out of :meth:`SinglePoint<AMSWorker.SinglePoint>` to match the behavior of |AMSJob| in similar situations.
However, other types of exceptions derived from |AMSPipeError| may also occur in |AMSWorker|.
These correspond to other errors defined by the pipe protocol and will propagate normally, because they represent programming and logic errors, protocol incompatibilities, or unsupported features.
In any case, |AMSWorker| will be ready to handle another call to :meth:`SinglePoint<AMSWorker.SinglePoint>` after an error.

PLAMS also provides the |AMSWorkerPool| class, which represents a pool of running |AMSWorker| instances, which dynamically pick tasks from a queue of calculations to be performed. This is useful for workflows that require the execution of many trivially parallel simple tasks. Using the |AMSWorkerPool| we could write the above example as:

.. code-block:: python

      molecules = read_molecules('folder/with/xyz/files')

      sett = Settings()
      sett.input.dftb.Model = 'GFN1-xTB'
      sett.runscript.nproc = 1 # every worker is a serial process now

      with AMSWorkerPool(sett, num_workers=4) as pool:
           results = pool.SinglePoints(molecules.items())
      for r in results:
         print('Energy of {} = {}'.format(r.name, r.get_energy()))



AMSWorker API
~~~~~~~~~~~~~

.. autoclass:: AMSWorker
    :exclude-members: __init__, __weakref__
    :no-private-members:



AMSWorkerResults API
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: AMSWorkerResults
    :exclude-members: __init__, __weakref__
    :no-private-members:



AMSWorkerPool API
~~~~~~~~~~~~~~~~~

.. autoclass:: AMSWorkerPool
    :exclude-members: __init__, __weakref__
    :no-private-members:


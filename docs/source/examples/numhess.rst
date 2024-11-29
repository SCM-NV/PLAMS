Numerical Hessian
--------------------

This module implements a simple scheme for calculating a numerical Hessian matrix.
We define a new job type ``NumHessJob`` by extending |MultiJob|.
The constructor (``__init__``) of this new job accepts several new arguments and simply stores them.
These new arguments define the initial |Molecule|, the type of job used for single point calculations (``jobtype``), the size and unit of displacement step and the way of extracting gradients from single point results.

Then the |prerun| method takes the given |Molecule| and creates multiple copies of it, each one with one atom displaced along one axis.
For each of these molecules an instance of single point job is created and stored in the ``children`` dictionary.
Settings of ``NumHessJob`` are directly passed to children jobs, so creating a ``NumHessJob`` strongly resembles creating a regular single point job.

The dedicated |Results| subclass for ``NumHessJob`` takes care of extracting the Hessian from results of all single point jobs.
The returned Hessian can optionally be mass weighted.

The source code of the whole module with both aforementioned classes:

.. literalinclude:: ../../../recipes/numhess.py

An example usage::

    mol = Molecule('methanol.xyz')

    s = Settings()
    s.input.ams.task = 'SinglePoint'
    s.input.ams.Properties.Gradients = 'Yes'
    s.input.adf.basis.type = 'DZP'
    s.input.adf.symmetry = 'NOSYM'
    s.input.adf.xc.gga = 'PW91'
    s.runscript.nproc = 1

    j = NumHessJob(name='test', molecule=mol, settings=s, jobtype=AMSJob,
                   gradient = lambda x: x.get_gradients().reshape(-1))
    r = j.run(jobrunner=JobRunner(parallel=True, maxjobs=8))
    print(r.get_hessian(mass_weighted=True))


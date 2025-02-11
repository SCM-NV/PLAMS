Worked Example
--------------

Logging in PLAMS
~~~~~~~~~~~~~~~~

PLAMS has built-in logging which aims to simplify tracking the progress and status of jobs. This consists of progress logging to stdout and a logfile, and writing job summaries to CSV files. Each of these is explained below.

Progress Logger
~~~~~~~~~~~~~~~

PLAMS writes job progress to stdout and a plain text logfile. The location of this logfile is determined by the working directory of the default job manager, and is called ``logfile``.

Users can also write logs to the same locations using the ``log`` function. This takes a ``level`` argument. By convention in PLAMS, the level should be between 0-7, with 0 being the most and 7 the least important logging.

The level of logging that is written to stdout and the logfile can be changed through the ``config.LogSettings``.

.. code:: ipython3

   from scm.plams import Settings, AMSJob, from_smiles, log, config, init

   # this line is not required in AMS2025+
   init()

   counter = 0


   def get_test_job():
       global counter
       s = Settings()
       s.input.ams.Task = "SinglePoint"
       s.input.dftb
       counter += 1
       return AMSJob(name=f"test{counter}", molecule=from_smiles("C"), settings=s)

::

   PLAMS working folder: /path/plams/examples/Logging/plams_workdir.004

.. code:: ipython3

   config.log.stdout = 3
   config.log.file = 5
   config.jobmanager.hashing = None  # Force PLAMS to re-run identical test jobs

.. code:: ipython3

   job = get_test_job()
   job.run()
   log("Test job finished", 5)

::

   [10.02|15:33:22] JOB test1 STARTED
   [10.02|15:33:22] JOB test1 RUNNING
   [10.02|15:33:23] JOB test1 FINISHED
   [10.02|15:33:23] JOB test1 SUCCESSFUL

.. code:: ipython3

   with open(config.default_jobmanager.logfile, "r") as f:
       print(f.read())

::

   [10.02|15:33:22] JOB test1 STARTED
   [10.02|15:33:22] Starting test1.prerun()
   [10.02|15:33:22] test1.prerun() finished
   [10.02|15:33:22] JOB test1 RUNNING
   [10.02|15:33:22] Executing test1.run
   [10.02|15:33:23] Execution of test1.run finished with returncode 0
   [10.02|15:33:23] JOB test1 FINISHED
   [10.02|15:33:23] Starting test1.postrun()
   [10.02|15:33:23] test1.postrun() finished
   [10.02|15:33:23] JOB test1 SUCCESSFUL
   [10.02|15:33:23] Test job finished

Note that the logs from an AMS calculation can also be forwarded to the progress logs using the ``watch = True`` flag.

.. code:: ipython3

   job = get_test_job()
   job.run(watch=True);

::

   [10.02|15:33:23] JOB test2 STARTED
   [10.02|15:33:23] JOB test2 RUNNING
   [10.02|15:33:23] test2: AMS 2024.207  RunTime: Feb10-2025 15:33:23  ShM Nodes: 1  Procs: 6
   [10.02|15:33:24] test2: DFTB: SCC cycle
   [10.02|15:33:24] test2: cyc=  1 err=1.1E+00 method=1 nvec= 1 mix=0.075 e=    0.0000
   [10.02|15:33:24] test2: cyc=  2 err=1.1E+00 method=1 nvec= 1 mix=0.154 e=    0.0000
   [10.02|15:33:24] test2: cyc=  3 err=8.9E-01 method=1 nvec= 2 mix=0.201 e=    0.0000
   [10.02|15:33:24] test2: cyc=  4 err=1.7E-02 method=1 nvec= 3 mix=0.207 e=    0.0000
   [10.02|15:33:24] test2: cyc=  5 err=6.8E-03 method=1 nvec= 4 mix=0.213 e=    0.0000
   [10.02|15:33:24] test2: cyc=  6 err=2.6E-03 method=1 nvec= 5 mix=0.219 e=    0.0000
   [10.02|15:33:24] test2: cyc=  7 err=7.2E-05 method=1 nvec= 6 mix=0.226 e=    0.0000
   [10.02|15:33:24] test2: cyc=  8 err=6.8E-05 method=1 nvec= 1 mix=0.233 e=    0.0000
   [10.02|15:33:24] test2: cyc=  9 err=4.2E-05 method=1 nvec= 2 mix=0.240 e=    0.0000
   [10.02|15:33:24] test2: cyc= 10 err=6.2E-07 method=1 nvec= 3 mix=0.247 e=    0.0000
   [10.02|15:33:24] test2: cyc= 11 err=5.8E-08 method=1 nvec= 3 mix=0.254 e=    0.0000
   [10.02|15:33:24] test2: cyc= 12 err=3.6E-08 method=1 nvec= 4 mix=0.262 e=    0.0000
   [10.02|15:33:24] test2: cyc= 13 err=9.0E-11 method=1 nvec= 4 mix=0.270 e=    0.0000
   [10.02|15:33:24] test2: SCC cycle converged!
   [10.02|15:33:24] test2: NORMAL TERMINATION
   [10.02|15:33:24] JOB test2 FINISHED
   [10.02|15:33:24] JOB test2 SUCCESSFUL

Job Summary Logger
~~~~~~~~~~~~~~~~~~

For AMS2025+, PLAMS also writes summaries of jobs to a CSV file, the location of which by default is also determined by the job manager. It is called ``job_logfile.csv``.

.. code:: ipython3

   from scm.plams import MultiJob


   jobs = [get_test_job() for _ in range(3)]
   jobs[2].settings.input.ams.Task = "Not a task!"

   for job in jobs:
       job.run()

::

   [10.02|15:33:24] JOB test3 STARTED
   [10.02|15:33:24] JOB test3 RUNNING
   [10.02|15:33:25] JOB test3 FINISHED
   [10.02|15:33:25] JOB test3 SUCCESSFUL
   [10.02|15:33:25] JOB test4 STARTED
   [10.02|15:33:25] JOB test4 RUNNING
   [10.02|15:33:26] JOB test4 FINISHED
   [10.02|15:33:26] JOB test4 SUCCESSFUL
   [10.02|15:33:26] JOB test5 STARTED
   [10.02|15:33:26] JOB test5 RUNNING
   [10.02|15:33:34] WARNING: Job test5 finished with nonzero return code
   [10.02|15:33:34] WARNING: Main KF file ams.rkf not present in /path/plams/examples/Logging/plams_workdir.004/test5
   ... (PLAMS log lines truncated) ...
   [10.02|15:33:34] File ams.rkf not present in /path/plams/examples/Logging/plams_workdir.004/test5
   [10.02|15:33:34] Error message for job test5 was:
       Input error: value "Not a task!" found in line 1 for multiple choice key "Task" is not an allowed choice
   [10.02|15:33:34] File ams.rkf not present in /path/plams/examples/Logging/plams_workdir.004/test5
   [10.02|15:33:34] File ams.rkf not present in /path/plams/examples/Logging/plams_workdir.004/test5

These CSVs give overall information on the status of all jobs run by a given job manager.

.. code:: ipython3

   import csv

   try:
       with open(config.default_jobmanager.job_logger.logfile, newline="") as csvfile:
           reader = csv.DictReader(csvfile)
           for row in reader:
               print(f"{row['job_name']} {row['job_status']}: {row['job_get_errormsg']}")
   except AttributeError:
       pass

::

   test1 successful: 
   test2 successful: 
   test3 successful: 
   test4 successful: 
   test5 crashed: Input error: value "Not a task!" found in line 1 for multiple choice key "Task" is not an allowed choice

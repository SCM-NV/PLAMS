Worked Example
--------------

Logging in PLAMS
~~~~~~~~~~~~~~~~

PLAMS has built-in logging which aims to simplify tracking the progress
and status of jobs. This consists of progress logging to stdout and a
logfile, and writing job summaries to CSV files. Each of these is
explained below.

Progress Logger
~~~~~~~~~~~~~~~

PLAMS writes job progress to stdout and a plain text logfile. The
location of this logfile is determined by the working directory of the
default job manager, and is called ``logfile``.

Users can also write logs to the same locations using the ``log``
function. This takes a ``level`` argument. By convention in PLAMS, the
level should be between 0-7, with 0 being the most and 7 the least
important logging.

The level of logging that is written to stdout and the logfile can be
changed through the ``config.LogSettings``.

.. code:: ipython3

    from scm.plams import Settings, AMSJob, from_smiles, log, config
    
    counter = 0
    
    
    def get_test_job():
        global counter
        s = Settings()
        s.input.ams.Task = "SinglePoint"
        s.input.dftb
        counter += 1
        return AMSJob(name=f"test{counter}", molecule=from_smiles("C"), settings=s)

.. code:: ipython3

    config.log.stdout = 3
    config.log.file = 5
    config.jobmanager.hashing = None  # Force PLAMS to re-run identical test jobs

.. code:: ipython3

    job = get_test_job()
    job.run()
    log("Test job finished", 5)


.. parsed-literal::

    [29.11|09:41:26] JOB test1 STARTED
    [29.11|09:41:26] JOB test1 RUNNING
    [29.11|09:41:27] JOB test1 FINISHED
    [29.11|09:41:28] JOB test1 SUCCESSFUL


.. code:: ipython3

    with open(config.default_jobmanager.logfile, "r") as f:
        print(f.read())


.. parsed-literal::

    [29.11|09:41:26] Starting test1.prerun()
    [29.11|09:41:26] test1.prerun() finished
    [29.11|09:41:26] JOB test1 RUNNING
    [29.11|09:41:26] Executing test1.run
    [29.11|09:41:27] Execution of test1.run finished with returncode 0
    [29.11|09:41:27] JOB test1 FINISHED
    [29.11|09:41:27] Starting test1.postrun()
    [29.11|09:41:27] test1.postrun() finished
    [29.11|09:41:28] JOB test1 SUCCESSFUL
    [29.11|09:41:28] Test job finished
    


Note that the logs from an AMS calculation can also be forwarded to the
progress logs using the ``watch = True`` flag.

.. code:: ipython3

    job = get_test_job()
    job.run(watch=True);


.. parsed-literal::

    [29.11|09:41:28] JOB test2 STARTED
    [29.11|09:41:28] JOB test2 RUNNING
    [29.11|09:41:28] test2: AMS 2024.206  RunTime: Nov29-2024 09:41:28  ShM Nodes: 1  Procs: 6
    [29.11|09:41:28] test2: DFTB: SCC cycle
    [29.11|09:41:28] test2: cyc=  1 err=1.1E+00 method=1 nvec= 1 mix=0.075 e=    0.0000
    [29.11|09:41:28] test2: cyc=  2 err=1.1E+00 method=1 nvec= 1 mix=0.154 e=    0.0000
    [29.11|09:41:28] test2: cyc=  3 err=8.9E-01 method=1 nvec= 2 mix=0.201 e=    0.0000
    [29.11|09:41:28] test2: cyc=  4 err=1.7E-02 method=1 nvec= 3 mix=0.207 e=    0.0000
    [29.11|09:41:28] test2: cyc=  5 err=6.8E-03 method=1 nvec= 4 mix=0.213 e=    0.0000
    [29.11|09:41:28] test2: cyc=  6 err=2.6E-03 method=1 nvec= 5 mix=0.219 e=    0.0000
    [29.11|09:41:28] test2: cyc=  7 err=7.2E-05 method=1 nvec= 6 mix=0.226 e=    0.0000
    [29.11|09:41:28] test2: cyc=  8 err=6.8E-05 method=1 nvec= 1 mix=0.233 e=    0.0000
    [29.11|09:41:28] test2: cyc=  9 err=4.2E-05 method=1 nvec= 2 mix=0.240 e=    0.0000
    [29.11|09:41:28] test2: cyc= 10 err=6.2E-07 method=1 nvec= 3 mix=0.247 e=    0.0000
    [29.11|09:41:28] test2: cyc= 11 err=5.8E-08 method=1 nvec= 3 mix=0.254 e=    0.0000
    [29.11|09:41:28] test2: cyc= 12 err=3.6E-08 method=1 nvec= 4 mix=0.262 e=    0.0000
    [29.11|09:41:28] test2: cyc= 13 err=9.0E-11 method=1 nvec= 4 mix=0.270 e=    0.0000
    [29.11|09:41:28] test2: SCC cycle converged!
    [29.11|09:41:28] test2: NORMAL TERMINATION
    [29.11|09:41:28] JOB test2 FINISHED
    [29.11|09:41:28] JOB test2 SUCCESSFUL


Job Summary Logger
~~~~~~~~~~~~~~~~~~

PLAMS also writes summaries of jobs to a CSV file, the location of which
by default is also determined by the job manager. It is called
``job_logfile.csv``.

.. code:: ipython3

    from scm.plams import MultiJob
    
    
    jobs = [get_test_job() for _ in range(3)]
    jobs[2].settings.input.ams.Task = "Not a task!"
    
    for job in jobs:
        job.run()


.. parsed-literal::

    [29.11|09:41:28] JOB test3 STARTED
    [29.11|09:41:28] JOB test3 RUNNING
    [29.11|09:41:29] JOB test3 FINISHED
    [29.11|09:41:29] JOB test3 SUCCESSFUL
    [29.11|09:41:29] JOB test4 STARTED
    [29.11|09:41:29] JOB test4 RUNNING
    [29.11|09:41:30] JOB test4 FINISHED
    [29.11|09:41:30] JOB test4 SUCCESSFUL
    [29.11|09:41:30] JOB test5 STARTED
    [29.11|09:41:30] JOB test5 RUNNING
    [29.11|09:41:38] WARNING: Job test5 finished with nonzero return code
    [29.11|09:41:38] WARNING: Main KF file ams.rkf not present in /path/plams/examples/Logging/plams_workdir/test5
    [29.11|09:41:38] JOB test5 CRASHED
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] Obtaining results of test5 failed. Returned value is None
    [29.11|09:41:38] Obtaining results of test5 successful. However, no guarantee that they make sense
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] Obtaining results of test5 failed. Returned value is None
    [29.11|09:41:38] Obtaining results of test5 successful. However, no guarantee that they make sense
    [29.11|09:41:38] Could not read termination status from file None
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] Obtaining results of test5 failed. Returned value is None
    [29.11|09:41:38] Obtaining results of test5 successful. However, no guarantee that they make sense
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] Obtaining results of test5 failed. Returned value is None
    [29.11|09:41:38] Obtaining results of test5 successful. However, no guarantee that they make sense
    [29.11|09:41:38] Could not read termination status from file None
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] Obtaining results of test5 failed. Returned value is None
    [29.11|09:41:38] Obtaining results of test5 successful. However, no guarantee that they make sense
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] WARNING: Trying to obtain results of crashed or failed job test5
    [29.11|09:41:38] Obtaining results of test5 successful. However, no guarantee that they make sense
    [29.11|09:41:38] Obtaining results of test5 successful. However, no guarantee that they make sense


These CSVs give overall information on the status of all jobs run by a
given job manager.

.. code:: ipython3

    try:
        import pandas as pd
    
        df = pd.read_csv(config.default_jobmanager.job_logger.logfile)
        print(df[["job_name", "job_status", "job_ok", "job_get_errormsg"]])
    except ImportError:
        pass


.. parsed-literal::

      job_name  job_status  job_ok  \
    0    test1  successful    True   
    1    test2  successful    True   
    2    test3  successful    True   
    3    test4  successful    True   
    4    test5     crashed   False   
    
                                        job_get_errormsg  
    0                                                NaN  
    1                                                NaN  
    2                                                NaN  
    3                                                NaN  
    4  Input error: value "Not a task!" found in line...  


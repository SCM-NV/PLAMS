#!/usr/bin/env amspython
# coding: utf-8

# ## Logging in PLAMS

# PLAMS has built-in logging which aims to simplify tracking the progress and status of jobs. This consists of progress logging to stdout and a logfile, and writing job summaries to CSV files. Each of these is explained below.

# ### Progress Logger

# PLAMS writes job progress to stdout and a plain text logfile. The location of this logfile is determined by the working directory of the default job manager, and is called `logfile`.
#
# Users can also write logs to the same locations using the `log` function. This takes a `level` argument. By convention in PLAMS, the level should be between 0-7, with 0 being the most and 7 the least important logging.
#
# The level of logging that is written to stdout and the logfile can be changed through the `config.LogSettings`.

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


config.log.stdout = 3
config.log.file = 5
config.jobmanager.hashing = None  # Force PLAMS to re-run identical test jobs


job = get_test_job()
job.run()
log("Test job finished", 5)


with open(config.default_jobmanager.logfile, "r") as f:
    print(f.read())


# Note that the logs from an AMS calculation can also be forwarded to the progress logs using the `watch = True` flag.

job = get_test_job()
job.run(watch=True)


# ### Job Summary Logger

# For AMS2025+, PLAMS also writes summaries of jobs to a CSV file, the location of which by default is also determined by the job manager. It is called `job_logfile.csv`.

from scm.plams import MultiJob


jobs = [get_test_job() for _ in range(3)]
jobs[2].settings.input.ams.Task = "Not a task!"

for job in jobs:
    job.run()


# These CSVs give overall information on the status of all jobs run by a given job manager.

import csv

try:
    with open(config.default_jobmanager.job_logger.logfile, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            print(f"{row['job_name']} {row['job_status']}: {row['job_get_errormsg']}")
except AttributeError:
    pass

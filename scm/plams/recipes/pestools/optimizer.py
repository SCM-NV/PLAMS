#!/usr/bin/env python

import datetime
import os
import threading
import time

from scm.plams.core.errors import FileError, PlamsError
from scm.plams.core.functions import log
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.amsworker import AMSWorkerPool
from scm.plams.core.enums import JobStatus

__all__ = ["Optimizer"]


class Optimizer:
    """
    Machine that can optimize the geometry of a set of ACE molecules
    """

    def __init__(self, nproc=1, jobrunner=None):
        """
        Initiates an instance of the Optimizer class

        * ``nproc``           -- Number of processors to use in total for the optimizations
        """
        # The general job settings
        self.settings = Settings()
        self.settings.runscript.nproc = nproc

        # The settings for geometry optimizations
        self.go_settings = Settings()
        self.go_settings.input.ams.GeometryOptimization.Method = "Quasi-Newton"
        self.go_settings.input.ams.GeometryOptimization.CoordinateType = "Auto"  # Cartesian?
        # self.method = 'Quasi-Newton'
        # self.coordinate_type = 'Auto' # Cartesian?

        # Set the jobrunner
        if jobrunner is None:
            parallel = False
            if nproc > 1:
                parallel = True
            jobrunner = JobRunner(parallel=parallel, maxjobs=nproc)
        self.jobrunner = jobrunner

        # Wether or not to keep the files after running
        self.keep = None  # 'all' is the PLAMS default
        self.use_pipe = True
        self.last_timings = None
        self.watch = False

        # Private variables
        self.path = None
        self.nsteps = None

    def set_cleaning_preference(self, keep=None):
        """
        Sets the preferences for the files saved to disc

        * ``keep`` -- None or 'all'.
         The former will clean PLAMS folders of all files after runtime,
         the latter will leave all files as they are.
         See PLAMS for more details.
        """
        if keep is not None:
            self.use_pipe = False
        self.keep = keep

    def keep_ams_running(self, keep=True):
        """
        Technical setting for a set of (serial) geometry optimizations

        * ``keep`` -- Each of the parallel processes will most likely be passed multiple geometries.
         If `keep=True`, AMS will be kept running in the background, and geometries will be passed via a pipe.
         If `keep=False`, each geometry optimization will start a new AMS process.
        """
        self.use_pipe = keep

    def optimize_geometries(self, molecules, engine_settings=None, name="go"):
        """
        Optimize a set of molecules

        * ``engine_settings`` -- PLAMS Settings object:
                    engine_settings = Settings()
                    engine_settings.DFTB.Model = 'GFN1-xTB'
        """
        if not self.use_pipe or len(molecules) == 1:
            molecules, energies, resultlist = self._perform_tasks_separately(molecules, engine_settings, name)
        else:
            molecules, energies, resultlist = self._perform_tasks_using_pipe(molecules, engine_settings, name)
        remaining_energies = [en for en in energies if not en is None]
        log("%s out of %s optimizations succeeded" % (len(remaining_energies), len(energies)))
        return molecules, energies

    def compute_energies(self, molecules, engine_settings=None, name="sp"):
        """
        Compute the energies for a set of molecules

        * ``engine_settings`` -- PLAMS Settings object:
                    engine_settings = Settings()
                    engine_settings.DFTB.Model = 'GFN1-xTB'
        """
        if not self.use_pipe or len(molecules) == 1:
            molecules, energies, resultlist = self._perform_tasks_separately(
                molecules, engine_settings, name, "SinglePoint"
            )
        else:
            molecules, energies, resultlist = self._perform_tasks_using_pipe(
                molecules, engine_settings, name, "SinglePoint"
            )
        remaining_energies = [en for en in energies if not en is None]
        log("%s out of %s single point calculations succeeded" % (len(remaining_energies), len(energies)))
        return molecules, energies

    def _perform_tasks_separately(self, molecules, engine_settings=None, name="go", task="GeometryOptimization"):
        """
        Optimize a set of molecules
        """
        # Set the settings object
        if task == "GeometryOptimization":
            settings = self.get_go_settings()
        elif task == "SinglePoint":
            settings = self.get_sp_settings()
        else:
            raise PlamsError("Task not implemented: %s" % (task))

        # Add the engine settings
        if engine_settings is None:
            engine_settings = Settings()
            engine_settings.ForceField.Type = "UFF"
        engine_name = list(engine_settings.keys())[0]
        settings.input[engine_name] = engine_settings[engine_name]

        # Adjust the number of processes if thre is only one job!
        if len(molecules) > 0:
            settings.runscript.nproc = int(self.jobrunner.maxjobs / min(self.jobrunner.maxjobs, len(molecules)))

        # Run jobs
        taskmap = {"GeometryOptimization": "geometry optimization", "SinglePoint": "single point"}
        starttime = time.time()
        resultlist = []
        if self.watch:
            # Set up and start the progress monitor.
            done_event = threading.Event()
            args = (resultlist, len(molecules), done_event, 60)
            pmt = threading.Thread(target=progress_monitor, args=args)
            log(f"Running {len(molecules)} {taskmap[task]}s")
            pmt.start()

        # Run jobs
        for i, mol in enumerate(molecules):
            job = AMSJob(name="%s%i" % (name, i), molecule=mol, settings=settings)
            resultlist.append(job.run(jobrunner=self.jobrunner))

        # Read results
        if task == "GeometryOptimization":
            optimized_geometries, energies = self._read_results_go(resultlist)
        elif task == "SinglePoint":
            optimized_geometries, energies = self._read_results_sp(resultlist)
        endtime = time.time()
        self.last_timings = endtime - starttime

        # Stop the progreess monitor (all the jobs have finished)
        if self.watch:
            done_event.set()
            pmt.join()
            time_taken = datetime.timedelta(seconds=round(time.time() - starttime))
            txt = f"All {len(molecules)} {taskmap[task]}s done! Time taken: {time_taken}s"
            log(txt)

        # Place the new coordinates into the molecules
        for i, (crd, en) in enumerate(zip(optimized_geometries, energies)):
            if en is not None:
                if crd is not None:
                    molecules[i].from_array(crd)
            molecules[i].properties.energy = en

        # Remove files
        keep = self.keep
        # if len(molecules) == 1 or task == 'SinglePoint' : keep = 'all'
        if len(molecules) == 1:
            keep = "all"
        for i, r in enumerate(resultlist):
            r._clean(keep)
            if keep is None:
                r.job.jobmanager.remove_job(r.job)

        return molecules, energies, resultlist

    def _perform_tasks_using_pipe(self, molecules, engine_settings=None, name="go", task="GeometryOptimization"):
        """
        Optimize a set of ACE molecules
        """
        if len(molecules) == 0:
            return [], [], []

        # Add the engine settings
        settings = self.settings.copy()
        if engine_settings is None:
            engine_settings = Settings()
            engine_settings.ForceField.Type = "UFF"
        engine_name = list(engine_settings.keys())[0]
        settings.input[engine_name] = engine_settings[engine_name]

        # Lift out geometry optimization settings (now only 'Type')
        kwargs = {}
        kwargs["quiet"] = False

        for key in self.go_settings.input.ams.keys():
            if key.lower() == "geometryoptimization":
                go_sett = self.go_settings.input.ams[key]
                if "convergence" in go_sett:
                    convkeys = [k.lower() for k in go_sett.convergence.keys()]
                    if "Energy" in convkeys:
                        kwargs["convenergy"] = go_sett.Convergence.Energy
                    if "Gradients" in convkeys:
                        kwargs["convgradients"] = go_sett.Convergence.Gradients
                if "method" in go_sett:
                    kwargs["method"] = go_sett.method
                if "maxiterations" in go_sett:
                    kwargs["maxiterations"] = go_sett.maxiterations
                if "coordinatetype" in go_sett:
                    kwargs["coordinatetype"] = go_sett.coordinatetype
                if "pretendconverged" in go_sett:
                    kwargs["pretendconverged"] = go_sett.pretendconverged
        # Run the tasks
        molecule_list = []
        for i, mol in enumerate(molecules):
            molecule_list.append(("%s%i" % (name, i), mol.copy(), kwargs))
        # with AMSWorkerPool(settings, jobrunner=self.jobrunner) as job_pool :
        # maxjobs = self.jobrunner.maxjobs if self.jobrunner.parallel else 0
        # job_pool = AMSWorkerPool(settings, jobrunner=self.jobrunner)
        from scm.plams import config

        if not "default_jobmanager" in config:
            raise PlamsError("No default jobmanager found.")
        jobmanager = config.default_jobmanager
        settings.runscript.nproc = int(self.jobrunner.maxjobs / min(self.jobrunner.maxjobs, len(molecules)))
        nproc = min(self.jobrunner.maxjobs, len(molecules))
        job_pool = AMSWorkerPool(settings, nproc, workerdir_root=jobmanager.workdir)
        keepdir = False
        try:
            if task == "GeometryOptimization":
                resultlist = job_pool.GeometryOptimizations(molecule_list, watch=self.watch)
            elif task == "SinglePoint":
                resultlist = job_pool.SinglePoints(molecule_list, watch=self.watch)
            else:
                raise PlamsError("Task not implemented: %s" % (task))
        except Exception:
            keepdir = True
        finally:
            for worker in job_pool.workers:
                out, err = worker.stop(keepdir)
                # I am doing this because multiprocessing (use by ACE) also uses tempfile, and cleans the whole tempdir at the end
                worker._finalize.__call__()

        # Open the error file
        errfile = open(os.path.join(jobmanager.workdir, "amsworkerpool.err"), "a")

        # Read the results
        optimized_geometries = []
        energies = []
        for i, r in enumerate(resultlist):
            msg = "Worker failed"
            jobname = molecule_list[i][0]
            if r is not None:
                msg = r.get_errormsg()
                jobname = r.name
            if msg is not None:
                # print ('Error: ',i,msg)
                errfile.write("%-10s: %s\n" % (jobname, msg))
                optimized_geometries.append(None)
                energies.append(None)
                continue
            energy = r.get_energy()  # * Units.conversion_ratio('Hartree','kcal/mol')
            coords = r.get_main_molecule().as_array()
            optimized_geometries.append(coords)
            energies.append(energy)

        errfile.close()

        # Place the new coordinates into the molecules
        for i, (crd, en) in enumerate(zip(optimized_geometries, energies)):
            if crd is not None:
                molecules[i].from_array(crd)
            molecules[i].properties.energy = en

        return molecules, energies, resultlist

    def get_go_settings(self):
        """
        Get the settings for a geometry optimization job
        """
        # Set the settings object
        settings = self.settings.copy()
        settings.input.ams.Task = "GeometryOptimization"
        settings.update(self.go_settings)
        # settings.input.ams.GeometryOptimization.Method = self.method
        # settings.input.ams.GeometryOptimization.CoordinateType = self.coordinate_type
        return settings

    def get_sp_settings(self):
        """
        Get the settings for a geometry optimization job
        """
        # Set the settings object
        settings = self.settings.copy()
        settings.input.ams.Task = "SinglePoint"
        return settings

    def _read_results_go(self, resultlist):
        """
        Read the geometries and energies from the results objects
        """
        optimized_geometries = []
        energies = []
        for i, r in enumerate(resultlist):
            try:
                # I do it like this to make sure we hit the guardian, I think
                nsteps = r.readrkf("History", "nEntries")
            except FileError:
                pass
            # If the geometry was copied from another directory, then it is not new so we do not need it.
            if r.job.status == JobStatus.FAILED or not os.path.isfile(os.path.join(r.job.path, "ams.rkf")):
                optimized_geometries.append(None)
                energies.append(None)
                continue
            if nsteps is None:
                # Apparently it can still go wrong!
                nsteps = r.readrkf("History", "nEntries")
            if nsteps > 0:
                energy = r.get_property_at_step(nsteps, "Energy")  # * Units.conversion_ratio('Hartree','kcal/mol')
            else:
                try:
                    energy = r.get_energy()  # * Units.conversion_ratio('Hartree','kcal/mol')
                except FileError:
                    energy = None

            coords = r.get_main_molecule().as_array()
            optimized_geometries.append(coords)
            energies.append(energy)
            # The first time, set some variables
            if self.path is None:
                self.path = os.path.dirname(r.job.path)
                self.nsteps = nsteps
        return optimized_geometries, energies

    def _read_results_sp(self, resultlist):
        """
        Read the geometries and energies from the results objects
        """
        energies = []
        for i, r in enumerate(resultlist):
            try:
                energy = r.get_energy()  # * Units.conversion_ratio('Hartree','kcal/mol')
            except FileError:
                pass
            if r.job.status == JobStatus.FAILED or not os.path.isfile(os.path.join(r.job.path, "ams.rkf")):
                energies.append(None)
                continue
            energy = r.get_energy()  # * Units.conversion_ratio('Hartree','kcal/mol')
            energies.append(energy)
        return [None for en in energies], energies


def progress_monitor(results, num_jobs, event, t):
    """
    Provides regular updates on the progress of geometry optimizations in a different thread
    """
    num_done = 0
    starttime = time.time()
    width = len(str(num_jobs))
    while True:
        if event.wait(timeout=t):
            return
        try:
            num_done = len(
                [
                    True
                    for r in results
                    if r.job.status
                    not in [JobStatus.CREATED, JobStatus.STARTED, JobStatus.REGISTERED, JobStatus.RUNNING]
                ]
            )
        except RuntimeError:
            pass
        percent_done = 100.0 * num_done / num_jobs
        if percent_done > 5.0:
            dt = time.time() - starttime
            dtrem = dt / percent_done * (100 - percent_done)
            dtrem = datetime.timedelta(seconds=round(dtrem))
            trem = f", {dtrem}s remaining"
        else:
            trem = ""
        log(f"{str(num_done).rjust(width)} / {num_jobs} jobs finished:{percent_done:5.1f}%{trem}")

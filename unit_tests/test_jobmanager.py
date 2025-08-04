import pytest
import os
import uuid
from pathlib import Path

from scm.plams.core.basejob import MultiJob
from scm.plams.core.jobmanager import JobManager
from scm.plams.core.settings import JobManagerSettings
from scm.plams.core.errors import PlamsError
from scm.plams.core.functions import jobs_in_directory
from .test_basejob import DummySingleJob


class TestJobManager:

    def test_lazy_workdir(self):
        # Given job manager
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)

        # When first initialised
        # Then workdir does not exist
        assert not os.path.exists(job_manager._workdir)

        # When access workdir for the first time
        # Then workdir is created
        workdir = job_manager.workdir
        assert os.path.exists(workdir)
        assert os.path.exists(job_manager._workdir)

        # When access subsequent time
        # Then same workdir is returned
        assert job_manager.workdir == workdir

        os.rmdir(job_manager.workdir)

    def test_load_and_clean_do_not_create_workdir(self):
        # Given job manager
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)

        # When load job
        job = DummySingleJob()
        job.run()
        job.results.wait()
        job_manager.load_job(f"{job.path}/{job.name}.dill")

        # Then workdir not created
        assert not os.path.exists(job_manager._workdir)

        # When clean the jobmanager
        job_manager._clean()

        # Then workdir not created
        assert not os.path.exists(job_manager._workdir)

    def test_load_legacy_job_succeeds(self):
        # Given job run before additional properties were added
        job1 = DummySingleJob()
        job1.run()
        job1.results.wait()
        status = job1.status
        delattr(job1, "_status")
        delattr(job1, "_status_log")
        job1.__dict__["status"] = status
        job1.pickle()

        # When call load
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)
        job2 = job_manager.load_job(f"{job1.path}/{job1.name}.dill")

        # Then job loaded successfully
        assert job1.name == job2.name
        assert job1.id == job2.id
        assert job1.path == job2.path
        assert job1.settings == job2.settings
        assert job1._filenames == job2._filenames
        assert job2.status == "successful"
        assert job2.status_log == []
        assert job2.get_errormsg() is None

    def test_load_legacy_multijob_succeeds(self):
        # Given multi job run before additional properties were added
        job1 = DummySingleJob()
        mjob1 = MultiJob(children=[MultiJob(children=[job1])])
        mjob1.run()
        mjob1.results.wait()
        mstatus = mjob1.status
        status = job1.status
        delattr(mjob1, "_status")
        delattr(mjob1, "_status_log")
        delattr(job1, "_status")
        delattr(job1, "_status_log")
        mjob1.__dict__["status"] = mstatus
        job1.__dict__["status"] = status
        mjob1.pickle()

        # When call load
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)
        mjob2 = job_manager.load_job(f"{mjob1.path}/{mjob1.name}.dill")
        job2 = mjob2.children[0].children[0]

        # Then job loaded successfully
        assert mjob1.name == mjob2.name
        assert mjob1.path == mjob2.path
        assert mjob1.settings == mjob2.settings
        assert mjob2.status == "successful"
        assert mjob2.status_log == []
        assert mjob2.get_errormsg() is None
        assert job1.name == job2.name
        assert job1.id == job2.id
        assert job1.path == job2.path
        assert job1.settings == job2.settings
        assert job1._filenames == job2._filenames
        assert job2.status == "successful"
        assert job2.status_log == []
        assert job2.get_errormsg() is None

    @pytest.mark.parametrize(
        "path_exists,folder_exists,use_existing_folder,expected_workdir",
        [
            (True, False, False, "./{}/{}"),
            (True, True, False, "./{}/{}.002"),
            (True, False, True, "./{}/{}"),
            (True, True, True, "./{}/{}"),
            (False, False, False, None),
        ],
        ids=[
            "path_exists_new_folder",
            "path_exists_folder_renamed",
            "path_exists_new_folder_with_use_existing",
            "path_exists_reuse_folder_with_use_existing",
            "path_not_exists_errors",
        ],
    )
    def test_workdir_location(self, path_exists, folder_exists, use_existing_folder, expected_workdir):
        # Given path and folder which may already exist
        path = str(uuid.uuid4())
        folder = str(uuid.uuid4())
        expected_workdir = expected_workdir.format(path, folder) if expected_workdir else None
        if path_exists:
            os.mkdir(path)
            if folder_exists:
                os.mkdir(f"{path}/{folder}")

        if expected_workdir is None:
            # When create jobmanager where path does not exist
            # Then raises error
            with pytest.raises(PlamsError):
                job_manager = JobManager(
                    settings=JobManagerSettings(), path=path, folder=folder, use_existing_folder=use_existing_folder
                )
        else:
            # When create jobmanager where path and folder may exist
            job_manager = JobManager(
                settings=JobManagerSettings(), path=path, folder=folder, use_existing_folder=use_existing_folder
            )

            # Then workdir is created
            assert os.path.abspath(expected_workdir) == job_manager.workdir
            assert os.path.exists(job_manager.workdir)

            job_manager._clean()
            if os.path.exists(job_manager.workdir):
                os.rmdir(job_manager.workdir)

        if os.path.exists(f"{path}/{folder}"):
            os.rmdir(f"{path}/{folder}")
        if os.path.exists(path):
            os.rmdir(path)

    def test_current_rundir_location(self):
        # Given job manager
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)

        # When get the run dir inside and outside context
        rundir1 = job_manager.current_dir_for_jobs
        with jobs_in_directory("foo"):
            with jobs_in_directory("bar"):
                rundir2 = job_manager.current_dir_for_jobs
        rundir3 = job_manager.current_dir_for_jobs

        # Then run dir is the workdir outside the context, and the run dir inside the context
        assert rundir1 == Path(job_manager.workdir)
        assert rundir2 == Path(job_manager.workdir) / "foo" / "bar"
        assert rundir3 == Path(job_manager.workdir)

        os.rmdir(job_manager.workdir)

    def test_register(self):
        # Given job manager
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)

        # When register jobs
        base_name = "test_jobreg"
        job1 = DummySingleJob(name=base_name)
        job2 = DummySingleJob(name=base_name)
        job3 = DummySingleJob(name=base_name)
        job4 = DummySingleJob(name=base_name)
        jobs = [job1, job2, job3, job4]
        job_manager._register(job1)
        job_manager._register(job2)
        with jobs_in_directory("foo"):
            job_manager._register(job3)
            with jobs_in_directory("bar"):
                job_manager._register(job4)

        # Then jobs registered as expected
        def verify_job_registration(job, expected_name, expected_subdir=None):
            # Verify job manager set on job
            assert job.jobmanager == job_manager
            # Verify job name has postfix if duplicate run
            assert job.name == expected_name
            # Verify job path in workdir/subdir/job
            if expected_subdir:
                assert Path(job.path) == Path(job_manager.workdir, expected_subdir, expected_name)
            else:
                assert Path(job.path) == Path(job_manager.workdir, expected_name)
            # Verify job status is registered
            assert job.status == "registered"

        verify_job_registration(job1, base_name)
        verify_job_registration(job2, f"{base_name}.002")
        verify_job_registration(job3, base_name, expected_subdir="foo")
        verify_job_registration(job4, base_name, expected_subdir="foo/bar")
        assert job_manager.jobs == jobs
        assert job_manager.names == {
            "foo/bar/test_jobreg": 1,
            "foo/test_jobreg": 1,
            "test_jobreg": 2,
        }

        job_manager._clean()
        if os.path.exists(job_manager.workdir):
            os.rmdir(job_manager.workdir)

    def test_check_hash(self):
        # Given job manager
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)

        # When check_hash of jobs
        base_name = "test_jobreg"
        job1 = DummySingleJob(name=base_name)
        job2 = DummySingleJob(name=base_name, inp="foo")
        job3 = DummySingleJob(name=base_name, inp="foo")
        job4 = DummySingleJob(name=base_name, inp="foo2")
        jobs = [job1, job2, job3, job4]
        for job in jobs:
            job_manager._check_hash(job)

        assert len(job_manager.hashes) == 3

    def test_load_and_remove(self):
        # Given job manager and saved jobs
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)

        base_name = "test_jobreg"
        job1 = DummySingleJob(name=base_name)
        job2 = DummySingleJob(name=base_name)
        job3 = DummySingleJob(name=base_name)
        job4 = DummySingleJob(name=base_name)
        jobs = [job1, job2, job3, job4]
        job_manager._register(job1)
        job_manager._register(job2)
        with jobs_in_directory("foo"):
            job_manager._register(job3)
            with jobs_in_directory("bar"):
                job_manager._register(job4)
        for job in jobs:
            job.pickle()

        # When load jobs
        folder2 = str(uuid.uuid4())
        job_manager2 = JobManager(settings=JobManagerSettings(), folder=folder2)

        for job in jobs:
            loaded_job = job_manager2.load_job(f"{job.path}{os.path.sep}{job.name}.dill")

            # Then loaded job is the same as the saved job
            assert loaded_job.path == job.path
            assert loaded_job.name == job.name
            assert loaded_job.status == job.status
            assert loaded_job.hash() == job.hash()

        # When remove jobs
        for job in jobs:
            job_manager.remove_job(job)

        # Then jobs removed from job manager
        assert job_manager.jobs == []

        job_manager._clean()
        if os.path.exists(job_manager.workdir):
            os.rmdir(job_manager.workdir)

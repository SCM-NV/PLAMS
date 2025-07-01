import pytest
import os
import uuid

from scm.plams.core.basejob import MultiJob
from scm.plams.core.jobmanager import JobManager
from scm.plams.core.settings import JobManagerSettings
from scm.plams.core.errors import PlamsError
from scm.plams.unit_tests.test_basejob import DummySingleJob


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

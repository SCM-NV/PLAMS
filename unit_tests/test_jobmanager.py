import pytest
import os
import uuid

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

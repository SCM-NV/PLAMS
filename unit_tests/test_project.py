import os
import shutil
from pathlib import Path
import pytest
import uuid
import json
from typing import Dict

from scm.plams.core.errors import ProjectError
from scm.plams.core.project import Project

from .test_basejob import DummySingleJob


class TestProject:

    @pytest.fixture(scope="class", autouse=True)
    def tmp_test_project_dir(self):
        """
        Create new directory for test projects and clean it up after the tests are finished
        """
        test_project_dir = Path(__file__).parent.absolute() / "test_project"
        try:
            shutil.rmtree(test_project_dir)
        except FileNotFoundError:
            pass
        os.mkdir(test_project_dir)
        yield test_project_dir
        shutil.rmtree(test_project_dir)

    def get_random_name(self):
        """
        Get random unique name for project
        """
        return f"test_proj_{uuid.uuid4()}"

    def verify_metadata(self, project: Project, expected: Dict):
        """
        Verify project metadata as expected
        """
        with open(project._metadata_file) as f:
            metadata = json.load(f)
            assert metadata == expected

    def test_create_project(self, tmp_test_project_dir):
        # Given name
        name = self.get_random_name()

        # When create project
        project1 = Project.create(name, description="a test project", parent_dir=tmp_test_project_dir)

        # Then metadata file created in correct location
        assert project1._metadata_file == tmp_test_project_dir / name / ".project.json"
        self.verify_metadata(project1, {"name": name, "description": "a test project"})

        # When create project with same name without explicit parent dir
        project2 = Project.create(name, description="a test project")

        # Then created in current working directory
        assert project2._metadata_file == Path(os.getcwd()) / name / ".project.json"
        shutil.rmtree(project2.path)

        # When create project with same name in same directory
        # Then fails
        with pytest.raises(ProjectError, match=r"Project .* already exists"):
            Project.create(name, parent_dir=tmp_test_project_dir)

    def test_foo(self, tmp_test_project_dir):
        # Given name
        name = self.get_random_name()

        # When create project
        project1 = Project.create(name, description="a test project", parent_dir=tmp_test_project_dir)

        job = DummySingleJob()
        project1.run_job(job)

        assert 1 == 1

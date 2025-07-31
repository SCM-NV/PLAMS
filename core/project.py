import os
import psutil
from typing import Optional, Union, Dict, List
from pathlib import Path
import shutil
import json
from contextlib import contextmanager

from scm.plams.core.errors import ProjectError
from scm.plams.core.functions import log, get_config, config_context
from scm.plams.core.jobmanager import JobManager
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.basejob import SingleJob, Results
from scm.plams.core.enums import ProjectMode, JobMode


class Project:
    _available_physical_cores = psutil.cpu_count(logical=False)

    def __init__(
        self,
        name: str,
        description: Optional[str] = None,
        parent_dir: Optional[Union[str, os.PathLike]] = None,
        mode: ProjectMode = ProjectMode.CREATE,
    ):
        # Validate and set metadata
        if not name or not isinstance(name, str):
            raise ValueError("String value for 'name' must be specified, it cannot be an empty string or 'None'")
        if description and not isinstance(description, str):
            raise ValueError("Value for 'description' must be a string")
        self._name = name
        self._description = description
        self._path = Path(parent_dir if parent_dir else os.getcwd()) / name
        self._mode = mode
        self._jobs: List[SingleJob] = []

        # Initialise project directory
        if self._mode == ProjectMode.CREATE:
            if self._path.exists():
                raise ProjectError(
                    f"Project '{name}' already exists in the directory '{parent_dir}'. Modify the location to create a new project, or change the access mode to use an existing project."
                )
            self.path.mkdir(parents=True)
            log(f"Created project '{name}' in the directory '{parent_dir}'")
        # elif self._mode == "r":
        #     raise NotImplementedError("readonly mode not implemented")
        # elif self._mode == "r+":
        #     raise NotImplementedError("readwrite mode not implemented")
        elif mode == ProjectMode.OVERWRITE:
            if self.path.exists():
                shutil.rmtree(self.path)
                log(f"Deleted existing project '{name}' in the directory '{parent_dir}'")
            self.path.mkdir(parents=True)
            log(f"Created project '{name}' in the directory '{parent_dir}'")
        else:
            raise ValueError(
                f"Invalid access mode '{mode}', must be one of: 'x' (create), 'r' (read), 'r+' (readwrite) or 'w' (write)"
            )

        # Initialise job manager and job runners
        if self._mode != "r":
            self._save()
            self._job_manager = JobManager(
                settings=get_config().jobmanager, path=parent_dir, folder=name, use_existing_folder=True
            )
            self._serial_job_runner = JobRunner()

    @classmethod
    def create(
        cls, name: str, description: Optional[str] = None, parent_dir: Optional[Union[str, os.PathLike]] = None
    ) -> "Project":
        """
        Create a new project.
        """
        return cls(name, description=description, parent_dir=parent_dir)

    @classmethod
    def open(cls, name: str, parent_dir: Optional[Union[str, os.PathLike]] = None):
        raise NotImplementedError("open not implemented")

    # @classmethod
    # def delete(cls, project: "Project"):

    @property
    def name(self) -> str:
        """
        Name of the project. This is also the name of the project directory.

        :return: name of the project
        """
        return self._name

    @property
    def description(self) -> Optional[str]:
        """
        Description of the project and its contents.

        :return: description of the project
        """
        return self._description

    @description.setter
    def description(self, value: Optional[str]):
        self._description = value
        self._save()

    @property
    def path(self) -> Path:
        """
        Absolute path for the project directory.

        :return: full path of the project
        """
        return self._path.resolve()

    @property
    def _metadata_file(self) -> Path:
        return self.path / ".project.json"

    @property
    def _metadata(self) -> Dict:
        return {
            "name": self.name,
            "description": self.description,
            "jobs": [{
                "name": j.name,
                "path": j.path,
                "type": f"{j.__class__.__module__}.{j.__class__.__qualname__}"
            } for j in self._jobs]
        }

    def _save(self):
        metadata = self._metadata
        with open(self._metadata_file, "w") as f:
            json.dump(metadata, f, indent=4)

    @contextmanager
    def job_session(
        self, max_parallel_jobs: int = 1, max_parallel_threads: int = 256, cores_per_job: Optional[int] = None
    ):
        """
        Start a new session context for running jobs. All jobs run in this session context...

        :param max_parallel_jobs:
        :param max_parallel_threads:
        :param cores_per_job:
        :return:
        """
        # Validate session
        if max_parallel_jobs <= 0:
            raise ValueError(f"Value of 'max_parallel_jobs' must be greater than zero, but was {max_parallel_jobs}'")
        if max_parallel_threads <= 0:
            raise ValueError(
                f"Value of 'max_parallel_threads' must be greater than zero, but was {max_parallel_threads}'"
            )
        if cores_per_job is not None and cores_per_job <= 0:
            raise ValueError(f"Value of 'cores_per_job' must be greater than zero or None, but was {cores_per_job}'")

        # Set up parallelization settings
        parallel = max_parallel_jobs != 1 and max_parallel_threads != 1
        job_runner = JobRunner(parallel=parallel, maxjobs=max_parallel_jobs, maxthreads=max_parallel_threads)
        cores_per_job = (
            cores_per_job if cores_per_job is not None else (self._available_physical_cores if not parallel else 1)
        )

        # Enter job context
        log(
            f"Starting job session in project '{self.name}' with {max_parallel_jobs} in parallel and {cores_per_job} cores per job"
        )
        with config_context() as cfg:
            cfg.default_jobmanager = self._job_manager
            cfg.default_jobrunner = job_runner
            cfg.job.runscript.nproc = cores_per_job
            yield

    def run_job(self, job: SingleJob, mode: JobMode = JobMode.CREATE) -> Results:
        # get the job runner and manager from the current session
        # or use the defaults for this project if not in a session
        cfg = get_config()
        job_manager = cfg.default_jobmanager
        job_runner = cfg.default_jobrunner
        in_session = job_manager == self._job_manager
        if not in_session:
            job_manager = self._job_manager
            job_runner = self._serial_job_runner

        try:
            return job.run(job_manager=job_manager, jobrunner=job_runner)
        finally:
            self._jobs.append(job)
            self._save()




# %%
import copy
import re
import warnings
from collections import UserList
from pathlib import Path
from typing import Any, Dict, Hashable, Iterable, List, Optional, Tuple, Type, Union

from scm.plams import JobManager, SingleJob


class GroupedColNames:
    paths = "paths"
    input_settings = "input_settings"
    jobs = "jobs"
    names = "names"
    #########################
    molecule_info = ("n_atoms", "chemical_formula", "gyration_radius", "smiles")
    #########################
    jobs_info = ("ok", "check", "error")
    traj_status = "traj_status"
    amsjob_info = (
        "termination_status",
        "out_errors",
        "log_errors",
        "out_warnings",
        "log_warnings",
    )
    timings = ("CPUTime", "SysTime", "ElapsedTime")
    #########################
    labels = "labels"
    ########################


class JobsAnalysis:
    # _plot_type = ViewJobAnalysis
    _cols = GroupedColNames

    def __init__(
        self,
        paths: Optional[Union[List[str], List[Path]]] = None,
        jobs: Optional[List[SingleJob]] = None,
        data: Optional[Dict[str, List[Any]]] = None,
        extra_cols: Optional[Dict[str, List[Any]]] = None,
    ):
        """_summary_

        :param paths: it must be the path to the folder of each jobs, defaults to None
        :type paths: List[Union[str, Path]], optional
        :param jobs: _description_, defaults to None
        :type jobs: List[plams.SingleJob], optional
        :param input_settings: _description_, defaults to None
        :type input_settings: List[plams.Settings], optional
        :param data: _description_, defaults to None
        :type data: pd.DataFrame, optional
        """
        ################### initialization data ####################
        if data is not None:
            self.data = data
            if GroupedColNames.paths in self.data.keys():
                self.data[GroupedColNames.paths] = [Path(x) for x in self.data[GroupedColNames.paths]]
        else:
            self.data = {}

        if paths is not None:
            self.data[GroupedColNames.paths] = [Path(x) for x in paths]

        if jobs is not None:
            self.data[GroupedColNames.jobs] = jobs
            if self.names is not None:
                self.data[GroupedColNames.names] = self.names
            if self.paths is not None:
                self.data[GroupedColNames.paths] = self.paths

        if extra_cols is not None:
            self.data.update(extra_cols)
        ###################### Cache data ######################
        self.setting_paths_cols: List[str] = []
        self._cache_free_blocks = []
        self._params_cols = None
        self._sal_cols = None
        self._ir_check = None

    @property
    def paths(self):
        if GroupedColNames.jobs in self.data:
            return [Path(job.path) if job.status != "created" else None for job in self.data[GroupedColNames.jobs]]
        if GroupedColNames.paths in self.data:
            return self.data[GroupedColNames.paths]
        return None

    @property
    def names(self):
        if GroupedColNames.jobs in self.data:
            return [job.path for job in self.data[GroupedColNames.jobs]]
        if GroupedColNames.paths in self.data:
            return self.data[GroupedColNames.paths]
        return None

    @classmethod
    def load_paths_of_inputs(
        cls,
        base_path: str,
        pattern: str = "*/*.in",
    ):
        paths = [i.parent for i in Path(base_path).glob(pattern)]
        assert len(paths) > 0, "Any path found"
        return cls(paths=paths)

    @classmethod
    def load_only_job_failed(
        cls,
        logfile_path: Path,
    ):
        assert logfile_path.exists(), f"{logfile_path} does not exists"

        with open(logfile_path, "r") as log_file:
            log_data = log_file.read()
        regex_pattern = r"JOB (\S+) FAILED"
        matches = list(re.findall(regex_pattern, log_data))
        if len(matches) > 0:
            return cls(paths=matches)
        else:
            print("All ok")
            return None

    def load_jobs(
        self,
        job_loader,
        col_path_jobs=GroupedColNames.paths,
        suffix_out="",
    ):

        def load_job(path_folder, job_loader):
            """plams.load .dill file is very fast, always try that first"""
            if not isinstance(path_folder, (str, Path)):
                return None
            if not Path(path_folder).exists():
                return None
            path_files_in_job: Path = path_folder / (path_folder.name)
            is_dill = path_files_in_job.with_suffix(".dill")
            jm = JobManager({}, folder=Path.cwd(), use_existing_folder=True)
            if is_dill.exists():
                job = jm.load_job(is_dill)
            else:
                job = None
            if job is None:
                job = job_loader(path_folder)
            return job

        self.data[GroupedColNames.jobs + suffix_out] = [load_job(p, job_loader) for p in self.data[col_path_jobs]]
        self.data[GroupedColNames.names + suffix_out] = self.data[GroupedColNames.jobs].apply(lambda x: x.name)
        self.data[GroupedColNames.names + suffix_out] = self.data[GroupedColNames.names]

    def _check_jobs_types(self, job_type=SingleJob, jobs_col=GroupedColNames.jobs, raise_error=False):
        if len(self.data[jobs_col]) > 0:
            actual_job_type = self.get_col_type(col=jobs_col)
            type_check = isinstance(self.data[jobs_col][0], job_type)
        else:
            raise ValueError(f"{len(self.data[jobs_col])=} not higher than zero, something went wrong")

        if raise_error and not type_check:
            raise TypeError(f"The first job in {jobs_col} must be of type {job_type}, but found {actual_job_type}")
        else:
            return type_check

    def get_col_type(self, col: str = GroupedColNames.jobs):
        actual_job_type = type(self.data[col][0])
        return actual_job_type

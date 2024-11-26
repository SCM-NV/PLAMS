# %%
import re
import warnings
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Type, Union

from scm.plams import AMSJob, JobManager, Settings, SingleJob, to_smiles
from scm.plams.core.settings_analysis import compare_settings, print_in_table


class GroupedColNames:
    paths = "paths"
    jobs = "jobs"
    names = "names"
    #########################
    molecule_info = ("n_atoms", "chemical_formula", "gyration_radius", "smiles")
    #########################
    jobs_info = ("ok", "check", "error")
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
        paths: Optional[Union[List[Optional[str]], List[Optional[Path]], List[Path]]] = None,
        jobs: Optional[List[SingleJob]] = None,
        data: Optional[Dict[str, List[Any]]] = None,
        extra_cols: Optional[Dict[str, List[Any]]] = None,
    ):
        ################### initialization data ####################
        self.data = {}
        if data is not None:
            self.data = data
            if GroupedColNames.paths in self.data.keys():
                self.data[GroupedColNames.paths] = [Path(x) for x in self.data[GroupedColNames.paths]]

        if paths is not None:
            self.data[GroupedColNames.paths] = [Path(x) if x is not None else None for x in paths]

        if jobs is not None:
            self.data[GroupedColNames.jobs] = jobs
            self._generate_names_paths_from_jobs()

        if extra_cols is not None:
            self.data.update(extra_cols)
        ###################### Cache data ######################
        self.setting_paths_cols: List[str] = []
        self._cache_free_blocks = []
        self._params_cols = None
        self._sal_cols = None
        self._ir_check = None

    def _generate_names_paths_from_jobs(
        self, col_paths=GroupedColNames.paths, col_names=GroupedColNames.paths, col_jobs=GroupedColNames.jobs
    ):
        self.names(col_paths, col_jobs)
        self.paths(col_names, col_jobs)

    def paths(self, col_paths=GroupedColNames.paths, col_jobs=GroupedColNames.jobs):
        if col_paths in self.data:
            return self.data[col_paths]
        if col_jobs in self.data:
            self.data[col_paths] = [Path(job.path) if job.status != "created" else None for job in self.data[col_jobs]]
            return self.data[col_paths]
        return None

    def names(self, col_names=GroupedColNames.paths, col_jobs=GroupedColNames.jobs):
        if col_names in self.data:
            return self.data[col_names]
        if col_jobs in self.data:
            self.data[col_names] = [job.name for job in self.data[col_jobs]]
            return self.data[col_names]
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
            path_folder = Path(path_folder)
            if not Path(path_folder).exists():
                return None
            path_files_in_job = path_folder / (path_folder.name)
            is_dill = path_files_in_job.with_suffix(".dill")
            jm = JobManager({}, folder=Path.cwd(), use_existing_folder=True)
            job = None
            if is_dill.exists():
                job = jm.load_job(is_dill)
            if job is None:
                job = job_loader(path_folder)
            return job

        self.data[GroupedColNames.jobs + suffix_out] = [load_job(p, job_loader) for p in self.data[col_path_jobs]]
        self._generate_names_paths_from_jobs(
            col_jobs=GroupedColNames.jobs + suffix_out,
            col_names=GroupedColNames.names + suffix_out,
            col_paths=GroupedColNames.paths + suffix_out,
        )

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

    ############################################################################
    ####################          input analysis          ######################
    ############################################################################
    def get_nested_settings(
        self,
        settings_path: str,
        default_val: Any = None,
        default_type: Optional[Type] = None,
        job_col=GroupedColNames.jobs,
        col_suffix="",
    ):
        def get_value(job: SingleJob):
            value = job.settings.get_nested(settings_path, default=default_val)
            if default_type is not None and value is not None:
                return default_type(value)
            return value

        added_col = settings_path + col_suffix
        self.data[added_col] = [get_value(job) for job in self.data[job_col]]
        self.setting_paths_cols.append(added_col)
        self.setting_paths_cols.sort()
        return added_col

    def compare_settings(
        self,
        job_col=GroupedColNames.jobs,
        blocks_analysis: bool = True,
        keys_analysis: bool = False,
        default_settings: Optional[Settings] = None,
        flatten_list: bool = True,
        clean_singular_value: bool = True,
        col_suffix="",
    ):
        settings_list = [job.settings for job in self.data[job_col]]
        settings_summary = compare_settings(
            settings_list,
            blocks_analysis=blocks_analysis,
            keys_analysis=keys_analysis,
            default_settings=default_settings,
            flatten_list=flatten_list,
            clean_singular_value=clean_singular_value,
        )
        added_cols = []
        for k, v in settings_summary.items():
            k = ".".join(map(str, k)) + col_suffix
            self.data[k] = v
            added_cols.append(k)
            self.setting_paths_cols.append(k)
        self.setting_paths_cols.sort()
        return added_cols

    def generate_labels(
        self, select_cols: Union[List[str], Dict[str, str]], cols_separator="\n", col_label=GroupedColNames.labels
    ):
        """you have to give the cols of the data that you want"""
        col_names = list(self.data)
        if isinstance(select_cols, dict):
            cols = [col_names.index(v) for v in select_cols.values()]
            new_names = list(select_cols)
        else:
            cols = [col_names.index(v) for v in select_cols]
            new_names = select_cols

        def generate_label(row):
            support_list = []
            for col_i, name_i in zip(cols, new_names):
                support_list.append(f"{name_i}: {row[col_i]}")
            return f"{cols_separator}".join(support_list)

        self.data[col_label] = [generate_label(row) for row in zip(*self.data.values())]
        return col_label

    ############################################################################
    ##################          molecule analysis          #####################
    ############################################################################

    def get_molecules_infos(self, gyration_radius=False, smiles=False, jobs_col=GroupedColNames.jobs):

        if not self._check_jobs_types(job_type=SingleJob, raise_error=False, jobs_col=jobs_col):
            warnings.warn("the jobs are not of type: plams.AMSJob so no get_molecules_infos are present")
            return

        self.data["n_atoms"] = [len(x.molecule) for x in self.data[jobs_col]]
        self.data["chemical_formula"] = [x.molecule.get_formula() for x in self.data[jobs_col]]

        def _to_smiles(job):
            mol = job.molecule
            mol.delete_all_bonds()
            mol.guess_bonds()
            return to_smiles(mol, canonical=True)

        if smiles:
            self.data["smiles"] = [to_smiles(x) for x in self.data[jobs_col]]

        if gyration_radius:
            self.data["gyration_radius"] = [x.molecule.get_gyration_radius() for x in self.data[jobs_col]]

    ############################################################################
    ####################    error analysis and timings    ######################
    ############################################################################

    def get_job_info(
        self,
        jobs_col=GroupedColNames.jobs,
    ):
        def check_not_created(job) -> bool:
            if job.status not in ["created"]:
                return True
            return False

        def check_ams_status_not_none(job):
            try:
                status = job.results.readrkf("General", "termination status")
            except:
                return False
            if status is None:
                return False
            return True

        self.data["ok"] = [job.ok() if check_not_created(job) else None for job in self.data[jobs_col]]
        self.data["check"] = [
            job.check() if check_not_created(job) and check_ams_status_not_none(job) else None
            for job in self.data[jobs_col]
        ]
        self.data["error"] = [
            job.get_errormsg() if check_not_created(job) and check_ams_status_not_none(job) else None
            for job in self.data[jobs_col]
        ]

    def get_timings(self, col_jobs=GroupedColNames.jobs, custom_get_timings: Optional[Callable] = None):
        if self._check_jobs_types(job_type=AMSJob, raise_error=False):
            cols_needed = GroupedColNames.timings
            cols_needed_present = all([i in self.data for i in cols_needed])
            if not cols_needed_present:
                self.data["CPUTime"] = [job.results.readrkf("General", "CPUTime") for job in self.data[col_jobs]]
                self.data["SysTime"] = [job.results.readrkf("General", "SysTime") for job in self.data[col_jobs]]
                self.data["ElapsedTime"] = [
                    job.results.readrkf("General", "ElapsedTime") for job in self.data[col_jobs]
                ]
        # elif self._check_jobs_types(job_type=params.ParAMSJob, raise_error=False):

        #     def get_time(job: params.ParAMSJob):
        #         val = job.results.get_timings()
        #         if val is not None:
        #             return val.total_seconds() / 3600
        #         else:
        #             return None

        #     self.data["timings[h]"] = [get_time(job) for job in self.data[GroupedColNames.jobs]]

        elif custom_get_timings is not None:
            self.data["timings"] = [custom_get_timings(job) for job in self.data[col_jobs]]

    ############################################################################
    #####################          JobSummarize          #######################
    ############################################################################
    def get_job_summary(self, summarizer: Dict[str, Callable[[SingleJob], Any]], jobs_col=GroupedColNames.jobs):
        for k, v_call in summarizer.items():
            self.data[k] = [v_call(job) for job in self.data[jobs_col]]

    ############################################################################
    ##################            tests analysis           #####################
    ############################################################################
    def reasonable_checker(
        self, col_out: str, reasonable_checker: Callable[[SingleJob], bool], jobs_col=GroupedColNames.jobs
    ):
        self.data[col_out] = [reasonable_checker(job) for job in self.data[jobs_col]]

    def success_checker(
        self,
        ref_job_idx: int,
        indexes_to_check: List[int],
        success_checker: Callable[[SingleJob, SingleJob], Dict[str, Union[bool, float]]],
        jobs_col=GroupedColNames.jobs,
    ):
        for idx_i in indexes_to_check:
            res = success_checker(self.data[jobs_col][idx_i], self.data[jobs_col][ref_job_idx])
            for k, v in res:
                if k not in self.data:
                    self.data[k] = [None] * len(self.data[jobs_col])
                self.data[k][idx_i] = v

    def print_table(self):
        print_in_table(self.data)

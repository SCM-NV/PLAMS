# %%
import csv
import re
import warnings
from collections import UserList, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Type, Union

from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from scm.plams import (
    AMSJob,
    JobManager,
    JobManagerSettings,
    Settings,
    SingleJob,
    to_smiles,
)
from scm.plams.tools.settings_analysis import compare_settings
from scm.plams.tools.table_formatter import format_in_table


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


class HeaderStr(str):
    @property
    def is_settings(self):
        return self._is_settings

    @is_settings.setter
    def is_settings(self, val: bool):
        self._is_settings = val

    def __new__(cls, value, is_settings: bool = False):
        # Create a new instance of str
        obj = super(HeaderStr, cls).__new__(cls, value)
        # Add the extra property
        obj.is_settings = is_settings
        return obj


class SettingsCols:

    def __init__(self, __dict__: Dict[Union[str, HeaderStr], List[Any]]):
        self.__dict__ = __dict__

    @property
    def values(self):
        return [k for k in self.__dict__ if isinstance(k, HeaderStr) and k.is_settings]

    def __repr__(self) -> str:
        list_repr = str(sorted(self.values)).replace(",", ",\n\t")
        return f"SettingsCols({list_repr})"

    def rename(self, mapper: Dict[str, str]):
        for k_in, k_out in mapper.items():
            self.__dict__[HeaderStr(k_out, is_settings=True)] = self.__dict__.pop(k_in)


class JobsAnalysis:
    # _plot_type = ViewJobAnalysis
    _cols = GroupedColNames

    def __init__(
        self,
        paths: Optional[Union[List[Optional[str]], List[Optional[Path]], List[Path]]] = None,
        jobs: Optional[List[SingleJob]] = None,
        data: Optional[Dict[Union[str, HeaderStr], List[Any]]] = None,
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
        # self._params_cols = None
        # self._sal_cols = None
        # self._ir_check = None

        self.settings_cols = SettingsCols(__dict__=self.data)

    def apply(self, fn: Callable, col: str = GroupedColNames.jobs):
        assert col in self.data, f"{col} not in {self.data.keys()=}"
        return [fn(x) for x in self.data[col]]

    def _generate_names_paths_from_jobs(
        self, col_paths=GroupedColNames.paths, col_names=GroupedColNames.names, col_jobs=GroupedColNames.jobs
    ):
        self.names(col_names, col_jobs)
        self.paths(col_paths, col_jobs)

    def paths(self, col_paths=GroupedColNames.paths, col_jobs=GroupedColNames.jobs):
        if col_paths in self.data:
            return self.data[col_paths]
        if col_jobs in self.data:
            self.data[col_paths] = self.apply(
                fn=lambda job: Path(job.path) if job.status != "created" else None, col=col_jobs
            )
            return self.data[col_paths]
        return None

    def names(self, col_names=GroupedColNames.names, col_jobs=GroupedColNames.jobs):
        if col_names in self.data:
            return self.data[col_names]
        if col_jobs in self.data:
            self.data[col_names] = self.apply(fn=lambda job: job.name, col=col_jobs)
            return self.data[col_names]
        return None

    @classmethod
    def load_paths_of_inputs(
        cls,
        base_path: Union[str, Path],
        pattern: str = "*/*.in",
    ):
        paths = [i.parent for i in Path(base_path).glob(pattern)]
        assert len(paths) > 0, "Any path found"
        return cls(paths=paths)

    @classmethod
    def load_only_job_failed(
        cls,
        logfile_path: Union[str, Path],
    ):
        logfile_path = Path(logfile_path)
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

        def load_job(path_folder):
            """plams.load .dill file is very fast, always try that first"""
            if not isinstance(path_folder, (str, Path)):
                return None
            path_folder = Path(path_folder)
            if not Path(path_folder).exists():
                return None
            path_files_in_job = path_folder / (path_folder.name)
            is_dill = path_files_in_job.with_suffix(".dill")
            jm = JobManager(JobManagerSettings(), folder=Path.cwd(), use_existing_folder=True)
            job = None
            if is_dill.exists():
                job = jm.load_job(is_dill)
            if job is None:
                job = job_loader(path_folder)
            return job

        self.data[GroupedColNames.jobs + suffix_out] = self.apply(fn=load_job, col=col_path_jobs)
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
        self.data[HeaderStr(added_col, is_settings=True)] = self.apply(fn=get_value, col=job_col)
        return added_col

    def compare_settings(
        self,
        job_col=GroupedColNames.jobs,
        analyze_blocks: bool = True,
        analyze_keys: bool = False,
        default_settings: Optional[Settings] = None,
        flatten_list: bool = True,
        remove_unimportant_columns: bool = True,
        unimportant_variation_threshold: int = 1,
        none_is_unimportant: bool = False,
        col_suffix="",
    ):
        for k in self.settings_cols.values:
            self.data.pop(k)
        settings_list = [job.settings for job in self.data[job_col]]
        settings_summary = compare_settings(
            settings_list,
            analyze_blocks=analyze_blocks,
            analyze_keys=analyze_keys,
            default_settings=default_settings,
            flatten_list=flatten_list,
            remove_unimportant_columns=remove_unimportant_columns,
            unimportant_variation_threshold=unimportant_variation_threshold,
            none_is_unimportant=none_is_unimportant,
        )
        added_cols = []
        for k, v in settings_summary.items():
            k = ".".join(map(str, k)) + col_suffix
            self.data[HeaderStr(k, is_settings=True)] = v
            added_cols.append(k)
        return added_cols

    def generate_labels(
        self, select_cols: Union[List[str], Dict[str, str]], cols_separator="\n", col_label=GroupedColNames.labels
    ):
        """you have to give the cols of the data that you want"""
        # col_names = list(self.data)
        if isinstance(select_cols, dict):
            # cols = [col_names.index(v) for v in select_cols.values()]
            cols = list(select_cols.values())
            new_names = list(select_cols)
        else:
            # cols = [col_names.index(v) for v in select_cols]
            cols = select_cols
            new_names = select_cols

        def generate_label(row):
            support_list = []
            for col_i, name_i in zip(cols, new_names):
                support_list.append(f"{name_i}: {row[col_i]}")
            return f"{cols_separator}".join(support_list)

        self.data[col_label] = [generate_label(row) for row in self]
        return col_label

    def __getitem__(self, idx: int):
        row = {}
        for k in self.data:
            row[k] = self.data[k][idx]
        return row

    ############################################################################
    ##################          molecule analysis          #####################
    ############################################################################

    def get_molecules_infos(self, gyration_radius=False, smiles=False, col_jobs=GroupedColNames.jobs):

        if not self._check_jobs_types(job_type=SingleJob, raise_error=False, jobs_col=col_jobs):
            warnings.warn("the jobs are not of type: plams.AMSJob so no get_molecules_infos are present")
            return

        self.data["n_atoms"] = self.apply(fn=lambda x: len(x.molecule), col=col_jobs)
        self.data["chemical_formula"] = self.apply(fn=lambda x: x.molecule.get_formula(), col=col_jobs)

        def _to_smiles(job):
            mol = job.molecule
            mol.delete_all_bonds()
            mol.guess_bonds()
            return to_smiles(mol, canonical=True)

        if smiles:
            self.data["smiles"] = self.apply(fn=to_smiles, col=col_jobs)

        if gyration_radius:
            self.data["gyration_radius"] = self.apply(fn=lambda x: x.molecule.get_gyration_radius(), col=col_jobs)

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

        self.data["ok"] = self.apply(fn=lambda job: job.ok() if check_not_created(job) else None, col=jobs_col)
        self.data["check"] = self.apply(
            fn=lambda job: job.check() if check_not_created(job) and check_ams_status_not_none(job) else None,
            col=jobs_col,
        )
        self.data["error"] = self.apply(
            fn=lambda job: job.get_errormsg() if check_not_created(job) and check_ams_status_not_none(job) else None,
            col=jobs_col,
        )
        cols_added = ["ok", "check", "error"]
        return cols_added

    def get_timings(self, col_jobs=GroupedColNames.jobs, custom_get_timings: Optional[Callable] = None):
        cols_added = []
        if self._check_jobs_types(job_type=AMSJob, raise_error=False):
            cols_needed = GroupedColNames.timings
            cols_needed_present = all([i in self.data for i in cols_needed])
            if not cols_needed_present:
                self.data["CPUTime"] = self.apply(
                    fn=lambda job: job.results.readrkf("General", "CPUTime"), col=col_jobs
                )
                self.data["SysTime"] = self.apply(
                    fn=lambda job: job.results.readrkf("General", "SysTime"), col=col_jobs
                )
                self.data["ElapsedTime"] = self.apply(
                    fn=lambda job: job.results.readrkf("General", "ElapsedTime"), col=col_jobs
                )
                cols_added.extend(["CPUTime", "SysTime", "ElapsedTime"])
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
            cols_added.append("timings")
        return cols_added

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

    def success_checker_plotter(
        self,
        ref_job_idx: int,
        indexes_to_check: List[int],
        success_plot: Callable[[SingleJob, plt.Axes], plt.Axes],
        jobs_col=GroupedColNames.jobs,
        separate_plot=True,
        **plt_kwargs,
    ) -> Figure:
        """Note: success_plot should first plot the non-ref job and then the reference job, such that the legend is displaced well"""
        plt_kwargs.setdefault("layout", "tight")
        if separate_plot:
            plt_kwargs.setdefault("sharex", True)
            plt_kwargs.setdefault("sharey", True)
            plt_kwargs.setdefault(
                "ncols", 4 if "nrows" not in plt_kwargs else -(-len(indexes_to_check) // plt_kwargs["nrows"])
            )
            plt_kwargs.setdefault("nrows", -(-len(indexes_to_check) // plt_kwargs["ncols"]))
            fig, axes = plt.subplots(**plt_kwargs)
            axes = axes.ravel()
        else:
            if "ncols" in plt_kwargs:
                plt_kwargs.pop("ncols")
            if "nrows" in plt_kwargs:
                plt_kwargs.pop("nrows")
            fig, axes = plt.subplots(**plt_kwargs)
            axes = [axes] * len(indexes_to_check)

        text_labels = []
        for i, (idx_i, ax) in enumerate(zip(indexes_to_check, axes)):

            if GroupedColNames.labels in self.data:
                text_label = self.data[GroupedColNames.labels][idx_i]
            else:
                text_label = ""
                for col_i in self.settings_cols.values:
                    text_label += f"{col_i}: {self.data[col_i][idx_i]}\n"

            success_plot(self.data[jobs_col][idx_i], ax)
            if separate_plot:
                success_plot(self.data[jobs_col][ref_job_idx], ax)
                text_labels.append([text_label, "ref"])
            elif i == 0:
                # plot only once
                success_plot(self.data[jobs_col][ref_job_idx], ax)
                text_labels.append(text_label)
                text_labels.append("ref")
            else:
                text_labels.append(text_label)

        if separate_plot:
            for ax, text_label in zip(axes, text_labels):
                ax.set_title(text_label[0])
                ax.legend(text_label)
        else:
            axes[0].legend(text_labels)
        return fig

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

    ############################################################################
    ##################               visualize             #####################
    ############################################################################
    def view_table(
        self,
        indexes=None,
        ret_str: bool = False,
        print_on: bool = True,
        max_col_length: int = -1,
        max_rows_displayed: int = 30,
    ):

        if indexes is not None:
            data = defaultdict(list)
            for i in indexes:
                for k, v in self[i].items():
                    data[k].append(v)
        else:
            data = self.data
        return format_in_table(
            data,
            ret_str=ret_str,
            print_on=print_on,
            max_col_length=max_col_length,
            max_rows_displayed=max_rows_displayed,
        )

    def to_csv(self, data_path: Union[str, Path], pop_jobs_cols=GroupedColNames.jobs):
        headers = [str(key) for key in self.data.keys()]

        jobs_col_idx = None
        if pop_jobs_cols in headers:
            jobs_col_idx = headers.index(pop_jobs_cols)
            headers.pop(jobs_col_idx)

        with open(data_path, mode="w", newline="", encoding="utf-8") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(headers)
            for row in self:
                row_l = list(row.values())

                if jobs_col_idx is not None:
                    row_l.pop(jobs_col_idx)

                writer.writerow(row_l)

        print(f"Data saved to {data_path}")

    @classmethod
    def from_csv(cls, data_path: str):
        reconstructed_data: Dict[Union[str, HeaderStr], List[Any]] = {}

        with open(data_path, mode="r", newline="", encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            # Read headers
            headers = next(reader)
            # Initialize dictionary with headers as keys
            reconstructed_data = {header: [] for header in headers}
            # Read rows and populate the dictionary
            for row in reader:
                for header, value in zip(headers, row):
                    # Convert value back to original type (e.g., int if possible)
                    reconstructed_data[header].append(value if value != "" else None)
        return cls(data=reconstructed_data)

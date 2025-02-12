# %%
import csv
import re
import warnings
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Callable, Dict, Hashable, List, Optional, Tuple, Type, Union

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
from scm.plams.core.functions import get_logger, requires_optional_package
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

    def __iter__(self):
        return iter(self.values)


class JobsAnalysis:
    # _plot_type = ViewJobAnalysis
    _cols = GroupedColNames

    def __init__(
        self,
        paths: Optional[Union[List[Optional[str]], List[Optional[Path]], List[Path]]] = None,
        jobs: Optional[List[SingleJob]] = None,
        data: Optional[Dict[Union[str, HeaderStr], List[Any]]] = None,  # df_jobs.to_dict("list")
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

    def apply(self, fn: Callable, col: Optional[str] = GroupedColNames.jobs, indexes=None):
        if indexes is None:
            indexes = range(len(self))
        if col is not None:
            if col not in self.data:
                raise KeyError(f"{col} not in {self.data.keys()=}")
            return [fn(self.data[col][i]) for i in indexes]

        return [fn(self.data[i]) for i in indexes]

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
        use_dill=True,
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
            jm = JobManager(
                JobManagerSettings(),
                folder=Path.cwd(),
                use_existing_folder=True,
                job_logger=get_logger("none", fmt="csv"),
            )
            job = None
            if is_dill.exists() and use_dill:
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

    def __len__(self):
        first_key = list(self.data)[0]
        return len(self.data[first_key])

    ############################################################################
    ##################          molecule analysis          #####################
    ############################################################################

    def get_molecules_infos(self, gyration_radius=False, smiles=False, col_jobs=GroupedColNames.jobs):

        if not self._check_jobs_types(job_type=SingleJob, raise_error=False, jobs_col=col_jobs):
            warnings.warn("the jobs are not of type: plams.AMSJob so no get_molecules_infos are present")
            return

        def get_natoms(job):
            if job.molecule is None:
                return None
            return len(job.molecule)

        def get_formula(job):
            if job.molecule is None:
                return None
            return job.molecule.get_formula()

        self.data["n_atoms"] = self.apply(get_natoms, col=col_jobs)
        self.data["chemical_formula"] = self.apply(get_formula, col=col_jobs)

        def _to_smiles(job):
            mol = job.molecule
            mol.delete_all_bonds()
            mol.guess_bonds()
            return to_smiles(mol, canonical=True)

        if smiles:
            self.data["smiles"] = self.apply(fn=_to_smiles, col=col_jobs)

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

    def groupby(self, groupby: Callable[[Dict], Hashable], idxs=None):
        groups = defaultdict(list)
        if idxs is None:
            idxs = range(len(self))
        for i in idxs:
            groups[groupby(self[i])].append(i)
        return dict(groups)

    def assign_reference(self, is_ref: Callable[[Dict], bool], subset_idxs: Optional[List[int]] = None):
        if subset_idxs is None:
            subset_idxs = list(range(len(self)))
        groups_with_ref: Dict[str, Optional[Union[int, List[int]]]] = {"ref": None, "idxs": []}
        for idx_i in subset_idxs:
            if is_ref(self[idx_i]):
                if groups_with_ref["ref"] is not None:
                    ref_found = groups_with_ref["ref"]
                    raise ValueError(f"Collision of two refs found for the same group: {ref_found} and {idx_i}")
                groups_with_ref["ref"] = idx_i
            else:
                groups_with_ref["idxs"].append(idx_i)
        return groups_with_ref

    def job_plotter(
        self,
        success_plot: Callable[[SingleJob, str, plt.Axes], str],
        jobs_col=GroupedColNames.jobs,
        grouped_indexes: Union[Dict[str, List[int]], List[int]] = None,
        **plt_kwargs,
    ):
        if grouped_indexes is None or isinstance(grouped_indexes, list):
            if not isinstance(grouped_indexes, list):
                grouped_indexes = range(len(self))
            support = {}
            n_jobs = len(f"{len(self)}")
            for i in grouped_indexes:
                support["Job {:0{}}".format(i, n_jobs)] = [i]
            grouped_indexes = support

        print(grouped_indexes)

        n_axes = len(grouped_indexes)
        plt_kwargs.setdefault("layout", "tight")
        plt_kwargs.setdefault("sharex", True)
        plt_kwargs.setdefault("sharey", False)
        plt_kwargs.setdefault("ncols", 4 if "nrows" not in plt_kwargs else -(-n_axes // plt_kwargs["nrows"]))
        plt_kwargs.setdefault("nrows", -(-n_axes // plt_kwargs["ncols"]))
        fig, axes = plt.subplots(**plt_kwargs)
        axes = axes.ravel()

        for i, (title_i, ax) in enumerate(zip(grouped_indexes, axes)):
            for idx_i in grouped_indexes[title_i]:
                text_label = self.get_label(idx_i)
                success_plot(self.data[jobs_col][idx_i], text_label, ax)
            ax.legend()
            ax.set_title(f"{title_i}")
        return fig

    def success_checker_plotter(
        self,
        ref_job_idx: int,
        indexes_to_check: List[int],
        success_plot: Callable[[SingleJob, str, plt.Axes], plt.Axes],
        jobs_col=GroupedColNames.jobs,
        separate_plot=True,
        axes=None,
        **plt_kwargs,
    ) -> Figure:
        """Note: success_plot should first plot the non-ref job and then the reference job, such that the legend is displaced well"""

        n_axes = len(indexes_to_check)
        plt_kwargs.setdefault("layout", "tight")
        if separate_plot and axes is None:
            plt_kwargs.setdefault("sharex", True)
            plt_kwargs.setdefault("sharey", True)
            plt_kwargs.setdefault("ncols", 4 if "nrows" not in plt_kwargs else -(-n_axes // plt_kwargs["nrows"]))
            plt_kwargs.setdefault("nrows", -(-n_axes // plt_kwargs["ncols"]))
            fig, axes = plt.subplots(**plt_kwargs)
            axes = axes.ravel()
        else:
            if axes is None:
                if "ncols" in plt_kwargs:
                    plt_kwargs.pop("ncols")
                if "nrows" in plt_kwargs:
                    plt_kwargs.pop("nrows")
                fig, axes = plt.subplots(**plt_kwargs)
                axes = [axes] * n_axes
            else:
                axes = axes * n_axes
                fig = axes[0].get_figure()

        for i, (idx_i, ax) in enumerate(zip(indexes_to_check, axes)):
            if idx_i is None:
                continue
            ref_label = "REF\n" + self.get_label(ref_job_idx)
            text_label = self.get_label(idx_i)

            success_plot(self.data[jobs_col][idx_i], text_label, ax)
            if separate_plot:
                success_plot(self.data[jobs_col][ref_job_idx], ref_label, ax)
            elif i == 0:
                # plot only once
                success_plot(self.data[jobs_col][ref_job_idx], ref_label, ax)

        if separate_plot:
            for ax in axes:
                ax.legend()
        else:
            axes[0].legend()
        return fig

    def get_label(self, idx: int):
        if GroupedColNames.labels in self.data:
            text_label = self.data[GroupedColNames.labels][idx]
        else:
            text_label = ""
            for col_i in self.settings_cols.values:
                text_label += f"{col_i}: {self.data[col_i][idx]}\n"
        return text_label

    def success_checker_plotter_2(
        self,
        is_ref: Callable[[Dict], bool],
        groupby: Callable[[Dict], Hashable],
        success_plot: Callable[[SingleJob, plt.Axes], plt.Axes],
        jobs_col=GroupedColNames.jobs,
        **plt_kwargs,
    ) -> Figure:
        """Note: success_plot should first plot the non-ref job and then the reference job, such that the legend is displaced well"""
        groups = self.groupby(groupby)
        groups_refs = {}
        for k, vals in groups.items():
            groups_refs[k] = self.assign_reference(is_ref, subset_idxs=vals)

        fig = self.plot_groups(groups_refs, success_plot=success_plot, jobs_col=jobs_col, **plt_kwargs)
        return fig

    def plot_groups(
        self,
        groups_refs,
        success_plot: Callable[[SingleJob, plt.Axes], plt.Axes],
        jobs_col=GroupedColNames.jobs,
        **plt_kwargs,
    ):
        n_axes = len(groups_refs)
        plt_kwargs.setdefault("layout", "tight")
        plt_kwargs.setdefault("sharex", True)
        plt_kwargs.setdefault("sharey", False)
        plt_kwargs.setdefault("ncols", 4 if "nrows" not in plt_kwargs else -(-n_axes // plt_kwargs["nrows"]))
        plt_kwargs.setdefault("nrows", -(-n_axes // plt_kwargs["ncols"]))
        fig, axes = plt.subplots(**plt_kwargs)
        axes = axes.ravel()
        for (k, vals), ax in zip(groups_refs.items(), axes):
            fig = self.success_checker_plotter(
                ref_job_idx=vals["ref"],
                indexes_to_check=vals["idxs"],
                success_plot=success_plot,
                separate_plot=False,
                axes=[ax],
                jobs_col=jobs_col,
            )
            ax.set_title(f"Group {k}")
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
            for k, v in res.items():
                if k not in self.data:
                    self.data[k] = [None] * len(self.data[jobs_col])
                self.data[k][idx_i] = v

    @requires_optional_package("pigeon")
    def notebook_annotate(
        self,
        job_plotter: Optional[Callable[[SingleJob, plt.Axes], plt.Axes]] = None,
        ref_idx: Optional[int] = None,
        options: Optional[Union[List[str], Tuple[float, float]]] = None,
        shuffle: bool = False,
        include_skip: bool = True,
    ):
        """Wrap around pigeon

        Parameters
        ----------
        options: list(any) or tuple(start, end, [step]) or None
                if list: list of labels for binary classification task (Dropdown or Buttons)
                if tuple: range for regression task (IntSlider or FloatSlider)
                if None: arbitrary text input (TextArea)
        shuffle: bool, shuffle the examples before annotating
        include_skip: bool, include option to skip example while annotating
        job_plotter: func, function for plotting the job
        :param ref_idx: index for plotting two examples and this one is the ref, defaults to None
        :type ref_idx: Optional[int], optional
        """
        import io

        import matplotlib.pyplot as plt
        from IPython.display import Image, display
        from pigeon import annotate

        def plt_plot(idx, row, job_plotter, ref=None):
            fig, ax = plt.subplots()
            job_plotter(row["jobs"], ax)
            if ref:
                job_plotter(ref["jobs"], ax)
            ax.legend(["aaa", "REF"])
            return fig

        def plt_to_image(fig):
            buf = io.BytesIO()
            fig.savefig(buf, format="png")
            buf.seek(0)
            plt.close(fig)  # Close the figure to free resources
            return Image(data=buf.getvalue())

        if ref_idx is not None:
            ref = self[ref_idx]
        if job_plotter is not None:
            display_fn = lambda row: display(plt_to_image(plt_plot(*row, job_plotter, ref=ref)))

        annotations = annotate(
            enumerate(self),
            display_fn=display_fn,
            options=options,
            shuffle=shuffle,
            include_skip=include_skip,
        )
        self.add_annotations(
            annotations,
            annotation_col="annotations",
        )
        return annotations

    def add_annotations(
        self,
        annotations,
        annotation_col="annotations",
    ):
        pigeon_annotation = {x[0][0]: x[1] for x in annotations}
        annotations_all = []
        for i in range(len(self)):
            annotations_all.append(pigeon_annotation.get(i, None))
        self.data[annotation_col] = annotations_all
        return annotation_col

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


class AMSJobErrorChecker:

    @staticmethod
    def find_log_file(job: AMSJob):
        """get the ams.log path or if not exists returns None"""
        if job.path is None:
            return None
        log_path = Path(job.path) / "ams.log"
        if not log_path.exists():
            return None
        return log_path

    @staticmethod
    def find_out_file(job: SingleJob):
        if job.path is None:
            return None
        out_path = Path(job.path) / f"{job.name}.out"
        if not out_path.exists():
            return None
        return out_path

    @staticmethod
    def stop_calculation(job: AMSJob, file_name="interactive.in", reason="Stop") -> Any:
        if job.path is None:
            return None
        with open(Path(job.path) / file_name, "w") as f:
            f.write(reason)

    @staticmethod
    def stop_has_occurred(job: AMSJob, value="calculation interrupted by user.") -> bool:
        out_path = AMSJobErrorChecker.find_out_file(job)
        if out_path is None:
            return None
        with open(out_path, "r") as f:
            file = f.read()
        return value in file

    @staticmethod
    def status_rkf(job: AMSJob):
        if "ams" not in job.results.rkfs:
            if job.path is None:
                return None
            if (Path(job.path) / f"{job.name}.in").exists():
                return "Not Started"
            return None
        return job.results.rkfs["ams"].read(section="General", variable="termination status")

    @staticmethod
    def status_log(job: AMSJob):
        logfile_path = AMSJobErrorChecker.find_log_file(job)
        if logfile_path is None:
            return "not run"
        with open(logfile_path) as f:
            lines = f.readlines()
        if "NORMAL TERMINATION" in lines[-1]:
            return "normal"
        for line in lines[::-1]:
            if "*** MDStep" in line:
                return int(line.split("*** ")[-1].replace(" ***\n", ""))
        return "Not known"

    @staticmethod
    def status_log_md_step(job: AMSJob):
        logfile_path = AMSJobErrorChecker.find_log_file(job)
        if logfile_path is None:
            return None
        with open(logfile_path) as f:
            lines = f.readlines()
        for line in lines[::-1]:
            if "*** MDStep" in line:
                return int(line.split("*** ")[-1].replace(" ***\n", "").replace("MDStep", ""))
        return None

    @staticmethod
    def time_duration_log(job: AMSJob):
        log_file = AMSJobErrorChecker.find_log_file(job)
        if log_file is None:
            return None

        def timestamp_log(part: str):
            part = part.split(">  ", 1)[0]
            try:
                timestamp = datetime.strptime(part, "<%b%d-%Y> <%H:%M:%S")
                return timestamp
            except ValueError:
                return None

        with open(log_file) as f:
            lines = f.readlines()
            initial = timestamp_log(lines[0])
            final = timestamp_log(lines[-1])
            if initial is None or final is None:
                return None
            duration_seconds = (final - initial).total_seconds()
        return duration_seconds

    @staticmethod
    def time_expected_end_md(job: AMSJob):
        def last_md_step_log(logfile_path: Path):
            with open(logfile_path) as f:
                lines = f.readlines()
            for line in lines[::-1]:
                if "*** MDStep" in line:
                    return int(line.split("*** ")[-1].replace(" ***\n", "").replace("MDStep", ""))
            return None

        log_file = AMSJobErrorChecker.find_log_file(job)
        if log_file is None:
            return None
        step_status = last_md_step_log(log_file)
        if not isinstance(step_status, int):
            return None
        md_end = job.settings.get_nested("input.ams.MolecularDynamics.NSteps".split("."), None)
        if md_end is None:
            return md_end
        md_end = int(md_end)
        time_spent = AMSJobErrorChecker.time_duration_log(job)
        if time_spent is None:
            return None
        time_left = time_spent / step_status * (md_end - step_status)
        date_time = datetime.now() + timedelta(seconds=time_left)
        return date_time

    @staticmethod
    def errors_out(job: AMSJob) -> Any:
        out_path = AMSJobErrorChecker.find_out_file(job)
        if out_path is None:
            return None
        with open(out_path, "r") as file:
            # [line.strip() for line in job.results.grep_file(f'{job.name}.out', 'ERROR: ')] -> TBN MUCH slower!
            error_lines = [line.strip() for line in file if "ERROR" in line]
            error_out = "\n".join(error_lines)
            return error_out
        return None

    @staticmethod
    def errors_log(job: AMSJob) -> Any:
        log_path = AMSJobErrorChecker.find_log_file(job)
        if log_path is None:
            return None
        with open(log_path, "r") as file:
            error_lines = ["ERROR" + line.split("ERROR", 1)[-1].strip() for line in file if "ERROR" in line]
            error_log = "\n".join(error_lines)
            return error_log
        return None

    @staticmethod
    def warnings_out(job: AMSJob) -> Any:
        out_path = AMSJobErrorChecker.find_out_file(job)
        if out_path is None:
            return None
        with open(out_path, "r") as file:
            warning_lines = [line.strip() for line in file if "WARNING" in line]
            warning_out = "\n".join(warning_lines)
            return warning_out
        return None

    @staticmethod
    def warnings_log(job: AMSJob) -> Any:
        log_path = AMSJobErrorChecker.find_log_file(job)
        if log_path is None:
            return None
        with open(log_path, "r") as file:
            warning_lines = ["WARNING" + line.split("WARNING", 1)[-1].strip() for line in file if "WARNING" in line]
            warning_log = "\n".join(warning_lines)
            return warning_log
        return None

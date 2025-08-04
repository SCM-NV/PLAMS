import datetime
from typing import Optional, Sequence, Union, Dict, List, Callable, Any, Tuple, Hashable, Set, Literal, TYPE_CHECKING
import os
import csv
from pathlib import Path
import numpy as np
from numbers import Number
from dataclasses import dataclass, replace
from itertools import chain, islice

from scm.plams.core.basejob import Job, SingleJob
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.core.errors import PlamsError
from scm.plams.core.functions import requires_optional_package, load, config_context
from scm.plams.tools.table_formatter import format_in_table
from scm.plams.interfaces.molecule.rdkit import to_smiles
from scm.plams.mol.molecule import Molecule
from scm.plams.interfaces.adfsuite.inputparser import input_to_settings

try:
    from scm.libbase import UnifiedChemicalSystem as ChemicalSystem
    from scm.utils.conversions import chemsys_to_plams_molecule

    _has_scm_libbase = True
except ImportError:
    _has_scm_libbase = False

try:
    from scm.pisa.block import DriverBlock
    from scm.pisa.input_def import DRIVER_BLOCK_FILES, ENGINE_BLOCK_FILES

    _has_scm_pisa = True
except ImportError:
    _has_scm_pisa = False

if TYPE_CHECKING:
    from pandas import DataFrame


__all__ = ["JobAnalysis"]


class JobAnalysis:
    """
    Analysis tool for Jobs, which generates tables of data consisting of fields and their respective value for each job.

    The jobs and fields which are included in the analysis are customizable, to allow for flexible comparison.
    """

    @dataclass
    class _Field:
        key: str
        value_extractor: Callable[[Job], Any]
        display_name: Optional[str] = None
        fmt: Optional[str] = None
        from_settings: bool = False
        expansion_depth: int = 0

        def __post_init__(self):
            self.display_name = self.key if self.display_name is None else self.display_name

    _standard_fields = {
        "Path": _Field(key="Path", value_extractor=lambda j: j.path),
        "Name": _Field(key="Name", value_extractor=lambda j: j.name),
        "OK": _Field(key="OK", value_extractor=lambda j: j.ok()),
        "Check": _Field(key="Check", value_extractor=lambda j: j.check()),
        "ErrorMsg": _Field(key="ErrorMsg", value_extractor=lambda j: j.get_errormsg()),
        "ParentPath": _Field(key="ParentPath", value_extractor=lambda j: j.parent.path if j.parent else None),
        "ParentName": _Field(key="ParentName", value_extractor=lambda j: j.parent.name if j.parent else None),
        "Formula": _Field(
            key="Formula",
            value_extractor=lambda j: (
                JobAnalysis._mol_formula_extractor(j.molecule) if isinstance(j, SingleJob) else None
            ),
        ),
        "Smiles": _Field(
            key="Smiles",
            value_extractor=lambda j: (
                JobAnalysis._mol_smiles_extractor(j.molecule) if isinstance(j, SingleJob) else None
            ),
        ),
        "GyrationRadius": _Field(
            key="GyrationRadius",
            value_extractor=lambda j: (
                JobAnalysis._mol_gyration_radius_extractor(j.molecule) if isinstance(j, SingleJob) else None
            ),
            fmt=".4f",
        ),
        "CPUTime": _Field(
            key="CPUTime",
            value_extractor=lambda j: JobAnalysis._read_rkf(j, "General", "CPUTime"),
            fmt=".6f",
        ),
        "SysTime": _Field(
            key="SysTime",
            value_extractor=lambda j: JobAnalysis._read_rkf(j, "General", "SysTime"),
            fmt=".6f",
        ),
        "ElapsedTime": _Field(
            key="ElapsedTime",
            value_extractor=lambda j: JobAnalysis._read_rkf(j, "General", "ElapsedTime"),
            fmt=".6f",
        ),
    }

    StandardField = Literal[
        "Path",
        "Name",
        "OK",
        "Check",
        "ErrorMsg",
        "ParentPath",
        "ParentName",
        "Formula",
        "Smiles",
        "GyrationRadius",
        "CPUTime",
        "SysTime",
        "ElapsedTime",
    ]

    @staticmethod
    def _mol_formula_extractor(
        mol: Optional[Union[Molecule, Dict[str, Molecule], "ChemicalSystem", Dict[str, "ChemicalSystem"]]]
    ) -> Optional[str]:
        if isinstance(mol, dict):
            return ", ".join([f"{n}: {JobAnalysis._mol_formula_extractor(m)}" for n, m in mol.items()])
        elif isinstance(mol, Molecule):
            return mol.get_formula()
        elif _has_scm_libbase and isinstance(mol, ChemicalSystem):
            return mol.formula()
        return None

    @staticmethod
    def _mol_smiles_extractor(
        mol: Optional[Union[Molecule, Dict[str, Molecule], "ChemicalSystem", Dict[str, "ChemicalSystem"]]]
    ):
        if isinstance(mol, dict):
            return ", ".join([f"{n}: {JobAnalysis._mol_smiles_extractor(m)}" for n, m in mol.items()])
        elif isinstance(mol, Molecule):
            return to_smiles(mol)
        elif _has_scm_libbase and isinstance(mol, ChemicalSystem):
            return JobAnalysis._mol_smiles_extractor(chemsys_to_plams_molecule(mol))
        return None

    @staticmethod
    def _mol_gyration_radius_extractor(
        mol: Optional[Union[Molecule, Dict[str, Molecule], "ChemicalSystem", Dict[str, "ChemicalSystem"]]]
    ):
        if isinstance(mol, dict):
            return ", ".join([f"{n}: {JobAnalysis._mol_gyration_radius_extractor(m)}" for n, m in mol.items()])
        elif isinstance(mol, Molecule):
            return mol.get_gyration_radius()
        elif _has_scm_libbase and isinstance(mol, ChemicalSystem):
            return JobAnalysis._mol_gyration_radius_extractor(chemsys_to_plams_molecule(mol))
        return None

    _reserved_names = ["_jobs", "_fields", "StandardField", "_standard_fields", "_pisa_programs", "_Field"]

    def __init__(
        self,
        paths: Optional[Sequence[Union[str, os.PathLike]]] = None,
        jobs: Optional[Sequence[Job]] = None,
        loaders: Optional[Sequence[Callable[[str], Job]]] = None,
        standard_fields: Optional[Sequence["JobAnalysis.StandardField"]] = ("Path", "Name", "OK", "Check", "ErrorMsg"),
        await_results: bool = True,
    ):
        """
        Initialize new instance of |JobAnalysis| with a set of jobs.

        .. code:: python

            >>> ja = JobAnalysis(jobs=[job1, job2], standard_fields=["Name", "OK"])
            >>> ja

            | Name  | OK   |
            |-------|------|
            | job_1 | True |
            | job_2 | True |

        :param paths: one or more paths to folders from which to load jobs to add to the analysis
        :param jobs: one or more jobs to add to the analysis
        :param loaders: custom loading functions to generate jobs from a job folder
        :param standard_fields: keys of standard fields to include in analysis, defaults to ``("Path", "Name", "OK", "Check", "ErrorMsg")``
        :param await_results: whether to wait for the results of any passed jobs to finish, defaults to ``True``
        """
        self._jobs: Dict[str, Job] = {}
        self._fields: Dict[str, JobAnalysis._Field] = {}

        if _has_scm_pisa:
            self._pisa_programs = {value: key for key, value in ENGINE_BLOCK_FILES.items()}
            self._pisa_programs.update({value: key for key, value in DRIVER_BLOCK_FILES.items()})

        if jobs:
            for j in jobs:
                self._add_job(j)
                if await_results:
                    j.results.wait()

        if paths:
            for p in paths:
                self._load_job(p, loaders)

        if standard_fields:
            for sf in standard_fields:
                self._add_standard_field(sf)

    def copy(self) -> "JobAnalysis":
        """
        Produce a copy of this analysis with the same jobs and fields.

        .. code:: python

            >>> ja.copy()

            | Name  | OK   |
            |-------|------|
            | job_1 | True |
            | job_2 | True |

        :return: copy of the analysis
        """
        cpy = JobAnalysis()
        cpy._jobs = self.jobs
        cpy._fields = {k: replace(v) for k, v in self._fields.items()}
        return cpy

    @property
    def jobs(self) -> Dict[str, Job]:
        """
        Jobs currently included in analysis.

        .. code:: python

            >>> ja.jobs

            {
                '/path/job1': <scm.plams.interfaces.adfsuite.ams.AMSJob object at 0x1085a13d0>,
                '/path/job2': <scm.plams.interfaces.adfsuite.ams.AMSJob object at 0x15e389970>
            }

        :return: Dictionary of the job path and the |Job|
        """
        return {k: v for k, v in self._jobs.items()}

    @property
    def field_keys(self) -> List[str]:
        """
        Keys of current fields, as they appear in the analysis.

        .. code:: python

            >>> ja.field_keys

            ['Name', 'OK']

        :return: list of field keys
        """
        return [k for k in self._fields]

    def get_analysis(self) -> Dict[str, List]:
        """
        Gets analysis data. This is effectively a table in the form of a dictionary,
        where the keys are the field keys and the values are a list of data for each job.

        .. code:: python

            >>> ja.field_keys

            {
                'Name': ['job1', 'job2'],
                'OK': [True, True]
            }

        :return: analysis data as a dictionary of field keys/lists of job values
        """
        analysis = {col_name: self._get_field_analysis(col_name) for col_name in self._fields}
        analysis = self._expand_analysis(analysis)

        return analysis

    def _expand_analysis(self, analysis: Dict[str, Any]) -> Dict[str, Any]:
        """
        Expand analysis fields, converting individual job rows into multiple rows

        :param analysis: analysis data as a dictionary of field keys/lists of job values
        :return: analysis data as a dictionary of field keys/lists of (multiple) job values
        """

        if analysis.keys() and self._jobs:

            def expand(data, expand_fields):
                expanded_data: Dict[str, list] = {col_name: [] for col_name in data.keys()}
                for i in range(len(data[list(data.keys())[0]])):
                    job_data = {col_name: data[i] for col_name, data in data.items()}
                    valid_expand_fields = {
                        f
                        for f in expand_fields
                        if (
                            isinstance(job_data[f], Sequence)
                            or (isinstance(job_data[f], np.ndarray) and job_data[f].shape != ())
                        )
                        and not isinstance(job_data[f], str)
                    }
                    # Number of rows is the maximum expanded field
                    num_expanded_rows = max([len(job_data[f]) for f in valid_expand_fields], default=1)

                    # Convert multiple values to multiple rows of single values
                    for col_name in data:
                        expanded_data[col_name] += (
                            list(islice(chain(job_data[col_name], [None] * num_expanded_rows), num_expanded_rows))
                            if col_name in valid_expand_fields
                            else [job_data[col_name]] * num_expanded_rows
                        )
                return expanded_data

            # Recursively expand until complete
            depth = 1
            while expand_fields := {
                k for k, f in self._fields.items() if k in analysis.keys() and f.expansion_depth >= depth
            }:
                analysis = expand(analysis, expand_fields)
                depth += 1

        return analysis

    def _get_field_analysis(self, key) -> List:
        """
        Gets analysis data for field with a given key. This gives a list of data for the given field, with a value for each job.

        :param: key of the field
        :return: analysis data as list of job values
        """
        if key not in self._fields:
            raise KeyError(f"Field with key '{key}' is not part of the analysis.")

        value_extractor = self._fields[key].value_extractor

        def safe_value(job: Job):
            try:
                return value_extractor(job)
            except Exception as e:
                return f"ERROR: {str(e)}"

        with config_context() as cfg:
            cfg.log.stdout = 0
            cfg.log.file = 0
            return [safe_value(j) for j in self._jobs.values()]

    @requires_optional_package("pandas")
    def to_dataframe(self) -> "DataFrame":
        """
        Converts analysis data to a dataframe. The column names are the field keys and the column values are the values for each job.
        This method requires the `pandas <https://pandas.pydata.org/docs/index.html>`_ package.

        .. code:: python

            >>> print(ja.to_dataframe())

               Name    OK
            0  job1  True
            1  job2  True

        :return: analysis data as a dataframe
        """
        from pandas import DataFrame

        return DataFrame(self.get_analysis())

    def to_table(
        self,
        max_col_width: int = -1,
        max_rows: int = 30,
        fmt: Literal["markdown", "html", "rst"] = "markdown",
    ) -> str:
        """
        Converts analysis data to a pretty-printed table.

        .. code:: python

            >>> print(ja.to_table())

            | Name  | OK   |
            |-------|------|
            | job_1 | True |
            | job_2 | True |

        :param max_col_width: can be integer positive value or -1, defaults to -1 (no maximum width)
        :param max_rows: can be integer positive value or -1, defaults to 30
        :param fmt: format of the table, either markdown (default), html or rst
        :return: string representation of the table
        """

        def safe_format_value(v, vfmt):
            try:
                return format(v, vfmt)
            except (TypeError, ValueError, AttributeError):
                return str(v)

        def safe_format_values(f, vs):
            vfmt = self._fields[f].fmt
            return [safe_format_value(v, vfmt) for v in vs]

        data = {self._fields[f].display_name: safe_format_values(f, v) for f, v in self.get_analysis().items()}
        return format_in_table(data, max_col_width=max_col_width, max_rows=max_rows, fmt=fmt)

    @requires_optional_package("IPython")
    def display_table(
        self,
        max_col_width: int = -1,
        max_rows: int = 30,
        fmt: Literal["markdown", "html", "rst"] = "markdown",
    ) -> None:
        """
        Converts analysis data to a pretty-printed table which is then displayed using IPython.

        .. code:: python

            >>> ja.display_table()

            | Name  | OK   |
            |-------|------|
            | job_1 | True |
            | job_2 | True |

        :param max_col_width: can be integer positive value or -1, defaults to -1 (no maximum width)
        :param max_rows: can be integer positive value or -1, defaults to 30
        :param fmt: format of the table, either markdown (default), html or rst
        """
        from IPython.display import display, Markdown, HTML

        table = self.to_table(max_col_width=max_col_width, max_rows=max_rows, fmt=fmt)

        if fmt == "markdown":
            display(Markdown(table))
        elif fmt == "html":
            display(HTML(table))
        elif fmt == "rst":
            table = "\n".join([f"    {row}" for row in table.split("\n")])
            display(Markdown(table))

    def to_csv_file(self, path: Union[str, os.PathLike]) -> None:
        """
        Write the analysis to a csv file with the specified path.

        .. code:: python

            >>> ja.to_csv_file("./a.csv")
            >>> with open("./a.csv") as csv:
            >>>     print(csv.read())

            Name,OK
            job1,True
            job2,True

        :param path: path to save the csv file
        """
        data = self.get_analysis()
        keys = list(data.keys())
        num_rows = len(data[keys[0]]) if len(keys) > 0 else 0

        with open(path, mode="w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(keys)
            for i in range(num_rows):
                row = [data[k][i] for k in keys]
                writer.writerow(row)

    def get_timeline(self, max_intervals: int = 5, fmt: Literal["markdown", "html", "rst"] = "markdown") -> str:
        """
        Get depiction of timeline of jobs as they were run.
        Each job is represented as a horizontal bar of symbols, where each symbol indicates a different job status.

        These are as follows:

        * ``created``: ``.``
        * ``started``: ``-``
        * ``registered``: ``+``
        * ``running``: ``=``
        * ``finished``: ``*``
        * ``crashed``: ``x``
        * ``failed``: ``X``
        * ``successful``: ``>``
        * ``copied``: ``#``
        * ``preview``: ``~``
        * ``deleted``: ``!``

        e.g.

        .. code:: python

            >>> print(ja.get_timeline())

            | JobName    | ↓2025-02-03 15:16:52 | ↓2025-02-03 15:17:10 | ↓2025-02-03 15:17:28 | ↓2025-02-03 15:17:46 | ↓2025-02-03 15:18:03 | WaitDuration | RunDuration | TotalDuration |
            |------------|----------------------|----------------------|----------------------|----------------------|----------------------|--------------|-------------|---------------|
            | generate   | ==================== | ==================== | ==================== | ==========>          |                      | 0s           | 1m2s        | 1m2s          |
            | reoptimize |                      |                      |                      |           ====>      |                      | 0s           | 3s          | 3s            |
            | score      |                      |                      |                      |                ===>  |                      | 0s           | 2s          | 2s            |
            | filter     |                      |                      |                      |                   =* | >                    | 0s           | 1s          | 1s            |

        If multiple status changes occur within the same resolution period, the latest will be displayed.

        :param max_intervals: maximum number of datetime intervals to display i.e. the width and resolution of the timeline
        :param fmt: format of the table, either markdown (default) or html
        :return: string representation of timeline as a markdown (default), html or rst table
        """
        # Symbols for various job statuses
        status_symbols = {
            "created": ".",
            "started": "-",
            "registered": "+",
            "running": "=",
            "finished": "*",
            "crashed": "x",
            "failed": "X",
            "successful": ">",
            "copied": "#",
            "preview": "~",
            "deleted": "!",
        }

        # Order jobs by the start time
        job_statuses = [(j, j.status_log) for j in self._jobs.values()]
        ordered_job_statuses = sorted(job_statuses, key=lambda x: x[1][0] if x[1] else (datetime.datetime.max, None))

        # Calculate known start and end time and job durations
        def duration(start, end):
            if not start or not end:
                return "Unknown"
            dur = end - start
            d = dur.days
            h, r = divmod(dur.seconds, 3600)
            m, s = divmod(r, 60)
            dur_fmt = ""
            if d > 0:
                dur_fmt += f"{d}d"
            if h > 0 or d > 0:
                dur_fmt += f"{h}h"
            if m > 0 or h > 0 or d > 0:
                dur_fmt += f"{m}m"
            dur_fmt += f"{s}s"
            return dur_fmt

        # Calculate different durations
        start_time = min([s[0][0] for _, s in ordered_job_statuses if s], default=None)
        end_time = max([s[-1][0] for _, s in ordered_job_statuses if s], default=None)
        durations: Dict[str, List[str]] = {"WaitDuration": [], "RunDuration": [], "TotalDuration": []}
        for _, s in ordered_job_statuses:
            wait_duration = "Unknown"
            run_duration = "Unknown"
            total_duration = "Unknown"
            if s and len(s) > 1:
                statuses = [x[1] for x in s]
                if "created" in statuses:
                    created_idx = statuses.index("created")
                    if created_idx < len(statuses) - 1:
                        wait_duration = duration(s[created_idx][0], s[created_idx + 1][0])
                if "started" in statuses and "running" in statuses:
                    started_idx = statuses.index("started")
                    running_idx = statuses.index("running")
                    if started_idx < running_idx < len(statuses) - 1:
                        run_duration = duration(s[started_idx][0], s[running_idx + 1][0])
                total_duration = duration(s[0][0], s[-1][0])
            durations["WaitDuration"].append(wait_duration)
            durations["RunDuration"].append(run_duration)
            durations["TotalDuration"].append(total_duration)

        # Table data
        data = {}
        data["JobName"] = [j.name for j, _ in ordered_job_statuses]

        if start_time and end_time:
            # Calculate the column interval widths
            for num_intervals in range(max_intervals, 1, -1):
                interval = (end_time - start_time) / (num_intervals - 1)
                intervals = [start_time + i * interval for i in range(num_intervals)]
                str_intervals = [f"↓{intv.strftime('%Y-%m-%d %H:%M:%S')}" for intv in intervals]
                if len(set(str_intervals)) == len(str_intervals):
                    break

            num_positions = 20
            symbol_interval = interval / num_positions

            def get_col_and_position(status_time):
                col = int((status_time - start_time) // interval)
                pos = int((status_time - intervals[col]) // symbol_interval)
                if pos == num_positions:
                    col = col + 1
                    pos = 0
                return col, pos

            job_timelines = []
            use_html = fmt == "html"
            for job, statuses in ordered_job_statuses:
                job_timeline = [
                    ["&nbsp;" if use_html else " " for _ in range(num_positions)] for _ in range(num_intervals)
                ]
                if statuses:
                    for i, status in enumerate(statuses):
                        symbol = status_symbols.get(status[1], "?")
                        if i == (len(statuses) - 1):
                            col, pos = get_col_and_position(status[0])
                            job_timeline[col][pos] = symbol
                        else:
                            next_status = statuses[i + 1]
                            col_start, pos_start = get_col_and_position(status[0])
                            col_end, pos_end = get_col_and_position(next_status[0])
                            for col in range(col_start, col_end + 1):
                                for pos in range(
                                    pos_start if col == col_start else 0, pos_end if col == col_end else num_positions
                                ):
                                    job_timeline[col][pos] = symbol

                job_timelines.append(job_timeline)

            for i in range(num_intervals):
                data[str_intervals[i]] = ["".join(jt[i]) for jt in job_timelines]

        # Add durations
        data["WaitDuration"] = durations["WaitDuration"]
        data["RunDuration"] = durations["RunDuration"]
        data["TotalDuration"] = durations["TotalDuration"]

        return (
            format_in_table(data, max_rows=-1, fmt=fmt, monospace=True)
            .replace("<th>", '<th style="border-left: 1px solid black; border-right: 1px solid black;">')
            .replace("<td>", '<td style="border-left: 1px solid black; border-right: 1px solid black;">')
        )

    @requires_optional_package("IPython")
    def display_timeline(self, max_intervals: int = 5, fmt: Literal["markdown", "html", "rst"] = "markdown") -> None:
        """
        Get depiction of timeline of jobs as they were run and display using IPython.
        Each job is represented as a horizontal bar of symbols, where each symbol indicates a different job status.

        These are as follows:

        * ``created``: ``.``
        * ``started``: ``-``
        * ``registered``: ``+``
        * ``running``: ``=``
        * ``finished``: ``*``
        * ``crashed``: ``x``
        * ``failed``: ``X``
        * ``successful``: ``>``
        * ``copied``: ``#``
        * ``preview``: ``~``
        * ``deleted``: ``!``

        e.g.

        .. code:: python

            >>> ja.display_timeline()

            | JobName    | ↓2025-02-03 15:16:52 | ↓2025-02-03 15:17:10 | ↓2025-02-03 15:17:28 | ↓2025-02-03 15:17:46 | ↓2025-02-03 15:18:03 | WaitDuration | RunDuration | TotalDuration |
            |------------|----------------------|----------------------|----------------------|----------------------|----------------------|--------------|-------------|---------------|
            | generate   | ==================== | ==================== | ==================== | ==========>          |                      | 0s           | 1m2s        | 1m2s          |
            | reoptimize |                      |                      |                      |           ====>      |                      | 0s           | 3s          | 3s            |
            | score      |                      |                      |                      |                ===>  |                      | 0s           | 2s          | 2s            |
            | filter     |                      |                      |                      |                   =* | >                    | 0s           | 1s          | 1s            |

        If multiple status changes occur within the same resolution period, the latest will be displayed.
        :param max_intervals: maximum number of datetime intervals to display i.e. the width and resolution of the timeline
        :param fmt: format of the table, either markdown (default), html or rst
        """
        from IPython.display import display, Markdown, HTML

        table = self.get_timeline(max_intervals=max_intervals, fmt=fmt)

        if fmt == "markdown":
            display(Markdown(table))
        elif fmt == "html":
            display(HTML(table))
        elif fmt == "rst":
            table = "\n".join([f"    {row}" for row in table.split("\n")])
            display(Markdown(table))

    def _add_job(self, job: Job):
        """
        Add a job to this |JobAnalysis| instance.

        :param job: |Job| to add to the analysis
        """
        if job.path in self._jobs:
            raise KeyError(f"Job with path '{job.path}' has already been added to the analysis.")
        self._jobs[job.path] = job

    def add_job(self, job: Job) -> "JobAnalysis":
        """
        Add a job to the analysis. This adds a row to the analysis data.

        .. code:: python

            >>> ja.add_job(job3)

            | Name  | OK   |
            |-------|------|
            | job_1 | True |
            | job_2 | True |
            | job_3 | True |

        :param job: |Job| to add to the analysis
        :return: copy of |JobAnalysis| with the additional job
        """
        cpy = self.copy()
        cpy._add_job(job)
        return cpy

    def remove_job(self, job: Union[str, os.PathLike, Job]) -> "JobAnalysis":
        """
        Remove a job from the analysis. This removes a row from the analysis data.

        .. code:: python

            >>> ja.remove_job(job2)

            | Name  | OK   |
            |-------|------|
            | job_1 | True |

        :param job: |Job| or path to a job to remove from the analysis
        :return: copy of |JobAnalysis| with the job removed
        """
        path = job.path if isinstance(job, Job) else str(os.path.abspath(job))
        if path not in self._jobs:
            raise KeyError(f"Job with path '{path}' is not part of the analysis.")

        cpy = self.copy()
        cpy._jobs.pop(path)
        return cpy

    def _load_job(self, path: Union[str, os.PathLike], loaders: Optional[Sequence[Callable[[str], Job]]] = None):
        """
        Add job to this |JobAnalysis| instance by loading from a given path to the job folder.
        If no dill file is present in that location, or the dill unpickling fails, the loaders will be used to load the given job from the folder.

        :param path: path to folder from which to load the job
        :param loaders: functions to try and load jobs, defaults to :meth:`~scm.plams.interfaces.adfsuite.ams.AMSJob.load_external` followed by |load_external|
        """
        path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"Cannot find job file in location '{path}'")

        dill_file = path / f"{path.name}.dill"

        job = None
        loaders = (
            loaders
            if loaders
            else [
                AMSJob.load_external,
                SingleJob.load_external,
            ]
        )

        use_loaders = not dill_file.exists()
        if not use_loaders:
            try:
                job = load(dill_file)
            except Exception:
                use_loaders = True

        if use_loaders:
            for loader in loaders:
                try:
                    job = loader(str(path))
                    break
                except Exception:
                    pass

        if not job:
            raise PlamsError(f"Could not load job from path '{path}'")

        self._add_job(job)

    def load_job(
        self, path: Union[str, os.PathLike], loaders: Optional[Sequence[Callable[[str], Job]]] = None
    ) -> "JobAnalysis":
        """
        Add job to the analysis by loading from a given path to the job folder.
        If no dill file is present in that location, or the dill unpickling fails, the loaders will be used to load the given job from the folder.

        .. code:: python

            >>> ja.load_job("path/job3")

            | Name  | OK   |
            |-------|------|
            | job_1 | True |
            | job_2 | True |
            | job_3 | True |

        :param path: path to folder from which to load the job
        :param loaders: functions to try and load jobs, defaults to :meth:`~scm.plams.interfaces.adfsuite.ams.AMSJob.load_external` followed by |load_external|
        :return: copy of |JobAnalysis| with the additional loaded job
        """
        cpy = self.copy()
        cpy._load_job(path, loaders)
        return cpy

    def filter_jobs(self, predicate: Callable[[Dict[str, Any]], bool]) -> "JobAnalysis":
        """
        Retain jobs from the analysis where the given predicate for field values evaluates to ``True``.
        In other words, this removes rows(s) from the analysis data where the filter function evaluates to ``False`` given a dictionary of the row data.

        .. code:: python

            >>> ja

            | Name  | OK    |
            |-------|-------|
            | job_1 | True  |
            | job_2 | True  |
            | job_3 | False |

            >>> ja.filter_jobs(lambda data: not data["OK"])

            | Name  | OK    |
            |-------|-------|
            | job_3 | False |

        :param predicate: filter function which takes a dictionary of field keys and their values and evaluates to ``True``/``False``
        :return: copy of |JobAnalysis| with the job filter applied
        """
        # Need to make sure a path is associated with job row(s) as this is the job key
        cpy = self.copy()
        requires_path = "Path" not in cpy
        try:
            if requires_path:
                cpy = cpy.add_standard_field("Path")

            analysis = cpy.get_analysis()
            jobs_to_remove = set()
            for i, p in enumerate(analysis["Path"]):
                data = {k: v[i] for k, v in analysis.items()}
                if not predicate(data):
                    jobs_to_remove.add(p)

            for j in jobs_to_remove:
                cpy = cpy.remove_job(j)
        finally:
            if requires_path:
                cpy = cpy.remove_field("Path")

        return cpy

    def sort_jobs(
        self,
        field_keys: Optional[Sequence[str]] = None,
        sort_key: Optional[Callable[[Dict[str, Any]], Any]] = None,
        reverse: bool = False,
    ) -> "JobAnalysis":
        """
        Sort jobs according to a single or multiple fields. This is the order the rows will appear in the analysis data.

        Either one of ``field_keys`` or ``key`` must be provided.
        If ``field_keys`` is provided, the values from these field(s) will be used to sort, in the order they are specified.
        If ``sort_key`` is provided, the sorting function will be applied to all fields.

        .. code:: python

            >>> ja.sort_jobs(field_keys=["Name"], reverse=True)

            | Name  | OK    |
            |-------|-------|
            | job_2 | True  |
            | job_1 | True  |

        :param field_keys: field keys to sort by,
        :param sort_key: sorting function which takes a dictionary of field keys and their values
        :param reverse: reverse sort order, defaults to ``False``
        :return: copy of |JobAnalysis| with the job sorted
        """
        cpy = self.copy()
        analysis = cpy.get_analysis()
        sort_key = sort_key if sort_key else lambda data: tuple([str(v) for v in data.values()])
        key_set = set(field_keys) if field_keys else set(cpy.field_keys)

        def key(ik):
            i, _ = ik
            return sort_key({k: v[i] for k, v in analysis.items() if k in key_set})

        sorted_keys = sorted(enumerate(cpy._jobs.keys()), key=key, reverse=reverse)
        cpy._jobs = {k: cpy._jobs[k] for _, k in sorted_keys}
        return cpy

    def _add_field(
        self,
        key: str,
        value_extractor: Callable[[Job], Any],
        display_name: Optional[str] = None,
        fmt: Optional[str] = None,
        expansion_depth: int = 0,
    ):
        """
        Add a new field to this |JobAnalysis| instance. This adds a column to the analysis data.

        :param key: unique identifier for the field
        :param value_extractor: callable to extract the value for the field from a job
        :param display_name: name which will appear for the field when displayed in table
        :param fmt: string format for how field values are displayed in table
        :param expansion_depth: whether to expand field of multiple values into multiple rows, and recursively to what depth
        """
        if key in self._fields:
            raise KeyError(f"Field with key '{key}' has already been added to the analysis.")

        self._set_field(
            key=key,
            value_extractor=value_extractor,
            display_name=display_name,
            fmt=fmt,
            expansion_depth=expansion_depth,
        )

    def add_field(
        self,
        key: str,
        value_extractor: Callable[[Job], Any],
        display_name: Optional[str] = None,
        fmt: Optional[str] = None,
        expansion_depth: int = 0,
    ) -> "JobAnalysis":
        """
        Add a new field to the analysis. This adds a column to the analysis data.

        .. code:: python

            >>> ja.add_field("N", lambda j: len(j.molecule), display_name="Num Atoms")

            | Name  | OK    | Num Atoms |
            |-------|-------|-----------|
            | job_1 | True  | 4         |
            | job_2 | True  | 6         |

        :param key: unique identifier for the field
        :param value_extractor: callable to extract the value for the field from a job
        :param display_name: name which will appear for the field when displayed in table
        :param fmt: string format for how field values are displayed in table
        :param expansion_depth: whether to expand field of multiple values into multiple rows, and recursively to what depth
        :return: copy of |JobAnalysis| with the additional field
        """
        cpy = self.copy()
        cpy._add_field(
            key=key,
            value_extractor=value_extractor,
            display_name=display_name,
            fmt=fmt,
            expansion_depth=expansion_depth,
        )
        return cpy

    def _set_field(
        self,
        key: str,
        value_extractor: Callable[[Job], Any],
        display_name: Optional[str] = None,
        fmt: Optional[str] = None,
        expansion_depth: int = 0,
    ):
        """
        Set a field in this |JobAnalysis| instance. This adds or modifies a column in the analysis data.

        .. code:: python

            >>> ja.set_field("N", lambda j: len(j.molecule), display_name="Num Atoms")

            | Name  | OK    | Num Atoms |
            |-------|-------|-----------|
            | job_1 | True  | 4         |
            | job_2 | True  | 6         |

        :param key: unique identifier for the field
        :param value_extractor: callable to extract the value for the field from a job
        :param display_name: name which will appear for the field when displayed in table
        :param fmt: string format for how field values are displayed in table
        :param expansion_depth: whether to expand field of multiple values into multiple rows, and recursively to what depth
        """
        self._fields[key] = self._Field(
            key=key,
            value_extractor=value_extractor,
            display_name=display_name,
            fmt=fmt,
            expansion_depth=expansion_depth,
        )

    def set_field(
        self,
        key: str,
        value_extractor: Callable[[Job], Any],
        display_name: Optional[str] = None,
        fmt: Optional[str] = None,
        expansion_depth: int = 0,
    ) -> "JobAnalysis":
        """
        Set a field in the analysis. This adds or modifies a column in the analysis data.

        .. code:: python

            >>> ja.set_field("N", lambda j: len(j.molecule), display_name="Num Atoms")

            | Name  | OK    | Num Atoms |
            |-------|-------|-----------|
            | job_1 | True  | 4         |
            | job_2 | True  | 6         |

        :param key: unique identifier for the field
        :param value_extractor: callable to extract the value for the field from a job
        :param display_name: name which will appear for the field when displayed in table
        :param fmt: string format for how field values are displayed in table
        :param expansion_depth: whether to expand field of multiple values into multiple rows, and recursively to what depth
        :return: copy of |JobAnalysis| with the column set
        """
        cpy = self.copy()
        cpy._set_field(
            key=key,
            value_extractor=value_extractor,
            display_name=display_name,
            fmt=fmt,
            expansion_depth=expansion_depth,
        )
        return cpy

    def format_field(self, key: str, fmt: Optional[str] = None) -> "JobAnalysis":
        """
        Apply a string formatting to a given field. This will apply when ``to_table`` is called.

        .. code:: python

            >>> ja.format_field("N", "03.0f")

            | Name  | OK    | Num Atoms |
            |-------|-------|-----------|
            | job_1 | True  | 004       |
            | job_2 | True  | 006       |

        :param key: unique identifier of the field
        :param fmt: string format of the field e.g. ``.2f``
        :return: copy of |JobAnalysis| with field formatting applied
        """
        if key not in self._fields:
            raise KeyError(f"Field with key '{key}' is not part of the analysis.")

        cpy = self.copy()
        cpy._fields[key] = replace(cpy._fields[key], fmt=fmt)
        return cpy

    def rename_field(self, key: str, display_name: str) -> "JobAnalysis":
        """
        Give a display name to a field in the analysis. This is the header of the column in the analysis data.

        .. code:: python

            >>> ja.rename_field("N", "N Atoms")

            | Name  | OK    | N Atoms |
            |-------|-------|---------|
            | job_1 | True  | 004     |
            | job_2 | True  | 006     |

        :param key: unique identifier for the field
        :param display_name: name of the field
        :return: copy of |JobAnalysis| with the renamed field
        """

        if key not in self._fields:
            raise KeyError(f"Field with key '{key}' is not part of the analysis.")

        cpy = self.copy()
        cpy._fields[key] = replace(cpy._fields[key], display_name=display_name)
        return cpy

    def expand_field(self, key: str, depth: int = 1) -> "JobAnalysis":
        """
        Expand field of multiple values into multiple rows for each job.
        For nested values, the depth can be provided to determine the level of recursive expansion.

        .. code:: python

            >>> (ja
            >>>  .add_field("Step", lambda j: get_steps(j))
            >>>  .add_field("Energy", lambda j: get_energies(j)))

            | Name  | OK    | Step      | Energy             |
            |-------|-------|-----------|--------------------|
            | job_1 | True  | [1, 2, 3] | [42.1, 43.2, 42.5] |
            | job_2 | True  | [1, 2]    | [84.5, 112.2]      |

            >>> (ja
            >>>  .expand_field("Step")
            >>>  .expand_field("Energy"))

            | Name  | OK    | Step | Energy |
            |-------|-------|------|--------|
            | job_1 | True  | 1    | 42.1   |
            | job_1 | True  | 2    | 43.2   |
            | job_1 | True  | 3    | 42.5   |
            | job_2 | True  | 1    | 84.5   |
            | job_2 | True  | 1    | 112.2  |

        :param key: unique identifier of field to expand
        :param depth: depth of recursive expansion, defaults to 1
        :return: copy of |JobAnalysis| with the field expanded
        """
        if key not in self._fields:
            raise KeyError(f"Field with key '{key}' is not part of the analysis.")

        cpy = self.copy()
        cpy._fields[key].expansion_depth = depth
        return cpy

    def collapse_field(self, key: str) -> "JobAnalysis":
        """
        Collapse field of multiple rows into single row of multiple values for each job.

        .. code:: python

            >>> ja

            | Name  | OK    | Step | Energy |
            |-------|-------|------|--------|
            | job_1 | True  | 1    | 42.1   |
            | job_1 | True  | 2    | 43.2   |
            | job_1 | True  | 3    | 42.5   |
            | job_2 | True  | 1    | 84.5   |
            | job_2 | True  | 1    | 112.2  |

            >>> (ja
            >>>  .collapse_field("Step")
            >>>  .collapse_field("Energy"))

            | Name  | OK    | Step      | Energy             |
            |-------|-------|-----------|--------------------|
            | job_1 | True  | [1, 2, 3] | [42.1, 43.2, 42.5] |
            | job_2 | True  | [1, 2]    | [84.5, 112.2]      |

        :param key: unique identifier of field to collapse
        :return: copy of |JobAnalysis| with the field collapsed applied
        """
        if key not in self._fields:
            raise KeyError(f"Field with key '{key}' is not part of the analysis.")

        cpy = self.copy()
        cpy._fields[key].expansion_depth = 0
        return cpy

    def reorder_fields(self, order: Sequence[str]) -> "JobAnalysis":
        """
        Reorder fields based upon the given sequence of field keys. This is the order the columns will appear in the analysis data.

        Any specified fields will be placed first, with remaining fields placed after with their order unchanged.

        .. code:: python

            >>> ja.reorder_fields(["Name", "Step"])

            | Name  | Step | OK    | Energy |
            |-------|------|-------|--------|
            | job_1 | 1    | True  | 42.1   |
            | job_1 | 2    | True  | 43.2   |
            | job_1 | 3    | True  | 42.5   |
            | job_2 | 1    | True  | 84.5   |
            | job_2 | 1    | True  | 112.2  |

        :param order: sequence of fields to be placed at the start of the field ordering
        :return: copy of |JobAnalysis| with the fields reordered
        """

        def key(field_key):
            try:
                return order.index(field_key)
            except ValueError:
                return len(order)

        return self.sort_fields(sort_key=key)

    def sort_fields(self, sort_key: Callable[[str], Any], reverse: bool = False) -> "JobAnalysis":
        """
        Sort fields according to a sort key. This is the order the columns will appear in the analysis data.

        .. code:: python

            >>> ja.sort_fields(lambda k: len(k))

            | OK    | Name  | Step | Energy |
            |-------|-------|------|--------|
            | True  | job_1 | 1    | 42.1   |
            | True  | job_1 | 2    | 43.2   |
            | True  | job_1 | 3    | 42.5   |
            | True  | job_2 | 1    | 84.5   |
            | True  | job_2 | 1    | 112.2  |

        :param sort_key: sorting function which accepts the field key
        :param reverse: reverse sort order, defaults to ``False``
        :return: copy of |JobAnalysis| with the fields sorted
        """
        cpy = self.copy()
        sorted_keys = sorted(cpy._fields.keys(), key=sort_key, reverse=reverse)
        cpy._fields = {k: cpy._fields[k] for k in sorted_keys}
        return cpy

    def _remove_field(self, key: str):
        """
        Remove a field from this |JobAnalysis| instance. This removes a column from the analysis data.

        :param key: unique identifier of the field
        """
        if key not in self._fields:
            raise KeyError(f"Field with key '{key}' is not part of the analysis.")

        self._fields.pop(key)

    def remove_field(self, key: str) -> "JobAnalysis":
        """
        Remove a field from the analysis. This removes a column from the analysis data.

        .. code:: python

            >>> ja.remove_field("OK")

            | Name  | OK   |
            |-------|------|
            | job_1 | True |
            | job_2 | True |

        :param key: unique identifier of the field
        :return: copy of |JobAnalysis| with the field removed
        """
        cpy = self.copy()
        cpy._remove_field(key)
        return cpy

    def remove_fields(self, keys: Sequence[str]) -> "JobAnalysis":
        """
        Remove multiple fields from the analysis. This removes columns from the analysis data.

        .. code:: python

            >>> ja.remove_fields(["OK", "N"])

            | Name  |
            |-------|
            | job_1 |
            | job_2 |

        :param keys: unique identifiers of the fields
        :return: copy of |JobAnalysis| with the fields removed
        """
        cpy = self.copy()
        for key in keys:
            cpy._remove_field(key)
        return cpy

    def filter_fields(self, predicate: Callable[[List[Any]], bool]) -> "JobAnalysis":
        """
        Retain fields from the analysis where the given predicate evaluates to ``True`` given the field values.
        In other words, this removes column(s) from the analysis data where the filter function evaluates to ``False``
        given all the row values.

        .. code:: python

            >>> ja

            | OK    | Name  | Step | Energy |
            |-------|-------|------|--------|
            | True  | job_1 | 1    | 42.1   |
            | True  | job_1 | 2    | 43.2   |
            | True  | job_1 | 3    | 42.5   |
            | True  | job_2 | 1    | 84.5   |
            | True  | job_2 | 1    | 112.2  |

            >>> ja.filter_fields(lambda vals: all([not isinstance(v, int) or v > 50 for v in vals]))

            | Name  | Energy |
            |-------|--------|
            | job_1 | 42.1   |
            | job_1 | 43.2   |
            | job_1 | 42.5   |
            | job_2 | 84.5   |
            | job_2 | 112.2  |

        :param predicate: filter function which takes values and evaluates to ``True``/``False``
        :return: copy of |JobAnalysis| with the field filter applied
        """
        cpy = self.copy()
        for n, vals in cpy.get_analysis().items():
            if not predicate(vals):
                cpy._remove_field(n)
        return cpy

    def remove_empty_fields(self) -> "JobAnalysis":
        """
        Remove field(s) from the analysis which have ``None`` for all values. This removes column(s) from the analysis data,
        where all rows have empty values.

        .. code:: python

            >>> ja.add_standard_field("ParentName")

            | Name  | OK    | ParentName |
            |-------|-------|------------|
            | job_1 | True  | None       |
            | job_2 | True  | None       |

            >>> ja.remove_empty_fields()

            | Name  | OK    |
            |-------|-------|
            | job_1 | True  |
            | job_2 | True  |

        :return: copy of |JobAnalysis| with empty fields removed
        """
        return self.filter_fields(lambda vals: any([not self._is_empty_value(v) for v in vals]))

    @staticmethod
    def _is_empty_value(val) -> bool:
        """
        Check if a value is considered empty i.e. is ``None`` or has no value.
        """
        if val is None:
            return True

        if isinstance(val, np.ndarray):
            return val.shape == () and val.item() is None
        elif isinstance(val, (Sequence, Dict)):
            return len(val) == 0

        return False

    def remove_uniform_fields(self, tol: float = 1e-08, ignore_empty: bool = False) -> "JobAnalysis":
        """
        Remove field(s) from the analysis which evaluate the same for all values. This removes column(s) from the analysis data,
        where all rows have the same value.

        .. code:: python

            >>> ja.add_standard_field("ParentName")

            | Name  | OK    | ParentName |
            |-------|-------|------------|
            | job_1 | True  | None       |
            | job_2 | True  | None       |
            | job_3 | True  | p_job_4    |

            >>> ja.remove_uniform_fields()

            | Name  | ParentName |
            |-------|------------|
            | job_1 | None       |
            | job_2 | None       |
            | job_3 | p_job_4    |

            >>> ja.remove_uniform_fields(ignore_empty=True)

            | Name  |
            |-------|
            | job_1 |
            | job_2 |
            | job_3 |

        :param tol: absolute tolerance for numeric value comparison, all values must fall within this range
        :param ignore_empty: when ``True`` ignore ``None`` values and empty containers in comparison, defaults to ``False``
        :return: copy of |JobAnalysis| with uniform fields removed
        """

        def is_uniform(vals: List[Any]):
            """
            Check if a list of values is considered uniform
            """
            # Skip over None values if set to be ignored
            vals = [v for v in vals if not self._is_empty_value(v) or not ignore_empty]
            # Empty list is uniform
            if not vals:
                return True
            # Check if all numeric values, and if so evaluate range within tolerance
            if all([isinstance(v, Number) and not isinstance(v, bool) for v in vals]):
                return np.ptp(vals) <= tol

            # Check if all iterable values, and if so evaluate elements individually
            if all(
                [
                    (
                        (isinstance(v, Sequence) and not isinstance(v, str))
                        or (isinstance(v, np.ndarray) and v.shape != ())
                    )
                    for v in vals
                ]
            ):
                l = len(vals[0])
                if any(len(v) != l for v in vals):
                    return False
                return all([is_uniform([v[i] for v in vals]) for i in range(l)])

            # Check if dictionary, and if so evaluate keys and values are uniform
            if all([isinstance(v, dict) for v in vals]):
                ks = set(vals[0].keys())
                if any(set(v.keys()) != ks for v in vals):
                    return False
                for k in ks:
                    if not all([is_uniform([v[k] for v in vals])]):
                        return False
                return True

            try:
                return all([v == vals[0] for v in vals])
            except ValueError:
                return False

        return self.filter_fields(lambda vals: not is_uniform(vals))

    def add_standard_fields(self, keys: Sequence["JobAnalysis.StandardField"]) -> "JobAnalysis":
        """
        Adds multiple standard fields to the analysis.

        These are:

        * ``Path``: for |Job| attribute :attr:`~scm.plams.core.basejob.Job.path`
        * ``Name``: for |Job| attribute :attr:`~scm.plams.core.basejob.Job.name`
        * ``OK``: for |Job| method :meth:`~scm.plams.core.basejob.Job.ok`
        * ``Check``: for |Job| method :meth:`~scm.plams.core.basejob.Job.check`
        * ``ErrorMsg``: for |Job| method :meth:`~scm.plams.core.basejob.Job.get_errormsg`
        * ``ParentPath``: for attribute :attr:`~scm.plams.core.basejob.Job.path` of |Job| attribute :attr:`~scm.plams.core.basejob.Job.parent`
        * ``ParentName``: for attribute :attr:`~scm.plams.core.basejob.Job.name` of |Job| attribute :attr:`~scm.plams.core.basejob.Job.parent`
        * ``Formula``: for method :meth:`~scm.plams.mol.molecule.Molecule.get_formula` of |Job| attribute :attr:`~scm.plams.core.basejob.SingleJob.molecule`
        * ``Smiles``: for function :func:`~scm.plams.interfaces.molecule.rdkit.to_smiles` for |Job| attribute :attr:`~scm.plams.core.basejob.SingleJob.molecule`
        * ``GyrationRadius``: for function :meth:`~scm.plams.mol.molecule.Molecule.gyration_radius` for |Job| attribute :attr:`~scm.plams.core.basejob.SingleJob.molecule`
        * ``CPUTime``: for method :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.readrkf` with ``General/CPUTime`` for |Job| attribute :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results`
        * ``SysTime``: for method :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.readrkf` with ``General/SysTime`` for |Job| attribute :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results`
        * ``ElapsedTime``: for method :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.readrkf` with ``General/ElapsedTime`` for |Job| attribute :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results`

        .. code:: python

            >>> ja

            | Name  |
            |-------|
            | job_1 |
            | job_2 |

            >>> ja.add_standard_fields(["Path", "Smiles"])

            | Name  | Path        | Smiles |
            |-------|-------------|--------|
            | job_1 | /path/job_1 | N      |
            | job_2 | /path/job_2 | C=C    |

        :param keys: sequence of keys for the analysis fields
        :return: copy of |JobAnalysis| with the standard fields added
        """
        cpy = self.copy()
        for key in keys:
            cpy._add_standard_field(key)
        return cpy

    def _add_standard_field(self, key: "JobAnalysis.StandardField"):
        """
        Adds a standard field to this |JobAnalysis| instance.

        :param key: key for the analysis field
        """
        if key not in self._standard_fields:
            raise KeyError(f"'{key}' is not one of the standard fields: {', '.join(self._standard_fields)}.")

        if key in self._fields:
            raise KeyError(f"Field with key '{key}' has already been added to the analysis.")

        self._fields[key] = replace(self._standard_fields[key])

    def add_standard_field(self, key: "JobAnalysis.StandardField") -> "JobAnalysis":
        """
        Adds a standard field to the analysis.

        These are:

        * ``Path``: for |Job| attribute :attr:`~scm.plams.core.basejob.Job.path`
        * ``Name``: for |Job| attribute :attr:`~scm.plams.core.basejob.Job.name`
        * ``OK``: for |Job| method :meth:`~scm.plams.core.basejob.Job.ok`
        * ``Check``: for |Job| method :meth:`~scm.plams.core.basejob.Job.check`
        * ``ErrorMsg``: for |Job| method :meth:`~scm.plams.core.basejob.Job.get_errormsg`
        * ``ParentPath``: for attribute :attr:`~scm.plams.core.basejob.Job.path` of |Job| attribute :attr:`~scm.plams.core.basejob.Job.parent`
        * ``ParentName``: for attribute :attr:`~scm.plams.core.basejob.Job.name` of |Job| attribute :attr:`~scm.plams.core.basejob.Job.parent`
        * ``Formula``: for method :meth:`~scm.plams.mol.molecule.Molecule.get_formula` of |Job| attribute :attr:`~scm.plams.core.basejob.SingleJob.molecule`
        * ``Smiles``: for function :func:`~scm.plams.interfaces.molecule.rdkit.to_smiles` for |Job| attribute :attr:`~scm.plams.core.basejob.SingleJob.molecule`
        * ``GyrationRadius``: for function :meth:`~scm.plams.mol.molecule.Molecule.gyration_radius` for |Job| attribute :attr:`~scm.plams.core.basejob.SingleJob.molecule`
        * ``CPUTime``: for method :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.readrkf` with ``General/CPUTime`` for |Job| attribute :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results`
        * ``SysTime``: for method :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.readrkf` with ``General/SysTime`` for |Job| attribute :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results`
        * ``ElapsedTime``: for method :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.readrkf` with ``General/ElapsedTime`` for |Job| attribute :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results`

        .. code:: python

            >>> ja

            | Name  |
            |-------|
            | job_1 |
            | job_2 |

            >>> ja.add_standard_field("Path")

            | Name  | Path        |
            |-------|-------------|
            | job_1 | /path/job_1 |
            | job_2 | /path/job_2 |

        :param key: key for the analysis field
        :return: copy of |JobAnalysis| with the additional standard field
        """
        cpy = self.copy()
        cpy._add_standard_field(key)
        return cpy

    @staticmethod
    def get_settings_field_key(key_tuple: Tuple[Hashable, ...]) -> str:
        """
        Gets the key for a |JobAnalysis| settings field corresponding to a given nested |Settings| key.
        This will be a Pascal-case string of the settings nested key path e.g.``("input", "ams", "task")`` will give the key ``InputAmsTask``.

        .. code:: python

            >>> JobAnalysis.get_settings_field_key(("input", "ams", "task"))

            InputAmsTask

        :param key_tuple: nested tuple of keys in the |Settings| object
        :return: corresponding key used for the |JobAnalysis| field
        """
        return "".join([str(k).title() for k in key_tuple])

    def add_settings_field(
        self,
        key_tuple: Tuple[Hashable, ...],
        display_name: Optional[str] = None,
        fmt: Optional[str] = None,
        expansion_depth: int = 0,
    ) -> "JobAnalysis":
        """
        Add a field for a nested key from the job settings to the analysis.
        The key of the field will be a Pascal-case string of the settings nested key path e.g. ``("input", "ams", "task")`` will appear as field ``InputAmsTask``.

        .. code:: python

            >>> ja.add_settings_field(("input", "ams", "task"), display_name="Task")

            | Name  | Task        |
            |-------|-------------|
            | job_1 | SinglePoint |
            | job_2 | SinglePoint |

        :param key_tuple: nested tuple of keys in the settings object
        :param display_name: name which will appear for the field when displayed in table
        :param fmt: string format for how field values are displayed in table
        :param expansion_depth: whether to expand field of multiple values into multiple rows, and recursively to what depth
        :return: copy of |JobAnalysis| with the settings fields added
        """
        key = JobAnalysis.get_settings_field_key(key_tuple)
        cpy = self.copy()
        cpy = cpy.add_field(
            key,
            lambda j, k=key_tuple: cpy._get_job_settings(j).get_nested(k),  # type: ignore
            display_name=display_name,
            fmt=fmt,
            expansion_depth=expansion_depth,
        )
        cpy._fields[key].from_settings = True
        return cpy

    def add_settings_fields(
        self,
        predicate: Optional[Callable[[Tuple[Hashable, ...]], bool]] = None,
        flatten_list: bool = True,
    ) -> "JobAnalysis":
        """
        Add a field for all nested keys which satisfy the predicate from the job settings to the analysis.
        The key of the fields will be a Pascal-case string of the settings nested key path e.g. ("input", "ams", "task") will appear as field ``InputAmsTask``.

        .. code:: python

            >>> ja.add_settings_fields(lambda k: len(k) >= 3 and k[2].lower() == "xc")

            | Name  | InputAdfXcDispersion | InputAdfXcGga |
            |-------|----------------------|---------------|
            | job_1 | Grimme3              | PBE           |
            | job_2 | Grimme3              | PBE           |

        :param predicate: optional predicate which evaluates to ``True`` or ``False`` given a nested key, by default will be ``True`` for every key
        :param flatten_list: whether to flatten lists in settings objects
        :return: copy of |JobAnalysis| with settings fields added
        """

        all_blocks: Set[Tuple[Hashable, ...]] = set()
        all_keys: Dict[Tuple[Hashable, ...], None] = {}  # Use dict as a sorted set for keys
        cpy = self.copy()
        for job in cpy._jobs.values():
            settings = cpy._get_job_settings(job)
            blocks = set(settings.block_keys(flatten_list=flatten_list))
            keys = {k: None for k in settings.nested_keys(flatten_list=flatten_list)}
            all_blocks = all_blocks.union(blocks)
            all_keys.update(keys)

        predicate = predicate if predicate else lambda _: True
        for key in all_keys:
            # Take only final nested keys i.e. those which are not block keys and satisfy the predicate
            if key not in all_blocks and predicate(key):
                field_key = "".join([str(k).title() for k in key])
                field = cpy._Field(
                    key=field_key, value_extractor=lambda j, k=key: cpy._get_job_settings(j).get_nested(k), from_settings=True  # type: ignore
                )
                if field_key not in cpy._fields:
                    cpy._fields[field_key] = field
        return cpy

    @staticmethod
    def _read_rkf(job: Job, section: str, variable: str, file: str = "ams"):
        if isinstance(job, AMSJob) and job.results is not None:
            return job.results.readrkf(section=section, variable=variable, file=file)
        else:
            return None

    @staticmethod
    def get_rkf_field_key(
        section: str,
        variable: str,
        file: str = "ams",
    ) -> str:
        """
        Gets the key for a |JobAnalysis| rkf field corresponding to a given file, section and varible.
        This will be a Pascal-case string of the file, section and varible e.g. ``AmsGeneralTask``.

        .. code:: python

            >>> JobAnalysis.get_rkf_field_key("General", "task")

            AmsGeneralTask

        :param section: section of the rkf file to read
        :param variable: variable in the section of the rkf file to read
        :param file: rkf file to use (without the .rkf extension)
        :return: corresponding key used for the |JobAnalysis| field
        """
        return JobAnalysis.get_settings_field_key((file, section, variable))

    def add_rkf_field(
        self,
        section: str,
        variable: str,
        file: str = "ams",
        display_name: Optional[str] = None,
        fmt: Optional[str] = None,
        expansion_depth: int = 0,
    ) -> "JobAnalysis":
        """
        Add a field for given variable in a section of a rkf file.
        The key of the fields will be a Pascal-case string of the file, section and varible e.g. ``AmsGeneralTask``

        .. code:: python

            >>> ja.add_rkf_field("General", "task")

            | Name  | ams:General%task |
            |-------|------------------|
            | job_1 | singlepoint      |
            | job_2 | singlepoint      |

        :param section: section of the rkf file to read
        :param variable: variable in the section of the rkf file to read
        :param file: rkf file to use (without the .rkf extension)
        :param display_name: name which will appear for the field when displayed in table
        :param fmt: string format for how field values are displayed in table
        :param expansion_depth: whether to expand field of multiple values into multiple rows, and recursively to what depth
        :return: copy of |JobAnalysis| with the additional field
        """
        key = JobAnalysis.get_rkf_field_key(section, variable, file)
        return self.add_field(
            key,
            lambda j: JobAnalysis._read_rkf(job=j, section=section, variable=variable, file=file),
            display_name=display_name,
            fmt=fmt,
            expansion_depth=expansion_depth,
        )

    def _get_job_settings(self, job: Job) -> Settings:
        """
        Get job settings converting any PISA input block to a standard settings object.
        """
        # Convert any PISA settings blocks to standard
        settings = Settings()
        if job.settings is not None:
            if isinstance(job.settings, Settings):
                settings = job.settings.copy()
            if _has_scm_pisa:
                if hasattr(job.settings, "input") and isinstance(job.settings.input, DriverBlock):
                    # Note use own input parser facade here to use caching
                    program = self._pisa_programs[job.settings.input.name].name.split(".")[0]
                    settings.input = input_to_settings(job.settings.input.get_input_string(), program=program)
        return settings

    def add_settings_input_fields(self, include_system_block: bool = False, flatten_list: bool = True) -> "JobAnalysis":
        """
        Add a field for each input key in the :attr:`~scm.plams.core.basejob.Job.settings` object across all currently added jobs.

        .. code:: python

            >>> ja.add_settings_input_fields()

            | Name  | InputAdfBasisType | InputAdfXcDispersion | InputAdfXcGga | InputAmsTask |
            |-------|-------------------|----------------------|---------------|--------------|
            | job_1 | TZP               | Grimme3              | PBE           | SinglePoint  |
            | job_2 | TZP               | Grimme3              | PBE           | SinglePoint  |

        :param include_system_block: whether to include keys for the system block, defaults to ``False``
        :param flatten_list: whether to flatten lists in settings objects
        :return: copy of |JobAnalysis| with the settings input fields added
        """

        def predicate(key_tuple: Tuple[Hashable, ...]):
            if len(key_tuple) == 0 or str(key_tuple[0]).lower() != "input":
                return False

            return (
                len(key_tuple) < 3
                or str(key_tuple[1]).lower() != "ams"
                or str(key_tuple[2]).lower() != "system"
                or include_system_block
            )

        return self.add_settings_fields(predicate, flatten_list)

    def remove_settings_fields(self) -> "JobAnalysis":
        """
        Remove all fields which were added as settings fields.

        .. code:: python

            >>> ja.add_settings_input_fields().remove_settings_fields()

            | Name  |
            |-------|
            | job_1 |
            | job_2 |

        :return: copy of |JobAnalysis| with settings input fields removed
        """
        cpy = self.copy()
        keys = [k for k, f in cpy._fields.items() if f.from_settings]
        for k in keys:
            cpy = cpy.remove_field(k)
        return cpy

    def __str__(self) -> str:
        """
        Get string representation of analysis as Markdown table with a maximum of 5 rows and column width of 12.

        .. code:: python

            >>> str(ja)

            | Name  | OK   |
            |-------|------|
            | job_1 | True |
            | job_2 | True |

        :return: markdown table of analysis
        """
        return self.to_table(max_col_width=12, max_rows=5)

    def __repr__(self) -> str:
        """
        Get string representation of analysis as Markdown table with a maximum of 5 rows and column width of 12.

        .. code:: python

            >>> ja

            | Name  | OK   |
            |-------|------|
            | job_1 | True |
            | job_2 | True |

        :return: markdown table of analysis
        """
        return self.to_table()

    def _repr_html_(self) -> str:
        return self.to_table(fmt="html")

    def __getitem__(self, key: str) -> List[Any]:
        """
        Get analysis data for a given field.

        .. code:: python

            >>> ja["Name"]

            ['job_1', 'job_2']

        :param key: unique identifier for the field
        :return: list of values for each job
        """
        # Get single analysis field unless fields are expanded, then the whole analysis must be calculated
        if any([f.expansion_depth > 0 for f in self._fields.values()]):
            if key not in self._fields:
                raise KeyError(f"Field with key '{key}' is not part of the analysis.")

            return self.get_analysis()[key]
        else:
            return self._get_field_analysis(key)

    def __setitem__(self, key: str, value: Callable[[Job], Any]) -> None:
        """
        Set analysis for given field.

        .. code:: python

            >>> ja["N"] = lambda j: len(j.molecule)
            >>> ja["N"]

            [4, 6]

        :param key: unique identifier for the field
        :param value: callable to extract the value for the field from a job
        """
        if not callable(value):
            raise TypeError("To set a field, the value must be a callable which accepts a Job.")

        if key in self._fields:
            self._fields[key] = replace(self._fields[key], value_extractor=value)
        else:
            self._add_field(key=key, value_extractor=value)

    def __delitem__(self, key: str) -> None:
        """
        Delete analysis for given field.

        .. code:: python

            >>> del ja["OK"]
            >>> ja

            | Name  |
            |-------|
            | job_1 |
            | job_2 |

        :param key: unique identifier for the field
        """
        self._remove_field(key)

    def __getattr__(self, key: str) -> List[Any]:
        """
        Fallback to get analysis for given field when an attribute is not present.

        .. code:: python

            >>> ja.Name

            ['job_1', 'job_2']

        :param key: unique identifier for the field
        :return: list of values for each job
        """
        try:
            return self[key]
        except KeyError:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute or analysis field with key '{key}'"
            )

    def __setattr__(self, key, value) -> None:
        """
        Fallback to set analysis for given field.

        .. code:: python

            >>> ja.N = lambda j: len(j.molecule)
            >>> ja.N

            [4, 6]

        :param key: unique identifier for the field
        :param value: callable to extract the value for the field from a job
        """
        if key in self._reserved_names or hasattr(self.__class__, key):
            super().__setattr__(key, value)
        else:
            self[key] = value

    def __delattr__(self, key) -> None:
        """
        Fallback to set analysis for given field.

        .. code:: python

            >>> del ja.OK
            >>> ja

            | Name  |
            |-------|
            | job_1 |
            | job_2 |

        :param key: unique identifier for the field
        """
        if key in self._reserved_names or hasattr(self.__class__, key):
            super().__delattr__(key)
        else:
            try:
                del self[key]
            except KeyError:
                raise AttributeError(
                    f"'{self.__class__.__name__}' object has no attribute or analysis field with key '{key}'"
                )

    def __contains__(self, key) -> bool:
        """
        Check whether a field is part of the analysis

        .. code:: python

            >>> "Name" in ja

            True

        :param key: unique identifier for the field
        :return: boolean flag for presence of field in the analysis
        """
        return key in self.field_keys

    def __dir__(self):
        """
        Return standard attributes, plus dynamically added field keys which can be accessed via dot notation.
        """
        return [x for x in super().__dir__()] + [
            k for k in self._fields.keys() if isinstance(k, str) and k.isidentifier()
        ]

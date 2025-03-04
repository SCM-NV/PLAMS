import os
import numpy as np
import pytest
import shutil
from datetime import datetime, timedelta
from pathlib import Path

from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.core.jobmanager import JobManager
from scm.plams.core.enums import JobStatus
from scm.plams.unit_tests.test_basejob import DummySingleJob
from scm.plams.unit_tests.test_helpers import temp_file_path, skip_if_no_scm_pisa, skip_if_no_scm_libbase
from scm.plams.tools.job_analysis import JobAnalysis
from scm.plams.core.settings import Settings, JobManagerSettings


class TestJobAnalysis:

    @pytest.fixture(scope="class")
    def dummy_single_jobs(self):
        # Generate dummy jobs for a selection of molecules and input settings
        smiles = ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", "Sys"]
        jobs = []
        for i, s in enumerate(smiles):
            sett = Settings()
            sett.input.ams.task = "GeometryOptimization" if i % 2 else "SinglePoint"
            sett.input.ams.Properties.NormalModes = "True" if i % 3 else "False"
            if i < 5:
                sett.input.ADF.Basis.Type = "TZP"
                if i % 2:
                    sett.input.ADF.xc.gga = "pbe"
            else:
                sett.input.DFTB
            if s == "Sys":
                sett.input.ams.System.Atoms = [
                    "Ar 0.0000000000       0.0000000000       0.0000000000",
                    "Ar 1.6050000000       0.9266471820       2.6050000000",
                ]
                mol = None
            else:
                mol = from_smiles(s)
            jobs.append(DummySingleJob(wait=i / 100, molecule=mol, settings=sett, name="dummyjob"))

        for j in jobs:
            j.run()
            j.ok()

        yield jobs

    def test_init_with_jobs(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        assert len(ja.jobs) == 10

    def test_init_with_paths(self, dummy_single_jobs):
        ja = JobAnalysis(paths=[j.path for j in dummy_single_jobs])
        assert len(ja.jobs) == 10

    def test_init_with_standard_fields(self, dummy_single_jobs):
        ja = JobAnalysis(standard_fields=("Name", "Smiles"))
        assert len(ja.field_keys) == 2

    def test_standard_fields(self, dummy_single_jobs):
        ja = JobAnalysis(
            jobs=dummy_single_jobs,
            standard_fields=[
                "Name",
                "OK",
                "Check",
                "ErrorMsg",
                "Formula",
                "Smiles",
                "GyrationRadius",
                "CPUTime",
                "SysTime",
                "ElapsedTime",
            ],
        ).format_field(
            "GyrationRadius", ".0f"
        )  # For different rdkit versions on GitHub and amspython

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | GyrationRadius | CPUTime | SysTime | ElapsedTime |
|--------------|------|-------|----------|---------|--------|----------------|---------|---------|-------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | 1              | None    | None    | None        |
| dummyjob.002 | True | True  | None     | CH4     | C      | 1              | None    | None    | None        |
| dummyjob.003 | True | True  | None     | H2O     | O      | 0              | None    | None    | None        |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | 1              | None    | None    | None        |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | 1              | None    | None    | None        |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | 2              | None    | None    | None        |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | 2              | None    | None    | None        |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | 3              | None    | None    | None        |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | 2              | None    | None    | None        |
| dummyjob.010 | True | True  | None     | None    | None   | None           | None    | None    | None        |"""
        )

        (
            ja.remove_fields(["CPUTime", "SysTime", "ElapsedTime", "GyrationRadius"]).add_standard_fields(
                ["ParentPath", "ParentName"]
            )
        )

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | ParentPath | ParentName |
|--------------|------|-------|----------|---------|--------|------------|------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | None       | None       |
| dummyjob.002 | True | True  | None     | CH4     | C      | None       | None       |
| dummyjob.003 | True | True  | None     | H2O     | O      | None       | None       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | None       | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | None       | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None       | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None       | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None       | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None       | None       |
| dummyjob.010 | True | True  | None     | None    | None   | None       | None       |"""
        )

        ja.remove_fields(["ParentPath", "ParentName", "Smiles", "Formula", "ErrorMsg", "Check", "OK"])

        assert (
            ja.to_table()
            == """\
| Name         |
|--------------|
| dummyjob     |
| dummyjob.002 |
| dummyjob.003 |
| dummyjob.004 |
| dummyjob.005 |
| dummyjob.006 |
| dummyjob.007 |
| dummyjob.008 |
| dummyjob.009 |
| dummyjob.010 |"""
        )

    def test_add_set_remove_rename_filter_fields(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_fields(["Formula", "Smiles", "CPUTime", "SysTime", "ElapsedTime"])
            .add_field("Wait", lambda j: j.wait)
            .set_field("Output", lambda _: None, display_name="foo")
            .set_field("Output", lambda j: j.results.read_file("$JN.out")[:5])
            .add_field(
                "Empty",
                lambda j: (
                    [] if j.wait == 0.01 else () if j.wait == 0.02 else np.array(None) if j.wait == 0.03 else None
                ),
            )
            .remove_field("Path")
            .rename_field("ErrorMsg", "Err")
        )

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | Err  | Formula | Smiles | CPUTime | SysTime | ElapsedTime | Wait | Output | Empty |
|--------------|------|-------|------|---------|--------|---------|---------|-------------|------|--------|-------|
| dummyjob     | True | True  | None | C2H6    | CC     | None    | None    | None        | 0.0  | Dummy  | None  |
| dummyjob.002 | True | True  | None | CH4     | C      | None    | None    | None        | 0.01 | Dummy  | []    |
| dummyjob.003 | True | True  | None | H2O     | O      | None    | None    | None        | 0.02 | Dummy  | ()    |
| dummyjob.004 | True | True  | None | CH4O    | CO     | None    | None    | None        | 0.03 | Dummy  | None  |
| dummyjob.005 | True | True  | None | C3H8    | CCC    | None    | None    | None        | 0.04 | Dummy  | None  |
| dummyjob.006 | True | True  | None | C4H10   | CCCC   | None    | None    | None        | 0.05 | Dummy  | None  |
| dummyjob.007 | True | True  | None | C3H8O   | CCCO   | None    | None    | None        | 0.06 | Dummy  | None  |
| dummyjob.008 | True | True  | None | C6H14   | CCCCCC | None    | None    | None        | 0.07 | Dummy  | None  |
| dummyjob.009 | True | True  | None | C4H10O  | CCCOC  | None    | None    | None        | 0.08 | Dummy  | None  |
| dummyjob.010 | True | True  | None | None    | None   | None    | None    | None        | 0.09 | Dummy  | None  |"""
        )

        ja.remove_empty_fields()

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | Formula | Smiles | Wait | Output |
|--------------|------|-------|---------|--------|------|--------|
| dummyjob     | True | True  | C2H6    | CC     | 0.0  | Dummy  |
| dummyjob.002 | True | True  | CH4     | C      | 0.01 | Dummy  |
| dummyjob.003 | True | True  | H2O     | O      | 0.02 | Dummy  |
| dummyjob.004 | True | True  | CH4O    | CO     | 0.03 | Dummy  |
| dummyjob.005 | True | True  | C3H8    | CCC    | 0.04 | Dummy  |
| dummyjob.006 | True | True  | C4H10   | CCCC   | 0.05 | Dummy  |
| dummyjob.007 | True | True  | C3H8O   | CCCO   | 0.06 | Dummy  |
| dummyjob.008 | True | True  | C6H14   | CCCCCC | 0.07 | Dummy  |
| dummyjob.009 | True | True  | C4H10O  | CCCOC  | 0.08 | Dummy  |
| dummyjob.010 | True | True  | None    | None   | 0.09 | Dummy  |"""
        )

        ja.remove_uniform_fields()

        assert (
            ja.to_table()
            == """\
| Name         | Formula | Smiles | Wait |
|--------------|---------|--------|------|
| dummyjob     | C2H6    | CC     | 0.0  |
| dummyjob.002 | CH4     | C      | 0.01 |
| dummyjob.003 | H2O     | O      | 0.02 |
| dummyjob.004 | CH4O    | CO     | 0.03 |
| dummyjob.005 | C3H8    | CCC    | 0.04 |
| dummyjob.006 | C4H10   | CCCC   | 0.05 |
| dummyjob.007 | C3H8O   | CCCO   | 0.06 |
| dummyjob.008 | C6H14   | CCCCCC | 0.07 |
| dummyjob.009 | C4H10O  | CCCOC  | 0.08 |
| dummyjob.010 | None    | None   | 0.09 |"""
        )

        ja.remove_uniform_fields(tol=0.1, ignore_empty=True).filter_fields(
            lambda vals: all([not v or "H" in v for v in vals])
        )

        assert (
            ja.to_table()
            == """\
| Formula |
|---------|
| C2H6    |
| CH4     |
| H2O     |
| CH4O    |
| C3H8    |
| C4H10   |
| C3H8O   |
| C6H14   |
| C4H10O  |
| None    |"""
        )

        ja.add_field(
            "MultiValue",
            lambda _: [[np.array([i if i != 1 else None for i in range(3)]) for _ in range(4)] for _ in range(5)],
        )
        ja.add_field(
            "Dict",
            lambda j: (
                {"a": (1, 2), "b": None, "c": {"d": [1, 2, 3], "e": None}}
                if j.wait < 0.05
                else {"b": "f", "c": {"e": "g", "d": [1, 2, 3]}, "a": (1, 2)}
            ),
        )
        ja.remove_uniform_fields(ignore_empty=True)

        assert (
            ja.to_table()
            == """\
| Formula |
|---------|
| C2H6    |
| CH4     |
| H2O     |
| CH4O    |
| C3H8    |
| C4H10   |
| C3H8O   |
| C6H14   |
| C4H10O  |
| None    |"""
        )

    def test_format_field(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .add_field("WaitS", lambda j: j.wait, fmt="e")
            .add_field("WaitMs", lambda j: j.wait * 1000)
            .format_field("WaitMs", "04.0f")
            .add_field("Output", lambda j: j.results.read_file("$JN.out")[:5])
            .remove_field("Path")
            .rename_field("ErrorMsg", "Err")
        )

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | Err  | Formula | Smiles | WaitS        | WaitMs | Output |
|--------------|------|-------|------|---------|--------|--------------|--------|--------|
| dummyjob     | True | True  | None | C2H6    | CC     | 0.000000e+00 | 0000   | Dummy  |
| dummyjob.002 | True | True  | None | CH4     | C      | 1.000000e-02 | 0010   | Dummy  |
| dummyjob.003 | True | True  | None | H2O     | O      | 2.000000e-02 | 0020   | Dummy  |
| dummyjob.004 | True | True  | None | CH4O    | CO     | 3.000000e-02 | 0030   | Dummy  |
| dummyjob.005 | True | True  | None | C3H8    | CCC    | 4.000000e-02 | 0040   | Dummy  |
| dummyjob.006 | True | True  | None | C4H10   | CCCC   | 5.000000e-02 | 0050   | Dummy  |
| dummyjob.007 | True | True  | None | C3H8O   | CCCO   | 6.000000e-02 | 0060   | Dummy  |
| dummyjob.008 | True | True  | None | C6H14   | CCCCCC | 7.000000e-02 | 0070   | Dummy  |
| dummyjob.009 | True | True  | None | C4H10O  | CCCOC  | 8.000000e-02 | 0080   | Dummy  |
| dummyjob.010 | True | True  | None | None    | None   | 9.000000e-02 | 0090   | Dummy  |"""
        )

    def test_expand_collapse_fields(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .remove_field("Path")
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .add_field("Atoms", lambda j: [at.symbol for at in j.molecule], expansion_depth=1)
            .add_field("MultiValue", lambda j: [[[f"w{i}" for i in range(int(j.wait * 100) // 2)]]])
        )

        assert (
            ja.to_table(max_rows=1000)
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Atoms                                    | MultiValue                   |
|--------------|------|-------|----------|---------|--------|------------------------------------------|------------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | C                                        | [[[]]]                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | C                                        | [[[]]]                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | [[[]]]                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | [[[]]]                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | [[[]]]                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | [[[]]]                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | [[[]]]                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | [[[]]]                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | C                                        | [[[]]]                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | [[[]]]                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | [[[]]]                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | [[[]]]                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | [[[]]]                       |
| dummyjob.003 | True | True  | None     | H2O     | O      | O                                        | [[['w0']]]                   |
| dummyjob.003 | True | True  | None     | H2O     | O      | H                                        | [[['w0']]]                   |
| dummyjob.003 | True | True  | None     | H2O     | O      | H                                        | [[['w0']]]                   |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | C                                        | [[['w0']]]                   |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | O                                        | [[['w0']]]                   |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | [[['w0']]]                   |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | [[['w0']]]                   |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | [[['w0']]]                   |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | [[['w0']]]                   |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | [[['w0', 'w1']]]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | [[['w0', 'w1']]]             |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | O                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | [[['w0', 'w1', 'w2']]]       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | O                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable | [[['w0', 'w1', 'w2', 'w3']]] |"""
        )

        ja.expand_field("MultiValue")

        assert (
            ja.to_table(max_rows=1000)
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Atoms                                    | MultiValue                 |
|--------------|------|-------|----------|---------|--------|------------------------------------------|----------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | C                                        | [[]]                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | C                                        | None                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | C                                        | [[]]                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None                       |
| dummyjob.003 | True | True  | None     | H2O     | O      | O                                        | [['w0']]                   |
| dummyjob.003 | True | True  | None     | H2O     | O      | H                                        | None                       |
| dummyjob.003 | True | True  | None     | H2O     | O      | H                                        | None                       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | C                                        | [['w0']]                   |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | O                                        | None                       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None                       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None                       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None                       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | [['w0', 'w1']]             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | [['w0', 'w1']]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | [['w0', 'w1', 'w2']]       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | O                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | [['w0', 'w1', 'w2']]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | [['w0', 'w1', 'w2', 'w3']] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | O                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                       |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable | [['w0', 'w1', 'w2', 'w3']] |"""
        )

        ja.expand_field("MultiValue", depth=2)

        assert (
            ja.to_table(max_rows=1000)
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Atoms                                    | MultiValue               |
|--------------|------|-------|----------|---------|--------|------------------------------------------|--------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | C                                        | []                       |
| dummyjob     | True | True  | None     | C2H6    | CC     | C                                        | None                     |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                     |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                     |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                     |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                     |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                     |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None                     |
| dummyjob.002 | True | True  | None     | CH4     | C      | C                                        | []                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None                     |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None                     |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None                     |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None                     |
| dummyjob.003 | True | True  | None     | H2O     | O      | O                                        | ['w0']                   |
| dummyjob.003 | True | True  | None     | H2O     | O      | H                                        | None                     |
| dummyjob.003 | True | True  | None     | H2O     | O      | H                                        | None                     |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | C                                        | ['w0']                   |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | O                                        | None                     |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None                     |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None                     |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None                     |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | ['w0', 'w1']             |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                     |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | ['w0', 'w1']             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | ['w0', 'w1', 'w2']       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | O                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                     |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | ['w0', 'w1', 'w2']       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | ['w0', 'w1', 'w2', 'w3'] |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | O                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None                     |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable | ['w0', 'w1', 'w2', 'w3'] |"""
        )

        ja.expand_field("MultiValue", depth=100)

        assert (
            ja.to_table(max_rows=1000)
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Atoms                                    | MultiValue |
|--------------|------|-------|----------|---------|--------|------------------------------------------|------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | C                                        | None       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None       |
| dummyjob     | True | True  | None     | C2H6    | CC     | H                                        | None       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None       |
| dummyjob.002 | True | True  | None     | CH4     | C      | H                                        | None       |
| dummyjob.003 | True | True  | None     | H2O     | O      | O                                        | w0         |
| dummyjob.003 | True | True  | None     | H2O     | O      | H                                        | None       |
| dummyjob.003 | True | True  | None     | H2O     | O      | H                                        | None       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | C                                        | w0         |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | O                                        | None       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | H                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | w0         |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | w1         |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | C                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | w0         |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | w1         |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | C                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | H                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | w0         |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | w1         |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | w2         |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | C                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | O                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | w0         |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | w1         |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | w2         |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | C                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | w0         |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | w1         |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | w2         |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | w3         |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | O                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | C                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | H                                        | None       |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable | w0         |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable | w1         |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable | w2         |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable | w3         |"""
        )

        ja.collapse_field("Atoms")

        assert (
            ja.to_table(max_rows=1000)
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Atoms                                                                                                | MultiValue |
|--------------|------|-------|----------|---------|--------|------------------------------------------------------------------------------------------------------|------------|
| dummyjob.003 | True | True  | None     | H2O     | O      | ['O', 'H', 'H']                                                                                      | w0         |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | ['C', 'O', 'H', 'H', 'H', 'H']                                                                       | w0         |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                                              | w0         |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                                              | w1         |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | ['C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                               | w0         |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | ['C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                               | w1         |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | ['C', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                                         | w0         |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | ['C', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                                         | w1         |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | ['C', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                                         | w2         |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'] | w0         |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'] | w1         |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'] | w2         |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | ['C', 'C', 'C', 'O', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                          | w0         |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | ['C', 'C', 'C', 'O', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                          | w1         |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | ['C', 'C', 'C', 'O', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                          | w2         |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | ['C', 'C', 'C', 'O', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                          | w3         |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable                                                             | w0         |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable                                                             | w1         |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable                                                             | w2         |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable                                                             | w3         |"""
        )

        ja.collapse_field("MultiValue")

        assert (
            ja.to_table(max_rows=1000)
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Atoms                                                                                                | MultiValue                   |
|--------------|------|-------|----------|---------|--------|------------------------------------------------------------------------------------------------------|------------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']                                                             | [[[]]]                       |
| dummyjob.002 | True | True  | None     | CH4     | C      | ['C', 'H', 'H', 'H', 'H']                                                                            | [[[]]]                       |
| dummyjob.003 | True | True  | None     | H2O     | O      | ['O', 'H', 'H']                                                                                      | [[['w0']]]                   |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | ['C', 'O', 'H', 'H', 'H', 'H']                                                                       | [[['w0']]]                   |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                                              | [[['w0', 'w1']]]             |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | ['C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                               | [[['w0', 'w1']]]             |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | ['C', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                                         | [[['w0', 'w1', 'w2']]]       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'] | [[['w0', 'w1', 'w2']]]       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | ['C', 'C', 'C', 'O', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                          | [[['w0', 'w1', 'w2', 'w3']]] |
| dummyjob.010 | True | True  | None     | None    | None   | ERROR: 'NoneType' object is not iterable                                                             | [[['w0', 'w1', 'w2', 'w3']]] |"""
        )

    def test_add_job(self, dummy_single_jobs):
        ja = JobAnalysis()
        for j in dummy_single_jobs:
            ja.add_job(j)
        assert len(ja.jobs) == 10

    @pytest.mark.parametrize("use_custom_loader", [False, True])
    def test_load_job(self, use_custom_loader, dummy_single_jobs):
        ja = JobAnalysis(paths=[j.path for j in dummy_single_jobs])

        job = DummySingleJob(molecule=from_smiles("N"), name="dummyloadjob")
        job.run().wait()
        path = Path(job.path)
        dill_file = path / f"{path.name}.dill"
        os.remove(dill_file)

        if use_custom_loader:
            calls = []

            def loader(p):
                calls.append(p)
                return job

            ja.load_job(job.path, loaders=[loader])
            assert calls == [job.path]
        else:
            ja.load_job(job.path)

        assert len(ja.jobs) == 11

    def test_filter_jobs(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .add_field("Wait", lambda j: j.wait)
            .add_field("Output", lambda j: j.results.read_file("$JN.out")[:5])
            .remove_field("Path")
            .filter_jobs(lambda d: d["Formula"] is not None and "O" not in d["Smiles"])
        )

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Wait | Output |
|--------------|------|-------|----------|---------|--------|------|--------|
| dummyjob     | True | True  | None     | C2H6    | CC     | 0.0  | Dummy  |
| dummyjob.002 | True | True  | None     | CH4     | C      | 0.01 | Dummy  |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | 0.04 | Dummy  |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | 0.05 | Dummy  |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | 0.07 | Dummy  |"""
        )

        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .add_field("Wait", lambda j: j.wait)
            .add_field("Output", lambda j: j.results.read_file("$JN.out")[:5])
            .add_field("Atoms", lambda j: [at.symbol for at in j.molecule], expansion_depth=1)
            .remove_field("Path")
            .filter_jobs(lambda d: d["Formula"] is not None and "O" not in d["Smiles"])
            .collapse_field("Atoms")
        )

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Wait | Output | Atoms                                                                                                |
|--------------|------|-------|----------|---------|--------|------|--------|------------------------------------------------------------------------------------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | 0.0  | Dummy  | ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']                                                             |
| dummyjob.002 | True | True  | None     | CH4     | C      | 0.01 | Dummy  | ['C', 'H', 'H', 'H', 'H']                                                                            |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | 0.04 | Dummy  | ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                                              |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | 0.05 | Dummy  | ['C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                               |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | 0.07 | Dummy  | ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'] |"""
        )

        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .add_field("Wait", lambda j: j.wait)
            .add_field("Output", lambda j: j.results.read_file("$JN.out")[:5])
            .add_field("Atoms", lambda j: [at.symbol for at in j.molecule], expansion_depth=1)
            .remove_field("Path")
            .filter_jobs(lambda d: d["Formula"] is not None and d["Atoms"] != "O")
            .collapse_field("Atoms")
        )

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Wait | Output | Atoms                                                                                                |
|--------------|------|-------|----------|---------|--------|------|--------|------------------------------------------------------------------------------------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | 0.0  | Dummy  | ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']                                                             |
| dummyjob.002 | True | True  | None     | CH4     | C      | 0.01 | Dummy  | ['C', 'H', 'H', 'H', 'H']                                                                            |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | 0.04 | Dummy  | ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                                              |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | 0.05 | Dummy  | ['C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']                               |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | 0.07 | Dummy  | ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'] |"""
        )

    def test_reorder_fields(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .add_standard_field("CPUTime")
            .add_standard_field("SysTime")
            .add_standard_field("ElapsedTime")
            .add_field("Wait", lambda j: j.wait)
        )

        assert ja.field_keys == [
            "Path",
            "Name",
            "OK",
            "Check",
            "ErrorMsg",
            "Formula",
            "Smiles",
            "CPUTime",
            "SysTime",
            "ElapsedTime",
            "Wait",
        ]

        ja.reorder_fields(["Name", "Wait"])

        assert ja.field_keys == [
            "Name",
            "Wait",
            "Path",
            "OK",
            "Check",
            "ErrorMsg",
            "Formula",
            "Smiles",
            "CPUTime",
            "SysTime",
            "ElapsedTime",
        ]

        ja.sort_fields(lambda k: str(k))

        assert ja.field_keys == [
            "CPUTime",
            "Check",
            "ElapsedTime",
            "ErrorMsg",
            "Formula",
            "Name",
            "OK",
            "Path",
            "Smiles",
            "SysTime",
            "Wait",
        ]

        ja.sort_fields(lambda k: len(k), reverse=True)

        assert ja.field_keys == [
            "ElapsedTime",
            "ErrorMsg",
            "CPUTime",
            "Formula",
            "SysTime",
            "Smiles",
            "Check",
            "Name",
            "Path",
            "Wait",
            "OK",
        ]

    def test_contains_field(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)

        assert "Path" in ja
        assert "foo" not in ja

    def test_sort_jobs(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_field("Wait", lambda j: j.wait)
            .sort_jobs(field_keys=["Formula"])
        )

        assert [j.name for j in ja.jobs.values()] == [
            "dummyjob",
            "dummyjob.005",
            "dummyjob.007",
            "dummyjob.006",
            "dummyjob.009",
            "dummyjob.008",
            "dummyjob.002",
            "dummyjob.004",
            "dummyjob.003",
            "dummyjob.010",
        ]

        ja.sort_jobs(field_keys=["OK", "Wait"], reverse=True)

        assert [(j.name, j.wait) for j in ja.jobs.values()] == [
            ("dummyjob.010", 0.09),
            ("dummyjob.009", 0.08),
            ("dummyjob.008", 0.07),
            ("dummyjob.007", 0.06),
            ("dummyjob.006", 0.05),
            ("dummyjob.005", 0.04),
            ("dummyjob.004", 0.03),
            ("dummyjob.003", 0.02),
            ("dummyjob.002", 0.01),
            ("dummyjob", 0.0),
        ]

        ja.sort_jobs(sort_key=lambda data: data["Wait"] * 100 % 5)

        assert [(j.name, int(j.wait * 100 % 5)) for j in ja.jobs.values()] == [
            ("dummyjob.006", 0.0),
            ("dummyjob", 0.0),
            ("dummyjob.007", 1.0),
            ("dummyjob.002", 1.0),
            ("dummyjob.003", 2.0),
            ("dummyjob.008", 2.0),
            ("dummyjob.009", 3.0),
            ("dummyjob.004", 3.0),
            ("dummyjob.010", 4.0),
            ("dummyjob.005", 4.0),
        ]

    def test_settings_fields(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .add_settings_input_fields()
            .add_settings_field(("runscript", "shebang"))
            .remove_field("Path")
        )

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | InputAdfBasisType | InputAmsTask         | InputAmsPropertiesNormalmodes | InputAdfXcGga | RunscriptShebang |
|--------------|------|-------|----------|---------|--------|-------------------|----------------------|-------------------------------|---------------|------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | TZP               | SinglePoint          | False                         | None          | #!/bin/sh        |
| dummyjob.002 | True | True  | None     | CH4     | C      | TZP               | GeometryOptimization | True                          | pbe           | #!/bin/sh        |
| dummyjob.003 | True | True  | None     | H2O     | O      | TZP               | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | TZP               | GeometryOptimization | False                         | pbe           | #!/bin/sh        |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | TZP               | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None              | GeometryOptimization | True                          | None          | #!/bin/sh        |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None              | SinglePoint          | False                         | None          | #!/bin/sh        |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None              | GeometryOptimization | True                          | None          | #!/bin/sh        |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None              | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.010 | True | True  | None     | None    | None   | None              | GeometryOptimization | False                         | None          | #!/bin/sh        |"""
        )

        ja.remove_settings_fields().add_settings_input_fields(include_system_block=True)

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | InputAdfBasisType | InputAmsTask         | InputAmsPropertiesNormalmodes | InputAdfXcGga | InputAmsSystemAtoms0                                  | InputAmsSystemAtoms1                                  |
|--------------|------|-------|----------|---------|--------|-------------------|----------------------|-------------------------------|---------------|-------------------------------------------------------|-------------------------------------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | TZP               | SinglePoint          | False                         | None          | None                                                  | None                                                  |
| dummyjob.002 | True | True  | None     | CH4     | C      | TZP               | GeometryOptimization | True                          | pbe           | None                                                  | None                                                  |
| dummyjob.003 | True | True  | None     | H2O     | O      | TZP               | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | TZP               | GeometryOptimization | False                         | pbe           | None                                                  | None                                                  |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | TZP               | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None              | GeometryOptimization | True                          | None          | None                                                  | None                                                  |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None              | SinglePoint          | False                         | None          | None                                                  | None                                                  |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None              | GeometryOptimization | True                          | None          | None                                                  | None                                                  |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None              | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.010 | True | True  | None     | None    | None   | None              | GeometryOptimization | False                         | None          | Ar 0.0000000000       0.0000000000       0.0000000000 | Ar 1.6050000000       0.9266471820       2.6050000000 |"""
        )

    def test_get_set_del_item(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs).add_standard_field("Formula").add_standard_field("Smiles")
        del ja["Path"]
        del ja["Check"]
        ja["OK"] = lambda j: "Yes" if j.ok() else "No"
        ja["Id"] = lambda j: j.name.split(".")[-1]

        assert ja["Name"] == [
            "dummyjob",
            "dummyjob.002",
            "dummyjob.003",
            "dummyjob.004",
            "dummyjob.005",
            "dummyjob.006",
            "dummyjob.007",
            "dummyjob.008",
            "dummyjob.009",
            "dummyjob.010",
        ]
        assert ja["Smiles"] == ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", None]
        assert ja["OK"] == ["Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"]
        assert ja["Id"] == ["dummyjob", "002", "003", "004", "005", "006", "007", "008", "009", "010"]
        with pytest.raises(KeyError):
            ja["Foo"]
        with pytest.raises(KeyError):
            del ja["Bar"]

        ja.filter_jobs(lambda data: data["Smiles"] in ("C", "CC", "CO")).add_field(
            "Atoms", lambda j: [at.symbol for at in j.molecule]
        )

        assert ja["Smiles"] == ["CC", "C", "CO"]
        assert ja["Atoms"] == [
            ["C", "C", "H", "H", "H", "H", "H", "H"],
            ["C", "H", "H", "H", "H"],
            ["C", "O", "H", "H", "H", "H"],
        ]

        ja.expand_field("Atoms")
        assert ja["Smiles"] == [
            "CC",
            "CC",
            "CC",
            "CC",
            "CC",
            "CC",
            "CC",
            "CC",
            "C",
            "C",
            "C",
            "C",
            "C",
            "CO",
            "CO",
            "CO",
            "CO",
            "CO",
            "CO",
        ]
        assert ja["Atoms"] == [
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "C",
            "H",
            "H",
            "H",
            "H",
            "C",
            "O",
            "H",
            "H",
            "H",
            "H",
        ]

    def test_get_set_del_attributes(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs).add_standard_field("Formula").add_standard_field("Smiles")
        del ja.Path
        del ja.Check

        ja.OK = lambda j: "Yes" if j.ok() else "No"
        ja.Id = lambda j: j.name.split(".")[-1]

        assert ja.Name == [
            "dummyjob",
            "dummyjob.002",
            "dummyjob.003",
            "dummyjob.004",
            "dummyjob.005",
            "dummyjob.006",
            "dummyjob.007",
            "dummyjob.008",
            "dummyjob.009",
            "dummyjob.010",
        ]
        assert ja.Smiles == ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", None]
        assert ja.OK == ["Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"]
        assert ja.Id == ["dummyjob", "002", "003", "004", "005", "006", "007", "008", "009", "010"]
        with pytest.raises(AttributeError):
            ja.Foo
        with pytest.raises(AttributeError):
            del ja.Bar

    def test_to_table(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .remove_field("Path")
            .add_settings_field(("Input", "AMS", "Properties", "NormalModes"))
            .remove_empty_fields()
            .remove_uniform_fields()
        )

        assert (
            ja.to_table(max_col_width=10, max_rows=6)
            == """\
| Name          | Formula | Smiles | InputAmsPr... |
|---------------|---------|--------|---------------|
| dummyjob      | C2H6    | CC     | False         |
| dummyjob.0... | CH4     | C      | True          |
| dummyjob.0... | H2O     | O      | True          |
| ...           | ...     | ...    | ...           |
| dummyjob.0... | C6H14   | CCCCCC | True          |
| dummyjob.0... | C4H10O  | CCCOC  | True          |
| dummyjob.0... | None    | None   | False         |"""
        )

    def test_to_csv(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .remove_field("Path")
            .add_settings_field(("Input", "AMS", "Properties", "NormalModes"))
            .remove_empty_fields()
            .remove_uniform_fields()
        )

        with temp_file_path(".csv") as tfp:
            ja.to_csv_file(tfp)
            with open(tfp) as tf:
                csv = tf.read()

        assert (
            csv
            == """\
Name,Formula,Smiles,InputAmsPropertiesNormalmodes
dummyjob,C2H6,CC,False
dummyjob.002,CH4,C,True
dummyjob.003,H2O,O,True
dummyjob.004,CH4O,CO,False
dummyjob.005,C3H8,CCC,True
dummyjob.006,C4H10,CCCC,True
dummyjob.007,C3H8O,CCCO,False
dummyjob.008,C6H14,CCCCCC,True
dummyjob.009,C4H10O,CCCOC,True
dummyjob.010,,,False
"""
        )

    def test_to_dataframe(self, dummy_single_jobs):
        try:
            import pandas  # noqa F401
        except ImportError:
            pytest.skip("Skipping test as cannot find pandas package.")

        df = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .remove_field("Path")
            .add_settings_field(("Input", "AMS", "Properties", "NormalModes"))
            .remove_empty_fields()
            .remove_uniform_fields()
            .to_dataframe()
        )

        assert df.shape == (10, 4)
        assert df.columns.to_list() == [
            "Name",
            "Formula",
            "Smiles",
            "InputAmsPropertiesNormalmodes",
        ]
        assert df.Formula.to_list() == [
            "C2H6",
            "CH4",
            "H2O",
            "CH4O",
            "C3H8",
            "C4H10",
            "C3H8O",
            "C6H14",
            "C4H10O",
            None,
        ]

    def test_get_timeline(self):
        jm = JobManager(JobManagerSettings())

        # Function to create dummy job with given statuses and timeline
        start = datetime.strptime("2025-01-01 12:00:00", "%Y-%m-%d %H:%M:%S")

        def get_job_with_statuses(name, final_status, start_time, wait_time, run_time):
            job = DummySingleJob(name=name)
            job.run(jobmanager=jm)
            statuses = [
                JobStatus.CREATED,
                JobStatus.STARTED,
                JobStatus.REGISTERED,
                JobStatus.RUNNING,
                JobStatus.FINISHED,
                final_status,
            ]
            times = [
                start_time,
                start_time + wait_time,
                start_time + wait_time + timedelta(seconds=0.05),
                start_time + wait_time + timedelta(seconds=0.1),
                start_time + wait_time + run_time,
                start_time + wait_time + run_time + timedelta(seconds=0.05),
            ]
            job._status_log = [(t, s) for t, s in zip(times, statuses)]
            return job

        # Given job analyser with no jobs
        ja = JobAnalysis()

        # When get timeline
        # Then empty table produced
        timeline = ja.get_timeline()
        assert (
            timeline
            == """\
<pre>
| JobName | WaitDuration | RunDuration | TotalDuration |
|---------|--------------|-------------|---------------|
</pre>"""
        )

        # Add single job that takes seconds to run
        job_s = get_job_with_statuses("tls", JobStatus.SUCCESSFUL, start, timedelta(seconds=0.1), timedelta(seconds=3))
        ja.add_job(job_s)
        timeline = ja.get_timeline()
        assert (
            timeline
            == """\
<pre>
| JobName | ↓2025-01-01 12:00:00 | ↓2025-01-01 12:00:01 | ↓2025-01-01 12:00:02 | ↓2025-01-01 12:00:03 | WaitDuration | RunDuration | TotalDuration |
|---------|----------------------|----------------------|----------------------|----------------------|--------------|-------------|---------------|
| tls     | .-+================= | ==================== | ===================* | >                    | 0s           | 3s          | 3s            |
</pre>"""
        )

        # Add single job that takes minutes to run
        job_m = get_job_with_statuses("tlm", JobStatus.FAILED, start, timedelta(seconds=2), timedelta(minutes=2))
        ja.add_job(job_m)
        timeline = ja.get_timeline()
        assert (
            timeline
            == """\
<pre>
| JobName | ↓2025-01-01 12:00:00 | ↓2025-01-01 12:00:30 | ↓2025-01-01 12:01:01 | ↓2025-01-01 12:01:31 | ↓2025-01-01 12:02:02 | WaitDuration | RunDuration | TotalDuration |
|---------|----------------------|----------------------|----------------------|----------------------|----------------------|--------------|-------------|---------------|
| tls     | ==>                  |                      |                      |                      |                      | 0s           | 3s          | 3s            |
| tlm     | .=================== | ==================== | ==================== | ===================* | X                    | 2s           | 2m0s        | 2m2s          |
</pre>"""
        )

        # Add single job that takes hours to run
        job_m = get_job_with_statuses("tlh", JobStatus.CRASHED, start, timedelta(seconds=10), timedelta(hours=2))
        ja.add_job(job_m)
        timeline = ja.get_timeline()
        assert (
            timeline
            == """\
<pre>
| JobName | ↓2025-01-01 12:00:00 | ↓2025-01-01 12:30:02 | ↓2025-01-01 13:00:05 | ↓2025-01-01 13:30:07 | ↓2025-01-01 14:00:10 | WaitDuration | RunDuration | TotalDuration |
|---------|----------------------|----------------------|----------------------|----------------------|----------------------|--------------|-------------|---------------|
| tls     | >                    |                      |                      |                      |                      | 0s           | 3s          | 3s            |
| tlm     | =X                   |                      |                      |                      |                      | 2s           | 2m0s        | 2m2s          |
| tlh     | ==================== | ==================== | ==================== | ===================* | x                    | 10s          | 2h0m0s      | 2h0m10s       |
</pre>"""
        )

        # Add more single jobs with delayed start
        job_ds = get_job_with_statuses(
            "tlds", JobStatus.SUCCESSFUL, start + timedelta(hours=1), timedelta(seconds=1), timedelta(seconds=30)
        )
        job_dm = get_job_with_statuses(
            "tldm", JobStatus.SUCCESSFUL, start + timedelta(hours=1), timedelta(minutes=2), timedelta(minutes=10)
        )

        ja.add_job(job_ds)
        ja.add_job(job_dm)
        timeline = ja.get_timeline(max_intervals=5)
        assert (
            timeline
            == """\
<pre>
| JobName | ↓2025-01-01 12:00:00 | ↓2025-01-01 12:30:02 | ↓2025-01-01 13:00:05 | ↓2025-01-01 13:30:07 | ↓2025-01-01 14:00:10 | WaitDuration | RunDuration | TotalDuration |
|---------|----------------------|----------------------|----------------------|----------------------|----------------------|--------------|-------------|---------------|
| tls     | >                    |                      |                      |                      |                      | 0s           | 3s          | 3s            |
| tlm     | =X                   |                      |                      |                      |                      | 2s           | 2m0s        | 2m2s          |
| tlh     | ==================== | ==================== | ==================== | ===================* | x                    | 10s          | 2h0m0s      | 2h0m10s       |
| tlds    |                      |                    = | >                    |                      |                      | 1s           | 30s         | 31s           |
| tldm    |                      |                    . | .======>             |                      |                      | 2m0s         | 10m0s       | 12m0s         |
</pre>"""
        )

        # Add more single jobs with dependencies
        job_dps = get_job_with_statuses(
            "tldps", JobStatus.SUCCESSFUL, start + timedelta(hours=1), timedelta(minutes=12), timedelta(seconds=10)
        )
        job_dpm = get_job_with_statuses(
            "tldpm", JobStatus.SUCCESSFUL, start + timedelta(hours=1), timedelta(minutes=12), timedelta(minutes=10)
        )

        ja.add_job(job_dps)
        ja.add_job(job_dpm)
        timeline = ja.get_timeline(max_intervals=10)
        assert (
            timeline
            == """\
<pre>
| JobName | ↓2025-01-01 12:00:00 | ↓2025-01-01 12:13:21 | ↓2025-01-01 12:26:42 | ↓2025-01-01 12:40:03 | ↓2025-01-01 12:53:24 | ↓2025-01-01 13:06:45 | ↓2025-01-01 13:20:06 | ↓2025-01-01 13:33:27 | ↓2025-01-01 13:46:48 | ↓2025-01-01 14:00:10 | WaitDuration | RunDuration | TotalDuration |
|---------|----------------------|----------------------|----------------------|----------------------|----------------------|----------------------|----------------------|----------------------|----------------------|----------------------|--------------|-------------|---------------|
| tls     | >                    |                      |                      |                      |                      |                      |                      |                      |                      |                      | 0s           | 3s          | 3s            |
| tlm     | ===X                 |                      |                      |                      |                      |                      |                      |                      |                      |                      | 2s           | 2m0s        | 2m2s          |
| tlh     | ==================== | ==================== | ==================== | ==================== | ==================== | ==================== | ==================== | ==================== | ===================* | x                    | 10s          | 2h0m0s      | 2h0m10s       |
| tlds    |                      |                      |                      |                      |          =>          |                      |                      |                      |                      |                      | 1s           | 30s         | 31s           |
| tldm    |                      |                      |                      |                      |          ...======== | =======>             |                      |                      |                      |                      | 2m0s         | 10m0s       | 12m0s         |
| tldps   |                      |                      |                      |                      |          ........... | .......=>            |                      |                      |                      |                      | 12m0s        | 10s         | 12m10s        |
| tldpm   |                      |                      |                      |                      |          ........... | .......============= | ==>                  |                      |                      |                      | 12m0s        | 10m0s       | 22m0s         |
</pre>"""
        )

        # Add more single jobs with dependencies
        job_n = DummySingleJob(name="tln")
        job_n.run(jobmanager=jm)
        job_n.results.wait()
        job_n._status_log = None

        job_e = DummySingleJob(name="tle")
        job_e.run(jobmanager=jm)
        job_e.results.wait()
        job_e._status_log = []

        ja.add_job(job_n)
        ja.add_job(job_e)
        timeline = ja.get_timeline()
        assert (
            timeline
            == """\
<pre>
| JobName | ↓2025-01-01 12:00:00 | ↓2025-01-01 12:30:02 | ↓2025-01-01 13:00:05 | ↓2025-01-01 13:30:07 | ↓2025-01-01 14:00:10 | WaitDuration | RunDuration | TotalDuration |
|---------|----------------------|----------------------|----------------------|----------------------|----------------------|--------------|-------------|---------------|
| tls     | >                    |                      |                      |                      |                      | 0s           | 3s          | 3s            |
| tlm     | =X                   |                      |                      |                      |                      | 2s           | 2m0s        | 2m2s          |
| tlh     | ==================== | ==================== | ==================== | ===================* | x                    | 10s          | 2h0m0s      | 2h0m10s       |
| tlds    |                      |                    = | >                    |                      |                      | 1s           | 30s         | 31s           |
| tldm    |                      |                    . | .======>             |                      |                      | 2m0s         | 10m0s       | 12m0s         |
| tldps   |                      |                    . | .......=>            |                      |                      | 12m0s        | 10s         | 12m10s        |
| tldpm   |                      |                    . | .......=======>      |                      |                      | 12m0s        | 10m0s       | 22m0s         |
| tln     |                      |                      |                      |                      |                      | Unknown      | Unknown     | Unknown       |
| tle     |                      |                      |                      |                      |                      | Unknown      | Unknown     | Unknown       |
</pre>"""
        )

        timeline = ja.get_timeline(fmt="html")
        assert (
            timeline
            == """\
<div style="max-width: 100%; overflow-x: auto;">
<table border="1" style="border-collapse: collapse; width: auto; font-family: monospace; ">
<thead><tr><th style="border-left: 1px solid black; border-right: 1px solid black;">JobName<th style="border-left: 1px solid black; border-right: 1px solid black;">↓2025-01-01 12:00:00                                                                                                    <th style="border-left: 1px solid black; border-right: 1px solid black;">↓2025-01-01 12:30:02                                                                                                    <th style="border-left: 1px solid black; border-right: 1px solid black;">↓2025-01-01 13:00:05                                                                                                    <th style="border-left: 1px solid black; border-right: 1px solid black;">↓2025-01-01 13:30:07                                                                                                    <th style="border-left: 1px solid black; border-right: 1px solid black;">↓2025-01-01 14:00:10                                                                                                    <th style="border-left: 1px solid black; border-right: 1px solid black;">WaitDuration<th style="border-left: 1px solid black; border-right: 1px solid black;">RunDuration<th style="border-left: 1px solid black; border-right: 1px solid black;">TotalDuration</th></tr></thead>
<tbody>
<tr><td style="border-left: 1px solid black; border-right: 1px solid black;">tls    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">0s          </td><td style="border-left: 1px solid black; border-right: 1px solid black;">3s         </td><td style="border-left: 1px solid black; border-right: 1px solid black;">3s           </td></tr>
<tr><td style="border-left: 1px solid black; border-right: 1px solid black;">tlm    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">=X&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">2s          </td><td style="border-left: 1px solid black; border-right: 1px solid black;">2m0s       </td><td style="border-left: 1px solid black; border-right: 1px solid black;">2m2s         </td></tr>
<tr><td style="border-left: 1px solid black; border-right: 1px solid black;">tlh    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">====================                                                                                                    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">====================                                                                                                    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">====================                                                                                                    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">===================*                                                                                                    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">x&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">10s         </td><td style="border-left: 1px solid black; border-right: 1px solid black;">2h0m0s     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">2h0m10s      </td></tr>
<tr><td style="border-left: 1px solid black; border-right: 1px solid black;">tlds   </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;=     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">1s          </td><td style="border-left: 1px solid black; border-right: 1px solid black;">30s        </td><td style="border-left: 1px solid black; border-right: 1px solid black;">31s          </td></tr>
<tr><td style="border-left: 1px solid black; border-right: 1px solid black;">tldm   </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">.======>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                                        </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">2m0s        </td><td style="border-left: 1px solid black; border-right: 1px solid black;">10m0s      </td><td style="border-left: 1px solid black; border-right: 1px solid black;">12m0s        </td></tr>
<tr><td style="border-left: 1px solid black; border-right: 1px solid black;">tldps  </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">.......=>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                                             </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">12m0s       </td><td style="border-left: 1px solid black; border-right: 1px solid black;">10s        </td><td style="border-left: 1px solid black; border-right: 1px solid black;">12m10s       </td></tr>
<tr><td style="border-left: 1px solid black; border-right: 1px solid black;">tldpm  </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">.......=======>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                                                                           </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">12m0s       </td><td style="border-left: 1px solid black; border-right: 1px solid black;">10m0s      </td><td style="border-left: 1px solid black; border-right: 1px solid black;">22m0s        </td></tr>
<tr><td style="border-left: 1px solid black; border-right: 1px solid black;">tln    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">Unknown     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">Unknown    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">Unknown      </td></tr>
<tr><td style="border-left: 1px solid black; border-right: 1px solid black;">tle    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td style="border-left: 1px solid black; border-right: 1px solid black;">Unknown     </td><td style="border-left: 1px solid black; border-right: 1px solid black;">Unknown    </td><td style="border-left: 1px solid black; border-right: 1px solid black;">Unknown      </td></tr>
</tbody>
</table>
</div>"""
        )

        jm._clean()
        shutil.rmtree(jm.workdir)


class TestJobAnalysisWithPisa(TestJobAnalysis):

    @pytest.fixture(scope="class")
    def dummy_single_jobs(self):
        skip_if_no_scm_pisa()

        from scm.input_classes.drivers import AMS
        from scm.input_classes.engines import DFTB, ADF

        # Generate dummy jobs for a selection of molecules and input settings
        smiles = ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", "Sys"]
        jobs = []
        for i, s in enumerate(smiles):
            sett = AMS()
            sett.Task = "GeometryOptimization" if i % 2 else "SinglePoint"
            sett.Properties.NormalModes = "Yes" if i % 3 else "No"
            if i < 5:
                sett.Engine = ADF()
                sett.Engine.Basis.Type = "TZP"
                if i % 2:
                    sett.Engine.XC.GGA = "pbe"
            else:
                sett.Engine = DFTB()
            if s == "Sys":
                sett.System.Atoms = [
                    "Ar 0.0000000000       0.0000000000       0.0000000000",
                    "Ar 1.6050000000       0.9266471820       2.6050000000",
                ]
                mol = None
            else:
                mol = from_smiles(s)

            if i > 5:
                sett = Settings({"input": sett})

            jobs.append(DummySingleJob(wait=i / 100, molecule=mol, settings=sett, name="dummyjob"))

        jm = JobManager(JobManagerSettings())
        for j in jobs:
            j.run(jobmanager=jm)
            j.ok()
        yield jobs

        jm._clean()
        shutil.rmtree(jm.workdir)

    def test_settings_fields(self, dummy_single_jobs):
        ja = (
            JobAnalysis(jobs=dummy_single_jobs)
            .add_standard_field("Formula")
            .add_standard_field("Smiles")
            .add_settings_input_fields()
            .add_settings_field(("runscript", "shebang"))
            .remove_field("Path")
        )

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | InputAdfBasisType | InputAmsTask         | InputAmsPropertiesNormalmodes | InputAdfXcGga | RunscriptShebang |
|--------------|------|-------|----------|---------|--------|-------------------|----------------------|-------------------------------|---------------|------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | TZP               | SinglePoint          | False                         | None          | #!/bin/sh        |
| dummyjob.002 | True | True  | None     | CH4     | C      | TZP               | GeometryOptimization | True                          | pbe           | #!/bin/sh        |
| dummyjob.003 | True | True  | None     | H2O     | O      | TZP               | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | TZP               | GeometryOptimization | False                         | pbe           | #!/bin/sh        |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | TZP               | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None              | GeometryOptimization | True                          | None          | #!/bin/sh        |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None              | SinglePoint          | False                         | None          | #!/bin/sh        |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None              | GeometryOptimization | True                          | None          | #!/bin/sh        |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None              | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.010 | True | True  | None     | None    | None   | None              | GeometryOptimization | False                         | None          | #!/bin/sh        |"""
        )

        ja.remove_settings_fields().add_settings_input_fields(include_system_block=True)

        # N.B. small discrepancy in system block atoms column name
        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | InputAdfBasisType | InputAmsTask         | InputAmsPropertiesNormalmodes | InputAdfXcGga | InputAmsSystem0Atoms_10                               | InputAmsSystem0Atoms_11                               |
|--------------|------|-------|----------|---------|--------|-------------------|----------------------|-------------------------------|---------------|-------------------------------------------------------|-------------------------------------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | TZP               | SinglePoint          | False                         | None          | None                                                  | None                                                  |
| dummyjob.002 | True | True  | None     | CH4     | C      | TZP               | GeometryOptimization | True                          | pbe           | None                                                  | None                                                  |
| dummyjob.003 | True | True  | None     | H2O     | O      | TZP               | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | TZP               | GeometryOptimization | False                         | pbe           | None                                                  | None                                                  |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | TZP               | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None              | GeometryOptimization | True                          | None          | None                                                  | None                                                  |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None              | SinglePoint          | False                         | None          | None                                                  | None                                                  |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None              | GeometryOptimization | True                          | None          | None                                                  | None                                                  |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None              | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.010 | True | True  | None     | None    | None   | None              | GeometryOptimization | False                         | None          | Ar 0.0000000000       0.0000000000       0.0000000000 | Ar 1.6050000000       0.9266471820       2.6050000000 |"""
        )


class TestJobAnalysisWithChemicalSystem(TestJobAnalysis):

    @pytest.fixture(scope="class")
    def dummy_single_jobs(self):
        # Generate dummy jobs for a selection of molecules and input settings
        skip_if_no_scm_libbase()
        from scm.utils.conversions import plams_molecule_to_chemsys

        smiles = ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", "Sys"]
        jobs = []
        for i, s in enumerate(smiles):
            sett = Settings()
            sett.input.ams.task = "GeometryOptimization" if i % 2 else "SinglePoint"
            sett.input.ams.Properties.NormalModes = "True" if i % 3 else "False"
            if i < 5:
                sett.input.ADF.Basis.Type = "TZP"
                if i % 2:
                    sett.input.ADF.xc.gga = "pbe"
            else:
                sett.input.DFTB
            if s == "Sys":
                sett.input.ams.System.Atoms = [
                    "Ar 0.0000000000       0.0000000000       0.0000000000",
                    "Ar 1.6050000000       0.9266471820       2.6050000000",
                ]
                mol = None
            else:
                mol = from_smiles(s)
                mol = plams_molecule_to_chemsys(mol)

            jobs.append(DummySingleJob(wait=i / 100, molecule=mol, settings=sett, name="dummyjob"))

        jm = JobManager(JobManagerSettings())
        for j in jobs:
            j.run(jobmanager=jm)
            j.ok()

        yield jobs

        jm._clean()
        shutil.rmtree(jm.workdir)

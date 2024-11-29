import pytest
from pathlib import Path

from scm.plams.unit_tests.test_basejob import DummySingleJob
from scm.plams.core.errors import ResultsError


class TestResults:
    """
    Test suite for the Results class.
    Not truly independent as relies upon the job runner/manager and results components.
    But this suite focuses on testing the methods on the results class itself.
    """

    @pytest.fixture(scope="class")
    def dummy_job(self):
        """
        Simple dummy single job.
        """
        job = DummySingleJob(name="test_results", inp="# Start\nDummy input %ID%\n# End")
        job.run()
        return job

    def test_files_after_run_as_expected(self, dummy_job):
        # Given results, when get files, then all files (except .dill) are available on the results class
        assert set(dummy_job.results.files) == {
            "test_results.out",
            "test_results.in",
            "test_results.err",
            "test_results.run",
        }

    def test_get_item_returns_file_path(self, dummy_job):
        # Given results, when get item with file name, then returns full file path
        assert dummy_job.results["test_results.in"] == str(Path(dummy_job.path) / "test_results.in")
        assert dummy_job.results["$JN.out"] == str(Path(dummy_job.path) / "test_results.out")

    def test_rename_updates_files(self):
        # Given results
        job = DummySingleJob(name="test_results_rename")
        results = job.run()

        # When rename files
        results.rename("test_results_rename.in", "new_test_results_rename.in")
        results.rename("$JN.out", "new_$JN.out")

        # Then name change applied
        assert set(results.files) == {
            "test_results_rename.run",
            "test_results_rename.err",
            "new_test_results_rename.in",
            "new_test_results_rename.out",
        }

    @pytest.mark.parametrize(
        "pattern,options,expected_match",
        [
            ("output", "", True),
            ("input", "", False),
            ("#", "-v", True),
        ],
    )
    def test_grep_file_as_expected(self, dummy_job, pattern, options, expected_match):
        # Given results, when grep output file, then get expected line matches
        grep_file = dummy_job.results.grep_file("$JN.out", pattern, options)
        grep_output = dummy_job.results.grep_output(pattern, options)

        if expected_match:
            assert grep_file == grep_output == [f"Dummy output {dummy_job.id}"]
        else:
            assert grep_file == grep_output == []

    def test_read_file_as_expected(self, dummy_job):
        # Given results, when read output file, then get expected content
        assert (
            dummy_job.results.read_file("$JN.out").replace("\r\n", "\n")
            == f"# Start\nDummy output {dummy_job.id}\n# End"
        )
        with pytest.raises(ResultsError):
            dummy_job.results.read_file("$JN.foo")

    def test_regex_file_as_expected(self, dummy_job):
        # Given results, when regex output file, then get expected matches
        assert dummy_job.results.regex_file(
            "$JN.out", "[0-9a-f]{8}-[a-f0-9]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"
        ) == [str(dummy_job.id)]
        assert dummy_job.results.regex_file("$JN.out", "u[a-z]{1}") == ["um", "ut", "ut"]

    def test_awk_file_as_expected(self, dummy_job):
        # Given results, when awk output file, then get expected content
        assert dummy_job.results.awk_file("$JN.out", script="{print $3}") == ["", str(dummy_job.id), ""]
        assert dummy_job.results.awk_output(script="{print $3}") == ["", str(dummy_job.id), ""]

    def test_get_file_chunk_as_expected(self, dummy_job):
        # Given results, when get file chunk, then file chunk as expected
        assert dummy_job.results.get_file_chunk("$JN.out") == [
            "# Start",
            f"Dummy output {dummy_job.id}",
            "# End",
        ]
        assert dummy_job.results.get_file_chunk("$JN.out", begin="# Start") == [
            f"Dummy output {dummy_job.id}",
            "# End",
        ]
        assert dummy_job.results.get_file_chunk("$JN.out", end="# End") == ["# Start", f"Dummy output {dummy_job.id}"]
        assert dummy_job.results.get_file_chunk(
            "$JN.out",
            begin="# Start",
            end="# End",
            inc_begin=True,
            inc_end=True,
            process=lambda s: s.upper() if "#" in s else s,
        ) == ["# START", f"Dummy output {dummy_job.id}", "# END"]

        assert dummy_job.results.get_output_chunk() == [
            "# Start",
            f"Dummy output {dummy_job.id}",
            "# End",
        ]
        assert dummy_job.results.get_output_chunk(begin="# Start") == [
            f"Dummy output {dummy_job.id}",
            "# End",
        ]
        assert dummy_job.results.get_output_chunk(end="# End") == ["# Start", f"Dummy output {dummy_job.id}"]
        assert dummy_job.results.get_output_chunk(
            begin="# Start",
            end="# End",
            inc_begin=True,
            inc_end=True,
            process=lambda s: s.upper() if "#" in s else s,
        ) == ["# START", f"Dummy output {dummy_job.id}", "# END"]

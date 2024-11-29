import pytest
import numpy as np

from scm.plams.core.errors import FileError
from scm.plams.tools.kftools import KFReader, KFFile, KFHistory


@pytest.fixture
def water_optimization_rkf(rkf_folder):
    return rkf_folder / "water_optimization" / "ams.rkf"


class TestKFReader:
    """
    Test suite for kf file reader
    """

    def test_read_extracts_single_and_multiple_values_as_expected(self, water_optimization_rkf):
        reader = KFReader(water_optimization_rkf)

        # Happy
        assert reader.read("General", "title") == "water_optimization"  # single string value
        assert reader.read("General", "CPUTime") == 0.732007  # single float value
        assert reader.read("Molecule", "nAtoms") == 3  # single int value
        assert not reader.read("Molecule", "eeUseChargeBroadening")  # single bool value
        assert reader.read("Molecule", "AtomSymbols").split() == ["O", "H", "H"]  # multiple string values
        assert reader.read("Molecule", "AtomicNumbers") == [8, 1, 1]  # multiple int values
        assert reader.read("Molecule", "AtomMasses") == [15.99491400, 1.00782500, 1.00782500]  # multiple float values
        assert reader.read("Molecule", "Coords") == [
            0.12646245476495446,
            0.12646245476495452,
            0.0,
            1.9124814444820302,
            -0.14921777462121424,
            0.0,
            -0.14921777462121424,
            1.9124814444820302,
            0.0,
        ]  # multiple float values

        # Unhappy
        # Case-sensitive section
        with pytest.raises(KeyError):
            reader.read("general", "cputime")

        # Case-sensitive variable
        with pytest.raises(KeyError):
            reader.read("General", "cputime")

        # Missing section and variable
        with pytest.raises(KeyError):
            reader.read("Foo", "Bar")

    def test_variable_type_gets_integer_type_codes_as_expected(self, water_optimization_rkf):
        reader = KFReader(water_optimization_rkf)

        # Happy
        assert reader.variable_type("Molecule", "nAtoms") == 1
        assert reader.variable_type("Molecule", "AtomicNumbers") == 1
        assert reader.variable_type("General", "CPUTime") == 2
        assert reader.variable_type("Molecule", "AtomMasses") == 2
        assert reader.variable_type("Molecule", "Coords") == 2
        assert reader.variable_type("General", "title") == 3
        assert reader.variable_type("Molecule", "AtomSymbols") == 3
        assert reader.variable_type("Molecule", "eeUseChargeBroadening") == 4

        # Unhappy
        # Case-sensitive section
        with pytest.raises(KeyError):
            reader.variable_type("general", "cputime")

        # Case-sensitive variable
        with pytest.raises(KeyError):
            reader.variable_type("General", "cputime")

        # Missing section and variable
        with pytest.raises(KeyError):
            reader.variable_type("Foo", "Bar")


class TestKFFile:
    """
    Test suite for kf file reading and writing
    """

    @pytest.mark.parametrize("return_as_list", [False, True])
    def test_read_existing_kf_file_extracts_single_and_multiple_values_as_expected(
        self, water_optimization_rkf, return_as_list
    ):
        file = KFFile(water_optimization_rkf, autosave=False)

        # Happy
        assert file.read("General", "title", return_as_list) == "water_optimization"  # single string value
        assert file.read("General", "CPUTime", return_as_list) == (
            [0.732007] if return_as_list else 0.732007
        )  # single float value
        assert file.read("Molecule", "nAtoms", return_as_list) == ([3] if return_as_list else 3)  # single int value
        assert file.read("Molecule", "eeUseChargeBroadening", return_as_list) == (
            [False] if return_as_list else False
        )  # single bool value
        assert file.read("Molecule", "AtomSymbols", return_as_list).split() == [
            "O",
            "H",
            "H",
        ]  # multiple string values
        assert file.read("Molecule", "AtomicNumbers", return_as_list) == [8, 1, 1]  # multiple int values
        assert file.read("Molecule", "AtomMasses", return_as_list) == [
            15.99491400,
            1.00782500,
            1.00782500,
        ]  # multiple float values
        assert file.read("Molecule", "Coords", return_as_list) == [
            0.12646245476495446,
            0.12646245476495452,
            0.0,
            1.9124814444820302,
            -0.14921777462121424,
            0.0,
            -0.14921777462121424,
            1.9124814444820302,
            0.0,
        ]  # multiple float values

        # Unhappy
        # Case-sensitive section
        with pytest.raises(KeyError):
            file.read("general", "cputime", return_as_list)

        # Case-sensitive variable
        with pytest.raises(KeyError):
            file.read("General", "cputime", return_as_list)

        # Missing section and variable
        with pytest.raises(KeyError):
            file.read("Foo", "Bar", return_as_list)

    def test_read_non_existing_file_errors(self):
        file = KFFile("not-a-file", autosave=False)
        with pytest.raises(FileError):
            file.read("foo", "bar")

    def test_write_delete_read_sections_for_new_file(self, rkf_folder):
        file = KFFile(rkf_folder / "test_kffile_write.rkf", autosave=False)

        # When write set of initial values
        file.write("Single", "TestStr", "s")
        file.write("Single", "TestBool", True)
        file.write("Single", "TestInt", 42)
        file.write("Single", "TestFloat", 42.123)
        file.write("Multi", "TestStr", ["s1", "s2"])
        file.write("Multi", "TestBool", [True, False])
        file.write("Multi", "TestInt", [42, 43])
        file.write("Multi", "TestFloat", [42.123, 0.0000001])
        file.write("ToDelete", "TestStr", "_")

        # And when update some values
        file.write("Single", "TestBool", False)
        file.write("Multi", "TestInt", [42, 44, 45])

        # And when delete a section
        file.delete_section("ToDelete")

        # Then can read sections with written values
        assert file.read_section("Single") == {"TestBool": False, "TestFloat": 42.123, "TestInt": 42, "TestStr": "s"}
        assert file.read_section("Multi") == {
            "TestBool": [True, False],
            "TestFloat": [42.123, 1e-07],
            "TestInt": [42, 44, 45],
            "TestStr": ["s1", "s2"],  # Note this is a list until saved
        }

        # And listing sections, keys and skeleton does not include deleted values
        assert file.sections() == ["Multi", "Single"]
        assert file.keys() == {"Multi", "Single"}
        assert file.get_skeleton() == {
            "Multi": {"TestFloat", "TestBool", "TestStr", "TestInt"},
            "Single": {"TestFloat", "TestBool", "TestStr", "TestInt"},
        }

    def test_write_unhappy_values_for_new_files_errors(self, rkf_folder):
        file = KFFile(rkf_folder / "test_kffile_write_unhappy.rkf", autosave=False)

        with pytest.raises(ValueError):
            file.write("Test", "InvalidType", {})

        with pytest.raises(ValueError):
            file.write("Test", "EmptyList", [])

        with pytest.raises(ValueError):
            file.write("Test", "ListDifferentTypes", [1, "two"])

        with pytest.raises(ValueError):
            file.write("Test", "ListInvalidTypes", [{}, {}])

    def test_contains(self, water_optimization_rkf):
        happy_file = KFFile(water_optimization_rkf, autosave=False)
        unhappy_file = KFFile("not-a-file", autosave=False)

        # Section present
        assert "Molecule" in happy_file
        assert "General" in happy_file
        # Section and variable present
        assert ("Molecule", "nAtoms") in happy_file
        assert ("General", "CPUTime") in happy_file
        # Section not present
        assert "Foo" not in happy_file
        assert ("Foo", "Bar") not in happy_file
        # Variable not present
        assert ("Molecule", "Bar") not in happy_file
        # No file present
        assert "Molecule" not in unhappy_file
        assert "General" not in unhappy_file
        assert ("Molecule", "nAtoms") not in unhappy_file
        assert ("General", "CPUTime") not in unhappy_file
        assert "Foo" not in unhappy_file
        assert ("Foo", "Bar") not in unhappy_file
        assert ("Molecule", "Bar") not in unhappy_file

        with pytest.raises(TypeError):
            ("a", "b", "c") in happy_file

        with pytest.raises(TypeError):
            ("a", "b", "c") in unhappy_file


class TestKFHistory:
    """
    Test suite for kf fie history reading
    """

    def test_read_all_extracts_values_as_expected(self, water_optimization_rkf):
        reader = KFReader(water_optimization_rkf)
        history = KFHistory(reader, "History")

        assert np.allclose(history.read_all("Energy"), [-5.7580155, -5.76503956, -5.76626166, -5.76628814])
        assert np.allclose(
            history.read_all("Coords"),
            np.array(
                [
                    [0.0, 0.0, 0.0, 1.88972612, 0.0, 0.0, 0.0, 1.88972612, 0.0],
                    [0.07562039, 0.07562039, 0.0, 1.90375528, -0.08964955, 0.0, -0.08964955, 1.90375528, 0.0],
                    [0.12238454, 0.12238454, 0.0, 1.90712991, -0.13978833, 0.0, -0.13978833, 1.90712991, 0.0],
                    [0.12646245, 0.12646245, 0.0, 1.91248144, -0.14921777, 0.0, -0.14921777, 1.91248144, 0.0],
                ]
            ),
        )

    def test_read_all_errors_as_expected(self, water_optimization_rkf):
        reader = KFReader(water_optimization_rkf)
        history = KFHistory(reader, "History")

        # Case sensitive
        with pytest.raises(KeyError):
            _ = KFHistory(reader, "history")

        # Section does not exist
        with pytest.raises(KeyError):
            _ = KFHistory(reader, "does-not-exist")

        # Variable case sensitive
        with pytest.raises(KeyError):
            history.read_all("coords")

        # Variable does not exist
        with pytest.raises(KeyError):
            history.read_all("does-not-exist")

    def test_iter_gets_values_when_name_present_otherwise_errors(self, water_optimization_rkf):
        reader = KFReader(water_optimization_rkf)
        history = KFHistory(reader, "History")

        assert np.allclose([v for v in history.iter("Energy")], [-5.7580155, -5.76503956, -5.76626166, -5.76628814])

        with pytest.raises(KeyError):
            _ = [v for v in history.iter("Foo")]

    @pytest.mark.parametrize("default", [None, 42.0])
    def test_iter_optional_gets_values_when_name_present_otherwise_gives_default_until_break(
        self, water_optimization_rkf, default
    ):
        reader = KFReader(water_optimization_rkf)
        history = KFHistory(reader, "History")

        assert np.allclose(
            [v for v in history.iter_optional("Energy", default=default)],
            [-5.7580155, -5.76503956, -5.76626166, -5.76628814],
        )

        for i, v in enumerate(history.iter_optional("Foo", default=default)):
            assert v == default
            i += 1
            if i > 5:
                break

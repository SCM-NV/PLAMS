import pytest
import builtins
from importlib import reload

from scm.plams.core.settings import Settings
from scm.plams.unit_tests.test_helpers import get_mock_import_function


@pytest.fixture
def ams_text_inputs():
    """
    Set of example AMS text inputs and the corresponding expected settings objects.
    """

    geometry_optimisation_input = """
Task GeometryOptimization

Engine ADF
  NumericalQuality Basic
  XC
    GGA PBE
  End
  basis
    core None
    type DZ
  End
EndEngine
"""

    geometry_optimisation_settings = Settings()
    geometry_optimisation_settings.adf.basis.core = "None"
    geometry_optimisation_settings.adf.basis.type = "DZ"
    geometry_optimisation_settings.adf.numericalquality = "Basic"
    geometry_optimisation_settings.adf.xc.gga = "PBE"
    geometry_optimisation_settings.ams.task = "GeometryOptimization"

    yield [(geometry_optimisation_input, geometry_optimisation_settings)]

    # Tear-down: reload module without the patched import function
    import scm.plams.interfaces.adfsuite.inputparser as inputparser

    reload(inputparser)


@pytest.fixture
def system_text_inputs():
    """
    Set of example system input texts and the corresponding expected settings objects.
    """

    water_system_input = """
System
  Atoms
              O       0.0000000000       0.0000000000       0.0000000000
              H       1.0000000000       0.0000000000       0.0000000000
              H       0.0000000000       1.0000000000       0.0000000000
  End
End
    """

    water_atoms_settings = Settings()
    water_atoms_settings.atoms._1 = [
        "O       0.0000000000       " "0.0000000000       0.0000000000",
        "H       1.0000000000       " "0.0000000000       0.0000000000",
        "H       0.0000000000       " "1.0000000000       0.0000000000",
    ]
    water_system_settings = Settings()
    water_system_settings.system = [water_atoms_settings]

    yield [(water_system_input, water_system_settings)]

    # Tear-down: reload module without the patched import function
    import scm.plams.interfaces.adfsuite.inputparser as inputparser

    reload(inputparser)


def test_to_settings_without_scmlibbase_succeeds(ams_text_inputs, monkeypatch):
    input_parser = get_monkeypatched_input_parser(monkeypatch)
    to_settings_succeeds(ams_text_inputs, input_parser)


def test_to_settings_with_scmlibbase_succeeds(ams_text_inputs):
    input_parser = get_input_parser_or_skip()
    to_settings_succeeds(ams_text_inputs, input_parser)


def test_to_dict_without_scmlibbase_succeeds(system_text_inputs, monkeypatch):
    input_parser = get_monkeypatched_input_parser(monkeypatch)
    to_dict_succeeds(system_text_inputs, input_parser)


def test_to_dict_with_scmlibbase_succeeds(system_text_inputs):
    input_parser = get_input_parser_or_skip()
    to_dict_succeeds(system_text_inputs, input_parser)


def get_monkeypatched_input_parser(monkeypatch):
    # Mock scm.libbase import failing (even when present in the env)
    mock_import_function = get_mock_import_function("scm.libbase")
    monkeypatch.setattr(builtins, "__import__", mock_import_function)

    # Reload the module without scm.libbase
    import scm.plams.interfaces.adfsuite.inputparser as inputparser

    reload(inputparser)
    from scm.plams.interfaces.adfsuite.inputparser import InputParserFacade, InputParser

    # Get an instance of the input parser facade using the fallback input parser
    input_parser = InputParserFacade()
    assert isinstance(input_parser.parser, InputParser)

    return input_parser


def get_input_parser_or_skip():
    from scm.plams.interfaces.adfsuite.inputparser import InputParserFacade, InputParser

    # Get an instance of the input parser facade using the scm.libbase parser
    # otherwise skip the test if the package is not loaded
    input_parser = InputParserFacade()
    if input_parser._has_scm_libbase:
        from scm.libbase import InputParser as InputParserScmLibbase

        assert isinstance(input_parser.parser, InputParserScmLibbase)
        return input_parser
    else:
        assert isinstance(input_parser.parser, InputParser)
        pytest.skip("Skipping test because optional 'scm.libbase' package is not available")


def to_settings_succeeds(input_texts, input_parser):
    for input_text, expected_settings in input_texts:
        actual_settings = input_parser.to_settings("ams", input_text)
        assert actual_settings == expected_settings


def to_dict_succeeds(input_texts, input_parser):
    for input_text, expected_settings in input_texts:
        actual_dict = input_parser.to_dict("ams", input_text)
        assert actual_dict == expected_settings.as_dict()

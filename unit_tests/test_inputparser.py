import pytest
import builtins
from importlib import reload
from threading import Thread

from scm.plams.core.settings import Settings
from scm.plams.unit_tests.test_helpers import get_mock_import_function, skip_if_no_ams_installation


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

    md_input = """
MolecularDynamics
  Checkpoint
    Frequency 1000000
  End
  InitialVelocities
    Temperature 600.0
  End
  NSteps 100000
  Shake
    All Bonds H *
  End
  Thermostat
    Tau 100
    Temperature 600.0
    Type Berendsen
  End
  TimeStep 1.0
  Trajectory
    SamplingFreq 100
  End
End
"""

    md_settings = Settings()
    md_thermo_settings = Settings()
    md_settings.ams.moleculardynamics.checkpoint.frequency = "1000000"
    md_settings.ams.moleculardynamics.initialvelocities.temperature = "600.0"
    md_settings.ams.moleculardynamics.nsteps = "100000"
    md_settings.ams.moleculardynamics.shake.all = ["Bonds H *"]
    md_thermo_settings.tau = "100"
    md_thermo_settings.temperature = "600.0"
    md_thermo_settings.type = "Berendsen"
    md_settings.ams.moleculardynamics.thermostat = [md_thermo_settings]
    md_settings.ams.moleculardynamics.timestep = "1.0"
    md_settings.ams.moleculardynamics.trajectory.samplingfreq = "100"

    yield [(geometry_optimisation_input, geometry_optimisation_settings), (md_input, md_settings)]

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

    water_settings = Settings()
    water_settings.atoms._1 = [
        "O       0.0000000000       0.0000000000       0.0000000000",
        "H       1.0000000000       0.0000000000       0.0000000000",
        "H       0.0000000000       1.0000000000       0.0000000000",
    ]
    water_system_settings = Settings()
    water_system_settings.system = [water_settings]

    graphene_system_input = """
System
   Atoms
      C   -1.2300000000000000 -2.1304224999999999  0.0000000000000000
      C   -1.2300000000000000 -0.7101408299999999  0.0000000000000000
      C    0.0000000000000000  0.0000000000000000  0.0000000000000000
      C    0.0000000000000000  1.4202816700000001  0.0000000000000000
   End
   Lattice
        2.4600000000000000   0.0000000000000000   0.0000000000000000
        0.0000000000000000   4.2608449999999998   0.0000000000000000
        0.0000000000000000   0.0000000000000000  80.0000000000000000
   End
End
"""
    graphene_settings = Settings()
    graphene_settings.atoms._1 = [
        "C   -1.2300000000000000 -2.1304224999999999  0.0000000000000000",
        "C   -1.2300000000000000 -0.7101408299999999  0.0000000000000000",
        "C    0.0000000000000000  0.0000000000000000  0.0000000000000000",
        "C    0.0000000000000000  1.4202816700000001  0.0000000000000000",
    ]
    graphene_settings.lattice._1 = [
        "2.4600000000000000   0.0000000000000000   0.0000000000000000",
        "0.0000000000000000   4.2608449999999998   0.0000000000000000",
        "0.0000000000000000   0.0000000000000000  80.0000000000000000",
    ]
    graphene_system_settings = Settings()
    graphene_system_settings.system = [graphene_settings]

    yield [(water_system_input, water_system_settings), (graphene_system_input, graphene_system_settings)]

    # Tear-down: reload module without the patched import function
    import scm.plams.interfaces.adfsuite.inputparser as inputparser

    reload(inputparser)


@pytest.fixture
def large_system_input():
    """
    Input with large molecule.
    """
    large_system_input = """
    System
       Atoms
    """
    for i in range(100):
        for j in range(100):
            large_system_input += f"""        C    {i}    {j}    0
"""
    large_system_input += """
  End
End
"""

    return large_system_input


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
    # If there is no AMS installation the input parser will not run so skip test with a warning
    skip_if_no_ams_installation()

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
    # If there is no AMS installation the input parser will not run so skip test with a warning
    skip_if_no_ams_installation()

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


def test_libbase_parser_cached(ams_text_inputs):
    input_parser = get_input_parser_or_skip()

    assert "ams" not in input_parser.parser.parsers.keys()

    input_parser.to_dict("ams", ams_text_inputs[0][0])

    assert "ams" in input_parser.parser.parsers.keys()

    assert input_parser.parser == input_parser.parser


def test_thread_safe(large_system_input):
    get_input_parser_or_skip()

    errors = []

    def parse():
        try:
            input_parser = get_input_parser_or_skip()
            input_parser.to_settings("ams", large_system_input)
        except Exception as e:
            errors.append(e)

    threads = [Thread(target=parse) for _ in range(100)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert len(errors) == 0

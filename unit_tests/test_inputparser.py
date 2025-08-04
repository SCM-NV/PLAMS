import pytest
import builtins
from importlib import reload
from threading import Thread

from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Molecule
from scm.plams.mol.atom import Atom
from scm.plams.interfaces.adfsuite.inputparser import input_to_settings
from .test_helpers import get_mock_import_function, skip_if_no_ams_installation


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
    geometry_optimisation_settings.ADF.Basis.Core = "None"
    geometry_optimisation_settings.ADF.Basis.Type = "DZ"
    geometry_optimisation_settings.ADF.NumericalQuality = "Basic"
    geometry_optimisation_settings.ADF.XC.GGA = "PBE"
    geometry_optimisation_settings.ams.Task = "GeometryOptimization"

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
    md_settings.ams.MolecularDynamics.Checkpoint.Frequency = "1000000"
    md_settings.ams.MolecularDynamics.InitialVelocities.Temperature = "600.0"
    md_settings.ams.MolecularDynamics.NSteps = "100000"
    md_settings.ams.MolecularDynamics.Shake.All = ["Bonds H *"]
    md_thermo_settings.Tau = "100"
    md_thermo_settings.Temperature = "600.0"
    md_thermo_settings.Type = "Berendsen"
    md_settings.ams.MolecularDynamics.Thermostat = [md_thermo_settings]
    md_settings.ams.MolecularDynamics.TimeStep = "1.0"
    md_settings.ams.MolecularDynamics.Trajectory.SamplingFreq = "100"

    yield [(geometry_optimisation_input, geometry_optimisation_settings), (md_input, md_settings)]

    # Tear-down: reload module without the patched import function
    import scm.plams.interfaces.adfsuite.inputparser as inputparser

    reload(inputparser)


@pytest.fixture
def system_text_inputs():
    """
    Set of example system input texts and the corresponding expected settings and molecule objects.
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
    water_settings.Atoms._1 = [
        "O       0.0000000000       0.0000000000       0.0000000000",
        "H       1.0000000000       0.0000000000       0.0000000000",
        "H       0.0000000000       1.0000000000       0.0000000000",
    ]
    water_system_settings = Settings()
    water_system_settings.System = [water_settings]

    water_mol = Molecule()
    water_mol.add_atom(Atom(symbol="O", coords=(0.0, 0.0, 0.0)))
    water_mol.add_atom(Atom(symbol="H", coords=(1.0, 0.0, 0.0)))
    water_mol.add_atom(Atom(symbol="H", coords=(0.0, 1.0, 0.0)))

    graphene_system_input = """
System Substrate
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
    graphene_settings.Atoms._1 = [
        "C   -1.2300000000000000 -2.1304224999999999  0.0000000000000000",
        "C   -1.2300000000000000 -0.7101408299999999  0.0000000000000000",
        "C    0.0000000000000000  0.0000000000000000  0.0000000000000000",
        "C    0.0000000000000000  1.4202816700000001  0.0000000000000000",
    ]
    graphene_settings.Lattice._1 = [
        "2.4600000000000000   0.0000000000000000   0.0000000000000000",
        "0.0000000000000000   4.2608449999999998   0.0000000000000000",
        "0.0000000000000000   0.0000000000000000  80.0000000000000000",
    ]
    graphene_settings._h = "Substrate"
    graphene_system_settings = Settings()
    graphene_system_settings.System = [graphene_settings]

    graphene_mol = Molecule()
    graphene_mol.add_atom(Atom(symbol="C", coords=(-1.2300000000000000, -2.1304224999999999, 0.0)))
    graphene_mol.add_atom(Atom(symbol="C", coords=(-1.2300000000000000, -0.7101408299999999, 0.0)))
    graphene_mol.add_atom(Atom(symbol="C", coords=(0.0, 0.0, 0.0)))
    graphene_mol.add_atom(Atom(symbol="C", coords=(0.0, 1.4202816700000001, 0.0)))
    graphene_mol.lattice = [(2.46, 0.0, 0.0), (0.0, 4.2608449999999998, 0.0), (0.0, 0.0, 80.0)]

    yield [
        (water_system_input, water_system_settings, water_mol),
        (graphene_system_input, graphene_system_settings, graphene_mol),
    ]

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


def test_input_to_settings_without_scmlibbase_succeeds(ams_text_inputs, monkeypatch):
    input_parser = get_monkeypatched_input_parser(monkeypatch)
    input_to_settings_succeeds(ams_text_inputs, input_parser)


def test_input_to_settings_with_scmlibbase_succeeds(ams_text_inputs):
    input_parser = get_input_parser_or_skip()
    input_to_settings_succeeds(ams_text_inputs, input_parser)


def test_input_to_settings_using_default_parser_succeeds(ams_text_inputs):
    skip_if_no_ams_installation()
    input_to_settings_succeeds(ams_text_inputs, None)


def test_to_dict_without_scmlibbase_succeeds(system_text_inputs, monkeypatch):
    input_parser = get_monkeypatched_input_parser(monkeypatch)
    to_dict_succeeds(system_text_inputs, input_parser)


def test_to_dict_with_scmlibbase_succeeds(system_text_inputs):
    input_parser = get_input_parser_or_skip()
    to_dict_succeeds(system_text_inputs, input_parser)


def test_get_system_blocks_from_input_as_molecules(system_text_inputs):
    # If there is no AMS installation the input file will not be present so skip test with a warning
    skip_if_no_ams_installation()

    from scm.plams.interfaces.adfsuite.inputparser import get_system_blocks_as_molecules_from_input
    from scm.libbase import InputFile

    # Combine system blocks into one text
    text_input = "\n".join([s[0] for s in system_text_inputs])
    expected = {s[1].System[0]._h or "": s[2] for s in system_text_inputs}

    # Given input file
    input_file = InputFile(program="ams", text_input=text_input)

    # When get system blocks as molecules
    actual = get_system_blocks_as_molecules_from_input(input_file)

    # Then molecules are under key corresponding to header of the block
    # and match the expected geometry
    assert list(actual.keys()) == ["", "Substrate"]
    for k, v in actual.items():
        assert v.lattice == expected[k].lattice
        assert v.label(4) == expected[k].label(4)


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


def input_to_settings_succeeds(input_texts, input_parser):
    for input_text, expected_settings in input_texts:
        actual_settings = input_to_settings(input_text, parser=input_parser)
        assert actual_settings == expected_settings


def to_dict_succeeds(input_texts, input_parser):
    for input_text, expected_settings, _ in input_texts:
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
            input_to_settings(large_system_input, parser=input_parser)
        except Exception as e:
            errors.append(e)

    threads = [Thread(target=parse) for _ in range(100)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert len(errors) == 0

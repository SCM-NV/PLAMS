import os
import json
import threading
from typing import List, Union, Tuple, Optional, Dict, TYPE_CHECKING

from scm.plams.core.settings import Settings
from scm.plams.core.errors import PlamsError
from scm.plams.mol.molecule import Molecule
from scm.plams.core.functions import requires_optional_package
from scm.plams.interfaces.adfsuite.amsworker import AMSWorker, AMSWorkerError
from scm.plams.interfaces.adfsuite.ams import AMSJob

if TYPE_CHECKING:
    from scm.libbase import InputParser as InputParserLibbase
    from scm.libbase import InputFile as InputFileLibbase

__all__: List[str] = ["get_system_blocks_as_molecules_from_input", "input_to_settings"]


class InputParser:
    """
    A utility class for converting text input into JSON dictionaries and plams.Settings.

    This is a legacy implementation for environments without access to scm.libbase.
    """

    # !!!!!!!  DEPRECATED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # This class has been deprecated. The equivalent in scm.libbase
    # should be used instead where possible. The scm.libbase version does the input
    # parsing via direct calls into libscm_base, instead of spawning an AMSWorker
    # and then pushing the input through the pipe. This implementation exists here
    # only to remove the dependency on scm.libbase when running in python
    # environments without access to the base library.
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    def __init__(self):
        sett = Settings()
        sett.input.LennardJones = Settings()
        sett.runscript.nproc = 1
        sett.runscript.preamble_lines = [f'export AMS_INPUTREADER_ROOTPATH="{os.getcwd()}"']
        self.worker = AMSWorker(sett)

    def __enter__(self):
        return self

    def stop(self, keep_workerdir=False):
        if self.worker is not None:
            self.worker.stop(keep_workerdir)

    def __exit__(self, *args):
        self.stop()

    def to_dict(self, program, text_input, string_leafs=True):
        """Run a string of text input through the input parser and produce a Python dictionary representing the JSONified input."""
        try:
            json_input = self.worker.ParseInput(program, text_input, string_leafs)
        except AMSWorkerError as exc:
            raise ValueError(f"Input parsing failed. {exc.get_errormsg()}") from exc
        return json.loads(json_input)

    # Renamed for alignment with the libbase InputParser and yet to maintain backwards compatibility
    _run = to_dict

    def to_settings(self, program, text_input):
        """Transform a string with text input into a PLAMS Settings object."""
        return input_to_settings(text_input, program=program, parser=self)

    @staticmethod
    def separate_engine_lines(lines):
        """
        Separate the engine lines from other lines of AMS input text
        """
        return _separate_engine_lines(lines)


class InputParserFacade:
    """
    A utility class for converting text input into JSON dictionaries and plams.Settings.

    Uses the scm.libbase implementation of InputParser if available, and otherwise the legacy implementation
    of the parser which spawns an AMSWorker instance.
    """

    try:
        from scm.libbase import InputParser as InputParserLibbase

        # Cache a single instance of the parser to avoid having to repeatedly reload input file definition JSON
        # But for this need to make access to the parser thread-safe
        input_parser_scm_libbase = InputParserLibbase()
        input_parser_lock = threading.Lock()
        _has_scm_libbase = True
    except ImportError:
        _has_scm_libbase = False

    @property
    def parser(self):
        """
        Get instance of a parser used to convert text input.
        """
        if self._has_scm_libbase:
            return self.input_parser_scm_libbase
        else:
            return InputParser()

    def to_dict(self, program: str, text_input: str, string_leafs: bool = True):
        """
        Run a string of text input through the input parser and produce a Python dictionary representing the JSONified input.
        """
        if self._has_scm_libbase:
            with self.input_parser_lock:
                return self.parser.to_dict(program, text_input, string_leafs)
        else:
            with self.parser as parser:
                return parser.to_dict(program, text_input, string_leafs)


@requires_optional_package("scm.libbase")
def get_system_blocks_as_molecules_from_input(input_file: "InputFileLibbase") -> Dict[str, Molecule]:
    """
    Get a dictionary of mappings between the System blocks in the input to |Molecule| instances.

    The keys in the dictionary are the names of the systems from the headers of the blocks, e.g.::

        System initial
            ...
        End

    would be returned as dict["initial"].
    The main system without a name in the header will have the empty string as key.

    :param input_file: input file containing the text with the System blocks
    :return: dictionary of mappings between the system block names and molecules
    """
    if not input_file.is_defined("System") and input_file.is_defined("LoadSystem"):
        raise PlamsError(
            f"Attempting to read System blocks for program '{input_file.program}' that does not have them defined"
        )

    if input_file.number_of_entries("System") > 0:
        tmp = Settings()
        tmp.input.ams.System = Settings(json.loads(input_file.get_json())).System
        mols = AMSJob.settings_to_mol(tmp)
    else:
        mols = {}

    for i in range(input_file.number_of_entries("LoadSystem")):
        header = input_file.get_header(f"LoadSystem[{i}]")
        if header in mols:
            raise PlamsError(f"Input error: duplicate system headers found in input: {header}")
        mols[header] = Molecule()
        mols[header].readrkf(
            input_file.get_string(f"LoadSystem[{i}]%File"), input_file.get_string(f"LoadSystem[{i}]%Section")
        )

    return mols


def input_to_settings(
    text_input: str,
    program: str = "ams",
    string_leafs: bool = True,
    parser: Optional[Union["InputParserLibbase", InputParser, InputParserFacade]] = None,
) -> "Settings":
    """
    Transform a string with text input into a |Settings| object.

    For AMS driver input the root level input will be set under settings.ams,
    while the engine block will be parsed separately and returned under settings.%engine%,
    where %engine% is the name of the engine, e.g. adf or dftb.

    :param text_input: AMS text input to covert to settings
    :param program: name of the program, defaults to ``ams``
    :param string_leafs: include string leafs, defaults to ``True``
    :param parser: use specific parser or defaults to internal parser if ``None``
    :return: |Settings| object with the structure of the parsed AMS input
    """
    input_parser = parser or InputParserFacade()
    if program in ["ams", "acerxn"]:
        # Settings for the program are special:
        # * Root level input needs to go under settings.input.ams.
        # * Engine block needs to go  to settings.input.%engine% where
        #   %engine% is the name of the engine, e.g. adf.
        input_settings = Settings()
        lines = text_input.splitlines()

        # Find the lines corresponding to the engine block.
        while True:
            lines, engine_lines = _separate_engine_lines(lines)
            if engine_lines is None:
                break
            # We have found a separate engine block.
            engine_name = engine_lines[0].split()[1]
            if len(engine_lines) == 2:
                # If it's empty we already know the result of parsing it.
                input_settings[engine_name] = Settings()
            else:
                input_settings[engine_name] = Settings(
                    input_parser.to_dict(engine_name.lower(), "\n".join(engine_lines[1:-1]), string_leafs)
                )

        input_settings["ams"] = Settings(input_parser.to_dict(program, "\n".join(lines), string_leafs))

    else:
        input_settings = Settings(input_parser.to_dict(program, text_input, string_leafs))

    return input_settings


def _separate_engine_lines(lines: List[str]) -> Tuple[List[str], Optional[List[str]]]:
    """
    Separate the engine lines from other lines of AMS input text.

    :param lines: lines of AMS input text
    :return: tuple consisting of the original lines without the engine lines and the engine lines themselves
    """

    # Find the lines corresponding to the engine block.
    engine_lines = None
    engine_first, engine_last = None, None
    inner_engines = 0
    ends = 0
    for i, line in enumerate(lines):
        if line == "" or line.isspace():
            continue
        if line.split(maxsplit=1)[0].lower() == "engine":
            if engine_first is None:
                engine_first = i
                continue
            else:
                inner_engines += 1
        elif line.split(maxsplit=1)[0].lower() == "endengine":
            ends += 1
        if ends > inner_engines:
            engine_last = i
            break

    if engine_first is not None and engine_last is not None:
        # We have found a separate engine block.
        engine_lines = lines[engine_first : engine_last + 1]
        # We have dealt with the engine lines, let's remove them from the input.
        del lines[engine_first : engine_last + 1]

    return lines, engine_lines

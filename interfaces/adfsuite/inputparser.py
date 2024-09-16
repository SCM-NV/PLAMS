import os
import json
import threading

from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.amsworker import AMSWorker, AMSWorkerError

__all__ = []

# !!!!!!!  DEPRECATED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# The class in this module has been deprecated. The equivalent in scm.libbase
# should be used instead where possible. The scm.libbase version does the input
# parsing via direct calls into libscm_base, instead of spawning an AMSWorker
# and then pushing the input through the pipe. This implementation exists here
# only to remove the dependency on scm.libbase when running in python
# environments without access to the base library.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


class InputParser:
    """
    A utility class for converting text input into JSON dictionaries and plams.Settings.

    This is a legacy implementation for environments without access to scm.libbase.
    """

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

    def _run(self, program, text_input, string_leafs=True):
        """Run a string of text input through the input parser and produce a Python dictionary representing the JSONified input."""
        try:
            json_input = self.worker.ParseInput(program, text_input, string_leafs)
        except AMSWorkerError as exc:
            raise ValueError(f"Input parsing failed. {exc.get_errormsg()}") from exc
        return json.loads(json_input)

    def to_settings(self, program, text_input):
        """Transform a string with text input into a PLAMS Settings object."""

        if program in ["ams", "acerxn"]:
            input_settings = Settings()
            lines = text_input.splitlines()

            # Find the lines corresponding to the engine block.
            while 1:
                lines, engine_lines = self.separate_engine_lines(lines)
                if engine_lines is None:
                    break
                # We have found a separate engine block.
                engine_name = engine_lines[0].split()[1].lower()
                if len(engine_lines) == 2:
                    # If it's empty we already know the result of parsing it.
                    input_settings[engine_name] = Settings()
                else:
                    input_settings[engine_name] = Settings(self._run(engine_name, "\n".join(engine_lines[1:-1])))

            input_settings["ams"] = Settings(self._run(program, "\n".join(lines)))

        else:
            input_settings = Settings(self._run(program, text_input))

        return input_settings

    @staticmethod
    def separate_engine_lines(lines):
        """
        Separate the engine lines from other lines of AMS input text
        """
        # Find the lines corresponding to the engine block.
        engine_lines = None
        engine_first, engine_last = None, None
        inner_engines = 0
        ends = 0
        for i, l in enumerate(lines):
            if l == "" or l.isspace():
                continue
            if l.split(maxsplit=1)[0].lower() == "engine":
                if engine_first is None:
                    engine_first = i
                    continue
                else:
                    inner_engines += 1
            elif l.split(maxsplit=1)[0].lower() == "endengine":
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


class InputParserFacade:
    """
    A utility class for converting text input into JSON dictionaries and plams.Settings.

    Uses the scm.libbase implementation of InputParser if available, and otherwise the legacy implementation
    of the parser which spawns an AMSWorker instance.
    """

    try:
        from scm.libbase import InputParser as InputParserScmLibbase

        # Cache a single instance of the parser to avoid having to repeatedly reload input file definition JSON
        # But for this need to make access to the parser thread-safe
        input_parser_scm_libbase = InputParserScmLibbase()
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

    def to_settings(self, program: str, text_input: str):
        """Transform a string with text input into a PLAMS Settings object."""

        if self._has_scm_libbase:
            with self.input_parser_lock:
                return self.parser.to_settings(program, text_input)
        else:
            with self.parser as parser:
                return parser.to_settings(program, text_input)

    def to_dict(self, program: str, text_input: str, string_leafs: bool = True):
        """
        Run a string of text input through the input parser and produce
        a Python dictionary representing the JSONified input.
        """
        if self._has_scm_libbase:
            with self.input_parser_lock:
                return self.parser.to_dict(program, text_input, string_leafs)
        else:
            with self.parser as parser:
                return parser._run(program, text_input, string_leafs)

from typing import Dict, Union, Optional, KeysView, List, Any, Tuple, NoReturn, TYPE_CHECKING, cast
from typing_extensions import LiteralString

if TYPE_CHECKING:
    from scm.plams.tools.kftools import KFFile, TRead
    from scm.plams.mol.molecule import Atom

from scm.plams.core.errors import FileError, PlamsError
from scm.plams.interfaces.adfsuite.scmjob import SCMJob, SCMResults
from scm.plams.interfaces.adfsuite.inputparser import input_to_settings
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Molecule
from scm.plams.core.functions import log

__all__ = ["AMSAnalysisJob", "AMSAnalysisResults", "convert_to_unicode"]

try:
    from scm.pisa.block import DriverBlock

    _has_scm_pisa = True
except ImportError:
    _has_scm_pisa = False


class AMSAnalysisPlot:
    """
    Class representing a plot of 2D or higher

    * ``x``       -- A list of lists containing the values in each of the multiple x-axes
    * ``y``       -- A list containing the values along the y-axis
    * ``y_sigma`` -- A list containing the standard deviation of the values on the y-axis
    * ``name``    -- The name of the plot

    The most important method is the write method, which returns a string containing all the plot info,
    and can also write a corresponding file if a filename is provided as argument.
    This file can be read by e.g. gnuplot.
    """

    def __init__(self) -> None:
        """
        Initiate an instance of the plot class
        """
        self.x: List[List[float]] = []
        self.x_units: List[str] = []
        self.x_names: List[str] = []

        self.y: Optional[List[float]] = None
        self.y_units: Optional[str] = None
        self.y_name: Optional[str] = None
        self.y_sigma: Optional[List[float]] = None

        self.properties: Optional[Dict[str, TRead]] = None
        self.name: Optional[str] = None
        self.section: Optional[str] = None

    def read_data(self, kf: "KFFile", sec: str) -> None:
        """
        Read the xy data for a section from the kf file
        """
        # Read all the x-values. There can be multiple axes for ND plots (n=3,4,....)
        sections = kf.get_skeleton()
        xkeys = [k for k in sections[sec] if "x(" in k and ")-axis" in k]
        xnums = sorted([int(k.split("(")[1].split(")")[0]) for k in xkeys])
        xnums = sorted([xnum for xnum in set(xnums)])
        for i in xnums:
            xkey = f"x({i})-axis"
            self.x.append(kf.read_reals(sec, xkey))
            x_name = kf.read_string(sec, f"{xkey}(label)")
            self.x_names.append(convert_to_unicode(x_name))
            self.x_units.append(convert_to_unicode(kf.read_string(sec, f"{xkey}(units)")))

        # Read the y-values
        ykey = "y-axis"
        y_name = kf.read_string(sec, f"{ykey}(label)")
        self.y = kf.read_reals(sec, ykey)
        self.y_name = convert_to_unicode(y_name)
        self.y_units = convert_to_unicode(kf.read_string(sec, f"{ykey}(units)"))

        self.y_sigma = kf.read_reals(sec, "sigma")

        self.read_properties(kf, sec)
        self.section = sec.split("(")[0] + "_" + sec.split("(")[1].split(")")[0]
        self.name = self.section

    def read_properties(self, kf: "KFFile", sec: str) -> None:
        """
        Read properties from the KF file
        """
        counter = 0
        properties: Dict[str, TRead] = {}
        while 1:
            counter += 1
            try:
                propname = kf.read_string(sec, f"Property({counter})").strip()
            except:
                break
            prop = kf.read(sec, propname)
            if isinstance(prop, str):
                prop = convert_to_unicode(prop.strip())
            properties[propname] = prop

        # Now set the instance variables
        self.properties = properties
        if "Legend" in properties:
            self.name = cast(str, properties["Legend"])

    def get_dimensions(self) -> int:
        """
        Get the dimensonality of the plot
        """
        return len(self.x)

    def write(self, outfilename: Optional[str] = None) -> LiteralString:
        """
        Print this plot to a text file
        """
        # Place property string
        parts = []
        properties = self.properties if self.properties is not None else {}
        for propname, prop in properties.items():
            parts.append(f"{propname:<30} {prop}\n")

        # Place the string with the column names
        x_name = ""
        for xname, xunit in zip(self.x_names, self.x_units):
            x_str = f"{xname}({xunit})"
            x_name += f"{x_str:>30} "
        y_name = f"{self.y_name}({self.y_units})"
        parts.append(f"{x_name} {y_name:>30} {'sigma':>30}\n")

        # Determine the number of values per axis
        ndims = len(self.x)
        axis_length = int(len(self.x[0]) ** (1.0 / ndims))

        # Place the values
        value_lists = self.x + [self.y] + [self.y_sigma]
        for i, values in enumerate(zip(*value_lists)):
            v_str = ""
            for v in values:
                v_str += f"{v:30.10e} "
            v_str += "\n"
            if (i + 1) % axis_length == 0:
                v_str += "\n"
            parts.append(v_str)
        block = "".join(parts)

        if outfilename is not None:
            outfile = open(outfilename, "w", encoding="utf8")
            outfile.write(block)
            outfile.close()

        return block

    @classmethod
    def from_kf(cls, kf: "KFFile", section: str, i: int = 1) -> "AMSAnalysisPlot":
        xy = cls()

        # Find the correct section in the KF file
        sections = kf.sections()
        matches = [s for s in sections if s.lower() == f"{section.lower()}({i})"]
        if len(matches) == 0:
            print("Sections: ", list(sections))
            raise PlamsError(
                'AMSAnalysisResults.get_xy(section,i): section must be one of the above. You specified "{}"'.format(
                    section
                )
            )
        sec = matches[0]

        # Get the data
        xy.read_data(kf, sec)
        return xy


class AMSAnalysisResults(SCMResults):
    _kfext = ".kf"
    _rename_map = {"plot.kf": "$JN" + _kfext}

    def get_molecule(self, *args: Any, **kwargs: Any) -> NoReturn:  # type: ignore[override]
        raise PlamsError("AMSAnalysisResults does not support the get_molecule() method.")

    def get_sections(self) -> KeysView[str]:
        """
        Read the sections available to make xy plots
        """
        if not self._kfpresent():
            raise FileError("File {} not present in {}".format(self.job.name + self.__class__._kfext, self.job.path))
        if self._kf.reader._sections is None:  # type: ignore[union-attr]
            self._kf.reader._create_index()  # type: ignore[union-attr]
        return self._kf.reader._sections.keys()  # type: ignore[union-attr]

    def get_xy(self, section: str = "", i: int = 1) -> AMSAnalysisPlot:
        """
        Get the AMSAnalysisPlot object for a specific section of the plot KFFile
        """
        if isinstance(self.job.settings.input, Settings):
            task = self.job.settings.input.Task
        else:
            task = self.job.settings.input.Task.val
        if section == "":
            section = task

        if not self._kfpresent():
            raise FileError("File {} not present in {}".format(self.job.name + self.__class__._kfext, self.job.path))
        xy = AMSAnalysisPlot.from_kf(self._kf, section, i)
        return xy

    def get_all_plots(self) -> List[AMSAnalysisPlot]:
        """
        Get a list of all the plot objects created by the analysis jobs
        """
        sections = self.get_sections()
        plots: List[AMSAnalysisPlot] = []
        for section in sections:
            if section == "General":
                continue
            if "History" in section:
                continue
            name_part = section.split("(")[0]
            num_part = int(section.split("(")[1].split(")")[0])
            xy = self.get_xy(name_part, num_part)
            plots.append(xy)
        return plots

    def write_all_plots(self) -> None:
        """
        Write all the plots created by the analysis job to file
        """
        plots = self.get_all_plots()
        for xy in plots:
            xy.write(f"{xy.section}.dat")

    def get_D(self, i: int = 1) -> Tuple[Optional[float], Optional[str]]:
        """returns a 2-tuple (D, D_units) from the AutoCorrelation(i) section on the .kf file."""

        # If there are multiple, it will read the first one
        sections = [sec for sec in self.get_sections() if "Integral" in sec]
        if len(sections) < i:
            return None, None
        section = sections[i - 1]
        plot = self.get_xy(section.split("(")[0], i)
        if not plot.properties or "DiffusionCoefficient" not in plot.properties.keys():
            return None, None

        D = cast(float, plot.properties["DiffusionCoefficient"])
        D_units = plot.y_units
        return D, D_units

    def recreate_settings(self) -> Optional[Settings]:
        """Recreate the input |Settings| instance for the corresponding job based on files present in the job folder. This method is used by |load_external|.

        Extract user input from the kf file and parse it back to a |Settings| instance using ``scm.base`` module. Remove the ``system`` branch from that instance.
        """
        user_input = self._kf.read_string("General", "user input")
        try:
            inp = input_to_settings(user_input, program="analysis")
        except:
            log(
                "Failed to recreate input settings from {}".format(
                    str(self.job.get_path() / (self.job.name + self.__class__._kfext))
                )
            )
            return None
        s = Settings()
        s.input = inp
        return s

    def recreate_molecule(self) -> Union[None, Molecule, Dict[str, Molecule]]:
        """Recreate the input molecule(s) for the corresponding job based on files present in the job folder.

        This method is used by |load_external|.
        It extracts data from the ``InputMolecule`` and ``InputMolecule(*)`` sections.
        """
        from scm.plams import AMSJob

        if "system" not in self.job.settings.input:
            return None

        self.job.settings.input.ams.system = self.job.settings.input.system
        del self.job.settings.input.system
        molecule = AMSJob.settings_to_mol(self.job.settings)
        del self.job.settings.input.ams
        return molecule


class AMSAnalysisJob(SCMJob):
    """A class for analyzing molecular dynamics trajectories using the ``analysis`` program."""

    results: AMSAnalysisResults
    _result_type = AMSAnalysisResults
    _command = "analysis"
    _subblock_end = "end"

    def __init__(self, **kwargs: Any):
        SCMJob.__init__(self, **kwargs)

    def _serialize_mol(self) -> None:
        """
        Use the method from AMSJob to move the molecule to the settings object
        """
        from scm.plams import AMSJob

        systems = AMSJob._serialize_molecule(self)  # type: ignore[arg-type]
        if len(systems) > 0:
            if _has_scm_pisa and isinstance(self.settings.input, DriverBlock):
                self.settings.system = systems
            else:
                self.settings.input.system = systems

    def _remove_mol(self) -> None:
        """
        Remove the molecule from the system block again
        """
        if hasattr(self.settings.input, "system"):
            del self.settings.input.system
        elif "system" in self.settings:
            del self.settings.system

    @staticmethod
    def _atom_suffix(atom: "Atom") -> str:
        """
        Return the suffix of an atom.
        """
        from scm.plams import AMSJob

        return AMSJob._atom_suffix(atom)

    def check(self) -> bool:
        try:
            grep = self.results.grep_file("$JN.err", "NORMAL TERMINATION")
        except:
            return False
        return len(grep) > 0


def convert_to_unicode(k: str) -> str:
    """
    Convert a string with ascii symbols representing unicode symbols

    Example k: 'abc\\u03c9def'
    """
    parts = k.split("\\u")
    # Collect the hexadecimals
    symbols = [chr(int(part[:4], 16)) for part in parts[1:]]
    # Now repair the parts
    parts = parts[:1] + ["".join([s, part[4:]]) for s, part in zip(symbols, parts[1:])]
    key = "".join(parts)

    return key

import os
from typing import Dict, Union, Optional

from scm.plams.core.errors import FileError, PlamsError
from scm.plams.interfaces.adfsuite.scmjob import SCMJob, SCMResults
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Molecule
from scm.plams.core.functions import log

__all__ = ["AMSAnalysisJob", "AMSAnalysisResults", "convert_to_unicode"]


class AMSAnalysisPlot:
    """
    Class representing a plot of 2D or higher

    * ``x``       -- A list of lists containing the values in each of the multiple x-axes
    * ``y``       -- A list containing the values along the y-axis
    * ``y_sigma`` -- A list containing the standard deviation of the values onthe y-axis
    * ``name``    -- The name of the plot

    The most important method is the write method, which returns a string containing all the plot info,
    and can also write a corresponding file if a filename is provided as argument.
    This file can be read by e.g. gnuplot.
    """

    def __init__(self):
        """
        Initiate an instance of the plot class
        """
        self.x = []
        self.x_units = []
        self.x_names = []

        self.y = None
        self.y_units = None
        self.y_name = None
        self.y_sigma = None  # standard deviation for y_values

        self.properties: Optional[Dict] = None
        self.name = None
        self.section = None

    def read_data(self, kf, sec):
        """
        Read the xy data for a section from the kf file
        """
        # Read all the x-values. There can be multiple axes for ND plots (n=3,4,....)
        sections = kf.get_skeleton()
        xkeys = [k for k in sections[sec] if "x(" in k and ")-axis" in k]
        xnums = sorted([int(k.split("(")[1].split(")")[0]) for k in xkeys])
        xnums = sorted([xnum for xnum in set(xnums)])
        for i in xnums:
            xkey = "x(%i)-axis" % (i)
            self.x.append(kf.read(sec, xkey))
            x_name = kf.read(sec, "%s(label)" % (xkey))
            self.x_names.append(convert_to_unicode(x_name))
            self.x_units.append(convert_to_unicode(kf.read(sec, "%s(units)" % (xkey))))

        # Read the y-values
        ykey = "y-axis"
        y_name = kf.read(sec, "%s(label)" % (ykey))
        self.y = kf.read(sec, ykey)
        self.y_name = convert_to_unicode(y_name)
        self.y_units = convert_to_unicode(kf.read(sec, "%s(units)" % (ykey)))

        self.y_sigma = kf.read(sec, "sigma")

        self.read_properties(kf, sec)
        self.section = sec.split("(")[0] + "_" + sec.split("(")[1].split(")")[0]
        self.name = self.section

    def read_properties(self, kf, sec):
        """
        Read properties from the KF file
        """
        counter = 0
        properties = {}
        while 1:
            counter += 1
            try:
                propname = kf.read(sec, "Property(%i)" % (counter)).strip()
            except:
                break
            properties[propname] = kf.read(sec, propname)
            if isinstance(properties[propname], str):
                properties[propname] = properties[propname].strip()
                properties[propname] = convert_to_unicode(properties[propname])

        # Now set the instance variables
        self.properties = properties
        if "Legend" in properties:
            self.name = properties["Legend"]

    def get_dimensions(self):
        """
        Get the dimensonality of the plot
        """
        return len(self.x)

    def write(self, outfilename=None):
        """
        Print this plot to a text file
        """
        # Place property string
        parts = []
        properties = self.properties if self.properties is not None else {}
        for propname, prop in properties.items():
            parts.append("%-30s %s\n" % (propname, prop))

        # Place the string with the column names
        x_name = ""
        for xname, xunit in zip(self.x_names, self.x_units):
            x_str = "%s(%s)" % (xname, xunit)
            x_name += "%30s " % (x_str)
        y_name = "%s(%s)" % (self.y_name, self.y_units)
        parts.append("%s %30s %30s\n" % (x_name, y_name, "sigma"))

        # Determine the number of values per axis
        ndims = len(self.x)
        axis_length = int(len(self.x[0]) ** (1.0 / ndims))

        # Place the values
        value_lists = self.x + [self.y] + [self.y_sigma]
        for i, values in enumerate(zip(*value_lists)):
            v_str = ""
            for v in values:
                v_str += "%30.10e " % (v)
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
    def from_kf(cls, kf, section, i=1):
        xy = cls()

        # Find the correct section in the KF file
        sections = kf.sections()
        matches = [s for s in sections if s.lower() == section.lower() + "(%i)" % (i)]
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

    def get_molecule(self, *args, **kwargs):
        raise PlamsError("AMSAnalysisResults does not support the get_molecule() method.")

    def get_sections(self):
        """
        Read the sections available to make xy plots
        """
        if not self._kfpresent():
            raise FileError("File {} not present in {}".format(self.job.name + self.__class__._kfext, self.job.path))
        if self._kf.reader._sections is None:
            self._kf.reader._create_index()
        return self._kf.reader._sections.keys()  # type: ignore

    def get_xy(self, section="", i=1):
        """
        Get the AMSAnalysisPlot object for a specific section of the plot KFFile
        """
        task = self.job.settings.input.Task
        if section == "":
            section = task

        if not self._kfpresent():
            raise FileError("File {} not present in {}".format(self.job.name + self.__class__._kfext, self.job.path))
        xy = AMSAnalysisPlot.from_kf(self._kf, section, i)
        return xy

    def get_all_plots(self):
        """
        Get a list of all the plot objects created by the analysis jobs
        """
        sections = self.get_sections()
        plots = []
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

    def write_all_plots(self):
        """
        Write all the plots created by the analysis job to file
        """
        plots = self.get_all_plots()
        for xy in plots:
            xy.write("%s" % (xy.section + ".dat"))

    def get_D(self, i=1):
        """returns a 2-tuple (D, D_units) from the AutoCorrelation(i) section on the .kf file."""

        # If there are multiple, it will read the first one
        sections = [sec for sec in self.get_sections() if "Integral" in sec]
        if len(sections) < i:
            return None, None
        section = sections[i - 1]
        plot = self.get_xy(section.split("(")[0], i)
        if not "DiffusionCoefficient" in plot.properties.keys():
            return None, None

        D = plot.properties["DiffusionCoefficient"]
        D_units = plot.y_units
        return D, D_units

    def recreate_settings(self):
        """Recreate the input |Settings| instance for the corresponding job based on files present in the job folder. This method is used by |load_external|.

        Extract user input from the kf file and parse it back to a |Settings| instance using ``scm.libbase`` module. Remove the ``system`` branch from that instance.
        """
        user_input = self._kf.read("General", "user input")
        try:
            from scm.libbase import InputParser

            inp = InputParser().to_settings("analysis", user_input)
        except:
            log(
                "Failed to recreate input settings from {}".format(
                    os.path.join(self.job.path, "".join([self.job.name, self.__class__._kfext]))
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

    _result_type = AMSAnalysisResults
    _command = "analysis"
    _subblock_end = "end"

    def __init__(self, **kwargs):
        SCMJob.__init__(self, **kwargs)

    def _serialize_mol(self):
        """
        Use the method from AMSJob to move the molecule to the settings object
        """
        from scm.plams import AMSJob

        systems = AMSJob._serialize_molecule(self)
        if len(systems) > 0:
            self.settings.input.system = systems

    def _remove_mol(self):
        """
        Remove the molecule from the system block again
        """
        if "system" in self.settings.input:
            del self.settings.input.system

    @staticmethod
    def _atom_suffix(atom):
        """
        Return the suffix of an atom.
        """
        from scm.plams import AMSJob

        return AMSJob._atom_suffix(atom)

    def check(self):
        try:
            grep = self.results.grep_file("$JN.err", "NORMAL TERMINATION")
        except:
            return False
        return len(grep) > 0


def convert_to_unicode(k):
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

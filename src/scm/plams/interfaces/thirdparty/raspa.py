"""
Run RASPA2 with PLAMS
Contributed by Patrick Melix

See the Documentation for an Example
"""

import shutil
from os import symlink
from os.path import join as opj

from scm.plams.core.basejob import SingleJob
from scm.plams.core.errors import JobError, ResultsError
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings

__all__ = ["RaspaJob", "RaspaResults"]


class RaspaResults(Results):
    """A Class for handling RASPA Results.

    Note that the Output file is guessed in ``self.collect``. Do not rely on it to point to the
    file you need and but check ``self._filenames['out']``!
    """

    def collect(self):
        """Try to Guess the Main Output File.

        Set it to whatever comes first in the 'Output/System_0' folder.
        """
        Results.collect(self)
        if "out" not in self.job._filenames:
            files = [f for f in self.files if "Output/System_0" in f]
            if files:
                self.job._filenames["out"] = files[0]

    def get_block_value(self, search, file=None, get_std=False, get_unit=False):
        """Return a Block Value from Output File.

        Can be used to retrieve block values from the final section of the output files.
        Uses `Results.grep_file` to retrieve the values.
        Remember to escape special `grep` characters in the search string.
        Returns the average value over all blocks.

        search: string
        String labeling the value you desire.

        file: string
        Defaults to the automatically detected output file. See `self.collect()`.

        get_std: boolean
        Also return the standard deviation of the value

        get_unit: boolean
        Also return the unit of the value
        """
        if not file:
            file = self.job._filenames["out"]

        chunk = self.grep_file(file, pattern=search, options="-A 8")
        last_line = chunk[-1].strip()
        if not last_line.startswith("Average"):
            raise ResultsError("Unable to get the requested value from {:}".format(file))

        last_line = last_line.split()
        ret = [float(last_line[1])]

        if get_std:
            ret.append(float(last_line[4]))
        if get_unit:
            ret.append(last_line[-1][1:-1])

        if len(ret) > 1:
            return tuple(ret)
        else:
            return ret[0]

    def get_value(self, search, file=None, get_std=False, get_unit=False):
        """Return a Single Value from Output File.

        Can be used to retrieve single values from the final section of the output files.
        Uses `Results.grep_file` to retrieve the values.
        Remember to escape special `grep` characters in the search string.

        search: string
        String labeling the value you desire.

        file: string
        Defaults to the automatically detected output file. See `self.collect()`.

        get_std: boolean
        Also return the standard deviation of the value

        get_unit: boolean
        Also return the unit of the value
        """
        if not file:
            file = self.job._filenames["out"]

        line = self.grep_file(file, pattern=search)[-1].split()
        # Two cases, with standard deviation and without
        if "+/-" in line:
            idx = line.index("+/-") - 1
            ret = [float(line[idx])]
            if get_std:
                ret.append(float(line[idx + 2]))
        else:
            if get_std:
                raise ResultsError("Did not find +/- in {:}. Does the value you requested have one?!".format(file))
            # extract first float
            for s in line:
                try:
                    ret = [float(s)]
                    break
                except BaseException:
                    continue
        # unit always in [], sometimes at the end, sometimes in the middle.
        # always use the first [] to find the unit (sometimes it is given twice
        # and the second one is empty)
        if get_unit:
            line = " ".join(line)
            idx1 = line.index("[")
            idx2 = line.index("]")
            ret.append(line[idx1 + 1 : idx2])

        if len(ret) > 1:
            return tuple(ret)
        else:
            return ret[0]

    def get_from_all_files(self, *args, output_folder="Output/System_0/", method="get_value", **kwargs):
        """Wrapper to Execute Methods on all Output files.

        output_folder: string
        The subfolder containing the relevant output files. Default should be fine.

        method: function
        The function to call, needs to be a function of the |RaspaResults| class that takes a *file* argument.
        All *args* and *kwargs* are passed on.

        If the *method* returns tuples, the tuples are unpacked and added to lists.
        """
        call_method = getattr(self, method)
        relevant_files = [f for f in self.files if output_folder in f]
        ret = []
        for f in relevant_files:
            data = call_method(*args, file=f, **kwargs)
            if isinstance(data, tuple):
                # initialize ret list to proper length
                if len(ret) == 0:
                    ret = [[] for i in range(len(data))]
                for i, d in enumerate(data):
                    ret[i].append(d)
            else:
                ret.append(data)

        if isinstance(ret[0], list):
            return tuple(ret)
        else:
            return ret

    def get_isotherm(
        self,
        output_folder="Output/System_0/",
        search_x="Partial pressure:",
        search_y=r"Average loading excess \[cm^3 (STP)/gr",
        get_std=False,
        get_unit=False,
    ):
        """Try to automatically collect Isotherm Data.

        output_folder: string
        The subfolder containing the relevant output files. Default should be fine.

        search_x: string
        String labeling the desired x-Axis values in the output.

        search_y: string
        String labeling the desired y-Axis values in the output.

        get_std: boolean
        Return the standard deviation of the y values as a list.

        get_unit: boolean
        Return the units of the values as strings `x_unit`, `y_unit`.

        Returns two lists: `x_values`, `y_values`. Optionally after that `x_std`, `x_unit`, `y_unit`.
        """
        ret = []
        x_info = self.get_from_all_files(search_x, output_folder=output_folder, get_unit=get_unit)
        y_info = self.get_from_all_files(search_y, output_folder=output_folder, get_std=get_std, get_unit=get_unit)

        if not get_std and not get_unit:
            ret = [x_info, y_info]
        elif get_std and not get_unit:
            ret = [x_info, *y_info]
        elif get_std:
            ret = [x_info[0], *y_info[0:2], x_info[1][0], y_info[2][0]]  # only one unit
        else:
            ret = [x_info[0], y_info[0], x_info[1][0], y_info[1][0]]  # only one unit

        return tuple(ret)


class RaspaJob(SingleJob):
    """A single job with RASPA2 class.

    Requires the ASE Python package.

    Files that should be copied to the job directory can be passed as a dictionary using the ``copy``
    argument. The dict should have the following layout: *dict[<filename in jobdir>] = <path to file>*.
    If you prefer symlinks instead of copies use the ``symlink`` argument with the same layout.
    The path strings given to symlink are not altered. Remember that the job directory is not equal
    to the current work directory if using relative symlinks.

    **Molecule parsing is not yet supported!**

    The environment variable ``$RASPA_DIR`` needs to be set for the runscript to work.
    """

    _result_type = RaspaResults
    _filenames = {"inp": "simulation.input", "run": "$JN.run", "err": "$JN.err", "stdout": "$JN.out"}

    def __init__(self, copy=None, symlink=None, **kwargs):
        SingleJob.__init__(self, **kwargs)
        self.copy_files = copy
        self.symlink_files = symlink
        # Output Name is not known before running
        self.settings.runscript.stdout_redirect = True

    def _get_ready(self):
        """Copy files to execution dir if self.copy_files is set."""
        SingleJob._get_ready(self)
        if self.copy_files:
            for filename, path in self.copy_files.items():
                dest = opj(self.path, filename)
                shutil.copy(path, dest)
        if self.symlink_files:
            for filename, path in self.symlink_files.items():
                dest = opj(self.path, filename)
                symlink(path, dest)
        return

    def get_input(self):
        """Transform all contents of ``settings.input`` branch into string with blocks, keys and values.

        RASPA2 does not interpret indentation but rather certain keys automatically start a subsection.
        Subsections are never closed but *automatically* ended if the next key-value pair does not fit into it.
        So here we only need to care about readability of the input file. The logic is automatically preserved
        by using the nested PLAMS ``settings`` logic.

        Components should be built using either ``component['INT'] = Settings()`` or ``component._INT = Settings()``.
        It needs to have a key ``moleculename`` (case insensitive).
        """

        def parse(key, value, indent=""):
            # transform bools to yes and no
            if value is True:
                value = "yes"
            elif value is False:
                value = "no"

            ret = indent + key

            if key.lower() == "component":
                n_comp = len(value)
                for i in range(n_comp):
                    # get component number
                    if str(i) in value:
                        number = str(i)
                    elif "_{}".format(i) in value:
                        number = "_{}".format(i)
                    else:
                        raise JobError(
                            "Check that you have a valid component {:} in settings.input.component".format(i)
                        )
                    ret += " " + str(i)

                    subkeys = list(value[number].keys())
                    mol_name = subkeys[[sk.lower() for sk in subkeys].index("moleculename")]
                    ret += " MoleculeName " + str(value[number][mol_name])
                    ret += "\n"
                    tmp_value = value[number].copy()
                    del tmp_value[mol_name]  # exclude MoleculeName
                    for el in tmp_value:
                        ret += parse(el, tmp_value[el], indent + "  ")
                    ret += "\n" + indent + key

                return ret

            elif isinstance(value, Settings):
                if "_h" in value:
                    ret += " " + str(value["_h"])
                ret += "\n"

                for el in value:
                    if el != "_h":
                        ret += parse(el, value[el], indent + "  ")
            else:
                ret += " " + str(value)
            ret += "\n"
            return ret

        inp = ""
        for item in self.settings.input:
            inp += parse(item, self.settings.input[item])

        return inp

    def get_runscript(self):
        """$RASPA_DIR has to be set in the environment!"""
        ret = "export DYLD_LIBRARY_PATH=${RASPA_DIR}/lib"
        ret += "\nexport LD_LIBRARY_PATH=${RASPA_DIR}/lib"
        ret += "\n$RASPA_DIR/bin/simulate "
        if self.settings.runscript.stdout_redirect:
            ret += " >" + self._filename("stdout")
        ret += "\n\n"
        return ret

    def check(self):
        """Look for the normal termination signal in output."""
        s = self.results.grep_output("Simulation finished on")
        return bool(s)

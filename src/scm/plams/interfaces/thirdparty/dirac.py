import os
from os.path import join as opj

from scm.plams.core.basejob import SingleJob
from scm.plams.core.private import saferun
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings

__all__ = ["DiracJob", "DiracResults"]


class DiracResults(Results):
    """A class for result of computation done with DIRAC."""

    _rename_map = {
        "DFCOEF": "$JN.dfcoef",
        "GRIDOUT": "$JN.grid",
        "dirac.xml": "$JN.xml",
    }

    def collect(self):
        """After collecting the files produced by job execution with parent method :meth:`Results.collect<scm.plams.core.results.Results.collect>` append the ``pam`` output to the regular output file."""
        Results.collect(self)
        pamfile = self.job._filename("out")
        process = saferun(["grep", "output file", pamfile], cwd=self.job.path)
        output = process.stdout.decode()
        diracfile = output.split(":")[-1].strip()
        if diracfile in self.files:
            pampath = opj(self.job.path, pamfile)
            diracpath = opj(self.job.path, diracfile)
            with open(pampath, "r") as f:
                pamoutput = f.readlines()
            with open(diracpath, "r") as f:
                diracoutput = f.readlines()
            with open(pampath, "w") as f:
                f.writelines(diracoutput)
                f.write("\n\n   " + "*" * 74 + "\n")
                f.write("   " + "*" * 30 + "  pam output  " + "*" * 30 + "\n")
                f.write("   " + "*" * 74 + "\n\n")
                f.writelines(pamoutput)
            os.remove(diracpath)
        self.refresh()


class DiracJob(SingleJob):
    """A class representing a single computational job with DIRAC."""

    _result_type = DiracResults
    _top = ["dirac"]
    _filenames = {"inp": "$JN.inp", "run": "$JN.run", "out": "$JN.out", "err": "$JN.err"}

    def __init__(self, **kwargs):
        SingleJob.__init__(self, **kwargs)
        self.settings.runscript.pam.noarch = True
        self.settings.runscript.pam.get = ["DFCOEF", "GRIDOUT", "dirac.xml"]

    def _get_ready(self):
        """Before generating runscript and input with parent method :meth:`SingleJob._get_ready<scm.plams.core.basejob.SingleJob._get_ready>` add proper ``mol`` and ``inp`` entries to ``self.settings.runscript.pam``. If already present there, ``mol`` will not be added."""
        s = self.settings.runscript.pam
        if "mol" not in s:
            s.mol = self.name + ".xyz"
            with open(opj(self.path, self.name + ".xyz"), "w") as f:
                f.write(str(len(self.molecule)) + "\n\n")
                for atom in self.molecule:
                    suffix = "b={block}" if hasattr(atom, "block") else ""
                    f.write(atom.str(suffix=suffix) + "\n")
        s.inp = self._filename("inp")
        SingleJob._get_ready(self)

    def get_input(self):
        """Transform all contents of ``input`` branch of ``settings`` into string with blocks, subblocks, keys and values.

        On the highest level alphabetic order of iteration is modified: keys occuring in class attribute ``_top`` are printed first. See :ref:`dirac-input` for details.
        """
        is_empty = lambda x: isinstance(x, Settings) and len(x) == 0

        def parse_key(key, value):
            ret = "." + key.upper() + "\n"
            if not (value is True or is_empty(value)):
                if isinstance(value, list):
                    for i in value:
                        ret += str(i) + "\n"
                else:
                    ret += str(value) + "\n"
            return ret

        def parse_block(block):
            enabler = "_en"
            ret = "**" + block.upper() + "\n"
            s = self.settings.input[block]
            for k, v in s.items():
                if not isinstance(v, Settings) or is_empty(v):
                    ret += parse_key(k, v)
            for k, v in s.items():
                if isinstance(v, Settings) and enabler in v:
                    ret += parse_key(k, v[enabler])
            for k, v in s.items():
                if isinstance(v, Settings) and len(v) > 0:
                    ret += "*" + k.upper() + "\n"
                    for kk, vv in v.items():
                        if kk != enabler:
                            ret += parse_key(kk, vv)
            return ret

        inp = ""
        for block in self._top:
            if block in self.settings.input:
                inp += parse_block(block)
        for block in self.settings.input:
            if block not in self._top:
                inp += parse_block(block)
        inp += "*END OF INPUT\n"
        return inp

    def get_runscript(self):
        """Generate a runscript. Returned string is a ``pam`` call followed by option flags generated based on ``self.settings.runscript.pam`` contents. See :ref:`dirac-runscript` for details."""
        r = self.settings.runscript.pam
        ret = "pam"
        for k, v in r.items():
            ret += " --%s" % k
            if v is not True:
                if isinstance(v, list):
                    ret += '="%s"' % " ".join(v)
                else:
                    ret += "=" + str(v)

        if self.settings.runscript.stdout_redirect:
            ret += " >" + self._filename("out")
        ret += "\n\n"
        return ret

    def check(self):
        """Check if the calculation was successful by examining the last line of ``pam`` output."""
        s = self.results.grep_output("exit           :")[0]
        status = s.split(":")[-1].strip()
        return status == "normal"

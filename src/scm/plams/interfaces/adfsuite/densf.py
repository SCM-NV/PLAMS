from scm.plams.core.errors import PlamsError
from scm.plams.interfaces.adfsuite.scmjob import SCMJob, SCMResults

__all__ = ["DensfJob", "DensfResults"]


class DensfResults(SCMResults):
    _kfext = ".t41"
    _rename_map = {"TAPE41": "$JN" + _kfext}

    def get_molecule(self, *args, **kwargs):
        raise PlamsError("DensfResults do not support get_molecule() method. You can get molecule from inputjob")


class DensfJob(SCMJob):
    """A class representing calculation of molecular properties on a grid using ``densf`` program.

    A new attribute ``inputjob`` is introduced to supply KF file from previously run job. The value can either be a string with a path to KF file or an instance of any type of |SCMJob| or |SCMResults| (in this case the path to corresponding KF file will be extracted automatically). If the value of ``inputjob`` is ``None``, no automatic handling occurs and user needs to manually supply path to input job using ``INPUTFILE`` keyword placed in ``myjob.settings.input``.

    The resulting ``TAPE41`` file is renamed to ``jobname.t41``.
    """

    _result_type = DensfResults
    _command = "densf"
    _top = ["inputfile", "units"]

    def __init__(self, inputjob=None, **kwargs):
        SCMJob.__init__(self, **kwargs)
        self.inputjob = inputjob

    def _serialize_mol(self):
        self.settings.input.inputfile = self.inputjob

    def _remove_mol(self):
        if "inputfile" in self.settings.input:
            del self.settings.input.inputfile

    def check(self):
        try:
            grep = self.results.grep_file("$JN.err", "NORMAL TERMINATION")
        except:
            return False
        return len(grep) > 0

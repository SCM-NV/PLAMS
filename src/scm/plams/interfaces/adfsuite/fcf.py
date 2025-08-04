from scm.plams.core.errors import PlamsError
from scm.plams.interfaces.adfsuite.scmjob import SCMJob, SCMResults

__all__ = ["FCFJob", "FCFResults"]


class FCFResults(SCMResults):
    _kfext = ".t61"
    _rename_map = {"TAPE61": "$JN" + _kfext}

    def get_molecule(self, *args, **kwargs):
        raise PlamsError("FCFResults do not support get_molecule() method. You can get molecules from job1 or job2")


class FCFJob(SCMJob):
    """A class representing calculation of Franck-Condon factors using ``fcf`` program.

    Two new attributes are introduced: ``inputjob1`` and ``inputjob2``. They are used to supply KF files from previous runs to ``fcf`` program. The value can either be a string with a path to KF file or an instance of any type of |SCMJob| or |SCMResults| (in this case the path to corresponding KF file will be extracted automatically). If the value of ``inputjob1`` or ``inputjob2`` is ``None``, no automatic handling occurs and user needs to manually supply paths to input jobs using proper keywords placed in ``myjob.settings.input`` (``STATES`` or ``STATE1`` and ``STATE2``).


    The resulting ``TAPE61`` file is renamed to ``jobname.t61``.
    """

    _result_type = FCFResults
    _command = "fcf"
    _top = ["states", "state1", "state2"]

    def __init__(self, inputjob1=None, inputjob2=None, **kwargs):
        SCMJob.__init__(self, **kwargs)
        self.inputjob1 = inputjob1
        self.inputjob2 = inputjob2

    def _serialize_mol(self):
        self.settings.input.state1 = self.inputjob1
        self.settings.input.state2 = self.inputjob2

    def _remove_mol(self):
        if "state1" in self.settings.input:
            del self.settings.input.state1
        if "state2" in self.settings.input:
            del self.settings.input.state2

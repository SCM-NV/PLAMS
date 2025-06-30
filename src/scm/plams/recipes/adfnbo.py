import os

from scm.plams.core.functions import log
from scm.plams.interfaces.adfsuite.ams import AMSJob

__all__ = ["ADFNBOJob"]


class ADFNBOJob(AMSJob):

    def prerun(self):  # noqa F811
        s = self.settings.input.ADF
        s.fullfock = True
        s.aomat2file = True
        s.symmetry = "NoSym"
        s.basis.core = "None"
        if "save" in s:
            if isinstance(s.save, str):
                s.save += " TAPE15"
            elif isinstance(s.save, list):
                s.save.append("TAPE15")
            else:
                log(
                    "WARNING: 'SAVE TAPE15' could not be added to the input settings of {}. Make sure (thisjob).settings.input.save is a string or a list.".format(
                        self.name
                    ),
                    1,
                )
        else:
            s.save = "TAPE15"

        if isinstance(self.settings.adfnbo, list):
            adfnbo_input = self.settings.adfnbo
        else:
            adfnbo_input = ["write", "spherical", "fock"]
            log(
                "WARNING: (thisjob).settings.adfnbo should be a list. Using default settings: write, fock, spherical", 1
            )

        self.settings.runscript.post = (
            'cp "' + os.path.join(self.path, "adf.rkf") + '" TAPE21\n'
            "$AMSBIN/adfnbo <<eor\n" + "\n".join(adfnbo_input) + "\neor\n\n$AMSBIN/gennbo6 FILE47\n"
        )

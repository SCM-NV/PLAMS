from functools import reduce

from scm.plams.core.basejob import SingleJob
from scm.plams.core.settings import Settings

__all__ = ["GamessJob"]


# ======================<>===========================


class GamessJob(SingleJob):
    """
    A class representing a single computational job with
    `GAMESS-US <http://www.msg.ameslab.gov/gamess/>`.
    A typical gamess input looks like:

      $contrl scftyp=rohf mult=2 cctyp=eom-ccsd runtyp=energy ispher=1 $end
      $system mwords=3 $end
      $basis  gbasis=accd $end
      $guess  guess=huckel $end
      $ccinp  ccprp=.false. $end
      $eominp nstate(1)=1,0,0,0 iroot(1)=1,1 ccprpe=.false.
              minit=1 noact=3 nuact=5 $end
      $data
     NH2...aug-cc-pVDZ basis...excited state @ g.s. MP2 geometry
     Cnv 2

     NITROGEN    7.0   0.0   0.0           -0.0370366120
     HYDROGEN    1.0   0.0  -0.8055298238  -0.6814816939
      $end
    Using Plams Settings it is written like:

    j = GamessJob()

    # data
    header = "NH2...aug-cc-pVDZ basis...excited state @ g.s. MP2 geometry\n"
    symmetry = "Cnv 2\n\n"
    n = "NITROGEN    7.0   0.0   0.0           -0.0370366120\n"
    h = "HYDROGEN    1.0   0.0  -0.8055298238  -0.6814816939"

    j.settings.input.data = header + symmetry + n + h
    # basis
    j.settings.input.basis.gbasis = 'accd'
    # guess
    j.settings.input.guess.guess = 'huckel'
    # ccinp
    j.settings.input.ccinp.ccprp = '.false.'
    # eominp
    j.settings.input.eominp['nstate(1)'] = '1,0,0,0'
    j.settings.input.eominp['iroot(1)'] = '1,1'
    j.settings.input.eominp.ccprpe = '.false.'
    j.settings.input.eominp.minit = 1
    j.settings.input.eominp.noact = 3
    j.settings.input.eominp.nuact = 5
    # control
    j.settings.input.contrl.scftyp = 'rhf'
    j.settings.input.contrl.runtype = 'optimize'
    j.settings.input.contrl.pp = 'mcp'
    j.settings.input.contrl.ispher = 1
    j.settings.input.contrl.coord = 'zmt'
    j.settings.input.contrl.nzvar = 6
    """

    _filenames = {"inp": "$JN.inp", "run": " $JN.run", "out": "$JN.out", "err": "$JN.err"}

    def get_input(self):
        """
        Transform all contents of ``input`` branch of ``settings`` into string
        with blocks, subblocks, keys and values.
        """

        def parse(key, value, indent=""):
            ret = ""
            if isinstance(value, Settings):
                ret += " ${} ".format(key)
                for el in value:
                    ret += " {}={} ".format(el, value[el])
                ret += "$end\n"
            else:
                ret += " ${}\n{}\n $end\n".format(key, value)

            return ret

        inp = [parse(item, self.settings.input[item]) for item in self.settings.input]

        return self.print_molecule() + reduce(lambda x, y: x + y, inp)

    def print_molecule(self):
        """
        pretty print a molecule in the GAMESS format.
        """
        mol = self.molecule
        if mol:
            # Read Symmetry from properties otherwise use C1
            if "symmetry" in mol.properties.keys():
                sym = mol.properties["symmetry"]
            else:
                sym = "C1"
            if sym != "C1":
                sym += "\n"
            ret = " $data\ntitle\n{}\n".format(sym)
            for at in mol.atoms:
                ret += "{} {}   {}\n".format(at.symbol, at.atnum, at.str(symbol=False, space=14, decimal=10))
            ret += " $end\n"
            return ret
        else:
            return "\n"

    def get_runscript(self):
        """
        Run Gamess Using rungms.
        """
        return "rungms {} > {}".format(self._filename("inp"), self._filename("out"))

    def check(self):
        """
        Look for the normal termination signal in Cp2k output
        """
        s = self.results.grep_output("EXECUTION OF GAMESS TERMINATED NORMALLY")
        return len(s) > 0

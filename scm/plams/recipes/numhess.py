from collections import OrderedDict
from itertools import product

from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
import numpy as np
from scm.plams.tools.units import Units

__all__ = ["NumHessJob", "NumHessResults"]  # names exported to the main namespace


class NumHessResults(Results):

    def get_hessian(self, mass_weighted=False):
        j = self.job
        n = len(j.molecule)
        hessian = np.empty((3 * n, 3 * n))

        for atom, axis in product(range(1, n + 1), range(3)):
            v1 = np.array(j.gradient(j.children[(atom, axis, -1)].results))
            v2 = np.array(j.gradient(j.children[(atom, axis, 1)].results))
            hessian[3 * (atom - 1) + axis] = (v2 - v1) / (2 * Units.convert(j.step, "angstrom", "bohr"))

        if mass_weighted:
            masses = [at.mass for at in j.molecule]
            weights = np.empty((n, n))
            for i, j in product(range(n), repeat=2):
                weights[i, j] = masses[i] * masses[j]
            hessian /= np.kron(weights, np.ones((3, 3)))
        return (hessian + hessian.T) / 2


class NumHessJob(MultiJob):
    """A class for calculating numerical hessian.

    The length and unit of the geometry displacement step can be adjusted.
    *gradient* should be a function that takes results of *jobtype* and
    returns energy gradient in hartree/bohr.
    """

    _result_type = NumHessResults

    def __init__(self, molecule, step=0.01, unit="angstrom", jobtype=None, gradient=None, **kwargs):
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)
        self.molecule = molecule
        self.step = Units.convert(step, unit, "angstrom")
        self.jobtype = jobtype  # who is going to calculate single points
        self.gradient = gradient  # function extracting gradients from children

    def prerun(self):  # noqa F811
        for atom, axis, step in product(range(1, 1 + len(self.molecule)), range(3), [-1, 1]):
            vec = [0, 0, 0]
            vec[axis] = self.step * step
            newmol = self.molecule.copy()
            newmol[atom].translate(vec)
            newname = "{}_{}_{}".format(atom, axis, step)
            self.children[(atom, axis, step)] = self.jobtype(name=newname, molecule=newmol, settings=self.settings)

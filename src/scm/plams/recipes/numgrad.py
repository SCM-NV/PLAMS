from collections import OrderedDict
from itertools import product

from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results

__all__ = ["NumGradJob", "NumGradResults"]  # names exported to the main namespace


class NumGradResults(Results):

    def _get_gradient_component(self, atom, coord, extract):
        """Get gradient component for *atom* along *coord*.

        *atom* should be integer, *coord* a single character string with 'x', 'y' or 'z'.
        *extract* should be a function that takes results of a single point job and
        returns a single number.
        """
        energies = [extract(self.job.children[(atom, coord, i)].results) for i in self.job.steps]
        coeffs, denom = self.job._coeffs[self.job.npoints]
        return sum([c * e for c, e in zip(coeffs, energies)]) / (denom * self.job.step)

    def get_gradient(self, extract):
        """Get the full gradient vector. Returns a list of length 3N,
        where N is the number of atoms in the calculated molecule.

        *extract* should be a function that takes results of a single point job and
        returns a single number.
        """
        return [
            self._get_gradient_component(atom, coord, extract)
            for atom, coord in product(range(1, 1 + len(self.job.molecule)), "xyz")
        ]


class NumGradJob(MultiJob):
    """A class for calculating numerical gradients of energy
    (or some other single-valued property) with respect to
    cartesian displacements.

    The differentiation is done using the final difference method
    with 2, 4, 6 or 8 points. The length of the step can be adjusted.
    """

    _result_type = NumGradResults

    # finite difference coefficients for different number of points:
    _coeffs = {
        2: ([-1, 1], 2),
        4: ([1, -8, 8, -1], 12),
        6: ([-1, 9, -75, 75, -9, 1], 60),
        8: ([3, -32, 168, -672, 672, -168, 32, -3], 840),
    }

    def __init__(self, molecule, npoints=2, step=0.01, unit="angstrom", jobtype=None, **kwargs):
        # initialize parent class and store all additional constructor arguments
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)
        self.molecule = molecule
        self.npoints = npoints
        self.step = step
        self.unit = unit
        self.jobtype = jobtype  # who is going to calculate single points

    def prerun(self):  # noqa F811
        self.steps = list(range(-(self.npoints // 2), 0)) + list(range(1, self.npoints // 2 + 1))

        for atom, axis, i in product(range(1, 1 + len(self.molecule)), "xyz", self.steps):
            v = (
                self.step * i if axis == "x" else 0,
                self.step * i if axis == "y" else 0,
                self.step * i if axis == "z" else 0,
            )
            newmol = self.molecule.copy()
            newmol[atom].translate(v, self.unit)
            newname = "{}_{}_{}".format(atom, axis, i)
            # settings of NumGradJob are directly passed to children single point jobs
            self.children[(atom, axis, i)] = self.jobtype(name=newname, molecule=newmol, settings=self.settings)

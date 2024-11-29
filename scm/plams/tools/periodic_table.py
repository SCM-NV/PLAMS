from scm.plams.core.errors import PTError
from typing import List
import numpy

__all__ = ["PeriodicTable", "PT"]


class PeriodicTable:
    """A singleton class for the periodic table of elements.

    For each element the following properties are stored: atomic symbol, atomic mass, atomic radius and number of connectors.

    Atomic mass is, strictly speaking, atomic weight, as present in Mathematica's ElementData function.

    Atomic radius and number of connectors are used by :meth:`~scm.plams.mol.molecule.Molecule.guess_bonds`. Note that values of radii are neither atomic radii nor covalent radii. They are somewhat "emprically optimized" for the bond guessing algorithm.

    .. note::

        This class is visible in the main namespace as both ``PeriodicTable`` and ``PT``.
    """

    # Note: is_metallic and is_electronegative give "pragmatic" values for these quantities for the benefit of the
    #       bond guessing algorithm.
    data: List[List] = [
        # [symbol, mass, radius, connectors, is_metallic, is_electronegative]
        # atomic weights from: http://www.ciaaw.org/atomic-weights.htm
        ["Xx", 0.00000, 0.00, 0, 0, 0],
        ["H", 1.00798, 0.30, 1, 0, 0],
        ["He", 4.00260, 0.99, 0, 0, 0],
        ["Li", 6.96750, 1.52, 8, 1, 0],
        ["Be", 9.01218, 1.12, 8, 1, 0],
        ["B", 10.81350, 0.88, 6, 0, 0],
        ["C", 12.01060, 0.77, 4, 0, 1],
        ["N", 14.00685, 0.70, 3, 0, 1],
        ["O", 15.99940, 0.66, 2, 0, 1],
        ["F", 18.99840, 0.64, 1, 0, 1],
        ["Ne", 20.17970, 1.60, 0, 0, 0],
        ["Na", 22.98977, 1.86, 8, 1, 0],
        ["Mg", 24.30550, 1.60, 8, 1, 0],
        ["Al", 26.98154, 1.43, 8, 1, 0],
        ["Si", 28.08500, 1.17, 8, 1, 0],
        ["P", 30.97376, 1.10, 8, 0, 1],
        ["S", 32.06750, 1.04, 2, 0, 1],
        ["Cl", 35.45150, 0.99, 1, 0, 1],
        ["Ar", 39.94800, 1.92, 0, 0, 0],
        ["K", 39.09830, 2.31, 8, 0, 0],
        ["Ca", 40.07800, 1.97, 8, 1, 0],
        ["Sc", 44.95591, 1.60, 8, 1, 0],
        ["Ti", 47.86700, 1.46, 8, 1, 0],
        ["V", 50.94150, 1.31, 8, 1, 0],
        ["Cr", 51.99610, 1.25, 8, 1, 0],
        ["Mn", 54.93804, 1.29, 8, 1, 0],
        ["Fe", 55.84500, 1.26, 8, 1, 0],
        ["Co", 58.93319, 1.25, 8, 1, 0],
        ["Ni", 58.69340, 1.24, 8, 1, 0],
        ["Cu", 63.54600, 1.28, 8, 1, 0],
        ["Zn", 65.38000, 1.33, 8, 1, 0],
        ["Ga", 69.72300, 1.41, 8, 1, 0],
        ["Ge", 72.63000, 1.22, 8, 1, 0],
        ["As", 74.92159, 1.21, 8, 0, 0],
        ["Se", 78.97100, 1.17, 8, 0, 1],
        ["Br", 79.90400, 1.14, 1, 0, 1],
        ["Kr", 83.79800, 1.97, 0, 0, 0],
        ["Rb", 85.46780, 2.44, 8, 1, 0],
        ["Sr", 87.62000, 2.15, 8, 1, 0],
        ["Y", 88.90584, 1.80, 8, 1, 0],
        ["Zr", 91.22400, 1.57, 8, 1, 0],
        ["Nb", 92.90637, 1.41, 8, 1, 0],
        ["Mo", 95.95000, 1.36, 8, 1, 0],
        ["Tc", 98.00000, 1.35, 8, 1, 0],
        ["Ru", 101.07000, 1.33, 8, 1, 0],
        ["Rh", 102.90550, 1.34, 8, 1, 0],
        ["Pd", 106.42000, 1.38, 8, 1, 0],
        ["Ag", 107.86820, 1.44, 8, 1, 0],
        ["Cd", 112.41400, 1.49, 8, 1, 0],
        ["In", 114.81800, 1.66, 8, 1, 0],
        ["Sn", 118.71000, 1.62, 8, 1, 0],
        ["Sb", 121.76000, 1.41, 8, 1, 0],
        ["Te", 127.60000, 1.37, 8, 0, 1],
        ["I", 126.90447, 1.33, 1, 0, 1],
        ["Xe", 131.29300, 2.17, 0, 0, 0],
        ["Cs", 132.90545, 2.62, 8, 1, 0],
        ["Ba", 137.32700, 2.17, 8, 1, 0],
        ["La", 138.90547, 1.88, 8, 1, 0],
        ["Ce", 140.11600, 1.818, 8, 1, 0],
        ["Pr", 140.90766, 1.824, 8, 1, 0],
        ["Nd", 144.24200, 1.814, 8, 1, 0],
        ["Pm", 145.00000, 1.834, 8, 1, 0],
        ["Sm", 150.36000, 1.804, 8, 1, 0],
        ["Eu", 151.96400, 2.084, 8, 1, 0],
        ["Gd", 157.25000, 1.804, 8, 1, 0],
        ["Tb", 158.92535, 1.773, 8, 1, 0],
        ["Dy", 162.50000, 1.781, 8, 1, 0],
        ["Ho", 164.93033, 1.762, 8, 1, 0],
        ["Er", 167.25900, 1.761, 8, 1, 0],
        ["Tm", 168.93422, 1.759, 8, 1, 0],
        ["Yb", 173.04500, 1.922, 8, 1, 0],
        ["Lu", 174.96680, 1.738, 8, 1, 0],
        ["Hf", 178.49000, 1.57, 8, 1, 0],
        ["Ta", 180.94788, 1.43, 8, 1, 0],
        ["W", 183.84000, 1.37, 8, 1, 0],
        ["Re", 186.20700, 1.37, 8, 1, 0],
        ["Os", 190.23000, 1.34, 8, 1, 0],
        ["Ir", 192.21700, 1.35, 8, 1, 0],
        ["Pt", 195.08400, 1.38, 8, 1, 0],
        ["Au", 196.96657, 1.44, 8, 1, 0],
        ["Hg", 200.59200, 1.52, 8, 1, 0],
        ["Tl", 204.38350, 1.71, 8, 1, 0],
        ["Pb", 207.20000, 1.75, 8, 1, 0],
        ["Bi", 208.98040, 1.70, 8, 1, 0],
        ["Po", 209.00000, 1.40, 8, 1, 0],
        ["At", 210.00000, 1.40, 1, 0, 1],
        ["Rn", 222.00000, 2.40, 0, 0, 0],
        ["Fr", 223.00000, 2.70, 8, 1, 0],
        ["Ra", 226.00000, 2.20, 8, 1, 0],
        ["Ac", 227.00000, 2.00, 8, 1, 0],
        ["Th", 232.03770, 1.79, 8, 1, 0],
        ["Pa", 231.03588, 1.63, 8, 1, 0],
        ["U", 238.02891, 1.56, 8, 1, 0],
        ["Np", 237.00000, 1.55, 8, 1, 0],
        ["Pu", 244.00000, 1.59, 8, 1, 0],
        ["Am", 243.00000, 1.73, 8, 1, 0],
        ["Cm", 247.00000, 1.74, 8, 1, 0],
        ["Bk", 247.00000, 1.70, 8, 1, 0],
        ["Cf", 251.00000, 1.86, 8, 1, 0],
        ["Es", 252.00000, 1.86, 8, 1, 0],
        ["Fm", 257.00000, 2.00, 8, 1, 0],
        ["Md", 258.00000, 2.00, 8, 1, 0],
        ["No", 259.00000, 2.00, 8, 1, 0],
        ["Lr", 266.00000, 2.00, 8, 1, 0],
        ["Rf", 267.00000, 2.00, 8, 1, 0],
        ["Db", 268.00000, 2.00, 8, 1, 0],
        ["Sg", 269.00000, 2.00, 8, 1, 0],
        ["Bh", 270.00000, 2.00, 8, 1, 0],
        ["Hs", 277.00000, 2.00, 8, 1, 0],
        ["Mt", 278.00000, 2.00, 8, 1, 0],
        ["Ds", 281.00000, 2.00, 8, 1, 0],
        ["Rg", 282.00000, 2.00, 8, 1, 0],
        ["Cn", 285.00000, 2.00, 8, 1, 0],
        ["Nh", 286.00000, 2.00, 8, 1, 0],
        ["Fl", 289.00000, 2.00, 8, 1, 0],
        ["Mc", 290.00000, 2.00, 8, 1, 0],
        ["Lv", 293.00000, 2.00, 8, 1, 0],
        ["Ts", 294.00000, 2.00, 8, 1, 0],
        ["Og", 294.00000, 2.00, 8, 1, 0],
    ]

    symtonum = {d[0]: i for i, d in enumerate(data)}

    # Collection of symbols used for different kinds of dummy atoms:
    dummysymbols = ["Xx", "El", "Eh", "J"]

    def __init__(self):
        raise PTError("Instances of PeriodicTable cannot be created")

    @classmethod
    def get_atomic_number(cls, symbol):
        """Convert atomic symbol to atomic number."""
        if symbol.lower().capitalize() in cls.dummysymbols:
            return 0
        try:
            number = cls.symtonum[symbol.capitalize()]
        except KeyError:
            raise PTError("trying to convert incorrect atomic symbol")
        return number

    @classmethod
    def get_symbol(cls, atnum):
        """Convert atomic number to atomic symbol."""
        try:
            symbol = cls.data[atnum][0]
        except IndexError:
            raise PTError("trying to convert incorrect atomic number")
        return symbol

    @classmethod
    def get_mass(cls, arg):
        """Convert atomic symbol or atomic number to atomic mass."""
        if isinstance(arg, str) and arg.lower().capitalize() in ["El", "Eh"]:
            return cls.get_mass("H")
        else:
            return cls._get_property(arg, 1)

    @classmethod
    def get_radius(cls, arg):
        """Convert atomic symbol or atomic number to radius."""
        return cls._get_property(arg, 2)

    @classmethod
    def get_connectors(cls, arg):
        """Convert atomic symbol or atomic number to number of connectors."""
        return cls._get_property(arg, 3)

    @classmethod
    def get_metallic(cls, arg):
        """Convert atomic symbol or atomic number to number of connectors."""
        return cls._get_property(arg, 4)

    @classmethod
    def get_electronegative(cls, arg):
        """Convert atomic symbol or atomic number to number of connectors."""
        return cls._get_property(arg, 5)

    @classmethod
    def set_mass(cls, element, value):
        """Set the mass of *element* to *value*."""
        cls.data[cls.get_atomic_number(element)][1] = value

    @classmethod
    def set_radius(cls, element, value):
        """Set the radius of *element* to *value*."""
        cls.data[cls.get_atomic_number(element)][2] = value

    @classmethod
    def set_connectors(cls, element, value):
        """Set the number of connectors of *element* to *value*."""
        cls.data[cls.get_atomic_number(element)][3] = value

    @classmethod
    def _get_property(cls, arg, prop):
        """Get property of element described by either symbol or atomic number. Skeleton method for :meth:`get_radius`, :meth:`get_mass` and  :meth:`get_connectors`."""
        if isinstance(arg, str):
            pr = cls.data[cls.get_atomic_number(arg)][prop]
        elif isinstance(arg, (int, numpy.int64, numpy.int32)):
            try:
                pr = cls.data[arg][prop]
            except KeyError:
                raise PTError("trying to convert incorrect atomic number")
        return pr


PT = PeriodicTable

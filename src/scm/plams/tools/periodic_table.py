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
        # [symbol, mass, radius, connectors, is_metallic, is_electronegative, electron_affinity, ionization_potential]
        # atomic weights from: http://www.ciaaw.org/atomic-weights.htm
        # https://pubchem.ncbi.nlm.nih.gov/ptable/electron-affinity/
        # https://pubchem.ncbi.nlm.nih.gov/ptable/ionization-energy/
        ["Xx", 0.00000, 0.00, 0, 0, 0, None, None],
        ["H", 1.00798, 0.30, 1, 0, 0, None, None],
        ["He", 4.00260, 0.99, 0, 0, 0, None, None],
        ["Li", 6.96750, 1.52, 8, 1, 0, 0.0227110811, 0.1981523452],
        ["Be", 9.01218, 1.12, 8, 1, 0, None, None],
        ["B", 10.81350, 0.88, 6, 0, 0, None, None],
        ["C", 12.01060, 0.77, 4, 0, 1, None, None],
        ["N", 14.00685, 0.70, 3, 0, 1, None, None],
        ["O", 15.99940, 0.66, 2, 0, 1, None, None],
        ["F", 18.99840, 0.64, 1, 0, 1, None, None],
        ["Ne", 20.17970, 1.60, 0, 0, 0, None, None],
        ["Na", 22.98977, 1.86, 8, 1, 0, 0.0201386286, 0.1888547667],
        ["Mg", 24.30550, 1.60, 8, 1, 0, None, None],
        ["Al", 26.98154, 1.43, 8, 1, 0, 0.0162064511, 0.2199814425],
        ["Si", 28.08500, 1.17, 8, 1, 0, 0.0508978112, 0.2995804744],
        ["P", 30.97376, 1.10, 8, 0, 1, None, None],
        ["S", 32.06750, 1.04, 2, 0, 1, None, None],
        ["Cl", 35.45150, 0.99, 1, 0, 1, None, None],
        ["Ar", 39.94800, 1.92, 0, 0, 0, None, None],
        ["K", 39.09830, 2.31, 8, 0, 0, None, None],
        ["Ca", 40.07800, 1.97, 8, 1, 0, None, None],
        ["Sc", 44.95591, 1.60, 8, 1, 0, 0.0069088726, 0.2411123028],
        ["Ti", 47.86700, 1.46, 8, 1, 0, 0.0029031965, 0.2509243718],
        ["V", 50.94150, 1.31, 8, 1, 0, 0.0192933941, 0.2479109274],
        ["Cr", 51.99610, 1.25, 8, 1, 0, 0.0244750486, 0.2486826632],
        ["Mn", 54.93804, 1.29, 8, 1, 0, None, None],
        ["Fe", 55.84500, 1.26, 8, 1, 0, 0.0059901395, 0.2903931438],
        ["Co", 58.93319, 1.25, 8, 1, 0, 0.024291302, 0.2896214081],
        ["Ni", 58.69340, 1.24, 8, 1, 0, 0.0424822164, 0.2807648214],
        ["Cu", 63.54600, 1.28, 8, 1, 0, 0.0451281676, 0.2839252631],
        ["Zn", 65.38000, 1.33, 8, 1, 0, None, None],
        ["Ga", 69.72300, 1.41, 8, 1, 0, 0.0110247967, 0.2204591837],
        ["Ge", 72.63000, 1.22, 8, 1, 0, 0.0496115849, 0.2903196452],
        ["As", 74.92159, 1.21, 8, 0, 0, None, None],
        ["Se", 78.97100, 1.17, 8, 0, 1, None, None],
        ["Br", 79.90400, 1.14, 1, 0, 1, None, None],
        ["Kr", 83.79800, 1.97, 0, 0, 0, None, None],
        ["Rb", 85.46780, 2.44, 8, 1, 0, 0.0171986828, 0.1535019187],
        ["Sr", 87.62000, 2.15, 8, 1, 0, None, None],
        ["Y", 88.90584, 1.80, 8, 1, 0, 0.0112820419, 0.228470536],
        ["Zr", 91.22400, 1.57, 8, 1, 0, 0.0156552112, 0.2437950033],
        ["Nb", 92.90637, 1.41, 8, 1, 0, 0.0328171447, 0.2483886686],
        ["Mo", 95.95000, 1.36, 8, 1, 0, 0.0274149943, 0.2606261929],
        ["Tc", 98.00000, 1.35, 8, 1, 0, 0.0202121272, 0.2675350654],
        ["Ru", 101.07000, 1.33, 8, 1, 0, 0.0385867883, 0.2705117605],
        ["Rh", 102.90550, 1.34, 8, 1, 0, 0.0417839793, 0.2741131941],
        ["Pd", 106.42000, 1.38, 8, 1, 0, 0.0204693725, 0.306379099],
        ["Ag", 107.86820, 1.44, 8, 1, 0, 0.0478476175, 0.2784128648],
        ["Cd", 112.41400, 1.49, 8, 1, 0, None, None],
        ["In", 114.81800, 1.66, 8, 1, 0, 0.0110247967, 0.2126315781],
        ["Sn", 118.71000, 1.62, 8, 1, 0, 0.0440991866, 0.2698870221],
        ["Sb", 121.76000, 1.41, 8, 1, 0, 0.0393217747, 0.3175141436],
        ["Te", 127.60000, 1.37, 8, 0, 1, None, None],
        ["I", 126.90447, 1.33, 1, 0, 1, None, None],
        ["Xe", 131.29300, 2.17, 0, 0, 0, None, None],
        ["Cs", 132.90545, 2.62, 8, 1, 0, 0.0173456801, 0.1431018606],
        ["Ba", 137.32700, 2.17, 8, 1, 0, None, None],
        ["La", 138.90547, 1.88, 8, 1, 0, 0.0183746611, 0.2049509698],
        ["Ce", 140.11600, 1.82, 8, 1, 0, 0.0183746611, 0.2035544955],
        ["Pr", 140.90766, 1.82, 8, 1, 0, None, None],
        ["Nd", 144.24200, 1.81, 8, 1, 0, None, None],
        ["Pm", 145.00000, 1.83, 8, 1, 0, None, None],
        ["Sm", 150.36000, 1.80, 8, 1, 0, None, None],
        ["Eu", 151.96400, 2.08, 8, 1, 0, None, None],
        ["Gd", 157.25000, 1.80, 8, 1, 0, None, None],
        ["Tb", 158.92535, 1.77, 8, 1, 0, None, None],
        ["Dy", 162.50000, 1.78, 8, 1, 0, None, None],
        ["Ho", 164.93033, 1.76, 8, 1, 0, None, None],
        ["Er", 167.25900, 1.76, 8, 1, 0, None, None],
        ["Tm", 168.93422, 1.76, 8, 1, 0, None, None],
        ["Yb", 173.04500, 1.92, 8, 1, 0, None, None],
        ["Lu", 174.96680, 1.74, 8, 1, 0, None, None],
        ["Hf", 178.49000, 1.57, 8, 1, 0, None, None],
        ["Ta", 180.94788, 1.43, 8, 1, 0, 0.0118332817, 0.289952152],
        ["W", 183.84000, 1.37, 8, 1, 0, 0.0299506976, 0.293259591],
        ["Re", 186.20700, 1.37, 8, 1, 0, 0.0055123983, 0.2895846587],
        ["Os", 190.23000, 1.34, 8, 1, 0, 0.0404242544, 0.3197191029],
        ["Ir", 192.21700, 1.35, 8, 1, 0, 0.0575126892, 0.3344188318],
        ["Pt", 195.08400, 1.38, 8, 1, 0, None, None],
        ["Au", 196.96657, 1.44, 8, 1, 0, 0.0848541849, 0.3390492464],
        ["Hg", 200.59200, 1.52, 8, 1, 0, None, None],
        ["Tl", 204.38350, 1.71, 8, 1, 0, 0.0073498644, 0.2244648598],
        ["Pb", 207.20000, 1.75, 8, 1, 0, 0.013229756, 0.2725697226],
        ["Bi", 208.98040, 1.70, 8, 1, 0, 0.0347648588, 0.2678658093],
        ["Po", 209.00000, 1.40, 8, 1, 0, 0.0698237121, 0.3093190448],
        ["At", 210.00000, 1.40, 1, 0, 1, None, None],
        ["Rn", 222.00000, 2.40, 0, 0, 0, None, None],
        ["Fr", 223.00000, 2.70, 8, 1, 0, None, None],
        ["Ra", 226.00000, 2.20, 8, 1, 0, None, None],
        ["Ac", 227.00000, 2.00, 8, 1, 0, None, None],
        ["Th", 232.03770, 1.79, 8, 1, 0, None, None],
        ["Pa", 231.03588, 1.63, 8, 1, 0, None, None],
        ["U", 238.02891, 1.56, 8, 1, 0, None, None],
        ["Np", 237.00000, 1.55, 8, 1, 0, None, None],
        ["Pu", 244.00000, 1.59, 8, 1, 0, None, None],
        ["Am", 243.00000, 1.73, 8, 1, 0, None, None],
        ["Cm", 247.00000, 1.74, 8, 1, 0, None, None],
        ["Bk", 247.00000, 1.70, 8, 1, 0, None, None],
        ["Cf", 251.00000, 1.86, 8, 1, 0, None, None],
        ["Es", 252.00000, 1.86, 8, 1, 0, None, None],
        ["Fm", 257.00000, 2.00, 8, 1, 0, None, None],
        ["Md", 258.00000, 2.00, 8, 1, 0, None, None],
        ["No", 259.00000, 2.00, 8, 1, 0, None, None],
        ["Lr", 266.00000, 2.00, 8, 1, 0, None, None],
        ["Rf", 267.00000, 2.00, 8, 1, 0, None, None],
        ["Db", 268.00000, 2.00, 8, 1, 0, None, None],
        ["Sg", 269.00000, 2.00, 8, 1, 0, None, None],
        ["Bh", 270.00000, 2.00, 8, 1, 0, None, None],
        ["Hs", 277.00000, 2.00, 8, 1, 0, None, None],
        ["Mt", 278.00000, 2.00, 8, 1, 0, None, None],
        ["Ds", 281.00000, 2.00, 8, 1, 0, None, None],
        ["Rg", 282.00000, 2.00, 8, 1, 0, None, None],
        ["Cn", 285.00000, 2.00, 8, 1, 0, None, None],
        ["Nh", 286.00000, 2.00, 8, 1, 0, None, None],
        ["Fl", 289.00000, 2.00, 8, 1, 0, None, None],
        ["Mc", 290.00000, 2.00, 8, 1, 0, None, None],
        ["Lv", 293.00000, 2.00, 8, 1, 0, None, None],
        ["Ts", 294.00000, 2.00, 8, 1, 0, None, None],
        ["Og", 294.00000, 2.00, 8, 1, 0, None, None],
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

    @classmethod
    def get_electron_affinity(cls, arg):
        """
        Get the electron affinity of the metal element

        Note: Returns None if data is not available
        """
        return cls._get_property(arg, 6)

    @classmethod
    def get_ionization_energy(cls, arg):
        """
        Get the electron affinity of the metal element

        Note: Returns None if data is not available
        """
        return cls._get_property(arg, 7)


PT = PeriodicTable

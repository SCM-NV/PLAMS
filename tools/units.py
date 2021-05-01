import collections
import math
import numpy as np

from ..core.errors import UnitsError


__all__ = ['Units']


class Units:
    """A singleton class for unit converter.

    All values are based on `2014 CODATA recommended values <http://physics.nist.gov/cuu/Constants>`_.

    The following constants and units are supported:

    *   constants:

        -   ``speed_of_light`` (also ``c``)
        -   ``electron_charge`` (also ``e``)
        -   ``Avogadro_constant`` (also ``NA``)
        -   ``Bohr_radius``

    *   distance:

        -   ``Angstrom``, ``A``
        -   ``Bohr``, ``au``, ``a.u.``
        -   ``nm``
        -   ``pm``

    *   reciprocal distance:

        -   ``1/Angstrom``, ``1/A``, ``Angstrom^-1``, ``A^-1``,
        -   ``1/Bohr``, ``Bohr^-1``

    *   angle:

        -    ``degree``, ``deg``,
        -    ``radian``, ``rad``,
        -    ``grad``
        -    ``circle``

    *   energy:

        -   ``au``, ``a.u.``, ``Hartree``
        -   ``eV``
        -   ``kcal/mol``
        -   ``kJ/mol``
        -   ``cm^-1``, ``cm-1``
        -   ``K``, ``Kelvin``

    *   dipole moment:

        -   ``au``, ``a.u.``
        -   ``Cm``
        -   ``Debye``, ``D``

    Example::

        >>> print(Units.constants['speed_of_light'])
        299792458
        >>> print(Units.constants['e'])
        1.6021766208e-19
        >>> print(Units.convert(123, 'angstrom', 'bohr'))
        232.436313431
        >>> print(Units.convert([23.32, 145.0, -34.7], 'kJ/mol', 'kcal/mol'))
        [5.573613766730401, 34.655831739961755, -8.293499043977056]
        >>> print(Units.conversion_ratio('kcal/mol', 'kJ/mol'))
        4.184


    """

    constants = {}
    constants['Bohr_radius']                         =  0.529177210903  #http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
    constants['Avogadro_constant'] = constants['NA'] =  6.022140857e23   #http://physics.nist.gov/cgi-bin/cuu/Value?na
    constants['speed_of_light'] = constants['c']     =  299792458   #http://physics.nist.gov/cgi-bin/cuu/Value?c
    constants['electron_charge'] = constants['e']    =  1.6021766208e-19   #http://physics.nist.gov/cgi-bin/cuu/Value?e
    constants['Boltzmann'] = constants['k_B'] = 1.380649e-23 #J/K


    distance = {}
    distance['A'] = distance['Angstrom']                 =  1.0
    distance['Bohr'] = distance['a.u.'] = distance['au'] =  1.0 / constants['Bohr_radius']
    distance['nm']                                       = distance['A'] / 10.0
    distance['pm']                                       = distance['A'] * 100.0
    distance['m']                                        = distance['A'] * 1e-10

    rec_distance = {}
    rec_distance['1/A'] = rec_distance['1/Ang'] = rec_distance['1/Angstrom'] = rec_distance['A^-1'] = rec_distance['Ang^-1'] = rec_distance['Angstrom^-1'] = 1.0
    rec_distance['1/m'] = rec_distance['m^-1'] = 1e10
    rec_distance['1/Bohr'] = rec_distance['Bohr^-1'] = constants['Bohr_radius']

    energy = {}
    energy['au'] = energy['a.u.'] = energy['Hartree'] = energy['Ha'] =  1.0
    energy['eV']                                                     =  27.211386245988   #http://physics.nist.gov/cgi-bin/cuu/Value?hrev
    energy['kJ/mol']                                                 =  4.359744650e-21 * constants['NA']  #http://physics.nist.gov/cgi-bin/cuu/Value?hrj
    energy['J']                                                      =  4.359744650e-18
    energy['kcal/mol']                                               =  energy['kJ/mol'] / 4.184
    energy['cm^-1'] = energy['cm-1']                                 =  219474.6313702   #http://physics.nist.gov/cgi-bin/cuu/Value?hrminv
    energy['K'] = energy['J'] / constants['k_B']

    mass = {}
    mass['au'] = mass['a.u.'] = mass['amu'] = 1.0
    mass['kg'] = 1.66053906660e-27
    mass['g'] = mass['kg'] * 1e3

    angle = {}
    angle['degree'] =  angle['deg'] = 1.0
    angle['radian'] =  angle['rad'] = math.pi / 180.0
    angle['grad']   =  100.0 / 90.0
    angle['circle'] =  1.0 / 360.0

    dipole = {}
    dipole['au'] = dipole['a.u.'] =  1.0
    dipole['Cm']                  =  constants['e'] * constants['Bohr_radius'] * 1e-10
    dipole['Debye'] = dipole['D'] =  dipole['Cm'] * constants['c']* 1e21

    forces = {}
    hessian = {}
    stress = {}
    for k,v in energy.items():
        forces[k+'/Angstrom'] = forces[k+'/Ang'] = forces[k+'/A'] = v * rec_distance['1/Angstrom']
        hessian[k+'/Angstrom^2'] = hessian[k+'/Ang^2'] = hessian[k+'/A^2'] = v  * rec_distance['1/Angstrom']**2
        stress[k+'/Angstrom^3'] = stress[k+'/Ang^3'] = stress[k+'/A^3'] = v  * rec_distance['1/Angstrom']**3
        forces[k+'/bohr'] = forces[k+'/au'] = forces[k+'/a.u.'] = v * rec_distance['1/Bohr']
        hessian[k+'/bohr^2'] = hessian[k+'/au^2'] = hessian[k+'/a.u.^2'] = v * rec_distance['1/Bohr']**2
        stress[k+'/bohr^3'] = stress[k+'/au^3'] = stress[k+'/a.u.^3'] = v * rec_distance['1/Bohr']**3
        forces[k+'/m'] = v * rec_distance['1/m']
        hessian[k+'/m^2'] = v  * rec_distance['1/m']**2
        stress[k+'/m^3'] = v  * rec_distance['1/m']**3
    forces['au'] = forces['a.u.'] = forces['Ha/bohr']
    hessian['au'] = hessian['a.u.'] = hessian['Ha/bohr^2']
    stress['au'] = stress['a.u.'] = stress['Ha/bohr^3']
    stress['Pa'] = stress['J/m^3']
    stress['GPa'] = stress['Pa'] * 1e-9
    stress['bar'] = stress['Pa'] * 1e-5
    stress['atm'] = stress['bar'] / 1.01325



    dicts = {}
    dicts['distance'] = distance
    dicts['energy'] = energy
    dicts['mass'] = mass
    dicts['angle'] = angle
    dicts['dipole'] = dipole
    dicts['reciprocal distance'] = rec_distance
    dicts['forces'] = forces
    dicts['hessian'] = hessian
    dicts['stress'] = stress


    def __init__(self):
        raise UnitsError('Instances of Units cannot be created')


    @classmethod
    def find_unit(cls, unit):
        ret = {}
        for quantity in cls.dicts:
            for k in cls.dicts[quantity]:
                if k.lower() == unit.lower():
                    ret[quantity] = k
        return ret


    @classmethod
    def conversion_ratio(cls, inp, out):
        """Return conversion ratio from unit *inp* to *out*."""
        if inp == out:
            return 1.
        inps = cls.find_unit(inp)
        outs = cls.find_unit(out)
        common = set(inps.keys()) & set(outs.keys())
        if len(common) > 0:
            quantity = common.pop()
            d = cls.dicts[quantity]
            return d[outs[quantity]]/d[inps[quantity]]
        else:
            if len(inps) == 0 and len(outs) == 0:
                raise UnitsError("Unsupported units: '{}' and '{}'".format(inp, out))
            if len(inps) > 0 and len(outs) > 0:
                raise UnitsError("Invalid unit conversion: '{}' is a unit of {} and '{}' is a unit of {}".format(inp, ', '.join(list(inps.keys())), out, ', '.join(list(outs.keys()))))
            else: #exactly one of (inps,outs) empty
                invalid, nonempty = (out,inps) if len(inps) else (inp,outs)
                if len(nonempty) == 1:
                    quantity = list(nonempty.keys())[0]
                    raise UnitsError("Invalid unit conversion: {} is not supported. Supported units for {}: {}".format(invalid, quantity, ', '.join(list(cls.dicts[quantity].keys()))))
                else:
                    raise UnitsError("Invalid unit conversion: {} is not a supported unit for {}".format(invalid, ', '.join(list(nonempty.keys()))))




    @classmethod
    def convert(cls, value, inp, out):
        """Convert *value* from unit *inp* to *out*.

        *value* can be a single number or a container (list, tuple, numpy.array etc.). In the latter case a container of the same type and length is returned. Conversion happens recursively, so this method can be used to convert, for example, a list of lists of numbers, or any other hierarchical container structure. Conversion is applied on all levels, to all values that are numbers (also numpy number types). All other values (strings, bools etc.) remain unchanged.
        """
        if inp == out:
            return value
        if value is None or isinstance(value, (bool, str)):
            return value
        if isinstance(value, collections.abc.Iterable):
            t = type(value)
            if t == np.ndarray:
                t = np.array
            v = [cls.convert(i, inp, out) for i in value]
            return t(v)
        if isinstance(value, (int, float, np.generic)):
            return value * cls.conversion_ratio(inp,out)
        return value

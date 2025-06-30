import collections
import math
from typing import Dict

from scm.plams.core.errors import UnitsError
import numpy as np

__all__ = ["Units"]


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

        -   ``Angstrom``, ``A``, ``Ang``
        -   ``Bohr``, ``au``, ``a.u.``
        -   ``nm``
        -   ``pm``

    *   reciprocal distance:

        -   ``1/Angstrom``, ``1/A``, ``Angstrom^-1``, ``A^-1``,
        -   ``1/Bohr``, ``Bohr^-1``

    *   angle:

        -    ``degree``, ``deg``
        -    ``radian``, ``rad``
        -    ``grad``
        -    ``circle``

    *   charge:

        -    ``coulomb``, ``C``
        -    ``e``

    *   energy:

        -   ``au``, ``a.u.``, ``Hartree``
        -   ``eV``
        -   ``kcal/mol``
        -   ``kJ/mol``
        -   ``cm^-1``, ``cm-1``
        -   ``K``, ``Kelvin``
        -   ``Hz``, ``Hertz``
        -   ``THz``

    *   dipole moment:

        -   ``au``, ``a.u.``, ``e*bohr``
        -   ``Debye``, ``D``
        -  All charge units multiplied by distance units, for example
        -   ``eA``,  ``e*A``
        -   ``Cm``,  ``C*m``

    *   molecular polarizability:

        -   ``au``, ``a.u.``, ``(e*bohr)^2/hartree``
        -   ``e*A^2/V``
        -   ``C*m^2/V``
        -   ``cm^3``
        -   ``bohr^3``
        -   ``A^3``, ``angstrom^3``, ``Ang^3``

    *   forces:

        -   All energy units divided by angstrom or bohr, for example
        -   ``eV/angstrom``
        -   ``hartree/bohr``

    *   hessian:

        -   All energy units divided by angstrom^2 or bohr^2, for example
        -   ``eV/angstrom^2``
        -   ``hartree/bohr^2``

    *   pressure:

        -   All energy units divided by angstrom^3 or bohr^3, for example
        -   ``eV/angstrom^3``
        -   ``hartree/bohr^3``
        -   And some more:
        -   ``Pa``
        -   ``GPa``
        -   ``bar``
        -   ``atm``


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
    constants["Bohr_radius"] = 0.529177210903  # A     http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
    constants["Avogadro_constant"] = constants["NA"] = (
        6.022140857e23  # 1/mol http://physics.nist.gov/cgi-bin/cuu/Value?na
    )
    constants["speed_of_light"] = constants["c"] = 299792458  # m/s http://physics.nist.gov/cgi-bin/cuu/Value?c
    constants["electron_charge"] = constants["e"] = 1.6021766208e-19  # C http://physics.nist.gov/cgi-bin/cuu/Value?e
    constants["Boltzmann"] = constants["k_B"] = 1.380649e-23  # J/K
    constants["vacuum_electric_permittivity"] = (
        8.8541878128e-12  # F*m-1=C/(V*m) https://physics.nist.gov/cgi-bin/cuu/Value?ep0
    )

    distance = {}
    distance["A"] = distance["Angstrom"] = distance["Ang"] = 1.0
    distance["Bohr"] = distance["bohr"] = distance["a.u."] = distance["au"] = 1.0 / constants["Bohr_radius"]
    distance["nm"] = distance["A"] / 10.0
    distance["pm"] = distance["A"] * 100.0
    distance["m"] = distance["A"] * 1e-10

    rec_distance = {}
    rec_distance["1/A"] = rec_distance["1/Ang"] = rec_distance["1/Angstrom"] = rec_distance["A^-1"] = rec_distance[
        "Ang^-1"
    ] = rec_distance["Angstrom^-1"] = 1.0
    rec_distance["1/m"] = rec_distance["m^-1"] = 1e10
    rec_distance["1/Bohr"] = rec_distance["Bohr^-1"] = constants["Bohr_radius"]

    energy = {}
    energy["au"] = energy["a.u."] = energy["Hartree"] = energy["Ha"] = 1.0
    energy["eV"] = 27.211386245988  # http://physics.nist.gov/cgi-bin/cuu/Value?hrev
    energy["kJ/mol"] = 4.359744650e-21 * constants["NA"]  # http://physics.nist.gov/cgi-bin/cuu/Value?hrj
    energy["J"] = 4.359744650e-18
    energy["kcal/mol"] = energy["kJ/mol"] / 4.184
    energy["cm^-1"] = energy["cm-1"] = 219474.6313702  # http://physics.nist.gov/cgi-bin/cuu/Value?hrminv
    energy["K"] = energy["J"] / constants["k_B"]
    energy["Hz"] = energy["Hertz"] = (
        6.57968389898681e15  # https://physics.nist.gov/cgi-bin/cuu/Convert?exp=0&num=1&From=hr&To=hz&Action=Only+show+factor
    )
    energy["THz"] = energy["Hz"] / 1e12

    mass = {}
    mass["au"] = mass["a.u."] = mass["amu"] = 1.0
    mass["kg"] = 1.66053906660e-27
    mass["g"] = mass["kg"] * 1e3

    time = {}
    time["s"] = 1.0
    time["ms"] = time["s"] * 1e3
    time["us"] = time["s"] * 1e6
    time["ns"] = time["s"] * 1e9
    time["ps"] = time["s"] * 1e12
    time["fs"] = time["s"] * 1e15
    time["au"] = time["a.u."] = time["s"] / 2.4188843265857e-17  # https://physics.nist.gov/cgi-bin/cuu/Value?aut

    angle = {}
    angle["degree"] = angle["deg"] = 1.0
    angle["radian"] = angle["rad"] = math.pi / 180.0
    angle["grad"] = 100.0 / 90.0
    angle["circle"] = 1.0 / 360.0

    charge = {}
    charge["a.u."] = charge["au"] = charge["e"] = 1.0
    charge["C"] = charge["coulomb"] = constants["e"]

    dipole = {}
    for k, v in charge.items():
        if (k == "au") or (k == "a.u."):  # remove 'au','a.u.' options
            continue
        for k1, v1 in distance.items():
            if (k1 == "au") or (k1 == "a.u."):
                continue
            dipole[k + "*" + k1] = v * v1
            dipole[k + k1] = v * v1
    dipole["au"] = dipole["a.u."] = dipole["e*bohr"]
    dipole["debye"] = dipole["D"] = dipole["Cm"] * constants["c"] * 1e21

    # from support info https://doi.org/10.48550/arXiv.2310.13310 it is preferable to highlight that this is molecular polarizability,
    # it should be also /mol units but it is usually omitted and for consistency in the dipole units I removed, but both dipole and molecular_polarizability should have /mol
    molecular_polarizability = {}
    molecular_polarizability["au"] = molecular_polarizability["a.u."] = molecular_polarizability[
        "e^2*bohr^2/hartree"
    ] = molecular_polarizability["(e*bohr)^2/hartree"] = 1.0
    molecular_polarizability["e*A^2/V"] = molecular_polarizability["e^2*A^2/eV"] = molecular_polarizability[
        "(e*A)^2/eV"
    ] = (constants["Bohr_radius"] ** 2 / energy["eV"])
    molecular_polarizability["C*m^2/V"] = molecular_polarizability["e*A^2/V"] * 1e-20 * constants["e"]
    molecular_polarizability["cm^3"] = (
        molecular_polarizability["C*m^2/V"] / (4 * np.pi * constants["vacuum_electric_permittivity"]) * 1e6
    )  # form https://en.wikipedia.org/wiki/Polarizability that refs Atkins book
    molecular_polarizability["A^3"] = molecular_polarizability["Ang^3"] = molecular_polarizability["Angstrom^3"] = (
        molecular_polarizability["cm^3"] * 1e24
    )
    molecular_polarizability["bohr^3"] = molecular_polarizability["Ang^3"] / constants["Bohr_radius"] ** 3

    forces = {}
    hessian = {}
    stress = {}
    for k, v in energy.items():
        for k1, v1 in distance.items():
            forces[k + "/" + k1] = v / v1
            hessian[k + "/" + k1 + "^2"] = v / v1**2
            stress[k + "/" + k1 + "^3"] = v / v1**3
    forces["au"] = forces["a.u."] = forces["Ha/bohr"]
    hessian["au"] = hessian["a.u."] = hessian["Ha/bohr^2"]
    stress["au"] = stress["a.u."] = stress["Ha/bohr^3"]
    stress["Pa"] = stress["J/m^3"]
    stress["GPa"] = stress["Pa"] * 1e-9
    stress["bar"] = stress["Pa"] * 1e-5
    stress["atm"] = stress["bar"] / 1.01325

    dicts = {}
    dicts["distance"] = distance
    dicts["energy"] = energy
    dicts["mass"] = mass
    dicts["time"] = time
    dicts["angle"] = angle
    dicts["dipole"] = dipole
    dicts["reciprocal distance"] = rec_distance
    dicts["forces"] = forces
    dicts["hessian"] = hessian
    dicts["stress"] = stress
    dicts["charge"] = charge
    dicts["molecular_polarizability"] = molecular_polarizability

    # Precomputed a dict mapping lowercased unit names to quantityName:conversionFactor pairs
    quantities_for_unit: Dict[str, Dict[str, float]] = {}
    for quantity in dicts:
        for unit, factor in dicts[quantity].items():
            unit = unit.lower()
            if unit not in quantities_for_unit:
                quantities_for_unit[unit] = {}
            quantities_for_unit[unit][quantity] = factor

    def __init__(self):
        raise UnitsError("Instances of Units cannot be created")

    @classmethod
    def find_unit(cls, unit):
        ret = {}
        u = unit.lower()
        quantities = cls.quantities_for_unit.get(u, {})
        for quantity in quantities:
            for k in cls.dicts[quantity]:
                if k.lower() == u:
                    ret[quantity] = k
                    break
        return ret

    @classmethod
    def conversion_ratio(cls, inp, out) -> float:
        """Return conversion ratio from unit *inp* to *out*."""
        if inp == out:
            return 1.0
        inps = cls.quantities_for_unit.get(inp.lower(), {})
        outs = cls.quantities_for_unit.get(out.lower(), {})
        common = set(inps.keys()) & set(outs.keys())
        if len(common) > 0:
            quantity = common.pop()
            return outs[quantity] / inps[quantity]
        else:
            if len(inps) == 0 and len(outs) == 0:
                raise UnitsError("Unsupported units: '{}' and '{}'".format(inp, out))
            if len(inps) > 0 and len(outs) > 0:
                raise UnitsError(
                    "Invalid unit conversion: '{}' is a unit of {} and '{}' is a unit of {}".format(
                        inp, ", ".join(list(inps.keys())), out, ", ".join(list(outs.keys()))
                    )
                )
            else:  # exactly one of (inps,outs) empty
                invalid, nonempty = (out, inps) if len(inps) else (inp, outs)
                if len(nonempty) == 1:
                    quantity = list(nonempty.keys())[0]
                    raise UnitsError(
                        "Invalid unit conversion: {} is not supported. Supported units for {}: {}".format(
                            invalid, quantity, ", ".join(list(cls.dicts[quantity].keys()))
                        )
                    )
                else:
                    raise UnitsError(
                        "Invalid unit conversion: {} is not a supported unit for {}".format(
                            invalid, ", ".join(list(nonempty.keys()))
                        )
                    )

    @classmethod
    def convert(cls, value, inp, out):
        """Convert *value* from unit *inp* to *out*.

        *value* can be a single number or a container (list, tuple, numpy.array etc.). In the latter case a container of the same type and length is returned. Conversion happens recursively, so this method can be used to convert, for example, a list of lists of numbers, or any other hierarchical container structure. Conversion is applied on all levels, to all values that are numbers (also numpy number types). All other values (strings, bools etc.) remain unchanged.
        """
        if value is None or isinstance(value, (bool, str)) or inp == out:
            return value
        if isinstance(value, collections.abc.Iterable):
            t = type(value)
            if t == np.ndarray:
                t = np.array  # type: ignore
            v = [cls.convert(i, inp, out) for i in value]
            return t(v)  # type: ignore
        if isinstance(value, (int, float, np.generic)):
            return value * cls.conversion_ratio(inp, out)  # type: ignore
        return value

    @classmethod
    def ascii2unicode(cls, string):
        """
        Converts '^2' to '²' etc., for prettier printing of units.
        """
        if string is None:
            return ""
        ret = (
            string.replace("^-1", "⁻¹")
            .replace("angstrom", "Å")
            .replace("^2", "²")
            .replace("^3", "³")
            .replace("degree", "°")
            .replace("deg.", "°")
            .replace("Ang", "Å")
            .replace("*", "⋅")
        )
        return ret

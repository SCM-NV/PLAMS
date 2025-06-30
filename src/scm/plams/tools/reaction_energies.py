import os

from scm.plams.interfaces.adfsuite.ams import AMSJob
import numpy as np
from scm.plams.mol.molecule import Molecule
from scm.plams import ReactionEquation

__all__ = ["get_stoichiometry", "balance_equation", "reaction_energy"]


def get_stoichiometry(job_or_molecule_or_path, as_dict=True):
    r = job_or_molecule_or_path
    d = None
    if isinstance(r, AMSJob):
        d = r.molecule.get_formula(as_dict=as_dict)
    elif isinstance(r, Molecule):
        d = r.get_formula(as_dict=as_dict)
    elif isinstance(r, dict):
        d = r.copy()
    elif isinstance(r, str):
        if os.path.isdir(r):
            d = AMSJob.load_external(r).molecule.get_formula(as_dict=as_dict)
        elif os.path.exists(r):
            try:
                d = Molecule(r).get_formula(as_dict=as_dict)
            except:
                d = AMSJob.load_external(r).molecule.get_formula(as_dict=as_dict)
        else:
            raise ValueError("The path {} does not exist.".format(r))

    else:
        raise TypeError("expected type AMSJob or dict but received {}".format(type(r)))

    return d


def balance_equation_new(reactants, products, normalization="r0", normalization_value=1.0):
    """
    Calculate stoichiometric coefficients (meant to replace balance_equation method)

    reactants: a list of amsjobs, or a list of paths to ams.results folders or ams.rkf files or .xyz files, or a list of Molecules, or a list of stoichiometry dicts, or a list of Molecules
        The reactants

    products: a list of amsjobs, or a list of paths to ams.results folders or ams.rkf files, or a list of Molecules or .xyz files, or a list of stoichiometry dicts, or a list of Molecules
        The products

    Returns: a 2-tuple (coeffs_reactants, coeffs_products)
        coeffs_reactants is a list with length == len(reactants)
        coeffs_products is a list with length == len(products)

    normalization: str
        'r0' for the first reactant, 'r1' for the second reactant, etc.
        'p0' for the first product, 'p1' for the second product, etc.
        This normalizes the chemical equation such that the coefficient in front of the specified species is normalization_value

    normalization_value: float or None
        The value to which the compound defined with 'normalization' should be normalized
        if None, whatever integer value the ReactionEquation object produces should be used.
    """

    def get_formulas(list_of_jobs):
        """
        Convert the list of molecules to a list of molecular formulas
        """
        formulas = []
        for r in list_of_jobs:
            d = get_stoichiometry(r)
            formula = "".join(["%s%i" % (el, n) for el, n in d.items()])
            formulas.append(formula)
        return formulas

    def get_normalization_index(normalization):
        if normalization.startswith("r"):
            normalization_index = int(normalization.split("r")[1])
            if normalization_index >= num_reactants:
                raise ValueError(
                    "Reactant index {} specified, but max value allowed is {}".format(
                        normalization_index, num_reactants - 1
                    )
                )
        elif normalization.startswith("p"):
            normalization_index = int(normalization.split("p")[1])
            if normalization_index >= num_products:
                raise ValueError(
                    "Product index {} specified, but max value allowed is {}".format(
                        normalization_index, num_products - 1
                    )
                )
            normalization_index += num_reactants
        else:
            raise ValueError(
                "Unknown normalization: {}. Should be r0, r1, r2, ... (for reactants), p0, p1, p2 ... (for products)"
            )

        return normalization_index

    if len(reactants) == 0:
        raise ValueError("The reactants list is empty.")
    if len(products) == 0:
        raise ValueError("The products list is empty.")

    # Set up the input, which can be lists of formulas, or lists of PLAMS molecule objects
    num_reactants = len(reactants)
    reactants = get_formulas(reactants)
    num_products = len(products)
    products = get_formulas(products)

    # Set up the minimal numbers of the coefficients
    ind = get_normalization_index(normalization)
    min_coeffs = np.zeros(num_reactants + len(products))
    min_coeffs[ind] = 1

    # Solve
    reaction = ReactionEquation(reactants, products)
    coeffs = reaction.balance(min_coeffs)

    if coeffs is None:
        strings = ["Something went wrong when solving the system of linear equations."]
        strings += ["Verify that the chemical equation can be balanced at all."]
        text = " ".join(strings)
        raise RuntimeError(text)

    # Normalize in the specified molecule
    if normalization_value is not None:
        coeffs = coeffs.astype(np.float64)
        coeffs /= coeffs[ind]
        coeffs *= normalization_value

    return coeffs[:num_reactants], coeffs[num_reactants:]


def balance_equation(reactants, products, normalization="r0", normalization_value=1.0):
    """
    Calculate stoichiometric coefficients
    This only works if
    * number_of_chemical_elements == len(reactants)+len(products), OR
    * number_of_chemical_elements == len(reactants)+len(products)-1

    Returns: a 2-tuple (coeffs_reactants, coeffs_products)
        coeffs_reactants is a list with length == len(reactants)
        coeffs_products is a list with length == len(products)

    reactants: a list of amsjobs, or a list of paths to ams.results folders or ams.rkf files or .xyz files, or a list of Molecules, or a list of stoichiometry dicts, or a list of Molecules
        The reactants

    products: a list of amsjobs, or a list of paths to ams.results folders or ams.rkf files, or a list of Molecules or .xyz files, or a list of stoichiometry dicts, or a list of Molecules
        The products

    normalization: str
        'r0' for the first reactant, 'r1' for the second reactant, etc.
        'p0' for the first product, 'p1' for the second product, etc.
        This normalizes the chemical equation such that the coefficient in front of the specified species is normalization_value

    normalization_value : float
        The coefficient to normalize to

    EXAMPLE:

    .. code-block:: python

        balance_equation(
            reactants=[
                {'N': 2, 'H': 8, 'Cr': 2, 'O': 7}
            ],
            products=[
                {'Cr': 2, 'O': 3},
                {'N': 2},
                {'H': 2, 'O': 1}
            ])

    The above returns a tuple ``([1.0], [1.0, 1.0, 1.0, 4.0])``


    """

    def get_stoichiometries_and_elements(list_of_jobs):
        stoich = []
        elements = set()
        for r in list_of_jobs:
            d = get_stoichiometry(r)
            stoich.append(d)
            for k in d:
                elements.add(k)
        return stoich, elements

    def get_normalization_index(normalization):
        if normalization.startswith("r"):
            normalization_index = int(normalization.split("r")[1])
            if normalization_index >= num_reactants:
                raise ValueError(
                    "Reactant index {} specified, but max value allowed is {}".format(
                        normalization_index, num_reactants - 1
                    )
                )
        elif normalization.startswith("p"):
            normalization_index = int(normalization.split("p")[1])
            if normalization_index >= num_products:
                raise ValueError(
                    "Product index {} specified, but max value allowed is {}".format(
                        normalization_index, num_products - 1
                    )
                )
            normalization_index += num_reactants
        else:
            raise ValueError(
                "Unknown normalization: {}. Should be r0, r1, r2, ... (for reactants), p0, p1, p2 ... (for products)"
            )

        return normalization_index

    if len(reactants) == 0:
        raise ValueError("The reactants list is empty.")
    if len(products) == 0:
        raise ValueError("The products list is empty.")

    stoich_r, elements_r = get_stoichiometries_and_elements(reactants)
    stoich_p, elements_p = get_stoichiometries_and_elements(products)
    elements = elements_r
    elements.update(elements_p)  # hopefully not necessary
    elements = list(elements)
    num_reactants = len(stoich_r)
    num_products = len(stoich_p)

    # EXAMPLE:
    # aCH4 + bO2 --> cCO2 + dH2O
    # mat =
    # [[-1 0    1 0], #C
    #  [-4 0    0 2], #H
    #  [0 -2    2 1]] #O
    #  CH4 O2 CO2 H2O
    #
    # Al2(SO4)3 + Ca(OH)2 → Al(OH)3 + CaSO4
    # mat = np.array([[-2,-3,-12,0,0,0],[0,0,0,-1,-2,-2],[1,0,0,0,3,3],[0,1,4,1,0,0]]).T

    mat = np.array([[-s.get(e, 0) for s in stoich_r] + [s.get(e, 0) for s in stoich_p] for e in elements])

    if mat.shape[0] >= mat.shape[1]:
        u, s, vh = np.linalg.svd(mat)
        coeffs = np.compress(s <= 1e-12, vh, axis=0)
    elif mat.shape[0] == mat.shape[1] - 1:
        # e.g. 3x4 matrix, so add a [1,1,1,1] row to find a particular solution
        newmat = np.concatenate((mat, np.array([[1] * mat.shape[1]])), axis=0)
        b = np.array([0] * (newmat.shape[0] - 1) + [1.0])
        try:
            coeffs = np.linalg.solve(newmat, b)
        except Exception as e:
            raise RuntimeError(
                "Something went wrong when solving the system of linear equations. Verify that the chemical equation can be balanced at all, and that it can be balanced uniquely except for multiplication by a constant. {}\nA={}\nb={}".format(
                    e, newmat, b
                )
            )
    else:
        raise ValueError(
            "The number of chemical elements must equal the number of molecules, or (the number of molecules-1). You have {} chemical elements: {}, and {} molecules".format(
                len(elements), elements, num_reactants + num_products
            )
        )

    coeffs = coeffs.ravel()
    if len(coeffs) != num_reactants + num_products:
        raise RuntimeError(
            "Something went wrong when solving the system of linear equations. Verify that the chemical equation can be balanced at all, and that it can be balanced uniquely except for multiplication by a constant."
        )

    normalization_index = get_normalization_index(normalization)
    if np.abs(coeffs[normalization_index]) < 1e-12:
        raise RuntimeError(
            "Trying to normalize an extremely small coefficient (close to zero): {coeffs[normalization_index]}."
        )
    coeffs /= coeffs[normalization_index]
    coeffs *= normalization_value

    # double-check that the equation is balanced
    if abs(np.sum(mat @ coeffs.reshape(-1, 1))) > 1e-10:
        raise RuntimeError("Stoichiometry double-check failed. mat = {}, coeffs = {}".format(mat, coeffs))

    return list(coeffs[:num_reactants]), list(coeffs[num_reactants:])


def reaction_energy(reactants, products, normalization="r0", unit="hartree"):
    """

    Calculates a reaction energy from an unbalanced chemical equation (the equation is first balanced)

    reactants: a list of amsjobs or paths to ams results folders,
        The recatnts
    products: a list of amsjobs or paths to ams results folders
        The products
    normalization: str
        normalize the chemical equation by setting the corresponding coefficient to 1.

        * 'r0': first reactant
        * 'r1': second reactant, ...
        * 'p0: first product,
        * 'p1': second product, ...

    unit: str
        Unit of the reaction energy

    Returns: a 3-tuple (coeffs_reactants, coeffs_products, reaction_energy)

    """

    my_reactants = [AMSJob.load_external(x) for x in reactants]
    my_products = [AMSJob.load_external(x) for x in products]
    coeffs_r, coeffs_p = balance_equation(my_reactants, my_products, normalization)

    energies = np.array([job.results.get_energy(unit=unit) for job in my_reactants + my_products])

    coeffs = np.concatenate(([-x for x in coeffs_r], coeffs_p))
    reaction_energy = np.dot(coeffs, energies)

    return coeffs_r, coeffs_p, reaction_energy

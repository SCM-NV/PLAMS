import numpy

from scm.plams.mol.molecule import Molecule
from scm.plams.interfaces.molecule.rdkit import to_smiles
from scm.plams.core.functions import requires_optional_package
from scm.plams.core.errors import PlamsError
from dataclasses import dataclass

from typing import Union, Optional, Iterable, Sequence, List, Dict, Set, Mapping, Literal, Any
import numpy.typing
import numpy as np


class ReactionEquationError(PlamsError):
    """
    Error with the reaction equation
    """


class UninitializedReactionError(PlamsError):
    """
    Reaction equation state not properly initialized
    """


class ReactionEquation:
    """
    Class representing a reaction, which can then be balanced to any selected set of molecules
    """

    nullspace_methods = ["sympy", "plams"]

    def __init__(
        self, reactants: Sequence[Union[Molecule, str]], products: Sequence[Union[Molecule, str]], method: str = "plams"
    ):
        """
        Initiate the reaction

        * ``reactants`` -- List of PLAMS molecules or strings representing a molecular formula
        * ``products`` -- List of PLAMS molecules or strings representing a molecular formula
        """
        self.method = method
        self.print_as_smiles = False

        # Not to be changed by user
        self._rformulas = [m.get_formula() if isinstance(m, Molecule) else m for m in reactants]
        self._pformulas = [m.get_formula() if isinstance(m, Molecule) else m for m in products]
        rsmiles = [to_smiles(m) if isinstance(m, Molecule) else None for m in reactants]
        psmiles = [to_smiles(m) if isinstance(m, Molecule) else None for m in products]
        self.rsmiles: Optional[List[Union[str, None]]] = None if None in rsmiles else rsmiles  # type: ignore[assignment]
        self.psmiles: Optional[List[Union[str, None]]] = None if None in psmiles else psmiles  # type: ignore[assignment]

        # Set the charges
        self._rcharges = []
        for m in reactants:
            q = 0.0
            if isinstance(m, Molecule):
                if "charge" in m.properties:
                    q = m.properties.charge
            self._rcharges.append(q)
        self._pcharges = []
        for m in products:
            q = 0.0
            if isinstance(m, Molecule):
                if "charge" in m.properties:
                    q = m.properties.charge
            self._pcharges.append(q)

        # General settings for this system
        self.elements: Optional[List[str]] = None
        self.reactant_elements: Optional[List[Dict[str, int]]] = None
        self.product_elements: Optional[List[Dict[str, int]]] = None
        self.matrix: Optional[numpy.typing.NDArray] = None
        self.rank = None
        self.basis: Optional[numpy.typing.NDArray] = None

        # Settings that will change with each balance call (depend on min_coeffs)
        self.min_coeffs: Optional[numpy.typing.NDArray] = None
        self.model = None
        self.solver = None
        self.coeffs = None
        self.equivalent_coeffs = None
        self.message = "Unsolved"

        # Turn off the error writing to stdout
        import logging

        logging.getLogger("pyomo.core").setLevel(logging.ERROR)

    def prepare_state(self) -> None:
        """
        Do the time consuming stuff that needs to be done before we call balance
        """
        # Set the element data from the formulas
        self._set_elements()

        # Get the matrix for the homogenous set of linear equations
        self._set_matrix()

        # Get the nullspace basis
        self.basis = self.get_nullspace_basis()
        # print ('basis: ')
        # print (self.basis)

    def balance(self, min_coeffs: Optional[Iterable[int]] = None) -> Optional[numpy.typing.NDArray]:
        """
        Balance the equation to any set of molecules specified in min_coeffs

        * ``min_coeffs`` -- Vector representing the minimal allowable (integer) coefficient of each molecule
                        [0, 0, 1, 0, 0]
        """
        if min_coeffs is not None:
            min_coeffs = numpy.array(min_coeffs)

        self.coeffs = None
        self.equivalent_coeffs = None
        self.message = "Unsolved"

        # Set the matrix and the nullspace basis, if they were not already set
        if self.matrix is None:
            self.prepare_state()
        self.set_minimum_coefficients(min_coeffs)

        # Exclude basis vectors based on the provided min_coeffs
        basis = self.get_reduced_basis()
        if len(basis) == 0:
            self.message = "Empty nullspace"
            return None

        # If the basis containt no coefficients for either products or reactants, it also fails
        nreactants = len(self._rformulas)
        if abs(basis[:, :nreactants]).sum() == 0:
            self.message = "Empty nullspace"
            return None
        if abs(basis[:, nreactants:]).sum() == 0:
            self.message = "Empty nullspace"
            return None

        # Now we work with the basis to get the coefficients
        # We use ILP to do so
        self.setup_optimizer(basis)
        coeffs = self.optimize_coefficients()
        if coeffs is None:
            return None
        coeffs = self.refine_optimization()

        # Check the validity of the results
        null_vector = self.matrix @ coeffs
        if (null_vector != 0).any():
            self.message = "Optimization failed"
            self.coeffs = None
            return None

        return coeffs

    def set_minimum_coefficients(self, min_coeffs: Optional[numpy.typing.NDArray] = None) -> None:
        """
        Set the minimum coefficient values

        * ``min_coeffs`` -- Vector representing the minimal allowable (integer) coefficient of each molecule
                        [0, 0, 1, 0, 0]
        """
        if min_coeffs is not None:
            min_coeffs = numpy.array(min_coeffs)
            if len(min_coeffs[min_coeffs != 0]) == 0:
                msg = "At least one non-zero coefficient needs to be provided as the min_coeffs argument."
                raise ReactionEquationError(msg)

        nreactants = len(self._rformulas)
        if min_coeffs is None:
            if self.matrix is None:
                raise UninitializedReactionError("Uninitialized: Call prepare_state() method first.")
            min_coeffs = numpy.zeros(self.matrix.shape[1])
            min_coeffs[nreactants] = 1
        self.min_coeffs = min_coeffs

    def get_nullspace_basis(self) -> numpy.typing.NDArray:
        """
        Here we try to get the nullspace vectors
        """
        use_sympy = True
        try:
            import sympy
        except ImportError:
            use_sympy = False

        if self.method == "sympy" and use_sympy:
            matrix = sympy.Matrix(self.matrix)
            sol = matrix.nullspace()
            # Convert to array of floats
            basis = []
            for c in sol:
                c = [v[0] for v in c.tolist()]
                basis.append(c)
            basis_array = numpy.array([[float(x) for x in v] for v in basis])

        elif self.method == "plams" or self.method == "sympy":
            from scm.plams.tools.plams_matrix import PLAMSMatrix

            plams_matrix = PLAMSMatrix(self.matrix)
            basis_array = plams_matrix.nullspace()

        else:
            raise ReactionEquationError("Null space method not known")

        return basis_array

    def get_reduced_basis(self) -> numpy.typing.NDArray:
        """
        Select only the rows in basis that share the same block with the main product
        """

        def get_row_indices(ind: int) -> numpy.typing.NDArray:
            """
            Get the indices of the rows in the same block with compound ind
            """
            if self.basis is None:
                raise UninitializedReactionError("Uninitialized: prepare_state() should be called.")
            # First find a row that has a nonzero value at position ind
            rowmap = None
            for i, row in enumerate(self.basis):
                if row[ind] != 0:
                    rowmap = (row != 0) * numpy.ones(self.basis.shape)
                    break

            indexmap = numpy.ones(len(self.basis))
            if rowmap is not None:
                indexmap = (abs(self.basis * rowmap) > 0).sum(axis=1)
            return indexmap

        # Add the rows for the blocks representing each mandatory compound
        if self.basis is None:
            raise UninitializedReactionError("Uninitialized: Call prepare_state() method first.")
        if self.min_coeffs is None:
            raise UninitializedReactionError("Uninitialized: Call set_minimum_coefficients() method first.")
        nmols = len(self.min_coeffs)
        nonzero_indices = iter(numpy.arange(nmols)[self.min_coeffs > 0])
        indexmap = get_row_indices(next(nonzero_indices))
        for i in nonzero_indices:
            newmap = get_row_indices(i)
            # Check for independent reactions
            oldbasis = (abs(self.basis[indexmap > 0])).sum(axis=0) > 0
            newbasis = (abs(self.basis[newmap > 0])).sum(axis=0) > 0
            # They are only independent if these two bases have not overlap at all
            if (oldbasis * newbasis).sum() == 0:
                self.message = "(independent equations)"
            indexmap += newmap

        # Cut out the relevant rows basis
        basis = self.basis[indexmap > 0]
        return basis

    @requires_optional_package("pyomo")
    def setup_optimizer(self, basis: numpy.typing.NDArray) -> None:
        """
        Solve the problem using ILP.

        min(sum(x)), with basis@y = x, and x >= min_coeffs

        * ``basis`` -- m x n matrix, representing the null space vectors as rows
        """
        from pyomo.environ import ConcreteModel
        from pyomo.environ import Var
        from pyomo.environ import Objective
        from pyomo.environ import SolverFactory
        from pyomo.environ import ConstraintList
        from pyomo.environ import Constraint
        from pyomo.environ import NonNegativeIntegers
        from pyomo.environ import Integers
        from pyomo.environ import minimize

        m = basis.shape[0]
        n = basis.shape[1]
        mat = basis.transpose()

        # Initiate the pyomo model
        model = ConcreteModel()

        # Set the variables
        # Note: I can set upper an lower bonds of variables directly: model.x[0].setub(0)
        model.x = Var(range(n), domain=NonNegativeIntegers)
        model.y = Var(range(m), domain=Integers)

        # Set the objective function we want to minimize
        f = sum(model.x[i] for i in range(n))
        model.objective = Objective(expr=f, sense=minimize)

        # Set up the main constraint
        model.demands = ConstraintList()
        for i in range(n):
            expr = sum(mat[i, j] * model.y[j] for j in range(m)) == model.x[i]
            model.demands.add(expr)

        # Set up the secondary constraint
        model.constraints = ConstraintList()
        if self.min_coeffs is None:
            raise UninitializedReactionError("Uninitialized: Call set_minimum_coefficients() method first.")
        for i in range(n):
            expr = model.x[i] >= self.min_coeffs[i]
            model.constraints.add(expr)

        # Set up the charge costraint
        charges = self._rcharges + [-q for q in self._pcharges]
        if True in [abs(q) > 1e-10 for q in charges]:
            expr = sum(model.x[i] * charges[i] for i in range(n)) == 0.0
            model.charge_constraint = Constraint(expr=expr)

        # Set up the solver
        solver = SolverFactory("cbc")

        # Set the instance variables
        self.model = model
        self.solver = solver

    def optimize_coefficients(self) -> numpy.typing.NDArray:
        """
        Solve the problem using ILP.

        min(sum(x)), with basis@y = x, and x >= min_coeffs
        """
        if self.model is None:
            raise UninitializedReactionError("Model not properly initialized. Call setup_optimizer method.")
        n = len(self.model.x)

        # Solve the problem
        result = self.solver.solve(self.model)
        if result.solver.status != "ok":
            self.message = "Optimization failed"
            return None
        coeffs = numpy.array([int(self.model.x[j].value) for j in range(n)])

        self.coeffs = coeffs
        message = "Success"
        if self.message != "Unsolved":
            message = "%s %s" % (message, self.message)
        self.message = message
        return coeffs

    def refine_optimization(self) -> numpy.typing.NDArray:
        """
        Find all solutions with the same cost, and select the one with the lowest indices

        Note: This is done by repressing the highest coefficient in the previous solution,
              and optimizing again, untill no improvement is made.
        """
        if self.model is None:
            raise UninitializedReactionError("Model not properly initialized. Call balance() method directly.")
        n = len(self.model.x)
        coeffs = self.coeffs
        cost = self.model.objective()
        sum_indices = (numpy.arange(n) * coeffs).sum()

        all_solutions = [coeffs]
        new_coeffs = coeffs.copy()
        upper = {}
        while 1:
            # Add upper boundary to x
            largest_index = numpy.arange(n)[(new_coeffs > 0) * (self.min_coeffs == 0)][-1]
            c = new_coeffs[largest_index]
            self.model.x[largest_index].setub(c - 1)
            upper[largest_index] = c - 1

            # I could try warmstarting from previous x: self.solver.solve(self.model, warmstart=True)
            # I can set a time limit with: self.solver.solve(self.model, timelimit=30)
            # Info on options:
            # https://pyomo.readthedocs.io/en/stable/library_reference/appsi/appsi.solvers.cbc.html

            # Optimize ceofficients again
            result = self.solver.solve(self.model)
            new_coeffs = numpy.array([int(self.model.x[j].value) for j in range(n)])

            # Stop if no equivalent solution was found
            if result.solver.status != "ok":
                break
            if self.model.objective() > cost:
                break

            # Only store if the new solution has lower molecule indices
            new_sum_indices = (numpy.arange(n) * new_coeffs).sum()
            if new_sum_indices < sum_indices:
                coeffs = new_coeffs
                sum_indices = new_sum_indices

            all_solutions.append(new_coeffs)

        self.coeffs = coeffs
        self.equivalent_coeffs = all_solutions
        return coeffs

    @property
    def reactants(self) -> List[str]:
        """
        Get the list of reactant formulas
        """
        return self._rformulas

    @property
    def products(self) -> List[str]:
        """
        Get the list of product formulas
        """
        return self._pformulas

    def set_charges(self, reactant_charges: List[float], product_charges: List[float]) -> None:
        """
        Set the charges for reactant and product molecules

        * ``reactant_charges`` -- List of charges of the reactant molecules
        * ``product_charges`` -- List of charges of the product molecules
        """
        if len(reactant_charges) != len(self._rformulas):
            strings = ["Number of supplied reactant charges should be %i " % (len(self._rformulas))]
            strings += ["not %i." % (len(reactant_charges))]
            raise ReactionEquationError("".join(strings))
        if len(product_charges) != len(self._pformulas):
            strings = ["Number of supplied product charges should be %i " % (len(self._pformulas))]
            strings += ["not %i." % (len(product_charges))]
            raise ReactionEquationError("".join(strings))
        self._rcharges = reactant_charges
        self._pcharges = product_charges

    def set_coefficients(self, icoeffs: int) -> None:
        """
        Set one of the entries in self.equivalent_coefficients as the final result
        """
        if self.equivalent_coeffs is None:
            raise UninitializedReactionError("Model not properly initialized. Call balance() method first.")
        self.coeffs = self.equivalent_coeffs[icoeffs]

    def __format__(self, fmt: Optional[str] = None) -> str:
        """Formats the chemical reaction.

        fmt: str
            If 'smiles', will print the reaction using SMILES if possible.

        Example:

        >>> print(f"{reaction:smiles}")
        """
        if self.message == "Unsolved":
            return "Equation not yet balanced"

        if self.coeffs is None:
            return "Inconsistent system: No solution"

        nmol = len(self._rformulas) + len(self._pformulas)
        nreactants = len(self._rformulas)
        indices = numpy.arange(nmol)[self.coeffs > 0]
        rindices = indices[indices < nreactants]
        pindices = indices[indices >= nreactants] - nreactants
        pcoeffs = self.coeffs[nreactants:]

        reactants = self._rformulas
        products = self._pformulas
        if (fmt == "smiles" or self.print_as_smiles) and self.rsmiles is not None:
            reactants = self.rsmiles
        if (fmt == "smiles" or self.print_as_smiles) and self.psmiles is not None:
            products = self.psmiles

        def chargestring(arr, i) -> str:
            if not arr:
                return ""
            if fmt == "smiles" or self.print_as_smiles:
                return ""
            val = arr[i]
            if val == 0:
                return ""
            return f"[{val:+}]"

        ret = (
            " + ".join(f"{self.coeffs[i]} {reactants[i]}{chargestring(self._rcharges, i)}" for i in rindices)
            + " => "
            + " + ".join(f"{pcoeffs[i]} {products[i]}{chargestring(self._pcharges, i)}" for i in pindices)
        )
        return ret

    def __str__(self) -> str:
        """
        Write the balanced reaction
        """
        return f"{self}"

    @property
    def reaction_charge(self) -> Union[int, str]:
        """
        Write the balanced reaction
        """
        if self.message == "Unsolved":
            return "Equation not yet balanced"

        if self.coeffs is None:
            return "Inconsistent system: No solution"

        nmol = len(self._rformulas) + len(self._pformulas)
        nreactants = len(self._rformulas)
        indices = numpy.arange(nmol)[self.coeffs > 0]
        rindices = indices[indices < nreactants]
        pindices = indices[indices >= nreactants] - nreactants
        pcoeffs = self.coeffs[nreactants:]

        reaction_charge = sum([-self.coeffs[i] * self._rcharges[i] for i in rindices])
        reaction_charge += sum([pcoeffs[i] * self._pcharges[i] for i in pindices])
        return reaction_charge

    def matrix_as_string(self, mat: Optional[numpy.typing.NDArray] = None, space: int = 8) -> str:
        """
        Print a numpy matrix in nice format
        """
        if mat is None:
            mat = self.matrix
        if mat is None:
            raise UninitializedReactionError("Uninitialized. Call prepare_state() method first.")
        form = f"%{space}s"
        form2 = f"%{space}.1e"
        lines = []
        for row in mat:
            strings = [form % (str(v)) if len(str(v)) <= space else form2 % (v) for v in row]
            s = " ".join(strings)
            lines.append(s)
        return "\n".join(lines)

    #################
    # Private methods
    #################

    def _set_elements(self) -> None:
        """
        Set the element data from the formulas
        """
        # Get all elements, as well as for each molecule how many of each they contain
        elements: Set[str] = set()
        reactants: List[Dict[str, int]] = []
        for formula in self._rformulas:
            d = self._elements_from_formula(formula)
            reactants.append(d)
            for el in d.keys():
                elements.add(el)
        products: List[Dict[str, int]] = []
        for formula in self._pformulas:
            d = self._elements_from_formula(formula)
            products.append(d)
            for el in d.keys():
                elements.add(el)

        self.elements = sorted([el for el in elements])
        self.reactant_elements = reactants
        self.product_elements = products

    def _set_matrix(self) -> None:
        """
        Get the matrix from the molecular formulas
        """
        if self.elements is None:
            raise UninitializedReactionError("Uninitialized: prepare_state() method was not called.")
        if self.reactant_elements is None:
            raise UninitializedReactionError("Uninitialized: prepare_state() method was not called.")
        if self.product_elements is None:
            raise UninitializedReactionError("Uninitialized: prepare_state() method was not called.")

        mat = []
        for el in self.elements:
            row = []
            for d in self.reactant_elements:
                num = d[el] if el in d else 0
                row.append(num)
            for d in self.product_elements:
                num = d[el] if el in d else 0
                row.append(-num)
            mat.append(row)
        self.matrix = numpy.array(mat)
        self.rank = numpy.linalg.matrix_rank(self.matrix)

    @staticmethod
    def _elements_from_formula(formula: str) -> Dict[str, int]:
        """
        Get all elements from the formula
        """
        # Find the capital letters
        letter_indices = [i for i, s in enumerate(formula) if not s.isdigit()]
        letters = [formula[i] for i in letter_indices]
        capitals = [s == s.upper() for s in letters]

        # Use them to get the elements
        lower = [i for i, b in enumerate(capitals) if b]
        upper = lower[1:] + [len(formula)]
        elements = ["".join(letters[i:j]) for i, j in zip(lower, upper)]

        # Get the correponding numbers from in between the elements
        lengths = [len(el) for el in elements]
        indices = [letter_indices[i] for i, b in enumerate(capitals) if b]
        lower = [i + l for i, l in zip(indices, lengths)]
        upper = indices[1:] + [len(formula)]
        str_numbers = [formula[l:u] for l, u in zip(lower, upper)]
        numbers = [int(s) if s.isdigit() else 1 for s in str_numbers]

        # Create a dictionary
        element_numbers = {el: i for el, i in zip(elements, numbers)}
        return element_numbers


@dataclass
class Species:
    formula: str
    charge: int = 0
    min_coeff: int = 0
    smiles: Optional[str] = None

    @classmethod
    def from_molecule(cls, mol: Molecule, min_coeff: int = 0, smiles: Optional[str] = None) -> "Species":
        """Initialize a Species from a PLAMS Molecule."""
        return cls(
            formula=mol.get_formula(as_dict=False),
            charge=int(mol.properties.get("charge", 0) or 0),
            min_coeff=min_coeff,
            smiles=smiles,
        )

    @classmethod
    def from_smiles(cls, smiles: str, min_coeff: int = 0) -> "Species":
        """Initialize a Species from a SMILES string."""
        from scm.plams import from_smiles as plams_from_smiles

        return cls.from_molecule(plams_from_smiles(smiles), min_coeff=min_coeff, smiles=smiles)

    @classmethod
    def from_dict(
        cls, d: Mapping[str, int], charge: int = 0, min_coeff: int = 0, smiles: Optional[str] = None
    ) -> "Species":
        """Initialize a Species from a dictionary.

        Example:

        >>> x = Species.from_dict({"H": 4, "C": 1}, min_coeff=0)
        """

        formula = "".join(f"{k}{v}" for k, v in d.items())
        return cls(formula=formula, charge=charge, min_coeff=min_coeff, smiles=smiles)


def balance(
    reactants: Sequence[Union[str, Species]],
    products: Sequence[Union[str, Species]],
    method: Literal["plams", "sympy"] = "plams",
) -> ReactionEquation:
    """Method to find balanced reaction converting reactants into products.

    reactants: sequence of str | Species
        The reactants. If string, should be the summary formula like "C2H6O". Note: an element may only appear once in the string, "HCOOH" is not allowed.

        To set charges or min_coeff, use a Species.

    products: sequence of str | Species
        The products.

    method: "plams" or "sympy" (default "plams")
        Which method to use. "sympy" requires that sympy is installed.

    Returns: ReactionEquation
        A balanced reaction equation.

    Example:

    >>> from scm.plams.tools.reaction import balance, Species
    >>> balance(reactants=["CH4", "O2"], products=["CO2", "H2O"])
    >>> balance(reactants=[Species("H", charge=1), Species("OH", charge=-1)], products=[Species("H2O")])

    """

    def to_species(x: Any) -> Species:
        if isinstance(x, str):
            return Species(formula=x)
        if isinstance(x, Species):
            return x
        raise ValueError(f"Cannot handle {x} of type {type(x)}")

    r = [to_species(x) for x in reactants]
    p = [to_species(x) for x in products]

    min_coeffs: Optional[numpy.typing.NDArray] = np.array([x.min_coeff for x in r] + [x.min_coeff for x in p])
    if min_coeffs is not None and all(x == 0 for x in min_coeffs):
        min_coeffs = None
    ret = ReactionEquation([x.formula for x in r], [x.formula for x in p])
    ret.method = method
    if all(x.smiles for x in r) and all(x.smiles for x in p):
        ret.rsmiles = [x.smiles for x in r]
        ret.psmiles = [x.smiles for x in p]
    ret.set_charges([x.charge for x in r], [x.charge for x in p])
    ret.balance(min_coeffs=min_coeffs)
    return ret

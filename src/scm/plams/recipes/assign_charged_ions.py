import numpy
from scm.plams import PlamsError
from scm.plams import PT


def assign_charged_ions(mol, use_system_charges=True):
    """
    Identify the charged submolecules (ions) and assign their charges

    Note: The results are written in the atom properties (region and analysis.charge)
    """
    mol = mol.copy()
    ions, ioncharges = extract_charged_ions(mol, use_system_charges)

    for i, (k, q) in enumerate(ioncharges.items()):
        for i_ion, indices in enumerate(ions[k]):
            name = "ion%i_%i" % (i + 1, i_ion)
            for iat in indices:
                mol.atoms[iat].properties.region = name
                if not use_system_charges:
                    mol.atoms[iat].properties.charge = 0.0
            if not use_system_charges:
                mol.atoms[indices[0]].properties.analysis.charge = q
    return mol


def extract_charged_ions(mol, use_system_charges=True):
    """
    Main functionality (charge assignment to ions)
    """
    # Assign atomic charges
    if use_system_charges:
        charges = get_system_charges(mol)
    else:
        charges = mol.guess_atomic_charges(adjust_to_systemcharge=False)
        charges = add_metal_charges(mol, charges)

    # Find ions
    molindices = mol.get_molecule_indices()
    ions = {}
    ioncharges = {}
    nions = {}
    formulas = {}
    for imol, atoms in enumerate(molindices):
        submol = mol.get_fragment(atoms)
        label = submol.label()
        q = sum([charges[i] for i in atoms])
        if abs(q) > 1e-10:
            if not label in ions.keys():
                ions[label] = []
                ioncharges[label] = q
                nions[label] = 0
                formulas[label] = submol.get_formula()
            ions[label].append(atoms)
            nions[label] += 1

    return ions, ioncharges


def get_system_charges(mol):
    """
    Read the atomic charges from the system
    """
    charges = []
    for at in mol.atoms:
        charge = None
        for k in at.properties.keys():
            if "charge" in at.properties[k]:
                charge = at.properties[k].charge
        if charge is None:
            raise PlamsError("Not all atoms are assigned a charge")
        charges.append(charge)
    return charges


def add_metal_charges(mol, charges):
    """
    Assign charges to the metal ions using a simple EEQ scheme
    """
    if "charge" in mol.properties:
        system_charge = mol.properties.charge
    else:
        system_charge = 0.0
    if system_charge == sum(charges):
        return charges
    dq = system_charge - sum(charges)

    elements = [at.symbol for at in mol.atoms]
    metals = [i for i, at in enumerate(mol.atoms) if at.is_metallic]

    metal_elements = set([elements[iat] for iat in metals])
    electron_affinities = {el: PT.get_electron_affinity(el) for el in metal_elements}
    ionization_energies = {el: PT.get_ionization_energy(el) for el in metal_elements}
    chi0 = {iat: 0.5 * (electron_affinities[elements[iat]] + ionization_energies[elements[iat]]) for iat in metals}
    J0 = {iat: ionization_energies[elements[iat]] - electron_affinities[elements[iat]] for iat in metals}

    # Set up the matrix (A) and vector b, to then solve Ax = b
    # The last element of the vector (and last row of A) hold the constraint sum(x_i) = Q_tot
    matrix = numpy.zeros((len(metals) + 1, len(metals) + 1))
    b = numpy.ones(len(metals) + 1)
    for i, iat in enumerate(metals):
        matrix[i, i] = J0[iat]
        b[i] = chi0[iat]
    matrix[-1, :-1] = 1.0
    matrix[:-1, -1] = 1.0
    b[-1] = dq

    # Solve Ax+b
    mcharges = numpy.linalg.solve(matrix, b)
    mcharges_int = [round(q) for q in mcharges]
    if sum(mcharges_int) != dq:
        # print([(elements[iat], mcharges[i]) for i, iat in enumerate(metals)])
        raise Exception("Predicted charges non-integer!")
    for i, iat in enumerate(metals):
        charges[iat] = mcharges_int[i]
    return charges

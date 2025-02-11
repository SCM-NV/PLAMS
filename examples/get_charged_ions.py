import sys
import numpy
from scm.plams import AMSJob, Settings, init


def main(filename):
    """
    Main functionality (charge assignment to ions)
    """
    # this line is not required in AMS2025+
    init()
    
    # Get the molecular system
    job = AMSJob.load_external(filename)
    mol = job.molecule
    elements = [at.symbol for at in mol.atoms]

    # Assign atomic charges
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
            ions[label] += atoms
            nions[label] += 1
    for k, q in ioncharges.items():
        print(formulas[k], q, nions[k])

    return ions, ioncharges


def add_metal_charges(mol, charges):
    """
    Assign charges to the metal ions
    """
    system_charge = mol.properties.charge if not isinstance(mol.properties.charge, Settings) else 0.0
    if system_charge == sum(charges):
        return charges
    dq = system_charge - sum(charges)
    # print ('Charge difference: ',dq)

    elements = [at.symbol for at in mol.atoms]
    metals = [i for i, at in enumerate(mol.atoms) if at.is_metallic]
    # print ('Metals: ',metals)

    electron_affinities, ionization_energies = read_ea_ie()
    chi0 = {iat: 0.5 * (electron_affinities[elements[iat]] + ionization_energies[elements[iat]]) for iat in metals}
    J0 = {iat: ionization_energies[elements[iat]] - electron_affinities[elements[iat]] for iat in metals}
    Jc = 1.0

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
    # print ('Matrix:')
    # print (matrix)
    # print ('b: ',b)

    # Solve Ax+b
    mcharges = numpy.linalg.solve(matrix, b)
    mcharges_int = [round(q) for q in mcharges]
    if sum(mcharges_int) != dq:
        print([(elements[iat], mcharges[i]) for i, iat in enumerate(metals)])
        raise Exception("Predicted charges non-integer!")
    for i, iat in enumerate(metals):
        charges[iat] = mcharges_int[i]
    return charges


def read_ea_ie():
    """
    Read in the electron affinities and ionization energies from file
    """
    electron_affinities = {}
    ionization_energies = {}
    for el, l in data.items():
        electron_affinities[el] = data[el][0]
        ionization_energies[el] = data[el][1]
    return electron_affinities, ionization_energies


# Element:[     ElectronAffinity,     IonizationEnergy]
data = {
    "Li": [0.0227110811, 0.1981523452],
    "Na": [0.0201386286, 0.1888547667],
    "Al": [0.0162064511, 0.2199814425],
    "Si": [0.0508978112, 0.2995804744],
    "Sc": [0.0069088726, 0.2411123028],
    "Ti": [0.0029031965, 0.2509243718],
    "V": [0.0192933941, 0.2479109274],
    "Cr": [0.0244750486, 0.2486826632],
    "Fe": [0.0059901395, 0.2903931438],
    "Co": [0.0242913020, 0.2896214081],
    "Ni": [0.0424822164, 0.2807648214],
    "Cu": [0.0451281676, 0.2839252631],
    "Ga": [0.0110247967, 0.2204591837],
    "Ge": [0.0496115849, 0.2903196452],
    "Rb": [0.0171986828, 0.1535019187],
    "Y": [0.0112820419, 0.2284705360],
    "Zr": [0.0156552112, 0.2437950033],
    "Nb": [0.0328171447, 0.2483886686],
    "Mo": [0.0274149943, 0.2606261929],
    "Tc": [0.0202121272, 0.2675350654],
    "Ru": [0.0385867883, 0.2705117605],
    "Rh": [0.0417839793, 0.2741131941],
    "Pd": [0.0204693725, 0.3063790990],
    "Ag": [0.0478476175, 0.2784128648],
    "In": [0.0110247967, 0.2126315781],
    "Sn": [0.0440991866, 0.2698870221],
    "Sb": [0.0393217747, 0.3175141436],
    "Cs": [0.0173456801, 0.1431018606],
    "La": [0.0183746611, 0.2049509698],
    "Ce": [0.0183746611, 0.2035544955],
    "Ta": [0.0118332817, 0.2899521520],
    "W": [0.0299506976, 0.2932595910],
    "Re": [0.0055123983, 0.2895846587],
    "Os": [0.0404242544, 0.3197191029],
    "Ir": [0.0575126892, 0.3344188318],
    "Au": [0.0848541849, 0.3390492464],
    "Tl": [0.0073498644, 0.2244648598],
    "Pb": [0.0132297560, 0.2725697226],
    "Bi": [0.0347648588, 0.2678658093],
    "Po": [0.0698237121, 0.3093190448],
}

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage: amspython get_chargee_ions.py path/to/ams.rkf")
        sys.exit(0)
    filename = sys.argv[1]

    ions, ioncharges = main(filename)

    outfile = open("charges.in", "w")
    for k, q in ioncharges.items():
        outfile.write("%8.1f: %s\n" % (q, " ".join(["%i" % (iat) for iat in ions[k]])))
    outfile.close()

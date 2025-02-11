import sys
from io import StringIO
import numpy
from scm.plams import RKFTrajectoryFile, init
from scm.flexmd import PDBMolecule


def main():
    """
    Main sctipt
    """
    # this line is not required in AMS2025+
    init()

    if sys.argv[1] == "-h":
        raise SystemExit("amspython get_hbonds.py path/to/ams.rkf path/to/indices.txt")

    # Deal with the input
    filename = sys.argv[1]
    indexfile = sys.argv[2]
    indices, elements = read_indices(indexfile)

    # Open the trajectory file
    rkf = RKFTrajectoryFile(filename)
    mol = rkf.get_plamsmol()
    nsteps = len(rkf)
    print("NSteps: ", nsteps)
    for i, at in enumerate(mol.atoms):
        if at.symbol in elements:
            if not i in indices:
                indices.append(i)

    # Create the PDB object
    print("Creating the PDB object")
    f = StringIO()
    mol.writexyz(f)
    f.seek(0)
    text = f.read()
    f.close()
    pdb = PDBMolecule(xyzstring=text, bonddetectorflag=False)
    print("PDB created")

    # Get the actual number of HBonds per selected atom
    heavy_atoms = [i for i, at in enumerate(mol.atoms) if not at.symbol == "H"]
    hydrogens = [i for i, at in enumerate(mol.atoms) if at.symbol == "H"]
    values = []
    print("%8s %8s %s" % ("Step", "Atom", "Neighbors"))
    for istep in range(nsteps):
        crd, cell = rkf.read_frame(istep, molecule=mol)
        # Create neighborlists
        d_indices, boxlist = pdb.divide_into_cubes(range(len(mol)))
        pdb.coords = crd
        pdb.set_cellvectors(cell)
        for iat in indices:
            atomlists = (heavy_atoms, hydrogens)
            atoms, hs = pdb.find_neighbours_using_cubes(iat, d_indices, boxlist, atomlists)
            hbonds = pdb.get_hbonds(iat, atoms, hs)
            print("%8i %8i %s" % (istep, iat, str(hbonds)))
            values.append(len(hbonds))

    # Compute the histogram
    bins = [i for i in range(max(values) + 1)]
    yvalues, xvalues = numpy.histogram(values, bins=bins)

    # Write to output
    outfile = open("hist.txt", "w")
    for x, y in zip(xvalues, yvalues):
        outfile.write("%20.10f %20.10f\n" % (x, y / nsteps))
    outfile.close()


def read_indices(indexfilename):
    """
    Read atom indices from a file
    """
    infile = open(indexfilename)
    lines = infile.readlines()
    infile.close()
    indices = []
    elements = set()
    for line in lines:
        words = line.split()
        if len(words) == 0:
            continue
        digits = [w.isdigit() for w in words]
        for w in words:
            if w.isdigit():
                indices.append(int(w))
            else:
                elements.add(w)
    return indices, elements


if __name__ == "__main__":
    main()

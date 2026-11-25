import sys
from scm.plams import RKFTrajectoryFile, DCDTrajectoryFile
from scm.flexmd import pdb_from_plamsmol, PSFTopology
from pathlib import Path


def main():
    """
    Main script
    """
    args = sys.argv
    rkf_file = Path(args[1] if len(args) > 1 else "ams.rkf")
    if not rkf_file.is_file():
        raise ValueError(f"""Could not find file '{rkf_file}'. 
To use this script, specify the path to an rkf file as an argument e.g.

amspython rkf_to_wrapped_dcd.py my_dir/ams.rkf

otherwise './ams.rkf' will be used.
""")

    rkf = RKFTrajectoryFile(rkf_file)
    mol = rkf.get_plamsmol()
    print("NSteps: ", len(rkf))

    pdb = pdb_from_plamsmol(mol)
    psf = PSFTopology(pdb=pdb)
    psf.write_psf("ams.psf")

    dcd = DCDTrajectoryFile("ams.dcd", mode="wb")

    for i in range(len(rkf)):
        if i % 100 == 0:
            print(i)
        crd, cell = rkf.read_frame(i)
        mol.from_array(crd)
        mol.map_atoms_to_bonds()
        dcd.write_next(coords=mol.as_array(), cell=cell)

    print("Created files 'ams.dcd' and 'ams.psf'")


if __name__ == "__main__":
    main()

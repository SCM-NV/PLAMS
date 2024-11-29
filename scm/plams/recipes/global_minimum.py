__all__ = ["global_minimum"]

import sys

try:
    from rdkit.Chem import AllChem, rdForceFieldHelpers
    from scm.plams.interfaces.molecule.rdkit import from_rdmol, to_rdmol
except ImportError:
    pass

from scm.plams.core.functions import finish, init


def global_minimum(mol, n_scans=1, no_h=True, no_ring=True, bond_orders=[1.0], job_type=False, path=".", **kwarg):
    """
    Find the global minimum of the ligand (RDKit UFF or user-defined PLAMS |Job|) by systematically varying dihedral angles within the molecule.

    :param |Molecule| mol: The input molecule
    :param int n_scans:: The number of times the global minimum search should be repeated
    :param bool no_h: If hydrogen-containing bonds should ignored
    :param bool no_ring: If bonds in ring systems should ignored
    :param list bond_orders: A list of accepted bond orders (floats); a bond will be ignored if its bond order is not in *bond_orders*
    :param type or bool job_type: A type object of a class derived from |Job|. Set to ''False'' to use RDKit UFF
    :param str path: The path where the PLAMS working directory will be stored
    :param dict kwarg: Keyword arguments for *job_type*
    :return |Molecule|: A copy of *mol* with a newly optimized geometry
    """
    # Creates guess bonds if no bonds are present
    if len(mol.bonds) == 0:
        mol.guess_bonds()

    # Create a list of 2-tuples (i.e. atomic indices) representing (valid) bonds within the molecule
    bond_list = find_bond(mol, no_h=no_h, no_ring=no_ring, bond_orders=bond_orders)

    # Search for the global minimum with RDKit UFF or with PLAMS at an user-defined level of theory
    if not job_type:
        if "rdkit.Chem.AllChem" not in sys.modules or "rdkit.Chem.rdForceFieldHelpers" not in sys.modules:
            raise ImportError("rdkit.chem module not found, aborting RDKit UFF optimization")
        if not no_ring:
            raise TypeError("no_ring=False is not supported in combination with RDKit UFF")
        if not rdForceFieldHelpers.UFFHasAllMoleculeParams(to_rdmol(mol)):
            raise ValueError("No UFF parameters found for one or more atoms")

        for i in range(n_scans):
            for bond in bond_list:
                mol = global_minimum_scan_rdkit(mol, bond)

        # Optimize the molecule even if no valid bonds are found
        if not bond_list:
            rdmol = to_rdmol(mol)
            AllChem.UFFGetMoleculeForceField(rdmol).Minimize()
            mol = from_rdmol(rdmol)

    else:
        init(path=path)
        for i in range(n_scans):
            for bond in bond_list:
                mol = global_minimum_scan_plams(mol, bond, job_type, **kwarg)

        # Optimize the molecule even if no valid bonds are found
        if not bond_list:
            job = job_type(molecule=mol, **kwarg)
            results = job.run()
            mol = results.get_main_molecule()
        finish()

    return mol


def find_bond(mol, no_h=True, no_ring=True, bond_orders=[1.0]):
    """
    Create a list of bonds. Each entry is a tuple with indices of atoms forming a dihedral.
    Consider only diherdals with axis being a single bond, so that rotation is possible.

    :param |Molecule| mol: The input molecule
    :param bool no_h: If hydrogen-containing bonds should ignored
    :param bool no_ring: If bonds in ring systems should ignored
    :param list bond_orders: A list of accepted bond orders (floats); A bond will be ignored if its bond order is not in *bond_orders*
    :return list: A list of 2-tuples containing the atomic indices of valid bonds
    """
    mol.set_atoms_id()

    # Mark atoms that can form an "axis" of a diherdal, i.e atoms with more than one (non-hydrogen) neighbor
    for atom in mol:
        neighbors = mol.neighbors(atom)
        if no_h:
            neighbors = [at for at in mol.neighbors(atom) if at.atnum != 1]
        if no_ring:
            neighbors = [mol.in_ring(at) for at in neighbors]
        atom.mark = len(neighbors) > 1

    # For each bond with both ends marked add one bond to the list
    ret = []
    for i, bond in enumerate(mol.bonds):
        if bond.atom1.mark and bond.atom2.mark and bond.order in bond_orders:
            at1, at2 = bond.atom1, bond.atom2
            ret.append((at1.id - 1 + 1, at2.id - 1 + 1))

    # Clean up the molecule
    mol.unset_atoms_id()
    for atom in mol:
        del atom.mark

    return ret


def global_minimum_scan_plams(mol, bond_tuple, job_type, **kwarg):
    """
    Optimize the molecule (A PLAMS |Job|) with 3 different values for the given dihedral angle and find the lowest energy conformer.
    The matching PLAMS |Results| object must have access to the |get_energy()| and |get_main_molecule()| functions.
    If required, functions can be added manually to a class with the |add_to_class()| function.

    :param |Molecule| mol: The input molecule
    :param tuple bond_tuple: A 2-tuple containing the atomic indices of valid bonds
    :param type job_type: A type object of a class derived from |Job|
    :param dict kwarg: Keyword arguments for *job_type*
    :return |Molecule|: A copy of *mol* with a newly optimized geometry
    """
    # Define a number of variables and create 3 copies of the ligand
    angles = (-120, 0, 120)
    mol_list = [mol.copy() for i in range(3)]
    for angle, mol in zip(angles, mol_list):
        bond = mol[bond_tuple]
        atom = mol[bond_tuple[0]]
        mol.rotate_bond(bond, atom, angle, unit="degree")

    # Optimize the geometry for all dihedral angles in angle_list
    # The geometry that yields the minimum energy is returned
    energy_list = []
    for mol in mol_list:
        job = job_type(**kwarg)
        job.molecule = mol
        results = job.run()
        energy_list.append(results.get_energy())
        mol_new = results.get_main_molecule()
        for at, at_new in zip(mol, mol_new):
            at.coords = at_new.coords
    minimum = energy_list.index(min(energy_list))
    return mol_list[minimum]


def global_minimum_scan_rdkit(mol, bond_tuple):
    """
    Optimize the molecule (RDKit UFF) with 3 different values for the given dihedral angle and find the lowest energy conformer.

    :param |Molecule| mol: The input molecule
    :param tuple bond_tuple: A 2-tuples containing the atomic indices of valid bonds
    :return |Molecule|: A copy of *mol* with a newly optimized geometry
    """
    # Define a number of variables and create 3 copies of the ligand
    uff = AllChem.UFFGetMoleculeForceField
    angles = (-120, 0, 120)
    mol_list = [mol.copy() for i in range(3)]
    for angle, mol in zip(angles, mol_list):
        bond = mol[bond_tuple]
        atom = mol[bond_tuple[0]]
        mol.rotate_bond(bond, atom, angle, unit="degree")

    # Optimize the geometry for all dihedral angles in angle_list
    # The geometry that yields the minimum energy is returned
    mol_list = [to_rdmol(mol, properties=False) for mol in mol_list]
    for rdmol in mol_list:
        uff(rdmol).Minimize()
    energy_list = [uff(rdmol).CalcEnergy() for rdmol in mol_list]
    minimum = energy_list.index(min(energy_list))
    return from_rdmol(mol_list[minimum])

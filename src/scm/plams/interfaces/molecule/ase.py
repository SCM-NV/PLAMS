import numpy as np
from typing import Optional, TYPE_CHECKING, List

from scm.plams.core.functions import requires_optional_package
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Atom, Molecule

if TYPE_CHECKING:
    from ase.atoms import Atoms as ASEAtoms


__all__ = ["toASE", "fromASE"]


@requires_optional_package("ase")
def toASE(molecule: Molecule, set_atomic_charges: bool = False) -> "ASEAtoms":
    """Convert a PLAMS |Molecule| to an ASE molecule (``ase.Atoms`` instance). Translate coordinates, atomic numbers, and lattice vectors (if present). The order of atoms is preserved.


    set_atomic_charges: bool
        If True, set_initial_charges() will be called with the average atomic charge (taken from molecule.properties.charge). The purpose is to preserve the total charge, not to set any reasonable initial charges.
    """
    import ase

    # iterate over PLAMS atoms
    for atom in molecule:

        # check if coords only consists of floats or ints
        if not all(isinstance(x, (int, float)) for x in atom.coords):
            raise ValueError("Non-Number in Atomic Coordinates, not compatible with ASE")

    ase_mol = ase.Atoms(numbers=molecule.numbers, positions=molecule.as_array())

    # get lattice info if any
    lattice = np.zeros((3, 3))
    pbc = [False, False, False]
    for i, vec in enumerate(molecule.lattice):

        # check if lattice only consists of floats or ints
        if not all(isinstance(x, (int, float)) for x in vec):
            raise ValueError("Non-Number in Lattice Vectors, not compatible with ASE")

        pbc[i] = True
        lattice[i] = np.array(vec)

    # save lattice info to ase_mol
    if any(pbc):
        ase_mol.set_pbc(pbc)
        ase_mol.set_cell(lattice)

    if set_atomic_charges:
        charge = molecule.properties.get("charge", 0)
        if not charge:
            atomic_charges = [0.0] * len(molecule)
        else:
            atomic_charges = [float(charge)] + [0.0] * (len(molecule) - 1)

        ase_mol.set_initial_charges(atomic_charges)

    return ase_mol


def fromASE(molecule: "ASEAtoms", properties: Optional[Settings] = None, set_charge: bool = False) -> Molecule:
    """Convert an ASE molecule to a PLAMS |Molecule|. Translate coordinates, atomic numbers, and lattice vectors (if present). The order of atoms is preserved.

    Pass a |Settings| instance through the ``properties`` option to inherit them to the returned molecule.
    """
    plams_mol = Molecule()

    # iterate over ASE atoms
    for atom in molecule:
        # add atom to plams_mol
        plams_mol.add_atom(Atom(atnum=atom.number, coords=tuple(atom.position)))

    # add Lattice if any
    if any(molecule.get_pbc()):
        lattice: List[List[float]] = []
        # loop over three booleans
        for i, boolean in enumerate(molecule.get_pbc().tolist()):
            if boolean:
                lattice.append(list(molecule.get_cell()[i]))

        # write lattice to plams_mol
        plams_mol.lattice = lattice.copy()

    if properties:
        plams_mol.properties.update(properties)
    if (properties and "charge" not in properties or not properties) and set_charge:
        plams_mol.properties.charge = sum(molecule.get_initial_charges())
        if "charge" in molecule.info:
            plams_mol.properties.charge += molecule.info["charge"]
    return plams_mol

import numpy as np
from typing import Optional, TYPE_CHECKING

from scm.plams.core.functions import add_to_class, requires_optional_package
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Atom, Molecule

try:
    import ase

    _has_ase = True
except ImportError:

    _has_ase = False

if TYPE_CHECKING:
    from ase import atoms as ASEAtoms


__all__ = ["toASE", "fromASE"]


@add_to_class(Molecule)
@requires_optional_package("ase")
def readase(self, f, **other):
    """Read Molecule using ASE engine

    The ``read`` function of the |Molecule| class passes a file descriptor into here, so in this case you must specify the *format* to be read by ASE::

        mol = Molecule('file.cif', inputformat='ase', format='cif')

    The ASE Atoms object then gets converted to a PLAMS Molecule and returned.
    All *other* options are passed to ``ASE.io.read()``.
    See https://wiki.fysik.dtu.dk/ase/ase/io/io.html on how to use it.

    .. note::

        The nomenclature of PLAMS and ASE is incompatible for reading multiple geometries, make sure that you only read single geometries with ASE! Reading multiple geometries is not supported, each geometry needs to be read individually.

    """
    from ase import io

    ase_mol = io.read(f, **other)
    mol = fromASE(ase_mol)
    # update self with the molecule read without overwriting e.g. settings
    self += mol
    # lattice does not survive soft update
    self.lattice = mol.lattice
    return


@add_to_class(Molecule)
def writease(self, f, **other):
    """Write molecular coordinates using ASE engine.

    The ``write`` function of the |Molecule| class passes a file descriptor into here, so in this case you must specify the *format* to be written by ASE.
    All *other* options are passed to ``ASE.io.write()``.
    See https://wiki.fysik.dtu.dk/ase/ase/io/io.html on how to use it.

    These two write the same content to the respective files::

        molecule.write('filename.anyextension', outputformat='ase', format='gen')
        molecule.writease('filename.anyextension', format='gen')

    """
    ase_mol = toASE(self)
    ase_mol.write(f, **other)
    return


if _has_ase:
    Molecule._readformat["ase"] = Molecule.readase
    Molecule._writeformat["ase"] = Molecule.writease


@requires_optional_package("ase")
def toASE(molecule: Molecule, set_atomic_charges: bool = False) -> "ASEAtoms":
    """Convert a PLAMS |Molecule| to an ASE molecule (``ase.Atoms`` instance). Translate coordinates, atomic numbers, and lattice vectors (if present). The order of atoms is preserved.


    set_atomic_charges: bool
        If True, set_initial_charges() will be called with the average atomic charge (taken from molecule.properties.charge). The purpose is to preserve the total charge, not to set any reasonable initial charges.
    """

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
        lattice = []
        # loop over three booleans
        for i, boolean in enumerate(molecule.get_pbc().tolist()):
            if boolean:
                lattice.append(tuple(molecule.get_cell()[i]))

        # write lattice to plams_mol
        plams_mol.lattice = lattice.copy()

    if properties:
        plams_mol.properties.update(properties)
    if (properties and "charge" not in properties or not properties) and set_charge:
        plams_mol.properties.charge = sum(molecule.get_initial_charges())
        if "charge" in molecule.info:
            plams_mol.properties.charge += molecule.info["charge"]
    return plams_mol

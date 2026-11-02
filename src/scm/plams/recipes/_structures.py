from typing import TYPE_CHECKING, Union

from scm.plams.mol.molecule import Molecule

if TYPE_CHECKING:
    from scm.base import ChemicalSystem

Structure = Union[Molecule, "ChemicalSystem"]


def is_chemical_system(structure: object) -> bool:
    """Return whether *structure* is a ChemicalSystem when scm.base is available."""
    try:
        from scm.base import ChemicalSystem
    except ImportError:
        return False
    return isinstance(structure, ChemicalSystem)

from typing import Any, TYPE_CHECKING, Union

from scm.plams.core.basejob import Job
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.mol.atom import Atom
from scm.plams.mol.molecule import Molecule

if TYPE_CHECKING:
    from scm.base import ChemicalSystem
    from scm.inputs import AMS as AMSInputModel

RecipeSettings = Union[Settings, Job, "AMSInputModel"]
Structure = Union[Molecule, "ChemicalSystem"]

try:
    from scm.inputs import InputModel

    _has_scm_inputs = True
except ImportError:
    _has_scm_inputs = False


def to_plams_settings(settings: Any, recipe_name: str) -> Any:
    """Convert a complete typed AMS model to the mutable form used by legacy recipes."""
    source_settings = settings.settings if isinstance(settings, Job) else settings
    if not _has_scm_inputs:
        return source_settings

    input_model = source_settings.get("input") if isinstance(source_settings, Settings) else source_settings
    if not isinstance(input_model, InputModel):
        return source_settings

    # Loading the generated AMS model is expensive; keep ordinary PLAMS imports lightweight.
    from scm.inputs import AMS as AMSInputModel

    if not isinstance(input_model, AMSInputModel):
        raise TypeError(f"{recipe_name} settings must be PLAMS Settings or scm.inputs.AMS")

    from scm.plams.interfaces.adfsuite.inputparser import input_to_settings

    converted = source_settings.copy() if isinstance(source_settings, Settings) else Settings()
    converted.input = input_to_settings(input_model.to_input(), program=AMSJob._command)
    systems = converted.input.ams.get("System")
    # Recipes access the single unnamed System as a Settings block, not a repeated entry.
    if isinstance(systems, list) and len(systems) == 1 and not systems[0].get("_h", ""):
        converted.input.ams.System = systems[0]
    return converted


def is_chemical_system(structure: object) -> bool:
    """Return whether *structure* is a ChemicalSystem when scm.base is available."""
    try:
        from scm.base import ChemicalSystem
    except ImportError:
        return False
    return isinstance(structure, ChemicalSystem)


def get_main_structure(results: Any, use_chemical_system: bool) -> Any:
    """Read a result structure in the same representation as the recipe input."""
    return results.get_main_system() if use_chemical_system else results.get_main_molecule()


def get_structure_charge(structure: Structure) -> float:
    """Return the charge from either supported structure representation."""
    if isinstance(structure, Molecule):
        return structure.properties.get("charge", 0)
    return structure.charge


def get_atomic_number(atom: Any) -> int:
    """Return an atom's atomic number from either supported representation."""
    return atom.atnum if isinstance(atom, Atom) else atom.element.Z

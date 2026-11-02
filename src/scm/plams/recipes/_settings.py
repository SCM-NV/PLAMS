from typing import Any, TYPE_CHECKING, Union

from scm.plams.core.basejob import Job
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob

if TYPE_CHECKING:
    from scm.inputs import AMS as AMSInputModel

RecipeSettings = Union[Settings, Job, "AMSInputModel"]

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
    return converted

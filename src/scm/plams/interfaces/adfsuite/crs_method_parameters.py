"""Apply CRS parameter presets to PLAMS Settings or typed CRS inputs."""

from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Tuple

from scm.plams.core.functions import requires_optional_package
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.crs_definitions import get_crs_method_parameters_metadata

if TYPE_CHECKING:
    from scm.inputs import CRS

__all__ = ["apply_parameter_set_to_settings", "apply_parameter_set_to_inputs", "get_parameter_set_options"]

_PARAMETER_BLOCKS = ("CRSPARAMETERS", "SACPARAMETERS", "DISPERSION", "EPSILON")


def get_parameter_set_options() -> Dict[str, Tuple[str, ...]]:
    """Return available preset names grouped by method."""
    options: Dict[str, List[str]] = {}
    for parameter_set, preset in get_crs_method_parameters_metadata().items():
        options.setdefault(preset["method"], []).append(parameter_set)
    return {method: tuple(names) for method, names in options.items()}


def _get_preset(parameter_set: str) -> Mapping[str, Any]:
    presets = get_crs_method_parameters_metadata()
    key = parameter_set.strip().lower()
    if key not in presets:
        raise ValueError(
            f"Unknown parameter_set {parameter_set!r}. Available parameter sets: {', '.join(sorted(presets))}"
        )
    return deepcopy(presets[key])


def apply_parameter_set_to_settings(settings: Settings, parameter_set: str) -> Settings:
    """Apply a parameter preset in place and return the same job Settings.

    ``parameter_set`` must be one of the preset names returned by
    :meth:`get_parameter_set_options`.
    """
    if not isinstance(settings, Settings):
        raise TypeError("settings must be a PLAMS Settings object")
    if "input" in settings and not isinstance(settings["input"], Settings):
        raise TypeError("settings.input must be a PLAMS Settings object")
    preset = _get_preset(parameter_set)
    replacement = Settings()
    replacement.method = preset["method"]
    for name, parameters in preset["parameter_blocks"].items():
        if name in ("Dispersion", "Epsilon"):
            block = Settings()
            for index, (key, value) in enumerate(parameters.items(), start=1):
                block[f"_{index}"] = f"{key} {value}"
        else:
            block = Settings(parameters)
        replacement[name.lower()] = block

    for name in _PARAMETER_BLOCKS:
        settings.input.pop(name.lower(), None)
    settings.input.update(replacement)
    return settings


@requires_optional_package("scm.inputs")
def apply_parameter_set_to_inputs(crs: CRS, parameter_set: str) -> CRS:
    """Apply a parameter preset in place and return the same typed CRS model.

    ``parameter_set`` must be one of the preset names returned by
    :meth:`get_parameter_set_options`.
    """
    from scm.inputs import CRS

    if not isinstance(crs, CRS):
        raise TypeError("crs must be a scm.inputs.CRS instance")
    preset = _get_preset(parameter_set)
    values: Dict[str, Any] = {"METHOD": preset["method"]}
    for name, parameters in preset["parameter_blocks"].items():
        if name == "CRSParameters":
            values["CRSPARAMETERS"] = CRS.CRSPARAMETERSBlock(**parameters)
        elif name == "SACParameters":
            values["SACPARAMETERS"] = CRS.SACPARAMETERSBlock(**parameters)
        else:
            values[name.upper()] = [f"{key} {value}" for key, value in parameters.items()]
    replacement = CRS(**values)

    crs.METHOD = replacement.METHOD
    for name in _PARAMETER_BLOCKS:
        if name in values:
            setattr(crs, name, getattr(replacement, name))
        else:
            delattr(crs, name)
    return crs

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict


RESULT_TABLE_UNSUPPORTED_PROPERTIES = ("SIGMAPROFILE", "PURESIGMAPROFILE", "SIGMAPOTENTIAL", "PURESIGMAPOTENTIAL")
RESULT_TABLE_LLE_PROPERTIES = ("LLE", "STABILITY", "BINMIXCOEF", "TERNARYMIX")

RESULT_TABLE_COMPONENT_BASE_COLUMNS = ("property", "method", "mixture", "cid", "name",)
RESULT_TABLE_COMPONENT_DEFAULT_QUANTITIES = (
    "frac1",
    "frac2",
    "xI",
    "xII",
    "solvent fraction",
    "composition molar fraction",
    "gamma",
    "gammaI",
    "gammaII",
    "actI",
    "actII",
    "logp",
    "vapor pressure",
    "henryc",
    "henrycnodim",
    "deltag",
    "solubility molar fraction",
    "solubility massfrac",
    "solubility mol_per_L_solvent",
    "solubility g_per_L_solvent",
    "solubility mol_per_L_solution",
    "solubility g_per_L_solution",
)
RESULT_TABLE_COMPONENT_EXTRA_QUANTITIES = (
    "mu",
    "mu pure",
    "mu gas",
    "mu in solvent 1",
    "mu in solvent 2",
    "xI0",
    "xII0",
    "E gas",
    "G solute",
    "gamma_wf",
    "gamma_vf",
    "fh_chi",
    "poly fraction",
    "polyrepeats",
)
RESULT_TABLE_COMPONENT_KNOWN_QUANTITIES = tuple(
    dict.fromkeys(RESULT_TABLE_COMPONENT_DEFAULT_QUANTITIES + RESULT_TABLE_COMPONENT_EXTRA_QUANTITIES)
)
RESULT_TABLE_COMPONENT_EXCLUDED_BY_PROPERTY = {
    "LLE": ("gamma",),
    "STABILITY": ("gamma",),
    "LOGP": ("gamma",),
}
RESULT_TABLE_MIXTURE_DEFAULT_QUANTITIES = (
    "temperature",
    "pressure",
    "excess G",
    "excess H",
    "Gibbs energy of mixing",
    "Enthalpy of vaporization",
    "showmiscgap",
    "unstable",
    "converged",
    "llle_detected",
    "phiI",
    "phiII",
    "tpd_w",
    "isobar",
    "flashpoint",
)
RESULT_TABLE_MIXTURE_EXTRA_QUANTITIES = (
    "Gibbs energy",
    "status_msg",
    "xI_unstable",
    "xII_unstable",
    "tpd_w_I",
    "tpd_w_II",
    "llle_source_phase",
    "llle_status_msg",
    "conv_code",
    "L_conv_code",
    "L_status_msg",
)
RESULT_TABLE_MIXTURE_KNOWN_QUANTITIES = tuple(
    dict.fromkeys(RESULT_TABLE_MIXTURE_DEFAULT_QUANTITIES + RESULT_TABLE_MIXTURE_EXTRA_QUANTITIES)
)

RESULT_TABLE_LLE_BASE_COLUMNS = ("property", "method", "mixture", "tie_line", "cid", "name",)
RESULT_TABLE_LLE_COLUMN_SPECS = {
    "x": {"quantities": ("x",), "kind": "component"},
    "temperature": {"quantities": ("temperature", "xlle"), "kind": "mixture"},
    "xI": {"quantities": ("xI", "xlle", "xll"), "kind": "component"},
    "xII": {"quantities": ("xII", "xlle", "xll"), "kind": "component"},
    "gammaI": {"quantities": ("gammaI",), "kind": "component"},
    "gammaII": {"quantities": ("gammaII",), "kind": "component"},
    "actI": {"quantities": ("actI",), "kind": "component"},
    "actII": {"quantities": ("actII",), "kind": "component"},
    "phiI": {"quantities": ("phiI",), "kind": "mixture"},
    "phiII": {"quantities": ("phiII",), "kind": "mixture"},
    "converged": {"quantities": ("converged",), "kind": "mixture"},
    "act_interp": {"quantities": ("xlle", "actxll"), "kind": "derived"},
    "pressure": {"quantities": ("xlle", "pressure"), "kind": "derived"},
    "w_min": {"quantities": ("w_min",), "kind": "component", "extra": True},
    "status_msg": {"quantities": ("status_msg",), "kind": "mixture", "extra": True},
}
RESULT_TABLE_LLE_DEFAULT_QUANTITIES = tuple(
    dict.fromkeys(
        quantity
        for spec in RESULT_TABLE_LLE_COLUMN_SPECS.values()
        if not spec.get("extra")
        for quantity in spec["quantities"]
    )
)
RESULT_TABLE_LLE_EXTRA_QUANTITIES = tuple(
    dict.fromkeys(
        quantity
        for spec in RESULT_TABLE_LLE_COLUMN_SPECS.values()
        if spec.get("extra")
        for quantity in spec["quantities"]
    )
)
RESULT_TABLE_LLE_KNOWN_QUANTITIES = tuple(
    dict.fromkeys(RESULT_TABLE_LLE_DEFAULT_QUANTITIES + RESULT_TABLE_LLE_EXTRA_QUANTITIES)
)


def _strip_comment_period(comment: str) -> str:
    stripped = comment.strip()
    return stripped[:-1] if stripped.endswith(".") else stripped


def _extract_result_quantity_metadata_from_kf_def(data: Dict[str, Any]) -> Dict[str, Dict[str, str]]:
    result: Dict[str, Dict[str, str]] = {}

    for section in data.values():
        if not isinstance(section, dict):
            continue
        for key, value in section.items():
            if not isinstance(value, dict):
                continue

            value_type = value.get("_type")
            if value_type not in {"section", "free_section", "subsection"} and "_comment" in value:
                result[key] = {
                    "symbol": str(value.get("_symbol") or key),
                    "name": str(value.get("_gui_name") or ""), #_strip_comment_period(str(value.get("_comment", ""))) or key,
                    "comment": str(value.get("_comment") or ""),
                    "unit": str(value.get("_unit") or ""),
                }

    return result


def _load_crs_kf_def_metadata(json_path: Path) -> Dict[str, Dict[str, str]]:
    with json_path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"CRS kf definition file {json_path} must contain a JSON object")
    return _extract_result_quantity_metadata_from_kf_def(data)


def _build_result_quantity_metadata(
    kf_def_metadata: Dict[str, Dict[str, str]],
    overrides: Dict[str, Dict[str, str]],
) -> Dict[str, Dict[str, str]]:
    result = {key: value.copy() for key, value in kf_def_metadata.items()}
    for key, override in overrides.items():
        result[key] = {**result.get(key, {}), **override}
    return result

RESULT_TABLE_QUANTITY_METADATA_OVERRIDES = {
    "x": {"symbol": "x", "name": "Feed molar composition", "unit": "fraction"},
    "act_interp": {
        "symbol": "a*",
        "name": "Interpolated activity",
        "comment": "Estimated from interpolated miscibility-gap data; not from a full LLE calculation.",
    },
}

# module-level private metadata loading helpers
def _load_crs_input_block_metadata(json_path: Path) -> Dict[str, Dict[str, Any]]:
    with json_path.open(encoding="utf-8") as handle:
        return _extract_crs_input_block_metadata(json.load(handle))

def _extract_crs_input_block_metadata(data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    result: Dict[str, Dict[str, Any]] = {}

    for entry_name, entry_def in data.items():
        if not isinstance(entry_def, dict):
            continue
        if entry_name.lower() in [
            "dispersion",
            "coef_cations",
            "coef_anions",
            "epsilon",
            "crsparameters",
            "sacparameters",
        ]:
            continue

        category = entry_def.get("_category")
        if category == "key":
            result[entry_name.lower()] = _key_metadata(entry_name, entry_def)
        elif category == "block":
            result[entry_name.lower()] = _block_metadata(entry_name, entry_def)

    return result

def _key_metadata(name: str, key_def: Dict[str, Any]) -> Dict[str, Any]:
    metadata = {
        "kind": "key",
        "input_name": name,
        "type": key_def.get("_type"),
        "description": key_def.get("_comment", ""),
        "unit": key_def.get("_unit"),
        "choices": tuple(key_def.get("_choices", ())),
        "hint":  key_def.get("_hint", ""),
    }

    if "_default" in key_def:
        metadata["default"] = key_def["_default"]

    return metadata

def _block_metadata(name: str, block_def: Dict[str, Any]) -> Dict[str, Any]:
    keys: Dict[str, Any] = {}
    blocks: Dict[str, Any] = {}

    for child_name, child_def in block_def.items():
        if not isinstance(child_def, dict):
            continue

        child_category = child_def.get("_category")
        if child_category == "key":
            keys[child_name.lower()] = _key_metadata(child_name, child_def)
        elif child_category == "block":
            blocks[child_name.lower()] = _block_metadata(child_name, child_def)

    return {
        "kind": "block",
        "input_name": name,
        "description": block_def.get("_comment", ""),
        "type": block_def.get("_type"),
        "header": bool(block_def.get("_header", False)),
        "keys": keys,
        "blocks": blocks,
    }

def get_block(data: Dict[str, Any], *path: str) -> Dict[str, Any]:
    block = data[path[0].lower()]
    for name in path[1:]:
        block = block["blocks"][name.lower()]
    return block

def get_block_keys(data: Dict[str, Any], *path: str) -> Dict[str, Any]:
    return get_block(data, *path)["keys"]

def get_block_child_names(data: Dict[str, Any], *path: str) -> frozenset[str]:
    block = get_block(data, *path)
    return frozenset(block["keys"]) | frozenset(block["blocks"])

# module-level private metadata/constants used to build CRSJob
_crs_kf_def_path = Path(os.environ["AMSBIN"]) / "kf_def/crs.json"

RESULT_TABLE_QUANTITY_METADATA = _build_result_quantity_metadata(
    _load_crs_kf_def_metadata(_crs_kf_def_path.resolve()),
    RESULT_TABLE_QUANTITY_METADATA_OVERRIDES,
)

_crs_json_path = Path(os.environ["AMSBIN"]) / "../data/input_def/crs.json"

CRS_DATA = _load_crs_input_block_metadata(_crs_json_path.resolve())

_CRS_METHOD_PARAMETER_BLOCKS = frozenset(
    ("CRSParameters", "SACParameters", "Dispersion", "Epsilon")
)

def _load_crs_method_parameters_metadata(json_path: Path) -> Dict[str, Dict[str, Any]]:
    with json_path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    if data.get("schema_version") != 1:
        raise ValueError(f"Unsupported CRS method parameter schema_version in {json_path}")

    presets = data.get("presets")
    if not isinstance(presets, list):
        raise ValueError(f"CRS method parameter file {json_path} must contain a presets list")

    result: Dict[str, Dict[str, Any]] = {}

    for preset in presets:
        if not isinstance(preset, dict):
            raise ValueError(f"CRS method parameter preset in {json_path} must be an object")

        parameter_set = preset.get("parameter_set")
        if not isinstance(parameter_set, str) or not parameter_set.strip():
            raise ValueError(f"CRS method parameter preset in {json_path} must define parameter_set")

        parameter_set = parameter_set.strip().lower()
        if parameter_set in result:
            raise ValueError(f"Duplicate CRS {parameter_set=!r} in {json_path}")

        method = preset.get("method")
        if not isinstance(method, str) or not method.strip():
            raise ValueError(f"CRS {parameter_set=!r} must define method")

        parameter_blocks = preset.get("parameter_blocks")
        if not isinstance(parameter_blocks, dict):
            raise ValueError(f"CRS {parameter_set=!r} must define parameter_blocks")

        for block_name, block in parameter_blocks.items():
            if block_name not in _CRS_METHOD_PARAMETER_BLOCKS:
                allowed = ", ".join(sorted(_CRS_METHOD_PARAMETER_BLOCKS))
                raise ValueError(
                    f"Unsupported parameter block {block_name!r} in CRS {parameter_set=!r}. "
                    f"Allowed blocks: {allowed}"
                )
            if not isinstance(block, dict):
                raise ValueError(
                    f"{block_name} in CRS {parameter_set=!r} must be an object"
                )

        normalized_preset = preset.copy()
        normalized_preset["parameter_set"] = parameter_set
        normalized_preset["method"] = method
        result[parameter_set] = normalized_preset

    return result

CRS_METHODS = (
    "COSMO-RS",
    "COSMOSAC2013",
    "COSMOSAC2016",
    "COSMOSACDHB",
    "COSMOSACDHB-MESP",
)

_crs_method_parameters_path = Path(os.environ["AMSBIN"]) / "../data/crs/method_parameters.json"

if _crs_method_parameters_path.is_file():
    CRS_METHOD_PARAMETERS_METADATA = _load_crs_method_parameters_metadata(
        _crs_method_parameters_path.resolve(),
    )
else:
    CRS_METHOD_PARAMETERS_METADATA = {}

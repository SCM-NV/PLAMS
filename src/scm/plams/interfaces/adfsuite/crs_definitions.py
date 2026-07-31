from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Set, Tuple


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

def _property_metadata(
    *,
    description: str,
    system_scope: str,
    input_keys: Dict[str, Sequence[str]],
    required_keys: Sequence[str] = (),
    comments: Sequence[str] = (),
    hint_keys: Sequence[str] = (),
    builder: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Return a complete, immutable CRS property metadata dictionary."""
    normalized_input_keys = {
        "top_level": tuple(input_keys.get("top_level", ())),
        "property": tuple(input_keys.get("property", ())),
        "compound": tuple(input_keys.get("compound", ())),
    }
    return {
        "description": description,
        "system_scope": system_scope,
        "input_keys": normalized_input_keys,
        "required_keys": tuple(required_keys),
        "comments": tuple(comments),
        "hint_keys": tuple(hint_keys),
        "builder": {} if builder is None else builder,
    }


_VAPOR_PRESSURE_KEYS = ("pvap", "tvap", "vp_equation", "vp_params")
_FUSION_KEYS = ("meltingpoint", "hfusion", "cpfusion")
_VLE_SWEEP_PROPERTY_KEYS = ("nfrac", "isotherm", "isobar", "flashpoint")

_VLE_SWEEP_MODE = {
    "keys": ("isotherm", "isobar", "flashpoint"),
    "default": "isotherm",
    "set": {
        "isotherm": {"property": {"isotherm": True}},
        "isobar": {"property": {"isobar": True}},
        "flashpoint": {"property": {"flashpoint": True}},
    },
    "mode_hints": {
        "isotherm": (
            "LLE boundaries are interpolated from the composition sweep; use LLE for robust phase-boundary calculations.",
        ),
        "flashpoint": ("Uses pure-compound flashpoints; vapor-pressure inputs may improve results.",),
    },
}

_VLE_SWEEP_BUILDER = {
    "roles": ("compound",),
    "mode": _VLE_SWEEP_MODE,
}

_SOLUBILITY_MODE = {
    "keys": ("gas", "liquid", "solid"),
    "default": "solid",
    "set": {
        "gas": {"property": {"isobar": True}},
        "liquid": {},
        "solid": {},
    },
    "compound_required_keys": {
        "solid": {
            "solute": {
                "any_of": (("meltingpoint", "hfusion"),),
            },
        },
    },
    "mode_hints": {
        "gas": (
            "Top-level pressure is the solute partial pressure.",
        ),
        "liquid": (
            "Assumes the pure liquid solute as the coexisting phase; use LLE for miscible systems.",
        ),
        "solid": (
            "Requires solid-solute fusion-correction inputs.",
        ),
    },
}

_SOLUBILITY_BUILDER = {
    "roles": ("solvent", "solute"),
    "mode": _SOLUBILITY_MODE,
}

CRS_PROPERTY_TYPE_METADATA: Dict[str, Dict[str, Any]] = {
    "ACTIVITYCOEF": _property_metadata(
        description="Activity coefficients in a solvent mixture.",
        system_scope="solvent_mixture_with_solutes",
        input_keys={
            "top_level": ("temperature", "massfraction"),
            "property": ("densitysolvent",),
            "compound": ("frac1", "density") + _VAPOR_PRESSURE_KEYS,
        },
        required_keys=("temperature", "frac1"),
        hint_keys=("densitysolvent", "density"),
        builder={
            "roles": ("solvent", "solute"),
        },
    ),
    "LOGP": _property_metadata(
        description="Partition coefficients between two immiscible solvent phases.",
        system_scope="mixture",
        input_keys={
            "top_level": ("temperature", "massfraction"),
            "property": ("volumequotient",),
            "compound": ("frac1", "frac2", "density"),
        },
        required_keys=("temperature", "frac1", "frac2"),
        hint_keys=("volumequotient",),
        builder={
            "roles": ("solvent", "solute"),
        },
    ),
    "SOLUBILITY": _property_metadata(
        description="Solubility of solutes in a solvent mixture or under gas-pressure conditions.",
        system_scope="solvent_with_solutes",
        input_keys={
            "top_level": ("temperature", "pressure", "massfraction"),
            "property": ("densitysolvent", "isobar"),
            "compound": ("frac1", "density") + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS,
        },
        required_keys=("temperature", "frac1", "meltingpoint", "hfusion"),
        hint_keys=("densitysolvent", "density") + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS,
        builder={**_SOLUBILITY_BUILDER},
    ),
    "PURESOLUBILITY": _property_metadata(
        description="Solubility of a solute in pure solvents over a temperature range.",
        system_scope="pure_solvent_with_solute",
        input_keys={
            "top_level": ("temperature", "pressure"),
            "property": ("isobar",),
            "compound": ("frac1", "density") + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS,
        },
        required_keys=("temperature", "frac1", "meltingpoint", "hfusion"),
        hint_keys=("density",) + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS,
        builder={**_SOLUBILITY_BUILDER},
    ),
    "VAPORPRESSURE": _property_metadata(
        description="Vapor pressure of a mixture at fixed temperature.",
        system_scope="mixture",
        input_keys={
            "top_level": ("temperature", "massfraction"),
            "compound": ("frac1",) + _VAPOR_PRESSURE_KEYS,
        },
        required_keys=("temperature", "frac1"),
        hint_keys=("temperature",) + _VAPOR_PRESSURE_KEYS,
        builder={
            "roles": ("compound",),
        },
    ),
    "PUREVAPORPRESSURE": _property_metadata(
        description="Pure-compound vapor pressure over a temperature range.",
        system_scope="pure_compounds",
        input_keys={
            "top_level": ("temperature",),
            "compound": _VAPOR_PRESSURE_KEYS,
        },
        required_keys=("temperature",),
        hint_keys=("temperature",) + _VAPOR_PRESSURE_KEYS,
        builder={
            "roles": ("compound",),
        },
        comments=(
            "Multiple COMPOUND blocks are treated as independent pure compounds.",
        ),
    ),
    "BOILINGPOINT": _property_metadata(
        description="Boiling temperature of a mixture for a pressure range.",
        system_scope="mixture",
        input_keys={
            "top_level": ("pressure", "massfraction"),
            "compound": ("frac1",) + _VAPOR_PRESSURE_KEYS,
        },
        required_keys=("pressure",),
        hint_keys=("pressure",) + _VAPOR_PRESSURE_KEYS,
        builder={
            "roles": ("compound",),
        },
    ),
    "PUREBOILINGPOINT": _property_metadata(
        description="Pure-compound boiling point over a pressure range.",
        system_scope="pure_compounds",
        input_keys={
            "top_level": ("pressure",),
            "compound": _VAPOR_PRESSURE_KEYS,
        },
        required_keys=("pressure",),
        hint_keys=("pressure",) + _VAPOR_PRESSURE_KEYS,
        builder={
            "roles": ("compound",),
        },
        comments=(
            "Multiple COMPOUND blocks are treated as independent pure compounds.",
        ),

    ),
    "FLASHPOINT": _property_metadata(
        description="Flash point of a mixture using user-supplied pure-compound flash points.",
        system_scope="mixture",
        input_keys={
            "top_level": ("massfraction",),
            "compound": ("frac1", "flashpoint") + _VAPOR_PRESSURE_KEYS,
        },
        required_keys=("frac1",),
        hint_keys=("flashpoint",) + _VAPOR_PRESSURE_KEYS,
        builder={
            "roles": ("compound",),
        },
    ),
    "BINMIXCOEF": _property_metadata(
        description="Binary-mixture coefficients over a composition range.",
        system_scope="binary_mixture",
        input_keys={
            "top_level": ("temperature", "pressure", "massfraction"),
            "property": _VLE_SWEEP_PROPERTY_KEYS,
            "compound": ("frac1",) + _VAPOR_PRESSURE_KEYS + ("flashpoint",),
        },
        required_keys=("temperature",),
        hint_keys=("flashpoint",) + _VAPOR_PRESSURE_KEYS,
        builder={**_VLE_SWEEP_BUILDER,},
    ),
    "TERNARYMIX": _property_metadata(
        description="Ternary mixture property sweep over composition space.",
        system_scope="ternary_mixture",
        input_keys={
            "top_level": ("temperature", "pressure", "massfraction"),
            "property": _VLE_SWEEP_PROPERTY_KEYS,
            "compound": ("frac1",) + _VAPOR_PRESSURE_KEYS + ("flashpoint",),
        },
        required_keys=("temperature",),
        hint_keys=("flashpoint",) + _VAPOR_PRESSURE_KEYS,
        builder={**_VLE_SWEEP_BUILDER,},
    ),
    "COMPOSITIONLINE": _property_metadata(
        description="Composition-line calculation between two endpoint phase compositions.",
        system_scope="binary_mixture",
        input_keys={
            "top_level": ("temperature", "pressure", "massfraction"),
            "property": _VLE_SWEEP_PROPERTY_KEYS,
            "compound": ("frac1", "frac2") + _VAPOR_PRESSURE_KEYS + ("flashpoint",),
        },
        required_keys=("temperature",),
        hint_keys=("flashpoint",) + _VAPOR_PRESSURE_KEYS,
        builder={**_VLE_SWEEP_BUILDER, "roles": ("solvent",)},
        comments=(
            "frac1 and frac2 define two endpoint solutions mixed along the composition line.",
        ),
    ),
    "LLE": _property_metadata(
        description="Liquid-liquid equilibrium for a ternary mixture.",
        system_scope="mixture",
        input_keys={
            "top_level": ("temperature", "massfraction"),
            "compound": ("frac1",),
        },
        required_keys=("temperature", "frac1"),
        builder={
            "roles": ("compound",),
        },
    ),
    "STABILITY": _property_metadata(
        description="Michelsen tangent-plane-distance stability test for a feed composition.",
        system_scope="mixture",
        input_keys={
            "top_level": ("temperature", "massfraction"),
            "compound": ("frac1",),
        },
        required_keys=("temperature", "frac1"),
        builder={
            "roles": ("compound",),
        },
    ),
    "SIGMAPROFILE": _property_metadata(
        description="Sigma profile for a solvent mixture.",
        system_scope="mixture",
        input_keys={
            "top_level": ("temperature", "massfraction"),
            "property": ("nprofile", "sigmamax"),
            "compound": ("frac1",),
        },
        required_keys=("temperature", "frac1"),
        builder={
            "roles": ("compound",),
        },
    ),
    "PURESIGMAPROFILE": _property_metadata(
        description="Sigma profile for pure compounds.",
        system_scope="pure_compounds",
        input_keys={
            "property": ("nprofile", "sigmamax"),
        },
        builder={
            "roles": ("compound",),
        },
        comments=(
            "Multiple COMPOUND blocks are treated as independent pure compounds.",
        ),
    ),
    "SIGMAPOTENTIAL": _property_metadata(
        description="Sigma potential for a solvent mixture.",
        system_scope="mixture",
        input_keys={
            "top_level": ("temperature", "massfraction"),
            "property": ("nprofile", "sigmamax"),
            "compound": ("frac1",),
        },
        required_keys=("temperature", "frac1"),
        builder={
            "roles": ("compound",),
        },
    ),
    "PURESIGMAPOTENTIAL": _property_metadata(
        description="Sigma potential for pure compounds.",
        system_scope="pure_compounds",
        input_keys={
            "property": ("nprofile", "sigmamax"),
        },
        builder={
            "roles": ("compound",),
        },
        comments=(
            "Multiple COMPOUND blocks are treated as independent pure compounds.",
        ),
    ),
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
    "COSMORS",
    "COSMO-RS",
    "COSMOSAC",
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

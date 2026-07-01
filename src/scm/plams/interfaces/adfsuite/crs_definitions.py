from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple


RESULT_TABLE_UNSUPPORTED_PROPERTIES = {
    "SIGMAPROFILE",
    "PURESIGMAPROFILE",
    "SIGMAPOTENTIAL",
    "PURESIGMAPOTENTIAL",
}
RESULT_TABLE_LLE_PROPERTIES = {"LLE", "STABILITY", "BINMIXCOEF", "TERNARYMIX"}
RESULT_TABLE_COMPONENT_BASE_COLUMNS = (
    "property",
    "mixture",
    "cid",
    "name",
)
RESULT_TABLE_LLE_BASE_COLUMNS = (
    "property",
    "mixture",
    "tie_line",
    "cid",
    "name",
)

RESULT_TABLE_COMPONENT_DEFAULT_QUANTITIES = (
    "frac1",
    "frac2",
    "solvent fraction",
    "composition molar fraction",
    "gamma",
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
    "phiI",
    "phiII",
    "tpd_w",
    "isobar",
    "flashpoint",
)
RESULT_TABLE_MIXTURE_EXTRA_QUANTITIES = (
    "Gibbs energy",
    "status_msg",
    "llle_detected",
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
RESULT_TABLE_QUANTITY_METADATA = {
    "frac1": {"symbol": "x", "name": "Feed molar composition", "unit": "fraction"},
    "x": {"symbol": "x", "name": "Feed molar composition", "unit": "fraction"},
    "frac2": {"symbol": "x_2", "name": "Second feed molar composition", "unit": "fraction"},
    "composition molar fraction": {"symbol": "x_calc", "name": "Calculated molar composition", "unit": "fraction"},
    "solvent fraction": {"symbol": "x_solvent", "name": "Solvent fraction", "unit": "fraction"},
    "poly fraction": {"symbol": "x_poly", "name": "Polymer fraction", "unit": "fraction"},
    "gamma": {"symbol": "gamma", "name": "Activity coefficient", "unit": "dimensionless"},
    "gammaI": {"symbol": "gamma_I", "name": "Activity coefficient in phase I", "unit": "dimensionless"},
    "gammaII": {"symbol": "gamma_II", "name": "Activity coefficient in phase II", "unit": "dimensionless"},
    "gamma_wf": {"symbol": "gamma_wf", "name": "Weight-fraction activity coefficient", "unit": "dimensionless"},
    "gamma_vf": {"symbol": "gamma_vf", "name": "Volume-fraction activity coefficient", "unit": "dimensionless"},
    "logp": {"symbol": "logP", "name": "Log10 partition coefficient", "unit": "dimensionless"},
    "fh_chi": {"symbol": "chi_FH", "name": "Flory-Huggins chi", "unit": "dimensionless"},
    "mu": {
        "symbol": "mu",
        "name": "Pseudo-chemical potential",
        "unit": "kcal/mol",
        "note": "Not a true chemical potential and does not include the ideal mixing term.",
    },
    "mu pure": {"symbol": "mu_pure", "name": "Pure liquid-phase pseudo-chemical potential", "unit": "kcal/mol"},
    "mu gas": {
        "symbol": "mu_gas",
        "name": "Gas-phase pseudo-chemical potential",
        "unit": "kcal/mol",
        "note": "Gas-phase reference pseudo-chemical potential, optionally refined using input vapor pressure data.",
    },
    "mu in solvent 1": {
        "symbol": "mu_solv_1",
        "name": "Pseudo-chemical potential in solvent 1",
        "unit": "kcal/mol",
    },
    "mu in solvent 2": {
        "symbol": "mu_solv_2",
        "name": "Pseudo-chemical potential in solvent 2",
        "unit": "kcal/mol",
    },
    "E gas": {"symbol": "E_gas", "name": "Gas-phase energy", "unit": "kcal/mol"},
    "G solute": {
        "symbol": "G_solute",
        "name": "Solute pseudo-energy",
        "unit": "kcal/mol",
        "note": "Computed as gas-phase energy plus solvation free energy.",
    },
    "excess G": {"symbol": "G_excess", "name": "Excess Gibbs energy", "unit": "kcal/mol"},
    "excess H": {"symbol": "H_excess", "name": "Excess enthalpy", "unit": "kcal/mol"},
    "Gibbs energy": {
        "symbol": "G_CRS",
        "name": "Pseudo-Gibbs energy",
        "unit": "kcal/mol",
        "note": "Computed as a composition-weighted sum of pseudo-chemical potentials plus the ideal mixing term; not a true thermodynamic Gibbs energy.",
    },
    "Gibbs energy of mixing": {"symbol": "DeltaG_mix", "name": "Gibbs energy of mixing", "unit": "kcal/mol"},
    "Enthalpy of vaporization": {"symbol": "DeltaH_vap", "name": "Enthalpy of vaporization", "unit": "kcal/mol"},
    "deltag": {
        "symbol": "DeltaG_solv",
        "name": "Solvation free energy",
        "unit": "kcal/mol",
        "note": "Computed as mu - mu_gas plus a gas-to-solution standard-state Gibbs energy correction; affected by input density or solvent density.",
    },
    "vapor pressure": {"symbol": "p_vap", "name": "Vapor pressure", "unit": "bar"},
    "henryc": {
        "symbol": "H",
        "name": "Henry's law constant",
        "unit": "mol/(L atm)",
        "note": "Henry's law constant computed as exp((mu_gas - mu)/RT) divided by the solvent molar volume.",
    },
    "henrycnodim": {"symbol": "H_cc", "name": "Dimensionless Henry's law constant", "unit": "dimensionless"},
    "temperature": {"symbol": "T", "name": "Temperature", "unit": "K"},
    "pressure": {"symbol": "P", "name": "Pressure", "unit": "bar"},
    "isobar": {"symbol": "isobar", "name": "Isobaric", "unit": ""},
    "flashpoint": {"symbol": "flashpoint", "name": "Flash-point calculation", "unit": ""},
    "showmiscgap": {
        "symbol": "misc_gap",
        "name": "Miscibility gap detected",
        "unit": "",
        "note": "Estimated by interpolation.",
    },
    "unstable": {
        "symbol": "unstable",
        "name": "TPD unstable",
        "unit": "",
        "note": "Indicates instability from the tangent-plane distance test.",
    },
    "converged": {"symbol": "converged", "name": "LLE converged", "unit": ""},
    "status_msg": {"symbol": "status", "name": "Status message", "unit": ""},
    "w_min": {"symbol": "w_min", "name": "Trial-phase composition minimizing TPD(w)", "unit": "fraction"},
    "tpd_w": {"symbol": "TPD_min", "name": "Minimum tangent-plane distance", "unit": ""},
    "phiI": {"symbol": "phi_I", "name": "Phase I fraction", "unit": "fraction"},
    "phiII": {"symbol": "phi_II", "name": "Phase II fraction", "unit": "fraction"},
    "actI": {"symbol": "a_I", "name": "Activity in phase I", "unit": "dimensionless"},
    "actII": {"symbol": "a_II", "name": "Activity in phase II", "unit": "dimensionless"},
    "act_interp": {
        "symbol": "a*",
        "name": "Interpolated activity",
        "unit": "dimensionless",
        "note": "Estimated from interpolated miscibility-gap data; not from a full LLE calculation.",
    },
    "xI": {"symbol": "x_I", "name": "Molar composition in phase I", "unit": "fraction"},
    "xII": {"symbol": "x_II", "name": "Molar composition in phase II", "unit": "fraction"},
    "solution molar fraction": {
        "symbol": "x_sol",
        "name": "Solubility molar fraction",
        "unit": "fraction",
        "note": "Equilibrium molar-fraction solubility; density is only used for volume-based solubilities.",
    },
    "solubility mol_per_L_solvent": {
        "symbol": "S_mol_L_solvent",
        "name": "Solubility in moles per L solvent",
        "unit": "mol/L solvent",
    },
    "solubility g_per_L_solvent": {
        "symbol": "S_g_L_solvent",
        "name": "Solubility in grams per L solvent",
        "unit": "g/L solvent",
    },
    "solubility mol_per_L_solution": {
        "symbol": "S_mol_L_solution",
        "name": "Solubility in moles per L solution",
        "unit": "mol/L solution",
    },
    "solubility g_per_L_solution": {
        "symbol": "S_g_L_solution",
        "name": "Solubility in grams per L solution",
        "unit": "g/L solution",
    },
    "solubility massfrac": {"symbol": "w_solubility", "name": "Solubility mass fraction", "unit": "fraction"},
}


def _property_metadata(
    *,
    description: str,
    system_scope: str,
    input_keys: Dict[str, Sequence[str]],
    required_keys: Sequence[str] = (),
    notes: Sequence[str] = (),
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
        "notes": tuple(notes),
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
        notes=(
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
        notes=(
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
        notes=(
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
        notes=(
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
        notes=(
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
_crs_json_path = Path(os.environ["AMSBIN"]) / "../data/input_def/crs.json"

CRS_DATA = _load_crs_input_block_metadata(_crs_json_path.resolve())

_CRS_METHOD_PARAMETER_BLOCKS = frozenset(
    ("CRSParameters", "SACParameters", "Dispersion", "Epsilon")
)

def _load_crs_method_parameters_metadata(
    json_path: Path,
    exposed_methods: Sequence[str],
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    with json_path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    if data.get("schema_version") != 1:
        raise ValueError(f"Unsupported CRS method parameter schema_version in {json_path}")

    presets = data.get("presets")
    if not isinstance(presets, list):
        raise ValueError(f"CRS method parameter file {json_path} must contain a presets list")

    result: Dict[str, Dict[str, Dict[str, Any]]] = {}
    exposed_method_set = set(exposed_methods)

    for method in exposed_methods:
        candidates = [
            preset for preset in presets
            if preset.get("input_method") == method
        ]

        if not candidates:
            continue

        default_candidates = [
            preset for preset in candidates
            if preset.get("default", False)
        ]

        if len(default_candidates) == 1:
            selected = default_candidates[0]
        elif len(default_candidates) > 1:
            raise ValueError(f"Multiple default CRS method parameter presets for {method!r}")
        elif len(candidates) == 1:
            selected = candidates[0]
        else:
            raise ValueError(
                f"Multiple CRS method parameter presets for {method!r}, but none is marked default"
            )

        result[method] = _extract_crs_method_parameter_blocks(selected)

    missing_methods = exposed_method_set - set(result)
    if missing_methods:
        missing = ", ".join(sorted(missing_methods))
        raise ValueError(f"Missing CRS method parameter presets for exposed methods: {missing}")

    return result


def _extract_crs_method_parameter_blocks(
    preset: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    blocks: Dict[str, Dict[str, Any]] = {}

    for block_name in _CRS_METHOD_PARAMETER_BLOCKS:
        block = preset.get(block_name)
        if block is None:
            continue
        if not isinstance(block, dict):
            preset_id = preset.get("id", "<unknown>")
            raise ValueError(f"{block_name} in CRS method preset {preset_id!r} must be an object")
        blocks[block_name] = block.copy()

    return blocks

CRS_METHODS = (
    "COSMORS",
    "COSMO-RS",
    "COSMOSAC",
    "COSMOSAC2013",
    "COSMOSAC2016",
    "COSMOSACDHB",
    "COSMOSACDHB-MESP",
)

_CRS_METHOD_PARAMETER_METHODS = (
    "COSMO-RS",
    "COSMOSAC2013",
    "COSMOSAC2016",
    "COSMOSACDHB",
    "COSMOSACDHB-MESP",
)

_crs_method_parameters_path = Path(os.environ["AMSBIN"]) / "../data/crs/method_parameters.json"
CRS_METHOD_PARAMETERS_METADATA = _load_crs_method_parameters_metadata(
    _crs_method_parameters_path.resolve(),
    _CRS_METHOD_PARAMETER_METHODS,
)

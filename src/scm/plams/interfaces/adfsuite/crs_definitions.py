from __future__ import annotations

import json
import os
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Tuple


__all__ = [
    "CRS_METHODS",
    "CRSInputMetadata",
    "CRSMethodParameterMetadata",
    "get_block",
    "get_block_child_names",
    "get_block_keys",
    "get_crs_input_data",
    "get_crs_input_def_path",
    "get_crs_method_parameters_metadata",
    "get_crs_method_parameters_path",
]


CRSInputMetadata = Dict[str, Dict[str, Any]]
CRSMethodParameterMetadata = Dict[str, Dict[str, Any]]

CRS_METHODS = (
    "COSMO-RS",
    "COSMOSAC2013",
    "COSMOSAC2016",
    "COSMOSACDHB",
    "COSMOSACDHB-MESP",
)

_CRS_METHOD_PARAMETER_BLOCKS = frozenset(
    ("CRSParameters", "SACParameters", "Dispersion", "Epsilon")
)


def get_crs_input_def_path() -> Path:
    """Return the path to ``data/input_def/crs.json``.

    The lookup is lazy so importing this module does not require a configured AMS
    environment.
    """
    return _find_file(
        relative_parts=("data", "input_def", "crs.json"),
        description="data/input_def/crs.json",
    )


def get_crs_method_parameters_path() -> Path:
    """Return the path to ``data/crs/method_parameters.json``."""
    return _find_file(
        relative_parts=("data", "crs", "method_parameters.json"),
        description="data/crs/method_parameters.json",
    )


@lru_cache(maxsize=None)
def get_crs_input_data() -> CRSInputMetadata:
    """Return normalized metadata from ``data/input_def/crs.json``."""
    with get_crs_input_def_path().open(encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError(f"CRS input definition file {get_crs_input_def_path()} must contain a JSON object")

    return _extract_crs_input_block_metadata(data)


@lru_cache(maxsize=None)
def get_crs_method_parameters_metadata() -> CRSMethodParameterMetadata:
    """Return normalized metadata from ``data/crs/method_parameters.json``."""
    with get_crs_method_parameters_path().open(encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError(
            f"CRS method parameter file {get_crs_method_parameters_path()} must contain a JSON object"
        )

    return _extract_crs_method_parameters_metadata(data, get_crs_method_parameters_path())


def get_block(data: Mapping[str, Any], *path: str) -> Dict[str, Any]:
    """Return a block metadata object by case-insensitive nested path.

    Examples:
        get_block(data, "compound") returns the COMPOUND block.
        get_block(data, "compound", "required") returns the COMPOUND/REQUIRED block.
    """
    if not path:
        raise ValueError("At least one block name must be provided")

    block = data[_normalize_name(path[0])]
    for name in path[1:]:
        block = block["blocks"][_normalize_name(name)]
    return block


def get_block_keys(data: Mapping[str, Any], *path: str) -> Dict[str, Any]:
    """Return scalar key metadata directly defined in the selected block."""
    return dict(get_block(data, *path)["keys"])


def get_block_child_names(data: Mapping[str, Any], *path: str) -> Tuple[str, ...]:
    """Return key and child-block names directly accepted by the selected block."""
    block = get_block(data, *path)
    return tuple(dict.fromkeys((*block["keys"], *block["blocks"])))


@lru_cache(maxsize=None)
def _find_file(relative_parts: Tuple[str, ...], description: str) -> Path:
    for root in _candidate_ams_roots():
        candidate = root.joinpath(*relative_parts)
        if candidate.is_file():
            return candidate.resolve()

    roots = ", ".join(str(root) for root in _candidate_ams_roots())
    raise FileNotFoundError(
        f"Could not locate {description}. Set AMSHOME or AMSBIN, or run from an AMS source tree. "
        f"Searched roots: {roots}"
    )


def _candidate_ams_roots() -> Tuple[Path, ...]:
    roots = []

    amshome = os.environ.get("AMSHOME")
    if amshome:
        roots.append(Path(amshome))

    amsbin = os.environ.get("AMSBIN")
    if amsbin:
        roots.append(Path(amsbin).resolve().parent)

    source_tree_root = _source_tree_root()
    if source_tree_root is not None:
        roots.append(source_tree_root)

    return tuple(dict.fromkeys(root.resolve() for root in roots))


def _source_tree_root() -> Path:
    # src/scm/plams/interfaces/adfsuite/crs_definitions.py -> AMS source root
    return Path(__file__).resolve().parents[8]


def _extract_crs_input_block_metadata(data: Mapping[str, Any]) -> CRSInputMetadata:
    result: CRSInputMetadata = {}

    for entry_name, entry_def in data.items():
        if not isinstance(entry_def, dict):
            continue
        if _skip_top_level_input_block(entry_name):
            continue

        category = entry_def.get("_category")
        if category == "key":
            result[_normalize_name(entry_name)] = _key_metadata(entry_name, entry_def)
        elif category == "block":
            result[_normalize_name(entry_name)] = _block_metadata(entry_name, entry_def)

    return result


def _skip_top_level_input_block(name: str) -> bool:
    return _normalize_name(name) in {
        "dispersion",
        "coef_cations",
        "coef_anions",
        "epsilon",
        "crsparameters",
        "sacparameters",
    }


def _key_metadata(name: str, key_def: Mapping[str, Any]) -> Dict[str, Any]:
    metadata = {
        "kind": "key",
        "input_name": name,
        "type": key_def.get("_type"),
        "description": key_def.get("_comment", ""),
        "unit": key_def.get("_unit"),
        "choices": tuple(key_def.get("_choices", ())),
        "hint": key_def.get("_hint", ""),
    }

    if "_default" in key_def:
        metadata["default"] = key_def["_default"]

    return metadata


def _block_metadata(name: str, block_def: Mapping[str, Any]) -> Dict[str, Any]:
    keys: Dict[str, Any] = {}
    blocks: Dict[str, Any] = {}

    for child_name, child_def in block_def.items():
        if not isinstance(child_def, dict):
            continue

        child_category = child_def.get("_category")
        if child_category == "key":
            keys[_normalize_name(child_name)] = _key_metadata(child_name, child_def)
        elif child_category == "block":
            blocks[_normalize_name(child_name)] = _block_metadata(child_name, child_def)

    return {
        "kind": "block",
        "input_name": name,
        "description": block_def.get("_comment", ""),
        "type": block_def.get("_type"),
        "header": bool(block_def.get("_header", False)),
        "unique": bool(block_def.get("_unique", False)),
        "keys": keys,
        "blocks": blocks,
    }


def _extract_crs_method_parameters_metadata(
    data: Mapping[str, Any], json_path: Path
) -> CRSMethodParameterMetadata:
    if data.get("schema_version") != 1:
        raise ValueError(f"Unsupported CRS method parameter schema_version in {json_path}")

    presets = data.get("presets")
    if not isinstance(presets, list):
        raise ValueError(f"CRS method parameter file {json_path} must contain a presets list")

    result: CRSMethodParameterMetadata = {}

    for preset in presets:
        if not isinstance(preset, dict):
            raise ValueError(f"CRS method parameter preset in {json_path} must be an object")

        parameter_set = preset.get("parameter_set")
        if not isinstance(parameter_set, str) or not parameter_set.strip():
            raise ValueError(f"CRS method parameter preset in {json_path} must define parameter_set")

        normalized_parameter_set = _normalize_name(parameter_set)
        if normalized_parameter_set in result:
            raise ValueError(f"Duplicate CRS parameter_set={normalized_parameter_set!r} in {json_path}")

        method = preset.get("method")
        if not isinstance(method, str) or not method.strip():
            raise ValueError(f"CRS parameter_set={normalized_parameter_set!r} must define method")

        parameter_blocks = preset.get("parameter_blocks")
        if not isinstance(parameter_blocks, dict):
            raise ValueError(f"CRS parameter_set={normalized_parameter_set!r} must define parameter_blocks")

        _validate_parameter_blocks(parameter_blocks, normalized_parameter_set, json_path)

        normalized_preset = dict(preset)
        normalized_preset["parameter_set"] = normalized_parameter_set
        normalized_preset["method"] = method.strip()
        result[normalized_parameter_set] = normalized_preset

    return result


def _validate_parameter_blocks(
    parameter_blocks: Mapping[str, Any], parameter_set: str, json_path: Path
) -> None:
    for block_name, block in parameter_blocks.items():
        if block_name not in _CRS_METHOD_PARAMETER_BLOCKS:
            allowed = ", ".join(sorted(_CRS_METHOD_PARAMETER_BLOCKS))
            raise ValueError(
                f"Unsupported parameter block {block_name!r} in CRS parameter_set={parameter_set!r} "
                f"from {json_path}. Allowed blocks: {allowed}"
            )
        if not isinstance(block, dict):
            raise ValueError(f"{block_name} in CRS parameter_set={parameter_set!r} from {json_path} must be an object")


def _normalize_name(name: str) -> str:
    return name.strip().lower()

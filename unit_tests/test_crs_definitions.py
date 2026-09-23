from pathlib import Path
from typing import Any, Dict

import pytest

from scm.plams.interfaces.adfsuite.crs_definitions import (
    _extract_crs_input_block_metadata,
    _extract_crs_method_parameters_metadata,
    get_block,
    get_block_keys,
)


def test_input_metadata_and_nested_lookup() -> None:
    data = _extract_crs_input_block_metadata(
        {
            "COMPOUND": {
                "_category": "block",
                "Polymer": {"_category": "key", "_type": "bool", "_default": False},
                "FORM": {
                    "_category": "block",
                    "_header": True,
                    "Name": {"_category": "key", "_type": "string"},
                },
            }
        }
    )

    keys = get_block_keys(data, "Compound")
    assert set(keys) == {"polymer"}
    assert keys["polymer"]["input_name"] == "Polymer"
    assert keys["polymer"]["type"] == "bool"
    assert keys["polymer"]["default"] is False
    form = get_block(data, "compound", "Form")
    assert form["input_name"] == "FORM"
    assert form["header"] is True
    assert set(get_block_keys(data, "COMPOUND", "form")) == {"name"}
    assert "default" not in form["keys"]["name"]


def test_method_parameter_metadata() -> None:
    parameter_blocks = {"CRSParameters": {"chb": 1.5}, "Dispersion": {"enabled": False}}
    data = {
        "schema_version": 1,
        "presets": [
            {"parameter_set": " Example ", "method": " COSMO-RS ", "parameter_blocks": parameter_blocks}
        ],
    }

    result = _extract_crs_method_parameters_metadata(data, Path("method_parameters.json"))

    assert set(result) == {"example"}
    assert result["example"]["parameter_set"] == "example"
    assert result["example"]["method"] == "COSMO-RS"
    assert result["example"]["parameter_blocks"] == parameter_blocks


@pytest.mark.parametrize(
    "case,match",
    [("unsupported_schema", "schema_version"), ("duplicate_parameter_set", "Duplicate.*example")],
)
def test_invalid_method_parameter_metadata(case: str, match: str) -> None:
    data: Dict[str, Any] = {"schema_version": 1, "presets": []}
    if case == "unsupported_schema":
        data["schema_version"] = 2
    else:
        data["presets"] = [
            {"parameter_set": name, "method": "COSMO-RS", "parameter_blocks": {}}
            for name in ("Example", " example ")
        ]

    with pytest.raises(ValueError, match=match):
        _extract_crs_method_parameters_metadata(data, Path("method_parameters.json"))

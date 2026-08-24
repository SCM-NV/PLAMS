import importlib


def test_crs_definitions_import_without_ams_environment(monkeypatch):
    monkeypatch.delenv("AMSBIN", raising=False)
    monkeypatch.delenv("AMSHOME", raising=False)

    from scm.plams.interfaces.adfsuite import crs_definitions

    importlib.reload(crs_definitions)

    assert crs_definitions.CRS_METHODS


def test_crs_input_definition_metadata_helpers_support_nested_blocks():
    from scm.plams.interfaces.adfsuite import crs_definitions

    data = crs_definitions.get_crs_input_data()

    compound_keys = crs_definitions.get_block_keys(data, "compound")
    required_keys = crs_definitions.get_block_keys(data, "compound", "required")

    assert compound_keys["frac1"]["type"] == "float"
    assert set(required_keys) == {"name", "count"}


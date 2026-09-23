from typing import Any, Type

import pytest

from scm.plams import CRSJob, Settings

import json
from pathlib import Path

from scm.plams.interfaces.adfsuite import crs_input_builder
from scm.plams.interfaces.adfsuite.crs_definitions import (
    _extract_crs_input_block_metadata,
)

@pytest.fixture(autouse=True)
def minimal_crs_metadata(monkeypatch: pytest.MonkeyPatch) -> None:
    path = Path(__file__).parent / "fixtures" / "crs_minimal.json"
    with path.open(encoding="utf-8") as handle:
        metadata = _extract_crs_input_block_metadata(json.load(handle))

    monkeypatch.setattr(
        crs_input_builder,
        "get_crs_input_data",
        lambda: metadata,
    )

def test_to_settings() -> None:
    builder = CRSJob.input_builder("PURESIGMAPROFILE", nprofile=60, sigmamax=0.03)
    builder.add_compound("water.coskf")

    settings = builder.to_settings()

    assert settings.input.method == "COSMO-RS"
    assert settings.input.property._h == "PURESIGMAPROFILE"
    assert settings.input.property.nprofile == 60
    assert settings.input.property.sigmamax == 0.03
    assert len(settings.input.compound) == 1
    assert settings.input.compound[0]._h == "water.coskf"


@pytest.mark.parametrize(
    "key,value,error",
    [("unknown_key", 1, AttributeError), ("nprofile", "50", TypeError), ("nprofile", True, TypeError)],
)
def test_invalid_input(key: str, value: Any, error: Type[Exception]) -> None:
    with pytest.raises(error, match=key):
        CRSJob.input_builder("PURESIGMAPROFILE", **{key: value})


def test_required_input() -> None:
    builder = CRSJob.input_builder("SOLUBILITY")
    with pytest.raises(ValueError, match="temperature"):
        builder.to_settings(include_compounds=False)


@pytest.mark.parametrize("stage,missing", [(0, "solvent"), (1, "solute"), (2, "meltingpoint.*hfusion")])
def test_solubility_requirements(stage: int, missing: str) -> None:
    builder = CRSJob.input_builder("SOLUBILITY", mode="solid", temperature=298.15)
    if stage >= 1:
        builder.add_solvent("water.coskf", frac1=1.0)
    if stage >= 2:
        builder.add_solute("benzene.coskf")
    with pytest.raises(ValueError, match=missing):
        builder.to_settings()
    if stage == 2:
        builder.mode = "gas"
        assert builder.to_settings().input.property.isobar is True


@pytest.mark.parametrize("output_method", ["to_settings", "to_inputs"])
@pytest.mark.parametrize("has_compound", [False, True])
def test_exclude_compounds(output_method: str, has_compound: bool) -> None:
    if output_method == "to_inputs":
        pytest.importorskip("scm.inputs")
    builder = CRSJob.input_builder("SOLUBILITY", temperature=298.15)
    if has_compound:
        builder.add_solvent("water.coskf")  # Missing frac1 must also be ignored.
    output = getattr(builder, output_method)(include_compounds=False)
    if output_method == "to_settings":
        assert "compound" not in output.input
    else:
        assert output.COMPOUND == []


def test_to_inputs_conversion() -> None:
    pytest.importorskip("scm.inputs")
    compound = Settings()
    form = [Settings(), Settings()]
    form[0]._h = "conformer0.coskf"
    form[1]._h = "conformer1.coskf"
    compound.form = form
    builder = CRSJob.input_builder("PUREVAPORPRESSURE", temperature=[280.0, 300.0])
    builder.add_compound(compound)

    crs = builder.to_inputs()

    assert crs.TEMPERATURE == [280.0, 300.0]
    assert len(crs.COMPOUND) == 1
    assert [item.header for item in crs.COMPOUND[0].FORM] == ["conformer0.coskf", "conformer1.coskf"]

def test_to_job() -> None:
    builder = CRSJob.input_builder("PURESIGMAPROFILE", nprofile=60)
    builder.add_compound("water.coskf")
    job = builder.to_job(name="sigma_profile")

    assert isinstance(job, CRSJob)
    assert job.name == "sigma_profile"
    text = job.get_input().lower()
    assert "property puresigmaprofile" in text
    assert "nprofile 60" in text
    assert "compound water.coskf" in text

def test_modified_typed_inputs_in_job() -> None:
    pytest.importorskip("scm.inputs")

    builder = CRSJob.input_builder("LLE", temperature=298.15)
    builder.add_compound("Water.coskf", frac1=0.33)
    builder.add_compound("Ethanol.coskf", frac1=0.33)
    builder.add_compound("Benzene.coskf", frac1=0.34)

    crs = builder.to_inputs()
    crs.TECHNICAL.LLE.eps_g = 1.0e-5
    crs.TECHNICAL.LLE.debug = True

    job = CRSJob(settings=crs)
    text = job.get_input().lower()
    lines = [line.split() for line in text.splitlines()]

    assert ["technical"] in lines
    assert ["lle"] in lines
    assert ["debug"] in lines
    eps_g = next(fields[1] for fields in lines if fields and fields[0] == "eps_g")
    assert float(eps_g) == pytest.approx(1.0e-5)
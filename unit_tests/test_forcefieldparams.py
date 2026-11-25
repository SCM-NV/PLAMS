from scm.plams.interfaces.adfsuite.forcefieldparams import ForceFieldPatch


def test_fragment_serialization_preserves_parameters() -> None:
    parent = ForceFieldPatch.from_frcmod(
        "Magnesium ion parameters\n"
        "MASS\n"
        "MG 24.305 0.0\n\n"
        "BOND\n\n"
        "ANGLE\n\n"
        "DIHE\n\n"
        "IMPROPER\n\n"
        "NONBON\n"
        "  MG 1.360 0.01020237\n"
    )

    fragment = parent.get_fragment(["MG"])
    reparsed = ForceFieldPatch(str(fragment))

    assert reparsed.types == reparsed.ljtypes == ["MG"]
    assert reparsed.typelines == parent.typelines
    assert reparsed.ljlines == parent.ljlines

import re
from scm.plams.version import __version__


def test_version_pattern_correct():
    assert re.fullmatch(r"20\d{2}\.\d{3}", __version__)

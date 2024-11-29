import pytest

from scm.plams.core.enums import JobStatus


class TestStrEnum:
    """
    Check compatibility between Python 3.11 StrEnum and backwards compatible custom version
    """

    @pytest.mark.parametrize("value,expected", [("foo", None), ("created", JobStatus.CREATED), ("cREated", None)])
    def test_ctor(self, value, expected):
        if expected is None:
            with pytest.raises(ValueError):
                JobStatus(value)
        else:
            assert JobStatus(value) == expected

    def test_str_operations(self):
        assert JobStatus.CREATED == "created"
        assert JobStatus.CREATED != "CREATED"
        assert f"status={JobStatus.CREATED}" == "status=created"

import pytest
from unittest.mock import patch, Mock

from scm.plams.interfaces.adfsuite.errors import AMSBINEnvVarNotSetError, AMSExecutionError, AMSVersionError
from scm.plams.interfaces.adfsuite.utils import requires_ams


class TestRequiresAms:

    @requires_ams()
    def this_function_requires_any_ams(self):
        return True

    @requires_ams("2024.2")
    def this_function_requires_2025_ams(self):
        return True

    def test_requires_ams_when_no_amsbin_errors(self, monkeypatch):
        monkeypatch.setenv("AMSBIN", "")
        with pytest.raises(AMSBINEnvVarNotSetError):
            self.this_function_requires_any_ams()

    def test_requires_ams_when_no_exe_errors(self):
        mock_result = Mock()
        mock_result.stdout = ""
        mock_result.stderr = "Something went wrong"

        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(AMSExecutionError):
                self.this_function_requires_any_ams()

    def test_requires_ams_when_incompatible_version_errors(self):
        mock_result = Mock()
        mock_result.stdout = "release=2024.102"
        mock_result.stderr = ""

        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(AMSVersionError):
                self.this_function_requires_2025_ams()

    def test_requires_ams_when_compatible_version_returns(self):
        mock_result = Mock()
        mock_result.stdout = "release=2024.201"
        mock_result.stderr = ""

        with patch("subprocess.run", return_value=mock_result):
            assert self.this_function_requires_any_ams()
            assert self.this_function_requires_2025_ams()

        mock_result.stdout = "release=2025.101"
        with patch("subprocess.run", return_value=mock_result):
            assert self.this_function_requires_any_ams()
            assert self.this_function_requires_2025_ams()

import os
import pytest
import subprocess

from scm.plams.core.private import run_with_timeout


class TestRunWithTimeout:

    def test_happy_command(self):
        run_with_timeout(["echo", "1"])

    def test_non_existent_command(self):
        with pytest.raises(FileNotFoundError):
            run_with_timeout(["cmddoesnotexist"])

    def test_stderror_command(self):
        with pytest.raises(subprocess.CalledProcessError):
            run_with_timeout(["bash", "-c", "echo 'Something went wrong' >&2; exit 1"])

    def test_hanging_command(self):
        cmd = ["ping", "127.0.0.1", "-n", "2"] if os.name == "nt" else ["sleep", "1"]
        with pytest.raises(TimeoutError):
            run_with_timeout(cmd, timeout=0.5, poll_interval=0.1)
        run_with_timeout(cmd, timeout=2)

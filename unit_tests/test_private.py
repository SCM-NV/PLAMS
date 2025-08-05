import os

from scm.plams.core.private import run_with_timeout


class TestSafeSystemCall:

    def test_happy_command(self):
        assert run_with_timeout(["echo", "1"])

    def test_non_existent_command(self):
        assert not run_with_timeout(["cmddoesnotexist"])

    def test_hanging_command(self):
        cmd = ["ping", "127.0.0.1", "-n", "2"] if os.name == "nt" else ["sleep", "1"]
        assert not run_with_timeout(cmd, timeout=0.5, poll_interval=0.1)
        assert run_with_timeout(cmd, timeout=2)

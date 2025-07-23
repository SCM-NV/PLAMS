import os

from scm.plams.core.private import safe_system_call


class TestSafeSystemCall:

    def test_happy_command(self):
        assert safe_system_call("echo 1")

    def test_non_existent_command(self):
        assert not safe_system_call("cmddoesnotexist")

    def test_hanging_command(self):
        cmd = "timeout /T 1" if os.name == "nt" else "sleep 1"
        assert not safe_system_call(cmd, timeout=0.5, poll_interval=0.1)
        assert safe_system_call(cmd, timeout=2)

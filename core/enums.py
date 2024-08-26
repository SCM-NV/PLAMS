# Import StrEnum for Python >=3.11, otherwise use backwards compatible class
try:
    from enum import StrEnum
except ImportError:
    from enum import Enum

    class StrEnum(str, Enum):
        """
        Enum where all members are strings and so can be used for string comparison.
        """

        def __str__(self):
            return self.value


__all__ = ["JobStatus"]


class JobStatus(StrEnum):
    """
    Status of PLAMS job.
    """

    CREATED = "created"
    STARTED = "started"
    REGISTERED = "registered"
    RUNNING = "running"
    FINISHED = "finished"
    CRASHED = "crashed"
    FAILED = "failed"
    SUCCESSFUL = "successful"
    COPIED = "copied"
    PREVIEW = "preview"
    DELETED = "deleted"

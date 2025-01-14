from typing import Literal, Union

# Import StrEnum for Python >=3.11, otherwise use backwards compatible class
try:
    from enum import StrEnum  # type: ignore
except ImportError:
    from enum import Enum

    class StrEnum(str, Enum):  # type: ignore
        """
        Enum where all members are strings and so can be used for string comparison.
        """

        def __str__(self) -> str:
            return str(self.value)


__all__ = ["JobStatus", "JobStatusType"]


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


JobStatusType = Union[
    JobStatus,
    Literal[
        "created",
        "started",
        "registered",
        "running",
        "finished",
        "crashed",
        "failed",
        "successful",
        "copied",
        "preview",
        "deleted",
    ],
]

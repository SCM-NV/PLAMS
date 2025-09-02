from typing import Union, Optional


class AMSError(Exception):
    """
    General error relating to the AMS executable.
    """


class AMSBINEnvVarNotSetError(AMSError):
    """
    Error due to $AMSBIN environment variable not being set.
    """

    def __init__(self) -> None:
        super().__init__("$AMSBIN environment variable is not set")


class AMSExecutionError(AMSError):
    """
    Error due to the ams process not running successfully.
    """

    def __init__(self, command: str, error: Union[str, Exception]):
        super().__init__(f"Command: '{command}' did not run successfully.\nError was: '{error}'.")


class AMSVersionError(AMSError):
    """
    Error due to the version of AMS being incompatible for the required functionality.
    """

    def __init__(self, version: Optional[str], minimum_version: Optional[str]):
        super().__init__(
            f"Current version of AMS is '{version}', but this operation requires AMS '>={minimum_version}'."
        )

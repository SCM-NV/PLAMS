__all__ = [
    "PlamsError",
    "FileError",
    "ResultsError",
    "JobError",
    "PTError",
    "UnitsError",
    "MoleculeError",
    "TrajectoryError",
]


class PlamsError(Exception):
    """General PLAMS error."""


class FileError(PlamsError):
    """File or filesystem related error."""


class ResultsError(PlamsError):
    """|Results| related error."""


class JobError(PlamsError):
    """|Job| related error."""


class PTError(PlamsError):
    """:class:`Periodic table<scm.plams.utils.PeriodicTable>` error."""


class UnitsError(PlamsError):
    """:class:`Units converter<scm.plams.utils.Units>` error."""


class MoleculeError(PlamsError):
    """|Molecule| related error."""


class TrajectoryError(PlamsError):
    """:class:`Trajectory<scm.plams.trajectories.TrajectoryFile>` error."""

class ExampleError(PlamsError):
    pass

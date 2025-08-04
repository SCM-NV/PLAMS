__all__ = [
    "PlamsError",
    "FileError",
    "ResultsError",
    "JobError",
    "PTError",
    "UnitsError",
    "MoleculeError",
    "TrajectoryError",
    "MissingOptionalPackageError",
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


class MissingOptionalPackageError(PlamsError):
    """Missing optional package related error."""

    extras_install = {
        "rdkit": "chem",
        "ase": "chem",
        "psutil": "ams",
        "ubjson": "ams",
        "watchdog": "ams",
        "scipy": "analysis",
        "matplotlib": "analysis",
        "pandas": "analysis",
        "networkx": "analysis",
        "natsort": "analysis",
        "h5py": "analysis",
        "ipython": "analysis",
    }

    ams_install = {"scm.amspipe": "$AMSHOME/scripting/scm/amspipe"}

    def __init__(self, package_name: str):
        msg = f"The optional package '{package_name}' is required for this PLAMS functionality, but is not available. "
        if (extras_name := self.extras_install.get(package_name, None)) is not None:
            msg += f"It can be installed using the command: pip install 'plams[{extras_name}]'. "
        elif (ams_path := self.ams_install.get(package_name, None)) is not None:
            msg += f"It can be installed using the command: pip install {ams_path}. "
        msg += "Please install and try again."

        super().__init__(msg)

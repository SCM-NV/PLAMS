import os
import re

from typing import (
    Dict,
    List,
    Literal,
    Set,
    Tuple,
    Union,
    Optional,
    TYPE_CHECKING,
    Any,
    Sequence,
    Callable,
    TypeVar,
    Iterator,
    Type,
    cast,
)

import numpy as np

from scm.plams.core.basejob import SingleJob
from scm.plams.core.errors import FileError, JobError, MissingOptionalPackageError, PlamsError, PTError, ResultsError
from scm.plams.core.functions import get_config, log, parse_heredoc, requires_optional_package
from scm.plams.core.private import sha256
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.mol.atom import Atom
from scm.plams.mol.bond import Bond
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.converters import gaussian_output_to_ams, qe_output_to_ams, vasp_output_to_ams
from scm.plams.tools.kftools import KFFile, KFReader
from scm.plams.tools.units import Units

try:
    from scm.pisa.block import DriverBlock

    _has_scm_pisa = True
except ImportError:
    _has_scm_pisa = False

try:
    from scm.base import ChemicalSystem

    _has_scm_chemsys = True
except ImportError:
    _has_scm_chemsys = False

if TYPE_CHECKING:
    from scm.plams.core.jobrunner import JobRunner
    from scm.plams.core.jobmanager import JobManager
    from scm.plams.tools.kftools import TRead
    from scm.plams.interfaces.adfsuite.forcefieldparams import ForceFieldPatch
    from ase import Atoms as AseAtoms
    from watchdog.events import FileSystemEvent

T = TypeVar("T")


try:
    from watchdog.events import FileModifiedEvent, PatternMatchingEventHandler
    from watchdog.observers import Observer

    _has_watchdog = True

    class AMSJobLogTailHandler(PatternMatchingEventHandler):
        def __init__(self, job: "AMSJob", jobmanager: "JobManager"):
            super().__init__(
                patterns=[os.path.join(jobmanager.workdir, f"{job.name}*", "ams.log")],
                ignore_patterns=["*.rkf", "*.out"],
                ignore_directories=True,
                case_sensitive=True,
            )
            self._job = job
            self._seekto = 0

        def on_any_event(self, event: "FileSystemEvent") -> None:
            if (
                self._job.path is not None
                and event.src_path == os.path.join(self._job.path, "ams.log")
                and isinstance(event, FileModifiedEvent)
            ):
                try:
                    with open(event.src_path, "r") as f:
                        f.seek(self._seekto)
                        while True:
                            line = f.readline()
                            if not line:
                                break
                            log(f"{self._job.name}: " + line[25:-1])
                        self._seekto = f.tell()
                except FileNotFoundError:
                    self._seekto = 0

        def trigger(self) -> None:
            if self._job.path is None:
                return
            src_path = os.path.join(self._job.path, "ams.log")
            self.on_any_event(FileModifiedEvent(src_path))

except ImportError:
    _has_watchdog = False


__all__ = ["AMSJob", "AMSResults"]


class AMSResults(Results):
    """Results from an AMS calculation.

    Methods with a unit argument convert numerical results to that unit.

    Methods that expose RKF variables directly, including :meth:`readrkf`, :meth:`read_rkf_section`, and
    :meth:`get_history_property`, return values in their RKF storage units without conversion.

    When running with ``amspython``, use ``AKFReader`` to inspect the unit and metadata of an RKF variable::

        from scm.akfreader import AKFReader

        reader = AKFReader(results.rkfpath("engine"))
        print(reader.units("AMSResults%Energy"))
        print(reader.description("AMSResults%Energy"))
    """

    job: "AMSJob"

    def __init__(self, *args: Any, **kwargs: Any):
        Results.__init__(self, *args, **kwargs)
        self.rkfs: Dict[str, KFFile] = {}

    def collect(self) -> None:
        """Collect the job files and create an |KFFile| instance for each RKF file.

        Engine result files are discovered through ``ams.rkf`` and stored in :attr:`rkfs` under
        their filenames without the ``.rkf`` extension. This method is called automatically when
        job execution finishes.
        """
        Results.collect(self)
        self.collect_rkfs()

    def collect_rkfs(self) -> None:
        """Collect the RKF files referenced by the main ``ams.rkf`` file.

        Populate :attr:`rkfs` with the main file under ``ams`` and engine files under their filename
        without the ``.rkf`` extension. Log a warning when ``ams.rkf`` is unavailable.
        """
        rkfname = "ams.rkf"
        if rkfname in self.files:
            main = KFFile(str(self.job.get_path() / rkfname))
            n = main.read_int("EngineResults", "nEntries")
            for i in range(1, n + 1):
                files = main.read_string("EngineResults", "Files({})".format(i)).split("\x00")
                if files[0].endswith(".rkf"):
                    key = files[0][:-4]
                    self.rkfs[key] = KFFile(str(self.job.get_path() / files[0]))
            self.rkfs["ams"] = main

        else:
            log("WARNING: Main KF file {} not present in {}".format(rkfname, str(self.job.get_path())), 1)

    def _copy_to(self, newresults: "AMSResults") -> None:
        super()._copy_to(newresults)
        newresults.rkfs = {}
        newresults.collect_rkfs()

    def refresh(self) -> None:
        """Refresh the contents of ``files`` list.

        Also update the |KFFile| instances in :attr:`rkfs` if the job folder has moved since the
        results were pickled.
        """
        Results.refresh(self)
        to_remove = []
        for key, val in self.rkfs.items():
            if not os.path.isfile(val.path):
                if os.path.dirname(val.path) != self.job.path:
                    guessnewpath = str(self.job.get_path() / os.path.basename(val.path))
                    if os.path.isfile(guessnewpath):
                        self.rkfs[key] = KFFile(guessnewpath)
                    else:
                        to_remove.append(key)
                else:
                    to_remove.append(key)
        for i in to_remove:
            del self.rkfs[i]

    def engine_names(self) -> List[str]:
        """Return the identifiers of all engine-specific RKF files.

        :return: engine result identifiers, excluding the main ``ams.rkf`` file
        """
        self.refresh()
        ret = list(self.rkfs.keys())
        ret.remove("ams")
        return ret

    def get_main_engine_name(self) -> str:
        """Return the identifier of the main engine result file.

        Geometry-optimization step files and Hybrid subengine files are excluded. For molecular
        dynamics, the most recent ``MDStep`` file is selected.

        :return: main engine result identifier
        :raises ValueError: if no unique main engine result file can be determined
        """
        engine_names = self.engine_names()
        original_task = str(self.job.get_task()).lower()

        # if GO allows to save extra .rkf files
        if original_task == "geometryoptimization":
            engine_names = [x for x in engine_names if "GOStep" not in x]

        # remove hybrid engine sub engines
        engine_names = [x for x in engine_names if "hybrid-" not in x]

        # if MD find most recent MDStep
        if original_task == "moleculardynamics":
            engine_names = sorted(
                [x for x in engine_names if "term" not in x],
                key=lambda x: int(m.group(1)) if (m := re.match(r"[a-zA-Z]*(\d+)", x)) else -1,
            )
            engine_names = [engine_names[-1]]

        if len(engine_names) != 1:
            raise ValueError(
                f"Cannot get main engine name from {engine_names} for job in: {self.job.path} with {list(self.rkfs.keys())}"
            )
        return engine_names[0]

    def read_hybrid_term_rkf(self, section: str, variable: str, term: int, file: str = "engine") -> "TRead":
        """Read a variable from a Hybrid term's subengine RKF file.

        :param section: name of the RKF section
        :param variable: name of the RKF variable
        :param term: one-based Hybrid term index
        :param file: identifier of the Hybrid engine RKF file
        :return: value stored in the selected variable, without unit conversion; its type, shape, and unit depend on
            the variable
        """

        kf_file = cast(str, self.readrkf("EngineResults", f"Files({term})", file=file))
        kf = KFReader(str(self.job.get_path() / kf_file))
        return kf.read(section, variable)

    def rkfpath(self, file: str = "ams") -> str:
        """Return the absolute path of a chosen ``.rkf`` file.

        :param file: RKF file identifier. Use ``engine`` to select the unique engine result file.
        :return: absolute path to the selected RKF file
        """
        return self._access_rkf(lambda x: x.path, file)

    def readrkf(self, section: str, variable: str, file: str = "ams") -> "TRead":
        """Read a variable from a section of an RKF file.

        :param section: name of the RKF section; section names are case-sensitive
        :param variable: name of the RKF variable; variable names are case-sensitive
        :param file: RKF file identifier, defaults to ``ams``. For example, use ``something`` to read
            ``something.rkf``. If the job has a unique engine results file, use ``engine`` to select it.
        :return: value stored in the RKF variable, without unit conversion; its type, shape, and unit depend on the
            variable
        :raises KeyError: if the section or variable is not present in the selected RKF file
        """
        return self._access_rkf(lambda x: x.read(section, variable), file)

    def read_rkf_section(self, section: str, file: str = "ams") -> Dict[str, "TRead"]:
        """Return all variables from an RKF section.

        Section and variable names are case-sensitive. A missing section produces an empty dictionary.

        :param section: name of the RKF section
        :param file: RKF file identifier. Use ``engine`` to select the unique engine result file.
        :return: mapping from variable names to stored values, without unit conversion; each value's type, shape, and
            unit depend on the variable
        """
        return self._access_rkf(lambda x: x.read_section(section), file)

    def get_rkf_skeleton(self, file: str = "ams") -> Dict[str, Set[str]]:
        """Return the section and variable structure of an RKF file.

        :param file: RKF file identifier. Use ``engine`` to select the unique engine result file.
        :return: mapping from section names to sets of variable names
        """
        return self._access_rkf(lambda x: x.get_skeleton(), file)

    def get_molecule(self, section: str, file: str = "ams") -> Molecule:
        """Return a |Molecule| stored in an RKF section.

        All data used by this method is taken from the chosen ``.rkf`` file. The ``molecule`` attribute of the corresponding job is ignored.
        Coordinates and lattice vectors in the returned molecule are in angstrom.

        :param section: name of the RKF section containing the molecule
        :param file: RKF file identifier. Use ``engine`` to select the unique engine result file.
        :return: molecule reconstructed from the RKF section
        :raises ValueError: if the section does not contain a molecule
        """
        sectiondict = self.read_rkf_section(section, file)
        if sectiondict:
            return Molecule._mol_from_rkf_section(sectiondict)
        else:
            raise ValueError(f"Could not retrieve molecule from rkf section '{section}' for file '{file}.rkf'")

    @requires_optional_package("scm.base")
    def get_system(self, section: str, file: str = "ams") -> "ChemicalSystem":
        """Return a ``ChemicalSystem`` stored in an RKF section.

        All data used by this method is taken from the chosen ``.rkf`` file. The ``molecule`` attribute of the corresponding job is ignored.

        Note that ``ChemicalSystem`` is only available within AMS python. If unavailable, the call will raise an error.
        Coordinates and lattice vectors in the returned system are in angstrom.

        :param section: name of the RKF section containing the system
        :param file: RKF file identifier. Use ``engine`` to select the unique engine result file.
        :return: chemical system reconstructed from the RKF section
        """
        return ChemicalSystem.from_kf(self.rkfpath(file), section)

    @requires_optional_package("ase")
    def get_ase_atoms(self, section: str, file: str = "ams") -> "AseAtoms":
        """Return an ASE ``Atoms`` object stored in an RKF section.

        Positions and cell vectors are in angstrom.

        :param section: name of the RKF section containing the system
        :param file: RKF file identifier. Use ``engine`` to select the unique engine result file.
        :return: ASE atoms reconstructed from the RKF section
        """
        from ase import Atoms

        sectiondict = self.read_rkf_section(section, file)
        bohr2angstrom = Units.conversion_ratio("bohr", "angstrom")
        nLatticeVectors = cast(int, sectiondict.get("nLatticeVectors", 0))
        pbc = [True] * nLatticeVectors + [False] * (3 - nLatticeVectors)
        if nLatticeVectors > 0:
            cell = np.zeros((3, 3))
            lattice = np.array(sectiondict.get("LatticeVectors")).reshape(-1, 3)
            cell[: lattice.shape[0], : lattice.shape[1]] = lattice * bohr2angstrom
        else:
            cell = None
        atomsymbols = cast(str, sectiondict["AtomSymbols"]).split()
        positions = np.array(sectiondict["Coords"]).reshape(-1, 3) * bohr2angstrom
        return Atoms(symbols=atomsymbols, positions=positions, pbc=pbc, cell=cell)

    def get_input_molecule(self) -> Molecule:
        """Return a |Molecule| instance with the initial coordinates.

        All data used by this method is taken from ``ams.rkf`` file.
        The ``molecule`` attribute of the corresponding job is ignored.

        :return: input structure reconstructed from ``ams.rkf``, with coordinates and lattice vectors in angstrom
        """
        return self.get_molecule("InputMolecule", "ams")

    def get_input_system(self) -> "ChemicalSystem":
        """Return a ``ChemicalSystem`` instance with the initial coordinates.

        All data used by this method is taken from ``ams.rkf`` file.
        The ``molecule`` attribute of the corresponding job is ignored.

        Note that ``ChemicalSystem`` is only available within AMS python. If unavailable, the call will raise an error.

        :return: input system reconstructed from ``ams.rkf``, with coordinates and lattice vectors in angstrom
        """
        return self.get_system("InputMolecule", "ams")

    def get_input_molecules(self) -> Dict[str, Molecule]:
        """Return a dictionary mapping the name strings from the AMS input file to |Molecule| instances.

        The main molecule (aka the one that did not have a string in the block header in the input file) will be
        returned under the key of the empty string "". All data used by this method is taken from ``ams.rkf`` file.
        The ``molecule`` attribute of the corresponding job is ignored.

        :return: input molecules keyed by their System-block names, with coordinates and lattice vectors in angstrom;
            the unnamed main system uses an empty key
        """
        mols: Dict[str, Molecule] = {}
        skel = self.get_rkf_skeleton()
        if "InputMolecule" in skel:
            mols[""] = self.get_molecule("InputMolecule")
        if "InputMolecules" in skel:
            num_named_molecules = cast(int, self.readrkf("InputMolecules", "numNamedMolecules"))
            for imol in range(1, num_named_molecules + 1):
                mols[cast(str, self.readrkf("InputMolecules", f"Name({imol})"))] = self.get_molecule(
                    f"InputMolecule({imol})"
                )
        return mols

    def get_main_molecule(self) -> Molecule:
        """Return a |Molecule| instance with the final coordinates.

        All data used by this method is taken from ``ams.rkf`` file.
        The ``molecule`` attribute of the corresponding job is ignored.

        :return: final structure reconstructed from ``ams.rkf``, with coordinates and lattice vectors in angstrom
        """
        return self.get_molecule("Molecule", "ams")

    def get_main_system(self) -> "ChemicalSystem":
        """Return a ``ChemicalSystem`` instance with the final coordinates.

        All data used by this method is taken from ``ams.rkf`` file.
        The ``molecule`` attribute of the corresponding job is ignored.

        Note that ``ChemicalSystem`` is only available within AMS python. If unavailable, the call will raise an error.

        :return: final system reconstructed from ``ams.rkf``, with coordinates and lattice vectors in angstrom
        """
        return self.get_system("Molecule", "ams")

    @requires_optional_package("ase")
    def get_main_ase_atoms(self, get_results: bool = False) -> "AseAtoms":
        """Return an ASE ``Atoms`` instance with the final coordinates.

        Positions and cell vectors are in angstrom. When results are attached, ASE conventions are
        used: energy is in eV, forces in eV/angstrom, stress in eV/angstrom^3, and charges in
        elementary-charge units.

        :param get_results: attach available energy, forces, stress, and charges through an ASE
            ``SinglePointCalculator``
        :return: final structure as an ASE ``Atoms`` instance
        """
        from ase.calculators.singlepoint import SinglePointCalculator

        atoms = self.get_ase_atoms("Molecule", "ams")
        if get_results:
            atoms.set_calculator(SinglePointCalculator(atoms))
            try:
                energy = self.get_energy(unit="eV")
                atoms.calc.results["energy"] = energy
            except KeyError:
                pass
            try:
                forces = (-1) * np.array(self.get_gradients(dist_unit="angstrom", energy_unit="eV")).reshape(-1, 3)
                atoms.calc.results["forces"] = forces
            except KeyError:
                pass

            try:
                stress = np.array(self.get_stresstensor()).ravel() * Units.convert(1.0, "hartree/bohr^3", "eV/ang^3")
                if len(stress) == 9:
                    stress = stress.reshape(3, 3)
                    atoms.calc.results["stress"] = np.array(
                        [stress[0][0], stress[1][1], stress[2][2], stress[1][2], stress[0][2], stress[0][1]]
                    )
            except KeyError:
                pass

            try:
                charges = np.array(self.get_charges()).ravel()
                atoms.calc.results["charges"] = charges
            except KeyError:
                pass

        return atoms

    def get_history_molecule(self, step: int) -> Optional[Molecule]:
        """Return the molecule at a trajectory step from the ``History`` section of ``ams.rkf``.

        All data used by this method is taken from ``ams.rkf`` file. The ``molecule`` attribute of the corresponding job is ignored.
        Coordinates and lattice vectors in the returned molecule are in angstrom.

        :param step: one-based history step number
        :return: molecule at the requested step, or ``None`` if ``ams.rkf`` is unavailable
        :raises KeyError: if the requested step is outside the stored history
        :raises ResultsError: if the stored coordinates and molecular system have different sizes
        """
        if "ams" not in self.rkfs:
            return None

        main = self.rkfs["ams"]
        if not self.is_valid_stepnumber(main, step):
            return None
        flat_coords = main.read_reals("History", f"Coords({step})")
        coords = [flat_coords[i : i + 3] for i in range(0, len(flat_coords), 3)]
        if ("History", f"SystemVersion({step})") in main:
            system = self.get_system_version(main, step)
            mol = self.get_molecule(f"ChemicalSystem({system})")
            molsrc = f"ChemicalSystem({system})"
        else:
            mol = self.get_main_molecule()
            molsrc = "Molecule"
        assert mol is not None
        if len(mol) != len(coords):
            raise ResultsError(
                f'Coordinates taken from "History%Coords({step})" have incompatible length with molecule from {molsrc} section'
            )
        for at, c in zip(mol, coords):
            at.move_to(c, unit="bohr")

        if ("History", "LatticeVectors(" + str(step) + ")") in main:
            lattice = Units.convert(main.read_reals("History", "LatticeVectors(" + str(step) + ")"), "bohr", "angstrom")
            mol.lattice = [lattice[j : j + 3] for j in range(0, len(lattice), 3)]

        # Bonds from the reference molecule are probably outdated. Let us never use them ...
        mol.delete_all_bonds()
        # ... but instead use bonds from the history section if they are available:
        if all(
            ("History", i) in main for i in [f"Bonds.Index({step})", f"Bonds.Atoms({step})", f"Bonds.Orders({step})"]
        ):
            index = main.read_ints("History", f"Bonds.Index({step})")
            atoms = main.read_ints("History", f"Bonds.Atoms({step})")
            orders = main.read_reals("History", f"Bonds.Orders({step})")
            if not isinstance(orders, list):
                orders = [orders]
            for i in range(len(index) - 1):
                for j in range(index[i], index[i + 1]):
                    mol.add_bond(mol[i + 1], mol[atoms[j - 1]], orders[j - 1])
        if ("History", f"Bonds.CellShifts({step})") in main:
            assert mol.lattice
            cellShifts = main.read_ints("History", f"Bonds.CellShifts({step})")
            ndim = len(mol.lattice)
            for i, b in enumerate(mol.bonds):
                b.properties.suffix = " ".join([f"{cellShifts[ndim*i+j]}" for j in range(min(len(mol.lattice), ndim))])

        return mol

    def is_valid_stepnumber(self, main: KFFile, step: int) -> bool:
        """Check whether a one-based step number is present in an RKF history.

        :param main: RKF file containing the ``History`` section
        :param step: one-based history step number
        :return: ``True`` when the step is present
        :raises KeyError: if the history or requested step is absent
        """
        if "History" not in main:
            raise KeyError("'History' section not present in {}".format(main.path))
        n = main.read_int("History", "nEntries")
        if step > n or step <= 0:
            raise KeyError("Step {} not present in 'History' section of {}".format(step, main.path))
        return True

    def get_system_version(self, main: KFFile, step: int) -> Optional[int]:
        """Return the chemical-system section number used at a history step.

        :param main: RKF file containing the history
        :param step: one-based history step number
        :return: chemical-system section number, or ``None`` if no version is stored
        """
        if ("History", f"SystemVersion({step})") in main:
            version = main.read_int("History", f"SystemVersion({step})")
            if "SystemVersionHistory" in main:
                if ("SystemVersionHistory", "blockSize") in main:
                    blockSize = main.read_int("SystemVersionHistory", "blockSize")
                else:
                    blockSize = 1
                block = (version - 1) // blockSize + 1
                offset = (version - 1) % blockSize
                system = main.read_ints("SystemVersionHistory", f"SectionNum({block})")[offset]
            else:
                system = version
            return system
        else:
            return None

    def get_history_variables(self, history_section: str = "History") -> Optional[Set[str]]:
        """Return the variable names stored in an ``ams.rkf`` history section.

        :param history_section: history section name, such as ``History`` or ``MDHistory``
        :return: variable names, or ``None`` if ``ams.rkf`` is unavailable
        """
        if "ams" not in self.rkfs:
            return None
        main = self.rkfs["ams"]
        keylist = [var for sec, var in main if sec == history_section]
        # Now throw out all the last parts
        return set([key.split("(")[0] for key in keylist if len(key.split("(")) > 1])

    def get_history_length(self, history_section: str = "History") -> int:
        """Return the number of entries in an ``ams.rkf`` history section.

        :param history_section: history section name
        :return: value of ``nEntries`` in the selected section
        """
        return cast(int, self.readrkf(history_section, "nEntries"))

    def get_history_property(self, varname: str, history_section: str = "History") -> Optional[List]:
        """Return all saved values of a history variable.

        :param varname: variable name without the entry index
        :param history_section: history section name
        :return: stored values for all entries, without unit conversion, or ``None`` if the file or entry count is
            unavailable; each value's type, shape, and unit depend on the variable
        :raises KeyError: if the requested history section does not exist
        """
        if "ams" not in self.rkfs:
            return None
        main = self.rkfs["ams"]
        if history_section not in main:
            raise KeyError(f"The requested section '{history_section}' does not exist in {main.path}")
        if (history_section, "nScanCoord") in main:  # PESScan
            nentries = main.read_int(history_section, "nScanCoord")
            as_block = False
        elif (history_section, "nEntries") in main:
            nentries = main.read_int(history_section, "nEntries")
            as_block = self._values_stored_as_blocks(main, varname, history_section)
        else:
            return None
        if as_block:
            nblocks = main.read_int(history_section, "nBlocks")
            values = [
                main.read(history_section, f"{varname}({iblock})", return_as_list=True)
                for iblock in range(1, nblocks + 1)
            ]
            values = [val for blockvals in values for val in blockvals if isinstance(blockvals, list)]  # type: ignore[union-attr]
        else:
            values = [main.read(history_section, f"{varname}({step})") for step in range(1, nentries + 1)]  # type: ignore[misc]
        return values

    def get_property_at_step(self, step: int, varname: str, history_section: str = "History") -> Optional["TRead"]:
        """Return a history variable at one step.

        :param step: one-based history step number
        :param varname: variable name without the entry index
        :param history_section: history section name
        :return: stored value without unit conversion, or ``None`` if ``ams.rkf`` is unavailable; its type, shape,
            and unit depend on the variable
        """
        if "ams" not in self.rkfs:
            return None
        main = self.rkfs["ams"]
        as_block = self._values_stored_as_blocks(main, varname, history_section)
        if as_block:
            blocksize = main.read_int(history_section, "blockSize")
            iblock = int(np.ceil(step / blocksize))
            value = main.read(
                history_section, f"{varname}({iblock})"
            )  # this can return something that isn't a list, for example an int
            try:
                value = value[(step % blocksize) - 1]  # type: ignore[index]
            except TypeError:  # TypeError: 'int' object is not subscriptable
                pass
        else:
            value = main.read(history_section, f"{varname}({step})")
        return value

    def _values_stored_as_blocks(self, main: KFFile, varname: str, history_section: str) -> bool:
        """Return whether a trajectory variable is stored in blocks.

        :param main: trajectory RKF file
        :param varname: history variable name
        :param history_section: history section name
        :return: whether the values use block storage
        """
        nentries: int = main.read_int(history_section, "nEntries")
        as_block = False
        # This is extremely slow, because looping over main is very slow.
        # This is because the (sec,var) tuples are first stored in a set, then sorted, and only then yielded
        keylist = [var for sec, var in main if sec == history_section]
        if "nBlocks" in keylist:
            if not f"{varname}({nentries})" in keylist:
                as_block = True
        return as_block

    def get_atomic_temperatures_at_step(self, step: int, history_section: str = "MDHistory") -> Optional[np.ndarray]:
        """Return the atomic temperatures at a trajectory step.

        :param step: one-based trajectory step number
        :param history_section: section containing the velocities
        :return: array of shape ``(n_atoms,)`` containing temperatures in kelvin, or ``None`` if ``ams.rkf`` is
            unavailable
        """
        if not "ams" in self.rkfs:
            return None
        main = self.rkfs["ams"]
        if not self.is_valid_stepnumber(main, step):
            return None

        # Read the masses
        molname = "Molecule"
        if ("History", f"SystemVersion({step})") in main:
            system = self.get_system_version(main, step)
            molname = f"ChemicalSystem({system})"
        masses = np.array(main.read_reals(molname, "AtomMasses"))
        nats = len(masses)

        # Read the velocities
        velocities = np.array(self.get_property_at_step(step, "Velocities", history_section)).reshape((nats, 3))

        # Convert to SI units and compute temperatures
        m = masses * 1.0e-3 / Units.constants["NA"]
        vels = velocities * Units.conversion_ratio("Bohr", "Angstrom") * 1.0e5
        temperatures = (m.reshape((nats, 1)) * vels**2).sum(axis=1)
        temperatures /= 3 * Units.constants["k_B"]
        return temperatures

    def get_band_structure(
        self, bands: Optional[Sequence[int]] = None, unit: str = "hartree", only_high_symmetry_points: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], float]:
        """Return the electronic band structure from a DFTB, BAND, or Quantum ESPRESSO calculation.

        For unrestricted calculations, alternating stored bands belong to the spin-up and spin-down
        channels. The returned values can be passed to ``plot_band_structure``.

        :param bands: zero-based band indices to return, or ``None`` for all bands
        :param unit: energy unit for band and Fermi energies, defaults to ``hartree``
        :param only_high_symmetry_points: return only the first point of each path edge
        :return: path coordinates, spin-up energies, spin-down energies, point labels, and Fermi energy.
            Energy arrays have shape ``(len(x), len(bands))``; for restricted calculations the two
            arrays are identical. Labels are empty for points that are not high-symmetry points.
        """

        read_labels = True

        nBands = cast(int, self.readrkf("band_curves", "nBands", file="engine"))

        nEdges = cast(int, self.readrkf("band_curves", "nEdges", file="engine"))
        try:
            nSpin = cast(int, self.readrkf("band_curves", "nSpin", file="engine"))
        except KeyError:
            nSpin = 1

        if bands is None:
            bands = np.arange(nBands).tolist()

        spindown_bands = np.array(bands)
        if nSpin == 2:
            spindown_bands = np.array(bands) + nBands

        prevmaxx = 0
        complete_spinup_data = []
        complete_spindown_data = []
        labels = []
        x = []

        for i in range(nEdges):
            my_x = np.array(self.readrkf("band_curves", f"Edge_{i+1}_xFor1DPlotting", file="engine"))
            my_x += prevmaxx
            prevmaxx = np.max(my_x)

            if read_labels:
                my_labels = cast(str, self.readrkf("band_curves", f"Edge_{i+1}_labels", file="engine")).split()
                if len(my_labels) == 2:
                    if only_high_symmetry_points:
                        labels += [my_labels[0]]  # only the first point of the curve
                    else:
                        labels += [my_labels[0]] + [""] * (len(my_x) - 2) + [my_labels[1]]

            if only_high_symmetry_points:
                x.append(my_x[0:1])
            else:
                x.append(my_x)

            tA = cast(List[float], self.readrkf("band_curves", f"Edge_{i+1}_bands", file="engine"))
            A = np.array(tA).reshape(-1, nBands * nSpin)
            spinup_data = A[:, bands]
            spindown_data = A[:, spindown_bands]

            if only_high_symmetry_points:
                spinup_data = np.reshape(spinup_data[0, :], (-1, len(bands)))
                spindown_data = np.reshape(spindown_data[0, :], (-1, len(spindown_bands)))

            complete_spinup_data.append(spinup_data)
            complete_spindown_data.append(spindown_data)

        complete_spinup_data_ret = np.concatenate(complete_spinup_data)
        complete_spinup_data_ret = Units.convert(complete_spinup_data_ret, "hartree", unit)

        complete_spindown_data_ret = np.concatenate(complete_spindown_data)
        complete_spindown_data_ret = Units.convert(complete_spindown_data_ret, "hartree", unit)

        x_ret = np.concatenate(x).ravel()

        fermi_energy = cast(float, self.readrkf("BandStructure", "FermiEnergy", file="engine"))
        fermi_energy = Units.convert(fermi_energy, "hartree", unit)

        return x_ret, complete_spinup_data_ret, complete_spindown_data_ret, labels, fermi_energy  # type: ignore[return-value]

    def get_phonons_dos(
        self, unit: str = "hartree"
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """Return the phonon density of states (DOS).

        The returned values can be passed directly to
        :func:`~scm.plams.tools.plot.plot_phonons_dos`.

        :param unit: energy unit for the energy grid, defaults to ``hartree``
        :return: tuple containing the one-dimensional energy grid in ``unit``, the total DOS in ``1 / unit``,
            the DOS per species in ``1 / unit``, and the DOS per atom in ``1 / unit``. The decomposed DOS values are
            dictionaries of one-dimensional arrays. They are empty when the corresponding decomposition is absent.
            Atom keys contain the species and one-based atom index.
        """
        energyConv = Units.convert(1.0, "Hartree", unit)

        nEnergies = cast(int, self.readrkf("DOS_Phonons", "nEnergies", file="engine"))

        integrateDeltaE = cast(float, self.readrkf("DOS_Phonons", "IntegrateDeltaE", file="engine"))

        assert integrateDeltaE

        energy = np.array(cast(List[float], self.readrkf("DOS_Phonons", "Energies", file="engine"))) * energyConv
        dE = energy[1] - energy[0]

        total_dos = np.array(cast(List[float], self.readrkf("DOS_Phonons", "Total DOS", file="engine"))) / dE

        assert len(energy) == nEnergies
        assert len(total_dos) == nEnergies

        dos_per_species = {}
        dos_per_atom = {}

        nSpecies = None
        try:
            nSpecies = cast(int, self.readrkf("DOS_Phonons", "nSpecies", file="engine"))
        except KeyError:
            print(
                "Warning! The density of states (DOS) per atom and species is currently only supported by the Quantum Espresso engine."
            )
            print("         Only the total DOS is available.")
            pass

        if nSpecies:

            nAtoms = cast(int, self.readrkf("DOS_Phonons", "nAtoms", file="engine"))
            assert self.job.molecule is not None
            assert nAtoms == len(self.job.molecule)

            DOSperSpecies = np.array(
                cast(List[int], self.readrkf("DOS_Phonons", "DOS per species", file="engine"))
            ).reshape(nSpecies, -1)
            species = cast(str, self.readrkf("DOS_Phonons", "Species", file="engine")).split("\0")

            for i, s in enumerate(species):
                dos_per_species[s] = DOSperSpecies[i] / dE

            DOSperAtom = np.array(cast(List[int], self.readrkf("DOS_Phonons", "DOS per atom", file="engine"))).reshape(
                nAtoms, -1
            )
            atomToSpecies = np.array(cast(List[int], self.readrkf("DOS_Phonons", "Atom to Species", file="engine"))) - 1
            if isinstance(atomToSpecies, np.int64):
                atomToSpecies = atomToSpecies.reshape(1)

            for i in range(nAtoms):
                label = species[atomToSpecies[i]] + "(" + str(i + 1) + ")"
                dos_per_atom[label] = DOSperAtom[i] / dE

        return energy, total_dos, dos_per_species, dos_per_atom

    def get_phonons_band_structure(
        self, bands: Optional[Sequence[int]] = None, unit: str = "hartree", only_high_symmetry_points: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, List[str]]:
        """Return the phonon band structure.

        The returned values can be passed directly to
        :func:`~scm.plams.tools.plot.plot_phonons_band_structure`.

        :param bands: zero-based band indices to return, defaults to all bands
        :param unit: energy unit for the phonon bands, defaults to ``hartree``
        :param only_high_symmetry_points: return only the first point of each path edge, defaults to ``False``
        :return: tuple containing the one-dimensional reciprocal-path plotting coordinate, an array of shape
            ``(len(x), len(bands))`` with the phonon band energies in ``unit``, and the labels for the points in ``x``.
            Each column in the energy array is a band. Non-high-symmetry points have an empty label.
        """

        read_labels = True

        nBands = cast(int, self.readrkf("phonon_curves", "nBands", file="engine"))
        nEdges = cast(int, self.readrkf("phonon_curves", "nEdges", file="engine"))

        if bands is None:
            bands = np.arange(nBands).tolist()

        x = []
        y = []
        labels = []

        prevmaxx = 0
        for i in range(nEdges):
            my_x = np.array(
                cast(List[float], self.readrkf("phonon_curves", f"Edge_{i+1}_xFor1DPlotting", file="engine"))
            )
            my_x += prevmaxx
            prevmaxx = np.max(my_x)

            if read_labels:
                my_labels = cast(str, self.readrkf("phonon_curves", f"Edge_{i+1}_labels", file="engine")).split()
                if len(my_labels) == 2:
                    if only_high_symmetry_points:
                        labels += [my_labels[0]]  # only the first point of the curve
                    else:
                        labels += [my_labels[0]] + [""] * (len(my_x) - 2) + [my_labels[1]]

            if only_high_symmetry_points:
                x.append(my_x[0:1])
            else:
                x.append(my_x)

            tA = cast(List[float], self.readrkf("phonon_curves", f"Edge_{i+1}_bands", file="engine"))
            A = np.array(tA).reshape(-1, nBands)
            spinup_data = A[:, bands]

            if only_high_symmetry_points:
                spinup_data = np.reshape(spinup_data[0, :], (-1, len(bands)))

            y.append(spinup_data)

        y_ret = np.concatenate(y)
        y_ret = Units.convert(y_ret, "hartree", unit)

        x_ret = np.concatenate(x).ravel()

        return x_ret, y_ret, labels  # type: ignore[return-value]

    def get_phonons_thermodynamic_properties(
        self, temperature_unit: str = "K", properties_unit: Sequence[str] = ("hartree", "kB")
    ) -> Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, str]]:
        """Return phonon thermodynamic properties as functions of temperature.

        :param temperature_unit: unit for the returned temperatures, defaults to ``K``. This may be ``K`` or an
            energy-equivalent unit supported by :class:`~scm.plams.tools.units.Units`.
        :param properties_unit: two-element sequence containing the energy unit for internal and free energies,
            followed by the unit for specific heat and entropy, defaults to ``("hartree", "kB")``. The second value
            may be ``kB``, ``J/K``, ``J/K/mol``, ``eV/K``, or ``eV/K/mol``.
        :return: tuple containing a one-dimensional temperature array in ``temperature_unit``, a dictionary of
            property arrays, and a dictionary specifying the unit of each property. Possible property keys are
            ``Internal Energy``, ``Free Energy``, ``Specific Heat``, and ``Entropy``.
        """
        properties = {}
        units = {}

        energyConv = Units.convert(1.0, "Hartree", properties_unit[0])

        entropyConv = 1.0
        if properties_unit[1].lower() == "j/k":
            entropyConv = Units.constants["Boltzmann"]
        if properties_unit[1].lower() == "ev/k":
            entropyConv = Units.constants["Boltzmann"] / Units.constants["electron_charge"]
        elif properties_unit[1].lower() == "j/k/mol":
            entropyConv = Units.constants["Boltzmann"] * Units.constants["Avogadro_constant"]
        elif properties_unit[1].lower() == "ev/k/mol":
            entropyConv = (
                Units.constants["Boltzmann"] * Units.constants["Avogadro_constant"] / Units.constants["electron_charge"]
            )

        nPlots = cast(int, self.readrkf("Plot", "numPlots", file="engine"))

        temperature = np.array(
            cast(List[float], self.readrkf("Plot", "XValues(1)", file="engine"))
        )  # We assume that temperature is the same for all properties
        temperature = Units.convert(temperature, "K", temperature_unit)

        for iProp in range(nPlots):
            label = cast(str, self.readrkf("Plot", f"YLabel({iProp+1})", file="engine")).strip()

            if label not in ["Internal Energy", "Free Energy", "Specific Heat", "Entropy"]:
                continue

            values = np.array(cast(List[float], self.readrkf("Plot", f"YValues({iProp+1})", file="engine")))

            property_title = cast(str, self.readrkf("Plot", f"Title({iProp+1})", file="engine")).strip()
            stored_temperature_unit = cast(str, self.readrkf("Plot", f"XUnit({iProp+1})", file="engine")).strip()
            property_unit = cast(str, self.readrkf("Plot", f"YUnit({iProp+1})", file="engine")).strip()
            temperature_label = cast(str, self.readrkf("Plot", f"XLabel({iProp+1})", file="engine")).strip()

            assert property_title == label
            assert temperature_label == "Temperature"
            assert stored_temperature_unit == "Kelvin"
            if label in ["Internal Energy", "Free Energy"]:
                assert property_unit == "Hartree"
                units[label] = properties_unit[0]
                values *= energyConv
            if label in ["Specific Heat", "Entropy"]:
                assert property_unit == "kB"
                units[label] = properties_unit[1]
                values *= entropyConv

            properties[label] = values

        return temperature, properties, units

    def get_engine_results(self, engine: Optional[str] = None) -> Dict:
        """Return the ``AMSResults`` section of an engine RKF file.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: mapping from variable names to stored values, without unit conversion; each value's type, shape, and
            unit depend on the variable
        """
        return self._process_engine_results(lambda x: x.read_section("AMSResults"), engine)

    def get_engine_properties(self, engine: Optional[str] = None) -> Dict:
        """Return the ``Properties`` section of an engine RKF file.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: mapping from property names to stored values, without unit conversion; each value's type, shape, and
            unit depend on the property
        """

        def properties(kf: KFFile) -> Dict[str, "TRead"]:
            if not "Properties" in kf:
                return {}
            # There are two *kinds* of "Properties" sections:
            # - ADF simply has a set of variables in the properties section
            # - DFTB and some other engines have a *scheme* for storing "Properties". See the fortran 'RKFileModule' (RKFile.f90)
            #   for more details on this.
            ret = {}
            if ("Properties", "nEntries") in kf:
                # This is a 'RKFileModule' properties section:
                n = kf.read_int("Properties", "nEntries")
                for i in range(1, n + 1):
                    tp = kf.read_string("Properties", "Type({})".format(i)).strip()
                    stp = kf.read_string("Properties", "Subtype({})".format(i)).strip()
                    val = kf.read("Properties", "Value({})".format(i))
                    key = stp if stp.endswith(tp) else ("{} {}".format(stp, tp) if stp else tp)
                    ret[key] = val
            else:
                skeleton = kf.get_skeleton()
                for variable in skeleton["Properties"]:
                    ret[variable] = kf.read("Properties", variable)
            return ret

        return self._process_engine_results(properties, engine)

    def get_energy(self, unit: str = "hartree", engine: Optional[str] = None) -> float:
        """Return the final energy from ``AMSResults%Energy``.

        The meaning of the final energy depends on the engine; see the corresponding engine documentation.

        :param unit: energy unit of the returned value, defaults to ``hartree``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: final energy as a scalar in ``unit``
        """
        return self._process_engine_results(
            lambda x: x.read_real("AMSResults", "Energy"), engine
        ) * Units.conversion_ratio("au", unit)

    def get_energy_uncertainty(self, unit: str = "hartree", engine: Optional[str] = None) -> float:
        """Return the final energy uncertainty from ``AMSResults%EnergyUncertainty``.

        The meaning of the uncertainty depends on the engine; see the corresponding engine documentation.

        :param unit: energy unit of the returned value, defaults to ``hartree``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: final energy uncertainty as a scalar in ``unit``
        """
        return self._process_engine_results(
            lambda x: x.read_real("AMSResults", "EnergyUncertainty"), engine
        ) * Units.conversion_ratio("au", unit)

    def get_gradients(
        self, energy_unit: str = "hartree", dist_unit: str = "bohr", engine: Optional[str] = None
    ) -> np.ndarray:
        """Return the final nuclear energy gradients from ``AMSResults%Gradients``.

        These values are energy gradients, not forces. Forces are the negative gradients.

        :param energy_unit: energy unit of the returned gradients, defaults to ``hartree``
        :param dist_unit: distance unit of the returned gradients, defaults to ``bohr``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: gradient array of shape ``(n_atoms, 3)`` in ``energy_unit / dist_unit``
        """
        return (
            np.asarray(self._process_engine_results(lambda x: x.read_reals("AMSResults", "Gradients"), engine)).reshape(
                -1, 3
            )
            * Units.conversion_ratio("au", energy_unit)
            / Units.conversion_ratio("au", dist_unit)
        )

    def get_gradients_uncertainty(
        self, energy_unit: str = "hartree", dist_unit: str = "bohr", engine: Optional[str] = None
    ) -> np.ndarray:
        """Return the final-gradient uncertainties from ``AMSResults%GradientsUncertainty``.

        :param energy_unit: energy unit of the returned uncertainties, defaults to ``hartree``
        :param dist_unit: distance unit of the returned uncertainties, defaults to ``bohr``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: uncertainty array of shape ``(n_atoms, 3)`` in ``energy_unit / dist_unit``
        """
        return (
            np.asarray(
                self._process_engine_results(lambda x: x.read_reals("AMSResults", "GradientsUncertainty"), engine)
            ).reshape(-1, 3)
            * Units.conversion_ratio("au", energy_unit)
            / Units.conversion_ratio("au", dist_unit)
        )

    def get_gradients_magnitude_uncertainty(
        self, energy_unit: str = "hartree", dist_unit: str = "bohr", engine: Optional[str] = None
    ) -> np.ndarray:
        """Return the propagated uncertainty of each final-gradient magnitude.

        The values are read from ``AMSResults%GradientsMagnitudeUncertainty``.

        :param energy_unit: energy unit of the returned uncertainties, defaults to ``hartree``
        :param dist_unit: distance unit of the returned uncertainties, defaults to ``bohr``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: uncertainty array of shape ``(n_atoms,)`` in ``energy_unit / dist_unit``
        """
        return (
            np.asarray(
                self._process_engine_results(
                    lambda x: x.read_reals("AMSResults", "GradientsMagnitudeUncertainty"), engine
                )
            ).reshape(-1)
            * Units.conversion_ratio("au", energy_unit)
            / Units.conversion_ratio("au", dist_unit)
        )

    def get_stresstensor(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the clamped-ion Cartesian stress tensor from ``AMSResults%StressTensor``.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_lattice_vectors, n_lattice_vectors)`` in
            ``hartree/bohr^n_lattice_vectors``. This is energy per length for 1D periodicity, energy per area for 2D
            periodicity, and pressure for 3D periodicity.
        """
        return np.asarray(
            self._process_engine_results(lambda x: x.read_reals("AMSResults", "StressTensor"), engine)
        ).reshape(len(self.get_input_molecule().lattice), -1)

    def get_hessian(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the Cartesian nuclear Hessian from ``AMSResults%Hessian``.

        The matrix contains ordinary second derivatives of the total energy and is not mass-weighted.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(3 * n_atoms, 3 * n_atoms)`` in ``hartree/bohr^2``
        """
        return np.asarray(
            self._process_engine_results(lambda x: x.read_reals("AMSResults", "Hessian"), engine)
        ).reshape(3 * len(self.get_input_molecule()), -1)

    def get_elastictensor(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the elastic tensor in Voigt notation from ``AMSResults%ElasticTensor``.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array with shape ``(1, 1)``, ``(3, 3)``, or ``(6, 6)`` for 1D, 2D, or 3D periodicity,
            respectively, in ``hartree/bohr^n_lattice_vectors``
        """
        et_flat = np.asarray(
            self._process_engine_results(lambda x: x.read_reals("AMSResults", "ElasticTensor"), engine)
        )
        num_latvec = len(self.get_input_molecule().lattice)
        if num_latvec == 1:
            return et_flat.reshape(1, 1)
        elif num_latvec == 2:
            return et_flat.reshape(3, 3)
        else:
            return et_flat.reshape(6, 6)

    def get_frequencies(self, unit: str = "cm^-1", engine: Optional[str] = None) -> np.ndarray:
        """Return the vibrational frequencies.

        :param unit: unit for the returned frequencies, defaults to ``cm^-1``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_normal_modes,)`` containing frequencies in ``unit``
        """
        freqs = self.readrkf("Vibrations", "Frequencies[cm-1]", file=engine or "engine")
        freqs_arr = np.array(freqs).reshape((-1,))
        return freqs_arr * Units.conversion_ratio("cm^-1", unit)

    def get_frequency_spectrum(
        self,
        engine: Optional[str] = None,
        broadening_type: Literal["gaussian", "lorentzian"] = "gaussian",
        broadening_width: int = 40,
        min_x: int = 0,
        max_x: int = 4000,
        x_spacing: float = 0.5,
        post_process: Optional[Literal["max_to_1"]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return a broadened vibrational-frequency spectrum with unit mode intensities.

        Height-normalized profiles produce a dimensionless spectrum. Area-normalized profiles produce
        a spectrum in ``(cm^-1)^-1``. With ``max_to_1``, the spectrum is dimensionless.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :param broadening_type: line shape and whether unit weights define peak heights or integrated areas
        :param broadening_width: broadening width in ``cm^-1``
        :param min_x: lower frequency-grid bound in ``cm^-1``
        :param max_x: upper frequency-grid bound in ``cm^-1``
        :param x_spacing: frequency-grid spacing in ``cm^-1``
        :param post_process: normalize the spectrum maximum to one, or ``None``
        :return: frequency grid in ``cm^-1`` and broadened unit-weight intensities
        """
        from scm.plams.tools.postprocess_results import broaden_results

        frequencies = self.get_frequencies(engine=engine)
        intensities = frequencies * 0 + 1

        x_data, y_data = broaden_results(
            centers=frequencies,
            areas=intensities,
            broadening_width=broadening_width,
            broadening_type=broadening_type,
            x_data=(min_x, max_x, x_spacing),
            post_process=post_process,
        )
        return x_data, y_data

    def get_force_constants(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the vibrational force constants.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_normal_modes,)`` in ``hartree/bohr^2``
        """
        forceConstants = np.array(
            self._process_engine_results(lambda x: x.read_reals("Vibrations", "ForceConstants"), engine)
        )
        return forceConstants

    def get_pvdos(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the partial vibrational density of states (PVDOS).

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: dimensionless array of shape ``(n_normal_modes, n_atoms)`` with values between zero and one
        """
        pvdos = self._process_engine_results(lambda x: x.read_reals("Vibrations", "PVDOS"), engine)
        nNormalModes = self._process_engine_results(lambda x: x.read_int("Vibrations", "nNormalModes"), engine)
        nAtoms = len(self.get_main_molecule())
        return np.array(pvdos).reshape(nNormalModes, nAtoms)

    def get_reduced_masses(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the vibrational reduced masses.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_normal_modes,)`` in atomic mass units (amu)
        """
        reduced_masses = np.array(
            self._process_engine_results(lambda x: x.read_reals("Vibrations", "ReducedMasses"), engine)
        )
        return reduced_masses

    def get_normal_modes(
        self, engine: Optional[str] = None, mass_weighted_hessian_eigenvectors: Optional[bool] = False
    ) -> np.ndarray:
        """Return the normal-mode displacement vectors.

        By default, return the unweighted Cartesian displacement vectors stored in
        ``Vibrations%NoWeightNormalMode``. If ``mass_weighted_hessian_eigenvectors``
        is enabled, scale each atomic displacement by the square root of the ratio
        between its atomic mass and the mode's reduced mass.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :param mass_weighted_hessian_eigenvectors: return mass-weighted Hessian eigenvectors
        :return: dimensionless array of shape ``(n_normal_modes, n_atoms, 3)``
        """
        normal_modes_list = []
        num_normal_modes = self._process_engine_results(lambda x: x.read_int("Vibrations", "nNormalModes"), engine)
        for i in range(num_normal_modes):
            n_mode = np.array(
                self._process_engine_results(lambda x: x.read_reals("Vibrations", f"NoWeightNormalMode({i+1})"), engine)
            ).reshape(-1, 3)
            normal_modes_list.append(n_mode)
        normal_modes = np.array(normal_modes_list).reshape(num_normal_modes, -1, 3)
        if mass_weighted_hessian_eigenvectors:
            mol = self.get_main_molecule()
            masses = np.array(mol.get_masses()).reshape(1, -1, 1)
            reduced_masses = self.get_reduced_masses(engine=engine)
            reduced_masses = reduced_masses.reshape(-1, 1, 1)
            normal_modes_normalized = normal_modes * np.sqrt(masses) / np.sqrt(reduced_masses)
            return normal_modes_normalized
        return normal_modes

    def get_charges(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the atomic charges in elementary-charge units.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_atoms,)`` in elementary-charge units, ``e``
        """
        return np.asarray(self._process_engine_results(lambda x: x.read_reals("AMSResults", "Charges"), engine))

    def get_atom_types(self, engine: Optional[str] = None) -> List[str]:
        """Return the force-field atom type of each atom.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: atom types in system order
        """
        indices = self._process_engine_results(
            lambda x: x.read_ints("AMSResults", "AtomTyping.atomIndexToType"), engine
        )
        types = self._process_engine_results(
            lambda x: x.read_string("AMSResults", "AtomTyping.atomTypes"), engine
        ).split("\x00")

        return [types[i - 1] for i in indices]

    def get_dipolemoment(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the electric dipole moment in ``e*bohr``.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(3,)`` containing the Cartesian components in ``e*bohr``
        """
        return np.asarray(self._process_engine_results(lambda x: x.read_reals("AMSResults", "DipoleMoment"), engine))

    def get_dipolegradients(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the nuclear gradients of the electric dipole moment.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(3 * n_atoms, 3)`` in the atomic unit of dipole derivative, equivalent to ``e``
        """
        return np.asarray(
            self._process_engine_results(lambda x: x.read_reals("AMSResults", "DipoleGradients"), engine)
        ).reshape(-1, 3)

    def get_polarizability(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the polarizability in ``(e * bohr)^2 / hartree``.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: symmetric Cartesian tensor of shape ``(3, 3)``
        :raises ValueError: if the stored polarizability has an unsupported shape
        """
        p_components = np.asarray(
            self._process_engine_results(lambda x: x.read_reals("AMSResults", "Polarizability"), engine)
        )
        if p_components.shape == (6,):
            polarizability_matrix = np.array(
                [
                    [p_components[0], p_components[1], p_components[3]],
                    [p_components[1], p_components[2], p_components[4]],
                    [p_components[3], p_components[4], p_components[5]],
                ]
            )
        elif p_components.shape == (3, 3):  # type: ignore[comparison-overlap]
            polarizability_matrix = p_components
        else:
            raise ValueError(
                f"AMSResults-Polarizability shape is {p_components.shape} not in agreement with the option as inputs (6,) [xx,xy,yy,xz,zy,zz] or (3,3)"
            )
        return polarizability_matrix

    def get_zero_point_energy(self, unit: str = "hartree", engine: Optional[str] = None) -> float:
        """Return the vibrational zero-point energy from ``Vibrations%ZeroPointEnergy``.

        :param unit: energy unit of the returned value, defaults to ``hartree``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: zero-point energy as a scalar in ``unit``
        """
        return self._process_engine_results(
            lambda x: x.read_real("Vibrations", "ZeroPointEnergy"), engine
        ) * Units.conversion_ratio("au", unit)

    def get_ir_intensities(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the IR intensities in ``km/mol``.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_normal_modes,)`` in ``km/mol``
        """
        return np.asarray(self.readrkf("Vibrations", "Intensities[km/mol]", file=engine or "engine")).reshape((-1,))

    def get_raman_intensities(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the Raman intensities in ``angstrom^4/amu``.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_normal_modes,)`` in ``angstrom^4/amu``
        """
        return np.asarray(
            self._process_engine_results(lambda x: x.read_reals("Vibrations", "RamanIntens[A^4/amu]"), engine)
        ).reshape(
            -1,
        )

    def get_vcd_rotational_strength(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the vibrational rotational strengths in ``10^-44 esu^2 cm^2``.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_normal_modes,)`` in ``10^-44 esu^2 cm^2``
        """
        return np.asarray(
            self._process_engine_results(lambda x: x.read_reals("Vibrations", "RotationalStrength"), engine)
        ).reshape(
            -1,
        )

    def _get_ir_vcd_raman_spectrum(
        self,
        engine: Optional[str] = None,
        spectrum_type: Literal["ir", "raman", "vcd"] = "ir",
        broadening_type: Literal["gaussian", "lorentzian"] = "gaussian",
        broadening_width: int = 40,
        min_x: int = 0,
        max_x: int = 4000,
        x_spacing: float = 0.5,
        post_process: Optional[Literal["all_intensities_to_1", "max_to_1"]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return a broadened IR, Raman, or VCD spectrum.

        Height-normalized profiles preserve the unit of the line intensities. Area-normalized profiles
        return intensity per ``cm^-1``. With ``all_intensities_to_1``, unit line weights replace the
        physical intensities, giving dimensionless height profiles or area profiles in ``(cm^-1)^-1``.
        With ``max_to_1``, the returned spectrum is dimensionless.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :param spectrum_type: spectrum to calculate
        :param broadening_type: line-shape normalization and profile
        :param broadening_width: broadening width in ``cm^-1``
        :param min_x: lower frequency-grid bound in ``cm^-1``
        :param max_x: upper frequency-grid bound in ``cm^-1``
        :param x_spacing: frequency-grid spacing in ``cm^-1``
        :param post_process: replace line intensities with one, normalize the spectrum maximum, or ``None``
        :return: frequency grid in ``cm^-1`` and broadened intensities
        """
        from scm.plams.tools.postprocess_results import broaden_results

        frequencies = self.get_frequencies(engine=engine)
        if post_process == "all_intensities_to_1":
            intensities = frequencies * 0 + 1
        else:
            if spectrum_type == "ir":
                intensities = self.get_ir_intensities(engine=engine)
            elif spectrum_type == "raman":
                intensities = self.get_raman_intensities(engine=engine)
            elif spectrum_type == "vcd":
                intensities = self.get_vcd_rotational_strength(engine=engine)

        x_data, y_data = broaden_results(
            centers=frequencies,
            areas=intensities,
            broadening_width=broadening_width,
            broadening_type=broadening_type,
            x_data=(min_x, max_x, x_spacing),
            post_process=post_process if post_process == "max_to_1" else None,
        )

        return x_data, y_data

    def get_ir_spectrum(
        self,
        engine: Optional[str] = None,
        broadening_type: Literal["gaussian", "lorentzian"] = "gaussian",
        broadening_width: int = 40,
        min_x: int = 0,
        max_x: int = 4000,
        x_spacing: float = 0.5,
        post_process: Optional[Literal["all_intensities_to_1", "max_to_1"]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return a broadened IR spectrum.

        Frequencies are in ``cm^-1``. Height-normalized profiles produce intensities in ``km/mol``;
        area-normalized profiles produce intensities in ``(km/mol)/(cm^-1)``. With
        ``all_intensities_to_1``, the physical intensities are replaced by unit weights, giving a
        dimensionless height profile or an area profile in ``(cm^-1)^-1``. With ``max_to_1``, the
        spectrum is dimensionless.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :param broadening_type: line shape and whether line intensities define peak heights or integrated areas
        :param broadening_width: broadening width in ``cm^-1``
        :param min_x: lower frequency-grid bound in ``cm^-1``
        :param max_x: upper frequency-grid bound in ``cm^-1``
        :param x_spacing: frequency-grid spacing in ``cm^-1``
        :param post_process: replace line intensities with one, normalize the spectrum maximum, or ``None``
        :return: frequency grid in ``cm^-1`` and broadened intensities in the convention described above
        """
        data = self._get_ir_vcd_raman_spectrum(
            engine=engine,
            spectrum_type="ir",
            broadening_type=broadening_type,
            broadening_width=broadening_width,
            min_x=min_x,
            max_x=max_x,
            x_spacing=x_spacing,
            post_process=post_process,
        )

        return data

    def get_raman_spectrum(
        self,
        engine: Optional[str] = None,
        broadening_type: Literal["gaussian", "lorentzian"] = "gaussian",
        broadening_width: int = 40,
        min_x: int = 0,
        max_x: int = 4000,
        x_spacing: float = 0.5,
        post_process: Optional[Literal["all_intensities_to_1", "max_to_1"]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return a broadened Raman spectrum.

        Frequencies are in ``cm^-1``. Height-normalized profiles produce intensities in
        ``angstrom^4/amu``; area-normalized profiles produce intensities in
        ``(angstrom^4/amu)/(cm^-1)``. With ``all_intensities_to_1``, the physical intensities are
        replaced by unit weights, giving a dimensionless height profile or an area profile in
        ``(cm^-1)^-1``. With ``max_to_1``, the spectrum is dimensionless.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :param broadening_type: line shape and whether line intensities define peak heights or integrated areas
        :param broadening_width: broadening width in ``cm^-1``
        :param min_x: lower frequency-grid bound in ``cm^-1``
        :param max_x: upper frequency-grid bound in ``cm^-1``
        :param x_spacing: frequency-grid spacing in ``cm^-1``
        :param post_process: replace line intensities with one, normalize the spectrum maximum, or ``None``
        :return: frequency grid in ``cm^-1`` and broadened intensities in the convention described above
        """
        data = self._get_ir_vcd_raman_spectrum(
            engine=engine,
            spectrum_type="raman",
            broadening_type=broadening_type,
            broadening_width=broadening_width,
            min_x=min_x,
            max_x=max_x,
            x_spacing=x_spacing,
            post_process=post_process,
        )

        return data

    def get_vcd_spectrum(
        self,
        engine: Optional[str] = None,
        broadening_type: Literal["gaussian", "lorentzian"] = "gaussian",
        broadening_width: int = 40,
        min_x: int = 0,
        max_x: int = 4000,
        x_spacing: float = 0.5,
        post_process: Optional[Literal["all_intensities_to_1", "max_to_1"]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return a broadened VCD rotational-strength spectrum.

        Frequencies are in ``cm^-1``. Height-normalized profiles produce intensities in
        ``10^-44 esu^2 cm^2``; area-normalized profiles produce intensities in
        ``(10^-44 esu^2 cm^2)/(cm^-1)``. With ``all_intensities_to_1``, the physical intensities are
        replaced by unit weights, giving a dimensionless height profile or an area profile in
        ``(cm^-1)^-1``. With ``max_to_1``, the spectrum is dimensionless.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :param broadening_type: line shape and whether line intensities define peak heights or integrated areas
        :param broadening_width: broadening width in ``cm^-1``
        :param min_x: lower frequency-grid bound in ``cm^-1``
        :param max_x: upper frequency-grid bound in ``cm^-1``
        :param x_spacing: frequency-grid spacing in ``cm^-1``
        :param post_process: replace line intensities with one, normalize the spectrum maximum, or ``None``
        :return: frequency grid in ``cm^-1`` and broadened intensities in the convention described above
        """
        data = self._get_ir_vcd_raman_spectrum(
            engine=engine,
            spectrum_type="vcd",
            broadening_type=broadening_type,
            broadening_width=broadening_width,
            min_x=min_x,
            max_x=max_x,
            x_spacing=x_spacing,
            post_process=post_process,
        )

        return data

    def get_n_spin(self, engine: Optional[str] = None) -> int:
        """Return the number of stored spin channels.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: dimensionless integer; one for restricted or spin-orbit calculations and two for unrestricted
            calculations
        """
        return self._process_engine_results(lambda x: x.read_int("AMSResults", "nSpin"), engine)

    def get_orbital_energies(self, unit: str = "Hartree", engine: Optional[str] = None) -> np.ndarray:
        """Return the orbital energies for each spin channel.

        :param unit: energy unit for the returned values, defaults to ``Hartree``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_spin, n_orbitals)`` in ``unit``
        """
        return Units.convert(
            np.asarray(
                self._process_engine_results(lambda x: x.read_reals("AMSResults", "orbitalEnergies"), engine)
            ).reshape(self.get_n_spin(), -1),
            "Hartree",
            unit,
        )

    def get_orbital_occupations(self, engine: Optional[str] = None) -> np.ndarray:
        """Return the orbital occupations for each spin channel.

        Restricted occupations range from zero to two; unrestricted and spin-orbit occupations
        range from zero to one.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: dimensionless array of shape ``(n_spin, n_orbitals)``
        """
        return np.asarray(
            self._process_engine_results(lambda x: x.read_reals("AMSResults", "orbitalOccupations"), engine)
        ).reshape(self.get_n_spin(), -1)

    def get_homo_energies(self, unit: str = "Hartree", engine: Optional[str] = None) -> np.ndarray:
        """Return the HOMO energy of each spin channel.

        See also :meth:`are_orbitals_fractionally_occupied`.

        :param unit: energy unit for the returned values, defaults to ``Hartree``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_spin,)`` containing HOMO energies in ``unit``
        """
        return Units.convert(
            np.asarray(
                self._process_engine_results(lambda x: x.read_reals("AMSResults", "HOMOEnergy"), engine)
            ).reshape(-1),
            "Hartree",
            unit,
        )

    def get_lumo_energies(self, unit: str = "Hartree", engine: Optional[str] = None) -> np.ndarray:
        """Return the LUMO energy of each spin channel.

        See also :meth:`are_orbitals_fractionally_occupied`.

        :param unit: energy unit for the returned values, defaults to ``Hartree``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: array of shape ``(n_spin,)`` containing LUMO energies in ``unit``
        """
        return Units.convert(
            np.asarray(
                self._process_engine_results(lambda x: x.read_reals("AMSResults", "LUMOEnergy"), engine)
            ).reshape(-1),
            "Hartree",
            unit,
        )

    def get_smallest_homo_lumo_gap(self, unit: str = "Hartree", engine: Optional[str] = None) -> float:
        """Return the smallest HOMO-LUMO gap across spin channels.

        See also :meth:`are_orbitals_fractionally_occupied`.

        :param unit: energy unit for the returned gap, defaults to ``Hartree``
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: ``min(LUMO) - max(HOMO)`` in ``unit``
        """
        return Units.convert(
            self._process_engine_results(lambda x: x.read_real("AMSResults", "SmallestHOMOLUMOGap"), engine),
            "Hartree",
            unit,
        )

    def are_orbitals_fractionally_occupied(self, engine: Optional[str] = None) -> bool:
        """Return whether fractional orbital occupations were detected.

        In this case the HOMO and LUMO labels depend on an arbitrary occupied/empty demarcation.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: dimensionless boolean indicating whether fractional occupations are present
        """
        return self._process_engine_results(lambda x: x.read_logical("AMSResults", "fractionalOccupation"), engine)

    def get_timings(self) -> Dict[str, float]:
        """Return job timing statistics.

        :return: ``cpu``, ``system``, and ``elapsed`` times in seconds
        """
        ret: Dict[str, float] = {}
        try:
            # new AMS versions store timings on ams.rkf
            ret["elapsed"] = cast(float, self.readrkf("General", "ElapsedTime"))
            ret["system"] = cast(float, self.readrkf("General", "SysTime"))
            ret["cpu"] = cast(float, self.readrkf("General", "CPUTime"))
        except:
            # fall back to reading output, was needed for old AMS versions
            cpu = self.grep_output("Total cpu time:")
            system = self.grep_output("Total system time:")
            elapsed = self.grep_output("Total elapsed time:")
            ret["elapsed"] = float(elapsed[0].split()[-1])
            ret["system"] = float(system[0].split()[-1])
            ret["cpu"] = float(cpu[0].split()[-1])

        return ret

    def get_forcefield_params(
        self, engine: Optional[str] = None
    ) -> Tuple[List[float], List[str], Optional["ForceFieldPatch"]]:
        """Return force-field information from an engine RKF file.

        :param engine: engine RKF file identifier, defaults to the unique main engine
        :return: tuple containing the atomic charges in elementary-charge units, the force-field atom type for each
            atom, and the combined ``ForceFieldPatch``. The patch is ``None`` when the RKF file contains no
            patches.
        """
        from scm.plams.interfaces.adfsuite.forcefieldparams import forcefield_params_from_kf

        return self._process_engine_results(forcefield_params_from_kf, engine)

    def get_exit_condition_message(self) -> str:
        """Return the driver exit-condition message.

        :return: stored message, or an empty string if none is available
        """
        try:
            msg = self.readrkf("History", "ExitConditionMsg")
            assert isinstance(msg, str)
            return msg
        except:
            return ""

    def get_poissonratio(self, engine: Optional[str] = None) -> float:
        """Return the dimensionless Poisson ratio from ``AMSResults%PoissonRatio``.

        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: Poisson ratio as a dimensionless scalar
        """
        return self._process_engine_results(lambda x: x.read_real("AMSResults", "PoissonRatio"), engine)

    def get_youngmodulus(self, unit: str = "au", engine: Optional[str] = None) -> float:
        """Return Young's modulus from ``AMSResults%YoungModulus``.

        :param unit: pressure unit of the returned modulus, defaults to ``au`` (``hartree/bohr^3``); ``GPa`` is a
            commonly used alternative
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: Young's modulus as a scalar in ``unit``
        """
        return self._process_engine_results(
            lambda x: x.read_real("AMSResults", "YoungModulus"), engine
        ) * Units.conversion_ratio("au", unit)

    def get_shearmodulus(self, unit: str = "au", engine: Optional[str] = None) -> float:
        """Return the shear modulus from ``AMSResults%ShearModulus``.

        :param unit: pressure unit of the returned modulus, defaults to ``au`` (``hartree/bohr^3``); ``GPa`` is a
            commonly used alternative
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: shear modulus as a scalar in ``unit``
        """
        return self._process_engine_results(
            lambda x: x.read_real("AMSResults", "ShearModulus"), engine
        ) * Units.conversion_ratio("au", unit)

    def get_bulkmodulus(self, unit: str = "au", engine: Optional[str] = None) -> float:
        """Return the bulk modulus from ``AMSResults%BulkModulus``.

        :param unit: pressure unit of the returned modulus, defaults to ``au`` (``hartree/bohr^3``); ``GPa`` is a
            commonly used alternative
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: bulk modulus as a scalar in ``unit``
        """
        return self._process_engine_results(
            lambda x: x.read_real("AMSResults", "BulkModulus"), engine
        ) * Units.conversion_ratio("au", unit)

    @requires_optional_package("natsort")
    def get_pesscan_results(self, molecules: bool = True) -> Dict[str, Any]:
        """Return scan coordinates, energies, and convergence data for a PES scan.

        :param molecules: include ``Molecules`` with one structure per PES point
        :return: dictionary containing:

            - ``RaveledScanCoords`` and ``ScanCoords``: flat and grouped scan-coordinate names
            - ``nRaveledScanCoords`` and ``nScanCoords``: corresponding dimensionless counts
            - ``RaveledUnits`` and ``Units``: corresponding units, chosen from ``bohr``, ``bohr^2``,
              ``bohr^3``, and ``radian`` according to coordinate type
            - ``RaveledPESCoords``: values for every flat coordinate in its ``RaveledUnits`` entry
            - ``OrigScanCoords``: grouped names in their newline-separated RKF representation
            - ``nPESPoints``: dimensionless number of PES points
            - ``PES``: energies in hartree
            - ``Converged``: geometry-optimization convergence flags
            - ``HistoryIndices``: one-based indices into the ``History`` section
            - ``ConstrainedAtoms``: one-based indices of atoms occurring in scan coordinates
            - ``Properties``: engine properties in their RKF storage units; populated when
              ``CalcPropertiesAtPESPoints`` was enabled
            - ``Molecules``: structures with coordinates in angstrom, present only when ``molecules`` is ``True``
        """
        import re

        from natsort import natsorted

        def tolist(x: Any) -> List:
            if isinstance(x, list):
                return x
            else:
                return [x]

        nScanCoord = cast(int, self.readrkf("PESScan", "nScanCoord"))

        pes = tolist(self.readrkf("PESScan", "PES"))

        origscancoords: List[str] = tolist(self.get_history_property("ScanCoord", history_section="PESScan"))
        # one scan coordinate may have several variables
        scancoords = [x.split("\n") for x in origscancoords]
        pescoords = np.array(tolist(self.readrkf("PESScan", "PESCoords")))
        pescoords = pescoords.reshape(-1, sum(len(x) for x in scancoords))
        pescoords = np.transpose(pescoords)
        units: List[List[str]] = []
        for i in range(nScanCoord):
            units.append([])
            for j in range(len(scancoords[i])):
                if (
                    scancoords[i][j] in ["a", "b", "c"]
                    or "Dist" in scancoords[i][j]
                    or "Coordinate" in scancoords[i][j]
                ):
                    units[-1].append("bohr")
                elif "Volume" in scancoords[i][j]:
                    units[-1].append("bohr^3")
                elif "Area" in scancoords[i][j]:
                    units[-1].append("bohr^2")
                else:
                    units[-1].append("radian")

        converged = tolist(self.readrkf("PESScan", "GOConverged"))
        converged = [bool(x) for x in converged]
        historyindices = tolist(self.readrkf("PESScan", "HistoryIndices"))

        ret: Dict[str, Any] = {}

        raveled_scancoords = [x[i] for x in scancoords for i in range(len(x))]  # 1d list
        raveled_units = [x[i] for x in units for i in range(len(x))]  # 1d list

        constrained_atoms = []
        for sc in raveled_scancoords:
            constrained_atoms.extend(re.findall(r"\((\d+)\)", sc))

        try:
            constrained_atom_set = set(int(x) for x in constrained_atoms)
        except ValueError:
            constrained_atom_set = set()

        ret["RaveledScanCoords"] = raveled_scancoords
        ret["nRaveledScanCoords"] = len(raveled_scancoords)
        ret["ConstrainedAtoms"] = constrained_atom_set
        ret["ScanCoords"] = scancoords
        ret["nScanCoords"] = len(scancoords)
        ret["OrigScanCoords"] = origscancoords  # newline delimiter for joint scan coordinates
        ret["RaveledPESCoords"] = pescoords.tolist()
        ret["Units"] = units
        ret["RaveledUnits"] = raveled_units
        ret["Converged"] = converged
        ret["PES"] = pes
        ret["nPESPoints"] = len(pes)
        ret["HistoryIndices"] = historyindices

        if molecules:
            mols = []
            for ind in historyindices:
                mols.append(self.get_history_molecule(ind))
            ret["Molecules"] = mols

        ret["Properties"] = []
        for key in natsorted(title for title in self.rkfs.keys() if title.startswith("PESPoint")):
            amsresults = self.rkfs[key].read_section("AMSResults")  # type: ignore[index] # Python3.8 only - can be removed when support dropped
            ret["Properties"].append(amsresults)

        return ret

    def get_neb_results(self, molecules: bool = True, unit: str = "au") -> Dict[str, Any]:
        """Return results from a nudged elastic band calculation.

        :param molecules: include ``Molecules`` with all structures, including endpoints
        :param unit: energy unit for ``Energies``, barriers, and reaction energy, defaults to ``au``
        :return: dictionary containing dimensionless ``nImages`` and ``nIterations`` counts; the ``Climbing`` flag;
            ``HighestIndex``; ``Energies``, ``LeftBarrier``, ``RightBarrier``, and ``ReactionEnergy`` in ``unit``;
            one-based ``HistoryIndices``; and, when ``molecules`` is ``True``, ``Molecules`` with coordinates in
            angstrom. ``nImages`` excludes the endpoints, while ``Energies`` and ``Molecules`` include them.
        """

        def tolist(x: Any) -> List:
            if isinstance(x, list):
                return x
            else:
                return [x]

        ret: Dict[str, Any] = {}
        conversion_ratio = Units.conversion_ratio("au", unit)
        ret["nImages"] = self.readrkf("NEB", "nebImages")
        ret["nIterations"] = self.readrkf("NEB", "nebIterations")
        ret["Climbing"] = bool(self.readrkf("NEB", "climbing"))
        ret["HighestIndex"] = self.readrkf("NEB", "highestIndex")
        ret["LeftBarrier"] = cast(float, self.readrkf("NEB", "LeftBarrier")) * conversion_ratio
        ret["RightBarrier"] = cast(float, self.readrkf("NEB", "RightBarrier")) * conversion_ratio
        ret["ReactionEnergy"] = cast(float, self.readrkf("NEB", "ReactionEnergy")) * conversion_ratio
        history_dim = tolist(self.readrkf("NEB", "historyIndex@dim"))  # nimages, randombign
        history_dim.reverse()  # randombign, nimages
        history_indices_matrix = np.array(
            tolist(self.readrkf("NEB", "historyIndex"))
        )  # this matrix is padded with -1 values
        history_indices_matrix = history_indices_matrix.reshape(history_dim)
        history_indices = np.max(history_indices_matrix, axis=0, keepdims=False).tolist()
        if any(x == -1 for x in history_indices):
            raise ValueError("Found -1 in the 'converged' part of historyIndex. This should not happen!")
        ret["HistoryIndices"] = history_indices
        ret["Energies"] = [
            cast(float, self.get_property_at_step(ind, "Energy")) * conversion_ratio for ind in ret["HistoryIndices"]
        ]
        if molecules:
            ret["Molecules"] = [self.get_history_molecule(ind) for ind in cast(Sequence, ret["HistoryIndices"])]

        return ret

    def get_irc_results(self, molecules: bool = True, unit: str = "au") -> Dict[str, Any]:
        """Return the converged points from an intrinsic reaction coordinate calculation.

        ``IRCGradMax`` and ``IRCGradRms`` are returned in their RKF storage convention. AMS uses
        mass-weighted atomic units during IRC steps and ordinary atomic units during final minimization.

        :param molecules: include ``Molecules`` with the converged structures
        :param unit: energy unit for energies and barriers, defaults to ``au``
        :return: dictionary containing ``Energies``, ``RelativeEnergies``, ``LeftBarrier``, and ``RightBarrier`` in
            ``unit``; ``PathLength`` and ``ArcLength`` in angstrom; dimensionless ``IRCIteration`` and one-based
            ``HistoryIndices``; ``IRCDirection`` values; ``IRCGradMax`` and ``IRCGradRms`` in the convention above;
            and, when ``molecules`` is ``True``, ``Molecules`` with coordinates in angstrom
        """
        from itertools import compress

        def tolist(x: Any) -> List:
            if isinstance(x, list):
                return x
            else:
                return [x]

        sec = "History"
        items = [
            "IRCDirection",
            "Energy",
            "PathLength",
            "IRCIteration",
            "IRCGradMax",
            "IRCGradRms",
            "ArcLength",
            "HistoryIndices",
            "Molecules",
        ]
        d: Dict[str, List] = {}
        forw = {}
        back = {}
        reformed = {}
        forw_mask = None
        back_mask = None
        converged = self.get_history_property("Converged", history_section=sec) or []
        converged = tolist(converged)
        converged_mask = [x != 0 for x in converged]
        history_indices = [i for i, x in enumerate(converged_mask, 1) if x]  # raw, not rearranged
        for k in items:
            # first half is forward direction, second half is backward direction
            # second half should be reversed, path length made negative.
            if k == "Molecules":
                if not molecules:
                    continue
                d[k] = [self.get_history_molecule(ind) for ind in history_indices]  # rearrangement happens later
            elif k == "HistoryIndices":
                d[k] = history_indices  # rearrangement happens later
            else:
                d[k] = self.get_history_property(k, history_section=sec) or []
                d[k] = tolist(d[k])
                d[k] = list(compress(d[k], converged_mask))
            if k == "IRCDirection":
                forw_mask = [x == 1 for x in d[k]]
                back_mask = [x != 1 for x in d[k]]
                d[k] = ["Forward" if x == 1 else "Backward" if x == 2 else x for x in d[k]]

            forw[k] = list(compress(d[k], forw_mask))  # type: ignore[arg-type]
            back[k] = list(compress(d[k], back_mask))  # type: ignore[arg-type]
            back[k].reverse()
            if k == "PathLength":
                # print backwards direction as negative numbers
                back[k] = [-x for x in back[k]]
            reformed[k] = back[k] + forw[k]

        conversion_ratio = Units.convert(1.0, "au", unit)
        reformed["Energies"] = [x * conversion_ratio for x in reformed["Energy"]]
        del reformed["Energy"]
        max_energy = max(reformed["Energies"])
        if "Forward" in reformed["IRCDirection"]:
            reformed["LeftBarrier"] = max_energy - reformed["Energies"][0]
        if "Backward" in reformed["IRCDirection"]:
            reformed["RightBarrier"] = max_energy - reformed["Energies"][-1]
        reformed["RelativeEnergies"] = [x - reformed["Energies"][0] for x in reformed["Energies"]]

        return reformed

    def get_time_step(self, history_section: Literal["BinLog", "MDHistory"] = "MDHistory") -> float:
        """Return the time between adjacent saved MD frames in femtoseconds.

        This is the MD integration time step multiplied by the trajectory
        sampling frequency.

        :param history_section: history section containing the saved times, defaults to ``MDHistory``
        :return: time between the first two saved frames in femtoseconds
        :raises PlamsError: if the time cannot be determined from the first two saved frames
        """
        time1 = cast(Optional[float], self.get_property_at_step(1, "Time", history_section=history_section))
        time2 = cast(Optional[float], self.get_property_at_step(2, "Time", history_section=history_section))

        if time1 is None or time2 is None:
            raise PlamsError(f"Cannot determine time step from 'Time' variable of history section '{history_section}'")

        time_step = time2 - time1
        return time_step

    def _get_integer_start_end_every_max(
        self,
        start_fs: float,
        end_fs: Optional[float],
        every_fs: Optional[float],
        max_dt_fs: Optional[float],
        time_step: Optional[float] = None,
    ) -> Tuple[int, Optional[int], int, Optional[int]]:
        """Convert time-based trajectory bounds to integer indices.

        :param start_fs: start time in femtoseconds
        :param end_fs: end time in femtoseconds
        :param every_fs: sampling interval in femtoseconds
        :param max_dt_fs: maximum correlation time in femtoseconds
        :param time_step: interval between stored frames, determined automatically if omitted
        :return: start index, optional end index, sampling stride, and optional maximum correlation lag
        """
        if time_step is None:
            time_step = self.get_time_step()
        start_step = round(start_fs / time_step)
        end_step = round(end_fs / time_step) if end_fs is not None else None
        every = round(every_fs / time_step) if every_fs is not None else 1
        max_dt = round((max_dt_fs / time_step) / every) if max_dt_fs is not None else None
        return start_step, end_step, every, max_dt

    def get_velocity_acf(
        self,
        start_fs: float = 0,
        end_fs: Optional[float] = None,
        every_fs: Optional[float] = None,
        max_dt_fs: Optional[float] = None,
        atom_indices: Optional[Sequence[int]] = None,
        x: bool = True,
        y: bool = True,
        z: bool = True,
        normalize: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the velocity autocorrelation function for a fixed-composition trajectory.

        :param start_fs: start time in femtoseconds
        :param end_fs: end time in femtoseconds, or ``None`` for the end of the trajectory
        :param every_fs: sampling interval in femtoseconds, or ``None`` for every frame
        :param max_dt_fs: maximum correlation time in femtoseconds
        :param atom_indices: one-based atom indices, or ``None`` for all atoms
        :param x: include the x velocity component
        :param y: include the y velocity component
        :param z: include the z velocity component
        :param normalize: normalize the autocorrelation function to one at time zero
        :return: correlation times in femtoseconds and the autocorrelation function. The latter is
            dimensionless when normalized and otherwise in ``angstrom^2/fs^2``.
        """
        from scm.plams.trajectories.analysis import autocorrelation

        nEntries = cast(int, self.readrkf("MDHistory", "nEntries"))

        time_step = self.get_time_step()

        start_step, end_step, every, max_dt = self._get_integer_start_end_every_max(
            start_fs, end_fs, every_fs, max_dt_fs
        )

        data = np.array(self.get_history_property("Velocities", history_section="MDHistory")).reshape(nEntries, -1, 3)[
            start_step:end_step:every
        ]
        if atom_indices is not None:
            zero_based_atom_indices = [x - 1 for x in atom_indices]
            data = data[:, zero_based_atom_indices, :]

        n_dimensions = 3
        if not x or not y or not z:
            components = []
            if x:
                components += [0]
            if y:
                components += [1]
            if z:
                components += [2]

            data = data[:, :, components]
            n_dimensions = len(components)

        data *= Units.convert(1.0, "bohr", "angstrom")  # convert from bohr/fs to ang/fs

        vacf = (
            autocorrelation(data, max_dt=max_dt, normalize=normalize) * n_dimensions
        )  # multiply by n_dimensions to undo the averaging per component and instead average per atom

        times = np.arange(len(vacf)) * time_step * every

        return times, vacf

    def get_dipole_history(self, dipole_unit: str = "e*bohr") -> np.ndarray:
        """Return the molecular dipole moment for each saved frame.

        :param dipole_unit: dipole-moment unit of the returned values, defaults to ``e*bohr``
        :return: array of shape ``(n_frames, 3)`` in ``dipole_unit``
        """
        dipole_x = self.get_history_property(history_section="BinLog", varname="DipoleMoment_x")
        dipole_y = self.get_history_property(history_section="BinLog", varname="DipoleMoment_y")
        dipole_z = self.get_history_property(history_section="BinLog", varname="DipoleMoment_z")
        data = np.column_stack((dipole_x, dipole_y, dipole_z))  # type: ignore[arg-type]
        data *= Units.convert(1.0, "e*bohr", dipole_unit)
        return data

    def get_dipole_derivatives_acf(
        self,
        start_fs: float = 0,
        end_fs: Optional[float] = None,
        every_fs: Optional[float] = None,
        max_dt_fs: Optional[float] = None,
        x: bool = True,
        y: bool = True,
        z: bool = True,
        normalize: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate the autocorrelation function of the dipole derivative.

        :param start_fs: start time in femtoseconds, defaults to ``0``
        :param end_fs: end time in femtoseconds, defaults to the end of the trajectory
        :param every_fs: sampling interval in femtoseconds, defaults to every frame
        :param max_dt_fs: maximum correlation time in femtoseconds
        :param x: include the x component, defaults to ``True``
        :param y: include the y component, defaults to ``True``
        :param z: include the z component, defaults to ``True``
        :param normalize: normalize the autocorrelation function to one at time zero, defaults to ``False``
        :return: tuple containing the correlation times in femtoseconds and the dipole-derivative autocorrelation
            function. A normalized autocorrelation function is dimensionless; otherwise its unit is
            ``(e*bohr/fs)^2``.
        """
        from scm.plams.trajectories.analysis import autocorrelation

        # nEntries = self.readrkf('BinLog', 'nEntries')
        time_step = self.get_time_step(history_section="BinLog")

        data = self.get_dipole_history()

        start_step, end_step, every, max_dt = self._get_integer_start_end_every_max(
            start_fs, end_fs, every_fs, max_dt_fs, time_step
        )
        data = data[start_step:end_step:every, :]

        n_dimensions = 3
        if not x or not y or not z:
            components = []
            if x:
                components += [0]
            if y:
                components += [1]
            if z:
                components += [2]

            data = data[:, components]
            n_dimensions = len(components)

        data_deriv = np.diff(data, axis=0) / time_step
        if max_dt is None:  # default value is 5000 or num_timesteps // 2
            num_timesteps = data_deriv.shape[0]
            max_dt = num_timesteps // 2
            if max_dt > 5000:
                max_dt = 5000
        dipole_deriv_acf = (
            autocorrelation(data_deriv, max_dt=max_dt, normalize=normalize) * n_dimensions
        )  # multiply by n_dimensions to undo the averaging per component and instead average per atom

        times = np.arange(len(dipole_deriv_acf)) * time_step * every

        return times, dipole_deriv_acf

    @requires_optional_package("scipy")
    def get_diffusion_coefficient_from_velocity_acf(
        self, times: Optional[np.ndarray] = None, acf: Optional[np.ndarray] = None, n_dimensions: int = 3
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Integrate a velocity autocorrelation function to obtain a diffusion coefficient.

        If ``times`` or ``acf`` is None, then a default velocity autocorrelation function will be calculated.

        :param times: correlation times in femtoseconds
        :param acf: velocity autocorrelation in ``angstrom^2/fs^2``
        :param n_dimensions: number of spatial dimensions included in the autocorrelation
        :return: times in femtoseconds and cumulative diffusion coefficients in ``m^2/s``

        """

        from scipy.integrate import cumtrapz

        if times is None or acf is None:
            times, acf = self.get_velocity_acf(normalize=False)

        diffusion_coefficient = (1.0 / n_dimensions) * cumtrapz(acf, times)  # ang^2/fs
        diffusion_coefficient *= 1e-20 / 1e-15
        new_times = times[:-1]

        return np.array(new_times), diffusion_coefficient

    def get_power_spectrum(
        self,
        times: Optional[np.ndarray] = None,
        acf: Optional[np.ndarray] = None,
        max_dt_fs: Optional[float] = None,
        max_freq: Optional[float] = None,
        number_of_points: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate a power spectrum from the velocity autocorrelation function.

        If ``times`` or ``acf`` is None, then a default velocity autocorrelation function will be calculated.

        :param times: correlation times in femtoseconds, as returned by :meth:`get_velocity_acf`
        :param acf: velocity autocorrelation function
        :param max_dt_fs: maximum correlation time in femtoseconds used for the default autocorrelation function
        :param max_freq: maximum returned frequency in ``cm^-1``, defaults to ``5000``
        :param number_of_points: number of spectrum points, defaults to approximately ``1 cm^-1`` spacing
        :return: frequencies in ``cm^-1`` and corresponding intensities. The internally generated normalized
            autocorrelation function gives arbitrary intensity units; for a supplied autocorrelation function, the
            numerical scale follows that function.
        """
        from scm.plams.trajectories.analysis import power_spectrum

        if times is None or acf is None:
            times, acf = self.get_velocity_acf(max_dt_fs=max_dt_fs, normalize=True)

        return power_spectrum(times, acf, max_freq=max_freq, number_of_points=number_of_points)

    def get_ir_spectrum_md(
        self,
        times: Optional[np.ndarray] = None,
        acf: Optional[np.ndarray] = None,
        max_dt_fs: Optional[float] = None,
        max_freq: Optional[float] = 5000,
        number_of_points: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate an IR spectrum from a dipole-derivative autocorrelation function.

        If ``times`` or ``acf`` is omitted, a normalized dipole-derivative autocorrelation function is calculated.

        :param times: correlation times in femtoseconds, such as those returned by
            :meth:`get_dipole_derivatives_acf`
        :param acf: dipole-derivative autocorrelation function
        :param max_dt_fs: maximum correlation time in femtoseconds used when calculating the default autocorrelation
            function
        :param max_freq: maximum returned frequency in ``cm^-1``, defaults to ``5000``
        :param number_of_points: number of spectrum points, defaults to a spacing of approximately ``1 cm^-1``
        :return: tuple containing the frequencies in ``cm^-1`` and the IR intensities. The internally generated
            normalized autocorrelation function gives arbitrary intensity units. When an autocorrelation function is
            supplied, the numerical scale follows that function.
        """
        from scm.plams.trajectories.analysis import power_spectrum

        if times is None or acf is None:
            times, acf = self.get_dipole_derivatives_acf(max_dt_fs=max_dt_fs, normalize=True)

        return power_spectrum(times, acf, max_freq=max_freq, number_of_points=number_of_points)

    @staticmethod
    @requires_optional_package("scipy")
    def _get_green_kubo_viscosity(
        pressuretensor: np.ndarray,
        time_step: float,
        max_dt: float,
        volume: float,
        temperature: float,
        xy: bool = True,
        yz: bool = True,
        xz: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        from scipy.integrate import cumtrapz

        from scm.plams.tools.units import Units
        from scm.plams.trajectories.analysis import autocorrelation

        data = np.array(pressuretensor)

        components = []
        if yz:
            components += [3]
        if xz:
            components += [4]
        if xy:
            components += [5]

        data = data[:, components]
        acf = autocorrelation(data, max_dt=max_dt)
        times = np.arange(len(acf)) * time_step

        integrated_acf = cumtrapz(acf, times, initial=0)
        integrated_times = np.arange(len(integrated_acf)) * time_step

        k_B = Units.constants["k_B"]

        au2Pa = Units.convert(1.0, "hartree/bohr^3", "Pa")

        viscosity = (volume * 1e-30) / (k_B * temperature) * integrated_acf * au2Pa**2 * 1e-15 * 1e3

        return integrated_times, viscosity

    @requires_optional_package("scipy")
    def get_green_kubo_viscosity(
        self,
        start_fs: float = 0,
        end_fs: Optional[float] = None,
        every_fs: Optional[float] = None,
        max_dt_fs: Optional[float] = None,
        xy: bool = True,
        yz: bool = True,
        xz: bool = True,
        pressuretensor: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate viscosity from the off-diagonal pressure-tensor autocorrelation function.

        :param start_fs: start time in femtoseconds
        :param end_fs: end time in femtoseconds, or ``None`` for the end of the trajectory
        :param every_fs: sampling interval in femtoseconds, or ``None`` for every frame
        :param max_dt_fs: maximum correlation time in femtoseconds
        :param xy: include the xy pressure-tensor component
        :param yz: include the yz pressure-tensor component
        :param xz: include the xz pressure-tensor component
        :param pressuretensor: pressure tensors in ``hartree/bohr^3`` with shape ``(n_frames, 6)``, or
            ``None`` to read them from ``ams.rkf``
        :return: integration times in femtoseconds and cumulative viscosity in ``mPa*s``

        """
        from scipy.integrate import cumtrapz

        from scm.plams.tools.units import Units
        from scm.plams.trajectories.analysis import autocorrelation

        time_step: float = self.get_time_step()
        start_step, end_step, every, max_dt = self._get_integer_start_end_every_max(
            start_fs, end_fs, every_fs, max_dt_fs
        )
        if pressuretensor is None:
            rkf_pressuretensor = self.get_history_property("PressureTensor", "MDHistory") or []
            data = np.array(
                [x for x in rkf_pressuretensor if x is not None]
            )  # None might appear in currently running trajectories if the job was loaded with load_external
        else:
            data = pressuretensor
        data = np.array(data)[start_step:end_step:every]

        components = []
        if yz:
            components += [3]
        if xz:
            components += [4]
        if xy:
            components += [5]

        data = data[:, components]
        acf = autocorrelation(data, max_dt=max_dt)
        times = np.arange(len(acf)) * time_step * every

        integrated_acf = cumtrapz(acf, times)
        integrated_times = np.arange(len(integrated_acf)) * time_step * every

        V = self.get_main_molecule().unit_cell_volume()

        try:
            T = np.array(self.get_history_property("Temperature", "MDHistory"))[start_step:end_step:every]
            mean_T = np.mean(T)
        except KeyError:  # might be triggered for currently running trajectories, then just use the first temperature
            mean_T = cast(float, self.get_property_at_step(1, "Temperature", "MDHistory"))

        k_B = Units.constants["k_B"]

        au2Pa = Units.convert(1.0, "hartree/bohr^3", "Pa")

        viscosity = (V * 1e-30) / (k_B * mean_T) * integrated_acf * au2Pa**2 * 1e-15 * 1e3

        return integrated_times, viscosity

    def get_density_along_axis(
        self,
        axis: str = "z",
        density_type: str = "mass",
        start_fs: float = 0,
        end_fs: Optional[float] = None,
        every_fs: Optional[float] = None,
        bin_width: float = 0.1,
        atom_indices: Optional[Sequence[int]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the atomic density along a Cartesian coordinate axis.

        This only works if the axis is perpendicular to the other two axes. The
        system must be 3D-periodic and the number of atoms cannot change during
        the trajectory.

        :param axis: Cartesian axis, one of ``x``, ``y``, or ``z``
        :param density_type: ``mass`` for ``g/cm^3``, ``number`` for ``angstrom^-3``, or
            ``histogram`` for the mean number of atoms per frame and bin
        :param start_fs: start time in femtoseconds
        :param end_fs: end time in femtoseconds, or ``None`` for the end of the trajectory
        :param every_fs: sampling interval in femtoseconds, or ``None`` for every frame
        :param bin_width: target bin width in angstrom
        :param atom_indices: one-based atom indices, or ``None`` for all atoms
        :return: bin-center coordinates in angstrom and the selected density

        """
        assert axis in ["x", "y", "z"], f"Unknown axis: {axis}. Must be 'x', 'y', or 'z'"
        assert density_type in [
            "mass",
            "number",
            "histogram",
        ], f"Unknown density_type: {density_type}. Must be 'mass' or 'number'"

        main_mol = self.get_main_molecule()

        histogram = density_type == "histogram"

        bohr2ang = Units.convert(1.0, "bohr", "angstrom")

        start_step, end_step, every, _ = self._get_integer_start_end_every_max(start_fs, end_fs, every_fs, None)
        nEntries = cast(int, self.readrkf("History", "nEntries"))
        history_coords = np.array(self.get_history_property("Coords")).reshape(nEntries, -1, 3)
        coords = history_coords[start_step:end_step:every]
        nEntries = len(coords)

        axis2index = {"x": 0, "y": 1, "z": 2}
        index = axis2index[axis]

        if not histogram:
            other_indices = [0, 1, 2]
            other_indices.remove(index)

            assert (
                len(main_mol.lattice) == 3
            ), f"get_density_along_axis with density_type='mass' or 'number' can only be called for 3d-periodic systems. Current periodicity: {len(main_mol.lattice)}. Use density_type='histogram' instead."
            assert (
                np.abs(np.dot(main_mol.lattice[other_indices[0]], main_mol.lattice[index])) < 1e-6
            ), f"Axis {axis} must be perpendicular to the other two axes"
            assert (
                np.abs(np.dot(main_mol.lattice[other_indices[1]], main_mol.lattice[index])) < 1e-6
            ), f"Axis {axis} must be perpendicular to the other two axes"
            assert (
                np.abs(main_mol.lattice[index][other_indices[0]]) < 1e-6
            ), f"Density along {axis} requires that lattice vector {index+1} has only 1 non-zero component (along {axis}). The vector is {main_mol.lattice[index]}"
            assert (
                np.abs(main_mol.lattice[index][other_indices[1]]) < 1e-6
            ), f"Density along {axis} requires that lattice vector {index+1} has only 1 non-zero component (along {axis}). The vector is {main_mol.lattice[index]}"

            lattice_vectors = np.array(self.get_history_property("LatticeVectors")).reshape(-1, 3, 3)
            lattice_vectors = lattice_vectors[start_step:end_step:every] * bohr2ang

            vec1s = lattice_vectors[:, other_indices[0], :]
            vec2s = lattice_vectors[:, other_indices[1], :]
            cross_products = np.cross(vec1s, vec2s)
            areas = np.linalg.norm(cross_products, axis=1)
            mean_area = np.mean(areas)

        coords = coords[:, :, index]

        if atom_indices is not None:
            zero_based_atom_indices = [x - 1 for x in atom_indices]
            coords = coords[:, zero_based_atom_indices]

        coords *= bohr2ang

        avogadro = Units.constants["NA"]

        min_c = np.min(coords)
        max_c = np.max(coords)
        num_bins = round((max_c - min_c) / bin_width + 1)
        bins = np.linspace(min_c, max_c, num_bins, endpoint=True)

        if density_type == "mass":
            masses = np.array(self.readrkf("Molecule", "AtomMasses"))
            if atom_indices is not None:
                masses = masses[zero_based_atom_indices]
            masses_broadcasted = np.broadcast_to(masses, coords.shape)
            hist, bin_edges = np.histogram(coords, weights=masses_broadcasted, bins=bins)
            volume_slice_ang3 = mean_area * (bin_edges[1] - bin_edges[0])
            volume_slice_cm3 = volume_slice_ang3 * 1e-24
            density = (1.0 / nEntries) * (hist / avogadro) / volume_slice_cm3
        elif density_type == "number":
            hist, bin_edges = np.histogram(coords, bins=bins)
            volume_slice_ang3 = mean_area * (bin_edges[1] - bin_edges[0])
            density = (1.0 / nEntries) * hist / volume_slice_ang3
        elif density_type == "histogram":
            hist, bin_edges = np.histogram(coords, bins=bins)
            density = (1.0 / nEntries) * hist

        z = (bin_edges[:-1] + bin_edges[1:]) / 2.0
        return z, density

    def get_work_function_results(
        self, energy_unit: str = "hartree", dist_unit: str = "bohr"
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float, Tuple[float, float], Tuple[float, float]]:
        """Return the results of a work-function calculation.

        The two vacuum potentials and work functions correspond to the left and right sides of the
        slab, respectively. For Quantum ESPRESSO calculations using a dipole correction, the
        potential is evaluated immediately before the correction position.

        :param energy_unit: unit for potentials, Fermi energy, and work functions, defaults to ``hartree``
        :param dist_unit: unit for the coordinate perpendicular to the surface, defaults to ``bohr``
        :return: coordinate in ``dist_unit``; planar-average and macroscopic-average potentials, Fermi energy, and
            bulk potential in ``energy_unit``; left and right vacuum potentials in ``energy_unit``; and left and
            right work functions in ``energy_unit``
        """

        to_eunit = Units.conversion_ratio("au", energy_unit)
        to_dunit = Units.conversion_ratio("au", dist_unit)

        coordinate = np.array(self.readrkf("WorkFunction", "coordinate", file="engine")) * to_dunit
        planarAverage = np.array(self.readrkf("WorkFunction", "planarAverage", file="engine")) * to_eunit
        macroscopicAverage = np.array(self.readrkf("WorkFunction", "macroscopicAverage", file="engine")) * to_eunit
        Efermi = cast(float, self.readrkf("WorkFunction", "fermiEnergy", file="engine")) * to_eunit
        Vbulk = cast(float, self.readrkf("WorkFunction", "minMacroscopicAverPotential", file="engine")) * to_eunit
        leftVvacuum = cast(float, self.readrkf("WorkFunction", "leftVacuumPotential", file="engine")) * to_eunit
        rightVvacuum = cast(float, self.readrkf("WorkFunction", "rightVacuumPotential", file="engine")) * to_eunit
        leftWF = cast(float, self.readrkf("WorkFunction", "leftWorkFunction", file="engine")) * to_eunit
        rightWF = cast(float, self.readrkf("WorkFunction", "rightWorkFunction", file="engine")) * to_eunit

        return (
            coordinate,
            planarAverage,
            macroscopicAverage,
            Efermi,
            Vbulk,
            (leftVvacuum, rightVvacuum),
            (leftWF, rightWF),
        )

    def recreate_molecule(self) -> Union[None, Molecule, Dict[str, Molecule]]:
        """Recreate the input molecule(s) for the corresponding job based on files present in the job folder.

        This method is used by |load_external|. If ``ams.rkf`` is present in the job folder,
        extract data from the ``InputMolecule`` and ``InputMolecule(*)`` sections.

        :return: the input molecule or named input molecules, with coordinates and lattice vectors in angstrom; or
            ``None`` if they cannot be reconstructed
        """
        if "ams" in self.rkfs:
            mols = self.get_input_molecules()
            if len(mols) == 0:
                return None
            elif len(mols) == 1 and "" in mols:
                return mols[""]
            else:
                return mols
        return None

    def recreate_settings(self) -> Optional[Settings]:
        """Recreate the input |Settings| instance for the corresponding job based on files present in the job folder. This method is used by |load_external|.

        If ``ams.rkf`` is present in the job folder, extract user input and parse it back to a |Settings| instance using ``scm.base`` module. Remove the ``system`` branch from that instance.

        :return: reconstructed settings, or ``None`` if they cannot be reconstructed
        """
        if "ams" in self.rkfs:
            user_input = cast(str, self.readrkf("General", "user input"))
            try:
                from scm.plams.interfaces.adfsuite.inputparser import input_to_settings

                inp = input_to_settings(user_input)
            except:
                log("Failed to recreate input settings from {}".format(self.rkfs["ams"].path))
                return None
            s = Settings()
            s.input = inp
            if "system" in s.input.ams:
                del s.input.ams.system
            s.soft_update(get_config().job)
            return s
        return None

    def ok(self) -> bool:
        """Check if the execution of the associated :attr:`job` was successful or not.
        See :meth:`Job.ok<scm.plams.core.basejob.Job.ok>` for more information.

        :return: whether the associated job completed successfully
        """
        return self.job.ok()

    def get_errormsg(self) -> Optional[str]:
        """Return an error message for the associated job.

        See :meth:`Job.get_errormsg<scm.plams.core.basejob.Job.errormsg>` for more information.

        :return: error message, or ``None`` for a successful job
        """
        return self.job.get_errormsg()

    @property
    def name(self) -> str:
        """Return the name of the associated job.

        :return: :attr:`job.name`
        """
        return self.job.name

    class EnergyLandscape:
        """Stationary points, fragments, and connectivity from a PES exploration.

        Integer indexing uses one-based state identifiers; iteration follows landscape order.

        :param results: AMS results containing an ``EnergyLandscape`` section, or ``None`` to create
            an empty landscape
        """

        class State:
            """A local minimum or transition state in an energy landscape.

            ``energy`` is in hartree, and ``molecule`` coordinates and lattice vectors are in angstrom.
            ``id``, ``display_id``, ``reactantsID``, and ``productsID`` are one-based identifiers.

            :param landscape: containing energy landscape
            :param engfile: engine RKF identifier for this state
            :param energy: state energy in hartree
            :param mol: state geometry
            :param count: number of times the state was found
            :param isTS: whether the state is a transition state
            :param reactantsID: one-based reactant-state identifier
            :param productsID: one-based product-state identifier
            :param prefactorsFromReactant: forward rate prefactor in ``s^-1``
            :param prefactorsFromProduct: backward rate prefactor in ``s^-1``
            :param originalID: identifier displayed for a state copied to a sub-landscape
            """

            def __init__(
                self,
                landscape: "AMSResults.EnergyLandscape",
                engfile: str,
                energy: float,
                mol: Molecule,
                count: int,
                isTS: bool,
                reactantsID: Optional[int] = None,
                productsID: Optional[int] = None,
                prefactorsFromReactant: Optional[float] = None,
                prefactorsFromProduct: Optional[float] = None,
            ):
                self._landscape = landscape
                self.engfile = engfile
                self.energy = energy
                self.molecule = mol
                self.count = count
                self.isTS = isTS
                self.reactantsID = reactantsID
                self.productsID = productsID
                self.prefactorsFromReactant = prefactorsFromReactant
                self.prefactorsFromProduct = prefactorsFromProduct

            @property
            def id(self) -> int:
                """Return the one-based state identifier."""
                return self._landscape._states.index(self) + 1

            @property
            def reactants(self) -> Optional["AMSResults.EnergyLandscape.State"]:
                """Return the reactant state connected to this transition state."""
                return self._landscape._states[self.reactantsID - 1] if self.reactantsID is not None else None

            @property
            def products(self) -> Optional["AMSResults.EnergyLandscape.State"]:
                """Return the product state connected to this transition state."""
                return self._landscape._states[self.productsID - 1] if self.productsID is not None else None

            def __str__(self) -> str:
                if self.isTS:
                    lines = [
                        f"State {self.id}: {self.molecule.get_formula(False)} transition state @ {self.energy:.8f} Hartree (found {self.count} times"
                        + (f", results on {self.engfile})" if self.engfile is not None else ")")
                    ]
                    if self.reactantsID is not None:
                        lines += [f"  +- Reactants: {self.reactants}"]
                    if self.productsID is not None:
                        lines += [f"     Products:  {self.products}"]

                    if self.reactantsID is not None and self.productsID is not None:
                        eV = Units.convert(1.0, "eV", "hartree")
                        lines += [
                            f"     Prefactors: {self.prefactorsFromReactant:.3E}:{self.prefactorsFromProduct:.3E} s^-1"
                        ]
                        dEforward = (self.energy - self._landscape._states[self.reactantsID - 1].energy) / eV
                        dEbackward = (self.energy - self._landscape._states[self.productsID - 1].energy) / eV
                        lines += [f"     Barriers: {dEforward:.3f}:{dEbackward:.3f} eV"]
                    elif self.productsID is not None:
                        lines += [f"     Prefactors: {self.prefactorsFromReactant:.3E}:?"]
                else:
                    lines = [
                        f"State {self.id}: {self.molecule.get_formula(False)} local minimum @ {self.energy:.8f} Hartree (found {self.count} times"
                        + (f", results on {self.engfile})" if self.engfile is not None else ")")
                    ]
                return "\n".join(lines)

        class Fragment:
            """An isolated fragment represented in an energy landscape.

            ``energy`` is in hartree, and ``molecule`` coordinates and lattice vectors are in angstrom.
            The ``id`` property is a one-based identifier.

            :param landscape: containing energy landscape
            :param engfile: engine RKF identifier for this fragment
            :param energy: fragment energy in hartree
            :param mol: fragment geometry
            """

            def __init__(self, landscape: "AMSResults.EnergyLandscape", engfile: str, energy: float, mol: Molecule):
                self._landscape = landscape
                self.engfile = engfile
                self.energy = energy
                self.molecule = mol

            @property
            def id(self) -> int:
                """Return the one-based fragment identifier."""
                return self._landscape._fragments.index(self) + 1

            def __str__(self) -> str:
                lines = [
                    f"Fragment {self.id}: {self.molecule.get_formula(False)} local minimum @ {self.energy:.8f} Hartree (results on {self.engfile})"
                ]
                return "\n".join(lines)

        class FragmentedState:
            """A state composed of multiple isolated fragments.

            ``energy`` is in hartree, and the ``id`` property is a one-based identifier. Composition and
            connection values are stored internally as zero-based indices. Adsorption prefactors are in
            ``(bar*s)^-1`` and desorption prefactors are in ``s^-1``.

            :param landscape: containing energy landscape
            :param energy: fragmented-state energy in hartree
            :param composition: zero-based fragment indices
            :param connections: zero-based connected-state indices
            :param adsorptionPrefactors: adsorption prefactors in ``(bar*s)^-1`` for the connections
            :param desorptionPrefactors: desorption prefactors in ``s^-1`` for the connections
            """

            def __init__(
                self,
                landscape: "AMSResults.EnergyLandscape",
                energy: float,
                composition: Sequence[int],
                connections: Optional[Sequence[int]] = None,
                adsorptionPrefactors: Optional[Sequence[float]] = None,
                desorptionPrefactors: Optional[Sequence[float]] = None,
            ):
                self._landscape = landscape
                self.energy = energy
                self.composition = composition
                self.connections = connections
                self.adsorptionPrefactors = adsorptionPrefactors
                self.desorptionPrefactors = desorptionPrefactors

            @property
            def id(self) -> int:
                """Return the one-based fragmented-state identifier."""
                return self._landscape._fstates.index(self) + 1

            @property
            def fragments(self) -> List["AMSResults.EnergyLandscape.Fragment"]:
                """Return the fragments that compose this state."""
                return [self._landscape._fragments[i] for i in self.composition]

            def __str__(self) -> str:
                formula = ""
                for i, id in enumerate(self.composition):
                    formula += self._landscape._fragments[id].molecule.get_formula(False)
                    if i != len(self.composition) - 1:
                        formula += "+"
                lines = [
                    f"FragmentedState {self.id}: {formula} local minimum @ {self.energy:.8f} Hartree (fragments {[i+1 for i in self.composition]})"
                ]
                if self.connections is not None:
                    for i, iState in enumerate(self.connections):
                        lines += [f"  +- {self._landscape._states[iState]}"]

                        if self.adsorptionPrefactors is not None and self.desorptionPrefactors is not None:
                            if i == len(self.connections) - 1:
                                lines += [
                                    f"     Prefactors: {self.adsorptionPrefactors[i]:.3E}:{self.desorptionPrefactors[i]:.3E}"
                                ]
                            else:
                                lines += [
                                    f"  |  Prefactors: {self.adsorptionPrefactors[i]:.3E}:{self.desorptionPrefactors[i]:.3E}"
                                ]
                return "\n".join(lines)

        def __init__(self, results: "AMSResults"):
            self._states: List["AMSResults.EnergyLandscape.State"] = []
            self._fragments: List["AMSResults.EnergyLandscape.Fragment"] = []
            self._fstates: List["AMSResults.EnergyLandscape.FragmentedState"] = []
            if results is None:
                return

            if "EnergyLandscape" not in results.get_rkf_skeleton():
                return

            sec = results.read_rkf_section("EnergyLandscape")

            # If there is only 1 state in the 'EnergyLandscape' section, some variables that are normally lists are insted be build-in types (e.g. a 'float' instead of a 'list of floats').
            # For convenience here we make sure that the following variables are always 'lists':
            for var in ["energies", "counts", "isTS", "reactants", "products"]:
                if not isinstance(sec[var], list):
                    sec[var] = [sec[var]]  # type: ignore[list-item]

            nStates = cast(int, sec["nStates"])

            for iState in range(nStates):
                energy = cast(List[float], sec["energies"])[iState]
                resfile = os.path.splitext(cast(str, sec["fileNames"]).split("\0")[iState])[0]
                mol = results.get_molecule("Molecule", file=resfile)
                count = cast(List[int], sec["counts"])[iState]
                if not cast(List[bool], sec["isTS"])[iState]:
                    self._states.append(AMSResults.EnergyLandscape.State(self, resfile, energy, mol, count, False))
                else:
                    reactantsID = (
                        cast(List[int], sec["reactants"])[iState]
                        if cast(List[int], sec["reactants"])[iState] > 0
                        else None
                    )
                    productsID = (
                        cast(List[int], sec["products"])[iState]
                        if cast(List[int], sec["products"])[iState] > 0
                        else None
                    )
                    prefactorsFromReactant = (
                        cast(List[int], sec["prefactorsFromReactant"])[iState]
                        if cast(List[int], sec["products"])[iState] > 0
                        else None
                    )
                    prefactorsFromProduct = (
                        cast(List[float], sec["prefactorsFromProduct"])[iState]
                        if cast(List[bool], sec["products"])[iState] > 0
                        else None
                    )
                    self._states.append(
                        AMSResults.EnergyLandscape.State(
                            self,
                            resfile,
                            energy,
                            mol,
                            count,
                            True,
                            reactantsID,
                            productsID,
                            prefactorsFromReactant,
                            prefactorsFromProduct,
                        )
                    )

            if "nFragments" not in sec:
                return

            nFragments = cast(int, sec["nFragments"])

            for iFragment in range(nFragments):
                energy = cast(List[float], sec["fragmentsEnergies"])[iFragment]
                resfile = os.path.splitext(cast(str, sec["fragmentsFileNames"]).split("\0")[iFragment])[0]
                mol = results.get_molecule("Molecule", file=resfile)
                self._fragments.append(AMSResults.EnergyLandscape.Fragment(self, resfile, energy, mol))

            if "nFStates" not in sec:
                return

            nFragmentedStates = cast(int, sec["nFStates"])

            for iFState in range(nFragmentedStates):
                iEnergy = cast(float, sec["fStatesEnergy(" + str(iFState + 1) + ")"])

                value = cast(Union[List[int], int], sec["fStatesComposition(" + str(iFState + 1) + ")"])
                iComposition = [i - 1 for i in value] if isinstance(value, list) else [value - 1]

                value = cast(Union[List[int], int], sec["fStatesConnections(" + str(iFState + 1) + ")"])
                iConnections = [i - 1 for i in value] if isinstance(value, list) else [value - 1]

                value = cast(Union[List[int], int], sec["fStatesAdsorptionPrefactors(" + str(iFState + 1) + ")"])
                iAdsorptionPrefactors = value if isinstance(value, list) else [value]

                value = cast(Union[List[int], int], sec["fStatesDesorptionPrefactors(" + str(iFState + 1) + ")"])
                iDesorptionPrefactors = value if isinstance(value, list) else [value]

                self._fstates.append(
                    AMSResults.EnergyLandscape.FragmentedState(
                        self, iEnergy, iComposition, iConnections, iAdsorptionPrefactors, iDesorptionPrefactors
                    )
                )

        @property
        def minima(self) -> List["AMSResults.EnergyLandscape.State"]:
            """Return all local-minimum states.

            :return: local minima in landscape order
            """
            return [s for s in self._states if not s.isTS]

        @property
        def transition_states(self) -> List["AMSResults.EnergyLandscape.State"]:
            """Return all transition states.

            :return: transition states in landscape order
            """
            return [s for s in self._states if s.isTS]

        @property
        def fragments(self) -> List["AMSResults.EnergyLandscape.Fragment"]:
            """Return all isolated fragments.

            :return: isolated fragments in landscape order
            """
            return [f for f in self._fragments]

        @property
        def fragmented_states(self) -> List["AMSResults.EnergyLandscape.FragmentedState"]:
            """Return all fragmented states.

            :return: fragmented states in landscape order
            """
            return [fs for fs in self._fstates]

        def __str__(self) -> str:
            lines = ["All stationary points:"]
            lines += ["======================"]
            for s in self._states:
                lines += [str(s)]
            for f in self._fragments:
                lines += [str(f)]
            for fs in self._fstates:
                lines += [str(fs)]
            return "\n".join(lines)

        def __getitem__(self, i: int) -> "AMSResults.EnergyLandscape.State":
            return self._states[i - 1]

        def __iter__(self) -> Iterator["AMSResults.EnergyLandscape.State"]:
            return iter(self._states)

        def __len__(self) -> int:
            return len(self._states)

    def get_energy_landscape(self) -> "AMSResults.EnergyLandscape":
        """Return the energy landscape from an AMS PES exploration.

        The returned object is of the type ``AMSResults.EnergyLandscape`` and offers convenient access to the states' energies, geometries as well as the information on which transition states connect which minima.
        State and fragment energies are in hartree, while molecular coordinates are in angstrom.

        .. code-block:: python

            el = results.get_energy_landscape()
            print(el)

            for state in el:
                print(f"Energy = {state.energy} hartree")
                print(f"Is transition state = {state.isTS}")
                print("Geometry:", state.molecule)
                if state.isTS:
                    print(f"Forward  barrier: {state.energy - state.reactants.energy} hartree")
                    print(f"Backward barrier: {state.energy - state.products.energy} hartree")

        :return: stationary points, fragments, and their connectivity
        """
        return AMSResults.EnergyLandscape(self)

    # =========================================================================

    def _access_rkf(self, func: Callable[[KFFile], T], file: str = "ams") -> T:
        """Apply a function to a selected RKF file.

        :param func: callable receiving the selected |KFFile|
        :param file: RKF identifier. Use ``engine`` to select the unique engine result file.
        :return: value returned by ``func``
        :raises FileError: if the requested RKF file is absent
        :raises ValueError: if ``engine`` does not identify a unique file
        """
        # Try unique engine:
        if file == "engine":
            names = self.engine_names()
            if len(names) == 1:
                return func(self.rkfs[names[0]])
            else:
                raise ValueError(
                    "You cannot use 'engine' as 'file' argument if the engine results file is not unique. Please use the real name of the file you wish to read"
                )

        # Try:
        if file in self.rkfs:
            return func(self.rkfs[file])

        # Try harder:
        filename = file + ".rkf"
        self.refresh()
        if filename in self.files:
            self.rkfs[file] = KFFile(str(self.job.get_path() / filename))
            return func(self.rkfs[file])

        # Surrender:
        raise FileError("File {} not present in {}".format(filename, self.job.get_path()))

    def _process_engine_results(self, func: Callable[[KFFile], T], engine: Optional[str] = None) -> T:
        """Apply a function to an engine result RKF file.

        :param func: callable receiving the selected |KFFile|
        :param engine: engine RKF identifier, or ``None`` to select the unique engine result file
        :return: value returned by ``func``
        :raises FileError: if the requested engine result file is absent
        :raises ValueError: if no unique engine result file can be selected
        """
        names = self.engine_names()
        if engine is not None:
            if engine in names:
                return func(self.rkfs[engine])
            else:
                raise FileError(
                    f"File {engine}.rkf not present in {self.job.get_path()}\n engine names found are: {names}"
                )
        else:
            if len(names) == 1:
                return func(self.rkfs[names[0]])
            elif len(names) == 0:
                raise FileError("There is no engine .rkf present in {}".format(self.job.get_path()))
            else:
                raise ValueError(
                    "You need to specify the 'engine' argument when there are multiple engine result files present in the job folder"
                )


# ===========================================================================
# ===========================================================================
# ===========================================================================


class AMSJob(SingleJob):
    """A single computation with the AMS driver.

    :param molecule: input molecule, chemical system, or mapping of named systems
    :param args: positional arguments forwarded to |SingleJob|
    :param kwargs: keyword arguments forwarded to |SingleJob|
    """

    results: AMSResults
    molecule: Optional[Union[Molecule, Dict[str, Molecule], "ChemicalSystem", Dict[str, "ChemicalSystem"]]]

    _result_type = AMSResults
    _command = "ams"

    def __init__(
        self,
        molecule: Optional[Union[Molecule, Dict[str, Molecule], "ChemicalSystem", Dict[str, "ChemicalSystem"]]] = None,
        *args: Any,
        **kwargs: Any,
    ):
        def copy_mol(mol: T) -> T:
            if isinstance(mol, Molecule):
                return mol.copy()  # type: ignore[return-value]
            elif _has_scm_chemsys and isinstance(mol, ChemicalSystem):
                return mol.copy()
            else:
                return mol

        if molecule is not None:
            molecule = (
                {k: copy_mol(m) for k, m in molecule.items()} if isinstance(molecule, dict) else copy_mol(molecule)  # type: ignore[arg-type]
            )
        super().__init__(molecule, *args, **kwargs)

    def run(
        self,
        jobrunner: Optional["JobRunner"] = None,
        jobmanager: Optional["JobManager"] = None,
        watch: bool = False,
        **kwargs: Any,
    ) -> AMSResults:
        """Run the job using *jobmanager* and *jobrunner* (or defaults, if ``None``).

        If *watch* is set to ``True``, the contents of the AMS driver logfile will be forwarded line by line to the PLAMS logfile (and stdout), allowing for an easier monitoring of the running job.
        Note that forwarding the AMS driver logfile makes this method block until execution finishes,
        even when using a parallel |JobRunner|.

        Other keyword arguments (*\\*\\*kwargs*) are stored in ``run`` branch of job's settings.

        :param jobrunner: runner used to execute the job, or ``None`` for the configured default
        :param jobmanager: manager used for the job, or ``None`` for the configured default
        :param watch: forward new AMS log lines while the job runs
        :param kwargs: values stored in the ``run`` branch of the job settings
        :return: results associated with this job
        """

        if _has_watchdog and watch:
            jobmanager = get_config().default_jobmanager
            observer = Observer()
            event_handler = AMSJobLogTailHandler(self, jobmanager)
            observer.schedule(event_handler, jobmanager.workdir, recursive=True)
            observer.start()

            try:
                results = cast(AMSResults, super().run(jobrunner=jobrunner, jobmanager=jobmanager, **kwargs))
                results.wait()
                event_handler.trigger()
            finally:
                observer.stop()
                observer.join()

        else:
            results = cast(AMSResults, super().run(jobrunner=jobrunner, jobmanager=jobmanager, **kwargs))

        return results

    def get_input(self) -> str:
        """Generate the input file. This method is just a wrapper around :meth:`_serialize_input`.

        Each instance of |AMSJob| or |AMSResults| present as a value in ``settings.input`` branch is replaced with an absolute path to ``ams.rkf`` file of that job.

        If you need to use a path to some engine specific ``.rkf`` file rather than the main ``ams.rkf`` file, you can to it by supplying a tuple ``(x, name)`` where ``x`` is an instance of |AMSJob| or |AMSResults| and ``name`` is a string with the name of the ``.rkf`` file you want. For example, ``(myjob, 'dftb')`` will transform to the absolute path to ``dftb.rkf`` file in ``myjob``'s folder, if such a file is present.

        Instances of |KFFile| are replaced with absolute paths to corresponding files.

        :return: serialized AMS input
        """

        special = {
            AMSJob: lambda x: x.results.rkfpath(),
            AMSResults: lambda x: x.rkfpath(),
            KFFile: lambda x: x.path,
            tuple: lambda x: AMSJob._tuple2rkf(x),
        }
        return self._serialize_input(special)

    def get_runscript(self) -> str:
        """Generate the runscript. Returned string is of the form::

            unset AMS_SWITCH_LOGFILE_AND_STDOUT
            AMS_JOBNAME=jobname AMS_RESULTSDIR=. $AMSBIN/ams [-n nproc] --input=jobname.in [>jobname.out]

        ``-n`` flag is added if ``settings.runscript.nproc`` exists. ``[>jobname.out]`` is used based on ``settings.runscript.stdout_redirect``. If ``settings.runscript.preamble_lines`` exists, those lines will be added to the runscript verbatim before the execution of AMS. If ``settings.runscript.postamble_lines`` exists, those lines will be added to the runscript verbatim after the execution of AMS.

        :return: shell runscript for the job
        """
        ret = "unset AMS_SWITCH_LOGFILE_AND_STDOUT\n"
        ret += "unset SCM_LOGFILE\n"
        config_slurm = get_config().slurm
        if "nnode" in self.settings.runscript and config_slurm:
            # Running as a SLURM job step and user specified the number of nodes explicitly.
            ret += f'export SCM_SRUN_OPTIONS="$SCM_SRUN_OPTIONS -N {self.settings.runscript.nnode}"\n'
        elif "nproc" in self.settings.runscript and config_slurm:
            # Running as a SLURM job step and user asked for a specific number of tasks.
            # Make sure to use as few nodes as possible to avoid distributing jobs needlessly across nodes.
            # See: https://stackoverflow.com/questions/71382578
            nnode = 1
            for nnode in range(1, len(config_slurm.tasks_per_node) + 1):
                if sum(config_slurm.tasks_per_node[0:nnode]) >= self.settings.runscript.nproc:
                    break
            nn_flag = f"1-{nnode}" if nnode > 1 else f"{nnode}"
            ret += f'export SCM_SRUN_OPTIONS="$SCM_SRUN_OPTIONS -N {nn_flag}"\n'
        if _has_scm_pisa and isinstance(self.settings.input, DriverBlock):
            if self.settings.input.Engine.name in ("QuantumESPRESSO", "VASP"):
                ret += "export SCM_DISABLE_MPI=1\n"
        else:
            if "QuantumEspresso" in self.settings.input or "VASP" in self.settings.input:
                ret += "export SCM_DISABLE_MPI=1\n"
        if "preamble_lines" in self.settings.runscript:
            for line in self.settings.runscript.preamble_lines:
                ret += f"{line}\n"
        ret += 'AMS_JOBNAME="{}" AMS_RESULTSDIR=. $AMSBIN/ams'.format(self.name)
        if "nproc" in self.settings.runscript:
            ret += " -n {}".format(self.settings.runscript.nproc)
        ret += ' --input="{}" < /dev/null'.format(self._filename("inp"))
        if self.settings.runscript.stdout_redirect:
            ret += ' >"{}"'.format(self._filename("out"))
        ret += "\n\n"
        if "postamble_lines" in self.settings.runscript:
            for line in self.settings.runscript.postamble_lines:
                ret += f"{line}\n"
            ret += "\n"
        return ret

    def check(self) -> bool:
        """Check whether the main RKF file reports normal termination.

        :return: whether the job terminated normally without reported errors
        """
        try:
            status = cast(str, self.results.readrkf("General", "termination status"))
        except (FileError, KeyError) as e:
            log(str(e), 1)
            return False
        except:
            log(f"Could not read termination status from file {self.results.rkfpath()}", 1)
            return False

        if status is None:
            log(f"Could not read termination status from file {self.results.rkfpath()}", 1)
            return False

        if "NORMAL TERMINATION" in status:
            if "errors" in status:
                log("Job {} reported errors. Please check the output".format(self._full_name()), 1)
                return False
            if "warnings" in status:
                log("Job {} reported warnings. Please check the output".format(self._full_name()), 1)
            return True
        return False

    def get_errormsg(self) -> Optional[str]:
        """Return an error message for a failed job.

        :return: error message, or ``None`` for a successful job
        """
        if self.check():
            return None
        else:
            # Check if there is an error captured during the job process, or a previously cached error
            if self._error_msg:
                return self._error_msg

            default_msg = "Could not determine error message. Please check the output manually."
            msg = None
            try:
                # If not, the first place to check is the termination status on the ams.rkf.
                # If the AMS driver stopped with a known error (called StopIt in the Fortran code), the error will be in there.
                # Status can be:
                # - NORMAL TERMINATION with errors: find the error from the ams log file
                # - IN PROGRESS: probably means AMS was shut down hard from the outside
                #                e.g. it got SIGKILL from the scheduler for exceeding some resource limit
                #                find the last error from the stderr
                # Note AMS can crash before even creating an rkf, then can just check the output and error files.
                try:
                    termination_status = cast(str, self.results.readrkf("General", "termination status"))
                except FileError:
                    termination_status = None

                # First look for the last error in the logfile
                try:
                    log_err_lines = self.results.grep_file("ams.log", "ERROR: ")
                    if log_err_lines:
                        self._error_msg: Optional[str] = log_err_lines[-1].partition("ERROR: ")[2]
                        return self._error_msg
                except FileError:
                    pass

                # Then for a licensing issue, check the output logs directly
                try:
                    license_err_lines = self.results.get_output_chunk(
                        begin="LICENSE INVALID",
                        end="License file",
                        inc_begin=True,
                        inc_end=True,
                        match=1,
                    )
                    if license_err_lines:
                        self._error_msg = str.join("\n", license_err_lines)
                        return self._error_msg
                except FileError:
                    pass

                # For any other issue fall back to the error file directly
                if "$JN.err" in self.results:
                    # If the status is still "IN PROGRESS", that probably means AMS was shut down hard from the outside.
                    # E.g. it got SIGKILL from the scheduler for exceeding some resource limit.
                    # In this case useful information may be found on stderr.
                    with open(self.results["$JN.err"]) as err:
                        errlines = err.read().splitlines()
                    for el in reversed(errlines):
                        if el != "" and not el.isspace():
                            msg = (
                                f"Termination status: {termination_status}. Message: {el} . "
                                f"Check the files in {self.path} for more details."
                            )
                            break
            except:
                pass

            # Cache error message if called again
            self._error_msg = msg if msg else default_msg

            return self._error_msg

    def hash_input(self) -> str:
        """Calculate the hash of the input file.

        All instances of |AMSJob| or |AMSResults| present as values in ``settings.input`` branch are replaced with hashes of corresponding job's inputs. Instances of |KFFile| are replaced with absolute paths to corresponding files.

        :return: SHA-256 hash of the serialized input
        """
        special = {
            AMSJob: lambda x: x.hash_input(),
            AMSResults: lambda x: x.job.hash_input(),
            KFFile: lambda x: x.path,
            tuple: lambda x: AMSJob._tuple2rkf(x),
        }
        return sha256(self._serialize_input(special))

    # =========================================================================

    def _serialize_input(self, special: Dict[Type, Callable[[Any], str]]) -> str:
        """Transform the contents of ``settings.input`` branch into string with blocks, keys and values.

        First, the contents of ``settings.input`` are extended with entries returned by :meth:`_serialize_molecule`. Then the contents of ``settings.input.ams`` are used to generate AMS text input. Finally, every other (than ``ams``) entry in ``settings.input`` is used to generate engine specific input.

        Special values can be indicated with *special* argument, which should be a dictionary having types of objects as keys and functions translating these types to strings as values.

        :param special: mapping from special value types to string conversion functions
        :return: serialized AMS input
        """

        def unspec(value: T) -> Union[T, str]:
            """Convert a registered special value to a string.

            :param value: value to inspect
            :return: converted string or the original value
            """
            for spec_type in special:
                if isinstance(value, spec_type):
                    return special[spec_type](value)
            return value

        def serialize(key: str, value: Any, indent: int, end: str = "End") -> str:
            """Serialize one key-value pair from an input |Settings| instance.

            If the value is a nested |Settings| instance, use recursive calls to build the snippet for the entire block. Indent the result with *indent* spaces.

            :param key: input keyword or block name
            :param value: value or nested block to serialize
            :param indent: indentation width in spaces
            :param end: closing keyword for a block
            :return: serialized input snippet
            """
            ret = ""
            if isinstance(value, Settings):
                # Open block ...
                ret += " " * indent + key
                # ... with potential header
                if "_h" in value:
                    ret += " " + unspec(value["_h"])
                ret += "\n"

                # Free block, or explicitly placed entries with _1, _2, ...
                i = 1
                while f"_{i}" in value:
                    if isinstance(value[f"_{i}"], Settings):
                        for ckey in value[f"_{i}"]:
                            ret += serialize(ckey, value[f"_{i}"][ckey], indent + 2)
                    else:
                        ret += serialize("", value["_" + str(i)], indent + 2)
                    i += 1

                # Figure out the order in which we should serialize the entries in the block
                if "_o" in value:
                    # Ordered block: we need to serialize the contents in the order given by "_o"
                    children = []
                    for v in value["_o"]:
                        ckey, _, idx = v.partition("[")
                        if idx:
                            # Child is a recurring entry, e.g. SubBlock[5]
                            idx = int(idx.rstrip("]"))
                            children.append((ckey, value[ckey][idx - 1]))
                        else:
                            # Child is a unique entry
                            children.append((ckey, value[ckey]))
                else:
                    # Unordered block: normal Python order when iterating over a dict/Settings
                    children = [
                        (ckey, value[ckey]) for ckey in value if isinstance(ckey, str) and not ckey.startswith("_")
                    ]

                # Serialize all children
                for ckey, cvalue in children:
                    split_key = key.lower().split()
                    split_ckey = ckey.lower().split()
                    if len(split_key) > 0 and split_key[0] == "engine" and ckey.lower() == "input":
                        ret += serialize(ckey, cvalue, indent + 2, "EndInput")
                    # REB: For the hybrid engine. How to deal with the space in ckey (Engine DFTB)? Replace by underscore?
                    elif len(split_ckey) > 0 and split_ckey[0] == "engine":
                        engine = " ".join(ckey.split("_"))
                        ret += serialize(engine, cvalue, indent + 2, end="EndEngine") + "\n"
                    else:
                        ret += serialize(ckey, cvalue, indent + 2)

                # Close block
                if key.lower() == "input":
                    end = "endinput"
                ret += " " * indent + end + "\n"

            elif isinstance(value, list):
                for el in value:
                    ret += serialize(key, el, indent, end)
            elif value == "" or value is True:
                ret += " " * indent + key + "\n"
            elif value is False or value is None:
                pass
            else:
                ret += " " * indent + key + " " + str(unspec(value)) + "\n"
            return ret

        def merge_system_blocks(
            input_settings: Settings, systems: Dict[str, Union[Molecule, "ChemicalSystem"]]
        ) -> Dict[str, str]:
            """Combine parsed System blocks with explicitly supplied systems.

            Given a ``settings.input`` |Settings| object, combine the System blocks with any systems explicitly specified
            via the input molecule. These can either be PLAMS Molecules or Chemical Systems, but in both cases
            they are first serialised to a Settings instance, then to text.

            :param input_settings: input settings containing optional System blocks
            :param systems: explicit systems keyed by System-block header
            :return: merged System blocks keyed by header
            """
            ams_key = input_settings.find_case("ams")
            system_key = input_settings[ams_key].find_case("system")

            input_systems_settings: Dict[str, Settings] = {}
            for item in input_settings[ams_key]:
                if item == system_key:
                    input_system_settings = input_settings[ams_key][item]
                    input_system_settings_list = (
                        input_system_settings if isinstance(input_system_settings, list) else [input_system_settings]
                    )
                    input_systems_settings = {
                        s._h if "_h" in s else "": s for s in input_system_settings_list if isinstance(s, Settings)
                    }

            systems_settings = {n: self._serialize_single_molecule(n, m) for n, m in systems.items()}

            merged_systems: Dict[str, str] = {}
            for name in set(input_systems_settings.keys()).union(set(systems_settings.keys())):
                input_system_settings = input_systems_settings.get(name, Settings())
                system_settings = systems_settings.get(name, Settings())
                merged_system_settings = input_system_settings + system_settings
                merged_system = serialize("System", merged_system_settings, 0)
                merged_systems[name] = merged_system

            return merged_systems

        systems: Dict[str, Union[Molecule, "ChemicalSystem"]]
        if isinstance(self.molecule, Molecule) or (_has_scm_chemsys and isinstance(self.molecule, ChemicalSystem)):
            systems = {"": self.molecule}
        elif isinstance(self.molecule, dict):
            systems = self.molecule
        elif self.molecule is None:
            systems = {}
        else:
            raise JobError(
                f"Incorrect 'molecule' attribute of job {self._full_name()}. 'molecule' should be a Molecule, a ChemicalSystem, a dictionary or None, and not {type(self.molecule).__name__}"
            )

        if _has_scm_pisa and isinstance(self.settings.input, DriverBlock):
            # AMS specific way of writing input files:
            #   self.settings.input is an input class that knows how to serialize itself to text input.

            # Generate initial input text
            input_class: DriverBlock = self.settings.input
            txtinp = input_class.get_input_string()
            has_input_systems = hasattr(input_class, "System") and input_class.System.value_changed

            # Add/update any systems using the input molecules
            if not has_input_systems:
                # if there are no system blocks in the input settings there is an optimisation whereby any
                # chemical systems can be serialised straight to text
                system_blocks: Dict[str, str] = {}
                for name, system in systems.items():
                    if _has_scm_chemsys and isinstance(system, ChemicalSystem):
                        system_input = system.__format__(f'in:name="{name}"' if name else "in")
                    else:
                        system_settings = AMSJob._serialize_single_molecule(name, system)
                        system_input = serialize("System", system_settings, 0)
                    system_blocks[name] = system_input
                    txtinp += "\n" + system_input
            elif systems:
                # otherwise have to go first via serialisation to settings, merge, then serialise to text
                # to avoid duplication replace the original text system block with the updated text
                sys_text = str(input_class.System)
                sys_text_updated = ""
                system_blocks = merge_system_blocks(input_class.to_settings(), systems)
                for system_block in system_blocks.values():
                    sys_text_updated += "\n" + system_block
                txtinp = txtinp.replace(sys_text, sys_text_updated)

        else:
            # Open-source PLAMS way of writing input files:
            #    self.settings.input is a Settings object that we need to serialize to text input.

            input_settings: Settings = self.settings.input.copy()

            ams = input_settings.find_case("ams")
            syst = input_settings[ams].find_case("system")

            txtinp = ""
            has_input_systems = False

            # contents of the 'ams' block (AMS input) go first, excluding the system blocks
            for item in self.settings.input[ams]:
                if item == syst:
                    has_input_systems = True
                    continue
                txtinp += serialize(item, input_settings[ams][item], 0) + "\n"

            # then process the system blocks
            if not has_input_systems:
                # if there are no system blocks in the input settings there is an optimisation whereby any
                # chemical systems can be serialised straight to text
                system_blocks = {}
                for name, system in systems.items():
                    if _has_scm_chemsys and isinstance(system, ChemicalSystem):
                        system_input = str(system) + "\n"
                        if name:
                            system_input = system_input.replace("System", f"System {name}", 1)
                    else:
                        system_settings = AMSJob._serialize_single_molecule(name, system)
                        system_input = serialize("System", system_settings, 0)
                    system_blocks[name] = system_input
            else:
                # otherwise have to go first via serialisation to settings, merge, then serialise to text
                system_blocks = merge_system_blocks(input_settings, systems)

            for system_block in system_blocks.values():
                txtinp += system_block

            # and then engines
            for engine in input_settings:
                if engine != ams:
                    txtinp += "\n" + serialize(f"Engine {engine}", input_settings[engine], 0, end="EndEngine") + "\n"

        return txtinp

    @staticmethod
    def _serialize_single_molecule(name: str, mol: Union[Molecule, "ChemicalSystem"]) -> Settings:
        def serialize_chemsys_to_settings(mol: "ChemicalSystem") -> Settings:
            from scm.plams.interfaces.adfsuite.inputparser import input_to_settings

            sett = input_to_settings(str(mol), program=AMSJob._command)
            return sett.ams.system[0]

        def serialize_molecule_to_settings(mol: Molecule, name: Optional[str] = None) -> Settings:
            sett = Settings()

            if len(mol.lattice) in [1, 2] and mol.align_lattice():
                log(
                    "The lattice of {} Molecule supplied did not follow the convention required by AMS. "
                    "I rotated the whole system for you. You're welcome".format(name if name else "main"),
                    3,
                )

            sett.Atoms._1 = [
                atom.str(symbol=AMSJob._atom_symbol(atom), space=18, decimal=10, suffix=AMSJob._atom_suffix(atom))
                for atom in mol
            ]

            if mol.lattice:
                sett.Lattice._1 = ["{:16.10f} {:16.10f} {:16.10f}".format(*vec) for vec in mol.lattice]

            if len(mol.bonds) > 0:
                lines = ["{} {} {}".format(mol.index(b.atom1), mol.index(b.atom2), b.order) for b in mol.bonds]
                # Add bond properties if they are defined
                sett.BondOrders._1 = [
                    (
                        "{} {}".format(text, mol.bonds[i].properties.suffix)
                        if "suffix" in mol.bonds[i].properties
                        else text
                    )
                    for i, text in enumerate(lines)
                ]

            if "charge" in mol.properties:
                sett.Charge = mol.properties.charge

            if "regions" in mol.properties:
                # sett.Region = [Settings({'_h':name, '_1':['%s=%s'%(k,str(v)) for k,v in data.items()]}) for name, data in molecule.properties.regions.items()]
                sett.Region = [
                    Settings({"_h": name, "Properties": Settings({"_1": [line for line in data]})})
                    for name, data in mol.properties.regions.items()
                    if isinstance(data, list)
                ]

            return sett

        if isinstance(mol, Molecule):
            sett = serialize_molecule_to_settings(mol, name)
        elif _has_scm_chemsys and isinstance(mol, ChemicalSystem):
            sett = serialize_chemsys_to_settings(mol)
        else:
            raise PlamsError(f"Cannot serialize molecule of type {type(mol).__name__} to settings.")

        if name:
            sett._h = name

        return sett

    def _serialize_molecule(self) -> List[Settings]:
        """Return a list of |Settings| instances containing the information about one or more |Molecule| instances stored in the ``molecule`` attribute.

        Molecular charge is taken from ``molecule.properties.charge``, if present. Additional, atom-specific information to be put in ``atoms`` block after XYZ coordinates can be supplied with ``atom.properties.suffix``.

        If the ``molecule`` attribute is a dictionary, the returned list is of the same length as the size of the dictionary. Keys from the dictionary are used as headers of returned ``system`` blocks.

        :return: settings representing the job's System blocks
        """

        if self.molecule is None:
            return []

        if isinstance(self.molecule, Molecule) or (_has_scm_chemsys and isinstance(self.molecule, ChemicalSystem)):
            moldict = {"": self.molecule}
        elif isinstance(self.molecule, dict):
            moldict = self.molecule
        else:
            raise JobError(
                f"Incorrect 'molecule' attribute of job {self._full_name()}. 'molecule' should be a Molecule, a ChemicalSystem, a dictionary or None, and not {type(self.molecule).__name__}"
            )

        ret = [AMSJob._serialize_single_molecule(name, molecule) for name, molecule in moldict.items()]

        return ret

    # =========================================================================

    def get_task(self) -> Optional[str]:
        """Return the AMS task from the job settings.

        :return: task name, or ``None`` if it is not defined
        """
        if (
            isinstance(self.settings, Settings)
            and "input" in self.settings
            and "ams" in self.settings.input
            and "task" in self.settings.input.ams
        ):
            return self.settings.input.ams.task
        return None

    @staticmethod
    def _atom_suffix(atom: Atom) -> str:
        """Return atomic properties formatted for an AMS Atoms block.

        :param atom: atom whose properties are serialized
        :return: suffix containing the serialized properties
        :raises ValueError: if a property cannot be represented safely in AMS input
        """

        # Build key-value dictionary from properties
        keyval_dict = {}

        if _has_scm_chemsys:
            from scm.base import AtomAttributes

            allowed_attributes = AtomAttributes.Groups + ["mass", "region"]

            def skip_attribute(prefix: str, key: str) -> bool:
                return prefix == "" and key.lower() not in allowed_attributes

        else:

            def skip_attribute(prefix: str, key: str) -> bool:
                """
                Skip special atomic properties that are handled by _atom_symbol() already (handled explicitly below)
                or internal PLAMS properties which are not accepted as valid atom properties in AMS, and so can be pruned out.

                :param prefix: nested property prefix
                :param key: property key
                :return: whether the property should be omitted
                """
                return prefix == "" and key.lower() in ["suffix", "ghost", "name", "supercell", "rdkit"]

        def serialize(sett: Settings, prefix: str = "") -> None:
            for key, val in sett.items():
                if skip_attribute(prefix, key):
                    continue
                if isinstance(val, Settings):
                    # Recursively serialize nested Settings object
                    serialize(val, prefix + key + ".")
                    continue
                elif isinstance(val, list):
                    # Lists are converted to a comma separated string of elements
                    val = ",".join(str(v) for v in val)
                elif isinstance(val, set):
                    # Sets are converted to a *sorted* comma separated string of elements
                    val = ",".join(str(v) for v in sorted(val))
                key = prefix + key
                if not isinstance(val, str):
                    val = str(val)
                if "\n" in val:
                    raise ValueError(f"String representation of atomic property {key} may not include line breaks")
                if "{" in val or "}" in val:
                    raise ValueError(
                        f"String representation of atomic property {key} may not include curly brackets: {val}"
                    )
                if " " in val:
                    # Ensure that space containing values are quoted
                    if '"' not in val:
                        val = '"' + val + '"'
                    elif "'" not in val:
                        val = "'" + val + "'"
                    else:
                        raise ValueError(f"Atomic property {key} can not include both single and double quotes: {val}")
                keyval_dict[key] = str(val)

        serialize(atom.properties)

        # Assemble and return complete suffix string
        ret = " ".join(f"{key}={val}" for key, val in keyval_dict.items())
        if "suffix" in atom.properties:
            ret += " " + str(atom.properties.suffix)
        return ret.strip()

    @staticmethod
    def _atom_symbol(atom: Atom) -> str:
        """Return an atomic symbol formatted for AMS input.

        :param atom: atom whose symbol is formatted
        :return: symbol including optional ghost and custom-name markers
        """
        smb = atom.symbol
        if "ghost" in atom.properties and atom.properties.ghost:
            smb = ("Gh." + smb).rstrip(".")
        if "name" in atom.properties:
            smb = (smb + "." + str(atom.properties.name)).lstrip(".")
        return smb

    @staticmethod
    def _tuple2rkf(arg: Tuple[Union["AMSJob", "AMSResults"], str]) -> str:
        """Resolve a job/result and RKF identifier pair to an absolute path.

        :param arg: pair containing an |AMSJob| or |AMSResults| and an RKF identifier
        :return: absolute RKF path, or the string representation of an invalid pair
        """
        if len(arg) == 2 and isinstance(arg[1], str):
            if isinstance(arg[0], AMSJob):
                return arg[0].results.rkfpath(arg[1])
            if isinstance(arg[0], AMSResults):
                return arg[0].rkfpath(arg[1])
        return str(arg)

    @classmethod
    def load_external(
        cls,
        path: str,
        settings: Optional[Settings] = None,
        molecule: Optional[Molecule] = None,
        finalize: bool = False,
        fmt: str = "ams",  # type: ignore[override]
    ) -> "AMSJob":
        """Load an external job from *path*.

        In this context an "external job" is an execution of some external binary that was not managed by PLAMS, and hence does not have a ``.dill`` file. It can also be used in situations where the execution was started with PLAMS, but the Python process was terminated before the execution finished, resulting in steps 9-12 of :ref:`job-life-cycle` not happening.

        All the files produced by your computation should be placed in one folder and *path* should be the path to this folder or a file in this folder. The name of the folder is used as a job name. Input, output, error and runscript files, if present, should have names defined in ``_filenames`` class attribute (usually ``[jobname].in``, ``[jobname].out``, ``[jobname].err`` and ``[jobname].run``). It is not required to supply all these files, but in most cases one would like to use at least the output file, in order to use methods like :meth:`~scm.plams.core.results.Results.grep_output` or :meth:`~scm.plams.core.results.Results.get_output_chunk`. If *path* is an instance of an AMSJob, that instance is returned.

        This method is a class method, so it is called via class object and it returns an instance of that class::

            >>> a = AMSJob.load_external(path='some/path/jobname')
            >>> type(a)
            scm.plams.interfaces.adfsuite.ams.AMSJob

        You can supply |Settings| and |Molecule| instances as *settings* and *molecule* parameters, they will end up attached to the returned job instance. If you don't do this, PLAMS will try to recreate them automatically using methods :meth:`~scm.plams.core.results.Results.recreate_settings` and :meth:`~scm.plams.core.results.Results.recreate_molecule` of the corresponding |Results| subclass. If no |Settings| instance is obtained in either way, the defaults from ``config.job`` are copied.

        You can set the *finalize* parameter to ``True`` if you wish to run the whole :meth:`~Job._finalize` on the newly created job. In that case PLAMS will perform the usual :meth:`~Job.check` to determine the job status (*successful* or *failed*), followed by cleaning of the job folder (|cleaning|), |postrun| and pickling (|pickling|). If *finalize* is ``False``, the status of the returned job is *copied*.

        :param path: results directory or a file within it
        :param settings: settings to attach, or ``None`` to reconstruct them when possible
        :param molecule: molecule to attach, or ``None`` to reconstruct it when possible
        :param finalize: run the normal job finalization steps after loading
        :param fmt: input format: ``ams``, ``qe``, ``gaussian``, ``vasp``, or ``any``
        :return: job representing the external results

        The available formats behave as follows:

        * ``ams`` loads a finished AMS job.
        * ``qe`` converts Quantum ESPRESSO output to ``ams.rkf`` and ``qe.rkf``.
        * ``gaussian`` converts Gaussian output to AMS-compatible RKF files.
        * ``vasp`` converts a VASP ``OUTCAR`` to ``ams.rkf`` and ``vasp.rkf``.
        * ``any`` detects the format automatically.

        This method can also be used to convert a finished VASP job to an AMSJob. If you supply the path to a folder containing OUTCAR, then a subdirectory will be created in this folder called AMSJob. In the AMSJob subdirectory, two files will be created: ams.rkf and vasp.rkf, that contain some of the results from the VASP calculation. If the AMSJob subdirectory already exists, the existing ams.rkf and vasp.rkf files will be reused. NOTE: the purpose of loading VASP data this way is to let you call for example job.results.get_energy() etc., not to run new VASP calculations!
        """
        if isinstance(path, cls):
            return path

        preferred_name = None

        fmt = fmt.lower()

        # first check if path is a VASP output
        # in which case call cls._vasp_to_ams, which will return a NEW path
        # containing ams.rkf and vasp.rkf (the .rkf files will be created if they do not exist)
        if os.path.isdir(path):
            preferred_name = os.path.basename(os.path.abspath(path))
            if (
                not os.path.exists(os.path.join(path, "ams.rkf"))
                and os.path.exists(os.path.join(path, "OUTCAR"))
                and (fmt == "vasp" or fmt == "any")
            ):
                path = vasp_output_to_ams(path, overwrite=False)
        elif os.path.exists(path) and os.path.basename(path) == "OUTCAR" and (fmt == "vasp" or fmt == "any"):
            preferred_name = os.path.basename(os.path.dirname(os.path.abspath(path)))
            path = vasp_output_to_ams(os.path.dirname(path), overwrite=False)
        elif os.path.exists(path):
            try:
                from ase.io.formats import filetype
            except ImportError:
                raise MissingOptionalPackageError("ase")

            try:
                ft = filetype(path)
                # check if path is a Quantum ESPRESSO .out file
                # qe_output_to_ams returns a NEW path containng the ams.rkf and qe.rkf files (will be created if
                # they do not exist)
                if fmt == "qe":
                    assert ft == "espresso-out", f"The file {path} does not seem to be a Quantum ESPRESSO output file"
                if fmt == "gaussian":
                    assert ft == "gaussian-out", f"The file {path} does not seem to be a Gaussian output file"
                if ft == "espresso-out" and (fmt == "qe" or (fmt == "any" and not path.endswith("rkf"))):
                    path = qe_output_to_ams(path, overwrite=False)
                elif ft == "gaussian-out" and (fmt == "gaussian" or (fmt == "any" and not path.endswith("rkf"))):
                    path = gaussian_output_to_ams(path, overwrite=False)
            except Exception:
                # several types of exceptions can be raised, e.g. StopIteration and UnicodeDecodeError
                # assume that any exception means that the file was not a QE output file
                pass

        if not os.path.isdir(path):
            if os.path.exists(path):
                path = os.path.dirname(os.path.abspath(path))
            elif os.path.isdir(path + ".results"):
                path = path + ".results"
            elif os.path.isdir(path + "results"):
                path = path + "results"
            else:
                raise FileError("Path {} does not exist, cannot load from it.".format(path))

        job = cast(AMSJob, super(AMSJob, cls).load_external(path, settings, molecule, finalize))

        if preferred_name is not None:
            job.name = preferred_name

        if job.name.endswith(".results") and len(job.name) > 8:
            job.name = job.name[:-8]

        return job

    @classmethod
    def from_input(cls, text_input: str, **kwargs: Any) -> "AMSJob":
        """
        Create an AMS job from AMS-style text input.

        :param text_input: complete AMS input
        :param kwargs: additional |AMSJob| constructor arguments
        :return: job containing the parsed settings and systems
        :raises ImportError: if the SCM input parser is unavailable
        :raises JobError: if molecules are supplied both in the input and as a keyword argument

        Example::

           text = '''
           Task GeometryOptimization
           Engine DFTB
               Model GFN1-xTB
           EndEngine
           '''

           job = AMSJob.from_input(text)

        .. note::

            If *molecule* is included in the keyword arguments to this method, the *text_input* may not contain any System blocks. In other words, the molecules to be used either need to come from the *text_input*, or the keyword argument, but not both.

            If *settings* is included in the keyword arguments to this method, the |Settings| created from the *text_input* will be soft updated with the settings from the keyword argument. In other words, the *text_input* takes precedence over the *settings* keyword argument.
        """
        from scm.plams.interfaces.adfsuite.inputparser import input_to_settings

        sett = Settings()
        sett.input = input_to_settings(text_input, program=cls._command)
        mol = cls.settings_to_mol(sett)
        if mol:
            if "molecule" in kwargs:
                raise JobError("AMSJob.from_input(): molecule passed in both text_input and as keyword argument")
        else:
            mol = kwargs.pop("molecule", None)
        if "settings" in kwargs:
            sett.soft_update(kwargs.pop("settings"))
        return cls(molecule=mol, settings=sett, **kwargs)

    @classmethod
    def from_inputfile(cls, filename: str, heredoc_delimit: str = "eor", **kwargs: Any) -> "AMSJob":
        """Construct an :class:`AMSJob` instance from an AMS input or run file.

        If a runscript is provided, this method attempts to extract the input based
        on the heredoc delimiter (see *heredoc_delimit*).

        :param filename: path to the input or run file
        :param heredoc_delimit: heredoc delimiter used by a run file
        :param kwargs: additional arguments forwarded to :meth:`from_input`
        :return: job containing the parsed input
        """
        with open(filename, "r") as f:
            inp_file = parse_heredoc(f.read(), heredoc_delimit)
        return cls.from_input(inp_file, **kwargs)

    @staticmethod
    def settings_to_mol(s: Settings) -> Optional[Dict[str, Molecule]]:
        """Remove the ``s.input.ams.system`` block and convert it to molecules.

        The provided settings should be in the same style as the ones produced by the SCM input parser.
        Dictionary keys are taken from the header of each system block.
        The existing `s.input.ams.system` block is removed in the process, assuming it was present in the first place.

        :param s: settings produced by the SCM input parser
        :return: molecules keyed by System-block header, or ``None`` if no System block is present
        :raises KeyError: if multiple System blocks use the same header
        """
        from scm.plams.tools.units import Units

        def get_list(s: Settings) -> List:
            if "_1" in s and not "_2" in s and isinstance(s._1, list):
                return s._1
            else:
                i = 1
                l = []
                while "_" + str(i) in s:
                    l.append(s["_" + str(i)])
                    i += 1
                return l

        def read_mol(settings_block: Settings) -> Optional[Molecule]:
            """Return the molecule represented by one System block.

            :param settings_block: parsed System block
            :return: reconstructed molecule
            """
            if "geometryfile" in settings_block:
                mol = Molecule(settings_block.geometryfile)
            else:
                mol = Molecule()

                if "_h" in settings_block.atoms:
                    conv = Units.conversion_ratio(settings_block.atoms._h.strip("[]"), "Angstrom")
                else:
                    conv = 1.0
                for atom in get_list(settings_block.atoms):
                    # Extract arguments for Atom()
                    symbol, x, y, z, *suffix = atom.split(maxsplit=4)
                    coords = conv * float(x), conv * float(y), conv * float(z)

                    try:
                        at = Atom(symbol=symbol, coords=coords)
                    except PTError:  # It's either a ghost atom and/or an atom with a custom name
                        kwargs = {}
                        if symbol.startswith("Gh."):  # Ghost atom
                            kwargs["ghost"] = True
                            _, symbol = symbol.split(".", maxsplit=1)
                        if "." in symbol:  # Atom with a custom name
                            symbol, kwargs["name"] = symbol.split(".", maxsplit=1)
                        at = Atom(symbol=symbol, coords=coords, **kwargs)  # type: ignore[arg-type]
                    if suffix:
                        at.properties.soft_update(AMSJob._atom_suffix_to_settings(suffix[0]))
                    mol.add_atom(at)

                # Set the lattice vector if applicable
                if get_list(settings_block.lattice):
                    if "_h" in settings_block.lattice:
                        conv = Units.conversion_ratio(settings_block.lattice._h.strip("[]"), "Angstrom")
                    else:
                        conv = 1.0
                    mol.lattice = [[conv * float(j) for j in i.split()] for i in get_list(settings_block.lattice)]

            # Add bonds
            for bond in get_list(settings_block.bondorders):
                _at1, _at2, _order, *suffix = bond.split(maxsplit=3)
                at1, at2, order = mol[int(_at1)], mol[int(_at2)], float(_order)
                plams_bond = Bond(at1, at2, order=order)
                if suffix:
                    plams_bond.properties.suffix = suffix[0]
                mol.add_bond(plams_bond)

            # Set the molecular charge as a numeric value
            if settings_block.charge:
                charge = float(settings_block.charge)
                mol.properties.charge = int(charge) if charge.is_integer() else charge

            # Set the region info (used in ACErxn)
            if settings_block.region:
                for s_reg in settings_block.region:
                    if not "_h" in s_reg.keys():
                        raise JobError("Region block requires a header!")
                    if s_reg.properties:
                        mol.properties.regions[s_reg._h] = s_reg.properties._1

            # Apply the supercell keyword
            if "supercell" in settings_block:
                if isinstance(settings_block.supercell, str):
                    mol = mol.supercell(*[int(d) for d in settings_block.supercell.split()])
                else:
                    mol = mol.supercell(*settings_block.supercell)
                for at in mol:
                    del at.properties.supercell

            mol.properties.name = str(settings_block.get("_h", ""))
            return mol

        # Raises a KeyError if the `system` key is absent
        with s.suppress_missing():
            try:
                settings_list = s.input.ams.system
                if not isinstance(settings_list, list):
                    settings_list = [settings_list]
            except KeyError:  # The block s.input.ams.system is absent
                return None

        # Create a new dictionary with system headers as keys and molecules as values
        moldict: Dict[str, Optional[Molecule]] = {}
        for settings_block in settings_list:
            key = (
                str(settings_block._h) if ("_h" in settings_block) else ""
            )  # Empty string used as default system name.
            if key in moldict:
                raise KeyError(f"Duplicate system headers found in s.input.ams.system: {repr(key)}")
            moldict[key] = read_mol(settings_block)
            used_properties = ["atoms", "bondorders", "geometryfile", "lattice", "charge", "region", "supercell"]
            for used_property in used_properties:
                if used_property in settings_block:
                    settings_block.pop(used_property)
        s.input.ams.System = [system for system in settings_block if system]
        if not len(s.input.ams.System):
            del s.input.ams.System

        return {k: v for k, v in moldict.items() if v is not None}

    @staticmethod
    def _atom_suffix_to_settings(suffix: str) -> Settings:
        def is_int(s: str) -> bool:
            if "." in s:
                return False
            try:
                return int(s) == float(s)
            except ValueError:
                return False

        def is_float(s: str) -> bool:
            try:
                float(s)
                return True
            except ValueError:
                return False

        import shlex

        tokens = shlex.split(suffix)

        # Remove possible spaces around = signs
        i = 0
        while i < len(tokens):
            if tokens[i].endswith("="):
                tokens[i] += tokens.pop(i + 1)
            elif tokens[i].startswith("="):
                tokens[i - 1] += tokens.pop(i)
            else:
                i += 1

        properties = Settings()
        for i, t in enumerate(tokens):
            val: Any
            key, _, val = t.partition("=")
            if not val:
                # Token that is not a key=value pair? Let's accumulate those in the plain text suffix.
                if "suffix" in properties:
                    properties.suffix += " " + key
                else:
                    properties.suffix = key
                continue
            else:
                # We have a value. Let's make an educated guess for its type.
                if val.lower() == "true":
                    val = True
                elif val.lower() == "false":
                    val = False
                elif is_int(val):
                    val = int(val)
                elif is_float(val):
                    val = float(val)
                elif "," in val:
                    elem = val.split(",")
                    if key.lower() == "region":
                        # For regions we will output a set instead of a list
                        val = set(elem)
                    elif all(is_int(i) for i in elem):
                        val = [int(i) for i in elem]
                    elif all(is_float(f) for f in elem):
                        val = [float(f) for f in elem]
                    else:
                        # Anything else will come out as a list of strings
                        val = elem
            properties.set_nested(tuple(key.split(".")), val)

        return properties

    @staticmethod
    def _add_region(atom: Atom, name: str) -> None:
        """Add a region name to an atom.

        :param atom: atom to update
        :param name: region name; if ``None``, only normalize the existing region value to a set
        """
        if "region" not in atom.properties:
            atom.properties.region = set()
        if isinstance(atom.properties.region, str):
            atom.properties.region = set([atom.properties.region])
        if not isinstance(atom.properties.region, set):
            atom.properties.region = set(atom.properties.region)
        if name is None:
            return
        atom.properties.region.add(name)


def extract_engine_settings(settings: Settings) -> Settings:
    """Extract the first engine branch from a settings object.

    :param settings: complete AMS job settings
    :return: settings containing only the first engine branch

    Example::

        s.input.ams.Task = "singlepoint"
        s.runscript.nproc = 1
        s.input.ForceField.Type = "UFF"
        extract_engine_settings(s)
    """

    ret = Settings()

    if "input" in settings:
        for e in settings.input:
            if e != "ams":
                ret.input[e] = settings.input[e].copy()
                break

    return ret


def hybrid_committee_engine_settings(settings_list: List[Settings]) -> Settings:
    """Create settings for a Hybrid committee that averages its subengines.

    :param settings_list: top-level settings for each subengine
    :return: top-level settings for the Hybrid committee engine
    """

    def get_partial_settings(x: int) -> Settings:
        original_settings = settings_list[x]
        pure_settings = None
        pure_engine_name = ""
        for e in original_settings.input:
            if e != "ams":
                pure_engine_name = e
                pure_settings = original_settings.input[e]
                break

        assert pure_settings is not None

        pure_settings._h = f"{pure_engine_name} Engine{x+1}"

        return pure_settings

    N = len(settings_list)
    s = Settings()
    s.input.Hybrid.Engine = [get_partial_settings(x) for x in range(N)]
    s.input.Hybrid.Energy.Term = [f"Factor={1/N} Region=* UseCappingAtoms=No EngineID=Engine{x+1}" for x in range(N)]
    s.input.Hybrid.Committee.Enabled = "Yes"
    s.input.Hybrid.TweakRequestForSubEngines = (
        "No"  # If this is yes, it explicitly turns off any request for anything other than Energies/forces/charges
    )
    s.runscript.preamble_lines = ["export OMP_NUM_THREADS=1"]
    for ss in settings_list:
        if "runscript" in ss and "preamble_lines" in ss.runscript:
            s.runscript.preamble_lines += ss.runscript.preamble_lines

    return s

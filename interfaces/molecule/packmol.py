import os
import subprocess
import tempfile
from typing import Any, Dict, List, Literal, Optional, Tuple, Union, overload, Sequence, TYPE_CHECKING
from collections import Counter

import numpy as np
from scm.plams.core.errors import MoleculeError
from scm.plams.core.private import saferun
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.periodic_table import PeriodicTable
from scm.plams.tools.units import Units
from scm.plams.interfaces.molecule.rdkit import readpdb, writepdb
from scm.plams.core.functions import requires_optional_package, log
from scm.plams.core.settings import Settings
from scm.plams.core.jobmanager import JobManager

if TYPE_CHECKING:
    try:
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem
    except ImportError:
        pass

__all__ = [
    "packmol",
    "packmol_around",
    "packmol_in_void",
    "packmol_on_slab",
    "packmol_microsolvation",
    "PackMolError",
]


def tolist(x):
    return x if isinstance(x, list) else [x]


class PackMolError(MoleculeError):
    pass


class PackMolStructure:
    def __init__(
        self,
        molecule: Molecule,
        n_molecules: Optional[int] = None,
        n_atoms: Optional[int] = None,
        box_bounds: Optional[List[float]] = None,
        density: Optional[float] = None,
        fixed: bool = False,
        sphere: bool = False,
    ):
        """

        Class representing a packmol structure.

        molecule: |Molecule|
            The molecule

        n_molecules: int
            The number of molecules to insert

        n_atoms: int
            An approximate number of atoms to insert

        box_bounds: list of float
            [xmin, ymin, zmin, xmax, ymax, zmax] in angstrom. The min values should all be 0, i.e. [0., 0., 0., xmax, ymax, zmax]

        density: float
            Density in g/cm^3

        fixed: bool
            Whether the structure should be fixed at its original coordinates.

        sphere: bool
            Whether the molecules should be packed in a sphere. The radius is determined by getting the volume from the box bounds! Cannot be combined with ``fixed`` (``fixed`` takes precedence).

        """
        self.molecule = molecule
        if fixed:
            assert n_molecules is None or n_molecules == 1
            assert density is None
            self.n_molecules = 1
            if molecule.lattice and len(molecule.lattice) == 3:
                vecs = np.array(molecule.lattice)
                positive_mask = vecs >= 0
                negative_mask = vecs <= 0
                minxyz = np.sum(vecs * negative_mask, axis=0)
                maxxyz = np.sum(vecs * positive_mask, axis=0)
                self.box_bounds = minxyz.tolist() + maxxyz.tolist()
                # self.box_bounds = [
                #     0.0,
                #     0.0,
                #     0.0,
                #     molecule.lattice[0][0],
                #     molecule.lattice[1][1],
                #     molecule.lattice[2][2],
                # ]
            else:
                self.box_bounds = None
            self.fixed = True
            self.sphere = False
        else:
            if box_bounds and density:
                if n_molecules or n_atoms:
                    raise ValueError("Cannot set all n_molecules or n_atoms together with (box_bounds AND density)")
                n_molecules = self._get_n_molecules_from_density_and_box_bounds(self.molecule, box_bounds, density)
            assert n_molecules is not None or n_atoms is not None
            if n_molecules is None:
                assert n_atoms is not None
                self.n_molecules = self._get_n_molecules(self.molecule, n_atoms)
            else:
                self.n_molecules = n_molecules
            assert box_bounds or density
            self.box_bounds = box_bounds or self._get_box_bounds(self.molecule, self.n_molecules, density)
            self.fixed = False
            self.sphere = sphere

    def _get_n_molecules_from_density_and_box_bounds(
        self, molecule: Molecule, box_bounds: List[float], density: float
    ) -> int:
        """density in g/cm^3"""
        molecule_mass = molecule.get_mass(unit="g")
        volume_ang3 = self.get_volume(box_bounds)
        volume_cm3 = volume_ang3 * 1e-24
        n_molecules = int(density * volume_cm3 / molecule_mass)
        return n_molecules

    def get_volume(self, box_bounds: Optional[Sequence[float]] = None) -> float:
        bb = box_bounds or self.box_bounds
        if bb is None:
            raise ValueError("Cannot call get_volume when box_bounds is None.")
        vol = (bb[3] - bb[0]) * (bb[4] - bb[1]) * (bb[5] - bb[2])
        return vol

    def _get_n_molecules(self, molecule: Molecule, n_atoms: int):
        return n_atoms // len(molecule)

    def _get_box_bounds(self, molecule: Molecule, n_molecules: int, density: float):
        mass = n_molecules * molecule.get_mass(unit="g")
        volume_cm3 = mass / density
        volume_ang3 = volume_cm3 * 1e24
        side_length = volume_ang3 ** (1 / 3.0)
        return [0.0, 0.0, 0.0, side_length, side_length, side_length]

    def get_input_block(self, fname, tolerance):
        if self.n_molecules == 0 and not self.fixed:
            return ""
        if self.fixed:
            ret = f"""
            structure {fname}
            number 1
            fixed 0. 0. 0. 0. 0. 0.
            avoid_overlap yes
            end structure
            """
        elif self.sphere:
            vol = self.get_volume()
            radius = np.cbrt(3 * vol / (4 * np.pi))
            ret = f"""
            structure {fname}
              number {self.n_molecules}
              inside sphere 0. 0. 0. {radius}
            end structure
            """
        else:
            box_string = f"{self.box_bounds[0]+tolerance/2} {self.box_bounds[1]+tolerance/2} {self.box_bounds[2]+tolerance/2} {self.box_bounds[3]-tolerance/2} {self.box_bounds[4]-tolerance/2} {self.box_bounds[5]-tolerance/2}"
            ret = f"""
            structure {fname}
              number {self.n_molecules}
              inside box {box_string}
            end structure

        """
        return ret


class PackMol:
    def __init__(
        self,
        tolerance=2.0,
        structures: Optional[List[PackMolStructure]] = None,
        filetype="xyz",
        seed: int = -1,
        executable=None,
    ):
        """
        Class for setting up and running packmol.

        tolerance: float
            The packmol tolerance (approximate minimum interatomic distance)

        structures: list of PackMolStructure
            Structures to insert

        filetype: str
            One of 'xyz' or 'pdb'. Specifies the file format to use with packmol. 'pdb' requires rdkit.

        executable: str
            Path to the packmol executable. If not specified, $AMSBIN/packmol.exe will be used.

        seed: int
            Random seed. If -1, the current time is used as a random seed in packmol.

        Note: users are not recommended to use this class directly, but
        instead use the ``packmol``, ``packmol_on_slab`` and ``packmol_microsolvation``
        functions.

        """
        self.tolerance = tolerance
        self.structures = structures or []
        self.filetype = filetype
        if seed == -1 and "SCM_PACKMOL_SEED" in os.environ:
            self.seed = int(os.environ["SCM_PACKMOL_SEED"])
        else:
            self.seed = seed
        self.executable = executable or os.path.join(os.path.expandvars("$AMSBIN"), "packmol.exe")
        if not os.path.exists(self.executable):
            raise RuntimeError("PackMol exectuable not found: " + self.executable)

    def add_structure(self, structure: PackMolStructure):
        self.structures.append(structure)

    def _get_complete_box_bounds(self) -> Tuple[float, float, float, float, float, float]:
        min_x = min(s.box_bounds[0] for s in self.structures if s.box_bounds is not None)
        min_y = min(s.box_bounds[1] for s in self.structures if s.box_bounds is not None)
        min_z = min(s.box_bounds[2] for s in self.structures if s.box_bounds is not None)
        max_x = min(s.box_bounds[3] for s in self.structures if s.box_bounds is not None)
        max_y = min(s.box_bounds[4] for s in self.structures if s.box_bounds is not None)
        max_z = min(s.box_bounds[5] for s in self.structures if s.box_bounds is not None)

        # return min_x, min_y, min_z, max_x+self.tolerance, max_y+self.tolerance, max_z+self.tolerance
        return min_x, min_y, min_z, max_x, max_y, max_z

    def _get_complete_lattice(self) -> List[List[float]]:
        """
        returns a 3x3 list using the smallest and largest x/y/z box_bounds for all structures
        """
        if any(s.sphere for s in self.structures):
            return []
        (
            min_x,
            min_y,
            min_z,
            max_x,
            max_y,
            max_z,
        ) = self._get_complete_box_bounds()
        return [
            [max_x - min_x, 0.0, 0.0],
            [0.0, max_y - min_y, 0.0],
            [0.0, 0.0, max_z - min_z],
        ]

    def _get_complete_radius(self) -> float:
        """
        Calculates radius of sphere with the same volume as the
        cuboid from the box bounds

        :return: Radius in angstrom
        :rtype: float
        """
        volume = self._get_complete_volume()
        radius = np.cbrt(3 * volume / (4 * np.pi))

        return radius

    def _get_complete_volume(self) -> float:
        """Returns volume based on box bounds in ang^3

        :return: Volume in ang^3
        :rtype: float
        """
        (
            min_x,
            min_y,
            min_z,
            max_x,
            max_y,
            max_z,
        ) = self._get_complete_box_bounds()

        volume = (max_x - min_x) * (max_y - min_y) * (max_z - min_z)

        return volume

    def run(self):
        """
        returns: a Molecule with the packed structures
        """

        output_molecule = Molecule()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_fname = os.path.join(tmpdir, "output.xyz")
            input_fname = os.path.join(tmpdir, "input.inp")
            with open(input_fname, "w") as input_file:
                input_file.write(f"tolerance {self.tolerance}\n")
                input_file.write(f"seed {self.seed}\n")
                input_file.write(f"filetype {self.filetype}\n")
                input_file.write(f"output {output_fname}\n")

                for i, structure in enumerate(self.structures):
                    structure_fname = os.path.join(tmpdir, f"structure{i}.{self.filetype}")
                    if self.filetype == "pdb":
                        with open(structure_fname, "w") as f:
                            writepdb(structure.molecule, f)
                    else:
                        structure.molecule.write(structure_fname)
                    input_file.write(structure.get_input_block(structure_fname, tolerance=self.tolerance))

            with open(input_fname) as my_input:
                saferun(self.executable, stdin=my_input, stdout=subprocess.DEVNULL)

            if not os.path.exists(output_fname):
                raise PackMolError("Packmol failed. It may work if you try a lower density.")

            if self.filetype == "pdb":
                with open(output_fname, "r") as f:
                    output_molecule = readpdb(f)
            else:
                output_molecule = Molecule(output_fname)  # without periodic boundary conditions

            output_molecule.lattice = self._get_complete_lattice()

        return output_molecule


def sum_of_atomic_volumes(molecule: Molecule) -> float:
    """Returns the sum of atomic volumes (calculated using vdW radii) in angstrom^3."""
    return (4 / 3) * np.pi * sum(at.radius**3 for at in molecule)


def guess_density(molecules: Sequence[Molecule], coeffs: Sequence[Union[int, float]]) -> float:
    """Guess a density for a liquid of the given molecules and stoichiometric coefficients.

    This density is most of the time lower than the experimental density and is NOT meant
    to be accurate. Always equilibrate the density with NPT MD after creating a box with a
    density guessed using this method.

    Returns the guessed density in g/cm^3.
    """
    if len(molecules) != len(coeffs):
        raise ValueError(f"Incompatible lengths: {len(molecules)=}, {len(coeffs)=}")

    bond_volume_decrease_factor, estimated_volume_multiplier = 0.15, 18
    tot_estimated_volume = 0
    sum_atomic_masses = 0
    for mol, coeff in zip(molecules, coeffs):
        sum_atomic_volumes = 0
        estimated_volume = 0
        bond_volume_decrease = 0
        for at in mol:
            radius = PeriodicTable.get_radius(at.symbol)
            if PeriodicTable.get_metallic(at.symbol):
                radius *= 0.5
            bond_volume_decrease += (
                coeff * bond_volume_decrease_factor * min(len(at.bonds), 0.3) ** 0.3 * (4 / 3) * np.pi * radius**3
            )
            sum_atomic_volumes += coeff * (4 / 3) * np.pi * radius**3  # ang^3
            sum_atomic_masses += coeff * PeriodicTable.get_mass(at.symbol)

        estimated_volume += sum_atomic_volumes - bond_volume_decrease
        estimated_volume *= estimated_volume_multiplier + 10 / max(4, len(mol))
        tot_estimated_volume += estimated_volume

    return sum_atomic_masses * Units.conversion_ratio("amu", "g") / (tot_estimated_volume * 1e-24)


@overload
def packmol(
    molecules: Union[List[Molecule], Molecule],
    mole_fractions: Optional[List[float]] = ...,
    density: Optional[float] = ...,
    n_atoms: Optional[int] = ...,
    box_bounds: Optional[List[float]] = ...,
    n_molecules: Union[List[int], int, None] = ...,
    sphere: bool = ...,
    fix_first: bool = ...,
    keep_bonds: bool = ...,
    keep_atom_properties: bool = ...,
    region_names: Union[List[str], str, None] = ...,
    return_details: Literal[False] = ...,
    tolerance: float = ...,
    seed: int = ...,
    executable: Optional[str] = ...,
    _return_only_details: bool = ...,
) -> Molecule: ...


@overload
def packmol(
    molecules: Union[List[Molecule], Molecule],
    mole_fractions: Optional[List[float]] = ...,
    density: Optional[float] = ...,
    n_atoms: Optional[int] = ...,
    box_bounds: Optional[List[float]] = ...,
    n_molecules: Union[List[int], int, None] = ...,
    sphere: bool = ...,
    fix_first: bool = ...,
    keep_bonds: bool = ...,
    keep_atom_properties: bool = ...,
    region_names: Union[List[str], str, None] = ...,
    return_details: Literal[True] = ...,
    tolerance: float = ...,
    seed: int = ...,
    executable: Optional[str] = ...,
    _return_only_details: bool = ...,
) -> Tuple[Molecule, Dict[str, Any]]: ...


def packmol(
    molecules: Union[List[Molecule], Molecule],
    mole_fractions: Optional[List[float]] = None,
    density: Optional[float] = None,
    n_atoms: Optional[int] = None,
    box_bounds: Optional[List[float]] = None,
    n_molecules: Union[List[int], int, None] = None,
    sphere: bool = False,
    fix_first: bool = False,
    keep_bonds: bool = True,
    keep_atom_properties: bool = True,
    region_names: Union[List[str], str, None] = None,
    return_details: bool = False,
    tolerance: float = 2.0,
    seed: int = -1,
    executable: Optional[str] = None,
    _return_only_details: bool = False,  # get values of n_molecules, n_atoms, and mole_fractions in the returned dictionary
) -> Union[Molecule, Tuple[Molecule, Dict[str, Any]]]:
    """
    Create a fluid of the given ``molecules``. The function will use the
    given input parameters and try to obtain good values for the others.

    It is *strongly recommended* to specify ``density`` and/or ``box_bounds``. Otherwise you will
    get a (very inaccurate) guessed density in a cubic box (experimental feature).

    molecules : |Molecule| or list of Molecule
        The molecules to pack

    mole_fractions : list of float
        The mole fractions (in the same order as ``molecules``). Cannot be combined with ``n_molecules``. If not given, an equal (molar) mixture of all components will be created.

    density: float
        The total density (in g/cm^3) of the fluid

    n_atoms: int
        The (approximate) number of atoms in the final mixture

    box_bounds: list of float (length 6)
        The box in which to pack the molecules. The box is orthorhombic and should be specified as [xmin, ymin, zmin, xmax, ymax, zmax]. The minimum values should all be set to 0, i.e. set box_bounds=[0., 0., 0., xmax, ymax, zmax]. If not specified, a cubic box of appropriate dimensions will be used.

    n_molecules : int or list of int
        The (exact) number of molecules for each component (in the same order as ``molecules``). Cannot be combined with ``mole_fractions``.

    sphere: bool
        Whether the molecules should be packed in a sphere. The radius is determined by getting the volume from the box bounds!

    fix_first: bool
        Whether to keep the first molecule fixed. This can only be used with ``n_molecules=[1, ..., ...]``. Defaults to False.

    keep_bonds : bool
        If True, the bonds from the constituent molecules will be kept in the returned Molecule

    keep_atom_properties : bool
        If True, the atom.properties (e.g. force-field atom types) of the constituent molecules will be kept in
        the returned Molecule

    region_names : str or list of str
        Populate the region information for each atom. Should have the same length and order as ``molecules``.
        By default the regions are named ``mol0``, ``mol1``, etc.

    tolerance: float
        The packmol tolerance (approximately the minimum intermolecular distance). When packing
        a periodic box, half the tolerance will be excluded from each face of the box.

    return_details : bool
        Return a 2-tuple (Molecule, dict) where the dict has keys like 'n_molecules', 'mole_fractions', 'density',
        etc. They contain the actual details of the returned molecule, which may differ slightly from
        the requested quantities.

        Returned keys:

        * 'n_molecules': list of integer with actually added number of molecules
        * 'mole_fractions': list of float with actually added mole fractions
        * 'density': float, gives the density in g/cm^3
        * 'n_atoms': int, the number of atoms in the returned molecule
        * 'molecule_type_indices': list of int of length n_atoms. For each atom, give an integer index for which TYPE of molecule it belongs to.
        * 'molecule_indices': list of int of length n_atoms. For each atom, give an integer index for which molecule it belongs to
        * 'atom_indices_in_molecule': list of int of length n_atoms. For each atom, give an integer index for which position in the molecule it is.
        * 'volume': float. The volume of the bounding box / packed sphere in ang^3.

    executable : str
        The path to the packmol executable. If not specified, ``$AMSBIN/packmol.exe`` will be used (which is the correct path for the Amsterdam Modeling Suite).

    Useful combinations:

    * ``mole_fractions``, ``density``, ``n_atoms``: Create a mixture with a given density and approximate number of atoms (a cubic box will be created)

    * ``mole_fractions``, ``density``, ``box_bounds``: Create a mixture with a given density inside a given box (the number of molecules will approximately match the density and mole fractions)

    * ``n_molecules``, ``density``: Create a mixture with the given number of molecules and density (a cubic box will be created)

    * ``n_molecules``, ``box_bounds``: Create a mixture with the given number of molecules inside the given box

    Example:

    .. code-block:: python

        packmol(molecules=[from_smiles('O'), from_smiles('C')],
                mole_fractions=[0.8, 0.2],
                density=0.8,
                n_atoms=100)

    Returns: a |Molecule| or tuple (Molecule, dict)
        If return_details=False, return a Molecule. If return_details=True, return a tuple.


    """
    # Input arguments allow for lots of combinations.
    # Let's try to check that the specified combination makes sense ...
    if n_atoms is None and n_molecules is None and density is None:
        raise ValueError("Illegal combination of arguments: must specify either n_atoms, n_molecules or density")
    if n_atoms is not None and n_molecules is not None:
        raise ValueError("Illegal combination of arguments: n_atoms and n_molecules are mutually exclusive")
    if n_atoms is not None and box_bounds is not None and density is not None:
        raise ValueError("Illegal combination of arguments: n_atoms, box_bounds and density specified at the same time")
    if n_molecules is not None and box_bounds is not None and density is not None:
        raise ValueError(
            "Illegal combination of arguments: n_molecules, box_bounds and density specified at the same time"
        )
    if mole_fractions is not None and n_molecules is not None:
        raise ValueError("Illegal combination of arguments: mole_fractions and n_molecules are mutually exclusive")
    if fix_first:
        if n_molecules is None or np.isscalar(n_molecules) or n_molecules[0] != 1:
            raise ValueError(
                f"Illegal combination of arguments: fix_first requires that n_molecules is a list where the first element is 1. Received n_molecules={n_molecules}"
            )
    if isinstance(molecules, list):
        if n_molecules is not None:
            if not isinstance(n_molecules, list):
                raise ValueError("Illegal combination of arguments: molecules is a list, but n_molecules is not")
            if len(n_molecules) != len(molecules):
                raise ValueError("Illegal combination of arguments: len(n_molecules) != len(molecules)")
        if mole_fractions is not None:
            if not isinstance(mole_fractions, list):
                raise ValueError("Illegal combination of arguments: molecules is a list, but mole_fractions is not")
            if len(mole_fractions) != len(molecules):
                raise ValueError("Illegal combination of arguments: len(mole_fractions) != len(molecules)")
        if region_names is not None:
            if not isinstance(region_names, list):
                raise ValueError("Illegal combination of arguments: molecules is a list, but region_names is not")
            if len(region_names) != len(molecules):
                raise ValueError("Illegal combination of arguments: len(region_names) != len(molecules)")
    else:
        if n_molecules is not None and isinstance(n_molecules, list):
            raise ValueError("Illegal combination of arguments: n_molecules is a list, when molecules is not")
        if mole_fractions is not None and isinstance(mole_fractions, list):
            raise ValueError("Illegal combination of arguments: mole_fractions is a list, when molecules is not")
        if region_names is not None and isinstance(region_names, list):
            raise ValueError("Illegal combination of arguments: region_names is a list, when molecules is not")

    if _return_only_details:
        return_details = True

    molecules = tolist(molecules)
    if mole_fractions is None:
        mole_fractions = [1.0 / len(molecules)] * len(molecules)

    if n_molecules:
        n_molecules = tolist(n_molecules)
        sum_n_molecules = np.sum(n_molecules)
        if np.isclose(sum_n_molecules, 0):
            raise ValueError(
                f"The sum of n_molecules is {sum_n_molecules}, which is very close to 0. "
                f"Specify larger numbers. n_molecules specified: {n_molecules}"
            )
        if any(x < 0 for x in n_molecules):
            raise ValueError(f"All n_molecules must be >= 0. " f"n_molecules specified: {n_molecules}")

    xs = np.array(mole_fractions)
    sum_xs = np.sum(xs)
    if np.isclose(sum_xs, 0):
        raise ValueError(
            f"The sum of mole fractions is {sum_xs}, which is very close to 0. "
            f"Specify larger numbers. Mole fractions specified: {xs}"
        )
    if np.any(xs < 0):
        raise ValueError(f"All mole fractions must be >= 0. " f"Mole fractions specified: {mole_fractions}")

    atoms_per_mol = np.array([len(a) for a in molecules])
    masses = np.array([m.get_mass(unit="g") for m in molecules])

    coeffs = None

    if n_molecules:
        coeffs = np.int_(n_molecules)
    elif n_atoms:
        coeff_0 = n_atoms / np.dot(xs, atoms_per_mol)
        coeffs_floats = xs * coeff_0
        coeffs = np.int_(np.round(coeffs_floats))

    if (n_atoms or n_molecules) and not box_bounds:
        mass = np.dot(coeffs, masses)
        if density is not None:
            volume_cm3 = mass / density
        else:
            volume_cm3 = mass / guess_density(molecules, coeffs)
        volume_ang3 = volume_cm3 * 1e24
        side_length = volume_ang3 ** (1 / 3.0)
        box_bounds = [0.0, 0.0, 0.0, side_length, side_length, side_length]
    elif box_bounds and density and not n_molecules:
        volume_cm3 = (
            (box_bounds[3] - box_bounds[0]) * (box_bounds[4] - box_bounds[1]) * (box_bounds[5] - box_bounds[2]) * 1e-24
        )
        mass_g = volume_cm3 * density
        coeffs = mass_g / np.dot(xs, masses)
        coeffs = xs * coeffs
        coeffs = np.int_(np.round(coeffs))

    if coeffs is None:
        raise ValueError(
            f"Illegal combination of arguments: n_atoms={n_atoms}, "
            f"n_molecules={n_molecules}, box_bounds={box_bounds}, density={density}"
        )

    pm = PackMol(executable=executable, tolerance=tolerance, seed=seed)
    if sphere and len(molecules) == 2 and n_molecules and n_molecules[0] == 1:
        # Special case used by packmol_microsolvation
        s1 = PackMolStructure(molecules[0], n_molecules[0], box_bounds=box_bounds, sphere=False, fixed=True)
        s2 = PackMolStructure(molecules[1], n_molecules[1], box_bounds=box_bounds, sphere=True, fixed=False)
        pm.add_structure(s1)
        pm.add_structure(s2)
    else:
        for i, (mol, n_mol) in enumerate(zip(molecules, coeffs)):
            if fix_first and i == 0:
                s1 = PackMolStructure(mol, n_molecules=n_mol, box_bounds=box_bounds, sphere=False, fixed=True)
                pm.add_structure(s1)
            else:
                s1 = PackMolStructure(mol, n_molecules=n_mol, box_bounds=box_bounds, sphere=sphere)
                pm.add_structure(s1)

    if _return_only_details:
        ret = {
            "n_molecules": coeffs.tolist(),
            "mole_fractions": (coeffs / np.sum(coeffs)).tolist() if np.sum(coeffs) > 0 else [0.0] * len(coeffs),
            "n_atoms": np.dot([len(x) for x in molecules], coeffs),
        }
        return None, ret
    out = pm.run()

    # packmol returns the molecules sorted
    molecule_type_indices = []  # [0,0,0,...,1,1,1] # two different molecules with 3 and 5 atoms
    molecule_indices = (
        []
    )  # [0,0,0,1,1,1,2,2,2,....,58,58,58,58,58,59,59,59,59,59] # two different molecules with 3 and 5 atoms
    atom_indices_in_molecule = []  # [0,1,2,0,1,2,...,0,1,2,3,4,0,1,2,3,4]
    current = 0
    for i, (mol, n_mol) in enumerate(zip(molecules, coeffs)):
        molecule_type_indices += [i] * n_mol * len(mol)
        atom_indices_in_molecule += list(range(len(mol))) * n_mol

        temp = list(range(current, current + n_mol))
        molecule_indices += list(np.repeat(temp, len(mol)))
        current += n_mol
    assert len(molecule_type_indices) == len(out)
    assert len(molecule_indices) == len(out)
    assert len(atom_indices_in_molecule) == len(out)

    try:
        volume = out.unit_cell_volume(unit="angstrom")
        density = out.get_density() * 1e-3  # g / cm^3
    except (ValueError, ZeroDivisionError):  # not periodic, presumably when sphere=True
        volume = pm._get_complete_volume()
        if volume == 0:
            density = 0
        else:
            mass = out.get_mass(unit="g")
            density = mass / (volume * 1e-24)  # g / cm^3
    details = {
        "n_molecules": coeffs.tolist(),
        "mole_fractions": (coeffs / np.sum(coeffs)).tolist() if np.sum(coeffs) > 0 else [0.0] * len(coeffs),
        "n_atoms": len(out),
        "molecule_type_indices": molecule_type_indices,  # for each atom, indicate which type of molecule it belongs to by an integer index (starts with 0)
        "molecule_indices": molecule_indices,  # for each atoms, indicate which molecule it belongs to by an integer index (starts with 0)
        "atom_indices_in_molecule": atom_indices_in_molecule,
        "volume": volume,
        "density": density,
    }

    if sphere:
        details["radius"] = np.cbrt(volume * 3 / (4 * np.pi))

    if keep_atom_properties:
        for at, molecule_type_index, atom_index_in_molecule in zip(
            out, molecule_type_indices, atom_indices_in_molecule
        ):
            at.properties = molecules[molecule_type_index][atom_index_in_molecule + 1].properties.copy()

    if keep_bonds:
        out.delete_all_bonds()
        for imol, mol in enumerate(molecules):
            for b in mol.bonds:
                i1, i2 = sorted(mol.index(b))  # 1-based
                for iout, (molecule_type, atom_index_molecule) in enumerate(
                    zip(molecule_type_indices, atom_indices_in_molecule)
                ):
                    if molecule_type != imol:
                        continue
                    if i1 != atom_index_molecule + 1:
                        continue
                    new_i1 = iout + 1  # iout 0-based
                    new_i2 = iout + 1 + i2 - i1  # iout 0-based
                    out.add_bond(out[new_i1], out[new_i2], order=b.order)

    if region_names:
        region_names = tolist(region_names)
    else:
        region_names = [f"mol{i}" for i in range(len(molecules))]

    for at, molindex in zip(out, molecule_type_indices):
        AMSJob._add_region(at, region_names[molindex])

    tot_charge = sum(int(mol.properties.get("charge", 0)) * c for mol, c in zip(molecules, coeffs))
    if tot_charge != 0:
        out.properties.charge = tot_charge

    if return_details:
        return out, details

    return out


def get_packmol_solid_liquid_box_bounds(slab: Molecule):
    slab_max_z = max(at.coords[2] for at in slab)
    slab_min_z = min(at.coords[2] for at in slab)
    liquid_min_z = slab_max_z
    liquid_max_z = liquid_min_z + slab.lattice[2][2] - (slab_max_z - slab_min_z)
    box_bounds = [
        0.0,
        0.0,
        liquid_min_z + 1.5,
        slab.lattice[0][0],
        slab.lattice[1][1],
        liquid_max_z - 1.5,
    ]
    return box_bounds


def packmol_in_void(
    host: Molecule,
    molecules: Union[List[Molecule], Molecule],
    n_molecules: Union[List[int], int],
    keep_bonds: bool = True,
    keep_atom_properties: bool = True,
    region_names: Optional[List[str]] = None,
    tolerance: float = 2.0,
    return_details: bool = False,
    executable: Optional[str] = None,
):
    """
    Pack molecules inside voids in a crystal.

    host: Molecule
        The host molecule. Must be 3D-periodic and the cell must be orthorhombic (all angles 90 degrees) with the lattice vectors parallel to the cartesian axes (all off-diagonal components must be 0).

    For the other arguments, see the ``packmol`` function.

    Note: ``region_names`` needs to have one more element than the list of
    ``molecules``. For example ``region_names=['host', 'guest1',
    'guest2']``.

    """
    if len(host.lattice) != 3:
        raise ValueError("host in packmol_in_void must be 3D periodic")
    if host.cell_angles() != [90.0, 90.0, 90.0]:
        raise ValueError("host in packmol_in_void must be have orthorhombic cell")

    my_host = host.copy()
    my_host.map_to_central_cell(around_origin=False)
    box_bounds = [0.0, 0.0, 0.0, my_host.lattice[0][0], my_host.lattice[1][1], my_host.lattice[2][2]]

    my_molecules = [my_host] + tolist(molecules)
    my_n_molecules = [1] + tolist(n_molecules)

    ret = packmol(
        molecules=my_molecules,
        n_molecules=my_n_molecules,
        box_bounds=box_bounds,
        keep_bonds=keep_bonds,
        keep_atom_properties=keep_atom_properties,
        region_names=region_names,
        return_details=return_details,
        fix_first=True,
        tolerance=tolerance,
        executable=executable,
    )

    return ret


def _run_uff_md(
    ucs: "ChemicalSystem",
    nsteps: int = 1000,
    vectors=None,
    fixed_atoms: Optional[Sequence[int]] = None,
) -> "ChemicalSystem":
    """
    Runs UFF MD with SHAKE all bonds, keeps ``fixed_atoms`` (0-based atom indices) fixed,
    if ``vectors`` is not None will transform into those vectors

    Returns: The final system from the MD simulation.

    Raises: PackmolError if something goes worng.
    """
    from scm.plams import config

    thermostatted_region = "PACKMOL_thermostatted"
    md_ucs = ucs.copy()
    md_ucs.set_atoms_in_region(
        [x for x in range(len(md_ucs)) if fixed_atoms is None or x not in fixed_atoms], thermostatted_region
    )

    s = Settings()
    s.input.ForceField.Type = "UFF"
    s.input.ams.Task = "MolecularDynamics"
    if fixed_atoms:
        s.input.ams.Constraints.AtomList = " ".join(str(x + 1) for x in fixed_atoms)
    s.input.ams.MolecularDynamics.NSteps = nsteps
    s.input.ams.MolecularDynamics.TimeStep = 0.5
    s.input.ams.MolecularDynamics.Shake.All = "bonds * *"
    s.input.ams.MolecularDynamics.InitialVelocities.Type = "Zero"
    # s.input.ams.MolecularDynamics.InitialVelocities.Temperature = 10
    s.input.ams.MolecularDynamics.Thermostat.Temperature = 5
    s.input.ams.MolecularDynamics.Thermostat.Region = thermostatted_region
    s.input.ams.MolecularDynamics.Thermostat.Tau = 2
    s.input.ams.MolecularDynamics.Thermostat.Type = "Berendsen"

    if vectors is not None:
        x = vectors
        target_lattice_str = f"""
            {x[0][0]} {x[0][1]} {x[0][2]}
            {x[1][0]} {x[1][1]} {x[1][2]}
            {x[2][0]} {x[2][1]} {x[2][2]}
        """

        s.input.ams.MolecularDynamics.Deformation.StartStep = 1
        s.input.ams.MolecularDynamics.Deformation.StopStep = (nsteps * 3) // 4
        s.input.ams.MolecularDynamics.Deformation.TargetLattice._1 = target_lattice_str

    previous_config = config.copy()
    try:
        config.job.pickle = False
        config.log.stdout = 0

        with tempfile.TemporaryDirectory() as tmp_dir:
            job_manager = JobManager(config.jobmanager, path=tmp_dir)
            job = AMSJob(settings=s, molecule=md_ucs, name="shakemd")
            job.run(jobmanager=job_manager)

            if not job.ok():
                error_msg = job.results.get_errormsg()
                job_manager._clean()
                raise PackMolError(
                    f"Try a lower density or a less skewed cell! Original file in {job.path}. {error_msg}"
                )

            my_packed = job.results.get_main_system()
            my_packed.remove_region("PACKMOL_thermostatted")
            job_manager._clean()

    finally:
        config.job.pickle = previous_config.job.pickle
        config.log.stdout = previous_config.log.stdout

    return my_packed


@requires_optional_package("scm.libbase")
def packmol_around(
    current: Union[Molecule, "ChemicalSystem"],
    molecules: Union[List[Molecule], Molecule],
    mole_fractions: Optional[List[float]] = None,
    density: Optional[float] = None,
    n_atoms: Optional[int] = None,
    n_molecules: Union[List[int], int, None] = None,
    keep_bonds: bool = True,
    keep_atom_properties: bool = True,
    region_names: Union[List[str], str, None] = None,
    return_details: bool = False,
    tolerance: float = 2.0,
    seed: int = -1,
    executable: Optional[str] = None,
) -> Union[Molecule, Tuple[Molecule, Dict[str, Any]]]:
    """Pack around the current molecule.

    ``current``: Molecule
        Must have a 3D lattice

    ``density``: float
        Density in g/cm^3 of the *added* molecules (excluding ``current``).
        Example: To pack liquid water around a metal slab, set density to 1.0.
        The density is *estimated* from the available free volume and may
        be inaccurate.

    ``mole_fractions``: list of float
        Mole fractions of the *added* molecules.

    ``n_atoms``: float
        Approximate number of *added* molecules.

    ``region_names``: list of str
        Region names for the *added* molecules.

    In general, the arguments refer to the *added* molecules. For all other arguments, see the ``packmol`` function.

    In the returned ``Molecule``, the system will be mapped to ``[0..1]``. It has the same lattice has ``current``.

    .. important::

        The results from this function are almost always approximate! The output
        system will not exactly match your request.

    .. important::

        For non-orthorhombic cells, the results are always approximate. Typically,
        the density will be lower than what you request.

    """
    from scm.libbase import UnifiedChemicalSystem as ChemicalSystem
    from scm.utils.conversions import plams_molecule_to_chemsys, chemsys_to_plams_molecule

    if isinstance(current, Molecule):
        original_ucs = plams_molecule_to_chemsys(current)
    else:
        original_ucs = current.copy()
    assert isinstance(original_ucs, ChemicalSystem)

    if original_ucs.lattice.num_vectors != 3:
        raise ValueError(f"Input molecule `current` must have 3D lattice, got: {current.lattice}")

    lattice_is_orthorhombic = all(
        np.isclose(original_ucs.lattice.vectors[i][j], 0) for i in range(3) for j in range(3) if i != j
    )

    # step 1: store info about original system
    original_ucs.map_atoms(0)
    original_volume = original_ucs.lattice.get_volume()

    # step 2, get remaining volume
    def get_details_for_remaining_volume(original_ucs, molecules, **kwargs):
        sum_r3 = np.sum(np.fromiter((at.element.radius**3 for at in original_ucs), dtype=np.float32))
        current_atomic_volume = (4 / 3) * np.pi * sum_r3
        current_atomic_volume /= 0.74  # use packing efficiency in ccp as example to take up more volume
        remaining_volume = original_volume - current_atomic_volume
        # temporary value to call the original packmol with
        temp_L = np.cbrt(remaining_volume)
        box_bounds_for_remaining_volume = [0.0, 0.0, 0.0, temp_L, temp_L, temp_L]
        # it is unnecessary to actually pack the molecules, this is just used to get the "details"
        # details will contain the correct number of molecules to pack in the combined system
        _, details = packmol(
            molecules=molecules,
            return_details=True,
            box_bounds=box_bounds_for_remaining_volume,
            _return_only_details=True,
            **kwargs,
        )

        details["current_atomic_volume"] = current_atomic_volume
        return details

    molecules = tolist(molecules)
    details = get_details_for_remaining_volume(
        original_ucs,
        molecules,
        mole_fractions=mole_fractions,
        n_atoms=n_atoms,
        n_molecules=n_molecules,
        density=density,
        executable=executable,
    )
    # find cuboid parallel along x/y/z that is guaranteed to encompass the original lattice
    n_molecules = [1] + details["n_molecules"]

    system_for_packing = original_ucs.copy()
    if lattice_is_orthorhombic:
        box_bounds = [0, 0, 0] + np.diag(original_ucs.lattice.vectors).tolist()
    else:
        positive_mask = original_ucs.lattice.vectors >= 0
        negative_mask = original_ucs.lattice.vectors <= 0
        minxyz = np.sum(original_ucs.lattice.vectors * negative_mask, axis=0)
        maxxyz = np.sum(original_ucs.lattice.vectors * positive_mask, axis=0)
        box_bounds = minxyz.tolist() + maxxyz.tolist()
        v_occ = details["current_atomic_volume"]
        v_orig_free = original_ucs.lattice.get_volume() - v_occ
        v_new_free = np.prod(maxxyz - minxyz) - v_occ
        volume_multiplier = v_new_free / v_orig_free
        new_n_atoms = int(np.round(details["n_atoms"] * volume_multiplier))
        new_details = get_details_for_remaining_volume(
            system_for_packing,
            molecules,
            n_atoms=new_n_atoms,
            mole_fractions=details["mole_fractions"],
            executable=executable,
        )
        n_molecules = [1] + new_details["n_molecules"]

    my_molecules = [chemsys_to_plams_molecule(system_for_packing)] + tolist(molecules)
    if region_names is not None:
        region_names = tolist(region_names)
        if len(region_names) == len(my_molecules) - 1:
            # insert a dummy region name, it will not be returned anyway
            region_names = ["current"] + region_names
    my_packed, details = packmol(
        molecules=my_molecules,
        n_molecules=n_molecules,
        fix_first=True,
        box_bounds=box_bounds,
        return_details=True,
        tolerance=tolerance,
        executable=executable,
        keep_bonds=keep_bonds,
        keep_atom_properties=keep_atom_properties,
        region_names=region_names,
        seed=seed,
    )

    # remove the original substrate
    my_packed = plams_molecule_to_chemsys(my_packed)
    my_packed.remove_atoms(range(len(original_ucs)))

    ### start removing molecules outside the unit cell for non-orthorhombic cells
    ### for orthorhombic cells the fractional coordiantes are all in (0,1)
    ### so nothing will happen here
    # packed_n_molecules = my_packed.num_molecules()
    my_packed.lattice = original_ucs.lattice.copy()
    fractional_coords = my_packed.get_fractional_coordinates()
    # TODO: work out some reasonable margin depending on lattice vector lengths
    mask = (fractional_coords < 0) | (fractional_coords >= 1)
    row_indices = np.any(mask, axis=1).nonzero()[0]
    my_packed.select_atoms(row_indices)
    # TODO: select_molecule assumes that there are bonds defined - perhaps
    # one should really go through the details dictionary to find the indices
    # of added molecules regardless of connectivity
    my_packed.select_molecule()
    removed_atoms_from_my_packed = my_packed.get_selected_atoms()
    my_packed.remove_atoms(my_packed.get_selected_atoms())
    # removed_n_molecules = packed_n_molecules - my_packed.num_molecules()

    mti = details["molecule_type_indices"]
    removed_molecules_types = [mti[i + len(original_ucs)] for i in removed_atoms_from_my_packed]

    counter = Counter(removed_molecules_types)
    ret_details = dict(n_molecules=[])
    for imol, nmol in enumerate(details["n_molecules"]):
        if imol == 0:  # skip "current"
            continue
        lenmol = len(my_molecules[imol])
        if lenmol == 0:
            new_nmol = nmol
        else:
            new_nmol = nmol - counter[imol] // lenmol
        ret_details["n_molecules"].append(new_nmol)
    sum_molecules = np.sum(ret_details["n_molecules"])
    if sum_molecules == 0:
        ret_details["mole_fractions"] = [0] * len(ret_details["n_molecules"])
    else:
        ret_details["mole_fractions"] = (np.array(ret_details["n_molecules"]) / sum_molecules).tolist()

    ret = original_ucs + my_packed
    ret = chemsys_to_plams_molecule(ret)
    if return_details:
        return ret, ret_details
    return ret


def packmol_on_slab(
    slab: Molecule,
    molecules: Union[List[Molecule], Molecule],
    mole_fractions: Optional[List[float]] = None,
    density: float = 1.0,
    keep_bonds: bool = True,
    keep_atom_properties: bool = True,
    region_names: Optional[List[str]] = None,
    executable: Optional[str] = None,
):
    """

    Creates a solid/liquid interface with an approximately correct density. The
    density is calculated for the volume not occupied by the slab (+ 1.5
    angstrom buffer at each side of the slab).

    Returns: a |Molecule|

    slab : |Molecule|
        The system must have a 3D lattice (including a vacuum gap along z) and be orthorhombic. The vacuum gap will be filled with the liquid.

    For the other arguments, see ``packmol``.

    Example:

    .. code-block:: python

        packmol_on_slab(slab=slab_3d_with_vacuum_gap,
                        molecules=[from_smiles('O'), from_smiles('C')],
                        mole_fractions=[0.8, 0.2],
                        density=0.8)

    """
    if len(slab.lattice) != 3:
        raise ValueError("slab in packmol_on_slab must be 3D periodic: slab in xy-plane with vacuum gap along z-axis")
    if not all(np.isclose(slab.cell_angles(), [90.0, 90.0, 90.0])):
        raise ValueError("slab in packmol_on_slab must be have orthorhombic cell")

    liquid = packmol(
        molecules=molecules,
        mole_fractions=mole_fractions,
        density=density,
        box_bounds=get_packmol_solid_liquid_box_bounds(slab),
        keep_bonds=keep_bonds,
        keep_atom_properties=keep_atom_properties,
        region_names=region_names,
        executable=executable,
    )

    # Map all liquid molecules to [0..1]
    # NOTE: We need to be using the lattice of the slab for this!
    #       The lattice of the liquid is different ...
    liquid.lattice = slab.lattice
    # If the slab has cell-shifts for the bonds, the liquid also needs to have
    # them. If would not have cell-shifts, they would not be updated in the
    # map_to_central_cell call, even though they would become significant when
    # combining with the slab that has them: minimum image convention is only
    # assumed if no bond has cell-shifts.
    if liquid.bonds and any(b.has_cell_shifts() for b in slab.bonds):
        for b in liquid.bonds:
            b.properties.suffix = "0 0 0"
    liquid.map_to_central_cell(around_origin=False)
    if liquid.bonds and any(b.has_cell_shifts() for b in slab.bonds):
        for b in liquid.bonds:
            if b.properties.suffix == "0 0 0":
                del b.properties.suffix

    # Shift liquid molecules (now in [0..1]) on top of the slab.
    # The slab could be anywhere, e.g. [-0.5..0.5] ...
    slab_center_x = (max(at.coords[0] for at in slab) + min(at.coords[0] for at in slab)) / 2
    slab_center_y = (max(at.coords[1] for at in slab) + min(at.coords[1] for at in slab)) / 2
    liquid_center_x = (max(at.coords[0] for at in liquid) + min(at.coords[0] for at in liquid)) / 2
    liquid_center_y = (max(at.coords[1] for at in liquid) + min(at.coords[1] for at in liquid)) / 2
    liquid.translate([-liquid_center_x + slab_center_x, -liquid_center_y + slab_center_y, 0.0])

    out = slab.copy()
    for at in out:
        AMSJob._add_region(at, "slab")
    out.add_molecule(liquid)
    return out


def get_n_from_density_and_box_bounds(molecule, box_bounds, density):
    molecule_mass = molecule.get_mass(unit="g")
    volume_ang3 = (box_bounds[3] - box_bounds[0]) * (box_bounds[4] - box_bounds[1]) * (box_bounds[5] - box_bounds[2])
    volume_cm3 = volume_ang3 * 1e-24
    n_molecules = int(density * volume_cm3 / molecule_mass)
    return n_molecules


def packmol_microsolvation(
    solute: Molecule,
    solvent: Molecule,
    density: float = 1.0,
    threshold: float = 3.0,
    keep_bonds: bool = True,
    keep_atom_properties: bool = True,
    region_names: List[str] = ["solute", "solvent"],
    executable: Optional[str] = None,
):
    """
    Microsolvation of a ``solute`` with a ``solvent`` with an approximate ``density``.

    solute: |Molecule|
        The solute to be surrounded by solvent molecules

    solvent: |Molecule|
        The solvent molecule

    density: float
        Approximate density in g/cm^3

    threshold: float
        Distance in angstrom. Any solvent molecule for which at least 1 atom is within this threshold to the solute molecule will be kept

    For the other arguments, see ``packmol``.

    """

    solute_coords = solute.as_array()
    com = np.mean(solute_coords, axis=0)
    plams_solute = solute.copy()
    plams_solute.translate(-com)
    solute_coords = plams_solute.as_array()
    box_bounds = [0.0, 0.0, 0.0] + list(np.max(solute_coords, axis=0) - np.min(solute_coords, axis=0) + 3 * threshold)

    n_solvent = get_n_from_density_and_box_bounds(solvent, box_bounds, density=density)

    plams_solvated = packmol(
        [plams_solute, solvent],
        n_molecules=[1, n_solvent],
        box_bounds=box_bounds,
        keep_bonds=keep_bonds,
        keep_atom_properties=keep_atom_properties,
        region_names=region_names,
        sphere=True,
        executable=executable,
    )

    solute_indices = [i for i, at in enumerate(plams_solvated, 1) if i <= len(solute)]
    newmolecule = plams_solvated.get_complete_molecules_within_threshold(solute_indices, threshold=threshold)

    return newmolecule


@requires_optional_package("scm.libbase")
def packmol_around_md(
    current: Union[Molecule, "ChemicalSystem"],
    molecules: Union[Molecule, List[Molecule]],
    return_details: bool = False,
    always_run_md: bool = False,
    **kwargs,
) -> Molecule:
    """Pack around the current molecule, relax structure with MD.

    Experimental feature.

    ``current``: Molecule
        Must have a 3D lattice

    ``always_run_md``: bool
        If True, will run UFF MD also for orthorhombic cells. For nonorthorhombic cells, MD is always run irrespective of this flag.

    For all other arguments, see the ``packmol`` function.

    In the returned ``Molecule``, the system will be mapped to ``[0..1]``. It has the same lattice has ``current``.
    """
    from scm.libbase import (
        UnifiedChemicalSystem as ChemicalSystem,
        UnifiedLattice as Lattice,
    )
    from scm.utils.conversions import plams_molecule_to_chemsys, chemsys_to_plams_molecule

    loglevel = 7

    if isinstance(current, Molecule):
        original_ucs = plams_molecule_to_chemsys(current)
    else:
        original_ucs = current.copy()
    assert isinstance(original_ucs, ChemicalSystem)
    original_ucs.map_atoms(0)

    # step 1: store info about original system
    if original_ucs.lattice.num_vectors != 3:
        raise ValueError(f"Input molecule `current` must have 3D lattice, got: {current.lattice}")
    original_frac_coords = original_ucs.get_fractional_coordinates()

    original_volume = original_ucs.lattice.get_volume()
    original_lattice = original_ucs.lattice.copy()

    # step 2, get remaining volume
    current_atomic_volume = (
        (4 / 3) * np.pi * np.sum(np.fromiter((at.element.radius for at in original_ucs), dtype=np.float32))
    )
    current_atomic_volume /= 0.74  # use packing efficiency in ccp as example to take up more volume
    remaining_volume = original_volume - current_atomic_volume
    # temporary value to call the original packmol with
    temp_L = np.cbrt(remaining_volume)
    box_bounds_for_remaining_volume = [0.0, 0.0, 0.0, temp_L, temp_L, temp_L]
    # it is unnecessary to actually pack the molecules, this is just used to get the "details"
    # details will contain the correct number of molecules to pack in the combined system
    # TODO: reorganize the packmol function so that one can get this info without calling packmol
    log(f"Initial packing to determine number of molecules: {molecules}, {box_bounds_for_remaining_volume}", loglevel)
    _, details = packmol(molecules=molecules, return_details=True, box_bounds=box_bounds_for_remaining_volume, **kwargs)

    # find cuboid parallel along x/y/z that is guaranteed to encompass the original lattice
    maxcomponents = np.max(original_ucs.lattice.vectors, axis=0) - np.min(original_ucs.lattice.vectors, axis=0)
    box_bounds = [0.0, 0.0, 0.0] + list(maxcomponents)

    will_run_uff_md = always_run_md or any(
        (not np.isclose(original_lattice.vectors[i][j], 0) for i in range(3) for j in range(3) if i != j)
    )
    log(f"will_run_uff_md: {will_run_uff_md}", loglevel)
    if will_run_uff_md:
        if np.linalg.det(original_lattice.vectors) < 0:
            raise PackMolError("packmol_around cannot handle lattice where the determinant of the vectors is negative.")

    target_lattice = np.diag(maxcomponents)

    system_for_packing_type = "supercell"  # "supercell" or "distorted"
    # system_for_packing_type = "distorted"  # "supercell" or "distorted"

    n_molecules = [1] + details["n_molecules"]
    if system_for_packing_type == "supercell":
        # Create a supercell that should encompass the target x/y/z lattice
        # this is used for the initial packing of molecules to ensure there is no overlap
        # with the original atoms
        trafo = np.linalg.inv(original_ucs.lattice.vectors) @ np.array(target_lattice)
        trafo = np.sign(trafo) * np.ceil(np.abs(trafo))
        trafo = np.int_(trafo)
        supercell = original_ucs.make_supercell_trafo(trafo)
        supercell.map_atoms(0)
        system_for_packing = supercell
        tolerance = kwargs.get("tolerance", 1.5)
    else:
        # now distort the original system to the target lattice
        distorted = original_ucs.copy()
        distorted.lattice.vectors = np.diag(maxcomponents)
        distorted.set_fractional_coordinates(original_frac_coords)
        system_for_packing = distorted
        # in general we need higher tolerance here since we may be expanding the original system,
        # and we do not want the added molecules to enter in artificial "voids"
        tolerance = kwargs.get("tolerance", 1.5) * 1.3  # should depend on distortion_vol_expansion_factor somehow

    log(f"{system_for_packing_type=}", loglevel)
    log(f"{n_molecules=}, {box_bounds=}, {tolerance=}", loglevel)
    my_packed, details = packmol(
        molecules=[chemsys_to_plams_molecule(system_for_packing)] + tolist(molecules),
        n_molecules=n_molecules,
        fix_first=True,
        box_bounds=box_bounds,
        return_details=True,
        tolerance=tolerance,
    )
    # remove the original substrate
    my_packed = plams_molecule_to_chemsys(my_packed)
    my_packed.remove_atoms(range(len(system_for_packing)))
    my_packed.map_atoms_continuous()
    my_packed.lattice = Lattice()  # so that we can add_other without having incompatible lattices

    # now create a distorted system
    distorted = original_ucs.copy()
    distorted.lattice.vectors = np.diag(maxcomponents)
    distorted.set_fractional_coordinates(original_frac_coords)

    # distortion_vol_expansion_factor will be used to modify the tolerance when doing the packing
    # distortion_vol_expansion_factor = (distorted.lattice.get_volume() / original_volume) ** (1 / 3.0)
    # remove bonds to be able to do "shake all bonds * *" for the remaining molecules
    if will_run_uff_md:
        distorted.bonds.clear_bonds()
    distorted.add_other(my_packed)

    if will_run_uff_md:
        log("Running UFF MD", loglevel)
        distorted = _run_uff_md(
            distorted,
            nsteps=1500,
            vectors=original_ucs.lattice.vectors,
            fixed_atoms=list(range(len(original_ucs))),
        )

    distorted.remove_atoms(range(len(original_ucs)))
    distorted.map_atoms_continuous()
    distorted.lattice = Lattice()  # so that we can add_other without having incompatible lattices

    # this ensures that the original UCS is exactly preserved (including bonds etc.)
    out_ucs = original_ucs.copy()
    out_ucs.add_other(distorted)
    out_ucs.map_atoms(0)

    out_mol = chemsys_to_plams_molecule(out_ucs)

    if return_details:
        return out_mol, details
    return out_mol

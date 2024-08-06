import os
import subprocess
import tempfile
from typing import Any, Dict, List, Literal, Optional, Tuple, Union, overload

import numpy as np
from scm.plams.core.errors import MoleculeError
from scm.plams.core.private import saferun
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.mol.molecule import Molecule

try:
    from scm.plams.interfaces.molecule.rdkit import readpdb, writepdb
except ImportError:
    pass

__all__ = [
    "packmol",
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
                self.box_bounds = [
                    0.0,
                    0.0,
                    0.0,
                    molecule.lattice[0][0],
                    molecule.lattice[1][1],
                    molecule.lattice[2][2],
                ]
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
                self.n_molecules = self._get_n_molecules(self.molecule, n_atoms)
            else:
                self.n_molecules = n_molecules
            assert box_bounds or density
            self.box_bounds = box_bounds or self._get_box_bounds(self.molecule, self.n_molecules, density)
            self.fixed = False
            self.sphere = sphere

    def _get_n_molecules_from_density_and_box_bounds(self, molecule: Molecule, box_bounds: List[float], density: float):
        """density in g/cm^3"""
        molecule_mass = molecule.get_mass(unit="g")
        volume_ang3 = self.get_volume(box_bounds)
        volume_cm3 = volume_ang3 * 1e-24
        n_molecules = int(density * volume_cm3 / molecule_mass)
        return n_molecules

    def get_volume(self, box_bounds=None):
        bb = box_bounds or self.box_bounds
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
            # vol = 4*pi*r^3 /3
            # radius = (3*vol/(4*pi))**0.33333
            radius = (3 * vol / (4 * 3.14159)) ** 0.3333
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

        Note: users are not recommended to use this class directly, but
        instead use the ``packmol``, ``packmol_on_slab`` and ``packmol_microsolvation``
        functions.

        """
        self.tolerance = tolerance
        self.structures = structures or []
        self.filetype = filetype
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
        radius = (3 * volume / (4 * 3.14159)) ** 0.3333

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

            my_input = open(input_fname, "r")
            saferun(self.executable, stdin=my_input, stdout=subprocess.DEVNULL)
            my_input.close()

            if not os.path.exists(output_fname):
                raise PackMolError("Packmol failed. It may work if you try a lower density.")

            if self.filetype == "pdb":
                with open(output_fname, "r") as f:
                    output_molecule = readpdb(f)
            else:
                output_molecule = Molecule(output_fname)  # without periodic boundary conditions

            output_molecule.lattice = self._get_complete_lattice()

        return output_molecule


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
    executable: Optional[str] = ...,
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
    executable: Optional[str] = ...,
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
    executable: Optional[str] = None,
) -> Union[Molecule, Tuple[Molecule, Dict[str, Any]]]:
    """
    Create a fluid of the given ``molecules``. The function will use the
    given input parameters and try to obtain good values for the others. You *must*
    specify ``density`` and/or ``box_bounds``.

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
    if density is None and box_bounds is None:
        raise ValueError("Illegal combination of arguments: must specify either density or box_bounds")
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

    if (n_atoms or n_molecules) and density and not box_bounds:
        mass = np.dot(coeffs, masses)
        volume_cm3 = mass / density
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
            f"Illegal combination of arguments: n_atoms={n_atoms}, n_molecules={n_molecules}, box_bounds={box_bounds}, density={density}"
        )

    pm = PackMol(executable=executable, tolerance=tolerance)
    if sphere and len(molecules) == 2 and n_molecules and n_molecules[0] == 1:
        # Special case used by packmol_microsolvation
        pm.add_structure(
            PackMolStructure(molecules[0], n_molecules[0], box_bounds=box_bounds, sphere=False, fixed=True)
        )
        pm.add_structure(
            PackMolStructure(molecules[1], n_molecules[1], box_bounds=box_bounds, sphere=True, fixed=False)
        )
    else:
        for i, (mol, n_mol) in enumerate(zip(molecules, coeffs)):
            if fix_first and i == 0:
                pm.add_structure(
                    PackMolStructure(mol, n_molecules=n_mol, box_bounds=box_bounds, sphere=False, fixed=True)
                )
            else:
                pm.add_structure(PackMolStructure(mol, n_molecules=n_mol, box_bounds=box_bounds, sphere=sphere))

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
    except ValueError:  # not periodic, presumably when sphere=True
        volume = pm._get_complete_volume()
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
        details["radius"] = (volume * 3 / (4 * 3.1415926535)) ** (1 / 3.0)

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
    # them. If if would not have cell-shifts, they would not be updated in the
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

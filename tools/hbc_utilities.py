from scm.plams.core.functions import requires_optional_package
from scm.plams.version import __version__
from typing import Optional, List, Tuple, Dict, Any
from pathlib import Path
import numpy as np
import os


@requires_optional_package("scm.libbase")
def view_HBC(rkf_path: str, xyz_file: Optional[str] = None):
    """
    visulize the hydrogen bond centers (HBC)

    Args:
        rkf_path (str) : A COSKF file used to visulize the HBC.
        xyz_file (optional str) : The name of the XYZ file to write the visualization to. If not provided, a temporary file "tmp-HBC.xyz" will be used.
    """
    from scm.libbase import KFFile

    if not isinstance(rkf_path, str):
        raise TypeError(f"Expected `rkf_path` to be a string but got {type(rkf_path).__name__}")
    elif not Path(rkf_path).exists():
        raise FileNotFoundError(rkf_path)

    xyz_filename = xyz_file if xyz_file is not None else "tmp-HBC.xyz"

    with KFFile(rkf_path) as rkf:
        natom = rkf.read("COSMO", "Number of Atoms")
        atom_symbols = rkf.read("COSMO", "Atom Type").split()
        atom_coords = np.array(rkf.read("COSMO", "Atom Coordinates"))
        atom_coords = np.reshape(atom_coords, (natom, 3))

        n_HBC = rkf.read("HBC", "Number of HBCs")

        total_number = natom + n_HBC
        XYZ = f"{total_number}\n\n"

        atom_id = 0
        for s, (x, y, z) in zip(atom_symbols, atom_coords):
            XYZ += f"{s} {x:.8f} {y:.8f} {z:.8f}\n"
            atom_id += 1

        if n_HBC != 0:
            coords = np.array(rkf.read("HBC", "HBC Coordinates"))
            coords = np.reshape(coords, (3, n_HBC))
            coords = np.transpose(coords)

            for x, y, z in coords:
                XYZ += f"Xx {x:.8f} {y:.8f} {z:.8f}\n"
                atom_id += 1

    with open(xyz_filename, "w") as file:
        file.write(XYZ)
    os.system(f"$AMSBIN/amsview {xyz_filename}")


@requires_optional_package("scm.libbase")
def write_HBC_to_COSKF(
    rkf_path: str, HBC_xyz: List[np.ndarray], HBC_atom: List[int], HBC_angle: List[float], HBC_info: Dict[str, Any]
):
    """
    Write the hydrogen bond centers (HBC) information into the COSKF file
    """
    from scm.libbase import KFFile

    if not isinstance(rkf_path, str):
        raise TypeError(f"Expected `rkf_path` to be a string but got {type(rkf_path).__name__}")
    elif not Path(rkf_path).exists():
        raise FileNotFoundError(rkf_path)

    with KFFile(rkf_path) as rkf:
        n_HBC = len(HBC_atom)
        for key, value in HBC_info.items():
            rkf.write("HBC", key, value)
        rkf.write("HBC", "Number of HBCs", n_HBC)
        if n_HBC != 0:
            rkf.write("HBC", "HBC Coordinates", np.transpose(HBC_xyz).flatten())
            rkf.write("HBC", "HBC Atom", HBC_atom)
            rkf.write("HBC", "HBC angles", HBC_angle)


@requires_optional_package("scm.libbase")
def parse_mesp(
    densf_path: str,
    rkf_path: str,
    skip_mesp_atoms: Optional[List[str]] = None,
    allowd_HB_atoms: Optional[List[str]] = None,
) -> Tuple[List[np.ndarray], List[int], List[float], Dict[str, Any]]:
    """
    Extract local minima of the molecular electrostatic potential (MEP) from the Densf output (TAPE41 file),
    convert them to hydrogen bond centers (HBC)

    Args:
        densf_path (str): The path to a TAPE41 file containing the molecular electrostatic potential data from Densf calculation.
        rkf_path (str): The path to a COSKF file where the hydrogen bond center data will be written.
        skip_mesp_atoms (Optional[List[str]]): A list of atom types to exclude when preparing HBCs from MESP. Defaults to ["C", "H"].
        allowd_HB_atoms (Optional[List[str]]): A list of atom types used to determine the HBC when a hydrogen atom is bonded to them. The HBC of a hydrogen atom is calculated using the bond vector. Defaults to ["O", "N", "F"].

    Returns:
        HBC_xyz (List[np.ndarray]): A list of numpy arrays representing the coordinates of the hydrogen bond centers.
        HBC_Atom (List[int]): A list of integers representing the atom indices of the hydrogen bond centers.
        HBC_angle (List[float]): A list of floats representing the angles between the position vector and the eigenvector.
        HBC_info (Dict[str, Any]) : A dictionary containing metadata, including the ADF version, density grid type, and HBC script version.

    """
    from scm.libbase import KFFile, Units

    rkf = KFFile(rkf_path)
    densf = KFFile(densf_path)

    adf_version = f"{__version__}"
    Grid_type = densf.read_string("Grid", "grid type").strip()
    HBC_version = "MESP.2025.v1"
    HBC_info = {"adf version": adf_version, "densf grid type": Grid_type, "hbc version": HBC_version}

    if skip_mesp_atoms is None:
        skip_mesp_atoms = ["C", "H"]
    if allowd_HB_atoms is None:
        allowd_HB_atoms = ["O", "N", "F"]

    print(
        "Prepare the Hydrogen Bond Center (HBC) from the local minima of the molecular electrostatic potential (MESP).\n"
        f" *Excluding the following atoms: {', '.join(skip_mesp_atoms)}.\n"
        f" *The HBC of hydrogen bonded to {', '.join(allowd_HB_atoms)} is determined using the bond vector."
    )
    HBC_xyz: List[np.ndarray] = []
    HBC_atom: List[int] = []
    HBC_angle: List[float] = []
    BOHR = Units.convert("bohr", "angstrom", 1.0)

    atom_COSMO_radius = None
    if rkf.var_exists("COSMO", "Atom COSMO Radii"):
        atom_COSMO_radius = rkf.read_reals_np("COSMO", "Atom COSMO Radii") * BOHR

    COSMO_radius_default = {
        "H": 1.30,
        "C": 2.00,
        "N": 1.83,
        "O": 1.72,
        "F": 1.72,
        "Si": 2.48,
        "P": 2.13,
        "S": 2.16,
        "Cl": 2.05,
        "Br": 2.16,
        "I": 2.32,
    }

    if not densf.section_exists("Coulpot minima"):
        return HBC_xyz, HBC_atom, HBC_angle, HBC_info

    n_MESP = densf.read_int("Coulpot minima", "Number of minima")

    if n_MESP == 0:
        return HBC_xyz, HBC_atom, HBC_angle, HBC_info

    natom = rkf.read_int("COSMO", "Number of Atoms")
    atom_symbol = rkf.read_string("COSMO", "Atom Type").split()
    atom_coord = np.reshape(rkf.read_reals_np("COSMO", "Atom Coordinates"), (natom, 3))

    # Determine HBC through the projection of the MESP onto COSMO Cavity
    for i in range(n_MESP):
        Eigenvalues = densf.read_reals_np("Coulpot minima", f"Eigenvalues {i+1}")
        Eigenvectors = np.reshape(densf.read_reals_np("Coulpot minima", f"Eigenvectors {i+1}"), (3, 3))
        EV = Eigenvectors[np.argmax(Eigenvalues)]

        MESP_coord = densf.read_reals_np("Coulpot minima", f"Coords {i+1}") * BOHR

        sum_of_squares = np.sum(np.power(atom_coord - MESP_coord, 2.0), axis=1)
        MESP_atom_id = np.argmin(sum_of_squares)

        if atom_symbol[MESP_atom_id] in skip_mesp_atoms:
            continue

        PV = (MESP_coord - atom_coord[MESP_atom_id]) / np.sqrt(sum_of_squares[MESP_atom_id])

        if atom_COSMO_radius is not None:
            HB_center = atom_coord[MESP_atom_id] + PV * atom_COSMO_radius[MESP_atom_id]
        else:
            HB_center = atom_coord[MESP_atom_id] + PV * COSMO_radius_default[atom_symbol[MESP_atom_id]]

        HBC_xyz.append(HB_center)
        HBC_atom.append(int(MESP_atom_id + 1))

        dot_product = np.dot(PV, EV)
        cos_theta = np.clip(dot_product, -1.0, 1.0)
        angle_radians = np.arccos(cos_theta)
        angle_degrees = np.degrees(angle_radians)
        if angle_degrees > 90:
            angle_degrees = 180 - angle_degrees
        HBC_angle.append(angle_degrees)

    # determine HBC for Hydrogen atom bounded to allowd_HB_atoms
    for idx, atom in enumerate(atom_symbol):
        if atom == "H":
            cur_coord = atom_coord[idx]
            sum_of_squares = np.sum(np.power(atom_coord - cur_coord, 2.0), axis=1)
            non_zero_indices = np.where(sum_of_squares != 0)[0]
            smallest_non_zero_index = non_zero_indices[np.argmin(sum_of_squares[non_zero_indices])]
            if atom_symbol[smallest_non_zero_index] in allowd_HB_atoms:
                PV = (cur_coord - atom_coord[smallest_non_zero_index]) / np.sqrt(
                    sum_of_squares[smallest_non_zero_index]
                )
                HB_center = cur_coord + PV * COSMO_radius_default["H"]
                HBC_xyz.append(HB_center)
                HBC_atom.append(idx + 1)
                HBC_angle.append(0.0)

    densf.close()
    rkf.close()

    return HBC_xyz, HBC_atom, HBC_angle, HBC_info

#!/usr/bin/env amspython
# coding: utf-8

# ## Complete guide to storing and converting PLAMS Molecules between Python libraries and file formats

import os
from os.path import expandvars
from pathlib import Path

# Make sure to source amsbashrc.sh before launching this example so that
# the AMSHOME environment variable is set. Requires AMS2025+ to run this example.

AMSHOME = os.environ["AMSHOME"]
cif_file = f"{AMSHOME}/atomicdata/Molecules/IZA-Zeolites/ABW.cif"
xyz_file = f"{AMSHOME}/scripting/scm/params/examples/benchmark/ISOL6/e_13.xyz"
badxyz_file = f"{AMSHOME}/scripting/scm/plams/unit_tests/xyz/reactant2.xyz"

assert Path(cif_file).exists(), f"{cif_file} does not exist."
assert Path(xyz_file).exists(), f"{xyz_file} does not exist."


def head(filename, n: int = 4):
    """Print the first ``n`` lines of a file"""
    with open(filename, "r") as f:
        lines = f.readlines()
        lines = lines[: min(n, len(lines))]
    print("".join(lines))


# ### SMILES

# #### Load PLAMS Molecule from SMILES string

from scm.plams import from_smiles, Molecule, plot_molecule

mol = from_smiles("CCCCO")
print(f"{type(mol)=}")
plot_molecule(mol)


# #### Convert PLAMS Molecule to SMILES string
#
# Note: This requires that bonds are defined in the PLAMS Molecule.

from scm.plams import to_smiles

smiles = to_smiles(mol)
print(smiles)


# ### .xyz files

# #### Load PLAMS Molecule from .xyz file

from scm.plams import Molecule, plot_molecule

mol = Molecule(xyz_file)
print(f"{type(mol)=}")
plot_molecule(mol)


# #### Write PLAMS Molecule to .xyz file

mol.properties.comment = "The comment line (2nd line after the number of atoms)"
mol.write("out.xyz")


head("out.xyz")


# ### .cif files

# #### Load PLAMS Molecule from .cif file
#
# PLAMS cannot natively read .cif files. Instead, go through another library, for example ASE or pymatgen.

from ase.io import read
from scm.plams import fromASE

mol: Molecule = fromASE(read(cif_file))
print(f"{type(mol)=}")
plot_molecule(mol)


# #### Write PLAMS Molecule to .cif file
#
# PLAMS cannot natively export to .cif files. Instead, go through another library, for example ASE or pymatgen.
#
# ASE can be used to write many file formats. See https://wiki.fysik.dtu.dk/ase/ase/io/io.html

from scm.plams import toASE

toASE(mol).write("out.cif")
head("out.cif")


# ### AMS .in system block format
#
# #### Write PLAMS Molecule to AMS .in system file

mol.write("ams_system_block.in")
head("ams_system_block.in")


# #### Load PLAMS Molecule from AMS .in system file

from scm.plams import Molecule

mol = Molecule("ams_system_block.in")
plot_molecule(mol)


# ### POSCAR/CONTCAR (VASP input format)

# #### Write PLAMS Molecule to POSCAR/CONTCAR (VASP input format)
#
# ASE can be used to write many file formats. See https://wiki.fysik.dtu.dk/ase/ase/io/io.html

from scm.plams import toASE

toASE(mol).write("POSCAR")
head("POSCAR", 10)


# #### Load PLAMS Molecule from POSCAR/CONTCAR (VASP input format)

from scm.plams import fromASE
from ase.io import read

mol: Molecule = fromASE(read("POSCAR"))

print(f"{type(mol)=}")
plot_molecule(mol)


# ### ASE Atoms Python class

# #### Convert PLAMS Molecule to ASE Atoms

from scm.plams import toASE
from ase import Atoms
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

print(f"{type(mol)=}")
print(f"{mol.get_formula()=}")

ase_atoms: Atoms = toASE(mol)
print(f"{type(ase_atoms)=}")
print(f"{ase_atoms.get_chemical_formula()=}")

_, ax = plt.subplots(figsize=(2, 2))
plot_atoms(ase_atoms, rotation="-85x,5y,0z", ax=ax)


# #### Convert ASE Atoms to PLAMS Molecule

from scm.plams import fromASE, plot_molecule, Molecule

mol: Molecule = fromASE(ase_atoms)
print(f"{type(mol)=}")
plot_molecule(mol, rotation="-85x,5y,0z")


# ### RDKit Mol Python class

# #### Convert PLAMS Molecule to RDKit Mol

from scm.plams import to_rdmol, Molecule
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.ipython_useSVG = True
IPythonConsole.molSize = 250, 250

plams_mol = Molecule(xyz_file)
# guess bonds, the bonds will be included in the RDKit molecule
plams_mol.guess_bonds()

rdkit_mol = to_rdmol(plams_mol)
print(f"{type(rdkit_mol)=}")
rdkit_mol


# #### Convert RDKit Mol to PLAMS Molecule

from scm.plams import from_rdmol, plot_molecule, Molecule

mol: Molecule = from_rdmol(rdkit_mol)

print(f"{type(rdkit_mol)=}")
print(f"{type(mol)=}")
plot_molecule(mol)


# #### Convert problematic PLAMS Molecule to RDKit Mol

mol = Molecule(badxyz_file)
mol.guess_bonds()
plot_molecule(mol)


# This molecule will fail to convert to an RDKit Mol object, because RDKit does not like the AMS assignment of double bonds.

try:
    rdkit_mol = to_rdmol(mol)
except ValueError as exc:
    print("Failed to convert")


# The problem can be fixed by passing the argument `presanitize` to the `to_rdmol` function.

rdkit_mol = to_rdmol(mol, presanitize=True)
rdkit_mol


# ### SCM libbase UnifiedChemicalSystem Python class
#
# #### Convert PLAMS Molecule to UnifiedChemicalSystem

from scm.utils.conversions import plams_molecule_to_chemsys, chemsys_to_plams_molecule
from scm.plams import Molecule
from scm.libbase import UnifiedChemicalSystem

mol = Molecule(xyz_file)
chemsys = plams_molecule_to_chemsys(mol)
print(f"{type(chemsys)=}")
print(chemsys)


# #### Convert UnifiedChemicalSystem to PLAMS Molecule

from scm.utils.conversions import plams_molecule_to_chemsys, chemsys_to_plams_molecule
from scm.plams import Molecule
from scm.libbase import UnifiedChemicalSystem

mol = chemsys_to_plams_molecule(chemsys)
print(f"{type(chemsys)=}")
print(f"{type(mol)=}")
plot_molecule(mol)


# ### pymatgen Structure and Molecule Python classes

# Note that for this part of the example, the `pymatgen` package needs to be installed. This can be done via `amspackages`.


# #### Convert PLAMS Molecule to pymatgen Structure (periodic)
#
# There is no builtin converter between PLAMS Molecule and pymatgen Structure (periodic crystal). Instead, you need to go through the ASE interface to both packages:

from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
import scm.plams
from scm.plams import fromASE, toASE, Molecule
from ase.io import read


def convert_plams_molecule_to_pymatgen_structure(mol: Molecule) -> Structure:
    return AseAtomsAdaptor().get_structure(toASE(mol))


mol: scm.plams.Molecule = fromASE(read(cif_file))

pymatgen_structure: Structure = convert_plams_molecule_to_pymatgen_structure(mol)

print(f"{type(mol)=}")
print(f"{type(pymatgen_structure)=}")
print(pymatgen_structure)


# #### Convert pymatgen Structure (periodic) to PLAMS Molecule
#
# Go through the ASE interface:

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.structure import Structure
from scm.plams import fromASE
from scm.plams import Molecule


def pymatgen_structure_to_plams_molecule(pymatgen_structure: Structure) -> Molecule:
    return fromASE(AseAtomsAdaptor().get_atoms(pymatgen_structure))


print(f"{type(pymatgen_structure)=}")

mol = pymatgen_structure_to_plams_molecule(pymatgen_structure)
print(f"{type(mol)=}")


# #### Convert PLAMS Molecule to pymatgen Molecule (non-periodic)
#
# pymatgen has a special ``Molecule`` class for non-periodic systems. In PLAMS, the ``Molecule`` class is used for both periodic and non-periodic systems.

import pymatgen.core.structure
import scm.plams
from pymatgen.io.ase import AseAtomsAdaptor
from scm.plams import toASE


def convert_plams_molecule_to_pymatgen_molecule(
    mol: scm.plams.Molecule,
) -> pymatgen.core.structure.Molecule:
    return AseAtomsAdaptor().get_molecule(toASE(mol))


plams_molecule = scm.plams.Molecule(xyz_file)

pymatgen_molecule: pymatgen.core.structure.Molecule = convert_plams_molecule_to_pymatgen_molecule(plams_molecule)

print(f"{type(plams_molecule)=}")
print(f"{type(pymatgen_molecule)=}")
print(pymatgen_molecule)


# #### Convert pymatgen Molecule (non-periodic) to PLAMS Molecule

from pymatgen.io.ase import AseAtomsAdaptor
import pymatgen.core.structure
from scm.plams import fromASE
from scm.plams import Molecule


def pymatgen_molecule_to_plams_molecule(
    pymatgen_molecule: pymatgen.core.structure.Molecule,
) -> scm.plams.Molecule:
    return fromASE(AseAtomsAdaptor().get_atoms(pymatgen_molecule))


print(f"{type(pymatgen_molecule)=}")

mol = pymatgen_molecule_to_plams_molecule(pymatgen_molecule)
print(f"{type(mol)=}")
plot_molecule(mol)

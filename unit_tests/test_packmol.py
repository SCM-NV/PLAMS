import re

import numpy as np
import pytest
from dataclasses import dataclass
from typing import Optional, List, Union

from scm.plams.mol.molecule import Molecule, Atom
from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.interfaces.molecule.packmol import (
    PackMolStructure,
    packmol,
    guess_density,
    packmol_around,
)
from .test_helpers import skip_if_no_ams_installation


class TestPackmolStructure:

    @dataclass
    class PackMolTestCase:
        """
        Dataclass holding inputs for PackMolStructure and the expected values for attributes after initialization.
        """

        molecule: Molecule
        n_molecules: Optional[int] = None
        n_atoms: Optional[int] = None
        box_bounds: Optional[List[float]] = None
        density: Optional[float] = None
        fixed: bool = False
        sphere: bool = False
        expected_n_molecules: Optional[int] = None
        expected_box_bounds: Optional[List[float]] = None
        expected_fixed: Optional[bool] = False
        expected_sphere: Optional[bool] = False
        expected_error: Optional[str] = None
        expected_input: Optional[str] = None

    water = from_smiles("O")
    water_in_box = water.copy()
    water_in_box.lattice = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]

    test_cases = [
        # Happy
        PackMolTestCase(
            water,
            fixed=True,
            expected_n_molecules=1,
            expected_fixed=True,
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        PackMolTestCase(
            water,
            fixed=True,
            n_molecules=1,
            expected_n_molecules=1,
            expected_fixed=True,
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        PackMolTestCase(
            water,
            fixed=True,
            box_bounds=[0, 0, 0, 5, 5, 5],
            expected_n_molecules=1,
            expected_fixed=True,
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        PackMolTestCase(
            water_in_box,
            fixed=True,
            expected_n_molecules=1,
            expected_fixed=True,
            expected_box_bounds=[0, 0, 0, 10, 10, 10],
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        PackMolTestCase(
            water,
            fixed=True,
            sphere=True,
            expected_n_molecules=1,
            expected_fixed=True,
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        PackMolTestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            density=0.9,
            expected_n_molecules=30,
            expected_box_bounds=[0, 0, 0, 10, 10, 10],
            expected_input="""\
structure pm
  number 30
  inside box 0.05 0.05 0.05 9.95 9.95 9.95
end structure
""",
        ),
        PackMolTestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            n_molecules=20,
            expected_n_molecules=20,
            expected_box_bounds=[0, 0, 0, 10, 10, 10],
            expected_input="""\
structure pm
  number 20
  inside box 0.05 0.05 0.05 9.95 9.95 9.95
end structure
""",
        ),
        PackMolTestCase(
            water,
            n_molecules=25,
            density=1.0,
            expected_n_molecules=25,
            expected_box_bounds=[0.0, 0.0, 0.0, 9.077035146666434, 9.077035146666434, 9.077035146666434],
            expected_input="""\
structure pm
  number 25
  inside box 0.05 0.05 0.05 9.02703.* 9.02703.* 9.02703.*
end structure
""",
        ),
        PackMolTestCase(
            water,
            n_atoms=75,
            density=1.0,
            expected_n_molecules=25,
            expected_box_bounds=[0.0, 0.0, 0.0, 9.077035146666434, 9.077035146666434, 9.077035146666434],
            expected_input="""\
structure pm
  number 25
  inside box 0.05 0.05 0.05 9.02703.* 9.02703.* 9.02703.*
end structure
""",
        ),
        PackMolTestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            n_atoms=60,
            sphere=True,
            expected_n_molecules=20,
            expected_box_bounds=[0, 0, 0, 10, 10, 10],
            expected_sphere=True,
            expected_input="""\
structure pm
  number 20
  inside sphere 0. 0. 0. 6.20350.*
end structure
""",
        ),
        PackMolTestCase(
            water,
            n_molecules=0,
            box_bounds=[0, 0, 0, 5, 5, 5],
            expected_n_molecules=0,
            expected_box_bounds=[0, 0, 0, 5, 5, 5],
            expected_input="",
        ),
        # Unhappy
        PackMolTestCase(water, fixed=True, n_molecules=2, expected_error="n_molecules must be 1"),
        PackMolTestCase(water, fixed=True, density=42, expected_error="density cannot be set"),
        PackMolTestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            density=0.9,
            n_molecules=10,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        PackMolTestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            density=0.9,
            n_atoms=30,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        PackMolTestCase(
            water,
            n_molecules=30,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        PackMolTestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        PackMolTestCase(
            water,
            density=1,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        PackMolTestCase(
            water,
            n_atoms=10,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
    ]

    def get_structure_from_test_case(self, test_case):
        return PackMolStructure(
            molecule=test_case.molecule,
            n_molecules=test_case.n_molecules,
            n_atoms=test_case.n_atoms,
            box_bounds=test_case.box_bounds,
            density=test_case.density,
            fixed=test_case.fixed,
            sphere=test_case.sphere,
        )

    @pytest.mark.parametrize("test_case", test_cases)
    def test_structure_as_expected_given_inputs(self, test_case):

        if test_case.expected_error:
            with pytest.raises(ValueError, match=test_case.expected_error):
                self.get_structure_from_test_case(test_case)
        else:
            structure = self.get_structure_from_test_case(test_case)
            assert structure.n_molecules == test_case.expected_n_molecules
            assert structure.box_bounds == test_case.expected_box_bounds
            assert structure.fixed == test_case.expected_fixed
            assert structure.sphere == test_case.expected_sphere

    @pytest.mark.parametrize("test_case", test_cases)
    def test_get_input_block(self, test_case):
        if test_case.expected_error:
            with pytest.raises(ValueError, match=test_case.expected_error):
                self.get_structure_from_test_case(test_case)
        else:
            structure = self.get_structure_from_test_case(test_case)
            input_block = structure.get_input_block("pm", 0.1)

            assert re.fullmatch(test_case.expected_input, input_block, re.MULTILINE)


class TestPackMol:

    @dataclass
    class UnhappyTestCase:
        """
        Dataclass holding inputs for packmol and the expected error
        """

        molecules: Union[List[Molecule], Molecule]
        mole_fractions: Optional[List[float]] = None
        density: Optional[float] = None
        n_atoms: Optional[int] = None
        box_bounds: Optional[List[float]] = None
        n_molecules: Optional[Union[List[Optional[int]], int]] = None
        sphere: bool = False
        fix_first: bool = False
        region_names: Union[List[str], str, None] = None
        expected_error: str = ".*"

    water = from_smiles("O")
    acetonitrile = from_smiles("CC#N")
    ammonium = from_smiles("[NH4+]")
    chloride = from_smiles("[Cl-]")

    unhappy_test_cases = [
        UnhappyTestCase(molecules=water, expected_error="must specify either n_atoms, n_molecules or density"),
        UnhappyTestCase(
            molecules=water,
            n_atoms=300,
            n_molecules=100,
            expected_error="n_atoms and n_molecules are mutually exclusive",
        ),
        UnhappyTestCase(
            molecules=water,
            n_atoms=300,
            box_bounds=[0, 0, 0, 5, 5, 5],
            density=0.9,
            expected_error="n_atoms, box_bounds and density specified at the same time",
        ),
        UnhappyTestCase(
            molecules=water,
            n_molecules=100,
            box_bounds=[0, 0, 0, 5, 5, 5],
            density=0.9,
            expected_error="n_molecules, box_bounds and density specified at the same time",
        ),
        UnhappyTestCase(
            molecules=water,
            n_molecules=100,
            box_bounds=[0, 0, 0],
            expected_error=r"box_bounds must be a list of the form \[xmin, ymin, zmin, xmax, ymax, zmax\]",
        ),
        UnhappyTestCase(
            molecules=water,
            n_molecules=100,
            box_bounds=[0, 0, 0, 5, 5, None],
            expected_error=r"box_bounds must be a list of the form \[xmin, ymin, zmin, xmax, ymax, zmax\]",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=[100, 100],
            mole_fractions=[0.5, 0.5],
            expected_error="mole_fractions and n_molecules are mutually exclusive",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=1,
            fix_first=True,
            expected_error="fix_first requires that n_molecules is a list where the first element is 1",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=[42, 1],
            fix_first=True,
            expected_error="fix_first requires that n_molecules is a list where the first element is 1",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=42,
            expected_error="molecules is a list, but n_molecules is not",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile], n_molecules=[42], expected_error=r"len\(n_molecules\) != len\(molecules\)"
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            density=0.9,
            mole_fractions=42,
            expected_error="molecules is a list, but mole_fractions is not",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            density=0.9,
            mole_fractions=[42],
            expected_error=r"len\(mole_fractions\) != len\(molecules\)",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=[1, 2],
            region_names="foo",
            expected_error="molecules is a list, but region_names is not",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=[1, 2],
            region_names=["foo"],
            expected_error=r"len\(region_names\) != len\(molecules\)",
        ),
        UnhappyTestCase(
            molecules=water, n_molecules=[42], expected_error="n_molecules is a list, when molecules is not"
        ),
        UnhappyTestCase(
            molecules=water,
            density=0.9,
            mole_fractions=[0.5],
            expected_error="mole_fractions is a list, when molecules is not",
        ),
        UnhappyTestCase(
            molecules=water,
            n_molecules=42,
            region_names=["foo"],
            expected_error="region_names is a list, when molecules is not",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=[0, 0],
            expected_error="The sum of n_molecules is 0, which is very close to 0",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=[20, -10],
            expected_error="All n_molecules must be >= 0",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=[20, None],
            expected_error="the n_molecules list can contain a None value only if exactly two of box_bounds, density, and n_atoms are specified",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=[2000, None],
            density=0.9,
            box_bounds=[0, 0, 0, 10, 10, 10],
            expected_error="calculated value for missing n_molecules value is negative",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            n_molecules=[None, None],
            expected_error="multiple values of n_molecules are None",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            density=0.9,
            mole_fractions=[0, 0],
            expected_error="The sum of mole fractions is 0, which is very close to 0",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            density=0.9,
            mole_fractions=[1.0, -0.5],
            expected_error="All mole fractions must be >= 0",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            density=0.9,
            mole_fractions=[1.0, None],
            expected_error="one or more values of mole_fractions is None",
        ),
        UnhappyTestCase(
            molecules=[water, acetonitrile],
            density=0.9,
            mole_fractions=[0.5, 0.5],
            expected_error="n_atoms=None, n_molecules=None, box_bounds=None, density=0.9",
        ),
        UnhappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[None, None, 5],
            density=0.9,
            box_bounds=[10, 10, 10],
            expected_error="n_molecules, box_bounds and density specified at the same time",
        ),
        UnhappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[None, None, 5],
            density=0.9,
            expected_error="multiple values of n_molecules are None",
        ),
        UnhappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[None, 5, 5],
            density=0.9,
            expected_error="the n_molecules list can contain a None value only if exactly two of box_bounds, density, and n_atoms are specified",
        ),
        UnhappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[None, 5, 5],
            box_bounds=[0, 0, 0, 10, 10, 10],
            expected_error="the n_molecules list can contain a None value only if exactly two of box_bounds, density, and n_atoms are specified",
        ),
        UnhappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[None, 30, 5],
            density=0.9,
            box_bounds=[0, 0, 0, 10, 10, 10],
            expected_error="calculated value for missing n_molecules value is negative",
        ),
    ]

    @pytest.mark.parametrize("test_case", unhappy_test_cases)
    def test_pack_mol_unhappy(self, test_case):

        with pytest.raises(ValueError, match=test_case.expected_error):
            packmol(
                molecules=test_case.molecules,
                mole_fractions=test_case.mole_fractions,
                density=test_case.density,
                n_atoms=test_case.n_atoms,
                box_bounds=test_case.box_bounds,
                n_molecules=test_case.n_molecules,
                sphere=test_case.sphere,
                fix_first=test_case.fix_first,
                region_names=test_case.region_names,
                executable="no-op",
            )

    @dataclass
    class HappyTestCase:
        """
        Dataclass holding inputs for packmol and the expected return values
        """

        molecules: Union[List[Molecule], Molecule]
        mole_fractions: Optional[List[float]] = None
        density: Optional[float] = None
        n_atoms: Optional[int] = None
        box_bounds: Optional[List[float]] = None
        n_molecules: Optional[Union[List[Optional[int]], int]] = None
        sphere: bool = False
        fix_first: bool = False
        keep_bonds: bool = True
        keep_atom_properties: bool = True
        region_names: Union[List[str], str, None] = None
        expected_n_atoms: Optional[int] = None
        expected_n_molecules: Optional[List[int]] = None
        expected_mole_fractions: Optional[List[float]] = None
        expected_volume: Optional[float] = None
        expected_density: Optional[float] = None
        expected_radius: Optional[float] = None
        expected_region_names: Optional[List[str]] = None

    happy_test_cases = [
        HappyTestCase(
            water,
            n_atoms=194,
            density=1.0,
            expected_n_atoms=195,
            expected_n_molecules=[65],
            expected_mole_fractions=[1.0],
            expected_volume=1944.4885901260927,
            expected_density=1.0,
        ),
        HappyTestCase(
            water,
            density=1.0,
            box_bounds=[0.0, 0.0, 0.0, 8.0, 12.0, 14.0],
            expected_n_atoms=135,
            expected_n_molecules=[45],
            expected_mole_fractions=[1.0],
            expected_volume=1344,
            expected_density=1.0016253039797856,
        ),
        HappyTestCase(
            water,
            n_molecules=64,
            density=1.0,
            expected_n_atoms=192,
            expected_n_molecules=[64],
            expected_mole_fractions=[1.0],
            expected_volume=1914.5733810472298,
            expected_density=1.0,
        ),
        HappyTestCase(
            water,
            n_molecules=64,
            box_bounds=[0.0, 0.0, 0.0, 12.0, 13.0, 14.0],
            expected_n_atoms=192,
            expected_n_molecules=[64],
            expected_mole_fractions=[1.0],
            expected_volume=2184,
            expected_density=0.8766361634831642,
        ),
        HappyTestCase(
            [water],
            n_atoms=100,
            expected_n_atoms=99,
            expected_n_molecules=[33],
            expected_mole_fractions=[1.0],
            expected_volume=973.603502283497,
            expected_density=1.0139670792957152,
        ),
        HappyTestCase(
            [water, acetonitrile],
            n_atoms=200,
            mole_fractions=[0.666, 0.334],
            density=0.92,
            expected_n_atoms=201,
            expected_n_molecules=[33, 17],
            expected_mole_fractions=[0.66, 0.34],
            expected_volume=2332.679633984687,
            expected_density=0.92,
        ),
        HappyTestCase(
            [water, acetonitrile],
            n_atoms=200,
            mole_fractions=[0.666, 0.334],
            box_bounds=[0, 0, 0, 13.2, 13.2, 13.2],
            expected_n_atoms=201,
            expected_n_molecules=[33, 17],
            expected_mole_fractions=[0.66, 0.34],
            expected_volume=2299.968,
            expected_density=0.9330848356437622,
        ),
        HappyTestCase(
            [water, acetonitrile],
            n_molecules=[33, 17],
            density=0.92,
            expected_n_atoms=201,
            expected_n_molecules=[33, 17],
            expected_mole_fractions=[0.66, 0.34],
            expected_volume=2332.679633984687,
            expected_density=0.92,
        ),
        HappyTestCase(
            [water, acetonitrile],
            n_molecules=[33, 17],
            box_bounds=[0, 0, 0, 13.2, 13.2, 13.2],
            expected_n_atoms=201,
            expected_n_molecules=[33, 17],
            expected_mole_fractions=[0.66, 0.34],
            expected_volume=2299.968,
            expected_density=0.9330848356437622,
        ),
        HappyTestCase(
            molecules=[water],
            n_molecules=[100],
            density=1.0,
            sphere=True,
            expected_n_atoms=300,
            expected_n_molecules=[100],
            expected_mole_fractions=[1.0],
            expected_volume=2991.5209078862968,
            expected_density=1.0,
            expected_radius=8.93856517377096,
        ),
        HappyTestCase(
            [water, acetonitrile],
            n_molecules=[33, 17],
            sphere=True,
            expected_n_atoms=201,
            expected_n_molecules=[33, 17],
            expected_mole_fractions=[0.66, 0.34],
            expected_volume=2544.596756672291,
            expected_density=0.8433812774612032,
            expected_radius=8.469220758642498,
        ),
        HappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[6, 3, 1],
            density=0.4,
            sphere=True,
            expected_n_atoms=34,
            expected_n_molecules=[6, 3, 1],
            expected_mole_fractions=[0.6, 0.3, 0.1],
            expected_volume=820.55525521996,
            expected_density=0.4,
            expected_radius=5.807729858286062,
        ),
        HappyTestCase(
            [water, acetonitrile],
            n_atoms=200,
            mole_fractions=[0.666, 0.334],
            density=0.92,
            keep_bonds=False,
            expected_n_atoms=201,
            expected_n_molecules=[33, 17],
            expected_mole_fractions=[0.66, 0.34],
            expected_volume=2332.679633984687,
            expected_density=0.92,
            expected_region_names=["mol0", "mol0", "mol1", "mol1"],
        ),
        HappyTestCase(
            [water, acetonitrile],
            n_atoms=200,
            mole_fractions=[0.666, 0.334],
            density=0.92,
            region_names=["water", "acetonitrile"],
            expected_n_atoms=201,
            expected_n_molecules=[33, 17],
            expected_mole_fractions=[0.66, 0.34],
            expected_volume=2332.679633984687,
            expected_density=0.92,
            expected_region_names=["water", "acetonitrile"],
        ),
        HappyTestCase(
            [water],
            n_molecules=[None],
            density=1.0,
            box_bounds=[0.0, 0.0, 0.0, 8.0, 12.0, 14.0],
            expected_n_atoms=135,
            expected_n_molecules=[45],
            expected_mole_fractions=[1.0],
            expected_volume=1344,
            expected_density=1.0016253039797856,
        ),
        HappyTestCase(
            [water, acetonitrile],
            n_molecules=[33, None],
            density=0.92,
            box_bounds=[0, 0, 0, 10, 10, 23.32],
            expected_n_atoms=201,
            expected_n_molecules=[33, 17],
            expected_mole_fractions=[0.66, 0.34],
            expected_volume=2332,
            expected_density=0.9202681231843521,
        ),
        HappyTestCase(
            [water, acetonitrile],
            n_molecules=[None, 17],
            density=0.92,
            box_bounds=[0, 0, 0, 10, 10, 23.32],
            expected_n_atoms=201,
            expected_n_molecules=[33, 17],
            expected_mole_fractions=[0.66, 0.34],
            expected_volume=2332,
            expected_density=0.9202681231843521,
        ),
        HappyTestCase(
            molecules=[water, ammonium, chloride],
            density=0.4,
            box_bounds=[0, 0, 0, 5, 5, 5],
            expected_n_atoms=0,
            expected_n_molecules=[0, 0, 0],
            expected_mole_fractions=[0.0, 0.0, 0.0],
            expected_volume=125,
            expected_density=0.0,
        ),
        HappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[None, 2, 2],
            density=0.38,
            box_bounds=[0, 0, 0, 20, 20, 20],
            expected_n_atoms=300,
            expected_n_molecules=[96, 2, 2],
            expected_mole_fractions=[0.96, 0.02, 0.02],
            expected_volume=8000,
            expected_density=0.3811881797008519,
        ),
        HappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[96, None, 2],
            density=0.38,
            box_bounds=[0, 0, 0, 20, 20, 20],
            expected_n_atoms=300,
            expected_n_molecules=[96, 2, 2],
            expected_mole_fractions=[0.96, 0.02, 0.02],
            expected_volume=8000,
            expected_density=0.3811881797008519,
        ),
        HappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[96, 2, None],
            density=0.38,
            box_bounds=[0, 0, 0, 20, 20, 20],
            expected_n_atoms=300,
            expected_n_molecules=[96, 2, 2],
            expected_mole_fractions=[0.96, 0.02, 0.02],
            expected_volume=8000,
            expected_density=0.3811881797008519,
        ),
        HappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[96, 2, None],
            density=0.3811881797008519,
            n_atoms=300,
            expected_n_atoms=300,
            expected_n_molecules=[96, 2, 2],
            expected_mole_fractions=[0.96, 0.02, 0.02],
            expected_volume=8000,
            expected_density=0.3811881797008519,
        ),
        HappyTestCase(
            molecules=[water, ammonium, chloride],
            n_molecules=[96, 2, None],
            box_bounds=[0, 0, 0, 20, 20, 20],
            n_atoms=300,
            expected_n_atoms=300,
            expected_n_molecules=[96, 2, 2],
            expected_mole_fractions=[0.96, 0.02, 0.02],
            expected_volume=8000,
            expected_density=0.3811881797008519,
        ),
    ]

    @pytest.mark.parametrize("test_case", happy_test_cases)
    def test_pack_mol_happy_only_details(self, test_case):
        _, details = packmol(
            molecules=test_case.molecules,
            mole_fractions=test_case.mole_fractions,
            density=test_case.density,
            n_atoms=test_case.n_atoms,
            box_bounds=test_case.box_bounds,
            n_molecules=test_case.n_molecules,
            sphere=test_case.sphere,
            fix_first=test_case.fix_first,
            keep_bonds=test_case.keep_bonds,
            keep_atom_properties=test_case.keep_atom_properties,
            region_names=test_case.region_names,
            executable=".",  # hack, this just needs to be any path that exists, but the exe is not run
            _return_only_details=True,
        )

        assert details["n_atoms"] == test_case.expected_n_atoms
        assert details["n_molecules"] == test_case.expected_n_molecules
        assert details["mole_fractions"] == test_case.expected_mole_fractions

    @pytest.mark.parametrize("test_case", happy_test_cases)
    def test_pack_mol_happy(self, test_case):
        import numpy as np

        skip_if_no_ams_installation()

        mol, details = packmol(
            molecules=test_case.molecules,
            mole_fractions=test_case.mole_fractions,
            density=test_case.density,
            n_atoms=test_case.n_atoms,
            box_bounds=test_case.box_bounds,
            n_molecules=test_case.n_molecules,
            sphere=test_case.sphere,
            fix_first=test_case.fix_first,
            keep_bonds=test_case.keep_bonds,
            keep_atom_properties=test_case.keep_atom_properties,
            region_names=test_case.region_names,
            return_details=True,
        )

        assert details["n_atoms"] == test_case.expected_n_atoms
        assert details["n_molecules"] == test_case.expected_n_molecules
        assert details["mole_fractions"] == test_case.expected_mole_fractions
        assert np.isclose(details["volume"], test_case.expected_volume)
        assert np.isclose(details["density"], test_case.expected_density)
        if test_case.expected_radius:
            assert np.isclose(details["radius"], test_case.expected_radius)
        else:
            assert "radius" not in details

        assert len(mol.atoms) == test_case.expected_n_atoms
        mols = {}
        for m in mol.separate():
            f = m.get_formula()
            if f not in mols:
                mols[f] = [m, 1]
            else:
                mols[f][1] += 1

        if test_case.keep_bonds:
            assert [c for _, c in mols.values()] == [n for n in test_case.expected_n_molecules if n > 0]
            for m, _ in mols.values():
                if len(m) > 1:
                    assert m.bonds
        else:
            for m, _ in mols.values():
                assert not m.bonds

        if test_case.expected_region_names:
            for i, (m, _) in enumerate(mols.values()):
                assert all(a.properties.region == {test_case.expected_region_names[i]} for a in m.atoms)
        else:
            for i, (m, _) in enumerate(mols.values()):
                assert all(a.properties.region == {f"mol{i}"} for a in m.atoms)


class TestGuessDensity:

    @pytest.mark.parametrize(
        "smiles, expected, prev_implementation, reference",
        [
            ["O", 1.0140, 1.1392, 1],
            ["CO", 0.9035, 0.8465, 0.791],
            ["CCO", 0.8665, 0.7832, 0.789],
            ["CCCO", 0.8479, 0.7542, 0.803],
            ["CC(=O)C", 0.8284, 0.7439, 0.784],
            ["CC(=O)O", 0.9492, 0.8653, 1.049],
            ["CCOCC", 0.8368, 0.7375, 0.713],
            ["Cc1ccccc1", 0.7786, 0.6405, 0.867],
            ["c1ccccc1", 0.7762, 0.6329, 0.877],
            ["Oc1ccccc1", 0.8509, 0.6961, 1.07],
            ["c1ccccc1N", 0.8208, 0.6734, 1.021],
            ["c1ccc(cc1)[N+](=O)[O-]", 0.9360, 0.7680, 1.20],
            ["ClC(Cl)Cl", 0.8111, 0.7784, 1.49],
            ["CS(=O)C", 0.8151, 0.7319, 1.1004],
            ["COC(=O)C", 0.9216, 0.8233, 0.933],
            ["CCOC(=O)c1ccccc1", 0.8740, 0.7264, 1.045],
            ["CCCCCCCCO", 0.8170, 0.7091, 0.83],
            ["COC(=O)c1ccccc1C(=O)OC", 0.9354, 0.7799, 1.19],
            ["CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C", 0.9945, 0.8883, 0.96],
            ["C[Al](C)C", 0.8548, 0.8775, 0.752],
            ["C[Hg]C", 3.2238, 3.529, 2.961],
            ["He", 0.1189, 0.07978, 0.125],
            ["Li", 0.5296, 0.30693, 0.512],
            ["Be", 1.7124, 0.9923, 1.690],
            ["B", 0.3924, 0.3919, 2.08],
            ["Na", 0.9537, 0.5527, 0.927],
            ["Mg", 1.584, 0.9180, 1.584],
            ["Al", 2.4632, 1.4274, 2.375],
            ["Si", 2.3406, 2.7128, 2.57],
            ["Ti", 4.1059, 2.3794, 4.506],
            ["Fe", 7.4526, 4.3189, 7.86],
            ["Cu", 8.0890, 4.6876, 8.96],
            ["Ge", 5.3388, 6.1878, 5.323],
            ["As", 1.2193, 0.8178, 5.323],
            ["Zr", 6.2928, 3.6468, 6.52],
            ["Sn", 7.4538, 4.3195, 7.265],
            ["Sb", 5.7976, 6.7196, 6.697],
            ["Te", 1.4307, 0.9596, 6.240],
            ["W", 19.0859, 11.0605, 19.25],
            ["Pt", 19.8161, 11.4836, 21.450],
            ["Pb", 10.3207, 5.9810, 11.34],
            ["Po", 10.1663, 11.7831, 9.196],
            ["Es", 10.4543, 6.0584, 8.84],
        ],
        ids=[
            "water",
            "methanol",
            "ethanol",
            "propanol",
            "acetone",
            "acetic acid",
            "diethyl ether",
            "toluene",
            "benzene",
            "phenol",
            "aniline",
            "nitrobenzene",
            "chloroform",
            "dmso",
            "methyl acetate",
            "ethyl benzoate",
            "1-octanol",
            "dimethyl phthalate",
            "titanium isopropoxide",
            "trimethylaluminum",
            "dimethylmercury",
            "liquid_He",
            "liquid_Li",
            "liquid_Be",
            "liquid_B",
            "liquid_Na",
            "liquid_Mg",
            "liquid_Al",
            "liquid_Si",
            "liquid_Ti",
            "liquid_Fe",
            "liquid_Cu",
            "liquid_Ge",
            "liquid_As",
            "liquid_Zr",
            "liquid_Sn",
            "liquid_Sb",
            "liquid_Te",
            "liquid_W",
            "liquid_Pt",
            "liquid_Pb",
            "liquid_Po",
            "liquid_Es",
        ],
    )
    def test_guess_density(self, smiles, expected, prev_implementation, reference):
        try:
            mol = from_smiles(smiles)
        except Exception:
            mol = Molecule()
            atom = Atom(symbol=smiles)
            mol.add_atom(atom)

        mols = [mol]
        try:
            import scm.libbase  # noqa F401
            from scm.utils.conversions import plams_molecule_to_chemsys

            mols.append(plams_molecule_to_chemsys(mols[0]))
        except ImportError:
            pass

        for mol in mols:
            actual = guess_density([mol], [1])

            # keep track of previous implementation in case of making further improvements
            # uncomment this to see how they compare
            # abs_diff = np.abs(actual - reference)
            # abs_diff_prev = np.abs(prev_implementation - reference)
            # assert abs_diff < abs_diff_prev

            assert np.isclose(actual, expected, atol=0.0001), print(f"{actual=:}, {expected=}")
            if smiles not in ["B", "As", "Te"]:
                assert np.isclose(actual, reference, rtol=0.5), print(f"{actual=:}, {reference=}")


class TestPackMolAround:

    def test_pack_water_iteratively(self):
        skip_if_no_ams_installation()

        water = from_smiles("O")

        empty = Molecule()
        empty.lattice = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]

        num_atoms = []
        box = empty
        for _ in range(6):
            box = packmol_around(box, water, density=0.5)
            num_atoms.append(len(box))

        assert num_atoms == [51, 75, 87, 93, 96, 99]

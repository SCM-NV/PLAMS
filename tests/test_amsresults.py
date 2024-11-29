import numpy as np
import pytest
from unittest.mock import MagicMock
from ase import Atoms as AseAtoms

from scm.plams.interfaces.adfsuite.ams import AMSJob, AMSResults
from scm.plams.tools.kftools import KFFile
from scm.plams.core.errors import FileError, MissingOptionalPackageError
from scm.plams.mol.molecule import Molecule
from scm.plams.unit_tests.test_helpers import skip_if_no_ams_installation

# ToDo: Add tests for other job types e.g. MD, BAND etc. to test other result functions


class TestWaterOptimizationAMSResults:
    """
    Test suite for AMSResults, using rkf files from water optimization job.
    """

    @pytest.fixture
    def water_opt_results(self, rkf_folder):
        job = MagicMock(spec=AMSJob)
        job.status = "successful"
        job.path = str(rkf_folder / "water_optimization")
        results = AMSResults(job=job)
        return results

    def assert_water_molecule(self, actual, expected_coords, expected_type=Molecule):
        """
        Verify water molecule has the expected atom symbols and coordinates,
        and that the type is correct (Molecule/ASE Atoms/ChemicalSystem).
        """
        if isinstance(expected_coords, str):
            if expected_coords == "input":
                expected_coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
            elif expected_coords == "main":
                expected_coords = [
                    [0.0669210490964654, 0.06692104909646543, 0.0],
                    [1.0120415966947414, -0.07896264579120661, 0.0],
                    [-0.07896264579120661, 1.0120415966947414, 0.0],
                ]

        assert isinstance(actual, expected_type)
        if isinstance(actual, Molecule):
            assert [a.symbol for a in actual.atoms] == ["O", "H", "H"]
            assert np.allclose([a.coords for a in actual.atoms], expected_coords)
        elif isinstance(actual, AseAtoms):
            np.array_equal(actual.symbols.numbers, [8, 1, 1])
            np.array_equal(actual.positions, expected_coords)

    def test_collect_populates_files_and_rkfs(self, water_opt_results):
        # Given results of water optimization job
        # When call collect
        assert len(water_opt_results.files) == 0
        assert len(water_opt_results.rkfs) == 0
        water_opt_results.collect()

        # Then files populated with rkfs
        assert sorted(water_opt_results.files) == sorted(["ams.rkf", "dftb.rkf"])
        assert sorted([k for k in water_opt_results.rkfs.keys()]) == sorted(["dftb", "ams"])

    def test_collect_rkfs_populates_engine_rkfs_from_ams_rkf(self, water_opt_results):
        # Given results of water optimization job, with an ams.rkf in the present in files
        # When call collect rkfs
        assert len(water_opt_results.files) == 0
        assert len(water_opt_results.rkfs) == 0
        water_opt_results.files = ["ams.rkf"]
        water_opt_results.collect_rkfs()

        # Then rkfs populated
        assert sorted([k for k in water_opt_results.rkfs.keys()]) == sorted(["dftb", "ams"])

    def test_refresh_updates_moved_rkfs_and_deletes_removed_rkfs(self, water_opt_results, rkf_folder):
        # Given results of water optimization job, with two rkf files which have been moved, and one which does not exist
        water_opt_results.files = ["ams.rkf", "dftb.rkf", "does_not_exist.rkf"]
        water_opt_results.rkfs = {
            f.replace(".rkf", ""): KFFile(str(rkf_folder / "foo" / f)) for f in water_opt_results.files
        }

        # When refresh
        water_opt_results.refresh()

        # Then existing rkfs have paths updated to match the job path, and the non-existing file is removed
        assert sorted(water_opt_results.files) == sorted(["ams.rkf", "dftb.rkf"])
        assert sorted([k for k in water_opt_results.rkfs.keys()]) == sorted(["ams", "dftb"])

    def test_engine_names_returns_engine_specific_rkfs(self, water_opt_results):
        # Given results of water optimization job containing dftb engine rkf
        water_opt_results.collect()

        # When call get  engine names
        # Then dftb engine name returned
        assert water_opt_results.engine_names() == ["dftb"]

    def test_rkfpath_returns_absolute_path(self, water_opt_results, rkf_folder):
        # Given water optimization results
        # When get rkfpath for file
        # Then expected absolute path returned
        assert water_opt_results.rkfpath() == str(rkf_folder / "water_optimization" / "ams.rkf")
        assert water_opt_results.rkfpath("dftb") == str(rkf_folder / "water_optimization" / "dftb.rkf")
        with pytest.raises(FileError):
            assert water_opt_results.rkfpath("does_not_exist")

    @pytest.mark.parametrize(
        "section,variable,expected",
        [
            ["General", "task", "geometryoptimization"],
            ["InputMolecule", "nAtoms", 3],
            ["General", "Task", None],
            ["Foo", "Bar", None],
        ],
    )
    def test_readrkf_extracts_data_as_expected(self, water_opt_results, section, variable, expected):
        # Given water optimization results
        # When read rkf
        # Then expected values extracted
        water_opt_results.collect()
        if expected is not None:
            assert water_opt_results.readrkf(section, variable) == expected
        else:
            with pytest.raises(KeyError):
                water_opt_results.readrkf(section, variable)

    @pytest.mark.parametrize(
        "section,expected",
        [
            [
                "EngineResults",
                {
                    "Description(1)": "Results from the final step of the geometry optimization.",
                    "Files(1)": "dftb.rkf",
                    "Title(1)": "dftb",
                    "nEntries": 1,
                },
            ],
            ["general", {}],
            ["foo", {}],
        ],
    )
    def test_readrkf_section_extracts_data_as_expected(self, water_opt_results, section, expected):
        # Given water optimization results
        # When read rkf section
        # Then expected values extracted
        water_opt_results.collect()
        assert water_opt_results.read_rkf_section(section) == expected

    @pytest.mark.parametrize(
        "file,section,expected",
        [
            [
                "ams",
                "General",
                {
                    "CPUTime",
                    "ElapsedTime",
                    "SysTime",
                    "account",
                    "engine",
                    "file-ident",
                    "jobid",
                    "program",
                    "release",
                    "task",
                    "termination status",
                    "title",
                    "uid",
                    "user input",
                    "version",
                },
            ],
            [
                "dftb",
                "Molecule",
                {
                    "AtomMasses",
                    "AtomSymbols",
                    "AtomicNumbers",
                    "Charge",
                    "Coords",
                    "eeAttachTo",
                    "eeChargeWidth",
                    "eeEField",
                    "eeMulti",
                    "eeNLatticeVectors",
                    "eeNMulti",
                    "eeNZlm",
                    "eeUseChargeBroadening",
                    "eeXYZ",
                    "nAtoms",
                    "nAtomsTypes",
                },
            ],
            ["foo", None, None],
        ],
    )
    def test_get_rkf_skeleton_returns_expected_structure(self, water_opt_results, file, section, expected):
        # Given water optimization results
        # When get rkf skeleton
        # Then expected structure returned
        water_opt_results.collect()
        if expected is not None:
            skeleton = water_opt_results.get_rkf_skeleton(file)
            assert skeleton[section] == expected
        else:
            with pytest.raises(FileError):
                water_opt_results.get_rkf_skeleton(file)

    @pytest.mark.parametrize(
        "file,section,expected_coords",
        [
            ["ams", "InputMolecule", "input"],
            ["ams", "Molecule", "main"],
            ["dftb", "Molecule", "main"],
        ],
    )
    def test_get_molecule_returns_requested_section_molecule(self, water_opt_results, file, section, expected_coords):
        # Given water optimization results
        # When get molecule from given file section
        water_opt_results.collect()
        molecule = water_opt_results.get_molecule(section, file)

        # Then molecule as expected
        self.assert_water_molecule(molecule, expected_coords)

    @pytest.mark.parametrize(
        "file,section,expected_coords",
        [
            ["ams", "InputMolecule", "input"],
            ["ams", "Molecule", "main"],
            ["dftb", "Molecule", "main"],
        ],
    )
    def test_get_system_returns_requested_section_chemical_system_with_libbase_otherwise_errors(
        self, water_opt_results, file, section, expected_coords
    ):
        # Given water optimization results
        water_opt_results.collect()

        # When get molecule from given file section
        try:
            from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

            # Then molecule as expected when chemical system present
            molecule = water_opt_results.get_system(section, file)
            self.assert_water_molecule(molecule, expected_coords, ChemicalSystem)
        except ImportError:
            # Otherwise errors
            with pytest.raises(MissingOptionalPackageError):
                water_opt_results.get_system(section, file)

    @pytest.mark.parametrize(
        "section,expected_coords",
        [
            ["InputMolecule", "input"],
            ["Molecule", "main"],
        ],
    )
    def test_get_ase_atoms_returns_requested_section_molecule(self, water_opt_results, section, expected_coords):
        # Given water optimization results
        # When get ase atoms from given file section
        water_opt_results.collect()
        ase_atoms = water_opt_results.get_ase_atoms(section)

        # Then molecule as expected
        self.assert_water_molecule(ase_atoms, expected_coords, AseAtoms)

    def test_get_input_molecule_as_expected(self, water_opt_results):
        # Given water optimization results
        # When get input molecule
        water_opt_results.collect()
        molecule = water_opt_results.get_input_molecule()

        # Then molecule as expected (pre optimization)
        self.assert_water_molecule(molecule, "input")

    def test_get_input_system_as_expected_with_libbase_otherwise_errors(self, water_opt_results):
        # Given water optimization results
        water_opt_results.collect()

        # When get input molecule
        try:
            from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

            # Then molecule as expected when chemical system present
            molecule = water_opt_results.get_input_system()
            self.assert_water_molecule(molecule, "input", ChemicalSystem)
        except ImportError:
            # Otherwise errors
            with pytest.raises(MissingOptionalPackageError):
                water_opt_results.get_input_system()

    def test_get_input_molecules_has_initial_molecule_under_empty_key(self, water_opt_results):
        # Given water optimization results
        # When get input molecules with no named molecules
        water_opt_results.collect()
        molecule = water_opt_results.get_input_molecules()[""]

        # Then input molecule as expected (pre optimization) for empty key
        self.assert_water_molecule(molecule, "input")

    def test_get_main_molecule_as_expected(self, water_opt_results):
        # Given water optimization results
        # When get main molecule
        water_opt_results.collect()
        molecule = water_opt_results.get_main_molecule()

        # Then molecule as expected (post optimization)
        self.assert_water_molecule(molecule, "main")

    def test_get_main_system_as_expected(self, water_opt_results):
        # Given water optimization results
        # When get main molecule
        water_opt_results.collect()

        # Then molecule as expected (post optimization)
        try:
            from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

            # Then molecule as expected when chemical system present
            molecule = water_opt_results.get_main_system()
            self.assert_water_molecule(molecule, "main", ChemicalSystem)
        except ImportError:
            # Otherwise errors
            with pytest.raises(MissingOptionalPackageError):
                water_opt_results.get_main_system()

    @pytest.mark.parametrize("get_results", [False, True])
    def test_get_main_ase_atoms_as_expected(self, water_opt_results, get_results):
        # Given water optimization results
        # When get main ase atoms
        water_opt_results.collect()
        ase_atoms = water_opt_results.get_main_ase_atoms(get_results=get_results)

        # Then molecule as expected (post optimization)
        self.assert_water_molecule(ase_atoms, "main", AseAtoms)

        # And if results then energy/forces/charges decorated
        if get_results:
            assert ase_atoms.calc.results["energy"] == -156.90869381238653
            assert np.allclose(
                ase_atoms.calc.results["forces"],
                [
                    [-3.72826101e-04, -3.72826101e-04, -9.95035280e-16],
                    [1.11391933e-02, -1.07663672e-02, 1.54012759e-15],
                    [-1.07663672e-02, 1.11391933e-02, -5.45092314e-16],
                ],
            )
            assert np.allclose(ase_atoms.calc.results["charges"], [-0.67408472, 0.33704236, 0.33704236])

    @pytest.mark.parametrize(
        "step,expected_coords",
        [
            [1, "input"],
            [
                2,
                [
                    (0.04001658754120097, 0.04001658754120099, 0.0),
                    (1.0074239106376137, -0.04744049817881458, 0.0),
                    (-0.04744049817881458, 1.0074239106376137, 0.0),
                ],
            ],
            [
                3,
                [
                    (0.06476311169362219, 0.06476311169362221, 0.0),
                    (1.009209684748302, -0.07397279644192423, 0.0),
                    (-0.07397279644192423, 1.009209684748302, 0.0),
                ],
            ],
            [4, "main"],
        ],
    )
    def test_get_history_molecule_has_coordinates_for_expected_step(self, water_opt_results, step, expected_coords):
        # Given water optimization results
        # When get history molecule for given step
        water_opt_results.collect()
        molecule = water_opt_results.get_history_molecule(step)

        # Then molecule as expected (post optimization)
        self.assert_water_molecule(molecule, expected_coords)

    def test_is_valid_stepnumber_true_only_within_history_range(self, water_opt_results):
        # Given water optimization results
        # When check if valid stepnumber
        water_opt_results.collect()
        for i in range(-1, 6, 1):

            # Then valid if in history range
            if 0 < i < 5:
                assert water_opt_results.is_valid_stepnumber(water_opt_results.rkfs["ams"], i)
            else:
                with pytest.raises(KeyError):
                    water_opt_results.is_valid_stepnumber(water_opt_results.rkfs["ams"], i)

    def test_get_system_version_returns_none_when_not_present(self, water_opt_results):
        # Given water optimization results with no version info
        # When check if valid system version
        # Then get None
        water_opt_results.collect()
        assert water_opt_results.get_system_version(water_opt_results.rkfs["ams"], 1) is None

    def test_get_history_variables_as_expected(self, water_opt_results):
        # Given water optimization results with history
        # When check history variables
        # Then get all expected variables
        water_opt_results.collect()
        assert water_opt_results.get_history_variables() == {
            "Coords",
            "Energy",
            "Gradients",
            "ItemName",
            "maxGrad",
            "maxStep",
            "rmsGrad",
            "rmsStep",
        }

    def test_get_history_length_as_expected(self, water_opt_results):
        # Given water optimization results with history
        # When check history length
        # Then length matches number of steps
        water_opt_results.collect()
        assert water_opt_results.get_history_length() == 4

    def test_get_history_property_extracts_property_values_from_entries(self, water_opt_results):
        # Given water optimization results with history
        # When get property values
        # Then values returned for each step
        water_opt_results.collect()
        assert water_opt_results.get_history_property("Energy") == [
            -5.758015498298246,
            -5.765039556660849,
            -5.766261662993847,
            -5.7662881410725975,
        ]

    def test_get_property_at_step_extracts_property_values_from_entries_at_single_step(self, water_opt_results):
        # Given water optimization results with history
        # When get property values
        # Then values returned for each step
        water_opt_results.collect()
        assert water_opt_results.get_property_at_step(2, "Energy") == -5.765039556660849

    @pytest.mark.parametrize("engine,key,expected", [[None, "nOrbitals", 8], ["dftb", "nSpin", 1]])
    def test_get_engine_results_returns_all_results_as_dict(self, water_opt_results, engine, key, expected):
        # Given water optimization results with dftb engine
        # When get engine results
        water_opt_results.collect()
        results = water_opt_results.get_engine_results(engine)

        # Then dict of results returned
        assert results[key] == expected

    @pytest.mark.parametrize(
        "engine,key,expected",
        [[None, "DFTB Final Energy", -5.766288141072596], ["dftb", "Repulsion Energy", 0.038471496782890316]],
    )
    def test_get_engine_properties_results_returns_all_results_as_dict(self, water_opt_results, engine, key, expected):
        # Given water optimization results with dftb engine
        # When get engine properties
        water_opt_results.collect()
        results = water_opt_results.get_engine_properties()

        # Then dict of results returned
        assert results[key] == expected

    def test_get_energy_returns_in_given_unit(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get energy
        # Then given in correct unit
        water_opt_results.collect()
        assert water_opt_results.get_energy() == -5.766288141072596
        assert water_opt_results.get_energy("kJ/mol", "dftb") == -15139.387428634436

    def test_get_gradients_returns_in_given_unit(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get gradients
        # Then given in correct unit
        water_opt_results.collect()
        assert np.allclose(
            water_opt_results.get_gradients(),
            [
                [7.25031333e-06, 7.25031332e-06, 1.93503554e-17],
                [-2.16622821e-04, 2.09372508e-04, -2.99507132e-17],
                [2.09372508e-04, -2.16622821e-04, 1.06003578e-17],
            ],
        )
        assert np.allclose(
            water_opt_results.get_gradients("kJ/mol", "A", "dftb"),
            [
                [3.59722502e-02, 3.59722501e-02, 9.60063094e-14],
                [-1.07476877e00, 1.03879652e00, -1.48599722e-13],
                [1.03879652e00, -1.07476877e00, 5.25934130e-14],
            ],
        )

    def test_get_hessian_as_expected(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get hessian
        # Then get correct 9x9 matrix
        water_opt_results.collect()
        assert np.allclose(
            water_opt_results.get_hessian(),
            [
                [
                    5.25200795e-01,
                    -1.32396213e-01,
                    -1.90558991e-15,
                    -4.79749994e-01,
                    4.24745992e-02,
                    3.97943900e-15,
                    -4.54509777e-02,
                    8.99215247e-02,
                    -9.08916911e-14,
                ],
                [
                    -1.32396213e-01,
                    5.25200795e-01,
                    -2.70176255e-14,
                    8.99215244e-02,
                    -4.54509777e-02,
                    -2.34558963e-14,
                    4.24745996e-02,
                    -4.79749994e-01,
                    -6.88754534e-14,
                ],
                [
                    -1.90558991e-15,
                    -2.70176255e-14,
                    -5.29384822e-06,
                    8.90266923e-16,
                    1.52214697e-15,
                    2.77642325e-06,
                    -2.62293456e-14,
                    3.06345275e-14,
                    2.77642323e-06,
                ],
                [
                    -4.79749994e-01,
                    8.99215244e-02,
                    8.90266923e-16,
                    4.85406342e-01,
                    -8.19295674e-02,
                    1.64643292e-14,
                    -5.65626438e-03,
                    -7.99197226e-03,
                    8.81165912e-14,
                ],
                [
                    4.24745992e-02,
                    -4.54509777e-02,
                    1.52214697e-15,
                    -8.19295674e-02,
                    5.11073348e-02,
                    4.43914670e-16,
                    3.94550724e-02,
                    -5.65626438e-03,
                    -2.28694795e-14,
                ],
                [
                    3.97943900e-15,
                    -2.34558963e-14,
                    2.77642325e-06,
                    1.64643292e-14,
                    4.43914670e-16,
                    -1.03813440e-04,
                    2.56764504e-15,
                    1.24878717e-14,
                    1.00907518e-04,
                ],
                [
                    -4.54509777e-02,
                    4.24745996e-02,
                    -2.62293456e-14,
                    -5.65626438e-03,
                    3.94550724e-02,
                    2.56764504e-15,
                    5.11073348e-02,
                    -8.19295678e-02,
                    2.06439324e-15,
                ],
                [
                    8.99215247e-02,
                    -4.79749994e-01,
                    3.06345275e-14,
                    -7.99197226e-03,
                    -5.65626438e-03,
                    1.24878717e-14,
                    -8.19295678e-02,
                    4.85406342e-01,
                    9.84310364e-14,
                ],
                [
                    -9.08916911e-14,
                    -6.88754534e-14,
                    2.77642323e-06,
                    8.81165912e-14,
                    -2.28694795e-14,
                    1.00907518e-04,
                    2.06439324e-15,
                    9.84310364e-14,
                    -1.03813440e-04,
                ],
            ],
        )

    def test_get_frequencies_returns_in_given_unit(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get frequencies
        # Then returned in the given unit
        water_opt_results.collect()
        assert np.allclose(water_opt_results.get_frequencies(), [1427.92373935, 3674.50690527, 3785.96039844])
        assert np.allclose(water_opt_results.get_frequencies("eV", "dftb"), [0.17703998, 0.45558079, 0.46939927])

    def test_get_force_constants_as_expected(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get force contants
        # Then returned as expected
        water_opt_results.collect()
        assert np.allclose(water_opt_results.get_force_constants(), [0.08366119, 0.53330315, 0.58846565])

    def test_get_normal_modes_as_expected(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get normal modes
        # Then returned as expected
        water_opt_results.collect()
        assert np.allclose(
            water_opt_results.get_normal_modes(),
            [
                [
                    [5.04886076e-02, 5.04886076e-02, -1.41263536e-22],
                    [-1.03645924e-01, -6.97644912e-01, -2.31906718e-17],
                    [-6.97644912e-01, -1.03645924e-01, -2.54429472e-29],
                ],
                [
                    [3.46056664e-02, 3.46056660e-02, 3.99118436e-24],
                    [-6.91731267e-01, 1.42514233e-01, 3.25036299e-18],
                    [1.42514232e-01, -6.91731261e-01, -6.01340612e-29],
                ],
                [
                    [-5.06985283e-02, 5.06985285e-02, 1.76560066e-22],
                    [6.97032293e-01, -1.07590134e-01, 1.13086149e-16],
                    [1.07590135e-01, -6.97032299e-01, 1.48847456e-30],
                ],
            ],
        )

    def test_get_charges_as_expected(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get charges
        # Then returned as expected
        water_opt_results.collect()
        assert np.allclose(water_opt_results.get_charges(), [-0.67408472, 0.33704236, 0.33704236])

    def test_get_dipolemoment_as_expected(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get dipolemoment
        # Then returned as expected
        water_opt_results.collect()
        assert np.allclose(water_opt_results.get_dipolemoment(), [0.50904814, 0.50904814, 0.0])

    def test_get_zero_point_energy_returns_in_given_unit(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get polarizability
        # Then returns in given unit
        water_opt_results.collect()
        assert water_opt_results.get_zero_point_energy() == 0.020249244724897086
        assert water_opt_results.get_zero_point_energy("kJ/mol", "dftb") == 53.164384700766135

    def test_get_ir_intensities_as_expected(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get ir intensitites
        # Then returned as expected
        water_opt_results.collect()
        assert np.allclose(water_opt_results.get_ir_intensities(), [126.33789359, 31.24055122, 78.6282858])

    def test_get_orbital_energies_and_occupations_and_homo_lumo_consistent(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get number of orbital energies, the homo/lumo and the gap
        water_opt_results.collect()
        energies = water_opt_results.get_orbital_energies()[0]
        occupations = water_opt_results.get_orbital_occupations()[0]
        homo = water_opt_results.get_homo_energies()[0]
        lumo = water_opt_results.get_lumo_energies()[0]
        gap = water_opt_results.get_smallest_homo_lumo_gap()

        # Then homo/lumo consistent with the highest and lowest (un)occupied energies
        assert water_opt_results.get_n_spin() == 1
        assert not water_opt_results.are_orbitals_fractionally_occupied()
        assert np.allclose(
            energies,
            [-0.75775416, -0.61195531, -0.54473148, -0.49953169, -0.15455526, -0.04850734, 0.34764644, 0.45289103],
        )
        assert np.allclose(occupations, [2, 2, 2, 2, 0, 0, 0, 0])
        assert homo == energies[3]
        assert lumo == energies[4]
        assert gap == (lumo - homo)

    def test_get_timings_as_expected(self, water_opt_results):
        # Given water optimization results with dftb engine
        # When get timings
        # Then as expected
        water_opt_results.collect()
        assert water_opt_results.get_timings() == {"cpu": 0.732007, "elapsed": 0.8396751880645752, "system": 0.042673}

    def test_recreate_molecule_returns_input_molecule(self, water_opt_results):
        # Given water optimization results
        # When get recreate molecule
        # Then returns a molecule the same as the input molecule
        water_opt_results.collect()
        self.assert_water_molecule(water_opt_results.recreate_molecule(), "input")

    def test_recreate_settings_returns_input_settings(self, water_opt_results):
        skip_if_no_ams_installation()
        # Given water optimization results
        # When get recreate job settings
        # Then returns the same settings
        water_opt_results.collect()
        assert water_opt_results.recreate_settings()["input"] == {
            "DFTB": {"Model": "GFN1-xTB"},
            "ams": {"Properties": {"NormalModes": "yes"}, "Task": "GeometryOptimization"},
        }

    def test_ok_return_value_reflects_job_ok(self, water_opt_results):
        # Given water optimization results
        # When get toggle job settings ok return value
        call_count = 0
        water_opt_results.job.ok.side_effect = lambda: call_count % 2 == 0

        # Then mirrored in the results return value
        assert water_opt_results.ok()
        call_count += 1
        assert not water_opt_results.ok()

    def test_get_errormsg_return_value_reflects_job_errormsg(self, water_opt_results):
        # Given water optimization results
        # When get toggle job settings ok return value
        call_count = 0
        water_opt_results.job.get_errormsg.side_effect = lambda: None if call_count % 2 == 0 else "error"

        # Then mirrored in the results return value
        assert water_opt_results.get_errormsg() is None
        call_count += 1
        assert water_opt_results.get_errormsg() == "error"

    def test_name_reflects_job_name(self, water_opt_results):
        # Given water optimization results
        # When set job name
        water_opt_results.job.name = "foo"

        # Then mirrored in the results return value
        assert water_opt_results.name == "foo"


class TestPropaneNitrileOptimizationAMSResults:
    """
    Test suite for AMSResults, using rkf files from propane nitrile optimization job.
    """

    @pytest.fixture
    def propane_nitrile_opt_results(self, rkf_folder):
        job = MagicMock(spec=AMSJob)
        job.status = "successful"
        job.path = str(rkf_folder / "propanenitrile")
        results = AMSResults(job=job)
        return results

    def test_get_atom_types_as_expected(self, propane_nitrile_opt_results):
        # Given propane nitrile optimization results with ff engine
        # When get atom types
        # Then returned as expected
        propane_nitrile_opt_results.collect()
        assert (
            propane_nitrile_opt_results.get_atom_types()
            == propane_nitrile_opt_results.get_atom_types("forcefield")
            == ["C_1", "C_3", "C_3", "N_1", "H_", "H_", "H_", "H_", "H_"]
        )

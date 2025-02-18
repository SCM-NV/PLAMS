import pytest

from scm.plams.mol.molecule import Molecule


class TestIdentify:
    """
    Test suite for identification methods for plams Molecule
    """

    @pytest.mark.parametrize(
        "mol1_name, mol2_name, max_expected",
        [
            ("benzene", "chlorophyl1", None),
            ("CO_6_1", "CO_6_3", None),
            ("CO_flat4_1", "CO_flat4_3", None),
            ("chlorophyl1", "chlorophyl2", 1),
            ("EZ1", "EZ2", 2),
            ("RS1", "RS2", 2),
            ("CO_6_1", "CO_6_2", 3),
            ("CO_6_3", "CO_6_4", 3),
            ("CO_6_3", "CO_6_5", 3),
            ("CO_6_4", "CO_6_5", 3),
            ("CO_flat4_1", "CO_flat4_2", 3),
            ("CO_flat4_3", "CO_flat4_4", 3),
            ("benzene", "benzene", 4),
        ],
    )
    def test_label_matches_as_expected_for_molecule_pair(self, mol1_name, mol2_name, max_expected, xyz_folder):
        mol1 = Molecule(xyz_folder / f"{mol1_name}.xyz")
        mol2 = Molecule(xyz_folder / f"{mol2_name}.xyz")
        mol1.guess_bonds()
        mol2.guess_bonds()
        for i in range(5):
            assert (mol1.label(i) == mol2.label(i)) == (False if max_expected is None else i <= max_expected), print(
                f"Failed for level {i}"
            )

    def test_reorder_as_expected(self, xyz_folder):
        chl1 = Molecule(xyz_folder / "chlorophyl1.xyz")
        chl2 = Molecule(xyz_folder / "chlorophyl2.xyz")
        benzene = Molecule(xyz_folder / "benzene.xyz")

        # Reorder with different mismatching molecules gives None
        assert chl1.reorder(benzene) is None

        # Reorder maintains the geometry of the first molecule, but with the atom ordering of the second
        chl3 = chl1.reorder(chl2)
        assert chl3.label(4) == chl1.label(4) != chl2.label(4)
        assert [at.symbol for at in chl3] == [at.symbol for at in chl2] != [at.symbol for at in chl1]

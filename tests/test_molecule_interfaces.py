import pytest
from unittest.mock import patch
import re

from scm.plams.mol.molecule import Molecule
from scm.plams.interfaces.molecule.ase import toASE, fromASE
from scm.plams.unit_tests.test_helpers import get_mock_find_spec, get_mock_open_function
from scm.plams.core.errors import MissingOptionalPackageError
from scm.plams.interfaces.molecule.rdkit import from_rdmol, to_rdmol, from_smiles, to_smiles, from_smarts
from scm.plams.interfaces.molecule.packmol import packmol


@pytest.fixture
def plams_mols(xyz_folder, pdb_folder, rkf_folder):
    """
    Selection of plams molecules loaded from various sources.
    """
    water_molecule = Molecule(rkf_folder / "water_optimization" / "ams.rkf")
    water_molecule_with_bonds = water_molecule.copy()
    water_molecule_with_bonds.guess_bonds()
    water_molecule_in_box = water_molecule.copy()
    water_molecule_in_box.lattice = [[100, 0, 0], [0, 100, 0], [0, 0, 100]]
    benzene = Molecule(xyz_folder / "benzene.xyz")
    chlorophyl = Molecule(xyz_folder / "chlorophyl1.xyz")
    chymotrypsin = Molecule(pdb_folder / "chymotrypsin.pdb")

    return {
        "water": water_molecule,
        "water_bonds": water_molecule_with_bonds,
        "water_box": water_molecule_in_box,
        "benzene": benzene,
        "chlorophyl": chlorophyl,
        "chymotrypsin": chymotrypsin,
    }


class TestASE:
    """
    Test suite for conversions to/from PLAMS molecule and ASE.
    """

    def test_to_ase_from_ase_roundtrip(self, plams_mols):
        for name, orig_mol in plams_mols.items():
            print(f"Testing roundtrip for molecule '{name}'")
            converted_mol = toASE(orig_mol)
            final_mol = fromASE(converted_mol)
            assert final_mol.label(4) == orig_mol.label(4)

    def test_to_ase_requires_ase_package(self, plams_mols):
        with get_mock_find_spec("scm.plams.core.functions", "ase"):
            with pytest.raises(MissingOptionalPackageError):
                toASE(plams_mols["water"])


class TestRDKit:
    """
    Test suite for conversions to/from PLAMS molecule and RDKit.
    """

    @pytest.fixture
    def smiles(self):
        return [
            "C(=O)(O)C(C)(C)Oc1ccc(cc1)CC(=O)Nc1cc(cc(c1)C)C",
            "CC(C)([C@@H](C(=O)NC)NC(=O)[C@@H](CN(C=O)O)CC(C)C)C",
            "O=C1N2[C@@H](C(=O)[O-])C(C)(C)[S@](=O)[C@@H]2[C@@H]1NC(=O)Cc1ccccc1",
            "c1cc(cc2c1CC[NH2+]C2)S(=O)(=O)N",
            "OC[C@H]1O[C@H]([C@@H]([C@@H]1O)O)n1cnc2c(N)nccc12",
            "C1Nc2nc(N)[nH]c(=O)c2N=C1CO",
            "n1c(C)[nH]c(=O)c2cc(ccc12)CN(c1sc(cc1)C(=O)N[C@@H](CCC(=O)O)C(=O)O)C",
            "c12c(cccc1)n(c(c2c1ccc(cc1)F)/C=C/[C@H](C[C@@H](O)CC(=O)O)O)C(C)C",
            "c12c(c(nc([nH+]1)N)N)c(ccc2)Sc1ccccc1",
            "Cc1nc(N)c(C[n+]2csc(CCO)c2C)cn1",
            "n1c([nH+]c(c(c1N)c1ccc(cc1)Cl)CC)N",
            "c1(nnc(NC(=O)C)s1)S(=O)(=O)[NH-]",
            "[C@H](C(=O)O)(Cc1ccccc1)[C@@H](C(=O)O)Cc1cc2OCOc2cc1",
            "c1(ccccc1)Cc1n(c(=O)[nH]c(=O)c1C(C)C)COCc1ccccc1",
            "O=C(O)C[C@@H](C(=O)O)NC(=O)Cc1c[nH]c2c1cccc2",
            "N1C(=O)/C(=C\\Nc2ccc(cc2)S(=O)(=O)NC)/c2ccccc12",
            "c1(O)cccc(c1C)C(=O)N[C@@H](Cc1ccccc1)[C@H](O)C(=O)N1CSC(C)(C)[C@H]1C(=O)NCc1c(cccc1)C",
            "[C@H]1(NC(=[NH2+])N)[C@H]([C@H](C(CC)CC)NC(=O)C)[C@H](O)[C@@H](C(=O)O)C1",
            "c1cc2c(c(c1)C)cc(n2Cc1cccc(c1)C(=[NH2+])N)C(=O)NCc1cc(cc(c1)Cl)Cl",
            "c12c(cccc1cccc2)CC(=O)O",
            "C1=CC(=O)C=C2CC[C@@H]3[C@@]([C@@]12C)([C@H](C[C@]1([C@H]3C[C@H]([C@@]1(C(=O)CO)O)C)C)O)F",
            "C1(=O)c2c(CO1)c(C)c(c(c2O)C/C=C(/CCC(=O)O)\\C)OC",
            "[NH3+][C@H](C(=O)O)CCCNC(=[NH2+])NCCC",
            "[C@@]1(c2cc(Oc3c(ccc(c3)[C@](c3n(cnc3)C)(C)[NH3+])C#N)ccc2)(C(=O)N(CCCC1)C)CC",
            "C1CN(CC1)C(=O)[C@H](C(C)C)[NH3+]",
            "O=C(O)[C@H](O)C(C)(C)CO",
            "[nH]1c(CCCC)nc2c(=O)[nH][nH]c(=O)c12",
            "O=c1[nH]c(=O)cnn1c1cc(C)c(c(c1)C)Oc1ccc(c(c1)C(C)C)O",
            "O=C(O)Cc1cc(c(Oc2ccc(O)c(c2)C(C)C)c(c1)Cl)Cl",
            "OC[C@@H]1[C@@H](O)C[C@]2(n3c(=O)[nH]c(=O)c(c3)C)[C@H]1C2",
            "[NH3+][C@H](C(=O)O)Cc1ccc(cc1)O",
            "OCc1cccc(Nc2ncc3cc(c(=O)n(c3n2)C)c2c(Cl)cccc2Cl)c1",
            "S(=O)(=O)(c1ccc(cc1)n1c(c2ccc(C)cc2)cc(C(F)(F)F)n1)[NH-]",
            "[NH2+]=C(N)c1ccc2cc(ccc2c1)C(=O)Nc1ccccc1",
            "O=C1N(Cc2ccc(F)cc2)C(=O)[C@H]2[C@H]3N(CCC3)[C@@H](c3ccc(cc3)C(=[NH2+])N)[C@@H]12",
            "n1cc(ccc1)[C@H]1[N@@H+](CCC1)C",
            "OC[C@H]1O[C@@H](n2ccc(nc2=O)N)C(F)(F)[C@@H]1O",
            "n1(c(nc(c1c1ccnc(n1)NC1CC1)c1cc(c(cc1)Cl)Cl)[C@H]1CC[N@@H+](CC1)C)CCC",
            "CSC[C@H]1[NH2+][C@H]([C@H](O)[C@@H]1O)c1c[nH]c2c(=O)[nH]cnc12",
            "c12/C(=N\\O)/C(=C\\3/C(=O)Nc4c3cccc4)/Nc1cccc2",
            "c1cc([C@@H](C(=O)O)C)ccc1c1ccccc1",
            "C[C@H]([NH3+])[P@](=O)(O)C[C@H](C(=O)N[C@@H](C)C(=O)O)Cc1ccc(cc1)c1ccccc1",
            "N(C(=O)[C@@H]([C@@H](C(=O)NO)O)CC(C)C)[C@@H](C(C)(C)C)C(=O)NC",
            "C(C)(C)SCC[C@@H](N)[C@H](O)C(=O)NNC(=O)c1cccc(c1)Cl",
            "c1cc(ccc1)c1ccc(cc1F)[C@H](C)C(=O)O",
            "c1([nH+]c2c(c(n1)N)C[C@@H](CC2)CN(c1cc(c(c(c1)OC)OC)OC)C)N",
            "c1(cc(cc(c1)/C=C/c1ccc(cc1)O)O)O",
            "[C@H]1([C@@H](Oc2ccc(O)cc2S1)c1ccc(cc1)OCC[NH+]1CCCCC1)c1ccc(O)cc1",
            "OCC(C)(C)[C@@H](O)C(=O)NCCC(=O)O",
            "c1c(c(ccc1F)C(=O)NCc1nc2c(s1)c(F)cc(F)c2F)OCC(=O)O",
            "c1cc(cnc1)c1ccnc(n1)Nc1c(ccc(c1)NC(=O)c1ccc(cc1)CN1CC[N@H+](C)CC1)C",
            "COc1nc(C)nc(N/C(=N/S(=O)(=O)c2ccccc2Cl)/O)n1",
            "O=C(O)CCCn1c2ccccc2c2c1cccc2",
            "[NH2+]1C[C@@H]([C@@H]([C@H]1C(=O)O)CC(=O)O)C(=C)C",
            "CC/C(=C(\\c1ccc(O)cc1)/CC)/c1ccc(O)cc1",
            "OCCOCn1c(=O)[nH]c(=O)c(c1)Cc1ccccc1",
            "O=C1NCC/C(=C\\2/NC(=NC2=O)N)/c2c1[nH]cc2",
            "c1n(cnc1C(=O)N)[C@@H](CO)CCn1ccc2c1cc(cc2)NC(=O)CCc1ccccc1",
            "OC[C@@H](CC)Nc1nc2n(C(C)C)cnc2c(n1)NCc1ccccc1",
            "Clc1c(CN2CCCC2=[NH2+])[nH]c(=O)[nH]c1=O",
            "n1c(nc2n(C(C)C)cnc2c1Nc1ccc(c(Cl)c1)C(=O)O)N[C@H](C(C)C)CO",
            "n1(cnc2c(=O)[nH]c(N)nc12)CCCCC(F)(F)P(=O)(O)O",
            "s1ccnc1NC(=O)c1c(cc(c(c1)Sc1n(ccn1)C)F)N",
            "c1c(c(cc(c1)C(=O)O)NC(CC)CC)N1[C@@](CCC1=O)(CO)C[NH3+]",
            "O=C1[C@@H]2CCCN2C(=O)CN1",
            "OC[C@H]1O[C@H](C[C@@H]1O)n1c(=O)[nH]c(=O)c(C)c1",
            "O=C1O[C@](CN1)(C)c1ccc(OC)c(c1)OCCC",
            "O=C(c1ccc(OC(F)F)c(OCC2CC2)c1)Nc1c(Cl)cncc1Cl",
            "c1ccc2c3c([nH]c2c1)[C@H](N1C(=O)CN(C(=O)[C@H]1C3)C)c1ccc2OCOc2c1",
            "O=S(=O)(NCC1CC1)c1ccc(c(c1)Nc1oc(cn1)c1cc(ccc1)c1c[nH+]ccc1)OC",
            "Oc1cc(N[C@@H](C(=O)NS(=O)(=O)c2cccc(N)c2)c2c(F)c(OCC)cc(OCC)c2)ccc1C(=[NH2+])N",
            "c1(cc(c(cc1)F)C)S(=O)(=O)N[C@@H](C(=O)NO)C1CCOCC1",
            "C(=C(\\NC(=O)c1ccccc1)/C(=O)O)\\c1ccc(cc1)Oc1c(cccc1)Br",
            "c1c(ccc(c1)F)c1c(c2nc(N[C@H](c3ccccc3)C)ncc2)n(n(C2CC[NH2+]CC2)c1=O)C",
            "c1cc(F)ccc1S(=O)(=O)C[C@@](O)(C)C(=O)Nc1cc(c(cc1)C#N)C(F)(F)F",
            "c1(cc(ccc1)C[NH3+])[C@H]1CCN(C(=O)c2cc(CCc3ccccc3)cnc2)CC1",
            "COc1ccc(cc1)c1c2c(NCCO)ncnc2oc1c1ccc(OC)cc1",
            "COc1ccc(cc1)c1c(C(=O)NCC)n[nH]c1c1cc(Cl)c(O)cc1O",
            "C(=O)(c1ccccc1)NNc1n[nH]c(=S)n1N",
            "Brc1n(c(C(C)C)c(n1)NC(=O)C)C",
            "c1(c(n(cc1NS(=O)(=O)c1ccccc1)C(C)(C)C)N(C)C)c1ccccc1",
            "Brc1ccc(cc1)/C(=C/c1ccccc1)/C1=CC(=O)N(C1)C(=O)C",
            "c12c(ccc(c1c(=O)cc(C(=O)[O-])o2)OCCC(C)C)CC=C",
            "S(=O)(=O)(c1ccc(N(=O)=O)cc1)Cc1ccc(N(C)C)cc1",
            "s1cc(n(/c/1=C/N(=O)=O)CCNC(=O)C)c1ccccc1",
            "c1cc(c2c(c3ccccc3)c(N(C(=O)C)C(=O)C)on2)ccc1",
            "c1(cc(c(cc1)NC(=O)C)SC)c1ccc(cc1)NC(=O)C",
            "n1cnc2c(c1N)ncn2CCC(=O)NCCc1c[nH]c2c1cccc2",
            "ClCS(=O)(=O)/C(=C(\\c1ccccc1)/N1CCOCC1)/C",
            "n1c(Cl)c(ccc1)NCc1nc(Cl)c(N)cc1",
            "c12c(=O)n(c3ccccc3)[nH]c1nc(Nc1ccccc1)s2",
            "s1c(nc(c1N(c1ccccc1)C)C)Nc1ccccc1",
            "s1c(nc(c1N(c1ccccc1)C)c1ccccc1)Nc1ccccc1",
            "c1cc(ccc1C(=C1CCCCC1)c1ccc(cc1)OC(=O)C)OC(=O)C",
            "n1cnc2c(c1N)ncn2CCC(=O)NCCc1ccc(cc1)O",
            "N(C)(C)c1ccc(cc1)C(=[NH2+])c1ccc(N(C)C)cc1",
            "c1(CCC(=O)O)cc(c2cccccc12)CCC(=O)O",
            "O1CC[NH+](CC1)CCNc1ccc(n[nH+]1)c1ccccc1",
            "n1c(nc(cc1)N)SSc1nccc(n1)N",
            "Fc1c(=O)[nH]c(=O)n(c1)c1nc(nc(c1)C)Cc1ccc(OC)cc1",
        ]

    def roundtrip_and_assert(self, molecules, from_mol, to_mol, level=4):
        for name, orig_mol in molecules.items():
            print(f"Testing roundtrip for molecule '{name}'")
            if level > 0:
                converted_mol = from_mol(orig_mol)
                final_mol = to_mol(converted_mol)
                assert final_mol.label(level) == orig_mol.label(level)
            else:
                with pytest.raises(Exception):
                    converted_mol = from_mol(orig_mol)
                    to_mol(converted_mol)

    def test_to_rdmol_from_rdmol_roundtrip(self, plams_mols):
        from_mol = lambda mol: to_rdmol(mol)
        to_mol = lambda mol: from_rdmol(mol)
        self.roundtrip_and_assert(plams_mols, from_mol, to_mol)

    @pytest.mark.parametrize("short_smiles,ff", [(True, None), (False, None), (False, "uff"), (False, "mmff")])
    def test_to_smiles_from_smiles_roundtrip(self, plams_mols, short_smiles, ff):
        from_mol = lambda mol: to_smiles(mol, short_smiles=short_smiles)
        to_mol = lambda mol: from_smiles(mol, forcefield=ff)

        # For benzene, we get Kekule bonds instead of original aromatic...
        lvl1_mols = {k: v for k, v in plams_mols.items() if k == "benzene"}

        # For chlorophyl we get a kekulize error...
        # For chymotrypsin we get a aromatic error...
        err_mols = {k: v for k, v in plams_mols.items() if k in ["chlorophyl", "chymotrypsin"]}

        # Everything else should match up to bond order level
        lvl2_mols = {k: v for k, v in plams_mols.items() if k not in lvl1_mols.keys() and k not in err_mols.keys()}

        self.roundtrip_and_assert(lvl2_mols, from_mol, to_mol, level=2)
        self.roundtrip_and_assert(lvl1_mols, from_mol, to_mol, level=1)
        self.roundtrip_and_assert(err_mols, from_mol, to_mol, level=-1)

    def test_from_smiles_to_smiles_from_smiles_roundtrip(self, smiles):
        from_mol = lambda mol: to_smiles(mol)
        to_mol = lambda mol: from_smiles(mol)
        smile_mols = {s: from_smiles(s) for s in smiles}
        self.roundtrip_and_assert(smile_mols, from_mol, to_mol, level=2)

    def test_to_smiles_as_expected(self, plams_mols):
        water = plams_mols["water"]
        benzene = plams_mols["benzene"]

        assert to_smiles(water) == "O"
        assert to_smiles(water, short_smiles=False) == "[H]O[H]"
        assert to_smiles(benzene) == "c1ccccc1"
        assert to_smiles(benzene, short_smiles=False) == "[H]c1c([H])c([H])c([H])c([H])c1[H]"

    def test_from_smiles_as_expected(self):
        hexane = "CCCCCC"

        assert from_smiles(hexane).get_formula() == "C6H14"
        assert len(from_smiles(hexane, nconfs=3, rms=100)) == 1
        assert len(from_smiles(hexane, nconfs=3, rms=0.1)) == 3

    def test_from_smarts_as_expected(self, smiles):
        assert from_smarts("[#6]:1:[#6]:[#6]:[#6]:[#6]:[#6]:1").get_formula() == "C6H6"

    def test_to_smiles_and_from_smiles_requires_rdkit_package(self, plams_mols, smiles):
        with get_mock_find_spec("scm.plams.core.functions", "rdkit"):
            with pytest.raises(MissingOptionalPackageError):
                to_smiles(plams_mols["water"])
            with pytest.raises(MissingOptionalPackageError):
                from_smiles(smiles[0])


class TestPackmol:
    """
    Test suite for Packmol interface
    """

    @pytest.fixture
    def five_water(self):
        return """15
Built with Packmol
     O       3.93094000       1.81434900       4.06953500
     H       4.37496800       0.94573600       3.97345900
     H       2.99438400       1.71747900       4.34826900
     O       0.94333700       1.99511300       3.59589200
     H       1.01251200       2.23480200       2.64792500
     H       1.00016500       1.02441800       3.73277600
     O       2.83697500       1.01430500       0.96475500
     H       2.35001800       1.86503900       0.96632300
     H       3.80955300       1.14961100       0.96427700
     O       3.01712100       4.36370600       3.89516700
     H       3.92366000       4.01359800       4.02358900
     H       2.40561000       3.67751500       3.54960800
     O       4.13810600       3.33002100       1.10528600
     H       3.85585300       3.08176500       2.01059300
     H       4.29633900       4.29567000       1.02343100
VEC1       5.30828236       0.00000000       0.00000000
VEC2       0.00000000       5.30828236       0.00000000
VEC3       0.00000000       0.00000000       5.30828236
"""

    def test_packmol(self, five_water, plams_mols):
        # Packmol requires calling AMS, so mock out the call to the executable and return an expected xyz
        with get_mock_open_function(lambda f: f.endswith("output.xyz"), five_water):
            with patch("os.path.exists", return_value=True):
                with patch("scm.plams.interfaces.molecule.packmol.saferun") as mock_saferun:
                    input_mol = plams_mols["water"]

                    # Packmol has a lot of internal function calls, so even returning successfully is pretty good
                    # But also sanity-check the temporary input files passed to AMS
                    def read_input_file(*args, **kwargs):
                        # Intercept the tmp input files and check they are as expected
                        input_file = kwargs["stdin"]
                        content = input_file.read()
                        assert int(re.search(r"number (\d+)", content).group(1)) == 5

                        structure_file_name = re.search(r"structure\s+(.+\.xyz)", content).group(1)
                        with open(structure_file_name) as structure_file:
                            s = structure_file.read()
                        assert (
                            s
                            == """3

         O       0.06692105       0.06692105       0.00000000
         H       1.01204160      -0.07896265       0.00000000
         H      -0.07896265       1.01204160       0.00000000
"""
                        )

                    mock_saferun.side_effect = read_input_file
                    mols = packmol(input_mol, n_atoms=15, density=1.0, executable="mock_exe")

                    # And assert on the number of returned molecules
                    assert len(mols) == 15

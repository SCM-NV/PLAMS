from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.interfaces.molecule.ase import toASE
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ase_calculator import AMSCalculator
from test_helpers import skip_if_no_ams_installation


def test_properties():
    s = Settings()
    s.input.ForceField
    s.input.ams.Task = "SinglePoint"
    s.runscript.nproc = 1
    job = AMSCalculator(s, name="Properties")
    assert "forces" not in job.implemented_properties
    job.ensure_property("forces")
    assert "forces" in job.implemented_properties


def test_ase_deepcopy_worker():
    # PLAMS Settings configuring the calculation
    skip_if_no_ams_installation()

    settings = Settings()
    settings.input.ams.Task = "SinglePoint"
    settings.input.ams.Properties.Gradients = "Yes"
    settings.input.ams.Properties.StressTensor = "Yes"
    settings.input.ForceField.Type = "UFF"
    settings.runscript.nproc = 1

    def get_atoms():
        mol = from_smiles("O")  # PLAMS Molecule
        mol.lattice = [
            [
                3.0,
                0.0,
                0.0,
            ],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 3.0],
        ]
        return toASE(mol)  # convert PLAMS Molecule to ASE Atoms

    # set up atoms for water in a box
    atoms1 = get_atoms()
    atoms2 = get_atoms()
    # set OH bonds for atoms1 and atoms2 to 0.9 and 1.0 A respectively
    atoms1.set_distance(0, 1, 0.9, fix=0)
    atoms2.set_distance(0, 1, 1.0, fix=0)
    atoms1.set_distance(0, 2, 0.9, fix=0)
    atoms2.set_distance(0, 2, 1.0, fix=0)

    with AMSCalculator(settings=settings, name="ASE_deepcopy", amsworker=True) as calc:
        from copy import deepcopy

        atoms1.calc = deepcopy(calc)
        atoms2.calc = deepcopy(calc)
        e1 = atoms1.get_potential_energy()
        e2 = atoms2.get_potential_energy()
        # if no deepcopy is made, this would result in a new calculation
        e3 = atoms1.get_potential_energy()
        assert calc._counter[calc.name] == 2
        assert e1 == e3
        assert e2 > e1

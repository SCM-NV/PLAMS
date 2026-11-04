from collections import OrderedDict
from typing import Any, Sequence, Set

from scm.plams.core.basejob import MultiJob
from scm.plams.core.functions import add_to_instance
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.mol.atom import Atom
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.units import Units
from scm.plams.recipes._utils import Structure, is_chemical_system, to_plams_settings

__all__ = ["CounterpoiseEnergyJob", "CounterpoiseEnergyResults"]


class CounterpoiseEnergyResults(Results):
    """Results class for counterpoise corrected interaction and binding energies."""

    def get_all_energies(self, unit="au"):
        """'A_geo_AB_basis_AB' should be read as "atoms A at the AB geometry with the AB basis set"
        BSSE is the basis set superposition error
        Eint is the interaction energy
        Ebind is the binding energy (the complex relative to relaxed and isolated A and B)
        _raw means no counterpoise correction
        _cp means counterpoise-corrected
        Edef is the deformation energy
        """
        E = {}

        E["AB_geo_AB_basis_AB"] = self.job.AB.results.get_energy(unit=unit)
        E["A_geo_A_basis_A"] = (
            Units.convert(self.job.A, "hartree", unit)
            if isinstance(self.job.A, float)
            else self.job.children["A_geo_A_basis_A"].results.get_energy(unit=unit)
        )
        E["B_geo_B_basis_B"] = (
            Units.convert(self.job.B, "hartree", unit)
            if isinstance(self.job.B, float)
            else self.job.children["B_geo_B_basis_B"].results.get_energy(unit=unit)
        )

        E["A_geo_AB_basis_AB"] = self.job.children["A_geo_AB_basis_AB"].results.get_energy(unit=unit)
        E["A_geo_AB_basis_A"] = self.job.children["A_geo_AB_basis_A"].results.get_energy(unit=unit)
        E["B_geo_AB_basis_AB"] = self.job.children["B_geo_AB_basis_AB"].results.get_energy(unit=unit)
        E["B_geo_AB_basis_B"] = self.job.children["B_geo_AB_basis_B"].results.get_energy(unit=unit)

        E["BSSE_A"] = E["A_geo_AB_basis_AB"] - E["A_geo_AB_basis_A"]
        E["BSSE_B"] = E["B_geo_AB_basis_AB"] - E["B_geo_AB_basis_B"]

        E["BSSE_tot"] = E["BSSE_A"] + E["BSSE_B"]

        E["Eint_raw"] = E["AB_geo_AB_basis_AB"] - E["A_geo_AB_basis_A"] - E["B_geo_AB_basis_B"]
        E["Eint_cp"] = E["AB_geo_AB_basis_AB"] - E["A_geo_AB_basis_AB"] - E["B_geo_AB_basis_AB"]

        E["Edef_A"] = E["A_geo_AB_basis_A"] - E["A_geo_A_basis_A"]
        E["Edef_B"] = E["B_geo_AB_basis_B"] - E["B_geo_B_basis_B"]

        E["Ebind_raw"] = E["AB_geo_AB_basis_AB"] - E["A_geo_A_basis_A"] - E["B_geo_B_basis_B"]
        E["Ebind_cp"] = E["Eint_cp"] + E["Edef_A"] + E["Edef_B"]

        return E


class CounterpoiseEnergyJob(MultiJob):
    """A class for calculating the counterpoise-corrected interaction and binding energies."""

    _result_type = CounterpoiseEnergyResults

    def __init__(
        self,
        AB,
        ids_state_A,
        settings_state_AB=None,
        A=None,
        B=None,
        settings_state_A=None,
        settings_state_B=None,
        **kwargs,
    ):
        """
        AB is a Molecule, ChemicalSystem, or a previously run AMSJob.
        If it is a structure, a job will be created with settings_state_AB (default task Geometry Optimization).
        If it is an AMSJob, the structure and energy will be read from that job. settings_state_AB will be ignored.

        ids_state_A is a list of one-based atom indices making up the first molecule A.
        All other atoms are assigned to belong to the second molecule B.

        Set A and B to the energies (in Hartree) of the isolated fragments.
        Alternatively, set A and B to AMSJobs instances or Molecules of previously run geometry optimizations for the isolated molecules, which will then be taken as starting points for the geometry optimization. If they are None, the structure from the complex AB will be taken as a starting point for the geometry optimization.

        settings_state_A: PLAMS Settings or a complete scm.inputs.AMS model for A. Defaults to settings_state_AB.
        settings_state_B: PLAMS Settings or a complete scm.inputs.AMS model for B. Defaults to settings_state_AB.

        The calculations involving ghost atoms are enforced to be single point calculations.

        kwargs: other options to be passed to the MultiJob constructor (for example the name)
        """
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)

        self.A = A
        self.B = B
        self.AB = AB

        recipe_name = type(self).__name__
        settings_state_AB = to_plams_settings(settings_state_AB, recipe_name)
        settings_state_A = to_plams_settings(settings_state_A, recipe_name)
        settings_state_B = to_plams_settings(settings_state_B, recipe_name)

        # copy the settings so that we wont modify the ones provided as input by the user
        settings_state_AB = settings_state_AB.copy() if settings_state_AB else Settings()
        if settings_state_A:
            settings_state_A = settings_state_A.copy()
        else:
            settings_state_A = settings_state_AB.copy()
            settings_state_A.input.ams.Task = "GeometryOptimization"
        if settings_state_B:
            settings_state_B = settings_state_B.copy()
        else:
            settings_state_B = settings_state_AB.copy()
            settings_state_B.input.ams.Task = "GeometryOptimization"

        # In case the charge key is not specified, explicitly set the value to 0.
        # This is to prevent the charge in molecule.properties.charge (set by get_main_molecule())
        # to be used in case of neutral systems
        for s in [settings_state_AB, settings_state_A, settings_state_B]:
            if not "charge" in s.input.ams.system:
                s.input.ams.system.charge = 0

        use_chemical_system = is_chemical_system(self.AB) or (
            isinstance(self.AB, AMSJob) and is_chemical_system(self.AB.molecule)
        )
        if isinstance(self.AB, Molecule) or is_chemical_system(self.AB):
            self.children["AB_geo_AB_basis_AB"] = AMSJob(
                molecule=self.AB, settings=settings_state_AB, name="AB_geo_AB_basis_AB"
            )
            main_child = self.children["AB_geo_AB_basis_AB"]
            self.AB = main_child
        elif isinstance(AB, AMSJob):
            main_child = self.AB
        else:
            raise ValueError(
                "Unknown type for argument AB: {} (must be Molecule, ChemicalSystem, or AMSJob)".format(type(self.AB))
            )

        def get_main_structure() -> Any:
            getter = main_child.results.get_main_system if use_chemical_system else main_child.results.get_main_molecule
            return getter()

        def set_ghost_atoms(structure: Structure, ghost_indices: Set[int]) -> None:
            if is_chemical_system(structure):
                for index, atom in enumerate(structure.atoms, start=1):
                    atom.is_ghost = index in ghost_indices
            else:
                for index in range(1, len(structure) + 1):
                    structure[index].properties.ghost = index in ghost_indices

        def extract_atoms(structure: Structure, indices: Sequence[int]) -> Structure:
            if is_chemical_system(structure):
                return structure.extract_atoms([index - 1 for index in indices])
            extracted = Molecule()
            extracted.lattice = structure.lattice
            for index in indices:
                atom = structure[index]
                extracted.add_atom(Atom(atnum=atom.atnum, coords=atom.coords))
            return extracted

        name = "A_geo_AB_basis_AB"
        s = settings_state_A.copy()
        s.input.ams.Task = "SinglePoint"
        job = AMSJob(settings=s, name=name)

        @add_to_instance(job)
        def prerun(self):  # noqa F811
            structure = get_main_structure()
            set_ghost_atoms(structure, set(range(1, len(structure) + 1)) - set(ids_state_A))
            self.molecule = structure

        self.children[name] = job

        name = "A_geo_AB_basis_A"
        s = settings_state_A.copy()
        s.input.ams.Task = "SinglePoint"
        job = AMSJob(settings=s, name=name)

        @add_to_instance(job)
        def prerun(self):  # noqa F811
            structure = get_main_structure()
            self.molecule = extract_atoms(structure, ids_state_A)

        self.children[name] = job

        name = "B_geo_AB_basis_AB"
        s = settings_state_B.copy()
        s.input.ams.Task = "SinglePoint"
        job = AMSJob(settings=s, name=name)

        @add_to_instance(job)
        def prerun(self):  # noqa F811
            structure = get_main_structure()
            set_ghost_atoms(structure, set(ids_state_A))
            self.molecule = structure

        self.children[name] = job

        name = "B_geo_AB_basis_B"
        s = settings_state_B.copy()
        s.input.ams.Task = "SinglePoint"
        job = AMSJob(settings=s, name=name)

        @add_to_instance(job)
        def prerun(self):  # noqa F811
            structure = get_main_structure()
            indices = [index for index in range(1, len(structure) + 1) if index not in ids_state_A]
            self.molecule = extract_atoms(structure, indices)

        self.children[name] = job

        def add_isolated_child(X, name, settings):
            s = settings.copy()
            job = AMSJob(settings=s, name=name)
            if isinstance(X, Molecule) or is_chemical_system(X):
                job.molecule = X
            else:

                @add_to_instance(job)
                def prerun(self):  # noqa F811
                    structure = get_main_structure()
                    indices = [
                        index
                        for index in range(1, len(structure) + 1)
                        if (name.startswith("A") and index in ids_state_A)
                        or (name.startswith("B") and index not in ids_state_A)
                    ]
                    self.molecule = extract_atoms(structure, indices)

            self.children[name] = job

        if not isinstance(self.A, float):
            add_isolated_child(self.A, "A_geo_A_basis_A", settings_state_A)
        if not isinstance(self.B, float):
            add_isolated_child(self.B, "B_geo_B_basis_B", settings_state_B)

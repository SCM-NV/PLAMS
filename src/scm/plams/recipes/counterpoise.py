from collections import OrderedDict

from scm.plams.core.basejob import MultiJob
from scm.plams.core.functions import add_to_instance
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.mol.atom import Atom
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.units import Units

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
        AB is either a Molecule or a previously run AMSJob (from which the molecule will be extracted with .results.get_main_molecule())
        If it is a Molecule, a job will be created with settings_state_AB (default task Geometry Optimization).
        If it is an AMSJob, the structure and energy will be read from that job. settings_state_AB will be ignored.

        ids_state_A is a list of one-based atom indices making up the first molecule A.
        All other atoms are assigned to belong to the second molecule B.

        Set A and B to the energies (in Hartree) of the isolated fragments.
        Alternatively, set A and B to AMSJobs instances or Molecules of previously run geometry optimizations for the isolated molecules, which will then be taken as starting points for the geometry optimization. If they are None, the structure from the complex AB will be taken as a starting point for the geometry optimization.

        settings_state_A: Setting object for molecule A. Defaults to the same as settings_state_AB. The Task will be applied only to the isolated molecule A, if such a job is run.
        settings_state_B: Setting object for molecule B. Defaults to the same as settings_state_AB. The Task will be applied only to the isolated molecule B, if such job is run.

        The calculations involving ghost atoms are enforced to be single point calculations.

        kwargs: other options to be passed to the MultiJob constructor (for example the name)
        """
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)

        self.A = A
        self.B = B
        self.AB = AB

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

        if isinstance(self.AB, Molecule):
            self.children["AB_geo_AB_basis_AB"] = AMSJob(
                molecule=self.AB, settings=settings_state_AB, name="AB_geo_AB_basis_AB"
            )
            main_child = self.children["AB_geo_AB_basis_AB"]
            self.AB = main_child
        elif isinstance(AB, AMSJob):
            main_child = self.AB
        else:
            raise ValueError("Unknown type for argument AB: {} (must be Molecule or AMSJob)".format(type(self.AB)))

        name = "A_geo_AB_basis_AB"
        s = settings_state_A.copy()
        s.input.ams.Task = "SinglePoint"
        job = AMSJob(settings=s, name=name)

        @add_to_instance(job)
        def prerun(self):  # noqa F811
            mol = main_child.results.get_main_molecule()
            for i in range(1, len(mol) + 1):
                mol[i].properties.ghost = not i in ids_state_A
            self.molecule = mol

        self.children[name] = job

        name = "A_geo_AB_basis_A"
        s = settings_state_A.copy()
        s.input.ams.Task = "SinglePoint"
        job = AMSJob(settings=s, name=name)

        @add_to_instance(job)
        def prerun(self):  # noqa F811
            mainmol = main_child.results.get_main_molecule()
            mol = Molecule()
            mol.lattice = mainmol.lattice
            for i in range(1, len(mainmol) + 1):
                if i in ids_state_A:
                    mol.add_atom(Atom(atnum=mainmol[i].atnum, coords=mainmol[i].coords))
            self.molecule = mol

        self.children[name] = job

        name = "B_geo_AB_basis_AB"
        s = settings_state_B.copy()
        s.input.ams.Task = "SinglePoint"
        job = AMSJob(settings=s, name=name)

        @add_to_instance(job)
        def prerun(self):  # noqa F811
            mol = main_child.results.get_main_molecule()
            for i in range(1, len(mol) + 1):
                mol[i].properties.ghost = i in ids_state_A
            self.molecule = mol

        self.children[name] = job

        name = "B_geo_AB_basis_B"
        s = settings_state_B.copy()
        s.input.ams.Task = "SinglePoint"
        job = AMSJob(settings=s, name=name)

        @add_to_instance(job)
        def prerun(self):  # noqa F811
            mainmol = main_child.results.get_main_molecule()
            mol = Molecule()
            mol.lattice = mainmol.lattice
            for i in range(1, len(mainmol) + 1):
                if i not in ids_state_A:
                    mol.add_atom(Atom(atnum=mainmol[i].atnum, coords=mainmol[i].coords))
            self.molecule = mol

        self.children[name] = job

        def add_isolated_child(X, name, settings):
            s = settings.copy()
            job = AMSJob(settings=s, name=name)
            if isinstance(X, Molecule):
                job.molecule = X
            else:

                @add_to_instance(job)
                def prerun(self):  # noqa F811
                    mainmol = main_child.results.get_main_molecule()
                    mol = Molecule()
                    mol.lattice = mainmol.lattice
                    for i in range(1, len(mainmol) + 1):
                        if (name.startswith("A") and i in ids_state_A) or (
                            name.startswith("B") and i not in ids_state_A
                        ):
                            mol.add_atom(Atom(atnum=mainmol[i].atnum, coords=mainmol[i].coords))
                    self.molecule = mol

            self.children[name] = job

        if not isinstance(self.A, float):
            add_isolated_child(self.A, "A_geo_A_basis_A", settings_state_A)
        if not isinstance(self.B, float):
            add_isolated_child(self.B, "B_geo_B_basis_B", settings_state_B)

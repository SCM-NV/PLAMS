#!/usr/bin/env amspython
from abc import ABC, abstractmethod
from collections import defaultdict, namedtuple
from dataclasses import dataclass
from natsort import natsorted
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import scm.plams as plams
import scm.reactmap
import warnings

__all__ = ["SMILESReaction", "XYZReaction", "MoleculeReaction", "run_ts_workflow"]


@dataclass
class BoostReaction(ABC):
    @abstractmethod
    def get_molecules_dict(self) -> Dict[str, plams.Molecule]:
        pass

    @staticmethod
    def mapping(reactant: plams.Molecule, product: plams.Molecule) -> plams.Molecule:
        """returns a new product with rearranged atom indices"""
        settings = scm.reactmap.Settings(print_progress=True, preprocess_map=True)
        reaction = scm.reactmap.Reaction(
            reactant=scm.reactmap.Molecule(plams_mol=reactant), product=scm.reactmap.Molecule(plams_mol=product)
        )
        if not reaction.is_feasible():
            raise ValueError(f"Unfeasible reaction: {reaction}")
        scm.reactmap.Map(reaction, settings=settings)

        plams.log("SCM Reaction Mapping enabled!")
        plams.log(reaction.mapping)
        for i, x in enumerate(reaction.mapping):
            if i != x:
                plams.log(f"Reactant atom {i+1} -> product atom {x+1}")

        temp_mol = product.get_fragment(reaction.mapping)  # get_fragment uses zero-based indices!
        return temp_mol


def reaction_examples(i: int) -> BoostReaction:
    reactions = {
        0: SMILESReaction(
            reactants_smiles=["C1=CC(=O)C(=O)C=C1", "[OH3+]"],
            # products_smiles = ["C1=C[C]([OH+])C(=O)C=C1", "O"],
            products_smiles=["C1=CC(=[OH+])C(=O)C=C1", "O"],
        ),
        1: SMILESReaction(
            reactants_smiles=["C1=CC(=[OH+])C(=O)C=C1", "O"],
            products_smiles=["C1=CC(O)C(=[OH+])C(O)=C1"],
        ),
        2: SMILESReaction(
            reactants_smiles=["C=C", "O"],
            products_smiles=["CCO"],
        ),
        3: SMILESReaction(
            reactants_smiles=["C1C=CC=C1", "C1=CC(=O)OC1=O"],
            products_smiles=["[H][C@@]12C3CC(C=C3)[C@]1([H])C(=O)OC2=O"],
        ),
    }
    return reactions[i]


def uff_settings() -> plams.Settings:
    """returns UFF engine settings"""
    s = plams.Settings()
    s.input.ForceField.Type = "UFF"
    s.runscript.nproc = 1
    return s


def run_ts_workflow(
    reaction: BoostReaction,
    engine_settings: plams.Settings,
    reaction_boost_equilibration_nsteps: int = 1000,
    reaction_boost_nsteps: int = 500,
    reaction_boost_strength: float = 0.7,
    reaction_boost_equilibration_uff: bool = True,
    reaction_boost_equilibration_strength: float = 0.05,
    do_ts_search: bool = True,
    do_irc: bool = True,
    run_equilibration: bool = True,
    temperature: float = 200,
):
    """Runs all the steps"""

    warnings.filterwarnings("ignore", "biadjacency_matrix will return a scipy.sparse array instead of a matrix")

    if reaction_boost_equilibration_uff:
        equilibration_engine_settings = uff_settings()
    else:
        equilibration_engine_settings = engine_settings.copy()

    molecules_dict = reaction.get_molecules_dict()
    if "final" in molecules_dict and False:
        molecules_dict["final"] = plams.preoptimize(molecules_dict["final"], settings=engine_settings, maxiterations=30)

    boost_job = run_reaction(
        molecules_dict,
        engine_settings=engine_settings,
        run_equilibration=run_equilibration,
        equilibration_nsteps=reaction_boost_equilibration_nsteps,
        equilibration_engine_settings=equilibration_engine_settings,
        equilibration_strength=reaction_boost_equilibration_strength,
        preoptimize_after_equilibration=True,
        nsteps=reaction_boost_nsteps,
        strength=reaction_boost_strength,
        temperature=temperature,
    )

    md_forward_barrier, md_backward_barrier = boost_job.get_barriers()
    plams.log(f"Forward MD: {md_forward_barrier*27.211:.3f}")
    plams.log(f"Backward MD: {md_backward_barrier*27.211:.3f}")

    ts_job = run_ts(
        boost_job,
        engine_settings=engine_settings,
    )

    if not ts_job.ok():
        plams.log("No TS found, perhaps barrierless?!")
        return

    irc_job = get_irc_job(
        ts_job,
        engine_settings=engine_settings,
    )
    irc_job.run()

    if not irc_job.ok():
        plams.log("IRC job failed!")
        return

    irc_results = irc_job.results.get_irc_results()
    irc_forward_barrier = irc_results["LeftBarrier"]
    irc_backward_barrier = irc_results["RightBarrier"]

    plams.log(f"Forward IRC: {irc_backward_barrier*27.211:.3f}")
    plams.log(f"Backward IRC: {irc_forward_barrier*27.211:.3f}")


@dataclass
class SMILESReaction(BoostReaction):
    reactants_smiles: List[str]
    products_smiles: List[str]
    name: str = "reaction"
    fix_stoichiometry: bool = True

    def get_molecules(self) -> Tuple[plams.Molecule, plams.Molecule]:
        if self.fix_stoichiometry:
            self.add_molecules_if_needed()

        reactants = self.get_combined(self.reactants_smiles)
        products = self.get_combined(self.products_smiles)
        products = self.mapping(reactants, products)
        return reactants, products

    def _fix_single_stoichoimetry(self, stoich: Dict[str, int]) -> bool:
        """Returns True if the stoichiometry was modified"""
        A = namedtuple("A", ["stoich", "smiles"])
        fixes = [
            A(stoich={"H": 2, "O": 1}, smiles="O"),
            A(stoich={"H": 1, "Cl": 1}, smiles="Cl"),
            A(stoich={"H": 1, "F": 1}, smiles="F"),
            A(stoich={"H": 1, "Br": 1}, smiles="Br"),
            A(stoich={"H": 2, "S": 1}, smiles="S"),
        ]

        for fix in fixes:
            if all(stoich[x] <= -y for x, y in fix.stoich.items()):
                self.products_smiles.append(fix.smiles)
                for k, v in fix.stoich.items():
                    stoich[k] += v
                return True
            if all(stoich[x] >= y for x, y in fix.stoich.items()):
                self.reactants_smiles.append(fix.smiles)
                for k, v in fix.stoich.items():
                    stoich[k] -= v
                return True
        return False

    def add_molecules_if_needed(self):
        """modifies self.reactants_smiles and self.products_smiles"""
        rmols = [plams.from_smiles(x) for x in self.reactants_smiles]
        pmols = [plams.from_smiles(y) for y in self.products_smiles]
        stoich = defaultdict(lambda: 0)
        for mol in rmols:
            for k, v in mol.get_formula(as_dict=True).items():
                stoich[k] -= v
        for mol in pmols:
            for k, v in mol.get_formula(as_dict=True).items():
                stoich[k] += v

        while self._fix_single_stoichoimetry(stoich):
            pass

    def get_molecules_dict(self) -> Dict[str, plams.Molecule]:
        r, p = self.get_molecules()
        ret = {"": r, "final": p}
        return ret

    @classmethod
    def get_combined(cls, smiles_list: List[str], margin: float = 1.0) -> plams.Molecule:
        ret = plams.Molecule()
        molecules = [plams.from_smiles(x) for x in smiles_list]
        tot_charge = sum(mol.properties.get("charge", 0) for mol in molecules)
        molecules = plams.preoptimize(molecules)
        for m in molecules:
            ret.add_molecule(m, margin=margin)
        ret.properties.charge = tot_charge
        return ret


@dataclass
class MoleculeReaction(BoostReaction):
    reactant: plams.Molecule
    product: plams.Molecule
    do_mapping: bool = False

    def get_molecules_dict(self) -> Dict[str, plams.Molecule]:
        """Returns a molecule dictionary"""
        if self.do_mapping:
            product = self.mapping(self.reactant, self.product)
        else:
            product = self.product
        ret = {"": self.reactant, "final": self.mapping(self.reactant, product)}
        return ret


@dataclass
class XYZReaction(BoostReaction):
    folder: Path  # a folder with .xyz files
    charge: float = 0

    def get_molecules_dict(self) -> Dict[str, plams.Molecule]:
        """
        read molecules from .xyz files in a directory ``self.folder``.
        Returns a dictionary suitable for AMSJob
        """
        # read all xyz files, dictionary key: Molecule
        molecules = plams.read_molecules(str(self.folder.resolve()))

        sorted_keys = list(natsorted(molecules.keys()))
        molecules[""] = molecules.pop(sorted_keys[0])  # the first molecule must have key ''
        sorted_keys.pop(0)

        if self.charge != 0:
            for mol in molecules.values():
                mol.properties.charge = self.charge

        # molecules[""].hydrogen_to_deuterium()  # slow down hydrogen motions
        # molecules[''].lattice = [[9, 0, 0], [0, 9, 0], [0, 0, 9]]  # optional, add lattice

        return molecules


class AMSReactionBoostJob(plams.AMSMDJob):
    def __init__(
        self,
        molecule,
        settings: plams.Settings = None,
        nsteps_per: int = None,
        nsteps: int = None,
        strength: float = 1.0,
        initial_fraction: float = 0.05,
        post_reaction_relaxation_steps: int = 0,
        **kwargs,
    ):
        if isinstance(molecule, plams.Molecule):
            molecule = {"": molecule}
        assert isinstance(molecule, dict)
        assert len(molecule) >= 2

        n_transitions = len(molecule) - 1
        if nsteps_per is None and nsteps is not None:
            nsteps_per = round(nsteps / n_transitions)
        elif nsteps is None and nsteps_per is not None:
            nsteps = nsteps_per * n_transitions + post_reaction_relaxation_steps

        settings = settings or dict()

        super().__init__(molecule=molecule, settings=settings, nsteps=nsteps, **kwargs)
        self.nsteps_per = nsteps_per
        self.strength = strength
        self.post_reaction_relaxation_steps = post_reaction_relaxation_steps
        self.n_transitions = n_transitions
        self.initial_fraction = initial_fraction

        self.settings += self.get_boost_settings(
            what="pair", strength=self.strength, initial_fraction=self.initial_fraction
        )

        self.sorted_keys = [x for x in molecule if x != ""]
        self.settings.input.ams.MolecularDynamics.ReactionBoost.TargetSystem = self.sorted_keys
        self.settings.input.ams.MolecularDynamics.ReactionBoost.NSteps = self.nsteps_per

    @classmethod
    def get_boost_settings(cls, what: str, strength: float = 1.0, initial_fraction: float = 0.05) -> plams.Settings:
        """Returns the bulk of the MolecularDynamics%ReactionBoost settings"""
        s = plams.Settings()
        s.input.ams.Log.Info = "ReactionBoost"
        if what.lower() == "rmsd":
            s.input.ams.MolecularDynamics.ReactionBoost.Type = "RMSD"
            s.input.ams.MolecularDynamics.ReactionBoost.RMSDRestraint.Type = "GaussianWell"
            # s.input.ams.MolecularDynamics.ReactionBoost.Region = "solute"
            s.input.ams.MolecularDynamics.ReactionBoost.Change = "TargetCoordinate"
            # s.input.ams.MolecularDynamics.ReactionBoost.InitialFraction = 0.01
            # s.input.ams.MolecularDynamics.ReactionBoost.RMSDRestraint.Erf.MaxForce = 0.1
        else:
            s.input.ams.MolecularDynamics.ReactionBoost.Type = "Pair"
            rb = s.input.ams.MolecularDynamics.ReactionBoost
            rb.Change = "LogForce"
            rb.InitialFraction = initial_fraction  # at the first time step

            rb.BondBreakingRestraints.Type = "Erf"  # "None" to disable
            # rb.BondBreakingRestraints.Type = "None"  # "None" to disable
            rb.BondBreakingRestraints.Erf.MaxForce = 0.05 * strength
            # rb.BondBreakingRestraints.Erf.ForceConstant = 0.1

            rb.BondMakingRestraints.Type = "Erf"  # "None" to disable
            rb.BondMakingRestraints.Erf.MaxForce = 0.2 * strength
            # rb.BondMakingRestraints.Erf.ForceConstant = 0.1

            rb.BondedRestraints.Type = "Harmonic"  # "None" to disable
            rb.BondedRestraints.Harmonic.ForceConstant = 0.001 * strength

            rb.NonBondedRestraints.Type = "Exponential"  # "None" to disable
            rb.NonBondedRestraints.Exponential.Epsilon = 1e-4 * strength

        return s

    def get_approximate_ts_index(self) -> int:
        """Returns 1-based index for frame with highest engine energy"""
        engine_energy_hartree = self.results.get_history_property("EngineEnergy", "History")
        index = np.argmax(engine_energy_hartree) + 1

        return index

    def get_barriers(self, index: int = None) -> Tuple[float, float]:
        """Returns the forward and backward potential energy barriers in hartree relative to the frame with index=index"""
        if index is None:
            index = self.get_approximate_ts_index()

        engine_energy_hartree = self.results.get_history_property("EngineEnergy")
        md_forward_barrier = np.max(engine_energy_hartree) - np.min(engine_energy_hartree[:index])
        md_backward_barrier = np.max(engine_energy_hartree) - np.min(engine_energy_hartree[index - 1 :])

        return md_forward_barrier, md_backward_barrier

    def get_reaction_coordinates(
        self,
    ) -> Tuple[
        List[Tuple[int, int, float]],  # bond making
        List[Tuple[int, int, float]],  # bond breaking
        List[Tuple[int, int, float]],  # bonded
    ]:
        """
        Reads the stdout log info and finds out which restraints were added

        NOTE: This assumes that n_transitions = 1!
        """

        bond_making_lines = self.results.grep_output("bond-making")
        bond_breaking_lines = self.results.grep_output("bond-breaking")
        bonded_lines = self.results.grep_output("bonded")

        bond_making_tuples = []
        bond_breaking_tuples = []
        bonded_tuples = []

        for line in bond_making_lines:
            splitline = line.split()
            target_value = float(splitline[5])  # bohr
            if target_value > 5:
                continue
            bond_making_tuples.append((int(splitline[2]), int(splitline[3]), 1.0))

        for line in bond_breaking_lines:
            splitline = line.split()  # missing space so use other indices
            initial_value = float(splitline[3])  # bohr
            if initial_value > 5:
                continue
            bond_breaking_tuples.append((int(splitline[1]), int(splitline[2]), -1.0))

        for line in bonded_lines:
            splitline = line.split()
            mean = 0.5 * (float(splitline[4]) + float(splitline[5])) * 0.529
            bonded_tuples.append((int(splitline[2]), int(splitline[3]), mean))

        return bond_making_tuples, bond_breaking_tuples, bonded_tuples

    def get_rc_atoms(self) -> List[int]:
        """
        Returns atom indices that are involved in bond-making or bond-breaking restraints
        """

        bm, bb, _ = self.get_reaction_coordinates()
        ret = set()
        for b in bm + bb:
            ret.add(b[0])
            ret.add(b[1])
        return list(ret)

    def get_reaction_coordinate_lines(self) -> Tuple[List[str], List[str], List[str]]:
        """
        Returns 3-tuple list of str: bond-making, bond-breaking, bonded
        Can be used e.g. to set
        s.input.ams.TransitionStateSearch.ReactionCoordinate.Distance = (
            bond_making_lines + bond_breaking_lines
        )
        s.input.ams.Constraints.Distance = bonded_lines
        """

        def rc_tuple_to_str_list(rc_tuple):
            if rc_tuple is None:
                return []
            ret = []
            for t in rc_tuple:
                ret.append(" ".join(str(x) for x in t))

            return ret

        bond_making, bond_breaking, bonded = self.get_reaction_coordinates()

        bond_making_lines = rc_tuple_to_str_list(bond_making)
        bond_breaking_lines = rc_tuple_to_str_list(bond_breaking)
        bonded_lines = rc_tuple_to_str_list(bonded)

        return bond_making_lines, bond_breaking_lines, bonded_lines


def run_reaction(
    molecules_dict: Dict[str, plams.Molecule],
    engine_settings: plams.Settings,
    run_equilibration: bool = True,
    equilibration_nsteps: int = 500,
    equilibration_strength: float = 0.05,
    equilibration_engine_settings: plams.Settings = None,
    nsteps: int = 1000,
    strength: float = 1.0,
    preoptimize_after_equilibration: bool = True,
    temperature: float = 200,
) -> AMSReactionBoostJob:
    """Runs equilibration and produciton boost MD. Returns the production job"""

    for i, at in enumerate(molecules_dict[""], 1):
        if (
            "Constraints" in engine_settings.input.ams
            and "Atom" in engine_settings.input.ams.Constraints
            and i not in engine_settings.input.ams.Constraints.Atom
        ):
            plams.AMSJob._add_region(at, "thermostatted")
        else:
            plams.AMSJob._add_region(at, "frozen")

    if run_equilibration:
        equilibration_engine_settings = equilibration_engine_settings or engine_settings.copy()
        eq_job = AMSReactionBoostJob(
            settings=equilibration_engine_settings,
            molecule=molecules_dict,
            name="eq",
            nsteps=equilibration_nsteps,
            strength=equilibration_strength,
            temperature=300,
            samplingfreq=20,
            tau=5,  # very short time constant (fs) to cool the system down
            timestep=1.0,  # fs
            thermostat="NHC",
            thermostat_region="thermostatted",
        )
        eq_job.run()

        # update the initial molecule
        molecules_dict[""] = eq_job.results.get_main_molecule()

        if preoptimize_after_equilibration:
            rc_atoms = eq_job.get_rc_atoms()
            preoptimization_job = run_preoptimization_job_fixing_active_atoms(
                engine_settings=engine_settings,
                molecule=molecules_dict[""],
                rc_atoms=rc_atoms,
                name="preoptimization_prod",
            )

            molecules_dict[""] = preoptimization_job.results.get_main_molecule()

    prod_job = AMSReactionBoostJob(
        settings=engine_settings,
        molecule=molecules_dict,
        name="prod",
        nsteps=nsteps,
        strength=strength,
        post_reaction_relaxation_steps=20,
        temperature=temperature,
        samplingfreq=10,
        tau=3.0,  # very short time constant (fs) to cool the system down
        timestep=0.5,
        thermostat="NHC",
        writeenginegradients=True,
        thermostat_region="thermostatted",
    )

    prod_job.run()

    return prod_job


def run_preoptimization_job_fixing_active_atoms(
    rc_atoms: List[int], molecule: plams.Molecule, engine_settings: plams.Settings, name: str = "preoptimization"
) -> plams.AMSJob:

    for at_i in rc_atoms:
        plams.AMSJob._add_region(molecule[at_i], "bondchange")

    s = plams.Settings()
    s += engine_settings
    s.input.ams.GeometryOptimization.Method = "Quasi-Newton"
    s.input.ams.GeometryOptimization.CoordinateType = "Auto"
    s.input.ams.GeometryOptimization.MaxIterations = 500
    s.input.ams.GeometryOptimization.PretendConverged = "Yes"

    s.input.ams.Task = "GeometryOptimization"
    s.input.ams.GeometryOptimization.Convergence.Quality = "Basic"
    s.input.ams.Constraints.Atom = rc_atoms
    preoptimization_job = plams.AMSJob(
        settings=s,
        molecule=molecule,
        name=name,
    )
    preoptimization_job.run()

    return preoptimization_job


def run_ts(
    boost_job: AMSReactionBoostJob,
    engine_settings: plams.Settings,
    coordinate_type: str = "Auto",
    preoptimization: bool = True,
) -> plams.AMSJob:
    """Returns a transition state search job based on a preiovus AMSReactionBoostJob"""

    index = boost_job.get_approximate_ts_index()
    approximate_ts_molecule = boost_job.results.get_history_molecule(index)

    if isinstance(boost_job.molecule, dict) and len(boost_job.molecule) == 2:
        bond_making_lines, bond_breaking_lines, bonded_lines = boost_job.get_reaction_coordinate_lines()
        rc_atoms = boost_job.get_rc_atoms()
        bonded_lines = []
    else:
        bond_making_lines, bond_breaking_lines, bonded_lines = [], [], []
        rc_atoms = []

    use_reaction_coordinates = False

    s = plams.Settings()
    s += engine_settings
    s.input.ams.GeometryOptimization.Method = "Quasi-Newton"
    s.input.ams.GeometryOptimization.CoordinateType = coordinate_type
    s.input.ams.GeometryOptimization.MaxIterations = 500
    s.input.ams.GeometryOptimization.PretendConverged = "Yes"

    if preoptimization and rc_atoms:
        preoptimization_job = run_preoptimization_job_fixing_active_atoms(
            molecule=approximate_ts_molecule,
            name="preoptimization_ts",
            engine_settings=engine_settings,
            rc_atoms=rc_atoms,
        )
        approximate_ts_molecule = preoptimization_job.results.get_main_molecule()

    s.input.ams.Task = "TransitionStateSearch"
    s.input.ams.GeometryOptimization.InitialHessian.Type = "Calculate"
    s.input.ams.GeometryOptimization.MaxIterations = 1000
    # s.input.ams.GeometryOptimization.Convergence.Gradients = 1e-4   # hartree/angstrom
    # s.input.ams.GeometryOptimization.Convergence.Energy = 1e-5   # hartree
    # s.input.ams.GeometryOptimization.Convergence.Quality = "Basic"

    if bond_making_lines or bond_breaking_lines and use_reaction_coordinates:
        s.input.ams.TransitionStateSearch.ReactionCoordinate.Distance = bond_making_lines + bond_breaking_lines

    if bonded_lines:
        s.input.ams.Constraints.Distance = bonded_lines

    s.input.ams.Properties.PESPointCharacter = "Yes"
    job = plams.AMSJob(settings=s, molecule=approximate_ts_molecule, name="ts")

    job.run()
    return job


def get_irc_job(ts_job: plams.AMSJob, engine_settings: plams.Settings) -> plams.AMSJob:
    s = plams.Settings()
    s += engine_settings
    s.input.ams.Task = "IRC"
    s.input.ams.IRC.MinEnergyProfile = "Yes"
    s.input.ams.IRC.Direction = "Both"
    s.input.ams.GeometryOptimization.PretendConverged = "Yes"
    s.input.ams.GeometryOptimization.Convergence.Quality = "Basic"

    job = plams.AMSJob(settings=s, molecule=ts_job.results.get_main_molecule(), name="irc")

    return job

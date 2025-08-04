from scm.plams.interfaces.molecule.rdkit import writepdb
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.periodic_table import PeriodicTable
from scm.plams.trajectories.rkfhistoryfile import molecules_to_rkf

try:
    from scm.plams.interfaces.molecule.ase import toASE
except ImportError:
    pass
import os
from typing import List, Tuple

from scm.plams.recipes.md.amsmdjob import AMSMDJob

"""


job = AMSMovingRestraintRMSDJob(
    visit_molecules=[mol1, mol2, mol3],
    kappa_scheduler=KappaExponentialScheduler(kappa=10000, n=4, base=10),
    nsteps=100000,
    ...
)

job = AMSMovingRestraintInternalJob(
    restraints=[
        BondRestraint(atom1, atom2, start=x, end=y, kappa_scheduler=KappaLinearScheduler(10000)),
        BondRestraint(atom1, atom2, start=x, end=y, kappa_scheduler=KappaLinearScheduler(10000)),
    ],
    consecutive=False
)

"""

__all__ = ["AMSMovingRestraintRMSDJob", "AMSMovingRestraintBondSwitchJob"]


def tolist(x):
    if not isinstance(x, list):
        x = [x]
    return x


def set_plumed_settings(settings, plumed_sett_list):
    """
    settings: Settings
        will be modified!

    plumed_sett_list: list of str
        The input to plumed. Should not include the final End
    """
    ret = ""
    indent = 6 * " "
    for line in plumed_sett_list:
        ret += f"{indent}{line}\n"
    ret += "    End"
    settings.input.ams.MolecularDynamics.Plumed.Input = ret


class AMSMovingRestraintBondSwitchJob(AMSMDJob):
    def __init__(
        self,
        visit_molecules: List[Molecule],
        kappa: List[float],
        start_step: int = 1,
        kappa_others: List[float] = None,
        **kwargs,
    ):
        AMSMDJob.__init__(self, **kwargs)
        self.visit_molecules = tolist(visit_molecules).copy()
        self.kappa = tolist(kappa)
        self.start_step = start_step
        self.nsteps = int(self.settings.input.ams.MolecularDynamics.NSteps)
        self.kappa_others = tolist(kappa_others) if kappa_others is not None else self.kappa

        previous_mol = self.molecule.copy()
        assert previous_mol is not None, "You must specify a molecule when initializing AMSMovingRestraintBondSwitchJob"

        if not previous_mol.bonds:
            previous_mol.guess_bonds()

        plumed_sett_list = [""]
        # plumed_sett_list.append('COORDINATION GROUPA=14 GROUPB=8,9,10,11,12,13,14,15 R_0=0.13 LABEL=c_14')
        # plumed_sett_list.append('PRINT ARG=c_14 FILE=c_14.txt STRIDE=1')
        counter = 0
        current_step = self.start_step
        nsteps_part = (self.nsteps - self.start_step) // len(self.visit_molecules)
        nsteps_per_kappa = nsteps_part // len(self.kappa)
        nsteps_per_kappa_others = nsteps_part // len(self.kappa_others)

        for counter, mol in enumerate(self.visit_molecules, 1):
            if not mol.bonds:
                mol.guess_bonds()

            bond_switches, other, nonbonded = self._get_bond_switches(previous_mol, mol)

            for oth in other:
                i1, i2, distance = oth
                label = f"bond{counter}_{i1}_{i2}"

                plumed_sett_list.append(f"DISTANCE ATOMS={i1},{i2} LABEL={label}")
                # now loop through all kappa values

                s = f"MOVINGRESTRAINT ARG={label}"
                for i, kap in enumerate(self.kappa_others, 0):
                    s += f" STEP{i}={current_step+nsteps_per_kappa_others*i} KAPPA{i}={kap} AT{i}={distance*0.1:.4f}"
                end_index = len(self.kappa_others)
                s += f" STEP{end_index}={current_step+nsteps_per_kappa_others*end_index} KAPPA{end_index}=0.0 AT{end_index}={distance*0.1:.4f}"

                plumed_sett_list.append(s)

            for bs in bond_switches:
                i1, i2, distance = bs
                label = f"bond{counter}_{i1}_{i2}"

                plumed_sett_list.append(f"DISTANCE ATOMS={i1},{i2} LABEL={label}")
                # now loop through all kappa values

                s = f"MOVINGRESTRAINT ARG={label}"
                for i, kap in enumerate(self.kappa, 0):
                    s += f" STEP{i}={current_step+nsteps_per_kappa*i} KAPPA{i}={kap} AT{i}={distance*0.1:.4f}"
                end_index = len(self.kappa)
                end_step = current_step + nsteps_per_kappa * end_index
                s += f" STEP{end_index}={end_step} KAPPA{end_index}=0.0 AT{end_index}={distance*0.1:.4f}"

                plumed_sett_list.append(s)

            for nb in nonbonded:
                i1, i2, distance = nb
                label = f"nonbond{counter}_{i1}_{i2}"

                plumed_sett_list.append(f"DISTANCE ATOMS={i1},{i2} LABEL={label}")
                # now loop through all kappa values

                s = f"MOVINGRESTRAINT ARG={label} VERSE=L"
                s += f" STEP0={current_step} KAPPA0=0 AT0={distance*0.1:.4f}"
                s += f" STEP1={end_step-1} KAPPA1=100000.0 AT1={distance*0.1:.4f}"
                s += f" STEP2={end_step} KAPPA2=0.0 AT2={distance*0.1:.4f}"

                plumed_sett_list.append(s)

            current_step += nsteps_per_kappa * end_index + 1

            previous_mol = mol

        set_plumed_settings(self.settings, plumed_sett_list)

    @staticmethod
    def _get_bond_switches(mol1: Molecule, mol2: Molecule) -> Tuple[List[Tuple[int, int, float]]]:
        assert len(mol1) == len(mol2), "mol1 and mol2 have different lengths!"
        bonds1 = set(tuple(sorted((mol1.index(b.atom1), mol1.index(b.atom2)))) for b in mol1.bonds)
        bonds2 = set(tuple(sorted((mol2.index(b.atom1), mol2.index(b.atom2)))) for b in mol2.bonds)

        all_reactive_atoms = set()

        atoms2 = toASE(mol2)

        # changed bonds
        targets = []
        t1 = bonds1 - bonds2
        t2 = bonds2 - bonds1
        for t in t1:
            d = atoms2.get_distance(t[0] - 1, t[1] - 1, mic=True)
            targets.append((t[0], t[1], d))
            all_reactive_atoms.add(t[0])
            all_reactive_atoms.add(t[1])
        for t in t2:
            d = atoms2.get_distance(t[0] - 1, t[1] - 1, mic=True)
            targets.append((t[0], t[1], d))
            all_reactive_atoms.add(t[0])
            all_reactive_atoms.add(t[1])

        all_non_reactive_atoms = set(range(1, len(mol1) + 1)) - all_reactive_atoms

        union_bonds = bonds1.union(bonds2)

        # other bonded
        other = []
        other_bonds = union_bonds - t1 - t2
        for ob in other_bonds:
            i1, i2 = ob
            # other.append( (i1, i2, 0.5*(mol1[i1].distance_to(mol1[i2]) + mol2[i1].distance_to(mol2[i2]))) )
            d = atoms2.get_distance(i1 - 1, i2 - 1, mic=True)
            other.append((i1, i2, d))

        # other nonbonded
        nonbonded = []
        covalent_radius_multiplier = 1.5

        for ireac in sorted(all_reactive_atoms):
            for inonreac in sorted(all_non_reactive_atoms):
                if (ireac, inonreac) in union_bonds or (inonreac, ireac) in union_bonds:
                    continue
                at1 = mol2[ireac]
                at2 = mol2[inonreac]
                r = PeriodicTable.get_radius(at1.symbol) + PeriodicTable.get_radius(at2.symbol)
                nonbonded.append((ireac, inonreac, covalent_radius_multiplier * r))
        # for i, at1 in enumerate(mol2, 1):
        #    for j, at2 in enumerate(mol2, 1):
        #        if j >= i:
        #            continue
        #        if (i, j) in bonds1 or (i, j) in bonds2:
        #            continue
        #        r = PeriodicTable.get_radius(at1.symbol) + PeriodicTable.get_radius(at2.symbol)
        #        nonbonded.append( (i, j, covalent_radius_multiplier*r) )

        return targets, other, nonbonded


class AMSMovingRestraintRMSDJob(AMSMDJob):
    def __init__(self, visit_molecules: List[Molecule], kappa: List[int], start_step: int = 1, **kwargs):
        AMSMDJob.__init__(self, **kwargs)
        self.visit_molecules = tolist(visit_molecules)
        self.kappa = tolist(kappa)
        self.start_step = start_step
        self.nsteps = int(self.settings.input.ams.MolecularDynamics.NSteps)
        self.pdb_files = self._get_pdb_file_names()

    def _get_pdb_file_names(self):
        return [f"mol{i}.pdb" for i in range(len(self.visit_molecules))]

    def prerun(self):  # noqa F811
        plumed_sett_list = [""]
        current_step = self.start_step
        nsteps_part = (self.nsteps - self.start_step) // len(self.visit_molecules)
        nsteps_per_kappa = nsteps_part // len(self.kappa)

        molecules_to_rkf(self.visit_molecules, os.path.join(self.path, "to_be_visited.rkf"))

        for pdb_rel, mol in zip(self.pdb_files, self.visit_molecules):
            pdb = os.path.join(self.path, pdb_rel)
            self._write_plumed_compatible_pdb(pdb, mol)
            name = os.path.splitext(pdb_rel)[0]

            label = f"rmsd_{name}"

            indent = 6 * " "
            plumed_sett_list.append(f"{indent}RMSD REFERENCE={pdb} LABEL={label}")
            # now loop through all kappa values

            s = f"{indent}MOVINGRESTRAINT ARG={label}"
            for i, kap in enumerate(self.kappa, 0):
                s += f" STEP{i}={current_step+nsteps_per_kappa*i} KAPPA{i}={kap} AT{i}=0.0"
            end_index = len(self.kappa)
            s += f" STEP{end_index}={current_step+nsteps_per_kappa*end_index} KAPPA{end_index}=0.0 AT{end_index}=0.0"

            current_step += nsteps_per_kappa * end_index + 1

            plumed_sett_list.append(s)

        plumed_sett_list.append("    End")

        self.settings.input.ams.MolecularDynamics.Plumed.Input = "\n".join(plumed_sett_list)

    @classmethod
    def _write_plumed_compatible_pdb(cls, pdb_file: str, molecule: Molecule):
        temporary_pdb = "temp.pdb"
        writepdb(molecule, temporary_pdb)
        char_index = 62
        try:
            xyzfile = os.path.splitext(pdb_file)[0] + ".xyz"
            molecule.write(xyzfile)
            with open(pdb_file, "w") as f:
                with open(temporary_pdb, "r") as temp:
                    for line in temp:
                        splitline = list(line)
                        if len(splitline) > char_index + 1 and splitline[char_index] == "0":
                            splitline[char_index] = "1"
                        f.write("".join(splitline) + "\n")
        finally:
            if os.path.exists(temporary_pdb):
                os.remove(temporary_pdb)

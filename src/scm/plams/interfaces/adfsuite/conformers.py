import math
import os
import re
from os.path import join as opj
from pathlib import Path
from typing import List

import numpy as np

from scm.plams.core.basejob import Job, Results, SingleJob
from scm.plams.core.functions import read_all_molecules_in_xyz_file, requires_optional_package
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob, AMSResults
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.kftools import KFFile
from scm.plams.tools.plot import plot_molecule
from scm.plams.tools.units import Units

__all__ = ["ConformersJob", "ConformersResults"]


class ConformersCollection:
    "in memory conformers collection"

    def __init__(self, path_in: str):
        path = Path(path_in)
        if path.is_dir():
            path = path / "conformers.xyz"
        if not path.name.endswith("conformers.xyz"):
            raise ValueError(f"{path=} not valid")

        assert path.exists()
        molecules = read_all_molecules_in_xyz_file(path)
        for m in molecules:
            m.properties = self.parse_comment(m.properties.pop("comment", ""))
        self.molecules = list(sorted(molecules, key=lambda m: m.properties["Energy in kcal/mol"], reverse=False))

    @property
    def energies(self):
        return list(map(lambda m: m.properties["Energy in kcal/mol"], self.molecules))

    @property
    def relative_energies(self):
        en = self.energies
        min_en = min(en)
        return [e - min_en for e in en]

    def parse_comment(self, c: str):
        pairs = re.findall(r"([^:]+):\s*([^\s]+)", c)
        result = {key.strip(): value for key, value in pairs}
        for k, v in result.items():
            if v.isdigit():
                result[k] = int(v)
            else:
                try:
                    result[k] = float(v)
                except ValueError:
                    pass  # leave as string
        return result

    def __str__(self):
        """
        Print conformer info
        """
        block = "%9s %20s\n" % ("Conformer", "Energy [kcal/mol]")
        min_energy = None
        if len(self.molecules) > 0:
            min_energy = min(self.energies)
        for i in range(len(self.molecules)):
            block += "%9i %20.6f\n" % (
                i,
                self.energies[i] - min_energy,
            )
        return block


class ConformersResults(Results):
    """
    A specialized |Results| subclass for accessing the results of ConformersJob.
    Conformers are sorted by energy, from lowest to highest.
    """

    rkfname = "conformers.rkf"

    def __init__(self, *args, **kwargs):
        Results.__init__(self, *args, **kwargs)
        self._conformers = None
        self.rkf = None

    def rkfpath(self) -> str:
        """Absolute path to the 'conformers.rkf' results file"""
        return opj(self.job.path, self.rkfname)

    def get_lowest_conformer(self) -> Molecule:
        """Return the conformer with the lowest energy"""
        return self._conformers[0]

    def get_conformers(self) -> List[Molecule]:
        """
        Return a list containing all conformers found.
        The conformers are sorted according to their energy, the first element being the lowest energy conformer.
        """
        return self._conformers.molecules

    def get_relative_energies(self, unit="au") -> List[float]:
        """
        Return the relative energies of the conformers i.e. the energy of the conformer minus the energy of the lowest conformer found.
        This list is sorted according to the energy of the conformers, the first element corresponding to the lowest energy conformer.
        So, by definition, the first element will have an energy of 0.
        """
        e = self._conformers.relative_energies
        return Units.convert(e, "kcal/mol", unit)

    def get_energies(self, unit="au") -> List[float]:
        """
        Return the energies of the conformers.
        This list is sorted according to the energy of the conformers, the first element corresponding to the lowest energy conformer.
        """
        e = self._conformers.energies
        return Units.convert(e, "kcal/mol", unit)

    def get_lowest_energy(self, unit="au") -> float:
        """Return the energy of the lowest-energy conformer."""
        return Units.convert(self._conformers.energies[0], "kcal/mol", unit)

    def get_boltzmann_distribution(self, temperature: float) -> List[float]:
        """
        Return the Boltzmann distribution at a given temperature: exp^(E_i/kB*temperature) / (sum_j exp^(E_j/kB*temperature)), where E_i is the energy of conformer i.
        This list is sorted according to the energy of the conformers, the first element corresponding to the lowest energy (and highest probability) conformer.
        The temperature is in Kelvin.
        """
        if temperature <= 0:
            raise ValueError(f"temperature ({temperature}) should be a positive number.")
        kB = 3.166819e-6  # Hartree / Kelvin
        weights = [math.exp(-e / (kB * temperature)) for e in self.get_relative_energies()]
        denominator = sum(weights)
        dist = [w / denominator for w in weights]
        return dist

    def __str__(self):
        return str(self._conformers)

    def get_energy_landscape(self):
        el = AMSResults.EnergyLandscape(None)
        for energy, mol in zip(self.get_energies(), self.get_conformers()):
            state = AMSResults.EnergyLandscape.State(el, None, energy, mol, 1, False)
            el._states.append(state)
        return el

    def collect(self):
        """Collect files present in the job folder.

        Use parent method from |Results| to get a list of all files in the results folder.
        Then instantiate ``self.rkf`` to be a |KFFile| instance for the main ``conformers.rkf`` output file.
        Also instantiate ``self._conformers`` to be a Conformers instance built from it.

        This method is called automatically during the final part of the job execution and there is no need to call it manually.
        """
        Results.collect(self)

        if self.rkfname in self.files:
            rkf_path = opj(self.job.path, self.rkfname)
            self.rkf = KFFile(rkf_path)
            self._conformers = ConformersCollection(self.job.path)

    @requires_optional_package("matplotlib")
    def plot_conformers(self, indices=None, temperature=298, unit="kcal/mol", lowest=True):
        """
        Function for plotting conformers

        indices: None, int or list of int
            If None, will plot at most 3 conformers.

            If int, will plot at most the given number of conformers.

            If list of int (zero-based indices), plot those conformers.

        temperature: float
            Temperature for relative population (printed above the figure).

        unit: str
            Unit for relative energies (printed above the figure)

        lowest: bool
            Only used if ``indices`` is an integer. If True, plot the N lowest energy conformers. If False, plot conformers evenly distributed from the most stable to the least stable.

        """
        import matplotlib.pyplot as plt

        molecules = self.get_conformers()
        energies = self.get_relative_energies(unit)
        populations = self.get_boltzmann_distribution(temperature)

        if isinstance(indices, int):
            N_plot = min(indices, len(energies))
            if lowest:
                indices = list(range(N_plot))
            else:
                indices = np.linspace(0, len(energies) - 1, N_plot, dtype=np.int32)
        if indices is None:
            indices = list(range(min(3, len(energies))))

        fig, axes = plt.subplots(1, len(indices), figsize=(12, 3))
        if len(indices) == 1:
            axes = [axes]

        for ax, i in zip(axes, indices):
            mol = molecules[i]
            E = energies[i]
            population = populations[i]

            plot_molecule(mol, ax=ax)
            ax.set_title(f"#{i+1}\nΔE = {E:.2f} {unit}\nPop.: {population:.3f} (T = {temperature} K)")

        return axes


class ConformersJob(SingleJob):

    _result_type = ConformersResults
    _command = "conformers"

    def __init__(self, name="conformers", molecule=None, **kwargs):
        Job.__init__(self, name=name, **kwargs)

        if molecule is None:
            self.molecule = None
        elif isinstance(molecule, Molecule):
            self.molecule = molecule.copy()
        else:
            raise NotImplementedError("TODO: Implement multiple molecules input.")

    def check(self):
        return True

    def get_input(self):
        sett = self.settings.copy()

        # If there are references to ConformersResults in the settings, expand them into full paths:
        def expand_results_into_paths(s):
            for k, v in s.items():
                if isinstance(v, Settings):
                    expand_results_into_paths(v)
                elif isinstance(v, ConformersResults):
                    s[k] = v.rkfpath()
                elif isinstance(v, list) and any([isinstance(x, ConformersResults) for x in v]):
                    s[k] = [x.rkfpath() if isinstance(x, ConformersResults) else x for x in v]

        expand_results_into_paths(sett)

        return AMSJob(settings=sett, molecule=self.molecule).get_input()

    def get_runscript(self):
        ret = ""
        if "preamble_lines" in self.settings.runscript:
            for line in self.settings.runscript.preamble_lines:
                ret += f"{line}\n"
        ret += 'AMS_JOBNAME="{}" AMS_RESULTSDIR=. $AMSBIN/conformers'.format(self.name)
        if "nproc" in self.settings.runscript:
            ret += " -n {}".format(self.settings.runscript.nproc)
        ret += ' <"{}"'.format(self._filename("inp"))
        if self.settings.runscript.stdout_redirect:
            ret += ' >"{}"'.format(self._filename("out"))
        ret += "\n\n"
        return ret

    @classmethod
    def load_external(cls, path, settings=None, molecule=None, finalize=False, jobname=None):
        """
        Load an external job from *path*.
        """
        if os.path.basename(path)[-4:] == ".rkf":
            path = os.path.dirname(path)
        # molecule = molecule or Molecule.read_userin(path + "/conformers.rkf", program="conformers")
        ret = super().load_external(path, settings, molecule, finalize, jobname)
        return ret

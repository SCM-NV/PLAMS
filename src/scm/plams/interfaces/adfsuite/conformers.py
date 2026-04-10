import math
import os
import re
from os.path import join as opj
from pathlib import Path
from typing import Any, Dict, List, Optional, Union, cast

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
        molecules = read_all_molecules_in_xyz_file(str(path))
        for m in molecules:
            c = m.properties.pop("comment", "")
            assert isinstance(c, str)
            m.properties = Settings(self.parse_comment(c))
        self.molecules = list(sorted(molecules, key=lambda m: m.properties["Energy in kcal/mol"], reverse=False))

    @property
    def energies(self) -> List[float]:
        return list(map(lambda m: m.properties["Energy in kcal/mol"], self.molecules))

    @property
    def relative_energies(self) -> List[float]:
        en = self.energies
        min_en = min(en)
        return [e - min_en for e in en]

    def parse_comment(self, c: str) -> Dict[str, Union[float, str, int]]:
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

    def __str__(self) -> str:
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

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        Results.__init__(self, *args, **kwargs)
        self._conformers: Optional[ConformersCollection] = None
        self.rkf: Optional[KFFile] = None

    def collect(self) -> None:
        """Collect files present in the job folder.

        Use parent method from |Results| to get a list of all files in the results folder.
        Then instantiate ``self.rkf`` to be a |KFFile| instance for the main ``conformers.rkf`` output file.
        Also instantiate ``self._conformers`` to be a Conformers instance built from it.

        This method is called automatically during the final part of the job execution and there is no need to call it manually.
        """
        Results.collect(self)

        if self.rkfname in self.files:
            p = self.get_path()
            self.rkf = KFFile(str(Path(p) / self.rkfname))
            self._conformers = ConformersCollection(path_in=p)

    def get_path(self) -> str:
        p = self.job.path
        if p is None:
            raise ValueError(f"self.job.path is None, maybe you did not run this job? {self.job=}|{self.job.name=}")
        return p

    @property
    def conformers(self) -> ConformersCollection:
        if self._conformers is None:
            self.collect()
        assert isinstance(self._conformers, ConformersCollection)
        return self._conformers

    def rkfpath(self) -> str:
        """Absolute path to the 'conformers.rkf' results file"""
        p = self.get_path()
        return opj(p, self.rkfname)

    def get_lowest_conformer(self) -> Molecule:
        """Return the conformer with the lowest energy"""
        return self.conformers.molecules[0]

    def get_conformers(self) -> List[Molecule]:
        """
        Return a list containing all conformers found.
        The conformers are sorted according to their energy, the first element being the lowest energy conformer.
        """
        return self.conformers.molecules

    def get_relative_energies(self, unit: str = "au") -> List[float]:
        """
        Return the relative energies of the conformers i.e. the energy of the conformer minus the energy of the lowest conformer found.
        This list is sorted according to the energy of the conformers, the first element corresponding to the lowest energy conformer.
        So, by definition, the first element will have an energy of 0.
        """
        e = self.conformers.relative_energies
        return Units.convert(e, "kcal/mol", unit)

    def get_energies(self, unit: str = "au") -> List[float]:
        """
        Return the energies of the conformers.
        This list is sorted according to the energy of the conformers, the first element corresponding to the lowest energy conformer.
        """
        e = self.conformers.energies
        return Units.convert(e, "kcal/mol", unit)

    def get_lowest_energy(self, unit: str = "au") -> float:
        """Return the energy of the lowest-energy conformer."""
        return Units.convert(self.conformers.energies[0], "kcal/mol", unit)

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

    def __str__(self) -> str:
        return str(self._conformers)

    def get_energy_landscape(self) -> AMSResults.EnergyLandscape:
        # Conformers results are not AMSResults, so initialize the landscape container manually.
        el = AMSResults.EnergyLandscape.__new__(AMSResults.EnergyLandscape)
        el._states = []
        el._fragments = []
        el._fstates = []
        for energy, mol in zip(self.get_energies(), self.get_conformers()):
            state = AMSResults.EnergyLandscape.State(el, "", energy, mol, 1, False)
            el._states.append(state)
        return el

    @requires_optional_package("matplotlib")
    def plot_conformers(
        self,
        indices: Optional[Union[int, List[int]]] = None,
        temperature: float = 298.0,
        unit: str = "kcal/mol",
        lowest: bool = True,
    ) -> Union[List[Any], Any]:
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

        returns:
            NDArray[Axes] or Axes depending how many conformers are present/requested to plot
        """
        import matplotlib.pyplot as plt

        molecules = self.get_conformers()
        energies = self.get_relative_energies(unit)
        populations = self.get_boltzmann_distribution(temperature)

        if isinstance(indices, int):
            N_plot = min(indices, len(energies))
            if lowest:
                indices_ok = list(range(N_plot))
            else:
                indices_ok = np.linspace(0, len(energies) - 1, N_plot, dtype=np.int32).tolist()
        else:
            if indices is None:
                indices_ok = list(range(min(3, len(energies))))
            else:
                indices_ok = indices

        fig, axes = plt.subplots(1, len(indices_ok), figsize=(12, 3))
        if len(indices_ok) == 1:
            axes = [axes]

        for ax, i in zip(axes, indices_ok):
            mol = molecules[i]
            E = energies[i]
            population = populations[i]

            plot_molecule(mol, ax=ax)
            ax.set_title(f"#{i+1}\nΔE = {E:.2f} {unit}\nPop.: {population:.3f} (T = {temperature} K)")

        return axes


class ConformersJob(SingleJob):

    _result_type = ConformersResults
    _command = "conformers"

    def __init__(self, name: str = "conformers", molecule: Optional[Molecule] = None, **kwargs: Any) -> None:
        Job.__init__(self, name=name, **kwargs)

        # TODO: include UCS
        # TODO: Implement multiple molecules input.
        if molecule is None:
            self.molecule = None
        elif isinstance(molecule, Molecule):
            self.molecule = molecule.copy()
        else:
            raise ValueError(f"Found: {type(molecule)=} but it must be Molecule or None")

    def check(self) -> bool:
        return True

    def get_input(self) -> str:
        sett = self.settings.copy()

        # If there are references to ConformersResults in the settings, expand them into full paths:
        def expand_results_into_paths(s: Settings) -> None:
            for k, v in s.items():
                if isinstance(v, Settings):
                    expand_results_into_paths(v)
                elif isinstance(v, ConformersResults):
                    s[k] = v.rkfpath()
                elif isinstance(v, list) and any([isinstance(x, ConformersResults) for x in v]):
                    s[k] = [x.rkfpath() if isinstance(x, ConformersResults) else x for x in v]

        expand_results_into_paths(sett)

        return AMSJob(settings=sett, molecule=self.molecule).get_input()

    def get_runscript(self) -> str:
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
    def load_external(
        cls,
        path: Union[str, Path],
        settings: Optional[Settings] = None,
        molecule: Optional[Molecule] = None,
        finalize: bool = False,
        jobname: Optional[str] = None,
    ) -> "ConformersJob":
        """
        Load an external job from *path*.
        """
        p = str(path)
        if os.path.basename(p)[-4:] == ".rkf":
            p = os.path.dirname(p)
        # molecule = molecule or Molecule.read_userin(path + "/conformers.rkf", program="conformers")
        ret = super().load_external(p, settings, molecule, finalize, jobname)
        return cast("ConformersJob", ret)
        #

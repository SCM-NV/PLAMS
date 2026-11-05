from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.recipes._utils import RecipeSettings, Structure, is_chemical_system, to_plams_settings
from scm.plams.tools.units import Units
from scm.plams.core.errors import FileError
from scm.plams.mol.molecule import Molecule
from scm.plams.core.enums import JobStatus

import numpy as np

from os.path import join as opj
from os.path import abspath, isdir, basename
from typing import Dict, List, Tuple, Union, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from scm.base import ChemicalSystem

__all__ = ["ADFFragmentJob", "ADFFragmentResults"]


class ADFFragmentResults(Results):

    def check(self):
        """Check if the calculation finished successfully.

        Overwriting the method of |MultiJob| because it only checks if the job
        finished, not how it finished.
        """
        return all(job.check() for job in self.job.children)

    def get_properties(self):
        """Redirect to |AMSResults| of the full calculation."""
        return self.job.full.results.get_properties()

    def get_main_molecule(self):
        """Redirect to |AMSResults| of the full calculation."""
        return self.job.full.results.get_main_molecule()

    def get_input_molecule(self):
        """Redirect to |AMSResults| of the full calculation."""
        return self.job.full.results.get_input_molecule()

    def get_main_system(self) -> "ChemicalSystem":
        """Redirect to |AMSResults| of the full calculation."""
        return self.job.full.results.get_main_system()

    def get_input_system(self) -> "ChemicalSystem":
        """Redirect to |AMSResults| of the full calculation."""
        return self.job.full.results.get_input_system()

    def get_energy(self, unit="au"):
        """Redirect to |AMSResults| of the full calculation."""
        return self.job.full.results.get_energy(unit)

    def get_dipole_vector(self, unit="au"):
        """Redirect to |AMSResults| of the full calculation."""
        return self.job.full.results.get_dipole_vector(unit)

    def get_energy_decomposition(self, unit="au") -> Dict[str, float]:
        """Get the energy decomposition of the fragment calculation.

        Args:
            unit (str, optional): The unit of the energy. Defaults to 'au'.

        Returns:
            Dict[str, float]: The energy decomposition.
        """
        energy_section = self.job.full.results.read_rkf_section("Energy", file="adf")
        ret = {}
        for k in [
            "Electrostatic Energy",
            "Kinetic Energy",
            "Elstat Interaction",
            "XC Energy",
        ]:
            ret[k] = Units.convert(energy_section[k], "au", unit)

        # most information not available from the KF file
        res = self.job.full.results
        res1 = self.job.f1.results
        res2 = self.job.f2.results

        # remaining information on the output file
        def try_grep_result(key, pattern, match_position=-1, split_position=4):
            match = res.grep_output(pattern)
            if match:
                match = Units.convert(float(match[match_position].split()[split_position]), "au", unit)
                ret[key] = match

        try_grep_result("E_int", "Total Bonding Energy:", match_position=-2, split_position=3)
        try_grep_result("E_int_disp", "Dispersion Energy:", split_position=2)
        try_grep_result("E_Pauli", "Pauli Repulsion (Delta", split_position=4)
        try_grep_result("E_elstat", "Electrostatic Interaction:", split_position=2)
        try_grep_result("E_orb", "Total Orbital Interactions:", split_position=3)

        ret["E_1"] = res1.get_energy(unit=unit)
        ret["E_2"] = res2.get_energy(unit=unit)

        if hasattr(self.job, "f1_opt"):
            ret["E_1_opt"] = self.job.f1_opt.results.get_energy(unit=unit)
            ret["E_prep_1"] = ret["E_1"] - ret["E_1_opt"]
        if hasattr(self.job, "f2_opt"):
            ret["E_2_opt"] = self.job.f2_opt.results.get_energy(unit=unit)
            ret["E_prep_2"] = ret["E_2"] - ret["E_2_opt"]

        if ("E_1_opt" in ret) and ("E_2_opt" in ret):
            ret["E_prep"] = ret["E_prep_1"] + ret["E_prep_2"]
            ret["E_bond"] = ret["E_int"] + ret["E_prep"]

        return ret

    def get_nocv_eigenvalues(self) -> Union[List[float], Tuple[List[float], List[float]]]:
        """Get the NOCV eigenvalues of the full calculation."""
        keys = self.job.full.results.get_rkf_skeleton(file="adf")
        if "NOCV" not in keys:
            raise KeyError("NOCV section not found in the rkf file")
        keys = keys["NOCV"]
        # restricted calculation
        if "NOCV_eigenvalues_restricted" in keys:
            return self.job.full.results.readrkf("NOCV", "NOCV_eigenvalues_restricted", file="adf")
        # unrestricted calculation
        elif "NOCV_eigenvalues_alpha" in keys:
            alpha = self.job.full.results.readrkf("NOCV", "NOCV_eigenvalues_alpha", file="adf")
            beta = self.job.full.results.readrkf("NOCV", "NOCV_eigenvalues_beta", file="adf")
            tpl = (
                [float(x) for x in alpha],
                [float(x) for x in beta],  #  make mypy happy by ensuring return is tuple of list of floats
            )
            return tpl
        else:
            raise KeyError("NOCV eigenvalues not found in the rkf file")

    def check_nocv_eigenvalues(self) -> bool:
        """Basic NOCV check.

        The eigenvalue #i should be equal to -eigenvalue #-i.

        Returns:
            bool: True if the check passes, False otherwise.
        """

        def _check_list(eigenvalues):
            for i in range(int(len(eigenvalues) / 2)):
                if not np.isclose(eigenvalues[i], -1 * eigenvalues[-i - 1], atol=1e-4):
                    return False
            return True

        NOCV_eigenvalues = self.job.results.get_nocv_eigenvalues()
        # spin-restricted calculation
        if isinstance(NOCV_eigenvalues, list):
            return _check_list(NOCV_eigenvalues)
        elif isinstance(NOCV_eigenvalues, tuple):
            for eigenvalues in NOCV_eigenvalues:
                if not _check_list(eigenvalues):
                    return False
            return True
        else:
            raise TypeError("Unexpected type for NOCV eigenvalues")

    def get_nocv_orbital_interaction(self, unit="kcal/mol") -> Union[List[float], Tuple[List[float]]]:
        """Get the NOCV orbital interactions of the full calculation."""

        def _calc_energies(oi: List[float]) -> List[float]:
            ret = []
            for i in range(int(len(oi) / 2)):
                a = oi[i]
                b = oi[len(oi) - i - 1]
                ret.append(a + b)
            return ret

        keys = self.job.full.results.get_rkf_skeleton(file="adf")
        if "NOCV" not in keys:
            raise KeyError("NOCV section not found in the rkf file")
        keys = keys["NOCV"]
        # restricted calculation
        if "NOCV_oi_restricted" in keys:
            oi = self.job.full.results.readrkf("NOCV", "NOCV_oi_restricted", file="adf")
        # unrestricted calculation
        elif "NOCV_oi_alpha" in keys:
            tpl = (
                self.job.full.results.readrkf("NOCV", "NOCV_oi_alpha", file="adf"),
                self.job.full.results.readrkf("NOCV", "NOCV_oi_beta", file="adf"),
            )
        else:
            raise KeyError("NOCV orbital interaction not found in rkf file")
        energies: Union[List[float], Tuple[List[float], List[float]]]  #  make mypy happy
        if isinstance(oi, tuple):
            energies = (_calc_energies(tpl[0]), _calc_energies(tpl[1]))
        elif isinstance(oi, list):
            energies = _calc_energies(oi)
        return Units.convert(energies, "kcal/mol", unit)


class ADFFragmentJob(MultiJob):
    """Subclass of |MultiJob| for ADF fragment calculations."""

    _result_type = ADFFragmentResults

    def __init__(
        self,
        fragment1: Optional[Structure] = None,
        fragment2: Optional[Structure] = None,
        fragment1_opt: Optional[Structure] = None,
        fragment2_opt: Optional[Structure] = None,
        full_settings: Optional[RecipeSettings] = None,
        frag1_settings: Optional[RecipeSettings] = None,
        frag2_settings: Optional[RecipeSettings] = None,
        **kwargs,
    ):
        """Create an ADFFragmentJob with the given fragments.

        The optimized fragment structures can be given as arguments.
        Single point calculations on these structures will then be performed
        to obtain the preparation energy.

        Mixed fragment representations are combined as a ChemicalSystem.
        Fragment-specific and full-system settings override shared settings.

        Subclass of |MultiJob|.

        Args:
            fragment1 (|Molecule| or ChemicalSystem): The first fragment.
            fragment2 (|Molecule| or ChemicalSystem): The second fragment.
            fragment1_opt (|Molecule| or ChemicalSystem, optional): The optimized first fragment.
            fragment2_opt (|Molecule| or ChemicalSystem, optional): The optimized second fragment.
            full_settings (Settings or scm.inputs.AMS): The settings for the full calculation.
            frag1_settings (Settings or scm.inputs.AMS): Settings for the first fragment calculation.
            frag2_settings (Settings or scm.inputs.AMS): Settings for the second fragment calculation.
            **kwargs: Further keyword arguments for |MultiJob|.
        """
        if "settings" in kwargs:
            kwargs["settings"] = to_plams_settings(kwargs["settings"], type(self).__name__)
        MultiJob.__init__(self, **kwargs)
        self.fragment1 = fragment1.copy() if fragment1 is not None else None
        self.frag1_settings = to_plams_settings(frag1_settings, type(self).__name__) or Settings()
        self.fragment2 = fragment2.copy() if fragment2 is not None else None
        if is_chemical_system(self.fragment1) or is_chemical_system(self.fragment2):
            from scm.base import ChemicalSystem

            if isinstance(self.fragment1, Molecule):
                self.fragment1 = ChemicalSystem(AMSJob(molecule=self.fragment1).get_input())
            if isinstance(self.fragment2, Molecule):
                self.fragment2 = ChemicalSystem(AMSJob(molecule=self.fragment2).get_input())
        self.frag2_settings = to_plams_settings(frag2_settings, type(self).__name__) or Settings()
        self.full_settings = to_plams_settings(full_settings, type(self).__name__) or Settings()

        if isinstance(fragment1_opt, Molecule) or is_chemical_system(fragment1_opt):
            self.fragment1_opt = fragment1_opt.copy()
        elif fragment1_opt is not None:
            raise TypeError("fragment1_opt should be a Molecule or ChemicalSystem object, or None")
        if isinstance(fragment2_opt, Molecule) or is_chemical_system(fragment2_opt):
            self.fragment2_opt = fragment2_opt.copy()
        elif fragment2_opt is not None:
            raise TypeError("fragment2_opt should be a Molecule or ChemicalSystem object, or None")

    def prerun(self):  # noqa F811
        """Prepare the fragments and the full calculation."""
        self.f1 = AMSJob(
            name="frag1",
            molecule=self.fragment1,
            settings=self.frag1_settings + self.settings,
        )
        self.f2 = AMSJob(
            name="frag2",
            molecule=self.fragment2,
            settings=self.frag2_settings + self.settings,
        )

        for fragment, name in (
            (self.fragment1, "subsystem1"),
            (self.fragment2, "subsystem2"),
        ):
            if isinstance(fragment, Molecule):
                for atom in fragment:
                    atom.properties.suffix = f"adf.f={name}"
            else:
                fragment.enable_atom_attributes("adf")
                for atom in fragment.atoms:
                    atom.adf.f = name

        self.full = AMSJob(
            name="full",
            molecule=self.fragment1 + self.fragment2,
            settings=self.full_settings + self.settings,
        )

        self.full.settings.input.adf.fragments.subsystem1 = (self.f1, "adf")
        self.full.settings.input.adf.fragments.subsystem2 = (self.f2, "adf")
        self.children = [self.f1, self.f2, self.full]

        if hasattr(self, "fragment1_opt"):
            self.f1_opt = AMSJob(
                name="frag1_opt",
                molecule=self.fragment1_opt,
                settings=self.frag1_settings + self.settings,
            )
            self.children.append(self.f1_opt)
        if hasattr(self, "fragment2_opt"):
            self.f2_opt = AMSJob(
                name="frag2_opt",
                molecule=self.fragment2_opt,
                settings=self.frag2_settings + self.settings,
            )
            self.children.append(self.f2_opt)

    @classmethod
    def load_external(cls, path: str, jobname: Optional[str] = None) -> "ADFFragmentJob":
        """Load the results of the ADFFragmentJob job from an external path.

        Args:
            path (str): The path to the job. It should at least have the
                        subfolders 'frag1', 'frag2' and 'full'.
            jobname (str, optional): The name of the job. Defaults to None.

        Returns:
            ADFFragmentJob: The job with the loaded results.
        """
        if not isdir(path):
            raise FileError("Path {} does not exist, cannot load from it.".format(path))
        path = abspath(path)
        jobname = basename(path) if jobname is None else str(jobname)

        job = cls(name=jobname)
        job.path = path
        job.status = JobStatus.COPIED

        job.f1 = AMSJob.load_external(opj(path, "frag1"))
        job.f2 = AMSJob.load_external(opj(path, "frag2"))
        job.full = AMSJob.load_external(opj(path, "full"))
        job.children = [job.f1, job.f2, job.full]

        if isdir(opj(path, "frag1_opt")):
            job.f1_opt = AMSJob.load_external(opj(path, "frag1_opt"))
            job.children.append(job.f1_opt)
        if isdir(opj(path, "frag2_opt")):
            job.f2_opt = AMSJob.load_external(opj(path, "frag2_opt"))
            job.children.append(job.f2_opt)

        return job

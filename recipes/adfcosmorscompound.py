import os
from collections import OrderedDict
from typing import List, Optional, Union, Dict, Mapping

from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.crs import CRSJob
from scm.plams.tools.kftools import KFFile
from scm.plams.tools.periodic_table import PeriodicTable
from scm.plams.mol.molecule import Molecule
from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.core.functions import add_to_instance
from scm.plams.interfaces.adfsuite.quickjobs import model_to_settings

__all__ = ["ADFCOSMORSCompoundJob", "ADFCOSMORSCompoundResults"]


class ADFCOSMORSCompoundResults(Results):
    """Results class for ADFCOSMORSCompoundJob"""

    def coskfpath(self):
        """
        Returns the path to the resulting .coskf
        """

        return os.path.join(self.job.path, self.job.name + ".coskf")

    def get_main_molecule(self):
        """
        Returns the optimized molecule
        """

        return self.job.children["solv"].results.get_main_molecule()

    def get_input_molecule(self):
        """
        Returns the input molecule
        """
        for job in self.job.children.values():
            return job.results.get_input_molecule()

    def get_sigma_profile(self, subsection: str = "profil"):
        """
        Returns the sigma profile of the molecule. For more details see `CRSResults.get_sigma_profile`.
        """
        return self.job.children["crs"].results.get_sigma_profile(subsection=subsection)


class ADFCOSMORSCompoundJob(MultiJob):
    """
    A class for performing the equivalent of Task: COSMO-RS Compound in the AMS GUI

    Args:
        molecule : a PLAMS |Molecule|

    Keyword Args:
        coskf_name  : A name for the generated .coskf file.  If nothing is specified, the name of the job will be used.
        coskf_dir  : The directory in which to place the generated .coskf file.  If nothing is specified, the file will be put in the plams directory corresponding to the job.
        preoptimization  : If None, do not preoptimize with a fast engine (then initial optimization is done with ADF). Otherwise, can be one of 'UFF', 'GAFF', 'GFNFF', 'GFN1-xTB', 'ANI-2x'. Note that you need valid licenses for ForceField or DFTB or MLPotential to use these preoptimizers.
        singlepoint (bool) :  Run a singlepoint in gasphase and with solvation to generate the .coskf file on the given Molecule. (no geometry optimization). Cannot be combined with ``preoptimization``.
        settings (Settings) : A |Settings| object.  settings.runscript.nproc, settings.input.adf.custom_options. If 'adf' is in settings.input it should be provided without the solvation block.
        name : an optional name for the calculation directory
        mol_info (dict) : an optional dictionary containing information will be written to the Compound Data section within the COSKF file.

    Example:

        .. code-block:: python

            mol = from_smiles('O')
            job = ADFCOSMORSCompoundJob(
                molecule = mol,
                preoptimization = 'UFF',
                coskf_dir = 'coskfs',
                coskf_name = 'Water',
                name = 'H2O',
                mol_info = {'CAS':'7732-18-5'}
            )
            job.run()
            print(job.results.coskfpath())

    """

    _result_type = ADFCOSMORSCompoundResults

    def __init__(
        self,
        molecule: Union[Molecule, None],
        coskf_name: Optional[str] = None,
        coskf_dir: Optional[str] = None,
        preoptimization: Optional[str] = None,
        singlepoint: bool = False,
        settings: Optional[Settings] = None,
        mol_info: Optional[Mapping[str, Union[float, int, str]]] = None,
        **kwargs,
    ):
        """

        Class for running the equivalent of "COSMO-RS Compound" in the AMS
        GUI. Note that these are ADF calculations, not COSMO-RS
        calculations!

        Initialize two or three jobs:

        (optional): Preoptimization with force field or semi-empirical method
        1. Gasphase optimization (BP86, DZP)
        2. Gasphase optimization (BP86, TZP, BeckeGrid Quality Good)
        3. Take optimized structure and run singlepoint with implicit solvation

        Access the result .coskf file with ``job.results.coskfpath()``.
        Note: this file will be called jobname.coskf, where jobname is the
        name of the ADFCOSMORSCompoundJob.

        """
        if preoptimization and singlepoint:
            raise ValueError("Cannot combine preoptimization with singlepoint")

        MultiJob.__init__(self, children=OrderedDict(), **kwargs)
        self.input_molecule = molecule

        self.mol_info = dict()
        if mol_info is not None:
            self.mol_info.update(mol_info)
        # self.mol_info = mol_info
        self.atomic_ion = False  # should be set when molecule is set if using a custom prerun() method
        if molecule is not None:
            self.mol_info["Molar Mass"] = molecule.get_mass()
            self.mol_info["Formula"] = molecule.get_formula()
            self.atomic_ion = len(molecule.atoms) == 1
            try:
                rings = molecule.locate_rings()
                flatten_atoms = [atom for subring in rings for atom in subring]
                nring = len(set(flatten_atoms))
                self.mol_info["Nring"] = int(nring)
            except:
                pass

        self.settings = settings or Settings()

        self.coskf_name = coskf_name
        self.coskf_dir = coskf_dir

        if self.coskf_dir is not None and not os.path.exists(self.coskf_dir):
            os.mkdir(self.coskf_dir)

        if self.coskf_name is not None and isinstance(self.coskf_name, str) and not self.coskf_name.endswith(".coskf"):
            self.coskf_name += ".coskf"

        gas_s = Settings()
        gas_s += self.adf_settings(solvation=False, settings=self.settings)
        gas_job = AMSJob(settings=gas_s, name="gas")

        if not singlepoint:
            if preoptimization:
                preoptimization_s = Settings()
                preoptimization_s.runscript.nproc = 1
                preoptimization_s.input.ams.Task = "GeometryOptimization"
                preoptimization_s += model_to_settings(preoptimization)
                preoptimization_job = AMSJob(
                    settings=preoptimization_s, name="preoptimization", molecule=self.input_molecule
                )
                self.children["preoptimization"] = preoptimization_job

            gas_s = Settings()
            gas_s.input.ams.Task = "GeometryOptimization"
            gas_s += self.adf_settings(solvation=False, settings=self.settings)
            gas_job = AMSJob(settings=gas_s, name="gas")

            if preoptimization:

                @add_to_instance(gas_job)
                def prerun(self):  # noqa: F811
                    self.molecule = self.parent.children["preoptimization"].results.get_main_molecule()

            else:
                gas_job.molecule = self.input_molecule

            self.children["gas"] = gas_job

        solv_s = Settings()
        solv_s.input.ams.Task = "SinglePoint"
        solv_job = AMSJob(settings=solv_s, name="solv")

        if singlepoint:

            gas_job.settings.input.ams.Task = "SinglePoint"

            @add_to_instance(gas_job)
            def prerun(self):  # noqa: F811
                self.molecule = self.parent.input_molecule

            self.children["gas"] = gas_job

        @add_to_instance(solv_job)
        def prerun(self):  # noqa: F811
            gas_job.results.wait()
            self.settings.input.ams.EngineRestart = "../gas/adf.rkf"
            self.settings.input.ams.LoadSystem.File = "../gas/ams.rkf"
            molecule_charge = gas_job.results.get_main_molecule().properties.get("charge", 0)
            self.settings.input.ams.LoadSystem._1 = f"# {self.parent.name}"
            self.settings.input.ams.LoadSystem._2 = f"# charge {molecule_charge}"
            self.settings += self.parent.adf_settings(
                solvation=True,
                settings=self.parent.settings,
                elements=list(set(at.symbol for at in self.parent.input_molecule)),
                atomic_ion=self.parent.atomic_ion,
            )
            # self.settings.input.ams.EngineRestart = self.parent.children['gas'].results.rkfpath(file='adf') # this doesn't work with PLAMS restart since the file will refer to the .res directory (so the job is rerun needlessly)
            # self.settings.input.ams.LoadSystem.File = self.parent.children['gas'].results.rkfpath(file='ams')
            # cannot copy to gasphase-ams.rkf etc. because that conflicts with PLAMS restarts
            # shutil.copyfile(gas_job.results.rkfpath(file='ams'), os.path.join(self.path, 'gasphase-ams.rkf'))
            # shutil.copyfile(gas_job.results.rkfpath(file='adf'), os.path.join(self.path, 'gasphase-adf.rkf'))

        @add_to_instance(solv_job)
        def postrun(self):
            self.parent.convert_to_coskf(
                self.results.rkfpath(file="adf"),
                os.path.join(
                    self.parent.coskf_dir if self.parent.coskf_dir is not None else self.parent.path,
                    self.parent.coskf_name if self.parent.coskf_name is not None else self.parent.name + ".coskf",
                ),
                self.parent.mol_info,
            )

        self.children["solv"] = solv_job

        sigma_s = Settings()
        sigma_s.input.property._h = "PURESIGMAPROFILE"

        compounds = [Settings()]
        sigma_s.input.compound = compounds
        crsjob = CRSJob(settings=sigma_s, name="sigma")

        @add_to_instance(crsjob)
        def prerun(self):  # noqa F811
            self.parent.children["solv"].results.wait()
            self.settings.input.compound[0]._h = os.path.join(
                self.parent.path if self.parent.coskf_dir is None else os.path.join(os.getcwd(), self.parent.coskf_dir),
                self.parent.coskf_name if self.parent.coskf_name is not None else self.parent.name + ".coskf",
            )

        self.children["crs"] = crsjob

    @staticmethod
    def _get_radii() -> Dict[str, float]:
        """Method to get the atomic radii from solvent.txt (for some elements the radii are instead the Klamt radii)"""
        with open(os.path.expandvars("$AMSHOME/data/gui/solvent.txt"), "r") as f:
            mod_allinger_radii = [float(x) for i, x in enumerate(f) if i > 0]
        radii = {PeriodicTable.get_symbol(i): r for i, r in enumerate(mod_allinger_radii, 1) if i <= 118}
        klamt_radii = {
            "H": 1.30,
            "C": 2.00,
            "N": 1.83,
            "O": 1.72,
            "F": 1.72,
            "Si": 2.48,
            "P": 2.13,
            "S": 2.16,
            "Cl": 2.05,
            "Br": 2.16,
            "I": 2.32,
        }
        radii.update(klamt_radii)

        return radii

    @staticmethod
    def solvation_settings(elements: Optional[List[str]] = None, atomic_ion: bool = False) -> Settings:
        sett = Settings()

        radii = {
            "H": 1.3,
            "He": 1.275,
            "Li": 2.125,
            "Be": 1.858,
            "B": 1.792,
            "C": 2.0,
            "N": 1.83,
            "O": 1.72,
            "F": 1.72,
            "Ne": 1.333,
            "Na": 2.25,
            "Mg": 2.025,
            "Al": 1.967,
            "Si": 2.48,
            "P": 2.13,
            "S": 2.16,
            "Cl": 2.05,
            "Ar": 1.658,
            "K": 2.575,
            "Ca": 2.342,
            "Sc": 2.175,
            "Ti": 1.992,
            "V": 1.908,
            "Cr": 1.875,
            "Mn": 1.867,
            "Fe": 1.858,
            "Co": 1.858,
            "Ni": 1.85,
            "Cu": 1.883,
            "Zn": 1.908,
            "Ga": 2.05,
            "Ge": 2.033,
            "As": 1.967,
            "Se": 1.908,
            "Br": 2.16,
            "Kr": 1.792,
            "Rb": 2.708,
            "Sr": 2.5,
            "Y": 2.258,
            "Zr": 2.117,
            "Nb": 2.025,
            "Mo": 1.992,
            "Tc": 1.967,
            "Ru": 1.95,
            "Rh": 1.95,
            "Pd": 1.975,
            "Ag": 2.025,
            "Cd": 2.083,
            "In": 2.2,
            "Sn": 2.158,
            "Sb": 2.1,
            "Te": 2.033,
            "I": 2.32,
            "Xe": 1.9,
            "Cs": 2.867,
            "Ba": 2.558,
            "La": 2.317,
            "Ce": 2.283,
            "Pr": 2.275,
            "Nd": 2.275,
            "Pm": 2.267,
            "Sm": 2.258,
            "Eu": 2.45,
            "Gd": 2.258,
            "Tb": 2.25,
            "Dy": 2.242,
            "Ho": 2.225,
            "Er": 2.225,
            "Tm": 2.225,
            "Yb": 2.325,
            "Lu": 2.208,
            "Hf": 2.108,
            "Ta": 2.025,
            "W": 1.992,
            "Re": 1.975,
            "Os": 1.958,
            "Ir": 1.967,
            "Pt": 1.992,
            "Au": 2.025,
            "Hg": 2.108,
            "Tl": 2.158,
            "Pb": 2.283,
            "Bi": 2.217,
            "Po": 2.158,
            "At": 2.092,
            "Rn": 2.025,
            "Fr": 3.033,
            "Ra": 2.725,
            "Ac": 2.567,
            "Th": 2.283,
            "Pa": 2.2,
            "U": 2.1,
            "Np": 2.1,
            "Pu": 2.1,
            "Am": 2.1,
            "Cm": 2.1,
            "Bk": 2.1,
            "Cf": 2.1,
            "Es": 2.1,
            "Fm": 2.1,
            "Md": 2.1,
            "No": 2.1,
            "Lr": 2.1,
            "Rf": 2.1,
            "Db": 2.1,
            "Sg": 2.1,
            "Bh": 2.1,
            "Hs": 2.1,
            "Mt": 2.1,
            "Ds": 2.1,
            "Rg": 2.1,
            "Cn": 2.1,
            "Nh": 2.1,
            "Fl": 2.1,
            "Mc": 2.1,
            "Lv": 2.1,
            "Ts": 2.1,
            "Og": 2.1,
        }  # from _get_radii()

        if elements:
            radii = {k: radii[k] for k in sorted(elements)}

        if atomic_ion:
            charge_method = "method=atom corr"
        else:
            charge_method = "method=Conj corr"

        sett.input.adf.solvation = {
            "surf": "Delley",
            "solv": "name=CRS cav0=0.0 cav1=0.0",
            "charged": charge_method,
            "c-mat": "Exact",
            "scf": "Var All",
            "radii": radii,
        }
        return sett

    @staticmethod
    def adf_settings(
        solvation: bool, settings=None, elements: Optional[List[str]] = None, atomic_ion: bool = False
    ) -> Settings:
        """
        Returns ADF settings with or without solvation

        If solvation == True, then also include the solvation block.
        """

        s = Settings()
        if settings:
            s = settings.copy()
        if "basis" not in s.input.adf and "xc" not in s.input.adf:
            s.input.adf.Basis.Type = "TZP"
            s.input.adf.Basis.Core = "Small"
            s.input.adf.XC.GGA = "BP86"
            s.input.adf.Symmetry = "NOSYM"
            s.input.adf.BeckeGrid.Quality = "Good"
        if solvation:
            s += ADFCOSMORSCompoundJob.solvation_settings(elements=elements, atomic_ion=atomic_ion)
        return s

    @staticmethod
    def convert_to_coskf(rkf: str, coskf: str, mol_info: Optional[Mapping[str, Union[float, int, str]]] = None):
        """
        rkf: str
            absolute path to adf.rkf

        coskf: str
            path to write out the resulting .coskf file

        mol_info: dict[str, float | int | str], optional
            Optional information to write out in the "Compound Data" section of the .coskf file

        """
        f = KFFile(rkf)
        cosmo = f.read_section("COSMO")
        coskf_file = KFFile(coskf, autosave=False)
        for k, v in cosmo.items():
            coskf_file.write("COSMO", k, v)
        if mol_info is not None:
            for key, value in mol_info.items():
                # print(f"write to coskf {key}: {value}")
                coskf_file.write("Compound Data", key, value)
        coskf_file.save()

from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.tools.kftools import KFFile
from scm.plams.mol.molecule import Molecule
from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.conformers import ConformersJob
from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob
import os


class ADFCOSMORSConfFilter:
    """
    This class allows the user to specify criteria on which to filter sets of conformers.

    Args:
        max_num_confs : The maximum number of conformers to bring forward to the next step.  These are automatically ordered by energy, so the first *n* lowest energy structures are carried carried forward to the next step.
        max_energy_range : The maximum allowable difference (in kcal/mol) between the lowest energy conformer and the highest energy conformer that will be carried forward.  For example, for a max_energy_range of 2 kcal/mol and a lowest energy structure of 1 kcal/mol, any structure with an energy above 3 (1+2) kcal/mol will be filtered out.
    """

    def __init__(self, max_num_confs=None, max_energy_range=None):
        self.max_num_confs = max_num_confs
        self.max_energy_range = max_energy_range


class ADFCOSMORSConfResults(Results):
    pass


class ADFCOSMORSConfJob(MultiJob):
    """
    This class allows for the user to implement a custom workflow for generating multiple conformers for use with COSMO-RS.  The class allows the user to input conformer generation strategies, job types to refine the energies/geometries, and filters to remove conformers between each calculation step.

    Args:
        molecule (|Molecule|) : A plams Molecule instance

    Keyword Args:
        conf_gen (ConformersJob) : A ConformersJob instance
        first_filter (ADFCOSMORSConfFilter) : The ADFCOSMORSConfFilter to apply to the initial set of sampled conformers
        additional : A list of (|Settings|, ADFCOSMORSConfFilter) tuples.  The elements of this list represent additional calculation steps (e.g., progressively higher levels of theory) and possibly different filters to apply at each step.
        final_filter (ADFCOSMORSConfFilter) : The ADFCOSMORSConfFilter to apply to the results of the gas phase ADF calculation for the conformers
        adf_singlepoint (bool) : A boolean indicating if the adf gas phase calculation in the cosmo-rs compound task should be a single point.  This defaults to False.
        initial_conformers : the (integer) number of initially sampled conformers.  This is only applied if the default conformer generation strategy is used.
        coskf_dir : a directory to put all the .coskf files generated for the conformers.  If this keyword is not specified, the .coskf files are put in the directory containing the adf gas phase calculation results
        coskf_name : a base name to be used with conformers.  All conformers will have the name *coskf_name*_i where i is the index of the unique conformer.  If not specified, the base name becomes simply *conformer*
        name : an optional name for the calculation directory.
        mol_info (dict) : an optional dictionary containing information will be written to the Compound Data section within the COSKF file.

    """

    _result_type = ADFCOSMORSConfResults

    def __init__(
        self,
        molecule: Molecule,
        conf_gen=None,
        first_filter=None,
        additional=None,
        final_filter=None,
        adf_singlepoint=False,
        initial_conformers=500,
        coskf_dir=None,
        coskf_name=None,
        mol_info={},
        **kwargs,
    ):

        super().__init__(children={}, **kwargs)

        self.job_count = 0

        self.mol = molecule
        mol_info["Molar Mass"] = molecule.get_mass()
        mol_info["Formula"] = molecule.get_formula()
        try:
            rings = molecule.locate_rings()
            flatten_atoms = [atom for subring in rings for atom in subring]
            nring = len(set(flatten_atoms))
            mol_info["Nring"] = int(nring)
        except:
            pass
        self.mol_info = mol_info

        self.adf_results = False
        self.cosmo_results = False

        self.adf_singlepoint = adf_singlepoint

        self.coskf_dir = coskf_dir
        self.coskf_name = coskf_name

        if self.coskf_dir is not None and not os.path.exists(self.coskf_dir):
            os.mkdir(self.coskf_dir)

        self.initial_conformers = initial_conformers

        if conf_gen is None:
            self.conf_gen = self.default_confgen()
        else:
            if not isinstance(conf_gen, ConformersJob):
                print(
                    "Wrong type for argument conf_gen.  Expected ConformersJob instance.  Using the default conformer generator."
                )
                self.conf_gen = self.default_confgen()
            else:
                self.conf_gen = conf_gen

        self.job_settings = [self.conf_gen.settings]
        self.filters = [None, first_filter]

        if additional is not None:
            for js, f in additional:
                self.job_settings.append(js)
                self.filters.append(f)

        self.final_filter = final_filter
        self.filters.append(final_filter)

        if not self.has_valid_filter_settings():
            pass

        self.children["job_0"] = self.conf_gen

    def default_confgen(self):

        sett = Settings()
        sett.input.AMS.Generator.RDKit
        sett.input.AMS.Generator.RDKit.InitialNConformers = self.initial_conformers
        return ConformersJob(name="conformers_uff", molecule=self.mol, settings=sett)

    def make_intermediate_job(self, settings):

        settings.input.AMS.InputConformersSet = self.children[f"job_{self.job_count-1}"].results
        self._add_filter(settings)
        return ConformersJob(name=f"additional_{self.job_count}", settings=settings)

    def make_adf_job(self):

        sett = ADFCOSMORSCompoundJob.adf_settings(False)
        if not self.adf_singlepoint:
            sett.input.AMS.Task = "Optimize"
            sett.input.AMS.GeometryOptimization.UseAMSWorker = "False"
            # sett.input.ams.GeometryOptimization.ConvergenceQuality = 'Good'
        else:
            sett.input.AMS.Task = "Score"
        sett.input.AMS.InputConformersSet = self.children[f"job_{self.job_count-1}"].results
        self._add_filter(sett)

        return ConformersJob(name="adf_conformers", settings=sett)

    def make_filter_job(self):
        sett = Settings()
        sett.input.AMS.Task = "Filter"
        sett.input.AMS.InputConformersSet = self.children["adf_job"].results
        self._add_filter(sett)

        return ConformersJob(name="adf_filter", settings=sett)

    def make_cosmo_job(self):

        sett = ADFCOSMORSCompoundJob.adf_settings(True, elements=list(set(at.symbol for at in self.mol)))
        sett.input.AMS.Task = "Replay"

        if self.final_filter is not None:
            previous_job = "filter_job"
        else:
            previous_job = "adf_job"

        self.children[previous_job].results.wait()
        sett.input.AMS.Replay.File = self.children[previous_job].results["conformers.rkf"]
        sett.input.AMS.Replay.StoreAllResultFiles = "True"

        return AMSJob(name="replay", settings=sett)

    def new_children(self):
        """Don't doc this

        :meta private:
        """
        self.job_count += 1
        if self.job_count < len(self.job_settings):
            settings = self.job_settings[self.job_count]
            new_job = self.make_intermediate_job(settings)
            return {f"job_{self.job_count}": new_job}

        if "adf_job" not in self.children:
            return {"adf_job": self.make_adf_job()}

        if "filter_job" not in self.children and self.final_filter is not None:
            return {"filter_job": self.make_filter_job()}

        if "cosmo_job" not in self.children:
            return {"cosmo_job": self.make_cosmo_job()}

        return None

    def postrun(self):
        self._make_coskfs()

    def _make_coskfs(self):

        base_name = self.coskf_name if self.coskf_name is not None else "conformer"
        if self.final_filter is not None:
            previous_job = "filter_job"
        else:
            previous_job = "adf_job"

        if self.coskf_dir is None:
            self.coskf_dir = self.children[previous_job].path
        for i, E in enumerate(self.children[previous_job].results.get_energies("Ha")):
            if f"Frame{i+1}.rkf" in self.children["cosmo_job"].results:
                cosmo_section = self.children["cosmo_job"].results.read_rkf_section("COSMO", f"Frame{i+1}")
                cosmo_section["Gas Phase Bond Energy"] = E
                name = f"{base_name}_{i}.coskf"

                fullname = os.path.join(self.coskf_dir, name)
                if os.path.exists(fullname):
                    os.remove(fullname)

                coskf = KFFile(os.path.join(self.coskf_dir, name), autosave=False)
                for key, val in cosmo_section.items():
                    coskf.write("COSMO", key, val)

                for key, value in self.mol_info.items():
                    # print(f"write to coskf {key}: {value}")
                    coskf.write("Compound Data", key, value)
                coskf.save()

    def _add_filter(self, sett):

        filt = self.filters[self.job_count]
        if filt is not None:
            if filt.max_num_confs is not None:
                sett.input.AMS.InputMaxConfs = filt.max_num_confs
            if filt.max_energy_range is not None:
                sett.input.AMS.InputMaxEnergy = filt.max_energy_range

    def has_valid_filter_settings(self):

        for js, f in zip(self.job_settings, self.filters):
            if js is not None and not isinstance(js, Settings):
                return False
            if f is not None and not isinstance(f, ADFCOSMORSConfFilter):
                return False
        return True

from scm.plams.core.functions import log
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.core.settings import Settings
from scm.plams.recipes.adffragment import ADFFragmentJob, ADFFragmentResults
from scm.plams.tools.units import Units
from scm.plams.core.errors import FileError
from scm.plams.core.functions import add_to_instance


from os.path import join as opj
from os.path import abspath, isdir, basename, relpath
from os import symlink
from typing import Dict, List, Union

__all__ = ["BANDFragmentJob", "BANDFragmentResults"]


class BANDFragmentResults(ADFFragmentResults):
    """Subclass of |ADFFragmentResults| for BAND calculations."""

    def get_energy_decomposition(self, unit="kJ/mol") -> Dict[str, float]:
        """Get the energy decomposition of the fragment calculation.

        Args:
            unit (str, optional): The unit of the energy. Defaults to 'kJ/mol'.

        Returns:
            Dict[str, float]: The energy decomposition.
        """
        res = self.job.full.results
        res1 = self.job.f1.results
        res2 = self.job.f2.results
        ret = {}
        pos = 2
        # E_int appears in a comment below the PEDA Table of the output
        ret["E_int"] = Units.convert(float(res.grep_output("E_int")[-2].split()[pos]), "au", unit)
        ret["E_int_disp"] = Units.convert(float(res.grep_output("E_disp")[-1].split()[pos]), "au", unit)
        ret["E_Pauli"] = Units.convert(float(res.grep_output("E_Pauli")[-1].split()[pos]), "au", unit)
        ret["E_elstat"] = Units.convert(float(res.grep_output("E_elstat")[-1].split()[pos]), "au", unit)
        ret["E_orb"] = Units.convert(float(res.grep_output("E_orb")[-1].split()[pos]), "au", unit)

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


class BANDFragmentJob(ADFFragmentJob):
    _result_type = BANDFragmentResults

    def create_mapping_setting(self) -> None:
        if "fragment" in self.full_settings.input.band:
            log(
                "Fragment already present in full_settings. Assuming that the user has already set up the mapping. Skipping the mapping setup.",
                level=1,
            )
            return
        # first fragment 1 then fragment 2
        set1 = Settings()
        set1.atommapping = {str(i + 1): str(i + 1) for i in range(len(self.fragment1))}
        set2 = Settings()
        set2.atommapping = {str(i + 1): str(i + 1 + len(self.fragment1)) for i in range(len(self.fragment2))}
        # get the correct restart files
        # working with relative paths does not work for unknown reasons
        # using symlinks instead
        set1.filename = f"{self.f1.name}.rkf"
        set2.filename = f"{self.f2.name}.rkf"
        self.full_settings.input.band.fragment = [set1, set2]

    def new_children(self) -> Union[None, List[AMSJob]]:
        """After the first round, add the full job to the children list."""
        if hasattr(self, "full"):
            return None
        else:
            # create the correct mapping settings for the full job
            self.create_mapping_setting()
            # create the full job
            self.full = AMSJob(
                name="full", molecule=self.fragment1 + self.fragment2, settings=self.settings + self.full_settings
            )
            # dependencies are optional, but let's set them up
            self.full.depend += [self.f1, self.f2]
            # save the fragment paths for the prerun of the full job
            self.full.frag_paths = []
            for job in [self.f1, self.f2]:
                self.full.frag_paths.append(job.path)

            # edit full prerun to create symlinks
            @add_to_instance(self.full)
            def prerun(self):
                """Create symlinks for the restart files."""
                for i, job in enumerate(["frag1", "frag2"]):
                    rel_path = relpath(self.frag_paths[i], self.path)
                    symlink(opj(rel_path, "band.rkf"), opj(self.path, f"{job}.rkf"))

            return [self.full]

    def prerun(self) -> None:
        """Creates the fragment jobs."""
        self.f1 = AMSJob(name="frag1", molecule=self.fragment1, settings=self.settings)
        self.f2 = AMSJob(name="frag2", molecule=self.fragment2, settings=self.settings)
        self.children += [self.f1, self.f2]

    @classmethod
    def load_external(cls, path: str, jobname: str = None) -> "BANDFragmentJob":
        """Load the results of the BANDFragmentJob job from an external path.

        Args:
            path (str): The path to the job. It should at least have the
                        subfolders 'frag1', 'frag2' and 'full'.
            jobname (str, optional): The name of the job. Defaults to None.

        Returns:
            BANDFragmentJob: The job with the loaded results.
        """
        if not isdir(path):
            raise FileError("Path {} does not exist, cannot load from it.".format(path))
        path = abspath(path)
        jobname = basename(path) if jobname is None else str(jobname)

        job = cls(name=jobname)
        job.path = path
        job.status = "copied"

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


class NOCVBandFragmentJob(BANDFragmentJob):
    _result_type = BANDFragmentResults

    def __init__(self, nocv_settings=None, **kwargs) -> None:
        """Create a BANDFragmentJob with final NOCV calculation.

        Args:
            nocv_settings (Settings, optional): Settings for the NOCV calculation. Defaults to None.

        """
        super().__init__(**kwargs)
        self.nocv_settings = nocv_settings or Settings()

    def new_children(self) -> Union[None, List[AMSJob]]:
        """After the first round, add the full job to the children list.
        After the second round, add the NOCV job to the children list."""
        # new_children of BANDFragmentJob creates the full job
        ret = super().new_children()
        if ret is not None:
            return ret
        if hasattr(self, "nocv"):
            return None
        else:
            # add NOCV run
            # settings for the NOCV calculation
            set = self.settings + self.full_settings + self.nocv_settings
            # set the restart file
            set.input.band.restart.file = "full.rkf"
            # copy the fragment settings
            set.input.band.fragment = [fset.copy() for fset in self.full.settings.input.band.fragment]

            self.nocv = AMSJob(name="NOCV", molecule=self.fragment1 + self.fragment2, settings=set)

            self.nocv.frag_paths = []
            for job in [self.f1, self.f2, self.full]:
                self.nocv.frag_paths.append(job.path)

            # edit NOCV prerun to create symlinks
            @add_to_instance(self.nocv)
            def prerun(self):
                """Create symlinks for the restart files."""
                for i, job in enumerate(["frag1", "frag2", "full"]):
                    rel_path = relpath(self.frag_paths[i], self.path)
                    symlink(opj(rel_path, "band.rkf"), opj(self.path, f"{job}.rkf"))

            return [self.nocv]

    @classmethod
    def load_external(cls, path: str, jobname: str = None) -> "NOCVBandFragmentJob":
        """Load the results of the BANDFragmentJob job from an external path.

        Args:
            path (str): The path to the job. It should at least have the
                        subfolders 'frag1', 'frag2', 'full' and 'NOCV'.
            jobname (str, optional): The name of the job. Defaults to None.

        Returns:
            NOCVBandFragmentJob: The job with the loaded results.
        """
        job = super().load_external(path, jobname)
        job.nocv = AMSJob.load_external(opj(path, "NOCV"))
        job.children.append(job.nocv)
        return job

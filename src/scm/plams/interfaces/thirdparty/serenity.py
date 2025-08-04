from scm.plams.core.basejob import SingleJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.core.functions import requires_optional_package
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.units import Units

from pathlib import Path
import numpy as np
import h5py
import shutil
import os

__all__ = ["SerenityJob", "SerenityResults", "SerenitySettings"]


class SerenitySettings(Settings):
    """A Subclass of the Settings class to set the iteration order through keys."""

    def __iter__(self):
        """Iteration through keys follows insertion order instead of lexicographical order."""
        return iter(self.keys())


class SerenityResults(Results):
    """A Subclass of the Results class to access the results given in the Serenity output. Can be extended as required."""

    def _get_energy_type(self, search="Total Energy", index=-1, unit="a.u."):
        s = self.grep_output(search)
        s = [item for item in s if not "..." in item]
        s = s[index]

        def convert_with_fallback(item):
            try:
                return float(item.split()[-1])
            except ValueError:
                try:
                    return float(item.split()[-2])
                except ValueError:
                    raise ValueError("Can not convert type to float.")

        if not isinstance(index, slice):
            return Units.convert(convert_with_fallback(s), "a.u.", unit)
        else:
            return [Units.convert(convert_with_fallback(x), "a.u.", unit) for x in s]

    def get_energy(self, index=-1, unit="a.u."):
        """Returns 'Total Energy (DFT)/(HF):' from the output file. Set ``index`` to choose the n-th occurence of the total energy in the output. Defaults to the last occurence."""
        return self._get_energy_type("Total Energy", index=index, unit=unit)

    def get_ccsd_energy_correction(self, index=-1, unit="a.u."):
        """Returns the 'CCSD Energy Correction' from the output file. Set ``index`` to choose the n-th occurrence of the CCSD energy correction in the output. Defaults to the last occurrence."""
        return self._get_energy_type("CCSD Energy Correction", index=index, unit=unit)

    def get_triples_energy_correction(self, index=-1, unit="a.u."):
        return self._get_energy_type("Triples Energy Correction", index=index, unit=unit)

    def get_local_ccsd_energy(self, index=-1, unit="a.u."):
        return self._get_energy_type("Total Local-CCSD Energy", index=index, unit=unit)


class SerenityJob(SingleJob):
    """A Subclass of the SingleJob class representing a computational job with Serenity."""

    _result_type = SerenityResults

    def _get_ready(self):
        """If molecule is defined, xyz files are generated accordingly. If there is data to transfer between ADF and Serenity, file formats are converted."""
        super()._get_ready()
        if isinstance(self.molecule, dict):
            for name, mol in self.molecule.items():
                mol.write(os.path.join(self.path, f"{name}.xyz"))
        elif isinstance(self.molecule, Molecule):
            self.molecule.write(os.path.join(self.path, "mol.xyz"))

        base_path = Path(self.path).parent
        ams_job_dirs = [
            d for d in base_path.iterdir() if d.is_dir() and any(f.startswith("SerenityAMS") for f in os.listdir(d))
        ]

        if ams_job_dirs:
            self._translate_AMSJob_to_Serenity(ams_job_dirs[0])  # Processes the first matching directory
        return

    @requires_optional_package("h5py")
    def _translate_AMSJob_to_Serenity(self, AMSJob_folder):
        def convert_to_HDF5(binary_file, hdf5_file, dataset_name, shape=None):
            hdf5_file_path = Path(self.path) / hdf5_file
            data = np.fromfile(str(binary_file), dtype=np.float64)
            if shape:
                data = data.reshape(shape)
            with h5py.File(str(hdf5_file_path), "w") as h5file:
                h5file.create_dataset(dataset_name, data=data, dtype="float64")

        n_basis_func_file_source = AMSJob_folder / "SerenityAMS.nBasisFunc.txt"
        n_basis_func_file_dest = Path(self.path) / "SerenityAMS.nBasisFunc.txt"
        shutil.copy(n_basis_func_file_source, n_basis_func_file_dest)
        with open(n_basis_func_file_dest, "r") as file:
            dim_size = int(file.readline().strip())
        convert_to_HDF5(
            AMSJob_folder / "SerenityAMS.ERIs.bin",
            "SerenityAMS.ERIs.h5",
            "ERIs",
            (dim_size, dim_size, dim_size, dim_size),
        )
        convert_to_HDF5(AMSJob_folder / "SerenityAMS.orbEnergies.bin", "SerenityAMS.orbEnergies.h5", "orbitalEnergies")

    def get_input(self):
        """Transforms all contents of the ``input`` branch of |Settings| into a string with blocks, subblocks, nested blocks, keys and values. Reserved keywords are handled in a specific manner."""
        currentworkdir = Path.cwd()
        _reserved_keywords = ["task", "system"]

        def parse(key, value, indent=""):
            ret = ""
            if isinstance(value, Settings):
                if not any(k == key for k in _reserved_keywords):
                    ret += "{}+{}\n".format(indent, key)
                    for el in value:
                        ret += parse(el, value[el], indent + "  ")
                    ret += "{}-{}\n".format(indent, key)

                elif "task" in key:
                    for el in value:
                        ret += "{}+{}  {}\n".format(indent, key, el)
                        for v in value[el]:
                            ret += parse(v, value[el][v], indent + "  ")
                        ret += "{}-{}\n".format(indent, key)

                elif "system" in key:
                    molecule_exists = hasattr(self, "molecule") and self.molecule is not None
                    if molecule_exists:
                        if isinstance(self.molecule, dict):
                            for el in value:
                                ret += "{}+{}\n".format(indent, key)
                                ret += "{}{}  name  {}\n".format(indent, indent, el)
                                ret += "{}{}  geometry  {}.xyz\n".format(indent, indent, el)
                                for v in value[el]:
                                    ret += parse(v, value[el][v], indent + "  ")
                                ret += "{}-{}\n".format(indent, key)
                        elif isinstance(self.molecule, Molecule):
                            for el in value:
                                ret += "{}+{}\n".format(indent, key)
                                ret += "{}{}  name  {}\n".format(indent, indent, el)
                                ret += "{}{}  geometry  mol.xyz\n".format(indent, indent)
                                for v in value[el]:
                                    ret += parse(v, value[el][v], indent + "  ")
                                ret += "{}-{}\n".format(indent, key)
                    else:
                        for el in value:
                            ret += "{}+{}\n".format(indent, key)
                            ret += "{}{}  name  {}\n".format(indent, indent, el)
                            for v in value[el]:
                                ret += parse(v, value[el][v], indent + "  ")
                            ret += "{}-{}\n".format(indent, key)

            elif isinstance(value, list):
                for el in value:
                    ret += parse(key, el, indent)

            elif value == "" or value is True:
                ret += "{}{}\n".format(indent, key)
            else:
                if key == "geometry":
                    geometry_path = currentworkdir / value
                    abs_geometry_path = geometry_path.resolve()
                    ret += "{}{}  {}\n".format(indent, key, abs_geometry_path)
                else:
                    ret += "{}{}  {}\n".format(indent, key, str(value))
            return ret

        inp = ""

        for item in self.settings.input:
            inp += parse(item, self.settings.input[item]) + "\n"

        return inp

    def get_runscript(self):
        """Returned runscript: ``serenity myinput.in |tee ser_myinput.out`` or ``serenity myinput.in``"""
        input_file = self._filename("inp")
        # output_file = self._filename("out")

        # return 'serenity {} |tee ser_{}'.format(input_file, output_file)
        return "serenity {}".format(input_file)

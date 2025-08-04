import inspect
import os
import subprocess
from itertools import cycle
from typing import Optional, List, Dict

import numpy as np

from scm.plams.interfaces.adfsuite.scmjob import SCMJob, SCMResults
from scm.plams.tools.units import Units
from scm.plams.core.functions import log

__all__ = ["CRSResults", "CRSJob"]


class CRSResults(SCMResults):
    """A |SCMResults| subclass for accessing results of |CRSJob|."""

    _kfext = ".crskf"
    _rename_map = {"CRSKF": "$JN.crskf"}

    @property
    def section(self) -> str:
        try:  # Return the cached value if possible
            return self._section  # type: ignore
        except AttributeError:
            try:
                self._section = self.job.settings.input.property._h.upper()
            except AttributeError:
                self._section = self.job.settings.input.t.upper()

            return self._section

    def get_energy(self, energy_type: str = "deltag", compound_idx: int = 0, unit: str = "kcal/mol") -> float:
        """Returns the solute solvation energy from an Activity Coefficients calculation."""
        E = self.readkf(self.section, energy_type)[compound_idx]
        return Units.convert(E, "kcal/mol", unit)

    def get_activity_coefficient(self, compound_idx: int = 0) -> float:
        """Return the solute activity coefficient from an Activity Coefficients calculation."""
        return self.readkf(self.section, "gamma")[compound_idx]

    def get_sigma_profile(self, subsection: str = "profil", as_df: bool = False) -> dict:
        r"""Grab all sigma profiles, returning a dictionary of Numpy Arrays.

        Values of :math:`\sigma` are stored under the ``"σ (e/A**2)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`__ package.

        .. _Pandas: https://pandas.pydata.org/

        """
        args = (subsection, "σ (e/A**2)", "chdval")
        try:
            return self._get_array_dict("SIGMAPROFILE", *args, as_df=as_df)
        except KeyError:
            return self._get_array_dict("PURESIGMAPROFILE", *args, as_df=as_df)

    def get_sigma_potential(self, subsection: str = "mu", unit: str = "kcal/mol", as_df: bool = False) -> dict:
        r"""Grab all sigma profiles, expressed in *unit*, and return a dictionary of Numpy Arrays.

        Values of :math:`\sigma` are stored under the ``"σ (e/A**2)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`__ package.

        .. _Pandas: https://pandas.pydata.org/

        """
        args = (subsection, "σ (e/A**2)", "chdval")
        try:
            return self._get_array_dict("SIGMAPOTENTIAL", *args, unit=unit, as_df=as_df)
        except KeyError:
            return self._get_array_dict("PURESIGMAPOTENTIAL", *args, unit=unit, as_df=as_df)

    def get_prop_names(self, section=None) -> list:
        r"""Read the section of the .crskf file and return a list of the properties that were calculated.  The section argument can be supplied to look at previously-calculated results.  If no section name is supplied, the function defaults to using the most recent property that was calculated."""
        if section is None:
            section = self.section
        try:
            return self._kf.get_skeleton()[section]
        except KeyError:
            raise KeyError("Cannot find section name: " + str(section))

    def get_results(self, section=None) -> dict:
        r"""Read the section from the most recent calculation type and return the result as a dictionary."""
        if section is None:
            section = self.section

        if hasattr(self, "_prop_dict") and self._prop_dict["section"] == section:
            return self._prop_dict

        props = self.get_prop_names()
        try:
            props.remove("ncomp")
            props.remove("nitems")
        except ValueError:
            raise ValueError("Results object is missing or incomplete.")

        # first get the two ranges for the indices
        ncomp = self.readkf(section, "ncomp")
        nitems = self.readkf(section, "nitems")
        try:
            nstruct = self.readkf(section, "nstruct")
        except:
            nstruct = ncomp

        np_dict = {"section": section}
        np_dict["ncomp"] = ncomp
        chunk_length = 160
        for prop in props:
            tmp = self.readkf(section, prop)
            if prop in ["filename", "name", "SMILES", "mol_filenames"]:
                if len(tmp) / ncomp == chunk_length:
                    np_dict[prop] = [tmp[i : i + chunk_length].strip() for i in range(0, len(tmp), chunk_length)]
                    continue
                else:
                    np_dict[prop] = tmp.split("\x00")
                    continue
            if prop == "struct names":
                if len(tmp) / nstruct == chunk_length:
                    np_dict[prop] = [tmp[i : i + chunk_length].strip() for i in range(0, len(tmp), chunk_length)]
                    continue
                else:
                    np_dict[prop] = tmp.split("\x00")
                    continue
            if not isinstance(tmp, list):
                np_dict[prop] = tmp
            else:
                np_dict[prop] = np.array(tmp)
                if len(tmp) == ncomp * nitems:
                    np_dict[prop].shape = (ncomp, nitems)

        setattr(self, "_prop_dict", np_dict)
        return np_dict

    def get_multispecies_dist(self):
        """
        This function returns multispecies distribution for each (compound,structure) pair.  The format is a list
        with indices corresponding to compound indices.  Each item in the list is a dictionary with a structure name : list pair, where the structure name corresponds to a structure the compound can be exist as and the list is the distribution of that compound in that structure over the number of points (mole fractions, temperatures, pressures).
        """
        res = self.get_results()
        property_name = res["property"].rstrip()

        if property_name == "LOGP":
            nPhase = 2
        else:
            nPhase = 1

        ncomp = self.readkf(self.section, "ncomp")
        struct_names = res["struct names"]
        num_points = self.readkf(self.section, "nitems")
        valid_structs: List[List[str]] = [[] for _ in range(ncomp)]
        comp_dist = res["comp distribution"].flatten()
        for i in range(len(struct_names)):
            for j in range(ncomp):
                if res["valid structs"][i * ncomp + j]:
                    valid_structs[j].append(struct_names[i])

        compositions: List[Dict[str, List[float]]] = [{vs: [] for vs in valid_structs[i]} for i in range(ncomp)]
        idx = 0
        for i in range(ncomp):
            for nfrac in range(num_points):
                for k in range(nPhase):
                    for j in range(len(valid_structs[i])):
                        compositions[i][valid_structs[i][j]].append(comp_dist[idx])
                        idx += 1

        return compositions

    def get_structure_energy(self, as_df: bool = False):
        """
        Retrieve the energy information for each structure in multispecies.
        If OUTPUT_ENERGY_COMPONENTS is set to True in the input file, this function returns:

        1. The energy of each structure in multispecies (units in kcal/mol).
        2. Information related to association with other compound, if any.

        Parameters:
            as_df (bool, optional): If True, returns the result as a list of Pandas DataFrames.
                                    If False, returns the result as a list of dictionaries.
                                    Default is False.
        Returns:
            List[dict] or List[pandas.DataFrame]: A list containing the energy data and association information for each structure.

        Energy Abbreviations:
            * s_idx: the index for each unique structure
            * CompIdx: the compound index in multispecies
            * FormIdx: the form index in multispecies
            * SpecIdx: the species index in multispecies
            * StrucIdx: the structure index in multispecies
            * z: the equilibrium concentration in multispecies
            * coskf: the corresponding coskf file for each s_idx
            * mu_res: the residual part of the pseudo-chemical potential
            * mu_comb: the combinatorial part of the pseudo-chemical potential
            * mu_disp: the energy contribution from the dispersive interaction
            * mu_pdh: the energy contribution from the Pitzer-Debye-Hückel term
            * mu_RTlnz: the energy contribution from the ideal mixing
            * mu_Ecosmo: the Ecosmo energy
            * mu_res_misfit : the electrostatic interaction in residual part of the pseudo-chemical potential
            * mu_res_hb : the hydrogen bond interaction in residual part of the pseudo-chemical potential
            * Assoc: True if the structure has any association with other compound
            * NumRepMonomer: the number of repeated monomers used for polymers
            * NumStrucPerComp: the number of structures per compound used for dimers, trimers

        Association Information Abbreviations:
            * ReqCompNameAssoc: the required compound name for the associating structure
            * ReqCompIdxAssoc: the required compound index (CompIdx) for the associating structure
            * NumReqCompAssoc: the number of the required compounds in the associating structure
        """

        section = "EnegyComponent"
        try:
            nspecies = self.readkf(section, "nspecies")
        except:
            log("The section of EnergyComponent is not found in the crskf file.")
            return None, None

        ms_index = self.readkf(section, "ms_index")
        ms_index = np.array(ms_index)
        ms_index = ms_index.reshape(nspecies, 4)

        mu_component = self.readkf(section, "mu_component")
        mu_component = np.array(mu_component)
        mu_component = mu_component.reshape(nspecies, int(len(mu_component) / nspecies))

        species_molfrac = self.readkf(section, "species_molfrac")
        species_coskf = self.readkf(section, "species_coskf").split("\x00")
        species_coskf = [os.path.basename(x.rstrip()) for x in species_coskf]

        Assoc = self.readkf(section, "Assoc")
        NumRepMonmer = self.readkf(section, "NumRepMonmer")
        NumStrucPerComp = self.readkf(section, "NumStrucPerComp")

        dict_species = {}
        dict_species["s_idx"] = [i + 1 for i in range(nspecies)]
        dict_species["CompIdx"] = ms_index[:, 0]
        dict_species["FormIdx"] = ms_index[:, 1]
        dict_species["SpecIdx"] = ms_index[:, 2]
        dict_species["StrucIdx"] = ms_index[:, 3]
        dict_species["z"] = species_molfrac
        dict_species["coskf"] = species_coskf
        dict_species["mu_res"] = mu_component[:, 0]
        dict_species["mu_comb"] = mu_component[:, 1]
        dict_species["mu_disp"] = mu_component[:, 2]
        dict_species["mu_pdh"] = mu_component[:, 3]
        dict_species["mu_RTlnz"] = mu_component[:, 4]
        dict_species["mu_Ecosmo"] = mu_component[:, 5]
        dict_species["mu_res_misfit"] = mu_component[:, 6]
        dict_species["mu_res_hb"] = mu_component[:, 7]
        dict_species["Assoc"] = Assoc
        dict_species["NumRepMonmer"] = NumRepMonmer
        dict_species["NumStrucPerComp"] = NumStrucPerComp

        if np.sum(Assoc) > 0:
            Assoc_s_idx = self.readkf(section, "Assoc_s_idx")
            ReqCompIdxAssoc = self.readkf(section, "ReqCompIdxAssoc")
            NumReqCompAssoc = self.readkf(section, "NumReqCompAssoc")
            ReqCompNameAssoc = self.readkf(section, "ReqCompNameAssoc").split("\x00")

            dict_Asson = {}
            dict_Asson["Assoc_s_idx"] = Assoc_s_idx
            dict_Asson["ReqCompIdxAssoc"] = ReqCompIdxAssoc
            dict_Asson["NumReqCompAssoc"] = NumReqCompAssoc
            dict_Asson["ReqCompNameAssoc"] = ReqCompNameAssoc
        else:
            dict_Asson = None

        if as_df:
            try:
                import pandas as pd

                return pd.DataFrame(dict_species), pd.DataFrame(dict_Asson)
            except ImportError:
                method = inspect.stack()[2][3]
                raise ImportError("{}: as_df=True requires the 'pandas' package".format(method))
        else:
            return dict_species, dict_Asson

    def plot(
        self,
        *arrays: "np.ndarray",
        x_axis: Optional[str] = None,
        plot_fig: bool = True,
        x_label: Optional[str] = None,
        y_label: Optional[str] = None,
    ):
        """Plot, show and return a series of COSMO-RS results as a matplotlib Figure instance.

        Accepts the output of, *e.g.*, :meth:`CRSResults.get_sigma_profile`:
        A dictionary of Numpy arrays or a Pandas DataFrame.

        Returns a matplotlib Figure_ instance which can be further modified to the users liking.
        Automatic plotting of the resulting figure can be disabled with the *plot_fig* argument.

        .. note::
            This method requires the `matplotlib <https://matplotlib.org/index.html>`__ package.

        .. note::
            The name of the dictionary/DataFrame key containing the index (*i.e.* the x-axis) can,
            and should, be manually specified in *x_axis* if a custom *x_axis* is passed
            to :meth:`CRSResults._get_array_dict`.
            This argument can be ignored otherwise.

        .. _Figure: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html#matplotlib.figure.Figure

        """  # noqa

        def get_x_axis(array, x_axis):
            """Find and return the index and its name."""
            if x_axis is None:
                return np.arange(array.shape[1])

            if isinstance(x_axis, str):
                ret = self._prop_dict[x_axis]
            else:
                ret = np.array(x_axis, copy=False)
            ret = ret.ravel()  # Flatten it
            return ret[: array.shape[1]]

        # Check running enviroment
        try:
            from IPython import get_ipython

            ipython = get_ipython()
            if ipython is not None:
                if "zmqshell" in str(type(ipython)):
                    terminal = "jupyter"
                else:
                    terminal = "interactive"
            else:
                terminal = "script"
        except ImportError:
            terminal = "script"

        # Check if matplotlib is installed
        try:
            import matplotlib

            if plot_fig:
                if terminal == "jupyter":
                    ipython.run_line_magic("matplotlib", "inline")
                else:
                    matplotlib.use("TkAgg")
            elif not plot_fig:
                matplotlib.use("Agg")

            import matplotlib.pyplot as plt
        except ImportError:
            method = self.__class__.__name__ + ".plot"
            raise ImportError("{}: this method requires the 'matplotlib' package".format(method))

        self.get_results()

        # Create a dictionary of 1d arrays
        array_dict = {}
        for array in arrays:
            name: Optional[str] = None
            if isinstance(array, str):  # Array refers to a section in the kf file
                name = array
                array = self._prop_dict[array]

            # Ensure it's a 2D array
            array = np.array(array, ndmin=2, dtype=float, copy=False)

            # Fill the array dict with 1d arrays
            base_key = "" if name is None else name + " "
            iterator = enumerate(array, 1) if array.shape[0] != 1 else zip(cycle(" "), array)
            for i, array_1d in iterator:
                key = f"{base_key}{i}"
                array_dict[key] = array_1d
        # Retrieve the index and its name
        index = get_x_axis(array, x_axis)
        # print ("INDEX::::", index)
        if x_label is None:
            if isinstance(x_axis, str):
                x_label = x_axis
            else:
                x_label = ""

        if y_label is None:
            y_label = ""

        # Assign various series to the plot
        fig, ax = plt.subplots()
        for k, v in array_dict.items():
            ax.plot(index, v, label=k)

        # Add the legend and x-label
        ax.legend()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        # Show and return
        if plot_fig:
            if terminal == "jupyter":
                pass
            elif terminal == "interactive":
                plt.show(block=False)
            else:
                plt.show()
        return fig

    def _get_array_dict(
        self,
        section: str,
        subsection: str,
        x_axis: str,
        index_subsection: str,
        unit: str = "kcal/mol",
        as_df: bool = False,
    ) -> dict:
        """Create dictionary or DataFrame containing all values in *section*/*subsection*.

        Takes the following arguments:
            * The *section*/*subsection* of the desired quantity.
            * The desired name of the index (*x_axis*).
            * The name of subsection containing the index (*index_subsection*).
            * The *unit* of the output quanty (ignore this keyword if not applicable).
            * If the result should be returned as Pandas DataFrame (*as_df*).

        """
        ret = self._construct_array_dict(section, subsection, unit)

        # Create the index
        index = self.readarray(section, index_subsection, dtype=float)
        if section in ("BINMIXCOEF", "COMPOSITIONLINE", "TERNARYMIX"):
            ncomponent = 3 if section == "TERNARYMIX" else 2
            index.shape = ncomponent, len(index) // ncomponent
            iterator = np.nditer(index.astype(str), flags=["external_loop"], order="F")
            ret[x_axis] = np.array([" / ".join(str(i) for i in item) for item in iterator])
        else:
            ret[x_axis] = index

        # Return a dictionary of arrays or a DataFrame
        if not as_df:
            return ret
        else:
            return self._dict_to_df(ret, section, x_axis)

    def _construct_array_dict(self, section: str, subsection: str, unit: str = "kcal/mol") -> dict:
        """Construct dictionary containing all values in *section*/*subsection*."""
        # Use filenames as keys
        _filenames = self.readkf(section, "filename").split()
        filenames = [_filenames] if not isinstance(_filenames, list) else _filenames

        # Grab the keys and the number of items per key
        keys = [os.path.basename(key) for key in filenames] + ["Total"]
        nitems = self.readkf(section, "nitems")

        # Use sigma profiles/potentials as values
        ratio = Units.conversion_ratio("kcal/mol", unit)
        values = ratio * self.readarray(section, subsection, dtype=float)
        values.shape = len(values) // nitems, nitems

        ret = dict(zip(keys, values))
        try:
            ret["Total"] = self.readarray(section, subsection + "tot", dtype=float)
        except KeyError:
            pass
        return ret

    @staticmethod
    def _dict_to_df(array_dict: dict, section: str, x_axis: str):
        """Attempt to convert a dictionary into a DataFrame."""
        try:
            import pandas as pd
        except ImportError:
            method = inspect.stack()[2][3]
            raise ImportError("{}: as_df=True requires the 'pandas' package".format(method))

        index = pd.Index(array_dict.pop(x_axis), name=x_axis)
        df = pd.DataFrame(array_dict, index=index)
        df.columns.name = section.lower()
        return df


class CRSJob(SCMJob):
    """A |SCMJob| subclass intended for running COSMO-RS jobs."""

    _command = "crs"
    _result_type = CRSResults
    _subblock_end = "end"

    def __init__(self, **kwargs) -> None:
        """Initialize a :class:`CRSJob` instance."""
        super().__init__(**kwargs)
        self.settings.ignore_molecule = True

    @staticmethod
    def database() -> str:
        database_path = os.path.join(os.environ["SCM_PKG_ADFCRSDIR"], "ADFCRS-2018")
        if not os.path.isdir(database_path):
            raise FileNotFoundError("The ADFCRS-2018 database does not seem to be installed")
        return database_path

    @staticmethod
    def coskf_from_database(name: str) -> str:
        if not name.endswith(".coskf"):
            name += ".coskf"
        return os.path.join(CRSJob.database(), name)

    @staticmethod
    def cos_to_coskf(filename: str) -> str:
        """Convert a .cos file into a .coskf file with the :code:`$AMSBIN/cosmo2kf` command.

        Returns the filename of the new .coskf file.

        """
        filename_out = filename + "kf"
        try:
            amsbin = os.environ["AMSBIN"]
        except KeyError:
            raise EnvironmentError(
                "cos_to_coskf: Failed to load 'cosmo2kf' from '$AMSBIN/'; "
                "the 'AMSBIN' environment variable has not been set"
            )

        args = [os.path.join(amsbin, "cosmo2kf"), filename, filename_out]
        subprocess.run(args)
        return filename_out

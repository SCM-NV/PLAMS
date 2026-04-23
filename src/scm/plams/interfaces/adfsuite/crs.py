import inspect
import os
import subprocess
from itertools import cycle
from typing import Optional, List, Dict, TYPE_CHECKING, Set, Union, Any, Tuple, Sequence, cast

import numpy as np

from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.scmjob import SCMJob, SCMResults
from scm.plams.tools.units import Units
from scm.plams.core.functions import log

if TYPE_CHECKING:
    import pandas as pd
    from matplotlib.figure import Figure

__all__ = ["CRSResults", "CRSJob"]


class CRSResults(SCMResults):
    """A |SCMResults| subclass for accessing results of |CRSJob|."""

    _kfext = ".crskf"
    _rename_map = {"CRSKF": "$JN.crskf"}

    @property
    def section(self) -> str:
        try:  # Return the cached value if possible
            return self._section  # type: ignore[has-type]
        except AttributeError:
            try:
                self._section = self.job.settings.input.property._h.upper()
            except AttributeError:
                self._section = self.job.settings.input.t.upper()

            return self._section

    def get_energy(self, energy_type: str = "deltag", compound_idx: int = 0, unit: str = "kcal/mol") -> float:
        """Returns the solute solvation energy from an Activity Coefficients calculation."""
        E = cast(List[float], self.readkf(self.section, energy_type))[compound_idx]
        return Units.convert(E, "kcal/mol", unit)

    def get_activity_coefficient(self, compound_idx: int = 0) -> float:
        """Return the solute activity coefficient from an Activity Coefficients calculation."""
        return cast(List[float], self.readkf(self.section, "gamma"))[compound_idx]

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

    def get_prop_names(self, section: Optional[str] = None) -> Set[str]:
        r"""Read the section of the .crskf file and return a list of the properties that were calculated.  The section argument can be supplied to look at previously-calculated results.  If no section name is supplied, the function defaults to using the most recent property that was calculated."""
        if section is None:
            section = self.section
        try:
            return self._kf.get_skeleton()[section]
        except KeyError:
            raise KeyError("Cannot find section name: " + str(section))

    def get_results(self, section: Optional[str] = None) -> dict:
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
        ncomp = cast(int, self.readkf(section, "ncomp"))
        nitems = cast(int, self.readkf(section, "nitems"))
        try:
            nstruct = cast(int, self.readkf(section, "nstruct"))
        except:
            nstruct = ncomp

        np_dict: Dict[str, Any] = {"section": section}
        np_dict["ncomp"] = ncomp
        chunk_length = 160
        for prop in props:
            tmp = self.readkf(section, prop)
            if prop in ["filename", "name", "SMILES", "mol_filenames"]:
                tmp = cast(str, tmp)
                if len(tmp) / ncomp == chunk_length:
                    np_dict[prop] = [tmp[i : i + chunk_length].strip() for i in range(0, len(tmp), chunk_length)]
                    continue
                else:
                    np_dict[prop] = tmp.split("\x00")
                    continue
            if prop == "struct names":
                tmp = cast(str, tmp)
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

    def get_multispecies_dist(self) -> List[Dict[str, List[float]]]:
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

        ncomp = cast(int, self.readkf(self.section, "ncomp"))
        struct_names = res["struct names"]
        num_points = cast(int, self.readkf(self.section, "nitems"))
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

    def get_structure_energy(self, as_df: bool = False) -> Tuple[Optional[Dict], Optional[Dict]]:
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
            nspecies = cast(int, self.readkf(section, "nspecies"))
        except:
            log("The section of EnergyComponent is not found in the crskf file.")
            return None, None

        ms_index = np.array(self.readkf(section, "ms_index"))
        ms_index = ms_index.reshape(nspecies, 4)

        mu_component = np.array(self.readkf(section, "mu_component"))
        mu_component = mu_component.reshape(nspecies, int(len(mu_component) / nspecies))

        species_molfrac = self.readkf(section, "species_molfrac")
        species_coskf = cast(str, self.readkf(section, "species_coskf")).split("\x00")
        species_coskf = [os.path.basename(x.rstrip()) for x in species_coskf]

        Assoc = self.readkf(section, "Assoc")
        NumRepMonmer = self.readkf(section, "NumRepMonmer")
        NumStrucPerComp = self.readkf(section, "NumStrucPerComp")

        dict_species: Dict[str, Any] = {}
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
            ReqCompNameAssoc = cast(str, self.readkf(section, "ReqCompNameAssoc")).split("\x00")

            dict_Asson: Dict[str, Any] = {}
            dict_Asson["Assoc_s_idx"] = Assoc_s_idx
            dict_Asson["ReqCompIdxAssoc"] = ReqCompIdxAssoc
            dict_Asson["NumReqCompAssoc"] = NumReqCompAssoc
            dict_Asson["ReqCompNameAssoc"] = ReqCompNameAssoc
        else:
            dict_Asson = None  # type: ignore[assignment]

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
    ) -> "Figure":
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

        def get_x_axis(array: np.ndarray, x_axis: Optional[Union[str, np.ndarray]]) -> np.ndarray:
            """Find and return the index and its name."""
            if x_axis is None:
                return np.arange(array.shape[1])

            if isinstance(x_axis, str):
                ret = self._prop_dict[x_axis]  # type: ignore[attr-defined]
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
        _filenames = cast(str, self.readkf(section, "filename")).split()
        filenames = [_filenames] if not isinstance(_filenames, list) else _filenames  # type: ignore[list-item]

        # Grab the keys and the number of items per key
        keys = [os.path.basename(key) for key in filenames] + ["Total"]
        nitems = cast(int, self.readkf(section, "nitems"))

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
    def _dict_to_df(array_dict: dict, section: str, x_axis: str) -> "pd.DataFrame":
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
    PROBLEM_TYPES: Tuple[str, ...] = (
        "ACTIVITYCOEF",
        "BINMIXCOEF",
        "BOILINGPOINT",
        "COMPOSITIONLINE",
        "FLASHPOINT",
        "LLE",
        "LOGP",
        "PUREBOILINGPOINT",
        "PURESIGMAPOTENTIAL",
        "PURESIGMAPROFILE",
        "PURESOLUBILITY",
        "PUREVAPORPRESSURE",
        "SIGMAPOTENTIAL",
        "SIGMAPROFILE",
        "SOLUBILITY",
        "STABILITY",
        "TERNARYMIX",
        "VAPORPRESSURE",
    )
    COMPOUND_KEYS: Dict[str, Dict[str, str]] = {
        "name": {
            "description": "Optional compound name label.",
        },
        "frac1": {
            "description": "Phase-1 mole fraction, or mass fraction when MASSFRACTION is used.",
            "unit": "fraction",
        },
        "frac2": {
            "description": "Phase-2 mole fraction, or mass fraction when MASSFRACTION is used.",
            "unit": "fraction",
        },
        "nring": {
            "description": "COSMO-RS ring-atom count parameter for the compound.",
            "unit": "count",
        },
        "meltingpoint": {
            "description": "Pure-compound melting point for solubility calculations.",
            "unit": "K",
        },
        "hfusion": {
            "description": "Pure-compound enthalpy of fusion for solubility calculations.",
            "unit": "kcal/mol",
        },
        "cpfusion": {
            "description": "Pure-compound heat capacity of fusion for solubility calculations.",
            "unit": "kcal/(mol K)",
        },
        "scalearea": {
            "description": "Expert scaling factor for the COSMO surface area.",
            "unit": "dimensionless",
        },
        "pvap": {
            "description": "Pure-compound vapor pressure used together with tvap.",
            "unit": "bar",
        },
        "tvap": {
            "description": "Temperature corresponding to pvap.",
            "unit": "K",
        },
        "vp_equation": {
            "description": "Vapor-pressure correlation name such as Antoine or VPM1.",
        },
        "vp_params": {
            "description": "Coefficients for the selected vapor-pressure correlation.",
        },
        "density": {
            "description": "Pure-compound density used for solvent-molecule volume calculations.",
            "unit": "kg/L",
        },
        "polymer": {
            "description": "Treat the compound as a polymer using monomer data from the COSMO result file.",
        },
        "averagemwpoly": {
            "description": "Average molecular weight for polymer compounds.",
            "unit": "g/mol",
        },
        "flashpoint": {
            "description": "Pure-compound flash point.",
            "unit": "K",
        },
        "dielectric_const": {
            "description": "Dielectric constant of the solvent.",
        },
        "drophbond": {
            "description": "Disable hydrogen-bond terms for this compound.",
        },
        "cosmofile": {
            "description": "Treat the file as an ASCII .cosmo file instead of a KF-based COSMO result file.",
        },
        "compkffile": {
            "description": "Treat the file as a FastSigma-generated .compkf file.",
        },
        "sigmafile": {
            "description": "Treat the file as an ASCII sigma profile file.",
        },
        "FORM": {
            "description": "Multiple-form settings for the compound, including conformers and associated/dissociated forms."
        },
    }
    _COMPOUND_TEMPLATE_KEYS = tuple(key for key in COMPOUND_KEYS if key != "FORM")
    _COMPOUND_TEMPLATE_KEY_SET = frozenset(_COMPOUND_TEMPLATE_KEYS)
    FORM_KEYS: Dict[str, Dict[str, str]] = {
        "name": {
            "description": "Optional form name label.",
        },
        "count": {
            "description": "Relative count or multiplicity of the form within the compound.",
        },
        "nring": {
            "description": "COSMO-RS ring-atom count parameter for the form.",
            "unit": "count",
        },
        "Hcorr": {
            "description": "Enthalpy correction for this form when modeling multiple forms of a compound.",
            "unit": "kcal/mol",
        },
        "Scorr": {
            "description": "Entropy correction for this form when modeling multiple forms of a compound.",
            "unit": "kcal/mol",
        },
        "drophbond": {
            "description": "Disable hydrogen-bond terms for this form.",
        },
    }
    _FORM_KEY_SET = frozenset(FORM_KEYS)

    def __init__(self, **kwargs: Any) -> None:
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
    def _default_database_coskf(name: str) -> str:
        """Return the intended ADFCRS-2018 path for *name* without requiring a configured AMS installation."""
        if not name.endswith(".coskf"):
            name += ".coskf"
        try:
            return CRSJob.coskf_from_database(name)
        except (FileNotFoundError, KeyError):
            return os.path.join("$SCM_PKG_ADFCRSDIR", "ADFCRS-2018", name)

    @staticmethod
    def compound_template(
        path: str,
        *,
        name: Optional[str] = None,
        frac1: Optional[float] = None,
        frac2: Optional[float] = None,
        nring: Optional[int] = None,
        meltingpoint: Optional[float] = None,
        hfusion: Optional[float] = None,
        cpfusion: Optional[float] = None,
        scalearea: Optional[float] = None,
        pvap: Optional[float] = None,
        tvap: Optional[float] = None,
        vp_equation: Optional[str] = None,
        vp_params: Optional[str] = None,
        density: Optional[float] = None,
        polymer: Optional[bool] = None,
        averagemwpoly: Optional[float] = None,
        flashpoint: Optional[float] = None,
        dielectric_const: Optional[float] = None,
        drophbond: Optional[bool] = None,
        cosmofile: Optional[bool] = None,
        compkffile: Optional[bool] = None,
        sigmafile: Optional[bool] = None,
    ) -> Settings:
        """Create a compound settings block with the given path and explicit COMPOUND keyword/value pairs."""
        compound = Settings()
        compound._h = path

        if nring is None and path.lower().endswith(".coskf"):
            try:
                nring = CRSJob._read_or_determine_nring(path)
            except Exception:
                nring = None

        compound_values = {
            "name": name,
            "frac1": frac1,
            "frac2": frac2,
            "nring": nring,
            "meltingpoint": meltingpoint,
            "hfusion": hfusion,
            "cpfusion": cpfusion,
            "scalearea": scalearea,
            "pvap": pvap,
            "tvap": tvap,
            "vp_equation": vp_equation,
            "vp_params": vp_params,
            "density": density,
            "polymer": polymer,
            "averagemwpoly": averagemwpoly,
            "flashpoint": flashpoint,
            "dielectric_const": dielectric_const,
            "drophbond": drophbond,
            "cosmofile": cosmofile,
            "compkffile": compkffile,
            "sigmafile": sigmafile,
        }
        for key, value in compound_values.items():
            if value is not None:
                compound[key] = value
        return compound

    @staticmethod
    def _compound_template(
        path: str,
        *,
        name: Optional[str] = None,
        frac1: Optional[float] = None,
        frac2: Optional[float] = None,
        nring: Optional[int] = None,
        meltingpoint: Optional[float] = None,
        hfusion: Optional[float] = None,
        cpfusion: Optional[float] = None,
        scalearea: Optional[float] = None,
        pvap: Optional[float] = None,
        tvap: Optional[float] = None,
        vp_equation: Optional[str] = None,
        vp_params: Optional[str] = None,
        density: Optional[float] = None,
        polymer: Optional[bool] = None,
        averagemwpoly: Optional[float] = None,
        flashpoint: Optional[float] = None,
        dielectric_const: Optional[float] = None,
        drophbond: Optional[bool] = None,
        cosmofile: Optional[bool] = None,
        compkffile: Optional[bool] = None,
        sigmafile: Optional[bool] = None,
    ) -> Settings:
        """Backward-compatible wrapper for :meth:`compound_template`."""
        return CRSJob.compound_template(
            path,
            name=name,
            frac1=frac1,
            frac2=frac2,
            nring=nring,
            meltingpoint=meltingpoint,
            hfusion=hfusion,
            cpfusion=cpfusion,
            scalearea=scalearea,
            pvap=pvap,
            tvap=tvap,
            vp_equation=vp_equation,
            vp_params=vp_params,
            density=density,
            polymer=polymer,
            averagemwpoly=averagemwpoly,
            flashpoint=flashpoint,
            dielectric_const=dielectric_const,
            drophbond=drophbond,
            cosmofile=cosmofile,
            compkffile=compkffile,
            sigmafile=sigmafile,
        )

    @staticmethod
    def _normalize_multispecies_form(form: Settings) -> Settings:
        """Validate and normalize a FORM block for use in multispecies compounds."""
        if not isinstance(form, Settings):
            raise TypeError(f"FORM entries must be Settings instances, got {type(form).__name__}")

        path = getattr(form, "_h", None)
        if not path:
            raise ValueError("FORM entries must define a non-empty _h header/path")

        invalid_keys = sorted(set(form.keys()) - CRSJob._FORM_KEY_SET - {"_h"})
        if invalid_keys:
            allowed = ", ".join(CRSJob.FORM_KEYS)
            invalid = ", ".join(invalid_keys)
            raise ValueError(f"Unsupported FORM key(s): {invalid}. Allowed keys: {allowed}")
        return form

    @staticmethod
    def _flatten_multispecies_forms(forms: Sequence[Any]) -> List[Any]:
        """Normalize variadic or list-based FORM input into a flat list of Settings entries."""
        if len(forms) == 1 and isinstance(forms[0], Sequence) and not isinstance(forms[0], Settings):
            flattened = list(forms[0])
        else:
            flattened = list(forms)

        if not flattened:
            raise ValueError("_multispecies_template requires at least one FORM entry")
        return flattened

    @staticmethod
    def multispecies_template(*forms: Any, **compound_kwargs: Any) -> Settings:
        """Create a compound settings block with nested FORM blocks for multispecies compounds."""
        if "FORM" in compound_kwargs:
            raise ValueError("FORM must be passed via the forms argument, not as a COMPOUND keyword")

        invalid_keys = sorted(set(compound_kwargs) - CRSJob._COMPOUND_TEMPLATE_KEY_SET)
        if invalid_keys:
            allowed = ", ".join(CRSJob._COMPOUND_TEMPLATE_KEYS)
            invalid = ", ".join(invalid_keys)
            raise ValueError(f"Unsupported COMPOUND key(s): {invalid}. Allowed keys: {allowed}")

        compound = Settings()
        for key, value in compound_kwargs.items():
            compound[key] = value

        normalized_forms = CRSJob._flatten_multispecies_forms(forms)
        compound.form = [CRSJob._normalize_multispecies_form(form) for form in normalized_forms]
        return compound

    @staticmethod
    def _multispecies_template(*forms: Any, **compound_kwargs: Any) -> Settings:
        """Backward-compatible wrapper for :meth:`multispecies_template`."""
        return CRSJob.multispecies_template(*forms, **compound_kwargs)

    @staticmethod
    def _read_or_determine_nring(coskf_file: str) -> int:
        """Read Nring from a COSKF file, or determine it from the molecular graph if missing."""
        from scm.plams.mol.molecule import Molecule
        from scm.plams.tools.kftools import KFFile

        kf = KFFile(coskf_file)
        try:
            compound_data = kf.read_section("Compound Data")
            nring = compound_data.get("Nring")
            if nring is not None:
                return int(nring)
        except Exception:
            pass

        mol = Molecule(coskf_file)
        rings = mol.locate_rings()
        flatten_atoms = [atom for subring in rings for atom in subring]
        return len(set(flatten_atoms))

    @staticmethod
    def conformers_template(
        path: str,
        *,
        name: Optional[str] = None,
        count: Optional[float] = None,
        nring: Optional[int] = None,
        Hcorr: Optional[float] = None,
        Scorr: Optional[float] = None,
        drophbond: Optional[bool] = None,
    ) -> Settings:
        """Create a single FORM settings block from one conformer COSMO file."""
        if not path:
            raise ValueError("_conformers_template requires a non-empty conformer path")
        form = Settings()
        form._h = path

        if nring is None:
            try:
                nring = CRSJob._read_or_determine_nring(path)
            except Exception:
                nring = None

        form_values = {
            "name": name,
            "count": count,
            "nring": nring,
            "Hcorr": Hcorr,
            "Scorr": Scorr,
            "drophbond": drophbond,
        }
        for key, value in form_values.items():
            if value is not None:
                form[key] = value
        return CRSJob._normalize_multispecies_form(form)

    @staticmethod
    def _conformers_template(
        path: str,
        *,
        name: Optional[str] = None,
        count: Optional[float] = None,
        nring: Optional[int] = None,
        Hcorr: Optional[float] = None,
        Scorr: Optional[float] = None,
        drophbond: Optional[bool] = None,
    ) -> Settings:
        """Backward-compatible wrapper for :meth:`conformers_template`."""
        return CRSJob.conformers_template(
            path,
            name=name,
            count=count,
            nring=nring,
            Hcorr=Hcorr,
            Scorr=Scorr,
            drophbond=drophbond,
        )

    @staticmethod
    def problem_types() -> Tuple[str, ...]:
        """Return the supported COSMO-RS problem types."""
        return CRSJob.PROBLEM_TYPES

    @staticmethod
    def compound_keys() -> Dict[str, Dict[str, str]]:
        """Return metadata for supported COSMO-RS COMPOUND subkeys used on ordinary compounds."""
        return {key: value.copy() for key, value in CRSJob.COMPOUND_KEYS.items()}

    @staticmethod
    def form_keys() -> Dict[str, Dict[str, str]]:
        """Return metadata for FORM subkeys used for multiple forms such as conformers of a compound."""
        return {key: value.copy() for key, value in CRSJob.FORM_KEYS.items()}

    @staticmethod
    def settings_template(problem_type: str) -> Settings:
        """Return a suggested :class:`~scm.plams.core.settings.Settings` template for a COSMO-RS problem type.

        The returned settings contain default ADFCRS-2018 compounds for common workflows and can be edited before use.

        Example:

        .. code:: python

            >>> s = CRSJob.settings_template("LLE")
            >>> job = CRSJob(settings=s)
        """
        if not isinstance(problem_type, str):
            raise TypeError(f"problem_type must be a string, got {type(problem_type).__name__}")

        water = CRSJob._default_database_coskf("Water")
        octanol = CRSJob._default_database_coskf("1-Octanol")
        hexanone = CRSJob._default_database_coskf("2-Hexanone")
        methanol = CRSJob._default_database_coskf("Methanol")
        ethanol = CRSJob._default_database_coskf("Ethanol")
        benzene = CRSJob._default_database_coskf("Benzene")

        normalized = problem_type.upper()
        s = Settings()
        s.input.property._h = normalized

        def compounds(*items: Tuple[str, Dict[str, Any]]) -> List[Settings]:
            return [CRSJob.compound_template(path, **kwargs) for path, kwargs in items]

        pure_types = {"PURESIGMAPROFILE", "PURESIGMAPOTENTIAL", "PUREVAPORPRESSURE", "PUREBOILINGPOINT"}
        ternary_types = {"TERNARYMIX", "STABILITY", "LLE"}
        supported_types = set(CRSJob.problem_types())

        if normalized not in supported_types:
            supported = ", ".join(CRSJob.problem_types())
            raise ValueError(f"Unsupported COSMO-RS problem type '{problem_type}'. Supported types: {supported}")

        if normalized in pure_types:
            s.input.compound = compounds((methanol, {"frac1": 1.0}))
            if normalized == "PUREVAPORPRESSURE":
                s.input.temperature = "273.15 373.15 10"
            elif normalized == "PUREBOILINGPOINT":
                s.input.pressure = "0.101325 1.01325 10"
            return s

        if normalized == "ACTIVITYCOEF":
            s.input.temperature = 298.15
            s.input.compound = compounds(
                (water, {"frac1": 1.0}),
                (benzene, {}),
                (ethanol, {}),
                (methanol, {}),
            )
            return s

        if normalized in {"SIGMAPROFILE", "SIGMAPOTENTIAL"}:
            s.input.temperature = 298.15
            s.input.compound = compounds((water, {"frac1": 0.5}), (ethanol, {"frac1": 0.5}))
            return s

        if normalized == "VAPORPRESSURE":
            s.input.temperature = 298.15
            s.input.compound = compounds(
                (water, {"frac1": 0.25, "pvap": 1.0, "tvap": 373.15}),
                (methanol, {"frac1": 0.25}),
                (ethanol, {"frac1": 0.25, "vp_equation": "Antoine", "vp_params": "5.37229 1670.409 -40.191 0.0 0.0"}),
                (
                    hexanone,
                    {
                        "frac1": 0.25,
                        "vp_equation": "VPM1",
                        "vp_params": "-6474.348470271438 -6.057589837807771 0.003390587477679571 51.07134238467479 0.0",
                    },
                ),
            )
            return s

        if normalized == "FLASHPOINT":
            s.input.massfraction = ""
            s.input.compound = compounds((ethanol, {"frac1": 0.442, "flashpoint": 286.0}), (water, {"frac1": 0.558}))
            return s

        if normalized == "BINMIXCOEF":
            s.input.temperature = 298.14
            s.input.property.Nfrac = 50
            s.input.property.isotherm = ""
            s.input.compound = compounds((water, {"frac1": 0.5}), (methanol, {"frac1": 0.5}))
            return s

        if normalized == "BOILINGPOINT":
            s.input.pressure = "0.101325 1.01325 10"
            s.input.compound = compounds(
                (water, {"frac1": 0.25, "pvap": 1.0, "tvap": 373.15}),
                (methanol, {"frac1": 0.25}),
                (ethanol, {"frac1": 0.25, "vp_equation": "Antoine", "vp_params": "5.37229 1670.409 -40.191 0.0 0.0"}),
                (
                    hexanone,
                    {
                        "frac1": 0.25,
                        "vp_equation": "VPM1",
                        "vp_params": "-6474.348470271438 -6.057589837807771 0.003390587477679571 51.07134238467479 0.0",
                    },
                ),
            )
            return s

        if normalized in ternary_types:
            s.input.temperature = 298.15
            if normalized == "TERNARYMIX":
                s.input.property.Nfrac = 20
                s.input.property.isobar = ""
            s.input.compound = compounds((water, {"frac1": 0.4}), (ethanol, {"frac1": 0.4}), (benzene, {"frac1": 0.2}))
            return s

        if normalized == "LOGP":
            s.input.temperature = 298.15
            s.input.property.volumequotient = 6.766
            s.input.compound = compounds(
                (octanol, {"frac1": 0.725, "frac2": 0.0}),
                (water, {"frac1": 0.275, "frac2": 1.0}),
                (benzene, {}),
                (ethanol, {}),
                (methanol, {}),
            )
            return s

        if normalized == "COMPOSITIONLINE":
            s.input.pressure = 1.01325
            s.input.property.Nfrac = 10
            s.input.property.isobar = ""
            s.input.compound = compounds(
                (water, {"frac1": 0.0, "frac2": 0.9}),
                (ethanol, {"frac1": 0.3, "frac2": 0.1}),
                (benzene, {"frac1": 0.7, "frac2": 0.0}),
            )
            return s

        if normalized == "SOLUBILITY":
            s.input.temperature = "273.15 283.15 10"
            s.input.property.DensitySolvent = 1.0
            s.input.compound = compounds((water, {"frac1": 1.0}), (benzene, {"frac1": 0.0, "meltingpoint": 278.7, "hfusion": 2.37}))
            return s

        s.input.temperature = "273.15 373.15 10"
        s.input.compound = compounds((water, {"frac1": 1.0}), (benzene, {"frac1": 0.0}))
        return s

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

import inspect
import os
import subprocess
from itertools import cycle
from typing import (
    Optional,
    List,
    Dict,
    TYPE_CHECKING,
    Set,
    Union,
    Any,
    Tuple,
    Sequence,
    cast,
    Literal,
    NamedTuple,
)
import numpy as np

from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.scmjob import SCMJob, SCMResults
from scm.plams.tools.units import Units
from scm.plams.core.functions import log
from scm.plams.tools.kftools import KFFile

from pathlib import Path

if TYPE_CHECKING:
    import pandas as pd
    from matplotlib.figure import Figure

__all__ = ["CRSResults", "CRSJob", "CRSResultTables"]

ProblemType = Literal[
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
]


class CRSResultTables(NamedTuple):
    """Container returned by :meth:`CRSResults.get_result_table` with ``split=True``."""

    component: "pd.DataFrame"
    mixture: Optional["pd.DataFrame"]
    lle: Optional["pd.DataFrame"]


class CRSResults(SCMResults):
    """A |SCMResults| subclass for accessing results of |CRSJob|."""

    _kfext = ".crskf"
    _rename_map = {"CRSKF": "$JN.crskf"}
    _RESULT_TABLE_UNSUPPORTED_PROPERTIES = {
        "SIGMAPROFILE", "PURESIGMAPROFILE", "SIGMAPOTENTIAL", "PURESIGMAPOTENTIAL",
    }
    _RESULT_TABLE_LLE_PROPERTIES = {"LLE", "STABILITY", "BINMIXCOEF", "TERNARYMIX"}
    _RESULT_TABLE_COMPONENT_BASE_COLUMNS = ("property", "mixture", "cid", "name",)
    _RESULT_TABLE_LLE_BASE_COLUMNS = ("property", "mixture", "tie_line", "cid", "name",)

    _RESULT_TABLE_COMPONENT_DEFAULT_QUANTITIES = (
        "frac1",
        "frac2",
        "solvent fraction",
        "composition molar fraction",
        "gamma",
        "logp",
        "vapor pressure",
        "henryc",
        "henrycnodim",
        "deltag",
        "solubility molar fraction",
        "solubility massfrac",
        "solubility mol_per_L_solvent",
        "solubility g_per_L_solvent",
        "solubility mol_per_L_solution",
        "solubility g_per_L_solution",
    )
    _RESULT_TABLE_COMPONENT_EXTRA_QUANTITIES = (
        "mu",
        "mu pure",
        "mu gas",
        "mu in solvent 1",
        "mu in solvent 2",
        "xI0",
        "xII0",
        "E gas",
        "G solute",
        "gamma_wf",
        "gamma_vf",
        "fh_chi",
        "poly fraction",
        "polyrepeats",
    )
    _RESULT_TABLE_COMPONENT_KNOWN_QUANTITIES = tuple(
        dict.fromkeys(
            _RESULT_TABLE_COMPONENT_DEFAULT_QUANTITIES
            + _RESULT_TABLE_COMPONENT_EXTRA_QUANTITIES
        )
    )
    _RESULT_TABLE_COMPONENT_EXCLUDED_BY_PROPERTY = {
        "LLE": ("gamma",),
        "STABILITY": ("gamma",),
    }
    _RESULT_TABLE_MIXTURE_DEFAULT_QUANTITIES = (
        "temperature",
        "pressure",
        "excess G",
        "excess H",
        "Gibbs energy of mixing",
        "Enthalpy of vaporization",
        "showmiscgap",
        "unstable",
        "converged",
        "phiI",
        "phiII",
        "tpd_w",
        "isobar",
        "flashpoint",
    )
    _RESULT_TABLE_MIXTURE_EXTRA_QUANTITIES = (
        "Gibbs energy",
        "status_msg",
        "llle_detected",
        "xI_unstable",
        "xII_unstable",
        "tpd_w_I",
        "tpd_w_II",
        "llle_source_phase",
        "llle_status_msg",
        "conv_code",
        "L_conv_code",
        "L_status_msg",
    )
    _RESULT_TABLE_MIXTURE_KNOWN_QUANTITIES = tuple(
        dict.fromkeys(
            _RESULT_TABLE_MIXTURE_DEFAULT_QUANTITIES
            + _RESULT_TABLE_MIXTURE_EXTRA_QUANTITIES
        )
    )

    _RESULT_TABLE_LLE_COLUMN_SPECS = {
        "x": {"quantities": ("x",), "kind": "component"},
        "temperature": {"quantities": ("temperature", "xlle"), "kind": "mixture"},
        "xI": {"quantities": ("xI", "xlle", "xll"), "kind": "component"},
        "xII": {"quantities": ("xII", "xlle", "xll"), "kind": "component"},
        "gammaI": {"quantities": ("gammaI",), "kind": "component"},
        "gammaII": {"quantities": ("gammaII",), "kind": "component"},
        "actI": {"quantities": ("actI",), "kind": "component"},
        "actII": {"quantities": ("actII",), "kind": "component"},
        "phiI": {"quantities": ("phiI",), "kind": "mixture"},
        "phiII": {"quantities": ("phiII",), "kind": "mixture"},
        "converged": {"quantities": ("converged",), "kind": "mixture"},
        "act_interp": {"quantities": ("xlle", "actxll"), "kind": "derived"},
        "pressure": {"quantities": ("xlle", "pressure"), "kind": "derived"},
        "w_min": {"quantities": ("w_min",), "kind": "component", "extra": True},
        "status_msg": {"quantities": ("status_msg",), "kind": "mixture", "extra": True},
    }
    _RESULT_TABLE_LLE_DEFAULT_QUANTITIES = tuple(
        dict.fromkeys(
            quantity
            for spec in _RESULT_TABLE_LLE_COLUMN_SPECS.values()
            if not spec.get("extra")
            for quantity in spec["quantities"]
        )
    )
    _RESULT_TABLE_LLE_EXTRA_QUANTITIES = tuple(
        dict.fromkeys(
            quantity
            for spec in _RESULT_TABLE_LLE_COLUMN_SPECS.values()
            if spec.get("extra")
            for quantity in spec["quantities"]
        )
    )
    _RESULT_TABLE_LLE_KNOWN_QUANTITIES = tuple(
        dict.fromkeys(
            _RESULT_TABLE_LLE_DEFAULT_QUANTITIES
            + _RESULT_TABLE_LLE_EXTRA_QUANTITIES
        )
    )
    _RESULT_TABLE_QUANTITY_METADATA = {
        "frac1": {"symbol": "x", "name": "Feed molar composition", "unit": "fraction"},
        "x": {"symbol": "x", "name": "Feed molar composition", "unit": "fraction"},
        "frac2": {"symbol": "x_2", "name": "Second feed molar composition", "unit": "fraction"},
        "composition molar fraction": {"symbol": "x_calc", "name": "Calculated molar composition", "unit": "fraction"},
        "solvent fraction": {"symbol": "x_solvent", "name": "Solvent fraction", "unit": "fraction"},
        "poly fraction": {"symbol": "x_poly", "name": "Polymer fraction", "unit": "fraction"},
        "gamma": {"symbol": "gamma", "name": "Activity coefficient", "unit": "dimensionless"},
        "gammaI": {"symbol": "gamma_I", "name": "Activity coefficient in phase I", "unit": "dimensionless"},
        "gammaII": {"symbol": "gamma_II", "name": "Activity coefficient in phase II", "unit": "dimensionless"},
        "gamma_wf": {"symbol": "gamma_wf", "name": "Weight-fraction activity coefficient", "unit": "dimensionless"},
        "gamma_vf": {"symbol": "gamma_vf", "name": "Volume-fraction activity coefficient", "unit": "dimensionless"},
        "logp": {"symbol": "logP", "name": "Log10 partition coefficient", "unit": "dimensionless"},
        "fh_chi": {"symbol": "chi_FH", "name": "Flory-Huggins chi", "unit": "dimensionless"},
        "mu": {"symbol": "mu", "name": "Pseudo-chemical potential", "unit": "kcal/mol", "note": "Not a true chemical potential and does not include the ideal mixing term."},
        "mu pure": {"symbol": "mu_pure", "name": "Pure liquid-phase pseudo-chemical potential", "unit": "kcal/mol"},
        "mu gas": {"symbol": "mu_gas", "name": "Gas-phase pseudo-chemical potential", "unit": "kcal/mol", "note": "Gas-phase reference pseudo-chemical potential, optionally refined using input vapor pressure data."},
        "mu in solvent 1": {"symbol": "mu_solv_1", "name": "Pseudo-chemical potential in solvent 1", "unit": "kcal/mol"},
        "mu in solvent 2": {"symbol": "mu_solv_2", "name": "Pseudo-chemical potential in solvent 2", "unit": "kcal/mol"},
        "E gas": {"symbol": "E_gas", "name": "Gas-phase energy", "unit": "kcal/mol"},
        "G solute": {"symbol": "G_solute", "name": "Solute pseudo-energy", "unit": "kcal/mol", "note": "Computed as gas-phase energy plus solvation free energy."},
        "excess G": {"symbol": "G_excess", "name": "Excess Gibbs energy", "unit": "kcal/mol"},
        "excess H": {"symbol": "H_excess", "name": "Excess enthalpy", "unit": "kcal/mol"},
        "Gibbs energy": {"symbol": "G_CRS", "name": "Pseudo-Gibbs energy", "unit": "kcal/mol", "note": "Computed as a composition-weighted sum of pseudo-chemical potentials plus the ideal mixing term; not a true thermodynamic Gibbs energy."},
        "Gibbs energy of mixing": {"symbol": "DeltaG_mix", "name": "Gibbs energy of mixing", "unit": "kcal/mol"},
        "Enthalpy of vaporization": {"symbol": "DeltaH_vap", "name": "Enthalpy of vaporization", "unit": "kcal/mol"},
        "deltag": {"symbol": "DeltaG_solv", "name": "Solvation free energy", "unit": "kcal/mol", "note": "Computed as mu - mu_gas plus a gas-to-solution standard-state Gibbs energy correction; affected by input density or solvent density."},
        "vapor pressure": {"symbol": "p_vap", "name": "Vapor pressure", "unit": "bar"},
        "henryc": {"symbol": "H", "name": "Henry's law constant", "unit": "mol/(L atm)", "note": "Henry's law constant computed as exp((mu_gas - mu)/RT) divided by the solvent molar volume."},
        "henrycnodim": {"symbol": "H_cc", "name": "Dimensionless Henry's law constant", "unit": "dimensionless"},
        "temperature": {"symbol": "T", "name": "Temperature", "unit": "K"},
        "pressure": {"symbol": "P", "name": "Pressure", "unit": "bar"},
        "isobar": {"symbol": "isobar", "name": "Isobaric", "unit": ""},
        "flashpoint": {"symbol": "flashpoint", "name": "Flash-point calculation", "unit": ""},
        "showmiscgap": {"symbol": "misc_gap", "name": "Miscibility gap detected", "unit": "", "note": "Estimated by interpolation."},
        "unstable": {"symbol": "unstable", "name": "TPD unstable", "unit": "", "note": "Indicates instability from the tangent-plane distance test."},
        "converged": {"symbol": "converged", "name": "LLE converged", "unit": ""},
        "status_msg": {"symbol": "status", "name": "Status message", "unit": ""},
        "w_min": {"symbol": "w_min", "name": "Trial-phase composition minimizing TPD(w)", "unit": "fraction"},
        "tpd_w": {"symbol": "TPD_min", "name": "Minimum tangent-plane distance", "unit": ""},
        "phiI": {"symbol": "phi_I", "name": "Phase I fraction", "unit": "fraction"},
        "phiII": {"symbol": "phi_II", "name": "Phase II fraction", "unit": "fraction"},
        "actI": {"symbol": "a_I", "name": "Activity in phase I", "unit": "dimensionless"},
        "actII": {"symbol": "a_II", "name": "Activity in phase II", "unit": "dimensionless"},
        "act_interp": {
            "symbol": "a*",
            "name": "Interpolated activity",
            "unit": "dimensionless",
            "note": "Estimated from interpolated miscibility-gap data; not from a full LLE calculation.",
        },
        "xI": {"symbol": "x_I", "name": "Molar composition in phase I", "unit": "fraction"},
        "xII": {"symbol": "x_II", "name": "Molar composition in phase II", "unit": "fraction"},
        "solution molar fraction": {"symbol": "x_sol", "name": "Solubility molar fraction", "unit": "fraction", "note": "Equilibrium molar-fraction solubility; density is only used for volume-based solubilities."},
        "solubility mol_per_L_solvent": {
            "symbol": "S_mol_L_solvent",
            "name": "Solubility in moles per L solvent",
            "unit": "mol/L solvent",
        },
        "solubility g_per_L_solvent": {
            "symbol": "S_g_L_solvent",
            "name": "Solubility in grams per L solvent",
            "unit": "g/L solvent",
        },
        "solubility mol_per_L_solution": {
            "symbol": "S_mol_L_solution",
            "name": "Solubility in moles per L solution",
            "unit": "mol/L solution",
        },
        "solubility g_per_L_solution": {
            "symbol": "S_g_L_solution",
            "name": "Solubility in grams per L solution",
            "unit": "g/L solution",
        },
        "solubility massfrac": {"symbol": "w_solubility", "name": "Solubility mass fraction", "unit": "fraction"},
    }

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
        np_dict["nitems"] = nitems
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

    def get_result_table(
        self,
        section: Optional[str] = None,
        quantities: Union[str, Sequence[str]] = "default",
        split: bool = False,
        column_labels: Literal["raw", "key", "symbol", "name", "name_unit"] = "name",
    ) -> Union["pd.DataFrame", CRSResultTables]:
        """Return CRS property results as user-friendly pandas dataframe(s).

        By default, this method returns one combined dataframe with one row per
        ``(mixture, cid)`` pair. With ``split=True`` it returns separate component,
        mixture, and LLE tables. Sigma profile and sigma potential sections are
        intentionally handled by :meth:`get_sigma_profile` and
        :meth:`get_sigma_potential` instead. Quantity selection always uses raw
        CRS result keys; ``column_labels`` only controls the returned dataframe
        column names.
        """
        pd = self._import_pandas(self.__class__.__name__ + ".get_result_table")

        column_labels = self._normalize_column_label_mode(column_labels)
        results = self.get_results(section)
        property_name = self._result_property_name(results, section)
        if property_name in self._RESULT_TABLE_UNSUPPORTED_PROPERTIES:
            raise NotImplementedError(
                "{} results are not supported by get_result_table(); use get_sigma_profile() "
                "or get_sigma_potential() instead.".format(property_name)
            )

        component_quantities, mixture_quantities, lle_quantities = self._select_result_table_quantities(
            results, property_name, quantities
        )
        component_names = self._resolve_component_names(results, int(results["ncomp"]))
        component = self._build_result_component_table(pd, results, property_name, component_quantities, component_names)
        mixture = self._build_result_mixture_table(pd, results, property_name, mixture_quantities)
        lle = self._build_result_lle_table(pd, results, property_name, lle_quantities, component_names)

        if split:
            return CRSResultTables(
                component=self._rename_column_labels(component, column_labels),
                mixture=self._rename_column_labels(mixture, column_labels),
                lle=self._rename_column_labels(lle, column_labels),
            )

        if mixture is None:
            return self._rename_column_labels(component, column_labels)
        combined = component.merge(mixture, on=["property", "mixture"], how="left")
        return self._rename_column_labels(combined, column_labels)

    def get_result_table_metadata(
        self,
        section: Optional[str] = None,
        quantities: Union[str, Sequence[str]] = "default",
        split: bool = False,
    ) -> Union["pd.DataFrame", CRSResultTables]:
        """Return metadata describing CRS result table quantity columns.

        The returned dataframe describes raw CRS quantity keys, their table type,
        symbol, descriptive name, unit, and optional note. Quantity selection uses
        the same raw CRS result keys and defaults as :meth:`get_result_table`.
        """
        pd = self._import_pandas(self.__class__.__name__ + ".get_result_table_metadata")

        results = self.get_results(section)
        property_name = self._result_property_name(results, section)
        if property_name in self._RESULT_TABLE_UNSUPPORTED_PROPERTIES:
            raise NotImplementedError(
                "{} results are not supported by get_result_table_metadata(); use get_sigma_profile() "
                "or get_sigma_potential() instead.".format(property_name)
            )

        component_quantities, mixture_quantities, lle_quantities = self._select_result_table_quantities(
            results, property_name, quantities
        )
        lle_columns = self._select_lle_columns(results, property_name, lle_quantities)
        component = self._build_quantity_metadata_table(pd, "component", component_quantities)
        mixture = self._build_quantity_metadata_table(pd, "mixture", mixture_quantities)
        lle = self._build_quantity_metadata_table(pd, "lle", lle_columns)

        if split:
            return CRSResultTables(component=component, mixture=mixture, lle=lle)

        return pd.concat([component, mixture, lle], ignore_index=True)

    @staticmethod
    def _normalize_column_label_mode(column_labels: str) -> str:
        if column_labels == "key":
            return "raw"
        allowed = {"raw", "symbol", "name", "name_unit"}
        if column_labels not in allowed:
            raise ValueError("column_labels must be one of: raw, symbol, name, name_unit")
        return column_labels

    @staticmethod
    def _result_property_name(results: dict, section: Optional[str]) -> str:
        property_value = results.get("property", section)
        if property_value is None:
            raise KeyError("Cannot determine CRS property name.")
        return str(property_value).rstrip().upper()

    def _format_column_label(self, key: str, column_labels: str) -> str:
        if column_labels == "raw" or key in {"property", "mixture", "cid", "name", "molmass", "tie_line"}:
            return key

        metadata = self._RESULT_TABLE_QUANTITY_METADATA.get(key)
        if metadata is None:
            return key

        if column_labels == "symbol":
            return metadata.get("symbol") or key
        if column_labels == "name":
            return metadata.get("name") or key

        name = metadata.get("name") or key
        unit = metadata.get("unit")
        return "{} [{}]".format(name, unit) if unit else name

    def _rename_column_labels(self, df: Optional["pd.DataFrame"], column_labels: str) -> Optional["pd.DataFrame"]:
        if df is None:
            return None
        return df.rename(columns={column: self._format_column_label(column, column_labels) for column in df.columns})

    def _build_quantity_metadata_table(self, pd: Any, table: str, quantities: Sequence[str]) -> "pd.DataFrame":
        columns = ["table", "quantity", "symbol", "name", "unit", "note"]
        rows = []
        for quantity in quantities:
            metadata = self._RESULT_TABLE_QUANTITY_METADATA.get(quantity, {})
            rows.append(
                {
                    "table": table,
                    "quantity": quantity,
                    "symbol": metadata.get("symbol") or quantity,
                    "name": metadata.get("name") or quantity,
                    "unit": metadata.get("unit") or "",
                    "note": metadata.get("note") or "",
                }
            )
        return pd.DataFrame(rows, columns=columns)

    def _filter_component_quantities(
        self, property_name: str, quantities: Sequence[str]
    ) -> Tuple[str, ...]:
        excluded = set(self._RESULT_TABLE_COMPONENT_EXCLUDED_BY_PROPERTY.get(property_name, ()))
        return tuple(q for q in quantities if q not in excluded)

    def _select_result_table_quantities(
        self, results: dict, property_name: str, quantities: Union[str, Sequence[str]]
    ) -> Tuple[Tuple[str, ...], Tuple[str, ...], Tuple[str, ...]]:
        available_component, available_mixture, available_lle = self._available_result_table_quantities(
            results, property_name
        )

        if quantities == "default":
            component = tuple(q for q in self._RESULT_TABLE_COMPONENT_DEFAULT_QUANTITIES if q in results)
            component = self._filter_component_quantities(property_name, component)
            mixture = tuple(q for q in self._RESULT_TABLE_MIXTURE_DEFAULT_QUANTITIES if q in results)
            lle = (
                tuple(q for q in self._RESULT_TABLE_LLE_DEFAULT_QUANTITIES if q in results)
                if property_name in self._RESULT_TABLE_LLE_PROPERTIES
                else ()
            )
            return component, mixture, lle

        if quantities == "all":
            component = tuple(q for q in self._RESULT_TABLE_COMPONENT_KNOWN_QUANTITIES if q in results)
            component = self._filter_component_quantities(property_name, component)
            mixture = tuple(q for q in self._RESULT_TABLE_MIXTURE_KNOWN_QUANTITIES if q in results)
            lle = (
                tuple(q for q in self._RESULT_TABLE_LLE_KNOWN_QUANTITIES if q in results)
                if property_name in self._RESULT_TABLE_LLE_PROPERTIES
                else ()
            )
            return component, mixture, lle

        if isinstance(quantities, str):
            raise ValueError("quantities must be 'default', 'all', or a sequence of quantity names.")

        requested = tuple(quantities)
        known = set(available_component) | set(available_mixture) | set(available_lle)
        unknown = sorted(set(requested) - known)
        if unknown:
            available = ", ".join(sorted(known))
            raise KeyError("Unknown CRS result quantity/quantities: {}. Available quantities: {}".format(unknown, available))

        component = tuple(q for q in requested if q in available_component)
        component = self._filter_component_quantities(property_name, component)
        mixture = tuple(q for q in requested if q in available_mixture)
        lle = tuple(q for q in requested if q in available_lle)
        return component, mixture, lle

    def _available_result_table_quantities(
        self, results: dict, property_name: str
    ) -> Tuple[Tuple[str, ...], Tuple[str, ...], Tuple[str, ...]]:
        non_quantity_keys = {
            "section",
            "property",
            "ncomp",
            "nitems",
            "name",
            "filename",
            "molmass",
            "usepolyunits",
            "mixture",
            "nstruct",
            "valid structs",
            "nvalid structs",
            "comp distribution",
            "struct names",
            "ntriangle",
            "triangle",
        }
        component: List[str] = []
        mixture: List[str] = []
        lle: List[str] = []

        for key in results:
            if key in non_quantity_keys:
                continue
            if (
                property_name in self._RESULT_TABLE_LLE_PROPERTIES
                and key in self._RESULT_TABLE_LLE_KNOWN_QUANTITIES
            ):
                lle.append(key)
            if key in self._RESULT_TABLE_MIXTURE_KNOWN_QUANTITIES:
                mixture.append(key)
            if key in self._RESULT_TABLE_COMPONENT_KNOWN_QUANTITIES:
                component.append(key)

        return tuple(component), tuple(mixture), tuple(lle)

    def _resolve_component_names(self, results: dict, ncomp: int) -> List[str]:
        name_cache: Dict[str, str] = {}
        names = []
        raw_names = results.get("name")
        for cid in range(ncomp):
            raw_name = raw_names[cid]
            cache_key = "" if raw_name is None else str(raw_name)
            if cache_key not in name_cache:
                name_cache[cache_key] = self._resolve_compound_name(cache_key)
            names.append(name_cache[cache_key])
        return names

    @staticmethod
    def _resolve_compound_name(name_or_path: str) -> str:
        if not os.path.isfile(name_or_path):
            return name_or_path

        rkf = KFFile(name_or_path)
        try:
            name = str(rkf.read("Compound Data", "IUPAC")).strip()
            if name:
                return name
        except KeyError:
            pass

        try:
            other_name = str(rkf.read("Compound Data", "Other Name")).strip()
            if other_name:
                return other_name.split(";")[0].strip() or os.path.basename(name_or_path)
        except KeyError:
            pass

        return os.path.basename(name_or_path)

    def _build_result_component_table(
        self,
        pd: Any,
        results: dict,
        property_name: str,
        quantities: Sequence[str],
        component_names: Sequence[str],
    ) -> "pd.DataFrame":
        ncomp = int(results["ncomp"])
        nitems = int(results["nitems"])
        molmass = np.asarray(results["molmass"]).reshape(-1)
        rows = []
        for mixture in range(nitems):
            for cid in range(ncomp):
                row: Dict[str, Any] = {
                    "property": property_name,
                    "mixture": mixture,
                    "cid": cid,
                    "name": component_names[cid],
                    "molmass": molmass[cid],
                }
                for quantity in quantities:
                    row[quantity] = self._get_component_quantity_value(results[quantity], cid, mixture, ncomp, nitems)
                rows.append(row)

        # for quantity in quantities:
        #     print("quantity", quantity)
        #     self._get_component_quantity_value(results[quantity], cid, mixture, ncomp, nitems, check=True)

        return pd.DataFrame(rows, columns=list(self._RESULT_TABLE_COMPONENT_BASE_COLUMNS) + ["molmass"] + list(quantities))

    def _build_result_mixture_table(
        self, pd: Any, results: dict, property_name: str, quantities: Sequence[str]
    ) -> Optional["pd.DataFrame"]:
        if not quantities:
            return None

        nitems = int(results["nitems"])
        lle_not_applicable = bool(results.get("isobar", False)) or bool(results.get("flashpoint", False))
        rows = []
        for mixture in range(nitems):
            row: Dict[str, Any] = {"property": property_name, "mixture": mixture}
            for quantity in quantities:
                if quantity == "showmiscgap" and lle_not_applicable:
                    row[quantity] = np.nan
                else:
                    row[quantity] = self._get_mixture_quantity_value(results[quantity], mixture, nitems)
            rows.append(row)

        # for quantity in quantities:
        #     print("quantity", quantity)
        #     self._get_mixture_quantity_value(results[quantity], mixture, nitems, check=True)

        return pd.DataFrame(rows, columns=["property", "mixture"] + list(quantities))

    def _select_lle_columns(
        self, results: dict, property_name: str, quantities: Sequence[str]
    ) -> Tuple[str, ...]:
        if property_name not in self._RESULT_TABLE_LLE_PROPERTIES or not quantities:
            return ()

        requested = set(quantities)
        columns = []
        for column, spec in self._RESULT_TABLE_LLE_COLUMN_SPECS.items():
            if any(q in requested and q in results for q in spec["quantities"]):
                columns.append(column)
        return tuple(columns)

    def _build_result_lle_table(
        self,
        pd: Any,
        results: dict,
        property_name: str,
        quantities: Sequence[str],
        component_names: Sequence[str],
    ) -> Optional["pd.DataFrame"]:
        if property_name not in self._RESULT_TABLE_LLE_PROPERTIES:
            return None

        columns = (
            self._RESULT_TABLE_LLE_BASE_COLUMNS
            + self._select_lle_columns(results, property_name, quantities)
        )

        if not quantities:
            return pd.DataFrame(columns=columns)

        if property_name in {"LLE", "STABILITY"}:
            ncomp = int(results["ncomp"])
            nitems = 1
            mixture = 0
            rows = []

            for cid in range(ncomp):
                row = {
                    "property": property_name,
                    "mixture": mixture,
                    "tie_line": 0,
                    "cid": cid,
                    "name": component_names[cid],
                }
                for column in columns:
                    if column in self._RESULT_TABLE_LLE_BASE_COLUMNS:
                        continue
                    kind = self._RESULT_TABLE_LLE_COLUMN_SPECS[column]["kind"]
                    if kind == "component":
                        row[column] = self._get_component_quantity_value(results[column], cid, mixture, ncomp, nitems)
                    elif kind == "mixture":
                        row[column] = self._get_mixture_quantity_value(results[column], mixture, nitems)
                rows.append(row)
            return pd.DataFrame(rows, columns=columns)

        if property_name == "BINMIXCOEF":
            showmiscgap = bool(results.get("showmiscgap"))
            lle_not_applicable = bool(results.get("isobar", False)) or bool(results.get("flashpoint", False))
            if not showmiscgap or lle_not_applicable:
                return pd.DataFrame(columns=columns)

            if int(results["ncomp"]) != 2:
                raise ValueError("BINMIXCOEF LLE table expects ncomp == 2.")

            xlle = np.asarray(results["xlle"]).ravel()
            xI = [xlle[0], 1.0 - xlle[0]]
            xII = [xlle[1], 1.0 - xlle[1]]
            act_interp = [xlle[6], xlle[7]]
            pressure, temperature = xlle[2], xlle[3]

            rows = [
                {
                    "property": property_name,
                    "mixture": np.nan,
                    "tie_line": 0,
                    "cid": cid,
                    "name": component_names[cid],
                    "temperature": temperature,
                    "pressure": pressure,
                    "xI": xI[cid],
                    "xII": xII[cid],
                    "act_interp": act_interp[cid],
                }
                for cid in range(2)
            ]
            return pd.DataFrame(rows, columns=columns)

        #case for "TERNARYMIX"
        showmiscgap = bool(results.get("showmiscgap"))
        lle_not_applicable = bool(results.get("isobar", False)) or bool(results.get("flashpoint", False))
        if not showmiscgap or lle_not_applicable:
            return pd.DataFrame(columns=columns)

        if int(results["ncomp"]) != 3:
            raise ValueError("TERNARYMIX LLE table expects ncomp == 3.")

        nxll = int(results["nxll"])
        nitems = int(results["nitems"])
        xll = np.asarray(results["xll"], dtype=float).ravel().reshape(nxll, 6)
        actxll = np.asarray(results["actxll"], dtype=float).ravel().reshape(nxll, 3)
        temperature = self._get_mixture_quantity_value(results["temperature"], 0, nitems)

        rows = [
            {
                "property": property_name,
                "mixture": np.nan,
                "tie_line": tie_line,
                "cid": cid,
                "name": component_names[cid],
                "temperature": temperature,
                "xI": xll[tie_line, cid],
                "xII": xll[tie_line, cid + 3],
                "act_interp": actxll[tie_line, cid],
            }
            for tie_line in range(nxll)
            for cid in range(3)
        ]
        return pd.DataFrame(rows, columns=columns)

    @staticmethod
    def _get_component_quantity_value(value: Any, cid: int, mixture: int, ncomp: int, nitems: int, check: bool = False) -> Any:
        array = np.asarray(value)
        if check: print("Component value array shape:", array.shape)
        if array.shape == (ncomp, nitems):
            return array[cid, mixture]
        if array.shape == (ncomp,):
            return array[cid]
        if array.shape == (ncomp * nitems,):
            return array.reshape(ncomp, nitems)[cid, mixture]
        return value

    @staticmethod
    def _get_mixture_quantity_value(value: Any, mixture: int, nitems: int, check: bool = False) -> Any:
        array = np.asarray(value)
        if check: print("Mixture value array shape:", array.shape)
        if array.shape == (nitems,):
            return array[mixture]
        if array.shape == (1,):
            return array[0]
        if array.shape == ():
            return value
        return value

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
            pd = self._import_pandas(inspect.stack()[1][3], "as_df=True requires the 'pandas' package")
            return pd.DataFrame(dict_species), pd.DataFrame(dict_Asson)
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
    def _import_pandas(method: str, requirement: str = "this method requires the 'pandas' package") -> Any:
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("{}: {}".format(method, requirement))
        return pd

    @staticmethod
    def _dict_to_df(array_dict: dict, section: str, x_axis: str) -> "pd.DataFrame":
        """Attempt to convert a dictionary into a DataFrame."""
        pd = CRSResults._import_pandas(inspect.stack()[2][3], "as_df=True requires the 'pandas' package")

        index = pd.Index(array_dict.pop(x_axis), name=x_axis)
        df = pd.DataFrame(array_dict, index=index)
        df.columns.name = section.lower()
        return df


class CRSJob(SCMJob):
    """A |SCMJob| subclass intended for running COSMO-RS jobs."""

    _command = "crs"
    _result_type = CRSResults
    _subblock_end = "end"
    _PROPERTY_TYPE_METADATA: Dict[str, Dict[str, Any]] = {
        "ACTIVITYCOEF": {
            "description": "Activity coefficients in a solvent mixture, with optional Henry-law related inputs.",
            "top_level_keys": ["temperature", "massfraction"],
            "property_keys": ["densitysolvent"],
            "compound_keys": ["frac1", "density", "pvap", "tvap", "vp_equation", "vp_params"],
            "minimal_required_keys": ["temperature", "frac1"],
            "system_scope": "solvent_mixture_with_solutes",
            "notes": [
                "Henry's-law constants depend on solvent molar volume and gas-phase pseudochemical potential.",
                "Solvent molar volume uses densitysolvent first, then pure-compound density/molar-mass data, then COSMO volume estimates.",
                "Vapor-pressure inputs affect Henry's-law constants through gas-phase pseudochemical potential corrections.",
            ],
        },
        "BINMIXCOEF": {
            "description": "Binary-mixture coefficients over a composition range.",
            "top_level_keys": ["temperature", "pressure", "massfraction"],
            "property_keys": ["nfrac", "isotherm", "isobar", "flashpoint"],
            "compound_keys": ["frac1", "pvap", "tvap", "vp_equation", "vp_params", "flashpoint"],
            "minimal_required_keys": ["temperature"],
            "system_scope": "binary_mixture",
            "notes": [
                "Use only one of isotherm, isobar, or flashpoint at a time.",
                "isotherm uses temperature; isobar uses pressure; flashpoint uses pure-compound flashpoint inputs.",
                "The binary mixture is evaluated at  (nfrac+5) compositions.",
                "Vapor-pressure inputs affect VLE/flash results through gas-phase pseudochemical potential corrections.",
            ],
        },
        "BOILINGPOINT": {
            "description": "Boiling temperature of a mixture for a pressure range.",
            "top_level_keys": ["pressure", "massfraction"],
            "property_keys": [],
            "compound_keys": ["frac1", "pvap", "tvap", "vp_equation", "vp_params"],
            "minimal_required_keys": ["temperature"],
            "system_scope": "mixture",
            "notes": [
                "pressure may be a single value or a range specified as '0.1 1.0 10'.",
                "Vapor-pressure inputs affect boiling-point results through gas-phase pseudochemical potential corrections.",
            ],
        },
        "COMPOSITIONLINE": {
            "description": "Composition-line calculation between two endpoint phase compositions.",
            "top_level_keys": ["temperature", "pressure", "massfraction"],
            "property_keys": ["nfrac", "isotherm", "isobar", "flashpoint"],
            "compound_keys": ["frac1", "frac2", "pvap", "tvap", "vp_equation", "vp_params", "flashpoint"],
            "minimal_required_keys": ["temperature"],
            "system_scope": "binary_mixture",
            "notes": [
                "frac1 and frac2 define the endpoint phase compositions.",
                "Use only one of isotherm, isobar, or flashpoint at a time.",
                "The mixture is evaluated at (nfrac+1) compositions.",
                "Vapor-pressure inputs affect VLE/flash results through gas-phase pseudochemical potential corrections.",
            ],
        },
        "FLASHPOINT": {
            "description": "Flash point of a mixture using user-supplied pure-compound flash points.",
            "top_level_keys": ["massfraction"],
            "property_keys": [],
            "compound_keys": ["frac1", "flashpoint", "pvap", "tvap", "vp_equation", "vp_params"],
            "minimal_required_keys": ["frac1"],
            "system_scope": "mixture",
            "notes": [
                "Use pure-compound flashpoint inputs for flammable components.",
                "Vapor-pressure inputs affect flash results through gas-phase pseudochemical potential corrections.",
            ],
        },
        "LLE": {
            "description": "Liquid-liquid equilibrium for a ternary mixture.",
            "top_level_keys": ["temperature", "massfraction"],
            "property_keys": [],
            "compound_keys": ["frac1"],
            "minimal_required_keys": ["temperature", "frac1"],
            "system_scope": "mixture",
        },
        "LOGP": {
            "description": "Partition coefficients between two immiscible solvent phases.",
            "top_level_keys": ["temperature", "massfraction"],
            "property_keys": ["volumequotient"],
            "compound_keys": ["frac1", "frac2", "density"],
            "minimal_required_keys": ["temperature", "frac1", "frac2"],
            "system_scope": "mixture",
            "notes": [
                "volumequotient sets the solvent-1/solvent-2 molar-volume ratio; otherwise it is estimated from density/molar-mass data or COSMO volumes.",
            ],
        },
        "PUREBOILINGPOINT": {
            "description": "Pure-compound boiling point over a pressure range.",
            "top_level_keys": ["pressure"],
            "property_keys": [],
            "compound_keys": ["pvap", "tvap", "vp_equation", "vp_params"],
            "minimal_required_keys": ["pressure"],
            "system_scope": "pure_compounds",
            "notes": [
                "Multiple COMPOUND blocks may be supplied; each compound is treated independently as a pure compound.",
                "pressure may be a single value or a range specified as '0.1 1.0 10'.",
                "Vapor-pressure inputs affect boiling-point results through gas-phase pseudochemical potential corrections.",
            ],
        },
        "PURESIGMAPOTENTIAL": {
            "description": "Sigma potential for pure compounds.",
            "top_level_keys": [],
            "property_keys": ["nprofile", "sigmamax"],
            "compound_keys": ["frac1"],
            "minimal_required_keys": [],
            "system_scope": "pure_compounds",
            "notes": [
                "Multiple COMPOUND blocks may be supplied; each compound is treated independently as a pure compound.",
            ],
        },
        "PURESIGMAPROFILE": {
            "description": "Sigma profile for pure compounds.",
            "top_level_keys": [],
            "property_keys": ["nprofile", "sigmamax"],
            "compound_keys": ["frac1"],
            "minimal_required_keys": [],
            "system_scope": "pure_compounds",
            "notes": [
                "Multiple COMPOUND blocks may be supplied; each compound is treated independently as a pure compound.",
            ],
        },
        "PURESOLUBILITY": {
            "description": "Solubility of a solute in pure solvents over a temperature range.",
            "top_level_keys": ["temperature", "pressure"],
            "property_keys": ["isobar"],
            "compound_keys": [
                "frac1",
                "meltingpoint",
                "hfusion",
                "cpfusion",
                "density",
                "pvap",
                "tvap",
                "vp_equation",
                "vp_params",
            ],
            "minimal_required_keys": ["temperature", "frac1", "meltingpoint", "hfusion"],
            "system_scope": "pure_solvent_with_solute",
            "notes": [
                "temperature may be a single value or a range specified as '273.15 373.15 10'.",
                "Compound density affects solubility reported in volume-based units.",
                "For solid solutes, provide meltingpoint, hfusion, and optionally cpfusion.",
                "For liquid solutes, LLE is usually more appropriate.",
                "For gas solutes, use isobar with pressure set to the solute partial vapor pressure.",
                "Vapor-pressure inputs affect gas-solubility results through gas-phase pseudochemical potential corrections.",
            ],
        },
        "PUREVAPORPRESSURE": {
            "description": "Pure-compound vapor pressure over a temperature range.",
            "top_level_keys": ["temperature"],
            "property_keys": [],
            "compound_keys": ["pvap", "tvap", "vp_equation", "vp_params"],
            "minimal_required_keys": ["temperature"],
            "system_scope": "pure_compounds",
            "notes": [
                "Multiple COMPOUND blocks may be supplied; each compound is treated independently as a pure compound.",
                "temperature may be a single value or a range specified as '273.15 373.15 10'.",
                "Vapor-pressure inputs affect vapor-pressure results through gas-phase pseudochemical potential corrections.",
            ],
        },
        "SIGMAPOTENTIAL": {
            "description": "Sigma potential for a solvent mixture.",
            "top_level_keys": ["temperature", "massfraction"],
            "property_keys": ["nprofile", "sigmamax"],
            "compound_keys": ["frac1"],
            "minimal_required_keys": ["temperature", "frac1"],
            "system_scope": "mixture",
        },
        "SIGMAPROFILE": {
            "description": "Sigma profile for a solvent mixture.",
            "top_level_keys": ["temperature", "massfraction"],
            "property_keys": ["nprofile", "sigmamax"],
            "compound_keys": ["frac1"],
            "minimal_required_keys": ["temperature", "frac1"],
            "system_scope": "mixture",
        },
        "SOLUBILITY": {
            "description": "Solubility of solutes in a solvent mixture or under gas-pressure conditions.",
            "top_level_keys": ["temperature", "pressure", "massfraction"],
            "property_keys": ["densitysolvent", "isobar"],
            "compound_keys": [
                "frac1",
                "meltingpoint",
                "hfusion",
                "cpfusion",
                "density",
                "pvap",
                "tvap",
                "vp_equation",
                "vp_params",
            ],
            "minimal_required_keys": ["temperature", "frac1", "meltingpoint", "hfusion"],
            "system_scope": "solvent_with_solutes",
            "notes": [
                "temperature may be a single value or a range specified as '273.15 373.15 10'.",
                "densitysolvent or compound density affects solubility reported in volume-based units.",
                "For solid solutes, provide meltingpoint, hfusion, and optionally cpfusion.",
                "For liquid solutes, LLE is usually more appropriate.",
                "For gas solutes, use isobar with pressure set to the solute partial vapor pressure.",
                "Vapor-pressure inputs affect gas-solubility results through gas-phase pseudochemical potential corrections.",
            ],
        },
        "STABILITY": {
            "description": "Michelsen tangent-plane-distance stability test for a feed composition.",
            "top_level_keys": ["temperature", "massfraction"],
            "property_keys": [],
            "compound_keys": ["frac1"],
            "minimal_required_keys": ["temperature", "frac1"],
            "system_scope": "mixture",
        },
        "TERNARYMIX": {
            "description": "Ternary mixture property sweep over composition space.",
            "top_level_keys": ["temperature", "pressure", "massfraction"],
            "property_keys": ["nfrac", "isotherm", "isobar", "flashpoint"],
            "compound_keys": ["frac1", "pvap", "tvap", "vp_equation", "vp_params", "flashpoint"],
            "minimal_required_keys": ["temperature"],
            "system_scope": "ternary_mixture",
            "notes": [
                "Use only one of isotherm, isobar, or flashpoint at a time.",
                "isotherm uses temperature; isobar uses pressure; flashpoint uses pure-compound flashpoint inputs.",
                "The ternary mixture is evaluated at (nfrac+1)*(nfrac+2)/2 compositions.",
                "Vapor-pressure inputs affect VLE/flash results through gas-phase pseudochemical potential corrections.",
            ],
        },
        "VAPORPRESSURE": {
            "description": "Vapor pressure of a mixture at fixed temperature.",
            "top_level_keys": ["temperature", "massfraction"],
            "property_keys": [],
            "compound_keys": ["frac1", "pvap", "tvap", "vp_equation", "vp_params"],
            "minimal_required_keys": ["temperature", "frac1"],
            "system_scope": "mixture",
            "notes": [
                "temperature may be a single value or a range specified as '273.15 373.15 10'.",
                "Vapor-pressure inputs affect mixture vapor-pressure results through gas-phase pseudochemical potential corrections.",
            ],
        },
    }
    _PROPERTY_KEY_DEFAULTS = {
        "isotherm": True,
        "isobar": False,
        "flashpoint": False,
        "nfrac": 10,
        "nprofile": 50,
        "sigmamax": 0.025,
    }
    _METHODS = (
        "COSMORS",
        "COSMO-RS",
        "COSMOSAC",
        "COSMOSAC2013",
        "COSMOSAC2016",
        "COSMOSACDHB",
        "COSMOSACDHB-MESP",
    )
    _METHOD_ALIASES = {"COSMORS": "COSMO-RS", "COSMOSAC": "COSMOSAC2013"}
    _DEFAULT_DISPERSION_PARAMETERS = {
        "H": -0.0340,
        "C": -0.0356,
        "N": -0.0224,
        "O": -0.0333,
        "F": -0.026,
        "Si": -0.04,
        "P": -0.045,
        "S": -0.052,
        "Cl": -0.0485,
        "Br": -0.055,
        "I": -0.062,
    }
    _METHOD_PARAMETERS_METADATA: Dict[str, Dict[str, Dict[str, Any]]] = {
        "COSMO-RS": {
            "CRSParameters": {
                "rav": 0.400,
                "aprime": 1510.0,
                "fcorr": 2.802,
                "chb": 8850.0,
                "sigmahbond": 0.00854,
                "aeff": 6.94,
                "lambda": 0.130,
                "omega": -0.212,
                "eta": -9.65,
                "chortf": 0.816,
                "hb_hnof": True,
                "hb_temp": True,
                "combi2005": True,
            },
            "Dispersion": _DEFAULT_DISPERSION_PARAMETERS,
        },
        "COSMOSAC2013": {
            "SACParameters": {
                "aeff": 6.4813,
                "sigma0": 0.01233,
                "qn": 79.532,
                "aes": 7877.13,
                "cohoh": 5786.72,
                "cotot": 2739.58,
                "cohot": 4707.75,
                "rav": 0.51,
                "qs": 0.57,
                "hb_hnof": True,
                "hb_notemp": True,
            },
            "Epsilon": {
                "H": 338.13,
                "C.sp3": 29160.92,
                "C.sp2": 30951.83,
                "C.sp": 20685.98,
                "N.sp3": 23488.54,
                "N.sp2": 22663.38,
                "N.sp": 6390.40,
                "O.sp3-H": 8527.06,
                "O.sp3": 8484.38,
                "O.sp2": 6736.85,
                "O.sp2-N": 12145.28,
                "F": 8435.13,
                "P": 82512.21,
                "S": 56067.81,
                "Cl": 45065.19,
                "Br": 62947.83,
                "I": 105910.88,
            },
        },
        "COSMOSAC2016": {
            "SACParameters": {
                "aeff": 5.8447,
                "fdecay": 3.57,
                "sigma0": 0.007,
                "rn": 66.69,
                "qn": 79.53,
                "aes": 5920.84,
                "bes": 1.3950e8,
                "cohoh": 3551.1,
                "cotot": 1077.26,
                "cohot": 3099.31,
                "omega": -0.212,
                "eta": -9.00,
                "hb_hnof": True,
                "hb_notemp": True,
            },
            "Dispersion": _DEFAULT_DISPERSION_PARAMETERS,
        },
        "COSMOSACDHB": {
            "SACParameters": {
                "aeff": 5.8447,
                "fdecay": 3.57,
                "sigma0": 0.0063,
                "rn": 66.69,
                "qn": 79.53,
                "aes": 5920.84,
                "bes": 1.3950e8,
                "cohoh": 33306.83,
                "cotot": 33306.83,
                "cohot": 33306.83,
                "rhbcut": 1.4432,
                "omega": -0.212,
                "eta": -9.00,
                "hb_hnof": True,
                "hb_notemp": True,
            },
            "Dispersion": _DEFAULT_DISPERSION_PARAMETERS,
        },
        "COSMOSACDHB-MESP": {
            "SACParameters": {
                "aeff": 5.8447,
                "fdecay": 3.57,
                "sigma0": 0.0063,
                "rn": 66.69,
                "qn": 79.53,
                "aes": 5920.84,
                "bes": 1.3950e8,
                "cohoh": 34234.17,
                "cotot": 34234.17,
                "cohot": 34234.17,
                "rhbcut": 1.3871,
                "omega": -0.212,
                "eta": -9.00,
                "hb_hnof": True,
                "hb_notemp": True,
            },
            "Dispersion": _DEFAULT_DISPERSION_PARAMETERS,
        },
    }

    _COMPOUND_KEY_METADATA: Dict[str, Dict[str, str]] = {
        "name": {
            "description": "Optional compound name label.",
            "type": "str",
        },
        "frac1": {
            "description": "Phase-1 mole fraction, or mass fraction when MASSFRACTION is used.",
            "type": "float",
            "unit": "fraction",
        },
        "frac2": {
            "description": "Phase-2 mole fraction, or mass fraction when MASSFRACTION is used.",
            "type": "float",
            "unit": "fraction",
        },
        "nring": {
            "description": "COSMO-RS ring-atom count parameter for the compound.",
            "type": "int",
            "unit": "count",
        },
        "meltingpoint": {
            "description": "Pure-compound melting point for solubility calculations.",
            "type": "float",
            "unit": "K",
        },
        "hfusion": {
            "description": "Pure-compound enthalpy of fusion for solubility calculations.",
            "type": "float",
            "unit": "kcal/mol",
        },
        "cpfusion": {
            "description": "Pure-compound heat capacity of fusion for solubility calculations.",
            "type": "float",
            "unit": "kcal/(mol K)",
        },
        "scalearea": {
            "description": "Expert scaling factor for the COSMO surface area.",
            "type": "float",
            "unit": "dimensionless",
        },
        "pvap": {
            "description": "Pure-compound vapor pressure used together with tvap.",
            "type": "float",
            "unit": "bar",
        },
        "tvap": {
            "description": "Temperature corresponding to pvap.",
            "type": "float",
            "unit": "K",
        },
        "vp_equation": {
            "description": "Vapor-pressure correlation name such as Antoine or VPM1.",
            "type": "str",
        },
        "vp_params": {
            "description": "Coefficients for the selected vapor-pressure correlation.",
            "type": "str",
        },
        "density": {
            "description": "Pure-compound density used for solvent-molecule volume calculations.",
            "type": "float",
            "unit": "kg/L",
        },
        "polymer": {
            "description": "Treat the compound as a polymer using monomer data from the COSMO result file.",
            "type": "bool",
        },
        "averagemwpoly": {
            "description": "Average molecular weight for polymer compounds.",
            "type": "float",
            "unit": "g/mol",
        },
        "flashpoint": {
            "description": "Pure-compound flash point.",
            "type": "float",
            "unit": "K",
        },
        "dielectric_const": {
            "description": "Dielectric constant of the solvent.",
            "type": "float",
        },
        "drophbond": {
            "description": "Disable hydrogen-bond terms for this compound.",
            "type": "bool",
        },
        "cosmofile": {
            "description": "Treat the file as an ASCII .cosmo file instead of a KF-based COSMO result file.",
            "type": "bool",
        },
        "compkffile": {
            "description": "Treat the file as a FastSigma-generated .compkf file.",
            "type": "bool",
        },
        "sigmafile": {
            "description": "Treat the file as an ASCII sigma profile file.",
            "type": "bool",
        },
        "form": {
            "description": "Multiple-form settings for the compound, including conformers and associated/dissociated forms.",
            "type": "block",
        },
    }
    _compound_block_KEYS = tuple(key for key in _COMPOUND_KEY_METADATA if key != "form")
    _compound_block_KEY_SET = frozenset(_compound_block_KEYS)
    _FORM_KEY_METADATA: Dict[str, Dict[str, str]] = {
        "name": {
            "description": "Optional form name label.",
            "type": "str",
        },
        "count": {
            "description": "Relative count of this form within the compound, e.g. form A: 1, form B: 1.",
            "type": "float",
        },
        "nring": {
            "description": "COSMO-RS ring-atom count parameter for the form.",
            "type": "int",
            "unit": "count",
        },
        "Hcorr": {
            "description": "Enthalpy correction for this form when modeling multiple forms of a compound.",
            "type": "float",
            "unit": "kcal/mol",
        },
        "Scorr": {
            "description": "Entropy correction for this form when modeling multiple forms of a compound.",
            "type": "float",
            "unit": "kcal/mol",
        },
        "drophbond": {
            "description": "Disable hydrogen-bond terms for this form.",
            "type": "bool",
        },
        "species": {
            "description": "Optional nested species list for multispecies forms.",
            "type": "block",
        },
    }
    _FORM_KEY_SET = frozenset(_FORM_KEY_METADATA)
    _SPECIES_KEY_METADATA: Dict[str, Dict[str, str]] = {
        "name": {
            "description": "Optional species name label.",
            "type": "str",
        },
        "count": {
            "description": "Stoichiometric count of this species within the form, e.g. Mg2+: 1, Cl-: 2.",
            "type": "float",
        },
        "nring": {
            "description": "COSMO-RS ring-atom count parameter for the species.",
            "type": "int",
            "unit": "count",
        },
        "Hcorr": {
            "description": "Enthalpy correction for this species.",
            "type": "float",
            "unit": "kcal/mol",
        },
        "Scorr": {
            "description": "Entropy correction for this species.",
            "type": "float",
            "unit": "kcal/mol",
        },
        "drophbond": {
            "description": "Disable hydrogen-bond terms for this species.",
            "type": "bool",
        },
        "structure": {
            "description": "One or more nested structure entries for this species.",
            "type": "block",
        },
    }
    _SPECIES_KEY_SET = frozenset(_SPECIES_KEY_METADATA)
    _STRUCTURE_KEY_METADATA: Dict[str, Dict[str, str]] = {
        "name": {
            "description": "Optional structure name label.",
            "type": "str",
        },
        "count": {
            "description": "Relative count of this structure within the species, e.g. structure A: 1, structure B: 1",
            "type": "float",
        },
        "nring": {
            "description": "COSMO-RS ring-atom count parameter for the structure.",
            "type": "int",
            "unit": "count",
        },
        "Hcorr": {
            "description": "Enthalpy correction for this structure.",
            "type": "float",
            "unit": "kcal/mol",
        },
        "Scorr": {
            "description": "Entropy correction for this structure.",
            "type": "float",
            "unit": "kcal/mol",
        },
        "drophbond": {
            "description": "Disable hydrogen-bond terms for this structure.",
            "type": "bool",
        },
    }
    _STRUCTURE_KEY_SET = frozenset(_STRUCTURE_KEY_METADATA)

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
    def property_types() -> Tuple[str, ...]:
        """Return the supported COSMO-RS problem types."""
        return tuple(CRSJob._PROPERTY_TYPE_METADATA)

    @staticmethod
    def _normalize_choice(
        value: str,
        *,
        value_name: str,
        allowed: Sequence[str],
        aliases: Optional[Dict[str, str]] = None,
        supported_values: Optional[Sequence[str]] = None,
    ) -> str:
        """Normalize and validate a string choice against supported values."""
        if not isinstance(value, str):
            raise TypeError(f"{value_name} must be a string, got {type(value).__name__}")

        normalized = value.upper()
        if aliases is not None:
            normalized = aliases.get(normalized, normalized)
        if normalized not in allowed:
            supported = ", ".join(supported_values if supported_values is not None else allowed)
            raise ValueError(f"Unsupported {value_name} {value!r}. Supported values: {supported}")
        return normalized

    @staticmethod
    def _normalize_property_type(property_type: str) -> str:
        """Normalize and validate a COSMO-RS problem type."""
        return CRSJob._normalize_choice(
            property_type,
            value_name="property_type",
            allowed=CRSJob._PROPERTY_TYPE_METADATA,
            supported_values=CRSJob.property_types(),
        )

    @staticmethod
    def property_type_metadata(
        property_type: ProblemType,
        as_summary: bool = False,
    ) -> Union[Dict[str, Any], Tuple[str, ...]]:
        """Return discoverability metadata for a COSMO-RS problem type."""
        normalized = CRSJob._normalize_property_type(property_type)
        data = CRSJob._PROPERTY_TYPE_METADATA[normalized]
        metadata = {
            key: value.copy() if isinstance(value, list) else value
            for key, value in data.items()
        }
        # metadata = copy.deepcopy(CRSJob._PROPERTY_TYPE_METADATA[normalized])
        if not as_summary:
            return metadata

        summary = [
            f"{normalized}: {metadata['description']}",
            f"system_scope: {metadata['system_scope']}",
            f"top_level_keys: {', '.join(metadata['top_level_keys']) or '-'}",
            f"property_keys: {', '.join(metadata['property_keys']) or '-'}",
            f"compound_keys: {', '.join(metadata['compound_keys']) or '-'}",
            f"minimal_required_keys: {', '.join(metadata['minimal_required_keys']) or '-'}",
        ]
        for note in metadata.get("notes", []):
            summary.append(f"note: {note}")
        return tuple(summary)

    @staticmethod
    def _property_type_keys(property_type: str) -> Set[str]:
        """Return PROPERTY-block keys supported by a normalized property type."""
        return set(CRSJob._PROPERTY_TYPE_METADATA[property_type]["property_keys"])

    @staticmethod
    def _log_unsupported_property_options(property_type: str, options: Dict[str, Any]) -> None:
        """Log provided property-block options that are not supported for *property_type*."""
        supported_keys = CRSJob._property_type_keys(property_type)
        unsupported = [
            option_name
            for option_name, value in options.items()
            if value is not None and option_name not in supported_keys
        ]
        if unsupported:
            log(f"problem type '{property_type}' ignores {', '.join(unsupported)}", level=3)

    @staticmethod
    def property_block(
        property_type: ProblemType,
        *,
        include_defaults: bool = False,
        volumequotient: Optional[float] = None,
        densitysolvent: Optional[float] = None,
        nfrac: Optional[int] = None,
        isotherm: Optional[bool] = None,
        isobar: Optional[bool] = None,
        flashpoint: Optional[bool] = None,
        nprofile: Optional[int] = None,
        sigmamax: Optional[float] = None,
    ) -> Settings:
        """Create a PROPERTY block for a COSMO-RS problem type.

        By default, only the property type and explicitly provided options are written.
        If an option is omitted, CRS uses the default from the input definition.
        Pass ``include_defaults=True`` to write the documented PROPERTY defaults explicitly.

        ``densitysolvent`` and ``volumequotient`` are presence-sensitive: if omitted,
        CRS estimates volume information from COSMO volumes, so they are written only
        when explicitly provided.
        """
        normalized = CRSJob._normalize_property_type(property_type)
        s = Settings()
        s.input.property._h = normalized

        options = {
            "volumequotient": volumequotient,
            "densitysolvent": densitysolvent,
            "nfrac": nfrac,
            "isotherm": isotherm,
            "isobar": isobar,
            "flashpoint": flashpoint,
            "nprofile": nprofile,
            "sigmamax": sigmamax,
        }
        CRSJob._log_unsupported_property_options(normalized, options)
        property_keys = CRSJob._property_type_keys(normalized)

        def get_value(key: str, value: Any) -> Any:
            if value is not None:
                return value
            if include_defaults:
                return CRSJob._PROPERTY_KEY_DEFAULTS.get(key)
            return None

        if "densitysolvent" in property_keys and densitysolvent is not None:
            s.input.property.densitysolvent = densitysolvent

        if "volumequotient" in property_keys:
            if volumequotient is not None:
                s.input.property.volumequotient = volumequotient
            else:
                log(
                    f"volumequotient is not provided for problem type '{property_type}'; "
                    "the phase volume ratio will be estimated automatically from density or COSMO volume data.",
                    level=3,
                )

        if "nfrac" in property_keys:
            nfrac_value = get_value("nfrac", nfrac)
            if nfrac_value is not None:
                s.input.property.nfrac = nfrac_value

            mode_values = {
                "isotherm": isotherm,
                "isobar": isobar,
                "flashpoint": flashpoint,
            }
            if include_defaults and all(value is None for value in mode_values.values()):
                mode_values["isotherm"] = CRSJob._PROPERTY_KEY_DEFAULTS["isotherm"]
            mode_flags = (
                ("isotherm", mode_values["isotherm"]),
                ("isobar", mode_values["isobar"]),
                ("flashpoint", mode_values["flashpoint"]),
            )
            selected_flags = [name for name, enabled in mode_flags if enabled]
            if len(selected_flags) > 1:
                raise ValueError("Use only one of isotherm, isobar, or flashpoint at a time")
            if selected_flags:
                s.input.property[selected_flags[0]] = ""

        if "nprofile" in property_keys and "sigmamax" in property_keys:
            nprofile_value = get_value("nprofile", nprofile)
            sigmamax_value = get_value("sigmamax", sigmamax)
            if nprofile_value is not None:
                s.input.property.nprofile = nprofile_value
            if sigmamax_value is not None:
                s.input.property.sigmamax = sigmamax_value

        if normalized in {"SOLUBILITY", "PURESOLUBILITY"} and isobar:
            s.input.property.isobar = ""

        return s


    @staticmethod
    def methods() -> Tuple[str, ...]:
        """Return the supported COSMO-RS/SAC methods."""
        return CRSJob._METHODS

    @staticmethod
    def _set_non_none(block: Settings, values: Dict[str, Any]) -> Settings:
        for key, value in values.items():
            if value is not None:
                block[key] = value
        return block

    @staticmethod
    def _normalize_method(method: str) -> str:
        """Normalize and validate a COSMO-RS/SAC method name."""
        return CRSJob._normalize_choice(
            method,
            value_name="method",
            allowed=CRSJob._METHODS,
            aliases=CRSJob._METHOD_ALIASES,
            supported_values=CRSJob.methods(),
        )

    @staticmethod
    def _set_method_parameters(
        block: Settings,
        parameters: Dict[str, Any],
        block_name: str,
        method_name: str,
    ) -> None:
        """Set validated method parameters on a CRS method subblock."""

        method_metadata = CRSJob._METHOD_PARAMETERS_METADATA.get(method_name, {})
        if block_name not in method_metadata:
            raise ValueError(f"{block_name} is not supported for method {method_name!r}")

        allowed_keys = set(method_metadata[block_name])

        if block_name in {"CRSParameters", "SACParameters"}:
            allowed_keys.update({"hb_all", "hb_hnof", "hb_temp", "hb_notemp"})

        if block_name == "CRSParameters":
            allowed_keys.update({"combi1998", "combi2005"})

        for key, value in parameters.items():
            normalized_key = key.lower()
            if normalized_key not in allowed_keys:
                allowed = ", ".join(sorted(allowed_keys))
                raise ValueError(f"Unsupported {block_name} key {key!r}. Allowed keys: {allowed}")
            if value is not None:
                block[normalized_key] = value


    @staticmethod
    def _copy_method_parameter_defaults(method: str) -> Dict[str, Dict[str, Any]]:
        """Return a copy of default method parameters for a normalized method."""
        defaults = CRSJob._METHOD_PARAMETERS_METADATA.get(method, {})
        return {block_name: values.copy() for block_name, values in defaults.items()}

    @staticmethod
    def method_block(
        method: str = "COSMO-RS",
        *,
        include_defaults: bool = False,
        crsparameters: Optional[Dict[str, Any]] = None,
        sacparameters: Optional[Dict[str, Any]] = None,
        dispersion: Optional[Dict[str, float]] = None,
        epsilon: Optional[Dict[str, float]] = None,
    ) -> Settings:
        """Return method-selection and method-parameter settings for COSMO-RS/SAC.
        The method name and validated parameter names are case-insensitive.
        """
        s = Settings()
        normalized_method = CRSJob._normalize_method(method)
        s.input.method = normalized_method

        method_parameter_defaults = CRSJob._copy_method_parameter_defaults(normalized_method) if include_defaults else {}

        for block_name, user_parameters in (
            ("CRSParameters", crsparameters),
            ("SACParameters", sacparameters),
        ):
            parameters = {
                **method_parameter_defaults.get(block_name, {}),
                **(user_parameters or {}),
            }
            if parameters:
                CRSJob._set_method_parameters(
                    s.input[block_name],
                    parameters,
                    block_name,
                    normalized_method,
                )

        for block_name, user_parameters in (
            ("Dispersion", dispersion),
            ("Epsilon", epsilon),
        ):
            parameters = {
                **method_parameter_defaults.get(block_name, {}),
                **(user_parameters or {}),
            }
            # parameters = {key: value for key, value in parameters.items() if value is not None}
            if parameters:
                CRSJob._set_non_none(s.input[block_name], parameters)
                # block = s.input[block_name]
                # for key, value in parameters.items():
                #     block[key] = value

        return s

    @staticmethod
    def _infer_nring_if_missing(path: str, nring: Optional[int]) -> Optional[int]:
        if nring is not None:
            return nring
        try:
            return CRSJob._read_or_determine_nring(path)
        except Exception:
            return None


    @staticmethod
    def compound_block(
        path: Optional[str] = None,
        *,
        forms: Optional[Union[Settings, Sequence[Settings]]] = None,
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
        """Create a COMPOUND settings block from a compound path or FORM entries.

        Pass either ``path`` for a single-file compound or ``forms`` for a compound
        with nested FORM blocks.

        Example::

            benzene = CRSJob.compound_block(
                "benzene.coskf",
                name="benzene",
                frac1=0.3,
                meltingpoint=278.7,
                hfusion=2.37,
            )

            form0 = CRSJob.form_block("conf_0.coskf", Hcorr=0.0)
            form1 = CRSJob.form_block("conf_1.coskf", Hcorr=0.5)
            compound = CRSJob.compound_block(
                forms=[form0, form1],
                name="solute",
                frac1=1.0,
            )
        """

        if path is not None and forms is not None:
            raise ValueError("compound_block accepts either a compound path or FORM entries, not both")
        if path is None and forms is None:
            raise ValueError("compound_block requires either a compound path or FORM entries")
        if path == "":
            raise ValueError("compound_block requires a non-empty compound path")

        compound = Settings()

        if path is not None:
            compound._h = path
            nring = CRSJob._infer_nring_if_missing(path, nring)
        else:
            compound.form = [
                CRSJob._normalize_form(form)
                for form in CRSJob._ensure_settings_list(forms, "forms")
            ]

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

        CRSJob._set_non_none(compound, compound_values)

        return compound

    @staticmethod
    def form_block(
        path: Optional[str] = None,
        *,
        species: Optional[Union[Settings, Sequence[Settings]]] = None,
        name: Optional[str] = None,
        count: Optional[float] = None,
        nring: Optional[int] = None,
        Hcorr: Optional[float] = None,
        Scorr: Optional[float] = None,
        drophbond: Optional[bool] = None,
    ) -> Settings:
        """Create a FORM settings block from a form path or SPECIES entries.

        Pass either ``path`` for a single-file form or ``species`` for a form with
        nested SPECIES blocks.

        Example::

            neutral = CRSJob.form_block("lowest_energy_conformer.coskf")

            mg = CRSJob.species_block("mg.coskf", name="Mg2+", count=1.0)
            cl = CRSJob.species_block("cl.coskf", name="Cl-", count=2.0)
            ion = CRSJob.form_block(
                species=[mg, cl],
                name="MgCl2 dissociated ion",
                count=1.0,
            )
        """
        if path is not None and species is not None:
            raise ValueError("form_block accepts either a form path or SPECIES entries, not both")
        if path is None and species is None:
            raise ValueError("form_block requires either a form path or SPECIES entries")
        if path == "":
            raise ValueError("form_block requires a non-empty form path")

        form = Settings()

        if path is not None:
            form._h = path
            nring = CRSJob._infer_nring_if_missing(path, nring)
        else:
            form.species = [
                CRSJob._normalize_species(species_entry)
                for species_entry in CRSJob._ensure_settings_list(species, "species")
            ]

        form_values = {
            "name": name,
            "count": count,
            "nring": nring,
            "Hcorr": Hcorr,
            "Scorr": Scorr,
            "drophbond": drophbond,
        }

        CRSJob._set_non_none(form, form_values)

        return CRSJob._normalize_form(form)

    @staticmethod
    def species_block(
        path: Optional[str] = None,
        *,
        structures: Optional[Union[Settings, Sequence[Settings]]] = None,
        name: Optional[str] = None,
        count: Optional[float] = None,
        nring: Optional[int] = None,
        Hcorr: Optional[float] = None,
        Scorr: Optional[float] = None,
        drophbond: Optional[bool] = None,
    ) -> Settings:
        """Create a SPECIES settings block from a species path or STRUCTURE entries.

        Pass either ``path`` for a single-file species or ``structures`` for a species
        with nested STRUCTURE blocks.

        Example::

            single_anion = CRSJob.species_block("poly_anion.coskf", name="poly_anion", count=1.0)

            structure_a = CRSJob.structure_block("poly_anion_a.coskf", Hcorr=0.0)
            structure_b = CRSJob.structure_block("poly_anion_b.coskf", Hcorr=0.4)
            multi_anion = CRSJob.species_block(
                structures=[structure_a, structure_b],
                name="poly_anion",
                count=1.0,
            )
        """
        if path is not None and structures is not None:
            raise ValueError("species_block accepts either a species path or STRUCTURE entries, not both")
        if path is None and structures is None:
            raise ValueError("species_block requires either a species path or STRUCTURE entries")
        if path == "":
            raise ValueError("species_block requires a non-empty species path")

        species = Settings()

        if path is not None:
            species._h = path
            nring = CRSJob._infer_nring_if_missing(path, nring)
        else:
            species.structure = [
                CRSJob._normalize_structure(structure_entry)
                for structure_entry in CRSJob._ensure_settings_list(structures, "structures")
            ]

        species_values = {
            "name": name,
            "count": count,
            "nring": nring,
            "Hcorr": Hcorr,
            "Scorr": Scorr,
            "drophbond": drophbond,
        }
        CRSJob._set_non_none(species, species_values)

        return CRSJob._normalize_species(species)

    @staticmethod
    def structure_block(
        path: str,
        *,
        name: Optional[str] = None,
        count: Optional[float] = None,
        nring: Optional[int] = None,
        Hcorr: Optional[float] = None,
        Scorr: Optional[float] = None,
        drophbond: Optional[bool] = None,
    ) -> Settings:
        """Create a STRUCTURE settings block from a structure path.

        Example::

            structure_a = CRSJob.structure_block(
                "poly_anion_a.coskf",
                name="poly_anion conformer A",
                count=1.0,
                Hcorr=0.0
            )
        """
        if not path:
            raise ValueError("structure_block requires a non-empty structure path")
        structure = Settings()
        structure._h = path

        structure_values = {
            "name": name,
            "count": count,
            "nring": CRSJob._infer_nring_if_missing(path, nring),
            "Hcorr": Hcorr,
            "Scorr": Scorr,
            "drophbond": drophbond,
        }

        CRSJob._set_non_none(structure, structure_values)

        return CRSJob._normalize_structure(structure)

    @staticmethod
    def _ensure_settings_list(value: Any, block_name: str) -> List[Settings]:
        """Normalize a single Settings object or a sequence of Settings objects into a non-empty list."""
        if isinstance(value, Settings):
            items = [value]
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            items = list(value)
        else:
            raise TypeError(f"{block_name} must be a Settings instance or a sequence of Settings instances")

        if not items:
            raise ValueError(f"{block_name} must contain at least one entry")

        invalid = [type(item).__name__ for item in items if not isinstance(item, Settings)]
        if invalid:
            raise TypeError(f"{block_name} entries must be Settings instances, got {', '.join(invalid)}")
        return items

    @staticmethod
    def _normalize_form(form: Settings) -> Settings:
        """Validate and normalize a FORM block for use in multispecies compounds."""
        if not isinstance(form, Settings):
            raise TypeError(f"FORM entries must be Settings instances, got {type(form).__name__}")

        path = dict.get(form, "_h", None)
        has_path = bool(path)
        has_species = "species" in form and form["species"] not in (None, False)
        if has_path == has_species:
            raise ValueError("FORM entries must define exactly one of _h or species")

        invalid_keys = sorted(set(form.keys()) - CRSJob._FORM_KEY_SET - {"_h"})
        if invalid_keys:
            allowed = ", ".join(CRSJob._FORM_KEY_METADATA)
            invalid = ", ".join(invalid_keys)
            raise ValueError(f"Unsupported FORM key(s): {invalid}. Allowed keys: {allowed}")

        if has_species:
            form.species = [CRSJob._normalize_species(species) for species in CRSJob._ensure_settings_list(form.species, "species")]
        return form

    @staticmethod
    def _normalize_species(species: Settings) -> Settings:
        """Validate and normalize a SPECIES block nested under FORM."""
        if not isinstance(species, Settings):
            raise TypeError(f"SPECIES entries must be Settings instances, got {type(species).__name__}")

        path = dict.get(species, "_h", None)
        has_path = bool(path)
        has_structure = "structure" in species and species["structure"] not in (None, False)
        if has_path == has_structure:
            raise ValueError("SPECIES entries must define exactly one of _h or structure")

        invalid_keys = sorted(set(species.keys()) - CRSJob._SPECIES_KEY_SET - {"_h"})
        if invalid_keys:
            allowed = ", ".join(CRSJob._SPECIES_KEY_METADATA)
            invalid = ", ".join(invalid_keys)
            raise ValueError(f"Unsupported SPECIES key(s): {invalid}. Allowed keys: {allowed}")

        if has_structure:
            species.structure = [CRSJob._normalize_structure(structure) for structure in CRSJob._ensure_settings_list(species.structure, "structure")]
        return species

    @staticmethod
    def _normalize_structure(structure: Settings) -> Settings:
        """Validate and normalize a STRUCTURE block nested under SPECIES."""
        path = getattr(structure, "_h", None)
        if not path:
            raise ValueError("STRUCTURE entries must define a non-empty _h header/path")

        invalid_keys = sorted(set(structure.keys()) - CRSJob._STRUCTURE_KEY_SET - {"_h"})
        if invalid_keys:
            allowed = ", ".join(CRSJob._STRUCTURE_KEY_METADATA)
            invalid = ", ".join(invalid_keys)
            raise ValueError(f"Unsupported STRUCTURE key(s): {invalid}. Allowed keys: {allowed}")
        return structure

    @staticmethod
    def _read_or_determine_nring(coskf_file: str) -> int:
        """Read Nring from a COSKF file, or determine it from the molecular graph if missing."""
        from scm.plams.mol.molecule import Molecule
        from scm.plams.tools.kftools import KFFile

        if coskf_file.lower().endswith(".coskf"):
            kf = KFFile(coskf_file)

            try:
                compound_data = kf.read_section("Compound Data")
                nring = compound_data.get("Nring")
                if nring is not None:
                    return int(nring)

            except Exception as exc:
                log(f"Could not read Nring from COSKF file {coskf_file}; determining it from the molecular graph instead: {exc}", level=3)

        elif coskf_file.lower().endswith(".cosmo"):
            cosmo_file = coskf_file
            coskf_file = str(Path(cosmo_file).with_suffix(".coskf"))
            coskf_file = CRSJob.cos_to_coskf(cosmo_file, coskf_file)

        else:
            log(f"Nring determination is supported only for COSKF and COSMO files: {coskf_file}", level=3)
            return 0

        try:
            mol = Molecule(coskf_file)
            rings = mol.locate_rings()
            flatten_atoms = [atom for subring in rings for atom in subring]
            return len(set(flatten_atoms))
        except Exception as exc:
            log(f"Failed to determine Nring from molecular graph for {coskf_file}: {exc}", level=3)
            return 0

    @staticmethod
    def _format_key_metadata(key: str, metadata: Dict[str, str]) -> str:
        """Return a compact display string for one metadata entry."""
        type_unit = metadata.get("type", "")
        if metadata.get("unit"):
            type_unit = f"{type_unit} [{metadata['unit']}]"
        return f"{key}: {type_unit} - {metadata['description']}"

    @staticmethod
    def compound_keys(as_summary: bool = False) -> Union[Dict[str, Dict[str, str]], Tuple[str, ...]]:
        """Return supported COMPOUND subkey metadata, or compact display strings."""
        metadata = {key: value.copy() for key, value in CRSJob._COMPOUND_KEY_METADATA.items()}
        if as_summary:
            return tuple(CRSJob._format_key_metadata(key, value) for key, value in metadata.items())
        return metadata

    @staticmethod
    def form_keys(as_summary: bool = False) -> Union[Dict[str, Dict[str, str]], Tuple[str, ...]]:
        """Return metadata for FORM subkeys used for multiple forms such as conformers of a compound."""
        metadata = {key: value.copy() for key, value in CRSJob._FORM_KEY_METADATA.items()}
        if as_summary:
            return tuple(CRSJob._format_key_metadata(key, value) for key, value in metadata.items())
        return metadata

    @staticmethod
    def species_keys(as_summary: bool = False) -> Union[Dict[str, Dict[str, str]], Tuple[str, ...]]:
        """Return metadata for SPECIES subkeys nested under FORM."""
        metadata = {key: value.copy() for key, value in CRSJob._SPECIES_KEY_METADATA.items()}
        if as_summary:
            return tuple(CRSJob._format_key_metadata(key, value) for key, value in metadata.items())
        return metadata

    @staticmethod
    def structure_keys(as_summary: bool = False) -> Union[Dict[str, Dict[str, str]], Tuple[str, ...]]:
        """Return metadata for STRUCTURE subkeys nested under SPECIES."""
        metadata = {key: value.copy() for key, value in CRSJob._STRUCTURE_KEY_METADATA.items()}
        if as_summary:
            return tuple(CRSJob._format_key_metadata(key, value) for key, value in metadata.items())
        return metadata

    @staticmethod
    def job_settings_template(property_type: ProblemType) -> Settings:
        """Return a full default :class:`~scm.plams.core.settings.Settings` object for a COSMO-RS problem type."""
        water = CRSJob._default_database_coskf("Water")
        octanol = CRSJob._default_database_coskf("1-Octanol")
        hexanone = CRSJob._default_database_coskf("2-Hexanone")
        methanol = CRSJob._default_database_coskf("Methanol")
        ethanol = CRSJob._default_database_coskf("Ethanol")
        benzene = CRSJob._default_database_coskf("Benzene")

        normalized = CRSJob._normalize_property_type(property_type)

        if normalized == "ACTIVITYCOEF":
            s = CRSJob.property_block(normalized)
            s.input.temperature = 298.15
            s.input.compound = [
                CRSJob.compound_block(water, frac1=1.0),
                CRSJob.compound_block(benzene),
                CRSJob.compound_block(ethanol),
                CRSJob.compound_block(methanol),
            ]
            return s

        if normalized == "LOGP":
            s = CRSJob.property_block(normalized, volumequotient=6.766)
            s.input.temperature = 298.15
            s.input.compound = [
                CRSJob.compound_block(octanol, frac1=0.725, frac2=0.0),
                CRSJob.compound_block(water, frac1=0.275, frac2=1.0),
                CRSJob.compound_block(benzene),
                CRSJob.compound_block(ethanol),
                CRSJob.compound_block(methanol),
            ]
            return s

        if normalized == "SOLUBILITY":
            s = CRSJob.property_block(normalized)
            s.input.temperature = "273.15 283.15 10"
            s.input.property.DensitySolvent = 1.0
            s.input.compound = [
                CRSJob.compound_block(water, frac1=1.0),
                CRSJob.compound_block(benzene, frac1=0.0, meltingpoint=278.7, hfusion=2.37),
            ]
            return s

        if normalized in {"VAPORPRESSURE", "BOILINGPOINT"}:
            s = CRSJob.property_block(normalized)
            if normalized == "VAPORPRESSURE":
                s.input.temperature = 298.15
            else:
                s.input.pressure = "0.101325 1.01325 10"
            s.input.compound = [
                CRSJob.compound_block(water, frac1=0.25, pvap=1.0, tvap=373.15),
                CRSJob.compound_block(methanol, frac1=0.25),
                CRSJob.compound_block(
                    ethanol,
                    frac1=0.25,
                    vp_equation="Antoine",
                    vp_params="5.37229 1670.409 -40.191 0.0 0.0",
                ),
                CRSJob.compound_block(
                    hexanone,
                    frac1=0.25,
                    vp_equation="VPM1",
                    vp_params="-6474.348470271438 -6.057589837807771 0.003390587477679571 51.07134238467479 0.0",
                ),
            ]
            return s

        if normalized == "FLASHPOINT":
            s = CRSJob.property_block(normalized)
            s.input.massfraction = ""
            s.input.compound = [
                CRSJob.compound_block(ethanol, frac1=0.442, flashpoint=286.0),
                CRSJob.compound_block(water, frac1=0.558),
            ]
            return s

        if normalized in {"STABILITY", "LLE"}:
            s = CRSJob.property_block(normalized)
            s.input.temperature = 298.15
            s.input.compound = [
                CRSJob.compound_block(water, frac1=0.4),
                CRSJob.compound_block(ethanol, frac1=0.4),
                CRSJob.compound_block(benzene, frac1=0.2),
            ]
            return s

        if normalized == "BINMIXCOEF":
            s = CRSJob.property_block(normalized, nfrac=50, isotherm=True)
            s.input.temperature = 298.14
            s.input.compound = [
                CRSJob.compound_block(water, frac1=0.5),
                CRSJob.compound_block(methanol, frac1=0.5),
            ]
            return s

        if normalized == "TERNARYMIX":
            s = CRSJob.property_block(normalized, nfrac=20)
            s.input.temperature = 298.15
            s.input.compound = [
                CRSJob.compound_block(water, frac1=0.4),
                CRSJob.compound_block(ethanol, frac1=0.4),
                CRSJob.compound_block(benzene, frac1=0.2),
            ]
            return s

        if normalized == "COMPOSITIONLINE":
            s = CRSJob.property_block(normalized, nfrac=10, isobar=True)
            s.input.pressure = 1.01325
            s.input.compound = [
                CRSJob.compound_block(water, frac1=0.0, frac2=0.9),
                CRSJob.compound_block(ethanol, frac1=0.3, frac2=0.1),
                CRSJob.compound_block(benzene, frac1=0.7, frac2=0.0),
            ]
            return s

        if normalized == "PURESOLUBILITY":
            s = CRSJob.property_block(normalized)
            s.input.temperature = "273.15 373.15 10"
            s.input.compound = [
                CRSJob.compound_block(water, frac1=1.0),
                CRSJob.compound_block(benzene, frac1=0.0),
            ]
            return s

        if normalized == "PUREVAPORPRESSURE":
            s = CRSJob.property_block(normalized)
            s.input.temperature = "273.15 373.15 10"
            s.input.compound = [CRSJob.compound_block(methanol, frac1=1.0)]
            return s

        if normalized == "PUREBOILINGPOINT":
            s = CRSJob.property_block(normalized)
            s.input.pressure = "0.101325 1.01325 10"
            s.input.compound = [CRSJob.compound_block(methanol, frac1=1.0)]
            return s

        if normalized in {"PURESIGMAPROFILE", "PURESIGMAPOTENTIAL"}:
            s = CRSJob.property_block(normalized)
            s.input.compound = [CRSJob.compound_block(methanol, frac1=1.0)]
            return s

        if normalized in {"SIGMAPROFILE", "SIGMAPOTENTIAL"}:
            s = CRSJob.property_block(normalized)
            s.input.temperature = 298.15
            s.input.compound = [
                CRSJob.compound_block(water, frac1=0.5),
                CRSJob.compound_block(ethanol, frac1=0.5),
            ]
            return s

        s = CRSJob.property_block(normalized)
        s.input.temperature = "273.15 373.15 10"
        s.input.compound = [
            CRSJob.compound_block(water, frac1=1.0),
            CRSJob.compound_block(benzene, frac1=0.0),
        ]
        return s

    @staticmethod
    def cos_to_coskf(filename: str, filename_out: Optional[str] = None) -> str:
        """Convert a .cos file into a .coskf file with the :code:`$AMSBIN/cosmo2kf` command.

        Returns the filename of the new .coskf file.

        """
        if filename_out is None:
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

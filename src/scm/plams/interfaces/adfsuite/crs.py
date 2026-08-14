from __future__ import annotations
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
    ClassVar,
    Callable,
    overload,
    TypeVar,
)
import numpy as np

from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.scmjob import SCMJob, SCMResults
from scm.plams.tools.units import Units
from scm.plams.core.functions import log
from scm.plams.tools.kftools import KFFile
from scm.plams.interfaces.adfsuite import crs_definitions as crs_defs
from scm.plams.interfaces.adfsuite.crs_definitions import (
    CRS_DATA,
    get_block_child_names,
    get_block_keys,
)

from pathlib import Path
from copy import deepcopy

if TYPE_CHECKING:
    import pandas as pd
    from matplotlib.figure import Figure


__all__ = [
    "CRSResults",
    "CRSJob",
    "CRSResultTables",
    "CRSMethodName",
    "CRSPropertyName",
    "CRSVaporPressureEquation",
    "CRSSolubilityMode",
    "CRSVLESweepMode",
    "CRSInputBuilder",
    "ACTIVITYCOEFInputBuilder",
    "LOGPInputBuilder",
    "SOLUBILITYInputBuilder",
    "PURESOLUBILITYInputBuilder",
    "VAPORPRESSUREInputBuilder",
    "PUREVAPORPRESSUREInputBuilder",
    "BOILINGPOINTInputBuilder",
    "PUREBOILINGPOINTInputBuilder",
    "FLASHPOINTInputBuilder",
    "BINMIXCOEFInputBuilder",
    "TERNARYMIXInputBuilder",
    "COMPOSITIONLINEInputBuilder",
    "LLEInputBuilder",
    "STABILITYInputBuilder",
    "SIGMAPROFILEInputBuilder",
    "PURESIGMAPROFILEInputBuilder",
    "SIGMAPOTENTIALInputBuilder",
    "PURESIGMAPOTENTIALInputBuilder",
]

PathLike = Union[str, os.PathLike]

CRSMethodName = Literal[
    "COSMO-RS",
    "COSMOSAC2013",
    "COSMOSAC2016",
    "COSMOSACDHB",
    "COSMOSACDHB-MESP",
]
CRSPropertyName = Literal[
    "ACTIVITYCOEF",
    "LOGP",
    "SOLUBILITY",
    "PURESOLUBILITY",
    "VAPORPRESSURE",
    "PUREVAPORPRESSURE",
    "BOILINGPOINT",
    "PUREBOILINGPOINT",
    "FLASHPOINT",
    "BINMIXCOEF",
    "TERNARYMIX",
    "COMPOSITIONLINE",
    "LLE",
    "STABILITY",
    "SIGMAPROFILE",
    "PURESIGMAPROFILE",
    "SIGMAPOTENTIAL",
    "PURESIGMAPOTENTIAL",
]
CRSVaporPressureEquation = Literal[
    "Antoine",
    "VPM1",
    "DIPPR101",
    "DIPPR115",
    "KDB",
]
CRSSolubilityMode = Literal["gas", "liquid", "solid"]
CRSVLESweepMode = Literal["isotherm", "isobar", "flashpoint"]
_CRSInputBuilderT = TypeVar("_CRSInputBuilderT", bound="CRSInputBuilder")

_VAPOR_PRESSURE_KEYS = ("pvap", "tvap", "vp_equation", "vp_params")
_FUSION_KEYS = ("meltingpoint", "hfusion", "cpfusion")
_VLE_SWEEP_PROPERTY_KEYS = ("nfrac", "isotherm", "isobar", "flashpoint")
_SIGMA_PROPERTY_KEYS = ("nprofile", "sigmamax")

_PURE_COMPOUND_COMMENTS = ("Multiple COMPOUND blocks are treated as independent pure compounds.",)

_VLE_SWEEP_MODE_CONFIG = {
    "options": ("isotherm", "isobar", "flashpoint"),
    "default": "isotherm",
    "input_mapping": {
        "isotherm": {"property": {"isotherm": True}},
        "isobar": {"property": {"isobar": True}},
        "flashpoint": {"property": {"flashpoint": True}},
    },
    "descriptions": {
        "isotherm": (
            "Composition scan at constant temperature.",
            "Phase boundaries are interpolated from the scan.",
            "Use LLE for robust phase-boundary calculations.",
        ),
        "isobar": ("Composition scan at constant pressure."),
        "flashpoint": ("Estimate flash point using pure-compound flash points. Vapor-pressure inputs may improve results.",),
    },
}

_SOLUBILITY_MODE_CONFIG = {
    "options": ("gas", "liquid", "solid"),
    "default": "solid",
    "input_mapping": {
        "gas": {"property": {"isobar": True}},
        "liquid": {},
        "solid": {},
    },
    "compound_required_keys": {
        "solid": {
            "solute": {
                "any_of": (("meltingpoint", "hfusion"),),
            },
        },
    },
    "descriptions": {
        "gas": (
            "Calculate gas solubility at fixed partial pressure.",
            "Top-level pressure is the solute partial pressure.",
        ),
        "liquid": (
            "Calculate liquid solubility assuming the pure liquid solute as the coexisting phase.",
            "Use LLE for miscible liquid systems.",
        ),
        "solid": (
            "Calculate solid solubility using solute fusion-correction inputs.",
        ),
    },
}


class CRSInputRoute(NamedTuple):
    """Route from a public CRS builder key to a PLAMS Settings input location."""

    scope: Literal["top_level", "property"]
    metadata: Dict[str, Any]


class CRSModeConfig:
    """Static mode configuration for one CRS property type."""

    options: Tuple[str, ...]
    default: Optional[str]
    input_mapping: Dict[str, Any]
    descriptions: Dict[str, Any]

    def __init__(self, metadata: Optional[Dict[str, Any]] = None) -> None:
        metadata = metadata or {}
        self.options = tuple(metadata.get("options", ()))
        self.default = metadata.get("default")
        self.input_mapping = metadata.get("input_mapping", {})
        self.descriptions = metadata.get("descriptions", {})



class CRSResultTables(NamedTuple):
    """Container returned by :meth:`CRSResults.get_result_table` with ``split=True``."""

    component: "pd.DataFrame"
    mixture: Optional["pd.DataFrame"]
    lle: Optional["pd.DataFrame"]


class CRSSigmaProfileColumn(NamedTuple):
    """Column mapping for sigma-profile KF arrays."""

    column: str
    compound_source: str
    mixture_total_source: Optional[str] = None
    hb_channel: Optional[int] = None


class CRSResults(SCMResults):
    """A |SCMResults| subclass for accessing results of |CRSJob|."""

    _kfext = ".crskf"
    _rename_map = {"CRSKF": "$JN.crskf"}
    _RESULT_TABLE_UNSUPPORTED_PROPERTIES = crs_defs.RESULT_TABLE_UNSUPPORTED_PROPERTIES
    _RESULT_TABLE_LLE_PROPERTIES = crs_defs.RESULT_TABLE_LLE_PROPERTIES
    _RESULT_TABLE_COMPONENT_BASE_COLUMNS = crs_defs.RESULT_TABLE_COMPONENT_BASE_COLUMNS
    _RESULT_TABLE_LLE_BASE_COLUMNS = crs_defs.RESULT_TABLE_LLE_BASE_COLUMNS
    _RESULT_TABLE_COMPONENT_DEFAULT_QUANTITIES = crs_defs.RESULT_TABLE_COMPONENT_DEFAULT_QUANTITIES
    _RESULT_TABLE_COMPONENT_EXTRA_QUANTITIES = crs_defs.RESULT_TABLE_COMPONENT_EXTRA_QUANTITIES
    _RESULT_TABLE_COMPONENT_KNOWN_QUANTITIES = crs_defs.RESULT_TABLE_COMPONENT_KNOWN_QUANTITIES
    _RESULT_TABLE_COMPONENT_EXCLUDED_BY_PROPERTY = crs_defs.RESULT_TABLE_COMPONENT_EXCLUDED_BY_PROPERTY
    _RESULT_TABLE_MIXTURE_DEFAULT_QUANTITIES = crs_defs.RESULT_TABLE_MIXTURE_DEFAULT_QUANTITIES
    _RESULT_TABLE_MIXTURE_EXTRA_QUANTITIES = crs_defs.RESULT_TABLE_MIXTURE_EXTRA_QUANTITIES
    _RESULT_TABLE_MIXTURE_KNOWN_QUANTITIES = crs_defs.RESULT_TABLE_MIXTURE_KNOWN_QUANTITIES
    _RESULT_TABLE_LLE_COLUMN_SPECS = crs_defs.RESULT_TABLE_LLE_COLUMN_SPECS
    _RESULT_TABLE_LLE_DEFAULT_QUANTITIES = crs_defs.RESULT_TABLE_LLE_DEFAULT_QUANTITIES
    _RESULT_TABLE_LLE_EXTRA_QUANTITIES = crs_defs.RESULT_TABLE_LLE_EXTRA_QUANTITIES
    _RESULT_TABLE_LLE_KNOWN_QUANTITIES = crs_defs.RESULT_TABLE_LLE_KNOWN_QUANTITIES
    _RESULT_TABLE_QUANTITY_METADATA = crs_defs.RESULT_TABLE_QUANTITY_METADATA
    _SIGMA_PROFILE_PROPERTIES = ("SIGMAPROFILE", "PURESIGMAPROFILE")
    _SIGMA_PROFILE_DEFAULT_COLUMNS = (
        CRSSigmaProfileColumn(column="profile", compound_source="profil", mixture_total_source="profiltot"),
        CRSSigmaProfileColumn(
            column="hbprofile",
            compound_source="hbprofil",
            mixture_total_source="hbprofiltot",
            hb_channel=0,
        ),
        CRSSigmaProfileColumn(
            column="ohprofile",
            compound_source="hbprofil",
            mixture_total_source="hbprofiltot",
            hb_channel=1,
        ),
        CRSSigmaProfileColumn(
            column="otprofile",
            compound_source="hbprofil",
            mixture_total_source="hbprofiltot",
            hb_channel=2,
        ),
    )
    _SIGMA_PROFILE_EXTRA_COLUMNS = (
        CRSSigmaProfileColumn(column="orthprofile", compound_source="chdorth"),
    )

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
        mixture, and LLE tables. Sigma profile sections return one row per
        ``(compound, sigma)`` pair. Sigma potential sections are intentionally
        handled by :meth:`get_sigma_potential` instead. Quantity selection always
        uses raw CRS result keys; ``column_labels`` only controls the returned
        dataframe column names.
        """
        pd = self._import_pandas(self.__class__.__name__ + ".get_result_table")

        column_labels = self._normalize_column_label_mode(column_labels)
        results = self.get_results(section)
        property_name = self._result_property_name(results, section)
        if property_name in self._SIGMA_PROFILE_PROPERTIES:
            if split:
                raise ValueError("split=True is not supported for sigma profile result tables.")
            table = self._build_sigma_profile_result_table(pd, results, property_name, quantities)
            return self._rename_column_labels(table, column_labels)
        if property_name in self._RESULT_TABLE_UNSUPPORTED_PROPERTIES:
            raise NotImplementedError(
                "{} results are not supported by get_result_table(); use get_sigma_profile() "
                "or get_sigma_potential() instead.".format(property_name)
            )

        component_quantities, mixture_quantities, lle_quantities = self._select_result_table_quantities(
            results, property_name, quantities
        )
        component_names = self._resolve_component_names(results, int(results["ncomp"]))
        component = self._build_result_component_table(
            pd, results, property_name, component_quantities, component_names
        )
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
        combined = component.merge(mixture, on=["property", "method", "mixture"], how="left")
        return self._rename_column_labels(combined, column_labels)

    def get_result_table_metadata(
        self,
        section: Optional[str] = None,
        quantities: Union[str, Sequence[str]] = "default",
        split: bool = False,
    ) -> Union["pd.DataFrame", CRSResultTables]:
        """Return metadata describing CRS result table quantity columns.

        The returned dataframe describes raw CRS quantity keys, their table type,
        symbol, descriptive name, unit, and optional comment. Quantity selection uses
        the same raw CRS result keys and defaults as :meth:`get_result_table`.
        """
        pd = self._import_pandas(self.__class__.__name__ + ".get_result_table_metadata")

        results = self.get_results(section)
        property_name = self._result_property_name(results, section)
        if property_name in self._SIGMA_PROFILE_PROPERTIES:
            if split:
                raise ValueError("split=True is not supported for sigma profile result table metadata.")
            quantities = tuple(column.column for column in self._select_sigma_profile_columns(results, quantities))
            return self._build_quantity_metadata_table(pd, "sigma_profile", ("frac1", "sigma") + quantities)
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
    def combine_result_tables(
        tables: Sequence[Union["pd.DataFrame", CRSResultTables]],
        *,
        table: Optional[Literal["component", "mixture", "lle", "all"]] = None,
        source_column: Optional[str] = "source",
        reindex_mixture: bool = False,
        require_same_property: bool = True,
    ) -> Union["pd.DataFrame", CRSResultTables]:
        """Combine result tables returned by :meth:`get_result_table`.

        Plain dataframes from ``split=False`` are accepted with ``table=None``.
        Split result containers from ``split=True`` require ``table`` to select
        one subtable, or ``table="all"`` to combine all subtables independently.
        Dataframe schemas must match exactly; this method does not add missing
        quantity columns.
        """
        pd = CRSResults._import_pandas(CRSResults.__name__ + ".combine_result_tables")
        if table not in {None, "component", "mixture", "lle", "all"}:
            raise ValueError("table must be one of: None, component, mixture, lle, all")
        if not tables:
            raise ValueError("No result tables were provided.")

        are_split_tables = [isinstance(item, CRSResultTables) for item in tables]
        are_dataframes = [isinstance(item, pd.DataFrame) for item in tables]
        if not all(split or dataframe for split, dataframe in zip(are_split_tables, are_dataframes)):
            raise TypeError("tables must contain only pandas DataFrames or CRSResultTables.")
        if any(are_split_tables) and any(are_dataframes):
            raise TypeError("Cannot combine pandas DataFrames and CRSResultTables in one call.")

        if all(are_dataframes):
            if table is not None:
                raise ValueError("table must be None when combining pandas DataFrames.")
            frames = [(index, cast("pd.DataFrame", item)) for index, item in enumerate(tables)]
            return CRSResults._combine_result_table_frames(
                frames,
                pd=pd,
                source_column=source_column,
                reindex_mixture=reindex_mixture,
                require_same_property=require_same_property,
                allow_empty=False,
            )

        if table is None:
            raise ValueError("table must be specified when combining CRSResultTables.")

        split_tables = [cast(CRSResultTables, item) for item in tables]
        if table == "all":
            component = CRSResults._combine_result_table_frames(
                [(index, item.component) for index, item in enumerate(split_tables)],
                pd=pd,
                source_column=source_column,
                reindex_mixture=reindex_mixture,
                require_same_property=require_same_property,
                allow_empty=True,
            )
            mixture = CRSResults._combine_result_table_frames(
                [(index, item.mixture) for index, item in enumerate(split_tables)],
                pd=pd,
                source_column=source_column,
                reindex_mixture=reindex_mixture,
                require_same_property=require_same_property,
                allow_empty=True,
            )
            lle = CRSResults._combine_result_table_frames(
                [(index, item.lle) for index, item in enumerate(split_tables)],
                pd=pd,
                source_column=source_column,
                reindex_mixture=reindex_mixture,
                require_same_property=require_same_property,
                allow_empty=True,
            )
            return CRSResultTables(component=component, mixture=mixture, lle=lle)  # type: ignore[arg-type]

        frames = [(index, getattr(item, table)) for index, item in enumerate(split_tables)]
        combined = CRSResults._combine_result_table_frames(
            frames,
            pd=pd,
            source_column=source_column,
            reindex_mixture=reindex_mixture,
            require_same_property=require_same_property,
            allow_empty=False,
        )
        assert combined is not None
        return combined

    @staticmethod
    def plot_sigma_profile_table(
        table: Any,
        *,
        y: Optional[Union[str, Sequence[str]]] = None,
        split: bool = False,
        ax: Optional[Any] = None,
        plot_fig: bool = True,
    ) -> "Figure":
        """Plot a sigma profile table returned by :meth:`get_result_table`.

        The table may use raw, symbol, name, or name-unit column labels. By default,
        the main sigma profile columns are plotted when present. Pass *y* to select
        explicit columns, using either the current dataframe labels or raw quantity
        keys such as ``"profile"`` and ``"hbprofile"``. With ``split=True``,
        each selected profile column is drawn on a separate subplot.
        """
        plt = CRSResults._import_matplotlib_pyplot("CRSResults.plot_sigma_profile_table")

        sigma_column = CRSResults._resolve_sigma_profile_table_column(table, "sigma")
        name_column = CRSResults._resolve_sigma_profile_table_column(table, "name")
        y_columns = CRSResults._resolve_sigma_profile_y_columns(table, y)
        if not y_columns:
            raise ValueError("Sigma profile table does not contain any plottable profile columns.")

        if split:
            fig, axes = CRSResults._sigma_profile_subplot_axes(plt, ax, len(y_columns))
        else:
            if ax is None:
                fig, ax = plt.subplots()
            else:
                fig = ax.figure
            axes = (ax,)

        groups = [(name, group.sort_values(sigma_column)) for name, group in table.groupby(name_column, sort=False)]
        color_by_name, linestyle_by_column = CRSResults._sigma_profile_plot_styles(plt, groups, y_columns)
        for axis_index, axis in enumerate(axes):
            axis_y_columns = (y_columns[axis_index],) if split else y_columns
            for name, group in groups:
                is_total = str(name) == "Total"
                for y_column in axis_y_columns:
                    axis.plot(
                        group[sigma_column],
                        group[y_column],
                        color=color_by_name[name],
                        linestyle=linestyle_by_column[y_column],
                        linewidth=2.0 if is_total else 1.2,
                        alpha=0.95 if is_total else 0.85,
                        label="{}: {}".format(name, y_column) if not split else str(name),
                    )
            axis.set_ylabel(axis_y_columns[0] if split else "sigma profile")
            axis.legend()

        axes[-1].set_xlabel(sigma_column)

        if plot_fig:
            plt.show()
        return fig

    @staticmethod
    def _sigma_profile_plot_styles(
        plt: Any,
        groups: Sequence[Tuple[Any, Any]],
        y_columns: Sequence[str],
    ) -> Tuple[Dict[Any, Any], Dict[str, str]]:
        color_cycle = cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        color_by_name = {name: next(color_cycle) for name, _ in groups}
        linestyles = ("-", "--", ":", "-.")
        linestyle_by_column = {
            column: linestyles[index % len(linestyles)]
            for index, column in enumerate(y_columns)
        }
        return color_by_name, linestyle_by_column

    @staticmethod
    def _sigma_profile_subplot_axes(plt: Any, ax: Optional[Any], nplots: int) -> Tuple[Any, Tuple[Any, ...]]:
        if ax is None:
            fig, axes = plt.subplots(nplots, 1, sharex=True, squeeze=False)
            return fig, tuple(axes[:, 0])

        if hasattr(ax, "ravel") and not hasattr(ax, "plot"):
            axes = tuple(ax.ravel())
        else:
            try:
                axes = tuple(ax)
            except TypeError:
                axes = (ax,)

        if len(axes) != nplots:
            raise ValueError("split=True requires {} axes; got {}".format(nplots, len(axes)))

        return axes[0].figure, axes

    @staticmethod
    def _default_sigma_profile_y_keys(table: Any) -> Tuple[str, ...]:
        has_oh = CRSResults._resolve_sigma_profile_table_column(
            table, "ohprofile", required=False
        ) is not None
        has_ot = CRSResults._resolve_sigma_profile_table_column(
            table, "otprofile", required=False
        ) is not None

        if has_oh or has_ot:
            return ("profile", "ohprofile", "otprofile")
        return ("profile", "hbprofile")

    @staticmethod
    def _resolve_sigma_profile_y_columns(
        table: Any,
        y: Optional[Union[str, Sequence[str]]],
    ) -> Tuple[str, ...]:
        requested = (
            CRSResults._default_sigma_profile_y_keys(table)
            if y is None
            else (y,) if isinstance(y, str)
            else tuple(y)
        )
        required = y is not None

        columns = tuple(
            CRSResults._resolve_sigma_profile_table_column(table, key, required=required)
            for key in requested
        )
        return tuple(column for column in columns if column is not None)

    @staticmethod
    def _resolve_sigma_profile_table_column(table: Any, key: str, *, required: bool = True) -> Optional[str]:
        for column in CRSResults._sigma_profile_table_column_candidates(key):
            if column in table.columns:
                return column

        if required:
            raise ValueError(
                "Sigma profile table is missing a {} column. Accepted labels: {}".format(
                    key, ", ".join(CRSResults._sigma_profile_table_column_candidates(key))
                )
            )
        return None

    @staticmethod
    def _sigma_profile_table_column_candidates(key: str) -> Tuple[str, ...]:
        metadata = CRSResults._RESULT_TABLE_QUANTITY_METADATA.get(key, {})
        candidates = [key]
        symbol = metadata.get("symbol")
        name = metadata.get("name")
        unit = metadata.get("unit")
        if symbol:
            candidates.append(symbol)
        if name:
            candidates.append(name)
            if unit:
                candidates.append("{} [{}]".format(name, unit))
        return tuple(dict.fromkeys(candidates))


    @staticmethod
    def plot_lle_phase_diagram(
        table: Any,
        *,
        experimental_phase_I: Optional[Any] = None,
        experimental_phase_II: Optional[Any] = None,
        experimental_temperature: Optional[Any] = None,
        temperature: Optional[float] = None,
        experimental_label: str = "exp",
        component_labels: Optional[Sequence[str]] = None,
        component_order: Optional[Sequence[int]] = None,
        plot_tielines: bool = True,
        plot_phase_boundaries: bool = False,
        plot_feed: bool = False,
        ax: Optional[Any] = None,
        plot_fig: bool = True,
    ) -> "Figure":
        """Plot T-x binary or isothermal ternary CRS phase diagrams.

        ``table`` must be a pandas.DataFrame prepared from one or more
        :meth:`get_result_table` outputs for one ``LLE``, ``STABILITY``,
        ``BINMIXCOEF``, or ``TERNARYMIX`` calculation.

        When present, ``converged`` and ``llle_detected`` are used to keep only
        converged, non-LLLE rows.

        Ternary systems require a single temperature slice; pass ``temperature`` when
        multiple temperatures are present.

        ``LLE``, ``BINMIXCOEF``, and ``TERNARYMIX`` tables draw phase-I and phase-II
        points. ``LLE`` may also draw feed compositions. ``STABILITY`` tables draw only
        stable and unstable feed compositions.

        Parameters
        ----------
        table
            CRS LLE result table.
        experimental_phase_I, experimental_phase_II
            Optional experimental phase compositions with shape ``(n_tie_lines, n_components)``.
        experimental_temperature
            Optional experimental temperatures for binary T-x plots.
        temperature
            Temperature to plot for ternary diagrams with multiple temperatures.
        experimental_label
            Legend label for experimental data.
        component_labels
            Axis labels.
        component_order
            Axis order by ``cid``. Experimental arrays are reordered from sorted
            ``cid`` order to this order.
        plot_tielines
            Whether to draw phase-I/phase-II tie-lines.
        plot_phase_boundaries
            Whether to draw phase-boundary curves.
        plot_feed
            Whether to draw feed compositions for ``LLE`` tables.
        ax
            Matplotlib axes to draw into. Ternary plots use ``python-ternary``.
        plot_fig
            Whether to show the figure.

        Returns
        -------
        matplotlib.figure.Figure
            The created or updated figure.
        """
        plt = CRSResults._import_matplotlib_pyplot("CRSResults.plot_lle_phase_diagram")
        pd = CRSResults._import_pandas("CRSResults.plot_lle_phase_diagram")

        if not isinstance(table, pd.DataFrame):
            raise TypeError("table must be a pandas DataFrame")

        method_color_cycle = cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        plot_style = {
            "marker_size": 36,
            "phase_marker": "o",
            "experimental_color": "#333333",

            "feed_marker": "^",
            "feed_label": "feed",
            "feed_color": "#B8860B",

            "stable_marker": "^",
            "stable_label": "feed (stable)",
            "stable_color": "#B8860B",
            "unstable_marker": "x",
            "unstable_label": "feed (unstable)",
            "unstable_color": "#B8860B",

            "tieline_linewidth": 1.0,
            "tieline_alpha": 0.30,
            "phase_boundary_linewidth": 1.5,
            "phase_boundary_alpha": 0.80,

            "grid_color": "#9A9A9A",
            "grid_linewidth": 0.5,
            "grid_alpha": 0.4,

            "axis_color": "#555555",
            "axis_boundary_linewidth": 1.2,
            "axis_label_fontsize": 14,
            "axis_labelpad": 12,
            "axis_label_offsets": {"bottom": 0.00, "right": 0.10, "left": 0.10},

            "tick_linewidth": 0.8,
            "tick_fontsize": 10,
            "tick_offset": 0.015,

            "binary_figsize": (7, 6),
            "ternary_figsize": (15, 13),
            "ternary_legend_loc": "upper right",
            "title_fontsize": 14,
            "legend_fontsize": 12,
        }

        # Find required result columns and validate table-level assumptions.
        def existing_column(candidates: Sequence[str]) -> Optional[str]:
            for candidate in candidates:
                if candidate in table.columns:
                    return candidate
            return None

        def required_column(candidates: Sequence[str]) -> str:
            column = existing_column(candidates)
            if column is None:
                raise ValueError(
                    "plot_lle_phase_diagram requires an LLE result table with a {} column. "
                    "Accepted column names: {}".format(candidates[-1], ", ".join(candidates))
                )
            return column

        def require_one_supported_property_name(property_column: str) -> str:
            property_names = tuple(str(value).upper() for value in table[property_column].dropna().unique())
            supported_properties = ("LLE", "STABILITY", "BINMIXCOEF", "TERNARYMIX")
            if len(property_names) != 1 or property_names[0] not in supported_properties:
                raise ValueError(
                    "plot_lle_phase_diagram only accepts one of {} result tables. Got {}.".format(
                        ", ".join(supported_properties), property_names or "no property"
                    )
                )
            return property_names[0]

        def require_single_system_table() -> None:
            names_per_cid = table[["cid", "name"]].dropna().drop_duplicates()
            grouped = names_per_cid.groupby("cid")["name"].nunique()

            if (grouped > 1).any():
                raise ValueError(
                    "plot_lle_phase_diagram only accepts results for one chemical system."
                )

        # Resolve all input columns before the plotting table is copied or filtered.
        required_column(("method",))
        required_column(("cid",))
        required_column(("name",))
        property_column = required_column(("property",))
        property_name = require_one_supported_property_name(property_column)
        require_single_system_table()

        xI_column = required_column(("xI", "x_I", "Phase I mole fraction"))
        xII_column = required_column(("xII", "x_II", "Phase II mole fraction"))
        temperature_column = required_column(("temperature", "T", "Temperature [K]", "Temperature"))
        converged_column = existing_column(("converged", "LLE flash solver converged"))
        llle_detected_column = existing_column(("llle_detected", "LLLE detected"))
        feed_column = None
        unstable_column = None
        if property_name in {"LLE", "STABILITY"}:
            feed_column = required_column(("frac1", "x", "Feed mole fraction"))
            unstable_column = required_column(("unstable", "Feed unstable"))

        # Build the filtered plot table and normalize component axis metadata.
        def filtered_plot_table(
            plot_table: Any,
            converged_column: Optional[str],
            llle_detected_column: Optional[str],
        ) -> Any:
            if converged_column is not None:
                plot_table = plot_table[plot_table[converged_column].astype(bool)].copy()

            if llle_detected_column is not None:
                plot_table = plot_table[~plot_table[llle_detected_column].astype(bool)].copy()

            if plot_table.empty:
                raise ValueError("No rows remain after filtering converged=True and llle_detected=False")

            return plot_table

        def normalized_component_plot_metadata(
            plot_table: Any,
            component_order: Optional[Sequence[int]],
            component_labels: Optional[Sequence[str]],
        ) -> Tuple[Tuple[int, ...], int, List[int], Tuple[str, ...]]:
            base_order = tuple(int(cid) for cid in sorted(plot_table["cid"].dropna().unique()))
            if component_order is None:
                normalized_order = base_order
            else:
                normalized_order = tuple(int(cid) for cid in component_order)

            if set(normalized_order) != set(base_order):
                raise ValueError(
                    "component_order must contain the same cid values as the table. "
                    "Expected {}, got {}.".format(base_order, normalized_order)
                )

            ncomp = len(normalized_order)
            if ncomp not in (2, 3):
                raise NotImplementedError("Only binary and ternary phase diagrams are supported")

            reorder_indices = [base_order.index(cid) for cid in normalized_order]

            if component_labels is None:
                name_table = plot_table[["cid", "name"]].drop_duplicates()
                labels = []
                for cid in normalized_order:
                    name = str(name_table[name_table["cid"] == cid]["name"].iloc[0])
                    labels.append("x{} ({})".format(len(labels) + 1, name))
                normalized_labels = tuple(labels)
            else:
                normalized_labels = tuple(component_labels)

            if len(normalized_labels) != ncomp:
                raise ValueError("component_labels must contain {} labels".format(ncomp))

            return normalized_order, ncomp, reorder_indices, normalized_labels

        # Filter plot rows and derive the component order used by all arrays.
        table = filtered_plot_table(table.copy(), converged_column, llle_detected_column)

        component_order, ncomp, reorder_indices, component_labels = normalized_component_plot_metadata(
            table,
            component_order,
            component_labels,
        )

        group_column = existing_column(("tie_line", "source"))
        if group_column is None:
            table["_lle_group"] = 0
            group_column = "_lle_group"

        # Convert grouped component rows into phase/feed arrays.
        def get_tie_line_groups(method_table: Any) -> List[Any]:
            groups = []
            expected_cids = set(component_order)

            for _, group in method_table.groupby(group_column, sort=True):
                values = group.set_index("cid")
                group_cids = set(int(cid) for cid in values.index)

                if group_cids != expected_cids or len(group_cids) != len(values.index):
                    raise ValueError(
                        "plot_lle_phase_diagram only accepts complete tie-lines for one binary or ternary system."
                    )

                groups.append(values)

            if not groups:
                raise ValueError("No LLE tie-lines found.")

            return groups

        def phase_array(groups: Sequence[Any], column: str) -> "np.ndarray":
            rows = [[float(values.loc[cid, column]) for cid in component_order] for values in groups]
            arr = np.asarray(rows, dtype=float)
            totals = arr.sum(axis=1, keepdims=True)
            totals[totals == 0.0] = 1.0
            return arr / totals

        def group_array(groups: Sequence[Any], column: str) -> "np.ndarray":
            values = []
            for group in groups:
                group_values = group[column].dropna()
                if group_values.empty:
                    raise ValueError("LLE group is missing a non-empty {} value".format(column))
                values.append(float(group_values.iloc[0]))
            return np.asarray(values, dtype=float)

        def calculated_series(method_table: Any) -> Tuple[Sequence[Any], "np.ndarray", "np.ndarray", "np.ndarray"]:
            groups = get_tie_line_groups(method_table)
            return (
                groups,
                phase_array(groups, xI_column),
                phase_array(groups, xII_column),
                group_array(groups, temperature_column),
            )

        def iter_method_tables(plot_table: Any) -> List[Tuple[str, Any]]:
            method_tables = [
                (str(method), method_table.copy())
                for method, method_table in plot_table.groupby("method", sort=True)
            ]
            if not method_tables:
                raise ValueError("plot_lle_phase_diagram requires at least one non-empty method group")
            return method_tables

        def first_feed_series(
            method_tables: Sequence[Tuple[str, Any]]
        ) -> Tuple["np.ndarray", "np.ndarray", "np.ndarray"]:
            assert feed_column is not None
            assert unstable_column is not None

            _, method_table = method_tables[0]
            groups = get_tie_line_groups(method_table)
            return (
                phase_array(groups, feed_column),
                group_array(groups, temperature_column),
                group_array(groups, unstable_column).astype(bool),
            )

        def require_single_method(method_tables: Sequence[Tuple[str, Any]]) -> Tuple[str, Any]:
            if len(method_tables) != 1:
                raise ValueError("STABILITY phase diagrams require a single method")
            return method_tables[0]

        # Normalize experimental arrays to the plot component order.
        def experimental_array(values: Optional[Any], name: str) -> Optional["np.ndarray"]:
            if values is None:
                return None
            arr = np.asarray(values, dtype=float)
            if arr.ndim != 2 or arr.shape[1] != ncomp:
                raise ValueError("{} must have shape (n_tie_lines, ncomp={})".format(name, ncomp))
            return arr[:, reorder_indices]

        # Convert optional experimental inputs once so plot helpers can assume arrays.
        exp_I = experimental_array(experimental_phase_I, "experimental_phase_I")
        exp_II = experimental_array(experimental_phase_II, "experimental_phase_II")
        exp_temperature = (
            None if experimental_temperature is None else np.asarray(experimental_temperature, dtype=float).ravel()
        )
        if exp_temperature is not None:
            for exp in (exp_I, exp_II):
                if exp is not None and len(exp_temperature) != len(exp):
                    raise ValueError("experimental_temperature must match the number of experimental tie-lines")

        # Define binary plot primitives and layer helpers.
        def plot_binary_phase_boundary(ax: Any, x: "np.ndarray", temperatures: "np.ndarray", color: str) -> None:
            order = np.argsort(temperatures)
            x_sort = x[order]
            temperature_sort = temperatures[order]
            ax.plot(
                x_sort,
                temperature_sort,
                color=color,
                alpha=plot_style["phase_boundary_alpha"],
                linewidth=plot_style["phase_boundary_linewidth"],
            )

        def scatter_binary_points(
            ax: Any,
            x: "np.ndarray",
            temperatures: "np.ndarray",
            marker: str,
            color: str,
            label: str,
            facecolors: Optional[str] = "none",
        ) -> None:
            scatter_kwargs = {
                "marker": marker,
                "color": color,
                "s": plot_style["marker_size"],
                "label": label,
            }
            if facecolors is not None:
                scatter_kwargs["facecolors"] = facecolors
            ax.scatter(x, temperatures, **scatter_kwargs)

        def plot_binary_calculated(
            ax: Any,
            phase_I: "np.ndarray",
            phase_II: "np.ndarray",
            temperatures: "np.ndarray",
            label: str,
            color: str,
        ) -> None:
            marker = plot_style["phase_marker"]

            scatter_binary_points(ax, phase_I[:, 0], temperatures, marker, color, label)
            scatter_binary_points(ax, phase_II[:, 0], temperatures, marker, color, "_nolegend_")

            if plot_phase_boundaries:
                plot_binary_phase_boundary(ax, phase_I[:, 0], temperatures, color)
                plot_binary_phase_boundary(ax, phase_II[:, 0], temperatures, color)

        def plot_binary_feed(ax: Any, feed: "np.ndarray", temperatures: "np.ndarray") -> None:
            marker = plot_style["feed_marker"]
            color = plot_style["feed_color"]
            label = plot_style["feed_label"]
            scatter_binary_points(ax, feed[:, 0], temperatures, marker, color, label)

        def plot_binary_feed_stability(
            ax: Any, feed: "np.ndarray", temperatures: "np.ndarray", unstable: "np.ndarray"
        ) -> None:
            stable = ~unstable
            marker = plot_style["unstable_marker"]
            color = plot_style["unstable_color"]
            label = plot_style["unstable_label"]
            if unstable.any():
                scatter_binary_points(
                    ax, feed[unstable, 0], temperatures[unstable], marker, color, label, facecolors=None
                )

            marker = plot_style["stable_marker"]
            color = plot_style["stable_color"]
            label = plot_style["stable_label"]
            if stable.any():
                scatter_binary_points(ax, feed[stable, 0], temperatures[stable], marker, color, label)

        def plot_binary_experimental(ax: Any, reference_temperature: "np.ndarray") -> None:
            marker = plot_style["phase_marker"]
            color = plot_style["experimental_color"]
            label = experimental_label
            if exp_I is not None:
                exp_T = exp_temperature if exp_temperature is not None else reference_temperature
                if len(exp_T) != len(exp_I):
                    raise ValueError(
                        "experimental_temperature is required when experimental_phase_I length "
                        "does not match the calculated temperature grid"
                    )
                scatter_binary_points(ax, exp_I[:, 0], exp_T, marker, color, label)

                if plot_phase_boundaries:
                    plot_binary_phase_boundary(ax, exp_I[:, 0], exp_T, color)

            if exp_II is not None:
                exp_T = exp_temperature if exp_temperature is not None else reference_temperature
                if len(exp_T) != len(exp_II):
                    raise ValueError(
                        "experimental_temperature is required when experimental_phase_II length "
                        "does not match the calculated temperature grid"
                    )
                if exp_I is None:
                    scatter_binary_points(ax, exp_II[:, 0], exp_T, marker, color, label)
                else:
                    scatter_binary_points(ax, exp_II[:, 0], exp_T, marker, color, "_nolegend_")

                if plot_phase_boundaries:
                    plot_binary_phase_boundary(ax, exp_II[:, 0], exp_T, color)

        def finish_figure(fig: "Figure") -> "Figure":
            fig.tight_layout()
            if plot_fig:
                plt.show()
            return fig

        # Draw binary T-x diagrams directly on a matplotlib axes.
        if ncomp == 2:
            if ax is None:
                fig, ax = plt.subplots(figsize=plot_style["binary_figsize"])
            else:
                fig = ax.figure

            method_tables = iter_method_tables(table)

            def plot_binary_stability(method_tables: Sequence[Tuple[str, Any]]) -> None:
                require_single_method(method_tables)
                feed, temperatures, unstable = first_feed_series(method_tables)
                plot_binary_feed_stability(ax, feed, temperatures, unstable)

            def plot_binary_lle(method_tables: Sequence[Tuple[str, Any]]) -> "np.ndarray":
                reference_temperature = None
                for label, method_table in method_tables:
                    color = next(method_color_cycle)
                    _, phase_I, phase_II, temperatures = calculated_series(method_table)
                    if reference_temperature is None:
                        reference_temperature = temperatures
                    plot_binary_calculated(ax, phase_I, phase_II, temperatures, label, color)
                assert reference_temperature is not None
                return reference_temperature

            def plot_binary_lle_feed(method_tables: Sequence[Tuple[str, Any]]) -> None:
                feed, feed_temperatures, unstable = first_feed_series(method_tables)
                if len(method_tables) == 1:
                    plot_binary_feed_stability(ax, feed, feed_temperatures, unstable)
                else:
                    plot_binary_feed(ax, feed, feed_temperatures)

            # Dispatch binary plot layers by property type.
            if property_name == "STABILITY":
                plot_binary_stability(method_tables)
            else:
                reference_temperature = plot_binary_lle(method_tables)

                if property_name == "LLE" and plot_feed:
                    plot_binary_lle_feed(method_tables)

                plot_binary_experimental(ax, reference_temperature)

            ax.set_xlabel(
                "{}".format(component_labels[0]),
                fontsize=plot_style["axis_label_fontsize"],
                labelpad=plot_style["axis_labelpad"],
            )
            ax.set_ylabel(
                "Temperature (K)",
                fontsize=plot_style["axis_label_fontsize"],
                labelpad=plot_style["axis_labelpad"],
            )
            ax.set_xlim(0.0, 1.0)
            ax.set_title("Binary {} phase diagram".format(property_name), fontsize=plot_style["title_fontsize"])
            ax.legend(fontsize=plot_style["legend_fontsize"])

            return finish_figure(fig)

        # Select one temperature slice before building a ternary plot.
        requested_temperature = temperature
        plot_table = table
        row_temperatures = np.asarray(plot_table[temperature_column], dtype=float)
        unique_temperatures = np.unique(np.round(row_temperatures, decimals=8))

        if requested_temperature is None:
            if len(unique_temperatures) != 1:
                raise ValueError(
                    "Ternary LLE plots require a single temperature. "
                    "Pass temperature=... to select one."
                )
            selected_temperature = float(unique_temperatures[0])
        else:
            selected_temperature = float(requested_temperature)
            mask = np.isclose(row_temperatures, selected_temperature)
            if not mask.any():
                raise ValueError("No ternary LLE data found at temperature {:g} K".format(selected_temperature))
            plot_table = plot_table.loc[mask].copy()

        try:
            import ternary
        except ImportError:
            raise ImportError(
                "CRSResults.plot_lle_phase_diagram: ternary LLE plots require the 'python-ternary' package"
            )

        # Define ternary plot primitives and layer helpers.
        def ternary_points(compositions: Any) -> List[Tuple[float, ...]]:
            comp = np.asarray(compositions, dtype=float)
            totals = comp.sum(axis=1, keepdims=True)
            totals[totals == 0.0] = 1.0
            return [tuple(row) for row in comp / totals]

        def scatter_ternary_points(
            tax: Any,
            compositions: "np.ndarray",
            marker: str,
            color: str,
            label: str,
            facecolors: Optional[str] = "none",
        ) -> None:
            scatter_kwargs = {
                "marker": marker,
                "color": color,
                "s": plot_style["marker_size"],
                "label": label,
            }
            if facecolors is not None:
                scatter_kwargs["facecolors"] = facecolors
            tax.scatter(ternary_points(compositions), **scatter_kwargs)

        def ordered_boundary_points(compositions: Any) -> List[Tuple[float, ...]]:
            points = ternary_points(compositions)

            if len(points) <= 2:
                return points

            comp = np.asarray(points, dtype=float)
            remaining = list(range(len(comp)))
            start = max(remaining, key=lambda i: comp[i, 0])
            ordered = [start]
            remaining.remove(start)

            while remaining:
                current = comp[ordered[-1]]
                distances = np.linalg.norm(comp[remaining] - current, axis=1)
                next_pos = int(np.argmin(distances))
                ordered.append(remaining.pop(next_pos))

            return [points[i] for i in ordered]

        def plot_ternary_phase_boundary(tax: Any, compositions: "np.ndarray", color: str) -> None:
            tax.plot(
                ordered_boundary_points(compositions),
                color=color,
                linewidth=plot_style["phase_boundary_linewidth"],
                alpha=plot_style["phase_boundary_alpha"],
            )

        def plot_ternary_tielines(
            tax: Any,
            phase_I: "np.ndarray",
            phase_II: "np.ndarray",
            color: str,
        ) -> None:
            for p1, p2 in zip(ternary_points(phase_I), ternary_points(phase_II)):
                tax.line(
                    p1,
                    p2,
                    color=color,
                    linewidth=plot_style["tieline_linewidth"],
                    alpha=plot_style["tieline_alpha"],
                )

        def plot_ternary_calculated(
            tax: Any, phase_I: "np.ndarray", phase_II: "np.ndarray", label: str, color: str
        ) -> None:
            marker = plot_style["phase_marker"]

            scatter_ternary_points(tax, phase_I, marker, color, label)
            scatter_ternary_points(tax, phase_II, marker, color, "_nolegend_")

            if plot_phase_boundaries:
                plot_ternary_phase_boundary(tax, phase_I, color)
                plot_ternary_phase_boundary(tax, phase_II, color)

            if plot_tielines:
                plot_ternary_tielines(tax, phase_I, phase_II, color)

        def plot_ternary_feed(tax: Any, feed: "np.ndarray") -> None:
            marker = plot_style["feed_marker"]
            color = plot_style["feed_color"]
            label = plot_style["feed_label"]
            scatter_ternary_points(tax, feed, marker, color, label)

        def plot_ternary_feed_stability(tax: Any, feed: "np.ndarray", unstable: "np.ndarray") -> None:
            stable = ~unstable
            marker = plot_style["unstable_marker"]
            color = plot_style["unstable_color"]
            label = plot_style["unstable_label"]
            if unstable.any():
                scatter_ternary_points(tax, feed[unstable], marker, color, label, facecolors=None)

            marker = plot_style["stable_marker"]
            color = plot_style["stable_color"]
            label = plot_style["stable_label"]
            if stable.any():
                scatter_ternary_points(tax, feed[stable], marker, color, label)

        def plot_ternary_experimental(tax: Any) -> None:
            if exp_I is not None:
                marker = plot_style["phase_marker"]
                color = plot_style["experimental_color"]
                label = experimental_label
                scatter_ternary_points(tax, exp_I, marker, color, label)
                if plot_phase_boundaries:
                    plot_ternary_phase_boundary(tax, exp_I, color)

            if exp_II is not None:
                marker = plot_style["phase_marker"]
                color = plot_style["experimental_color"]
                label = experimental_label if exp_I is None else "_nolegend_"
                scatter_ternary_points(tax, exp_II, marker, color, label)
                if plot_phase_boundaries:
                    plot_ternary_phase_boundary(tax, exp_II, color)

            if plot_tielines and exp_I is not None and exp_II is not None:
                plot_ternary_tielines(tax, exp_I, exp_II, plot_style["experimental_color"])

        # Draw ternary isothermal diagrams through python-ternary.
        fig, tax = ternary.figure(ax=ax, scale=1.0)
        fig.set_size_inches(*plot_style["ternary_figsize"])
        tax.gridlines(
            color=plot_style["grid_color"],
            multiple=0.1,
            linewidth=plot_style["grid_linewidth"],
            alpha=plot_style["grid_alpha"],
        )
        tax.boundary(
            linewidth=plot_style["axis_boundary_linewidth"],
            linestyle="-",
            zorder=10,
            axes_colors={axis: plot_style["axis_color"] for axis in ("l", "r", "b")},
        )

        method_tables = iter_method_tables(plot_table)

        def plot_ternary_lle(method_tables: Sequence[Tuple[str, Any]]) -> None:
            for label, method_table in method_tables:
                color = next(method_color_cycle)
                _, phase_I, phase_II, _ = calculated_series(method_table)
                plot_ternary_calculated(tax, phase_I, phase_II, label, color)

        def plot_ternary_lle_feed(method_tables: Sequence[Tuple[str, Any]]) -> None:
            feed, _, unstable = first_feed_series(method_tables)
            if len(method_tables) == 1:
                plot_ternary_feed_stability(tax, feed, unstable)
            else:
                plot_ternary_feed(tax, feed)

        def plot_ternary_stability(method_tables: Sequence[Tuple[str, Any]]) -> None:
            require_single_method(method_tables)
            feed, _, unstable = first_feed_series(method_tables)
            plot_ternary_feed_stability(tax, feed, unstable)

        # Dispatch ternary plot layers by property type.
        if property_name == "STABILITY":
            plot_ternary_stability(method_tables)
        else:
            plot_ternary_lle(method_tables)

            if property_name == "LLE" and plot_feed:
                plot_ternary_lle_feed(method_tables)

            plot_ternary_experimental(tax)

        tax.bottom_axis_label(
            component_labels[0],
            fontsize=plot_style["axis_label_fontsize"],
            offset=plot_style["axis_label_offsets"]["bottom"],
            color=plot_style["axis_color"],
        )
        tax.right_axis_label(
            component_labels[1],
            fontsize=plot_style["axis_label_fontsize"],
            offset=plot_style["axis_label_offsets"]["right"],
            color=plot_style["axis_color"],
        )
        tax.left_axis_label(
            component_labels[2],
            fontsize=plot_style["axis_label_fontsize"],
            offset=plot_style["axis_label_offsets"]["left"],
            color=plot_style["axis_color"],
        )
        tax.ticks(
            axis="lbr",
            multiple=0.1,
            linewidth=plot_style["tick_linewidth"],
            axes_colors={axis: plot_style["axis_color"] for axis in ("l", "r", "b")},
            tick_formats="%.1f",
            fontsize=plot_style["tick_fontsize"],
            offset=plot_style["tick_offset"],
        )

        ax_obj = tax.get_axes()
        ax_obj.set_aspect("equal")
        tax.clear_matplotlib_ticks()

        for spine in ax_obj.spines.values():
            spine.set_visible(False)

        ax_obj.legend(loc=plot_style["ternary_legend_loc"], fontsize=plot_style["legend_fontsize"])
        tax.set_title(
            "Ternary {} phase diagram at {:g} K".format(property_name, selected_temperature),
            fontsize=plot_style["title_fontsize"],
        )
        tax._redraw_labels()

        return finish_figure(fig)

    @staticmethod
    def _combine_result_table_frames(
        frames: Sequence[Tuple[int, Optional["pd.DataFrame"]]],
        *,
        pd: Any,
        source_column: Optional[str],
        reindex_mixture: bool,
        require_same_property: bool,
        allow_empty: bool,
    ) -> Optional["pd.DataFrame"]:
        usable_frames = [(source_index, frame) for source_index, frame in frames if frame is not None]
        if not usable_frames:
            if allow_empty:
                return None
            raise ValueError("No usable result tables were provided.")

        expected_columns = list(usable_frames[0][1].columns)
        for _, frame in usable_frames[1:]:
            columns = list(frame.columns)
            if columns != expected_columns:
                raise ValueError(
                    "All result tables must have identical columns in identical order. "
                    "Expected {}, got {}.".format(expected_columns, columns)
                )

        if source_column is not None and source_column in expected_columns:
            raise ValueError("source_column {!r} already exists in the result table.".format(source_column))
        if reindex_mixture and "mixture" not in expected_columns:
            raise ValueError("reindex_mixture=True requires a 'mixture' column.")

        if require_same_property:
            if "property" not in expected_columns:
                raise ValueError("require_same_property=True requires a 'property' column.")
            properties = []
            for _, frame in usable_frames:
                if frame.empty:
                    continue
                unique_properties = frame["property"].dropna().unique()
                if len(unique_properties) != 1:
                    raise ValueError("Each result table must contain exactly one non-null property value.")
                properties.append(unique_properties[0])
            if properties and len(set(properties)) != 1:
                raise ValueError("All result tables must have the same property.")

        combined_frames = []
        mixture_offset = 0
        for source_index, frame in usable_frames:
            combined_frame = frame.copy()
            if reindex_mixture:
                valid_mixture = combined_frame["mixture"].notna()
                if valid_mixture.any():
                    combined_frame.loc[valid_mixture, "mixture"] = (
                        combined_frame.loc[valid_mixture, "mixture"] + mixture_offset
                    )
                    mixture_offset = int(combined_frame.loc[valid_mixture, "mixture"].max()) + 1
            if source_column is not None:
                combined_frame[source_column] = source_index
            combined_frames.append(combined_frame)

        return pd.concat(combined_frames, ignore_index=True)

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
        if column_labels == "raw" or key in {"property", "method", "mixture", "cid", "name", "molmass", "tie_line"}:
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
        columns = ["table", "quantity", "symbol", "name", "unit", "comment"]
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
                    "comment": metadata.get("comment") or "",
                }
            )
        return pd.DataFrame(rows, columns=columns)

    def _filter_component_quantities(self, property_name: str, quantities: Sequence[str]) -> Tuple[str, ...]:
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
            raise KeyError(
                "Unknown CRS result quantity/quantities: {}. Available quantities: {}".format(unknown, available)
            )

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
            if property_name in self._RESULT_TABLE_LLE_PROPERTIES and key in self._RESULT_TABLE_LLE_KNOWN_QUANTITIES:
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

    def _select_sigma_profile_columns(
        self,
        results: dict,
        quantities: Union[str, Sequence[str]],
    ) -> Tuple[CRSSigmaProfileColumn, ...]:
        available_columns = self._available_sigma_profile_columns(results)
        available_by_name = {column.column: column for column in available_columns}

        if quantities == "default":
            return tuple(column for column in self._SIGMA_PROFILE_DEFAULT_COLUMNS if column.column in available_by_name)

        if quantities == "all":
            return available_columns

        if isinstance(quantities, str):
            raise ValueError("quantities must be 'default', 'all', or a sequence of quantity names.")

        requested = tuple(quantities)
        unknown = sorted(set(requested) - set(available_by_name))
        if unknown:
            raise KeyError(
                "Unknown sigma profile quantity/quantities: {}. Available quantities: {}".format(
                    unknown, ", ".join(sorted(available_by_name))
                )
            )
        return tuple(available_by_name[quantity] for quantity in requested)

    def _available_sigma_profile_columns(self, results: dict) -> Tuple[CRSSigmaProfileColumn, ...]:
        columns = []
        nhb = int(results.get("nhb", 0))
        for column in self._SIGMA_PROFILE_DEFAULT_COLUMNS + self._SIGMA_PROFILE_EXTRA_COLUMNS:
            if column.compound_source not in results:
                continue
            if column.hb_channel is not None and column.hb_channel >= nhb:
                continue
            columns.append(column)
        return tuple(columns)

    def _build_sigma_profile_result_table(
        self,
        pd: Any,
        results: dict,
        property_name: str,
        quantities: Union[str, Sequence[str]],
    ) -> "pd.DataFrame":
        columns = self._select_sigma_profile_columns(results, quantities)
        ncomp = int(results["ncomp"])
        nitems = int(results["nitems"])
        nhb = int(results["nhb"])
        sigma = np.asarray(results["chdval"], dtype=float).reshape(nitems)
        frac1 = self._sigma_profile_frac1_values(results, ncomp)
        method = results["method"]
        component_names = self._resolve_component_names(results, ncomp)

        rows = []
        for cid in range(ncomp):
            profile_values = {
                column.column: self._sigma_profile_component_values(results, column, cid, ncomp, nitems, nhb)
                for column in columns
            }
            for index, sigma_value in enumerate(sigma):
                row: Dict[str, Any] = {
                    "property": property_name,
                    "method": method,
                    "mixture": 0,
                    "cid": cid,
                    "name": component_names[cid],
                    "frac1": frac1[cid],
                    "sigma": sigma_value,
                }
                for column in columns:
                    row[column.column] = profile_values[column.column][index]
                rows.append(row)

        if property_name == "SIGMAPROFILE":
            total_values = self._sigma_profile_total_values(results, columns, nitems, nhb)
            if total_values:
                for index, sigma_value in enumerate(sigma):
                    row = {
                        "property": property_name,
                        "method": method,
                        "mixture": 0,
                        "cid": np.nan,
                        "name": "Total",
                        "frac1": np.nan,
                        "sigma": sigma_value,
                    }
                    for column in columns:
                        total_column_values = total_values.get(column.column)
                        row[column.column] = np.nan if total_column_values is None else total_column_values[index]
                    rows.append(row)

        return pd.DataFrame(
            rows,
            columns=[
                "property",
                "method",
                "mixture",
                "cid",
                "name",
                "frac1",
                "sigma",
                *(column.column for column in columns),
            ],
        )

    @staticmethod
    def _sigma_profile_frac1_values(results: dict, ncomp: int) -> "np.ndarray":
        frac1 = np.full(ncomp, np.nan)
        if "frac1" not in results:
            return frac1

        values = np.asarray(results["frac1"], dtype=float).reshape(-1)
        frac1[: min(ncomp, values.size)] = values[:ncomp]
        return frac1

    @staticmethod
    def _sigma_profile_component_values(
        results: dict,
        column: CRSSigmaProfileColumn,
        cid: int,
        ncomp: int,
        nitems: int,
        nhb: int,
    ) -> "np.ndarray":
        values = np.asarray(results[column.compound_source], dtype=float)
        if column.hb_channel is None:
            return values.reshape(ncomp, nitems)[cid]
        if nhb == 1 and values.shape == (ncomp, nitems):
            return values[cid]
        return values.reshape(nhb, ncomp, nitems)[column.hb_channel, cid]

    @staticmethod
    def _sigma_profile_total_values(
        results: dict,
        columns: Sequence[CRSSigmaProfileColumn],
        nitems: int,
        nhb: int,
    ) -> Dict[str, "np.ndarray"]:
        total_values = {}
        for column in columns:
            if column.mixture_total_source is None or column.mixture_total_source not in results:
                continue
            values = np.asarray(results[column.mixture_total_source], dtype=float)
            if column.hb_channel is None:
                total_values[column.column] = values.reshape(nitems)
            elif nhb == 1 and values.shape == (nitems,):
                total_values[column.column] = values
            else:
                total_values[column.column] = values.reshape(nhb, nitems)[column.hb_channel]
        return total_values

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
        method = results["method"]
        rows = []
        for mixture in range(nitems):
            for cid in range(ncomp):
                row: Dict[str, Any] = {
                    "property": property_name,
                    "method": method,
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

        return pd.DataFrame(
            rows, columns=list(self._RESULT_TABLE_COMPONENT_BASE_COLUMNS) + ["molmass"] + list(quantities)
        )

    def _build_result_mixture_table(
        self, pd: Any, results: dict, property_name: str, quantities: Sequence[str]
    ) -> Optional["pd.DataFrame"]:
        if not quantities:
            return None

        nitems = int(results["nitems"])
        method = results["method"]
        lle_not_applicable = bool(results.get("isobar", False)) or bool(results.get("flashpoint", False))
        rows = []
        for mixture in range(nitems):
            row: Dict[str, Any] = {"property": property_name, "method": method, "mixture": mixture}
            for quantity in quantities:
                if quantity == "showmiscgap" and lle_not_applicable:
                    row[quantity] = np.nan
                else:
                    row[quantity] = self._get_mixture_quantity_value(results[quantity], mixture, nitems)
            rows.append(row)

        return pd.DataFrame(rows, columns=["property", "method", "mixture"] + list(quantities))

    def _select_lle_columns(self, results: dict, property_name: str, quantities: Sequence[str]) -> Tuple[str, ...]:
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

        columns = self._RESULT_TABLE_LLE_BASE_COLUMNS + self._select_lle_columns(results, property_name, quantities)

        if not quantities:
            return pd.DataFrame(columns=columns)

        method = results["method"]
        if property_name in {"LLE", "STABILITY"}:
            ncomp = int(results["ncomp"])
            nitems = 1
            mixture = 0
            rows = []

            for cid in range(ncomp):
                row = {
                    "property": property_name,
                    "method": method,
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

            xlle = np.asarray(results["xlle"]).ravel()
            xI = [xlle[0], 1.0 - xlle[0]]
            xII = [xlle[1], 1.0 - xlle[1]]
            act_interp = [xlle[6], xlle[7]]
            pressure, temperature = xlle[2], xlle[3]

            rows = [
                {
                    "property": property_name,
                    "method": method,
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

        # case for "TERNARYMIX"
        showmiscgap = bool(results.get("showmiscgap"))
        lle_not_applicable = bool(results.get("isobar", False)) or bool(results.get("flashpoint", False))
        if not showmiscgap or lle_not_applicable:
            return pd.DataFrame(columns=columns)

        nxll = int(results["nxll"])
        nitems = int(results["nitems"])
        xll = np.asarray(results["xll"], dtype=float).ravel().reshape(nxll, 6)
        actxll = np.asarray(results["actxll"], dtype=float).ravel().reshape(nxll, 3)
        temperature = self._get_mixture_quantity_value(results["temperature"], 0, nitems)

        rows = [
            {
                "property": property_name,
                "method": method,
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
    def _get_component_quantity_value(
        value: Any, cid: int, mixture: int, ncomp: int, nitems: int, check: bool = False
    ) -> Any:
        array = np.asarray(value)
        if check:
            print("Component value array shape:", array.shape)
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
        if check:
            print("Mixture value array shape:", array.shape)
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

        if np.sum(Assoc) > 0:  # type: ignore[arg-type]
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
        # try:
        #     from IPython import get_ipython

        #     ipython = get_ipython()
        #     if ipython is not None:
        #         if "zmqshell" in str(type(ipython)):
        #             terminal = "jupyter"
        #         else:
        #             terminal = "interactive"
        #     else:
        #         terminal = "script"
        # except ImportError:
        #     terminal = "script"

        # Check if matplotlib is installed
        plt = CRSResults._import_matplotlib_pyplot(self.__class__.__name__ + ".plot")

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
    def _import_matplotlib_pyplot(method: str) -> Any:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(f"{method}: this method requires the 'matplotlib' package")
        return plt

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

    _METHODS = crs_defs.CRS_METHODS
    _METHOD_ALIASES = {"COSMORS": "COSMO-RS", "COSMOSAC": "COSMOSAC2013"}
    _METHOD_PARAMETERS_METADATA = crs_defs.CRS_METHOD_PARAMETERS_METADATA

    _COMPOUND_KEY_METADATA = get_block_keys(CRS_DATA, "compound")
    _COMPOUND_KEY_SET = get_block_child_names(CRS_DATA, "compound")

    _FORM_KEY_METADATA = get_block_keys(CRS_DATA, "compound", "form")
    _FORM_KEY_SET = get_block_child_names(CRS_DATA, "compound", "form")

    _SPECIES_KEY_METADATA = get_block_keys(CRS_DATA, "compound", "form", "species")
    _SPECIES_KEY_SET = get_block_child_names(CRS_DATA, "compound", "form", "species")

    _STRUCTURE_KEY_METADATA = get_block_keys(CRS_DATA, "compound", "form", "species", "structure")
    _STRUCTURE_KEY_SET = get_block_child_names(CRS_DATA, "compound", "form", "species", "structure")

    _REQUIRED_KEY_METADATA = get_block_keys(CRS_DATA, "compound", "required")
    _REQUIRED_KEY_SET = get_block_child_names(CRS_DATA, "compound", "required")

    def __init__(self, **kwargs: Any) -> None:
        """Initialize a :class:`CRSJob` instance."""
        super().__init__(**kwargs)
        self.settings.ignore_molecule = True

    @staticmethod
    def database() -> str:
        """Return the installed ADFCRS-2018 database directory."""
        database_path = Path(os.environ["SCM_PKG_ADFCRSDIR"]) / "ADFCRS-2018"
        if not database_path.is_dir():
            raise FileNotFoundError("The ADFCRS-2018 database does not seem to be installed")
        return os.fspath(database_path)

    @staticmethod
    def coskf_from_database(name: str) -> str:
        """Return an existing COSKF file path from the installed ADFCRS-2018 database."""
        if not name.endswith(".coskf"):
            name += ".coskf"

        path = Path(CRSJob.database()) / name
        if not path.is_file():
            raise FileNotFoundError(f"COSKF file not found in ADFCRS-2018 database: {name}")

        return os.fspath(path)

    @staticmethod
    def property_types() -> Tuple[CRSPropertyName, ...]:
        """Return the supported COSMO-RS property types."""
        return cast(Tuple[CRSPropertyName, ...], tuple(_CRS_INPUT_BUILDER_CLASSES))

    @staticmethod
    def _property_type_metadata() -> Dict[str, Dict[str, Any]]:
        """Return CRS property metadata derived from property-specific builders."""
        return {
            property_type: builder_class._metadata()
            for property_type, builder_class in _CRS_INPUT_BUILDER_CLASSES.items()
        }

    @staticmethod
    def _property_type_builder_class(property_type: str) -> type["CRSInputBuilder"]:
        normalized = CRSJob._normalize_property_type(property_type)
        return _CRS_INPUT_BUILDER_CLASSES[normalized]

    @staticmethod
    def property_type_metadata(
        property_type: str,
        as_summary: bool = True,
    ) -> Union[Dict[str, Any], Tuple[str, ...]]:
        """Return discoverability metadata for a COSMO-RS problem type."""
        metadata = CRSJob._property_type_builder_class(property_type)._metadata()
        if not as_summary:
            return deepcopy(metadata)

        summary = [
            f"{property_type}: {metadata['description']}",
            f"system_scope: {metadata['system_scope']}",
            f"top_level_keys: {', '.join(metadata['input_keys']['top_level']) or '-'}",
            f"property_keys: {', '.join(metadata['input_keys']['property']) or '-'}",
            f"compound_keys: {', '.join(metadata['input_keys']['compound']) or '-'}",
            f"required_keys: {', '.join(metadata['required_keys']) or '-'}",
        ]
        for comment in metadata.get("comments", []):
            summary.append(f"comment: {comment}")

        summary.extend(CRSJob._property_type_input_hint_descriptions(metadata))
        return tuple(summary)

    @staticmethod
    def _property_type_input_hint_descriptions(metadata: Dict[str, Any]) -> Tuple[str, ...]:
        """Return formatted input hints selected by the property-type metadata."""
        key_metadata = CRSJob._property_type_input_key_metadata(metadata)
        grouped: Dict[str, List[str]] = {}

        for key in metadata.get("hint_keys", ()):
            hint = key_metadata.get(key, {}).get("hint", "")
            if hint:
                grouped.setdefault(hint, []).append(key)

        return tuple(
            f"hint [{', '.join(keys)}]: {hint}"
            for hint, keys in grouped.items()
        )

    @staticmethod
    def _property_type_input_key_metadata(metadata: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """Return input-definition metadata for all keys declared by a property type."""
        key_metadata: Dict[str, Dict[str, Any]] = {}

        for key in metadata["input_keys"].get("top_level", ()):
            key_metadata[key] = CRS_DATA[key]

        property_key_metadata = get_block_keys(CRS_DATA, "property")
        for key in metadata["input_keys"].get("property", ()):
            key_metadata[key] = property_key_metadata[key]

        compound_key_metadata = get_block_keys(CRS_DATA, "compound")
        for key in metadata["input_keys"].get("compound", ()):
            key_metadata[key] = compound_key_metadata[key]

        return key_metadata

    @staticmethod
    def _format_key_metadata(key: str, metadata: Dict[str, str]) -> str:
        """Return a compact display string for one metadata entry."""
        type_unit = metadata.get("type", "")
        if metadata.get("unit"):
            type_unit = f"{type_unit} [{metadata['unit']}]"
        return f"{key}: {type_unit} - {metadata['description']}"

    @staticmethod
    def compound_keys(
        key: Optional[str] = None,
        *,
        as_summary: bool = True,
    ) -> Union[Dict[str, Any], Tuple[str, ...], str]:
        """Return COMPOUND key descriptions, or raw metadata with ``as_summary=False``."""
        if key is None:
            metadata = {key: value.copy() for key, value in CRSJob._COMPOUND_KEY_METADATA.items()}
            if as_summary:
                return tuple(CRSJob._format_key_metadata(key, value) for key, value in metadata.items())
            return metadata

        normalized = key.lower()
        if normalized not in CRSJob._COMPOUND_KEY_METADATA:
            allowed = ", ".join(sorted(CRSJob._COMPOUND_KEY_METADATA))
            raise KeyError(f"Unknown COMPOUND key {key!r}. Supported keys: {allowed}")

        metadata = CRSJob._COMPOUND_KEY_METADATA[normalized].copy()
        if as_summary:
            return CRSJob._format_key_metadata(normalized, metadata)
        return metadata

    @staticmethod
    def form_keys(as_summary: bool = True) -> Union[Dict[str, Dict[str, str]], Tuple[str, ...]]:
        """Return FORM key descriptions, or raw metadata with ``as_summary=False``."""
        metadata = {key: value.copy() for key, value in CRSJob._FORM_KEY_METADATA.items()}
        if as_summary:
            return tuple(CRSJob._format_key_metadata(key, value) for key, value in metadata.items())
        return metadata

    @staticmethod
    def species_keys(as_summary: bool = True) -> Union[Dict[str, Dict[str, str]], Tuple[str, ...]]:
        """Return SPECIES key descriptions, or raw metadata with ``as_summary=False``."""
        metadata = {key: value.copy() for key, value in CRSJob._SPECIES_KEY_METADATA.items()}
        if as_summary:
            return tuple(CRSJob._format_key_metadata(key, value) for key, value in metadata.items())
        return metadata

    @staticmethod
    def structure_keys(as_summary: bool = True) -> Union[Dict[str, Dict[str, str]], Tuple[str, ...]]:
        """Return STRUCTURE key descriptions, or raw metadata with ``as_summary=False``."""
        metadata = {key: value.copy() for key, value in CRSJob._STRUCTURE_KEY_METADATA.items()}
        if as_summary:
            return tuple(CRSJob._format_key_metadata(key, value) for key, value in metadata.items())
        return metadata

    @staticmethod
    def methods() -> Tuple[str, ...]:
        """Return the supported COSMO-RS/SAC methods."""
        return CRSJob._METHODS

    @staticmethod
    def _normalize_choice(
        value: str,
        *,
        value_name: str,
        allowed: Sequence[str],
        aliases: Optional[Dict[str, str]] = None,
    ) -> str:
        """Normalize and validate a string choice against supported values."""
        if not isinstance(value, str):
            raise TypeError(f"{value_name} must be a string, got {type(value).__name__}")

        normalized = value.upper()
        if aliases is not None:
            normalized = aliases.get(normalized, normalized)
        if normalized not in allowed:
            supported = ", ".join(allowed)
            raise ValueError(f"Unsupported {value_name} {value!r}. Supported values: {supported}")
        return normalized

    @staticmethod
    def _normalize_property_type(property_type: str) -> str:
        """Normalize and validate a COSMO-RS problem type."""
        return CRSJob._normalize_choice(
            property_type,
            value_name="property_type",
            allowed=CRSJob.property_types(),
        )

    @staticmethod
    def _set_non_none(block: Settings, values: Dict[str, Any]) -> Settings:
        for key, value in values.items():
            if value is not None:
                block[key] = value
        return block

    @staticmethod
    def _set_parameter_block(
        block: Settings,
        parameters: Dict[str, Any],
        lowercase_keys: bool = False,
    ) -> None:
        for key, value in parameters.items():
            if value is None:
                continue
            block[key.lower() if lowercase_keys else key] = value

    @staticmethod
    def method_block(
        method: str = "COSMO-RS",
        *,
        include_defaults: bool = False,
        parameter_set: Optional[str] = None,
        crsparameters: Optional[Dict[str, Any]] = None,
        sacparameters: Optional[Dict[str, Any]] = None,
        dispersion: Optional[Dict[str, float]] = None,
        epsilon: Optional[Dict[str, float]] = None,
    ) -> Settings:
        """Return method-selection and method-parameter settings for COSMO-RS/SAC.
        The method name and validated parameter names are case-insensitive.
        """
        sett = Settings()
        normalized_method = CRSJob._normalize_method(method)
        sett.input.method = normalized_method

        defaults = (
            CRSJob._copy_method_parameter_defaults(normalized_method, parameter_set)
            if include_defaults or parameter_set is not None
            else {}
        )

        user_blocks = {
            "CRSParameters": crsparameters,
            "SACParameters": sacparameters,
            "Dispersion": dispersion,
            "Epsilon": epsilon,
        }

        for block_name, user_parameters in user_blocks.items():
            if user_parameters is not None and not isinstance(user_parameters, dict):
                raise TypeError(f"{block_name} parameters must be a dictionary")

            parameters = {
                **defaults.get(block_name, {}),
                **(user_parameters or {}),
            }

            if parameters:
                CRSJob._set_parameter_block(
                    sett.input[block_name],
                    parameters,
                    lowercase_keys=block_name in {"CRSParameters", "SACParameters"},
                )
        return sett

    @staticmethod
    def _normalize_method(method: str) -> str:
        """Normalize and validate a COSMO-RS/SAC method name."""
        return CRSJob._normalize_choice(
            method,
            value_name="method",
            allowed=CRSJob._METHODS,
            aliases=CRSJob._METHOD_ALIASES,
        )

    @staticmethod
    def method_options(
        method: Optional[str] = None,
        *,
        exposed_only: bool = True,
        as_df: bool = False,
    ) -> Union[Tuple[Dict[str, Any], ...], "pd.DataFrame"]:
        metadata = CRSJob._METHOD_PARAMETERS_METADATA

        normalized_method = CRSJob._normalize_method(method) if method is not None else None

        rows = []
        for parameter_set, preset in metadata.items():
            if exposed_only and not preset.get("exposed", False):
                continue
            if normalized_method is not None and preset.get("method") != normalized_method:
                continue

            rows.append({
                "method": preset.get("method"),
                "parameter_set": parameter_set,
                "default": bool(preset.get("default", False)),
                "optimized_for_adf_sigma_profile": bool(
                    preset.get("optimized_for_adf_sigma_profile", False)
                ),
            })

        rows.sort(key=lambda row: (row["method"] or "", row["parameter_set"]))

        if as_df:
            pd = CRSResults._import_pandas(
                "parameter_set_table",
                "as_df=True requires the 'pandas' package",
            )
            return pd.DataFrame(rows)

        return tuple(rows)

    # @staticmethod
    # def _method_parameter_keys(method_name: str, block_name: str) -> Set[str]:
    #     """Return allowed parameter names for a method parameter block."""
    #     keys: Set[str] = set()

    #     for preset in CRSJob._METHOD_PARAMETERS_METADATA.values():
    #         if preset.get("method") != method_name:
    #             continue

    #         block = preset.get("parameter_blocks", {}).get(block_name)
    #         if block is None:
    #             continue

    #         keys.update(str(key).lower() for key in block)

    #     return keys

    # @staticmethod
    # def _set_method_parameters(
    #     block: Settings,
    #     parameters: Dict[str, Any],
    #     block_name: str,
    #     method_name: str,
    # ) -> None:
    #     """Set validated method parameters on a CRS method subblock."""

    #     allowed_keys = CRSJob._method_parameter_keys(method_name, block_name)
    #     for key, value in parameters.items():
    #         normalized_key = key.lower()
    #         if normalized_key not in allowed_keys:
    #             allowed = ", ".join(sorted(allowed_keys))
    #             raise ValueError(f"Unsupported {block_name} key {key!r}. Allowed keys: {allowed}")
    #         if value is not None:
    #             block[normalized_key] = value

    @staticmethod
    def _copy_method_parameter_defaults(
        method: str,
        parameter_set: Optional[str] = None,
    ) -> Dict[str, Dict[str, Any]]:
        metadata = CRSJob._METHOD_PARAMETERS_METADATA

        def copy_parameter_blocks(preset: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
            return {
                block_name: values.copy()
                for block_name, values in preset["parameter_blocks"].items()
            }

        if parameter_set is not None:
            parameter_set_key = parameter_set.strip().lower()
            preset = metadata.get(parameter_set_key)
            if preset is None:
                raise ValueError(f"Unknown CRS parameter_set {parameter_set_key!r}")
            if preset["method"] != method:
                raise ValueError(
                    f"CRS parameter_set {parameter_set_key!r} is for method "
                    f"{preset['method']!r}, not {method!r}"
                )
            return copy_parameter_blocks(preset)

        method_presets = [
            preset
            for preset in metadata.values()
            if preset["method"] == method
        ]

        if not method_presets:
            return {}

        if len(method_presets) == 1:
            return copy_parameter_blocks(method_presets[0])

        default_presets = [
            preset
            for preset in method_presets
            if preset.get("default", False)
        ]

        if len(default_presets) == 1:
            return copy_parameter_blocks(default_presets[0])

        raise ValueError(
            f"Multiple parameter sets are available for method {method!r}: "
            f"{', '.join(repr(preset['parameter_set']) for preset in method_presets)}. "
            "Specify one explicitly with parameter_set=..."
        )

    @staticmethod
    def compound_block(
        path: Optional[PathLike] = None,
        *,
        forms: Optional[Union[Settings, Sequence[Settings]]] = None,
        name: Optional[str] = None,
        frac1: Optional[float] = None,
        frac2: Optional[float] = None,
        meltingpoint: Optional[float] = None,
        hfusion: Optional[float] = None,
        cpfusion: Optional[float] = None,
        pvap: Optional[float] = None,
        tvap: Optional[float] = None,
        vp_equation: Optional[str] = None,
        vp_params: Optional[str] = None,
        density: Optional[float] = None,
        flashpoint: Optional[float] = None,
        dielectric_const: Optional[float] = None,
        polymer: Optional[bool] = None,
        averagemwpoly: Optional[float] = None,
        **kwargs: Any,
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
        path = CRSJob._validate_path_or_nested_entries("compound_block", path, forms)

        compound = Settings()

        if path is not None:
            compound._h = path
        else:
            compound.form = forms

        nring = CRSJob._infer_nring_if_missing(path) if path is not None else None

        compound_values = {
            "name": name,
            "frac1": frac1,
            "frac2": frac2,
            "nring": nring or None,
            "meltingpoint": meltingpoint,
            "hfusion": hfusion,
            "cpfusion": cpfusion,
            "pvap": pvap,
            "tvap": tvap,
            "vp_equation": vp_equation,
            "vp_params": vp_params,
            "density": density,
            "polymer": polymer,
            "averagemwpoly": averagemwpoly,
            "flashpoint": flashpoint,
            "dielectric_const": dielectric_const,
        }
        compound_values.update(kwargs)
        CRSJob._set_non_none(compound, compound_values)

        return CRSJob._normalize_compound(compound)

    @staticmethod
    def form_block(
        path: Optional[PathLike] = None,
        *,
        species: Optional[Union[Settings, Sequence[Settings]]] = None,
        name: Optional[str] = None,
        count: Optional[float] = None,
        **kwargs: Any,
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
        path = CRSJob._validate_path_or_nested_entries("form_block", path, species)

        form = Settings()

        if path is not None:
            form._h = path
        else:
            form.species = species

        nring = CRSJob._infer_nring_if_missing(path) if path is not None else None

        form_values = {
            "name": name,
            "count": count,
            "nring": nring or None,
        }
        form_values.update(kwargs)
        CRSJob._set_non_none(form, form_values)

        return CRSJob._normalize_form(form)

    @staticmethod
    def species_block(
        path: Optional[PathLike] = None,
        *,
        structures: Optional[Union[Settings, Sequence[Settings]]] = None,
        name: Optional[str] = None,
        count: Optional[float] = None,
        **kwargs: Any,
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
        path = CRSJob._validate_path_or_nested_entries("species_block", path, structures)

        species = Settings()

        if path is not None:
            species._h = path
        else:
            species.structure = structures

        nring = CRSJob._infer_nring_if_missing(path) if path is not None else None

        species_values = {
            "name": name,
            "count": count,
            "nring": nring or None,
        }
        species_values.update(kwargs)
        CRSJob._set_non_none(species, species_values)

        return CRSJob._normalize_species(species)

    @staticmethod
    def structure_block(
        path: PathLike,
        *,
        name: Optional[str] = None,
        count: Optional[float] = None,
        **kwargs: Any,
    ) -> Settings:
        """Create a STRUCTURE settings block from a structure path.

        Example::

            structure_a = CRSJob.structure_block(
                "poly_anion_a.coskf",
                name="poly_anion conformer A",
                count=1.0,
                Hcorr=0.0,
                Scorr=0.0,
            )
        """
        path = CRSJob._validate_path("structure_block", path)

        structure = Settings()
        structure._h = path

        nring = CRSJob._infer_nring_if_missing(path)

        structure_values = {
            "name": name,
            "count": count,
            "nring": nring or None,
        }
        structure_values.update(kwargs)
        CRSJob._set_non_none(structure, structure_values)

        return CRSJob._normalize_structure(structure)

    @staticmethod
    def _normalize_path(path: Optional[PathLike]) -> Optional[str]:
        return None if path is None else os.fspath(path)

    @staticmethod
    def _validate_path(method_name: str, path: PathLike) -> str:
        normalized_path = CRSJob._normalize_path(path)
        if not normalized_path:
            raise ValueError(f"{method_name} requires a non-empty path")
        absolute_path = Path(normalized_path).expanduser().resolve()
        if not absolute_path.is_file():
            raise FileNotFoundError(f"{method_name} path does not exist: {absolute_path}")
        return str(absolute_path)

    @staticmethod
    def _validate_path_or_nested_entries(
        method_name: str,
        path: Optional[PathLike],
        nested_entries: Any,
    ) -> Optional[str]:
        """Validate that exactly one of a path or nested entries was provided."""
        if path is not None and nested_entries is not None:
            raise ValueError(f"{method_name} accepts either path or nested entries, not both")
        if path is None and nested_entries is None:
            raise ValueError(f"{method_name} requires either path or nested entries")
        if path is not None:
            return CRSJob._validate_path(method_name, path)

        return None

    @staticmethod
    def _validate_keys(block_name: str, keys: Set[Any], allowed_keys: frozenset[str]) -> None:
        invalid_keys = sorted(
            str(key) for key in keys if key != "_h" and str(key).lower() not in allowed_keys
        )
        if invalid_keys:
            allowed = ", ".join(sorted(allowed_keys))
            invalid = ", ".join(invalid_keys)
            raise ValueError(f"Unsupported {block_name} key(s): {invalid}. Allowed keys: {allowed}")

    @staticmethod
    def _validate_settings_keys(block_name: str, block: Settings, allowed_keys: frozenset[str]) -> None:
        CRSJob._validate_keys(block_name, set(block.keys()), allowed_keys)

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
    def _normalize_nested_block(
        block: Settings,
        *,
        block_name: str,
        allowed_keys: frozenset[str],
        nested_key: Optional[str] = None,
        nested_input_name: Optional[str] = None,
        nested_normalizer: Optional[Callable[[Settings], Settings]] = None,
    ) -> Settings:
        """Validate a CRS nested block and recursively normalize child blocks."""
        if not isinstance(block, Settings):
            raise TypeError(f"{block_name} entries must be Settings instances, got {type(block).__name__}")

        path = block.get("_h")
        has_path = bool(path)

        if nested_key is None:
            if not has_path:
                raise ValueError(f"{block_name} entries must define a non-empty _h header/path")
        else:
            has_nested = nested_key in block and block[nested_key] not in (None, False)
            if has_path == has_nested:
                raise ValueError(f"{block_name} entries must define exactly one of _h or {nested_key}")

        if "_ipython_canary_method_should_not_exist_" in block:
            block.pop("_ipython_canary_method_should_not_exist_")

        CRSJob._validate_settings_keys(block_name, block, allowed_keys)

        if nested_key is not None and has_nested:
            if nested_normalizer is None:
                raise TypeError(f"{block_name} nested_normalizer is required when nested_key is set")
            block[nested_key] = [
                nested_normalizer(entry)
                for entry in CRSJob._ensure_settings_list(block[nested_key], nested_input_name or nested_key)
            ]
        return block

    @staticmethod
    def _normalize_compound(compound: Settings) -> Settings:
        """Validate and normalize a COMPOUND block."""
        return CRSJob._normalize_nested_block(
            compound,
            block_name="COMPOUND",
            allowed_keys=CRSJob._COMPOUND_KEY_SET,
            nested_key="form",
            nested_input_name="forms",
            nested_normalizer=CRSJob._normalize_form,
        )

    @staticmethod
    def _normalize_form(form: Settings) -> Settings:
        """Validate and normalize a FORM block for use in multispecies compounds."""
        return CRSJob._normalize_nested_block(
            form,
            block_name="FORM",
            allowed_keys=CRSJob._FORM_KEY_SET,
            nested_key="species",
            nested_input_name="species",
            nested_normalizer=CRSJob._normalize_species,
        )

    @staticmethod
    def _normalize_species(species: Settings) -> Settings:
        """Validate and normalize a SPECIES block nested under FORM."""
        return CRSJob._normalize_nested_block(
            species,
            block_name="SPECIES",
            allowed_keys=CRSJob._SPECIES_KEY_SET,
            nested_key="structure",
            nested_input_name="structure",
            nested_normalizer=CRSJob._normalize_structure,
        )

    @staticmethod
    def _normalize_structure(structure: Settings) -> Settings:
        """Validate and normalize a STRUCTURE block nested under SPECIES."""
        return CRSJob._normalize_nested_block(
            structure,
            block_name="STRUCTURE",
            allowed_keys=CRSJob._STRUCTURE_KEY_SET,
        )

    @staticmethod
    def _infer_nring_if_missing(path: str) -> Optional[int]:
        try:
            return CRSJob._read_or_determine_nring(path)
        except Exception:
            return None

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
                log(
                    f"Could not read Nring from COSKF file {coskf_file}; determining it from the molecular graph instead: {exc}",
                    level=3,
                )

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
    def _as_coskf_database(database: Any) -> Tuple[Any, bool]:
        """Return a ``COSKFDatabase`` instance and whether this helper opened it."""
        from pyCRS.Database import COSKFDatabase

        if isinstance(database, COSKFDatabase):
            return database, False

        if isinstance(database, (str, os.PathLike)):
            return COSKFDatabase(os.fspath(database)), True

        raise TypeError("database must be a COSKFDatabase instance or a database path")

    @staticmethod
    def _get_single_database_row(rows_by_identifier: Dict[str, List[Any]], identifier: str) -> Any:
        rows = rows_by_identifier.get(identifier, [])
        valid_rows = [row for row in rows if getattr(row, "compound_id", None) is not None]

        if not valid_rows:
            raise ValueError(f"No compound found in the COSKFDatabase for identifier {identifier!r}")
        if len(valid_rows) > 1:
            compound_ids = ", ".join(str(row.compound_id) for row in valid_rows)
            raise ValueError(
                f"Multiple compounds found in the COSKFDatabase for identifier {identifier!r}: {compound_ids}"
            )

        return valid_rows[0]

    @staticmethod
    def _database_compound_property_values(
        database: Any,
        compound_id: int,
        property_keys: Sequence[str],
        property_sources: Sequence[str],
    ) -> Dict[str, Any]:
        """Return CRS compound keyword values read from COSKFDatabase property tables."""
        db_to_crs_keys = {
            "dielectricconstant": "dielectric_const",
            "mn": "averagemwpoly",
        }
        crs_to_db_keys = {
            "dielectric_const": "dielectricconstant",
            "averagemwpoly": "Mn",
        }

        values: Dict[str, Any] = {}
        for source in property_sources:
            rows = database.get_physical_properties(compound_id=compound_id, source=source)
            if not rows:
                continue

            row = rows[0]
            for crs_key in property_keys:
                if crs_key in values:
                    continue

                db_key = crs_to_db_keys.get(crs_key, crs_key)
                if not hasattr(row, db_key):
                    continue

                value = getattr(row, db_key)
                if value is None:
                    continue

                normalized_key = db_to_crs_keys.get(db_key.lower(), crs_key)
                if normalized_key == "vp_params" and isinstance(value, str):
                    value = value.replace(",", " ")
                values[normalized_key] = value

        return values

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

    @staticmethod
    def compound_block_from_coskf_database(
        identifier: str,
        database: Any,
        *,
        conformer: bool = False,
        property_type: Optional[str] = None,
        property_sources: Sequence[str] = ("PhysicalProperty", "PropPred"),
        property_keys: Optional[Sequence[str]] = None,
        **compound_overrides: Any,
    ) -> Settings:
        """Create a COMPOUND block from a :class:`pyCRS.Database.COSKFDatabase` entry.

        ``database`` can be either a ``COSKFDatabase`` instance or a path to one.
        Imports from ``pyCRS`` are intentionally lazy to avoid creating an import
        cycle between PLAMS ``CRSJob`` and ``pyCRS.CRSManager``.

        Explicit keyword overrides are applied after database-derived values.
        """
        if not isinstance(identifier, str):
            raise TypeError(f"identifier must be a string, got {type(identifier).__name__}")

        invalid_override_keys = {"path", "forms"}
        invalid_overrides = sorted(invalid_override_keys.intersection(compound_overrides))
        if invalid_overrides:
            invalid = ", ".join(invalid_overrides)
            raise ValueError(f"compound_block_from_coskf_database does not accept override(s): {invalid}")

        normalized_property_type = None
        if property_type is not None:
            normalized_property_type = CRSJob._normalize_property_type(property_type)

        if property_keys is None:
            if normalized_property_type is None:
                requested_property_keys: Tuple[str, ...] = ()
            else:
                metadata = CRSJob._property_type_builder_class(normalized_property_type)._metadata()
                requested_property_keys = tuple(
                    key
                    for key in metadata["input_keys"]["compound"]
                    if key not in {"name", "frac1", "frac2", "nring"}
                )
        else:
            requested_property_keys = tuple(property_keys)

        db, close_database = CRSJob._as_coskf_database(database)
        try:
            if conformer:
                rows_by_identifier = db.get_conformers(identifier)
                rows = rows_by_identifier.get(identifier, [])
                valid_rows = [row for row in rows if getattr(row, "compound_id", None) is not None]
                if not valid_rows:
                    raise ValueError(f"No conformers found in the COSKFDatabase for identifier {identifier!r}")

                first_row = valid_rows[0]
                compound_values = CRSJob._database_compound_property_values(
                    db,
                    first_row.compound_id,
                    requested_property_keys,
                    property_sources,
                )
                compound_values.setdefault("name", first_row.name)
                compound_values.setdefault("nring", first_row.nring)
                compound_values.update(compound_overrides)

                forms = [
                    CRSJob.form_block(
                        row.get_full_coskf_path(),
                        name=row.name,
                        nring=row.nring,
                    )
                    for row in valid_rows
                ]
                return CRSJob.compound_block(forms=forms, **compound_values)

            rows_by_identifier = db.get_compounds(identifier)
            row = CRSJob._get_single_database_row(rows_by_identifier, identifier)
            compound_values = CRSJob._database_compound_property_values(
                db,
                row.compound_id,
                requested_property_keys,
                property_sources,
            )
            compound_values.setdefault("name", row.name)
            compound_values.setdefault("nring", row.nring)
            compound_values.update(compound_overrides)

            return CRSJob.compound_block(row.get_full_coskf_path(), **compound_values)

        finally:
            if close_database:
                db._close_connection()

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["ACTIVITYCOEF"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "ACTIVITYCOEFInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["LOGP"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "LOGPInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["SOLUBILITY"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[CRSSolubilityMode] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "SOLUBILITYInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["PURESOLUBILITY"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[CRSSolubilityMode] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "PURESOLUBILITYInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["VAPORPRESSURE"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "VAPORPRESSUREInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["PUREVAPORPRESSURE"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "PUREVAPORPRESSUREInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["BOILINGPOINT"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "BOILINGPOINTInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["PUREBOILINGPOINT"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "PUREBOILINGPOINTInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["FLASHPOINT"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "FLASHPOINTInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["BINMIXCOEF"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[CRSVLESweepMode] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "BINMIXCOEFInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["TERNARYMIX"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[CRSVLESweepMode] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "TERNARYMIXInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["COMPOSITIONLINE"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[CRSVLESweepMode] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "COMPOSITIONLINEInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["LLE"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "LLEInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["STABILITY"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "STABILITYInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["SIGMAPROFILE"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "SIGMAPROFILEInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["PURESIGMAPROFILE"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "PURESIGMAPROFILEInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["SIGMAPOTENTIAL"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "SIGMAPOTENTIALInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: Literal["PURESIGMAPOTENTIAL"],
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "PURESIGMAPOTENTIALInputBuilder": ...

    @staticmethod
    @overload
    def input_builder(
        property_type: CRSPropertyName,
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "CRSInputBuilder": ...

    @staticmethod
    def input_builder(
        property_type: CRSPropertyName,
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "CRSInputBuilder":
        """Return a property-specific CRS input builder.

        Example::

            builder = CRSJob.input_builder(
                "ACTIVITYCOEF",
                method="COSMOSAC2013",
                mode="gas",
                temperature="353.15 373.15 10",
            )

            builder.add_solvent(CRSJob.coskf_from_database("Water.coskf"), frac1=1.0)
            builder.add_solute(CRSJob.coskf_from_database("Benzene.coskf"))

            settings = builder.to_settings()

        Use ``builder.describe()`` to inspect supported keys, modes, and compound roles.
        """
        normalized_property_type = CRSJob._normalize_property_type(property_type)
        builder_class = _CRS_INPUT_BUILDER_CLASSES[normalized_property_type]
        return builder_class(
            method=method,
            mode=mode,
            use_defaults=use_defaults,
            **kwargs,
        )


class CRSInputBuilder:
    """Property-type builder for CRSJob input settings."""

    _PROPERTY_TYPE: ClassVar[Optional[str]] = None
    _DESCRIPTION: ClassVar[str] = ""
    _SYSTEM_SCOPE: ClassVar[str] = ""
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": (),
        "property": (),
        "compound": (),
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ()
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ()
    _COMMENTS: ClassVar[Tuple[str, ...]] = ()
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ()
    _COMPOUND_ROLE_CONFIG: ClassVar[Dict[str, Dict[str, Optional[int]]]] = {}
    _MODE_CONFIG: ClassVar[Dict[str, Any]] = {}

    _AUTO_DEFAULT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "pressure"),
        "property": ("nfrac", "nprofile", "sigmamax"),
    }

    _GLOBAL_TOP_LEVEL_KEYS: ClassVar[Tuple[str, ...]] = (
        "pdh_correction",
        "usepolycombiforpolymer",
    )

    _COMMON_DIR_ENTRIES: ClassVar[Tuple[str, ...]] = (
        "property_type",
        "method",
        "method_options",
        "set",
        "get",
        "describe",
        "apply_defaults",
        "to_settings",
        "to_job",
    )

    _INTERNAL_ATTRIBUTES: ClassVar[Set[str]] = {
        "_method",
        "_accepted_keys",
        "_mode",
        "_mode_config",
        "_active_mode_input_keys",
        "_compound_roles",
        "_compounds_by_role",
    }

    def __init__(
        self,
        property_type: Optional[CRSPropertyName] = None,
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> None:
        normalized = type(self)._resolve_property_type(property_type)

        metadata = type(self)._metadata()
        builder_metadata = metadata["builder"]

        accepted_keys = self._build_accepted_keys(metadata)
        mode_config = type(self)._mode_config()

        compound_roles = tuple(builder_metadata.get("roles", ()))
        compounds_by_role = {role: [] for role in compound_roles}

        object.__setattr__(self, "property_type", normalized)
        object.__setattr__(self, "metadata", metadata)
        object.__setattr__(self, "_accepted_keys", accepted_keys)
        object.__setattr__(self, "_method", CRSJob._normalize_method("COSMO-RS"))

        # User modes can set CRS input keys; for example, mode="gas" sets isobar=True.
        object.__setattr__(self, "_mode", None)
        object.__setattr__(self, "_mode_config", mode_config)
        object.__setattr__(self, "_active_mode_input_keys", set())

        object.__setattr__(self, "_compound_roles", compound_roles)
        object.__setattr__(self, "_compounds_by_role", compounds_by_role)

        self._reject_mode_controlled_keys(kwargs)
        self._set_method(method)
        self._set_mode_or_default(mode)

        for key, value in kwargs.items():
            self._set_input_value(key, value)

        if use_defaults:
            self.apply_defaults()

    @classmethod
    def _resolve_property_type(cls, property_type: Optional[CRSPropertyName]) -> str:
        if cls._PROPERTY_TYPE is None:
            raise TypeError("CRSInputBuilder is a base class; use CRSJob.input_builder()")
        if property_type is None:
            return cls._PROPERTY_TYPE

        normalized = CRSJob._normalize_property_type(property_type)
        if normalized != cls._PROPERTY_TYPE:
            raise ValueError(f"{cls.__name__} only supports property_type={cls._PROPERTY_TYPE!r}")
        return normalized

    @classmethod
    def _metadata(cls) -> Dict[str, Any]:
        """Return the property metadata declared by this property-specific builder class."""
        if cls._PROPERTY_TYPE is None:
            raise TypeError("CRSInputBuilder base class does not define property metadata")

        input_keys = {
            "top_level": tuple(cls._INPUT_KEYS.get("top_level", ())),
            "property": tuple(cls._INPUT_KEYS.get("property", ())),
            "compound": tuple(cls._INPUT_KEYS.get("compound", ())),
        }
        builder: Dict[str, Any] = {"roles": tuple(cls._COMPOUND_ROLES)}
        if cls._COMPOUND_ROLE_CONFIG:
            builder["role_config"] = deepcopy(cls._COMPOUND_ROLE_CONFIG)
        if cls._MODE_CONFIG:
            builder["mode"] = deepcopy(cls._MODE_CONFIG)

        return {
            "description": cls._DESCRIPTION,
            "system_scope": cls._SYSTEM_SCOPE,
            "input_keys": input_keys,
            "required_keys": tuple(cls._REQUIRED_KEYS),
            "comments": tuple(cls._COMMENTS),
            "hint_keys": tuple(cls._HINT_KEYS),
            "builder": builder,
        }

    @classmethod
    def _mode_config(cls) -> CRSModeConfig:
        return CRSModeConfig(cls._MODE_CONFIG)

    @classmethod
    def _build_accepted_keys(cls, metadata: Dict[str, Any]) -> Dict[str, CRSInputRoute]:
        """Build public builder-key routes for top-level and PROPERTY input keys."""
        accepted_keys = {
            key: CRSInputRoute(scope="top_level", metadata=CRS_DATA[key])
            for key in metadata["input_keys"].get("top_level", ()) + cls._GLOBAL_TOP_LEVEL_KEYS
        }

        property_key_metadata = get_block_keys(CRS_DATA, "property")
        accepted_keys.update(
            {
                key: CRSInputRoute(scope="property", metadata=property_key_metadata[key])
                for key in metadata["input_keys"].get("property", ())
            }
        )

        return accepted_keys

    def __dir__(self) -> List[str]:
        """Return property-type-specific completions for interactive use."""
        entries = set()
        entries.update(self._accepted_keys)
        entries.update(self._COMMON_DIR_ENTRIES)

        if self._mode_config.options:
            entries.difference_update(self._mode_config.options)
            entries.difference_update(self._active_mode_input_keys)
            entries.update({"mode", "mode_options"})

        if "compound" in self._compound_roles and hasattr(type(self), "add_compound"):
            entries.add("add_compound")
        if "solvent" in self._compound_roles and hasattr(type(self), "add_solvent"):
            entries.add("add_solvent")
        if "solute" in self._compound_roles and hasattr(type(self), "add_solute"):
            entries.add("add_solute")

        return sorted(entries)

    def __setattr__(self, key: str, value: Any) -> None:
        if key in self._INTERNAL_ATTRIBUTES:
            object.__setattr__(self, key, value)
            return

        self._set_input_value(key, value)

    @property
    def method(self) -> str:
        return self._method

    @method.setter
    def method(self, method: CRSMethodName) -> None:
        self._set_method(method)

    @property
    def method_options(self) -> Tuple[str, ...]:
        return CRSJob.methods()

    def get(self, key: str, default: Any = None) -> Any:
        """Return an input value if set, otherwise return default."""
        if key == "method":
            return self._method
        if key == "method_options":
            return CRSJob.methods()
        if key == "mode":
            return self._mode
        if key == "mode_options":
            return self._mode_config.options
        self._validate_input_key(key)
        return self._get_value(key) if self._has_value(key) else default

    def set(self: _CRSInputBuilderT, **kwargs: Any) -> _CRSInputBuilderT:
        """Set one or more property-type-aware CRS input values."""
        self._reject_mode_controlled_keys(kwargs)
        for key, value in kwargs.items():
            self._set_input_value(key, value)
        return self

    def apply_defaults(self) -> "CRSInputBuilder":
        """Set auto-default input values for unset supported keys."""
        for key, route in sorted(self._accepted_keys.items()):
            if self._has_value(key):
                continue

            has_default, value = self._auto_default_candidate(key, route)
            if not has_default:
                continue

            self._set_value(key, value)

        return self

    def describe(
        self,
        include_values: bool = False,
        include_compound_details: bool = False,
        include_guidance: bool = False,
    ) -> Tuple[str, ...]:
        """Return property-type-specific input keys and mode descriptions."""

        descriptions = [
            self._method_description(include_values=include_values),
            *self._input_descriptions(include_values=include_values),
            *self._compound_descriptions(include_details=include_compound_details),
            *self._mode_descriptions(),
        ]

        if include_guidance:
            descriptions.extend(self._comment_descriptions())
            descriptions.extend(self._input_hint_descriptions())

        return tuple(descriptions)

    def to_settings(self, *, include_compounds: bool = True) -> Settings:
        """Build PLAMS Settings, optionally omitting COMPOUND blocks."""
        self._validate_required_inputs()
        if include_compounds:
            self._validate_compound_role_counts(check="bounds")
            self._validate_compound_required_keys_by_mode()

        settings = Settings()
        settings.input.method = self._method
        settings.input.property._h = self.property_type

        for key, route in sorted(self._accepted_keys.items()):
            if not self._has_value(key):
                continue

            value = self._get_value(key)
            metadata = route.metadata

            if route.scope == "top_level":
                if metadata.get("type") == "bool" and value == metadata.get("default"):
                    continue
                settings.input[key] = value
            elif route.scope == "property":
                settings.input.property[key] = value

        if not include_compounds:
            return settings

        compounds = self._compound_blocks()
        if compounds:
            settings.input.compound = compounds

        return settings

    def to_job(
            self,
            name: Optional[str] = None,
            **kwargs,
    ):
        if name is None:
            return CRSJob(settings=self.to_settings(), **kwargs)
        return CRSJob(settings=self.to_settings(), name=name, **kwargs)

    # General value helpers
    def _validate_input_key(self, key: str) -> None:
        if key not in self._accepted_keys:
            raise AttributeError(self._invalid_key_message(key))

    def _set_input_value(self, key: str, value: Any) -> None:
        if key == "mode":
            self._set_mode(value)
            return
        if key == "method":
            self._set_method(value)
            return

        self._validate_input_key(key)
        self._set_value(key, value)

    def _set_value(self, key: str, value: Any) -> None:
        object.__setattr__(self, key, value)

    def _get_value(self, key: str) -> Any:
        return getattr(self, key)

    def _has_value(self, key: str) -> bool:
        return hasattr(self, key)

    def _set_method(self, method: CRSMethodName) -> None:
        object.__setattr__(self, "_method", CRSJob._normalize_method(method))


    # mode helpers
    def _set_mode_or_default(self, mode: Optional[str]) -> None:
        selected_mode = self._mode_config.default if mode is None else mode
        if selected_mode is not None:
            self._set_mode(selected_mode)

    def _mode_values(self, mode: str) -> Dict[str, Any]:
        if not self._mode_config.options:
            raise ValueError(f"{self.property_type} does not support mode")

        normalized_mode = mode.lower()
        if normalized_mode not in self._mode_config.options:
            allowed = ", ".join(self._mode_config.options)
            raise ValueError(f"Unsupported mode {mode!r} for {self.property_type}. Supported modes: {allowed}")

        values: Dict[str, Any] = {}
        mode_values = self._mode_config.input_mapping.get(normalized_mode, {})
        for _, scoped_values in mode_values.items():
            values.update(scoped_values)
        return values

    def _mode_controlled_property_keys(self) -> Set[str]:
        keys = set(self._mode_config.options)

        for mode in self._mode_config.options:
            property_values = self._mode_config.input_mapping.get(mode, {}).get("property", {})
            keys.update(property_values)

        return keys

    def _set_mode(self, mode: str) -> None:
        normalized_mode = mode.lower()
        values = self._mode_values(normalized_mode)

        for key in self._active_mode_input_keys:
            if hasattr(self, key):
                object.__delattr__(self, key)

        for key, value in values.items():
            object.__setattr__(self, key, value)

        object.__setattr__(self, "_mode", normalized_mode)
        object.__setattr__(self, "_active_mode_input_keys", set(values))

    def _reject_mode_controlled_keys(self, values: Dict[str, Any]) -> None:
        if not self._mode_config.options:
            return

        # Mode-controlled flags, such as isobar, must be selected through mode=.
        provided = self._mode_controlled_property_keys().intersection(values)
        if not provided:
            return

        provided_keys = ", ".join(sorted(provided))
        allowed_modes = ", ".join(self._mode_config.options)

        raise ValueError(
            f"{self.property_type} mode-controlled key(s) must be set with mode=, "
            f"not keyword arguments: {provided_keys}. Supported modes: {allowed_modes}"
        )

    # compound helpers
    def _add_compound_role(
        self, role: str,
        compound: Union[Settings, PathLike],
        *,
        overrides: Optional[Dict[str, Any]] = None,
    ) -> "CRSInputBuilder":
        """Normalize and store a COMPOUND block under a builder compound role."""
        if role not in self._compound_roles:
            allowed = ", ".join(self._compound_roles) or "none"
            raise ValueError(f"{self.property_type} does not support {role} compounds. Supported roles: {allowed}")

        self._validate_compound_role_counts(check="max", role=role)

        if isinstance(compound, Settings):
            normalized = CRSJob._normalize_compound(compound.copy())
        else:
            normalized = CRSJob.compound_block(compound)

        self._set_compound_overrides(normalized, overrides or {})

        if role == "solute":
            self._drop_solute_composition(normalized)

        self._compounds_by_role[role].append(normalized)
        return self

    def _drop_solute_composition(self, compound: Settings) -> None:
        """Remove composition keys from solute blocks; only solvent/generic compounds carry composition."""
        dropped = []

        for key in ("frac1", "frac2"):
            if key in compound and compound[key] is not None:
                dropped.append(key)
                del compound[key]

        if dropped:
            log(
                f"Ignoring solute composition key(s) for {self.property_type}: "
                f"{', '.join(dropped)}. Set composition only on solvent or generic compound roles.",
                level=3,
            )

    def _compound_input_keys(self) -> Set[str]:
        return set(self.metadata["input_keys"].get("compound", ()))

    def _validate_compound_override_key(self, key: str) -> str:
        allowed = self._compound_input_keys()

        if key not in allowed:
            allowed_keys = ", ".join(sorted(allowed)) or "none"
            raise ValueError(
                f"{key} is not valid for {self.property_type} compounds. "
                f"Supported compound keys: {allowed_keys}"
            )

        return key

    def _set_compound_overrides(
        self,
        compound: Settings,
        overrides: Dict[str, Any],
    ) -> None:
        """Apply validated compound-key overrides such as frac1 and frac2."""
        for key, value in overrides.items():
            if value is None:
                continue

            compound[self._validate_compound_override_key(key)] = value

    def _compound_blocks(self) -> List[Settings]:
        compounds: List[Settings] = []
        for role in self._compound_roles:
            compounds.extend(self._compounds_by_role[role])
        return compounds

    def _has_compound_key(self, compound: Settings, key: str) -> bool:
        return key in compound and compound[key] is not None

    def _active_compound_required_keys(self) -> Dict[str, Dict[str, Any]]:
        if self._mode is None:
            return {}

        mode_metadata = self.metadata["builder"].get("mode", {})
        requirements = mode_metadata.get("compound_required_keys", {})
        return requirements.get(self._mode, {})


    # validation helpers
    def _validate_required_inputs(self) -> None:
        """Validate required top-level and PROPERTY input keys."""
        required_keys = set(self.metadata["required_keys"])

        missing = sorted(
            key for key in required_keys
            if key in self._accepted_keys and not self._has_value(key)
        )

        if missing:
            raise ValueError(f"{self.property_type} missing required input key(s): {', '.join(missing)}")

    def _validate_compound_role_counts(
        self,
        *,
        check: Literal["max", "bounds"],
        role: Optional[str] = None,
    ) -> None:
        """Validate property-specific compound role cardinality constraints."""
        role_config = self.metadata["builder"].get("role_config", {})
        roles = (role,) if role is not None else self._compound_roles

        for current_role in roles:
            config = role_config.get(current_role, {})
            count = len(self._compounds_by_role.get(current_role, ()))
            min_count = config.get("min_count")
            max_count = config.get("max_count")

            if check == "max":
                if max_count is not None and count >= max_count:
                    raise ValueError(
                        f"{self.property_type} supports at most {max_count} compound block(s) "
                        f"for role {current_role!r}"
                    )
                continue

            if min_count is not None and max_count == min_count and count != min_count:
                raise ValueError(
                    f"{self.property_type} requires exactly {min_count} compound block(s) "
                    f"for role {current_role!r}; got {count}"
                )

            if min_count is not None and count < min_count:
                raise ValueError(
                    f"{self.property_type} requires at least {min_count} compound block(s) "
                    f"for role {current_role!r}; got {count}"
                )

            if max_count is not None and count > max_count:
                raise ValueError(
                    f"{self.property_type} supports at most {max_count} compound block(s) "
                    f"for role {current_role!r}; got {count}"
                )

    def _validate_compound_required_keys(
        self,
        *,
        role: str,
        compound: Settings,
        requirements: Dict[str, Any],
        index: int,
    ) -> None:
        missing = [
            key for key in requirements.get("all_of", ())
            if not self._has_compound_key(compound, key)
        ]
        if missing:
            raise ValueError(
                f"{self.property_type} {role} #{index} missing required key(s) "
                f"for mode={self._mode!r}: {', '.join(missing)}"
            )

        any_of = requirements.get("any_of", ())
        if any_of and not any(
            all(self._has_compound_key(compound, key) for key in required_group)
            for required_group in any_of
        ):
            alternatives = " or ".join(
                " + ".join(required_group)
                for required_group in any_of
            )
            raise ValueError(
                f"{self.property_type} {role} #{index} requires {alternatives} "
                f"for mode={self._mode!r}"
            )

    def _validate_compound_required_keys_by_mode(self) -> None:
        """Validate mode-specific compound requirements such as solid-solute fusion data."""
        requirements_by_role = self._active_compound_required_keys()

        for role, requirements in requirements_by_role.items():
            for index, compound in enumerate(self._compounds_by_role.get(role, ()), start=1):
                self._validate_compound_required_keys(
                    role=role,
                    compound=compound,
                    requirements=requirements,
                    index=index,
                )


    # formatting/default helpers
    def _auto_default_candidate(self, key: str, route: CRSInputRoute) -> Tuple[bool, Any]:
        if key not in self._AUTO_DEFAULT_KEYS.get(route.scope, ()):
            return False, None

        if key == "pressure" and not getattr(self, "isobar", False):
            return False, None

        if "default" not in route.metadata:
            return False, None

        return True, self._get_default_value(route.metadata)

    @staticmethod
    def _get_default_value(metadata: Dict[str, Any]) -> Any:
        value = metadata.get("default")
        if isinstance(value, list) and len(value) == 1:
            return value[0]
        return value

    def _input_value_state(
        self,
        key: str,
        route: CRSInputRoute,
        include_values: bool,
    ) -> Tuple[Optional[str], Any]:
        if not include_values:
            return None, None

        if self._has_value(key):
            return "set", self._get_value(key)

        has_default, value = self._auto_default_candidate(key, route)
        if has_default:
            return "auto-default", value

        return "unset", None

    def _invalid_key_message(self, key: str) -> str:
        allowed_keys = set(self._accepted_keys)
        allowed_keys.add("method")

        if self._mode_config.options:
            allowed_keys.add("mode")

        allowed = ", ".join(sorted(allowed_keys))
        return f"{key!r} is not valid for {self.property_type}. Allowed input keys: {allowed}"

    @staticmethod
    def _format_input_description(
        key: str,
        scope: str,
        metadata: Dict[str, Any],
        *,
        value: Any = None,
        value_label: Optional[str] = None,
    ) -> str:
        input_type = metadata.get("type", "")
        unit = metadata.get("unit")
        type_unit = f"{input_type} [{unit}]" if unit else input_type
        description = metadata.get("description", "")

        line = f"{key} [{scope}] {type_unit}: {description}"

        if value_label is not None:
            line = f"{line} {value_label}: {value!r}"

        return line

    def _method_description(self, *, include_values: bool) -> str:
        return self._format_input_description(
            "method",
            "top_level",
            {"type": "multiple_choice", "description": "COSMO-RS/SAC method."},
            value=self._method,
            value_label="set" if include_values else None,
        )

    def _input_descriptions(self, *, include_values: bool) -> Tuple[str, ...]:
        excluded_keys = set(self._GLOBAL_TOP_LEVEL_KEYS)
        excluded_keys.update(self._mode_controlled_property_keys())

        descriptions = []
        for key, route in sorted(self._accepted_keys.items()):
            if key in excluded_keys:
                continue

            value_label, value = self._input_value_state(key, route, include_values)
            descriptions.append(
                self._format_input_description(
                    key,
                    route.scope,
                    route.metadata,
                    value=value,
                    value_label=value_label,
                )
            )

        return tuple(descriptions)

    def _compound_descriptions(self, *, include_details: bool) -> Tuple[str, ...]:
        compound_keys = tuple(self.metadata["input_keys"].get("compound", ()))
        if not compound_keys:
            return ()

        if not include_details:
            return (f"compound_keys [compound]: {', '.join(compound_keys)}",)

        compound_key_metadata = get_block_keys(CRS_DATA, "compound")
        return tuple(
            self._format_input_description(key, "compound", compound_key_metadata[key])
            for key in compound_keys
        )

    def _mode_description_text(self, mode: str) -> str:
        descriptions = self._mode_config.descriptions.get(mode, ())
        if isinstance(descriptions, str):
            return descriptions
        return " ".join(descriptions)

    def _mode_state_label(self, mode: str) -> str:
        labels = []
        if mode == self._mode:
            labels.append("current")
        if mode == self._mode_config.default:
            labels.append("default")
        return f" ({', '.join(labels)})" if labels else ""

    def _mode_descriptions(self) -> Tuple[str, ...]:
        if not self._mode_config.options:
            return ()

        return tuple(
            f"mode [{mode}]{self._mode_state_label(mode)}: {self._mode_description_text(mode)}"
            for mode in self._mode_config.options
        )

    def _comment_descriptions(self) -> Tuple[str, ...]:
        return tuple(f"comment: {comment}" for comment in self.metadata.get("comments", ()))

    def _input_hint_descriptions(self) -> Tuple[str, ...]:
        return CRSJob._property_type_input_hint_descriptions(self.metadata)


class _CompoundRoleMixin:
    def add_compound(
        self: _CRSInputBuilderT,
        compound: Union[Settings, PathLike],
        **kwargs: Any,
    ) -> _CRSInputBuilderT:
        """Add a COMPOUND block for property types using generic compounds."""
        return self._add_compound_role("compound", compound, overrides=kwargs)


class _SolventRoleMixin:
    def add_solvent(
        self: _CRSInputBuilderT,
        compound: Union[Settings, PathLike],
        **kwargs: Any,
    ) -> _CRSInputBuilderT:
        """Add a COMPOUND block as a solvent."""
        return self._add_compound_role("solvent", compound, overrides=kwargs)


class _SoluteRoleMixin:
    def add_solute(
        self: _CRSInputBuilderT,
        compound: Union[Settings, PathLike],
        **kwargs: Any,
    ) -> _CRSInputBuilderT:
        """Add a COMPOUND block as a solute."""
        return self._add_compound_role("solute", compound, overrides=kwargs)

class _SigmaPropertyMixin:
    nprofile: int
    sigmamax: float

class ACTIVITYCOEFInputBuilder(_SolventRoleMixin, _SoluteRoleMixin, CRSInputBuilder):
    """Builder for ACTIVITYCOEF CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "ACTIVITYCOEF"
    _DESCRIPTION: ClassVar[str] = "Activity coefficients in a solvent mixture."
    _SYSTEM_SCOPE: ClassVar[str] = "solvent_mixture_with_solutes"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "massfraction"),
        "property": ("densitysolvent",),
        "compound": ("frac1", "density") + _VAPOR_PRESSURE_KEYS,
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "frac1")
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("densitysolvent", "density")
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("solvent", "solute")


class LOGPInputBuilder(_SolventRoleMixin, _SoluteRoleMixin, CRSInputBuilder):
    """Builder for LOGP CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "LOGP"
    _DESCRIPTION: ClassVar[str] = "Partition coefficients between two immiscible solvent phases."
    _SYSTEM_SCOPE: ClassVar[str] = "mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "massfraction"),
        "property": ("volumequotient",),
        "compound": ("frac1", "frac2", "density"),
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "frac1", "frac2")
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("volumequotient",)
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("solvent", "solute")
    _COMPOUND_ROLE_CONFIG: ClassVar[Dict[str, Dict[str, Optional[int]]]] = {
        "solvent": {"min_count": 2},
        "solute": {"min_count": 1},
    }


class _SolubilityModeMixin:
    @property
    def mode(self) -> Optional[CRSSolubilityMode]:
        return cast(Optional[CRSSolubilityMode], self._mode)

    @mode.setter
    def mode(self, mode: CRSSolubilityMode) -> None:
        self._set_mode(mode)

    @property
    def mode_options(self) -> Tuple[CRSSolubilityMode, ...]:
        return cast(Tuple[CRSSolubilityMode, ...], self._mode_config.options)


class _VLESweepModeMixin:
    @property
    def mode(self) -> Optional[CRSVLESweepMode]:
        return cast(Optional[CRSVLESweepMode], self._mode)

    @mode.setter
    def mode(self, mode: CRSVLESweepMode) -> None:
        self._set_mode(mode)

    @property
    def mode_options(self) -> Tuple[CRSVLESweepMode, ...]:
        return cast(Tuple[CRSVLESweepMode, ...], self._mode_config.options)


class SOLUBILITYInputBuilder(_SolubilityModeMixin, _SolventRoleMixin, _SoluteRoleMixin, CRSInputBuilder):
    """Builder for SOLUBILITY CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "SOLUBILITY"
    _DESCRIPTION: ClassVar[str] = "Solubility of solutes in a solvent mixture or under gas-pressure conditions."
    _SYSTEM_SCOPE: ClassVar[str] = "solvent_with_solutes"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "pressure", "massfraction"),
        "property": ("densitysolvent", "isobar"),
        "compound": ("frac1", "density") + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS,
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "frac1")
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("densitysolvent", "density") + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("solvent", "solute")
    _MODE_CONFIG: ClassVar[Dict[str, Any]] = _SOLUBILITY_MODE_CONFIG


class PURESOLUBILITYInputBuilder(_SolubilityModeMixin, _SolventRoleMixin, _SoluteRoleMixin, CRSInputBuilder):
    """Builder for PURESOLUBILITY CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "PURESOLUBILITY"
    _DESCRIPTION: ClassVar[str] = "Solubility of a solute in pure solvents over a temperature range."
    _SYSTEM_SCOPE: ClassVar[str] = "pure_solvent_with_solute"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "pressure"),
        "property": ("isobar",),
        "compound": ("frac1", "density") + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS,
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "frac1")
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("density",) + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("solvent", "solute")
    _COMPOUND_ROLE_CONFIG: ClassVar[Dict[str, Dict[str, Optional[int]]]] = {
        "solvent": {"min_count": 1},
        "solute": {"min_count": 1, "max_count": 1},
    }
    _MODE_CONFIG: ClassVar[Dict[str, Any]] = _SOLUBILITY_MODE_CONFIG


class VAPORPRESSUREInputBuilder(_CompoundRoleMixin, CRSInputBuilder):
    """Builder for VAPORPRESSURE CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "VAPORPRESSURE"
    _DESCRIPTION: ClassVar[str] = "Vapor pressure of a mixture at fixed temperature."
    _SYSTEM_SCOPE: ClassVar[str] = "mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "massfraction"),
        "compound": ("frac1",) + _VAPOR_PRESSURE_KEYS,
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "frac1")
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",) + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class PUREVAPORPRESSUREInputBuilder(_CompoundRoleMixin, CRSInputBuilder):
    """Builder for PUREVAPORPRESSURE CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "PUREVAPORPRESSURE"
    _DESCRIPTION: ClassVar[str] = "Pure-compound vapor pressure over a temperature range."
    _SYSTEM_SCOPE: ClassVar[str] = "pure_compounds"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature",),
        "compound": _VAPOR_PRESSURE_KEYS,
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",) + _VAPOR_PRESSURE_KEYS
    _COMMENTS: ClassVar[Tuple[str, ...]] = _PURE_COMPOUND_COMMENTS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class BOILINGPOINTInputBuilder(_CompoundRoleMixin, CRSInputBuilder):
    """Builder for BOILINGPOINT CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "BOILINGPOINT"
    _DESCRIPTION: ClassVar[str] = "Boiling temperature of a mixture for a pressure range."
    _SYSTEM_SCOPE: ClassVar[str] = "mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("pressure", "massfraction"),
        "compound": ("frac1",) + _VAPOR_PRESSURE_KEYS,
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("pressure",)
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("pressure",) + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class PUREBOILINGPOINTInputBuilder(_CompoundRoleMixin, CRSInputBuilder):
    """Builder for PUREBOILINGPOINT CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "PUREBOILINGPOINT"
    _DESCRIPTION: ClassVar[str] = "Pure-compound boiling point over a pressure range."
    _SYSTEM_SCOPE: ClassVar[str] = "pure_compounds"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("pressure",),
        "compound": _VAPOR_PRESSURE_KEYS,
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("pressure",)
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("pressure",) + _VAPOR_PRESSURE_KEYS
    _COMMENTS: ClassVar[Tuple[str, ...]] = _PURE_COMPOUND_COMMENTS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class FLASHPOINTInputBuilder(_CompoundRoleMixin, CRSInputBuilder):
    """Builder for FLASHPOINT CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "FLASHPOINT"
    _DESCRIPTION: ClassVar[str] = "Flash point of a mixture using user-supplied pure-compound flash points."
    _SYSTEM_SCOPE: ClassVar[str] = "mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("massfraction",),
        "compound": ("frac1", "flashpoint") + _VAPOR_PRESSURE_KEYS,
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",)
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("flashpoint",) + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class BINMIXCOEFInputBuilder(_VLESweepModeMixin, _CompoundRoleMixin, CRSInputBuilder):
    """Builder for BINMIXCOEF CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "BINMIXCOEF"
    _DESCRIPTION: ClassVar[str] = "Binary-mixture coefficients over a composition range."
    _SYSTEM_SCOPE: ClassVar[str] = "binary_mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "pressure", "massfraction"),
        "property": _VLE_SWEEP_PROPERTY_KEYS,
        "compound": ("frac1",) + _VAPOR_PRESSURE_KEYS + ("flashpoint",),
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("flashpoint",) + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Dict[str, Dict[str, Optional[int]]]] = {
        "compound": {"min_count": 2, "max_count": 2},
    }
    _MODE_CONFIG: ClassVar[Dict[str, Any]] = _VLE_SWEEP_MODE_CONFIG


class TERNARYMIXInputBuilder(_VLESweepModeMixin, _CompoundRoleMixin, CRSInputBuilder):
    """Builder for TERNARYMIX CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "TERNARYMIX"
    _DESCRIPTION: ClassVar[str] = "Ternary mixture property sweep over composition space."
    _SYSTEM_SCOPE: ClassVar[str] = "ternary_mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "pressure", "massfraction"),
        "property": _VLE_SWEEP_PROPERTY_KEYS,
        "compound": ("frac1",) + _VAPOR_PRESSURE_KEYS + ("flashpoint",),
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("flashpoint",) + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Dict[str, Dict[str, Optional[int]]]] = {
        "compound": {"min_count": 3, "max_count": 3},
    }
    _MODE_CONFIG: ClassVar[Dict[str, Any]] = _VLE_SWEEP_MODE_CONFIG


class COMPOSITIONLINEInputBuilder(_VLESweepModeMixin, _SolventRoleMixin, CRSInputBuilder):
    """Builder for COMPOSITIONLINE CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "COMPOSITIONLINE"
    _DESCRIPTION: ClassVar[str] = "Composition-line calculation between two endpoint phase compositions."
    _SYSTEM_SCOPE: ClassVar[str] = "binary_mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "pressure", "massfraction"),
        "property": _VLE_SWEEP_PROPERTY_KEYS,
        "compound": ("frac1", "frac2") + _VAPOR_PRESSURE_KEYS + ("flashpoint",),
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _HINT_KEYS: ClassVar[Tuple[str, ...]] = ("flashpoint",) + _VAPOR_PRESSURE_KEYS
    _COMMENTS: ClassVar[Tuple[str, ...]] = (
        "frac1 and frac2 define two endpoint solutions mixed along the composition line.",
    )
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("solvent",)
    _MODE_CONFIG: ClassVar[Dict[str, Any]] = _VLE_SWEEP_MODE_CONFIG


class LLEInputBuilder(_CompoundRoleMixin, CRSInputBuilder):
    """Builder for LLE CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "LLE"
    _DESCRIPTION: ClassVar[str] = "Liquid-liquid equilibrium for a ternary mixture."
    _SYSTEM_SCOPE: ClassVar[str] = "mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "massfraction"),
        "compound": ("frac1",),
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "frac1")
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class STABILITYInputBuilder(_CompoundRoleMixin, CRSInputBuilder):
    """Builder for STABILITY CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "STABILITY"
    _DESCRIPTION: ClassVar[str] = "Michelsen tangent-plane-distance stability test for a feed composition."
    _SYSTEM_SCOPE: ClassVar[str] = "mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "massfraction"),
        "compound": ("frac1",),
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "frac1")
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class SIGMAPROFILEInputBuilder(_SigmaPropertyMixin, _CompoundRoleMixin, CRSInputBuilder):
    """Builder for SIGMAPROFILE CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "SIGMAPROFILE"
    _DESCRIPTION: ClassVar[str] = "Sigma profile for a solvent mixture."
    _SYSTEM_SCOPE: ClassVar[str] = "mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("massfraction",),
        "property": _SIGMA_PROPERTY_KEYS,
        "compound": ("frac1",),
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",)
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class PURESIGMAPROFILEInputBuilder(_SigmaPropertyMixin, _CompoundRoleMixin, CRSInputBuilder):
    """Builder for PURESIGMAPROFILE CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "PURESIGMAPROFILE"
    _DESCRIPTION: ClassVar[str] = "Sigma profile for pure compounds."
    _SYSTEM_SCOPE: ClassVar[str] = "pure_compounds"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "property": _SIGMA_PROPERTY_KEYS,
    }
    _COMMENTS: ClassVar[Tuple[str, ...]] = _PURE_COMPOUND_COMMENTS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class SIGMAPOTENTIALInputBuilder(_SigmaPropertyMixin, _CompoundRoleMixin, CRSInputBuilder):
    """Builder for SIGMAPOTENTIAL CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "SIGMAPOTENTIAL"
    _DESCRIPTION: ClassVar[str] = "Sigma potential for a solvent mixture."
    _SYSTEM_SCOPE: ClassVar[str] = "mixture"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "top_level": ("temperature", "massfraction"),
        "property": _SIGMA_PROPERTY_KEYS,
        "compound": ("frac1",),
    }
    _REQUIRED_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "frac1")
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


class PURESIGMAPOTENTIALInputBuilder(_SigmaPropertyMixin, _CompoundRoleMixin, CRSInputBuilder):
    """Builder for PURESIGMAPOTENTIAL CRS input settings."""

    _PROPERTY_TYPE: ClassVar[str] = "PURESIGMAPOTENTIAL"
    _DESCRIPTION: ClassVar[str] = "Sigma potential for pure compounds."
    _SYSTEM_SCOPE: ClassVar[str] = "pure_compounds"
    _INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
        "property": _SIGMA_PROPERTY_KEYS,
    }
    _COMMENTS: ClassVar[Tuple[str, ...]] = _PURE_COMPOUND_COMMENTS
    _COMPOUND_ROLES: ClassVar[Tuple[str, ...]] = ("compound",)


_CRS_INPUT_BUILDER_CLASSES: Dict[str, type[CRSInputBuilder]] = {
    cls._PROPERTY_TYPE: cls
    for cls in (
        ACTIVITYCOEFInputBuilder,
        LOGPInputBuilder,
        SOLUBILITYInputBuilder,
        PURESOLUBILITYInputBuilder,
        VAPORPRESSUREInputBuilder,
        PUREVAPORPRESSUREInputBuilder,
        BOILINGPOINTInputBuilder,
        PUREBOILINGPOINTInputBuilder,
        FLASHPOINTInputBuilder,
        BINMIXCOEFInputBuilder,
        TERNARYMIXInputBuilder,
        COMPOSITIONLINEInputBuilder,
        LLEInputBuilder,
        STABILITYInputBuilder,
        SIGMAPROFILEInputBuilder,
        PURESIGMAPROFILEInputBuilder,
        SIGMAPOTENTIALInputBuilder,
        PURESIGMAPOTENTIALInputBuilder,
    )
}

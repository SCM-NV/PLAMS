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


__all__ = ["CRSResults", "CRSJob", "CRSResultTables", "CRSInputBuilder"]

PathLike = Union[str, os.PathLike]


class CRSResultTables(NamedTuple):
    """Container returned by :meth:`CRSResults.get_result_table` with ``split=True``."""

    component: "pd.DataFrame"
    mixture: Optional["pd.DataFrame"]
    lle: Optional["pd.DataFrame"]


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

        return pd.DataFrame(
            rows, columns=list(self._RESULT_TABLE_COMPONENT_BASE_COLUMNS) + ["molmass"] + list(quantities)
        )

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
    _PROPERTY_TYPE_METADATA = crs_defs.CRS_PROPERTY_TYPE_METADATA

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
    def property_types() -> Tuple[str, ...]:
        """Return the supported COSMO-RS problem types."""
        return tuple(CRSJob._PROPERTY_TYPE_METADATA)

    @staticmethod
    def property_type_metadata(
        property_type: str,
        as_summary: bool = True,
    ) -> Union[Dict[str, Any], Tuple[str, ...]]:
        """Return discoverability metadata for a COSMO-RS problem type."""
        normalized = CRSJob._normalize_property_type(property_type)
        data = CRSJob._PROPERTY_TYPE_METADATA[normalized]
        metadata = {key: value.copy() if isinstance(value, list) else value for key, value in data.items()}
        if not as_summary:
            return deepcopy(metadata)

        summary = [
            f"{normalized}: {metadata['description']}",
            f"system_scope: {metadata['system_scope']}",
            f"top_level_keys: {', '.join(metadata['input_keys']['top_level']) or '-'}",
            f"property_keys: {', '.join(metadata['input_keys']['property']) or '-'}",
            f"compound_keys: {', '.join(metadata['input_keys']['compound']) or '-'}",
            f"required_keys: {', '.join(metadata['required_keys']) or '-'}",
        ]
        for note in metadata.get("notes", []):
            summary.append(f"note: {note}")

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
            supported_values=CRSJob.methods(),
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
        if not Path(normalized_path).is_file():
            raise FileNotFoundError(f"{method_name} path does not exist: {normalized_path}")
        return normalized_path

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
    def compound_from_database(
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
            raise ValueError(f"compound_from_database does not accept override(s): {invalid}")

        normalized_property_type = None
        if property_type is not None:
            normalized_property_type = CRSJob._normalize_property_type(property_type)

        if property_keys is None:
            if normalized_property_type is None:
                requested_property_keys: Tuple[str, ...] = ()
            else:
                requested_property_keys = tuple(
                    key
                    for key in CRSJob._PROPERTY_TYPE_METADATA[normalized_property_type]["compound_keys"]
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
    def input_builder(
        property_type: str,
        *,
        method: str = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        **kwargs: Any,
    ) -> "CRSInputBuilder":
        """Return a property-type-aware builder for CRS input settings.
        ``method`` selects the COSMO-RS/SAC method.
        ``mode`` is supported only for selected property types.

        Example::

            crs_input = CRSJob.input_builder("SOLUBILITY", temperature="353.15 373.15 10", mode="gas")

            benzene = CRSJob.compound_block(CRSJob.coskf_from_database("Benzene.coskf"))
            benzene.pvap = 353.3
            benzene.tvap = 1.01325

            water = CRSJob.compound_block(CRSJob.coskf_from_database("Water.coskf"), frac1=1.0)
            water.vp_equation = "Antoine"
            water.vp_params = "5.40221 1838.675 -31.737 0.0 0.0"

            crs_input.add_solvent(water)
            crs_input.add_solute(benzene)

            settings = crs_input.to_settings()

        Use ``builder.describe()`` and ``CRSJob.property_type_metadata(...)`` for available keys.
        """
        return CRSInputBuilder(property_type)._configure(
            method=method,
            mode=mode,
            use_defaults=use_defaults,
            values=kwargs
        )


class CRSInputBuilder:
    """Property-type builder for CRSJob input settings."""

    _DEFAULT_INPUT_KEYS: ClassVar[Dict[str, Tuple[str, ...]]] = {
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
        "get",
        "describe",
        "apply_defaults",
        "to_settings",
    )

    _INTERNAL_ATTRIBUTES: ClassVar[Set[str]] = {
        "property_type",
        "metadata",
        "_method",
        "_top_level_keys",
        "_property_keys",
        "_flat_keys",
        "_mode",
        "_has_mode",
        "_mode_options",
        "_mode_default",
        "_mode_input_mapping",
        "_active_mode_input_keys",
        "_compound_roles",
        "_compounds_by_role",
    }

    def __init__(self, property_type: str) -> None:
        normalized = CRSJob._normalize_property_type(property_type)
        metadata = CRSJob._PROPERTY_TYPE_METADATA[normalized]
        mode_metadata = metadata["builder"].get("mode")
        builder_metadata = metadata["builder"]

        top_level_keys = set(metadata["input_keys"]["top_level"])
        property_keys = set(metadata["input_keys"]["property"])

        mode_options = tuple(mode_metadata.get("keys", ())) if mode_metadata else ()
        mode_default = mode_metadata.get("default") if mode_metadata else None
        mode_input_mapping = mode_metadata.get("set", {}) if mode_metadata else {}

        compound_roles = tuple(builder_metadata.get("roles", ()))
        compounds_by_role = {role: [] for role in compound_roles}

        object.__setattr__(self, "property_type", normalized)
        object.__setattr__(self, "metadata", metadata)
        object.__setattr__(self, "_top_level_keys", top_level_keys)
        object.__setattr__(self, "_property_keys", property_keys)
        object.__setattr__(self, "_flat_keys", top_level_keys | property_keys)
        object.__setattr__(self, "_method", CRSJob._normalize_method("COSMO-RS"))

        # User modes can set CRS input keys; for example, mode="gas" sets isobar=True.
        object.__setattr__(self, "_has_mode", mode_metadata is not None)
        object.__setattr__(self, "_mode", None)
        object.__setattr__(self, "_active_mode_input_keys", set())
        object.__setattr__(self, "_mode_options", mode_options)
        object.__setattr__(self, "_mode_default", mode_default)
        object.__setattr__(self, "_mode_input_mapping", mode_input_mapping)

        object.__setattr__(self, "_compound_roles", compound_roles)
        object.__setattr__(self, "_compounds_by_role", compounds_by_role)

    def __dir__(self) -> List[str]:
        """Return property-type-specific completions for interactive use."""
        entries = set()
        entries.update(self._flat_keys)
        entries.update(self._GLOBAL_TOP_LEVEL_KEYS)
        entries.update(self._COMMON_DIR_ENTRIES)

        if self._has_mode:
            entries.difference_update(self._mode_options)
            entries.difference_update(self._active_mode_input_keys)
            entries.update({"mode", "mode_options"})

        if "compound" in self._compound_roles:
            entries.add("add_compound")
        if "solvent" in self._compound_roles:
            entries.add("add_solvent")
        if "solute" in self._compound_roles:
            entries.add("add_solute")

        return sorted(entries)

    def __setattr__(self, key: str, value: Any) -> None:
        if key in self._INTERNAL_ATTRIBUTES:
            object.__setattr__(self, key, value)
            return

        if key in self._GLOBAL_TOP_LEVEL_KEYS:
            object.__setattr__(self, key, value)
            return

        if key in self._flat_keys:
            object.__setattr__(self, key, value)
            return

        if key == "mode":
            self._set_mode(value)
            return

        if key == "method":
            self._set_method(value)
            return

        raise AttributeError(self._invalid_key_message(key))

    @property
    def mode(self) -> Optional[str]:
        return self._mode

    @property
    def mode_options(self) -> Tuple[str, ...]:
        return self._mode_options

    def get(self, key: str, default: Any = None) -> Any:
        """Return an input value if set, otherwise return default."""
        if key == "method":
            return self._method
        if key == "method_options":
            return CRSJob.methods()
        if key == "mode":
            return self._mode
        if key == "mode_options":
            return self._mode_options
        if key in self._GLOBAL_TOP_LEVEL_KEYS:
            return getattr(self, key, default)
        if key not in self._flat_keys:
            raise AttributeError(self._invalid_key_message(key))
        return getattr(self, key, default)

    def add_compound(self, compound: Settings, *, frac1: Optional[float] = None) -> "CRSInputBuilder":
        """Add a COMPOUND block for property types using generic compounds."""
        return self._add_compound_role("compound", compound, overrides={"frac1": frac1})

    def add_solvent(
        self,
        compound: Settings,
        *,
        frac1: Optional[float] = None,
        frac2: Optional[float] = None,
    ) -> "CRSInputBuilder":
        """Add a COMPOUND block as a solvent."""
        return self._add_compound_role("solvent", compound, overrides={"frac1": frac1, "frac2": frac2})

    def add_solute(self, compound: Settings) -> "CRSInputBuilder":
        """Add a COMPOUND block as a solute."""
        return self._add_compound_role("solute", compound)

    def apply_defaults(self) -> "CRSInputBuilder":
        """Set default input values for unset supported keys."""
        for key in self._DEFAULT_INPUT_KEYS["top_level"]:
            if key not in self._top_level_keys or hasattr(self, key):
                continue
            if key == "pressure" and not getattr(self, "isobar", False):
                continue
            object.__setattr__(self, key, self._get_default_value(CRS_DATA[key]))

        property_key_metadata = get_block_keys(CRS_DATA, "property")
        for key in self._DEFAULT_INPUT_KEYS["property"]:
            if key not in self._property_keys or hasattr(self, key):
                continue
            object.__setattr__(self, key, self._get_default_value(property_key_metadata[key]))

        return self

    def describe(
        self,
        include_hints: bool = False,
        include_values: bool = False,
        include_compound_details: bool = False,
    ) -> Tuple[str, ...]:
        """Return property-type-specific input keys and mode descriptions."""
        descriptions: List[str] = [
            self._format_input_description(
                "method",
                "top_level",
                {"type": "multiple_choice", "description": "COSMO-RS/SAC method."},
                value=self._method,
                value_label="set" if include_values else None,
            )
        ]

        for key in sorted(self._top_level_keys):
            metadata = CRS_DATA[key]
            value_label, value = self._input_value_state(key, "top_level", metadata, include_values)
            descriptions.append(
                self._format_input_description(
                    key,
                    "top_level",
                    metadata,
                    value=value,
                    value_label=value_label,
                )
            )

        property_key_metadata = get_block_keys(CRS_DATA, "property")
        for key in sorted(self._visible_property_keys()):
            metadata = property_key_metadata[key]
            value_label, value = self._input_value_state(key, "property", metadata, include_values)

            descriptions.append(
                self._format_input_description(
                    key,
                    "property",
                    metadata,
                    value=value,
                    value_label=value_label,
                )
            )

        if include_compound_details:
            descriptions.extend(self._compound_key_detail_descriptions())
        else:
            compound_keys = self._compound_keys_description()
            if compound_keys is not None:
                descriptions.append(compound_keys)

        mode_description = self._mode_description()
        if mode_description is not None:
            descriptions.append(mode_description)
            descriptions.extend(self._mode_hint_descriptions())

        if include_hints:
            descriptions.extend(self._note_descriptions())
            descriptions.extend(self._input_hint_descriptions())

        return tuple(descriptions)

    def to_settings(self, *, include_compounds: bool = True) -> Settings:
        """Build PLAMS Settings, optionally omitting COMPOUND blocks."""
        self._validate_required_inputs()
        if include_compounds:
            self._validate_compound_required_keys_by_mode()

        settings = Settings()
        settings.input.method = self._method
        settings.input.property._h = self.property_type

        for key in self._top_level_keys:
            if not hasattr(self, key):
                continue

            value = getattr(self, key)
            metadata = CRS_DATA[key]
            if metadata.get("type") == "bool" and value == metadata.get("default"):
                continue

            settings.input[key] = value


        for key in self._GLOBAL_TOP_LEVEL_KEYS:
            if not hasattr(self, key):
                continue

            value = getattr(self, key)
            metadata = CRS_DATA[key]
            if value == metadata.get("default"):
                continue

            settings.input[key] = value

        for key in self._property_keys | self._active_mode_input_keys:
            if hasattr(self, key):
                settings.input.property[key] = getattr(self, key)

        if not include_compounds:
            return settings

        compounds = self._compound_blocks()
        if compounds:
            settings.input.compound = compounds

        return settings

    # compound helpers
    def _add_compound_role(
        self, role: str,
        compound: Settings,
        *,
        overrides: Optional[Dict[str, Any]] = None,
    ) -> "CRSInputBuilder":
        """Normalize and store a COMPOUND block under a builder compound role."""
        if role not in self._compound_roles:
            allowed = ", ".join(self._compound_roles) or "none"
            raise ValueError(f"{self.property_type} does not support {role} compounds. Supported roles: {allowed}")

        normalized = CRSJob._normalize_compound(compound.copy())
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

    def _set_compound_overrides(
        self,
        compound: Settings,
        overrides: Dict[str, Any],
    ) -> None:
        """Apply validated compound-key overrides such as frac1 and frac2."""
        allowed = self._compound_input_keys()

        for key, value in overrides.items():
            if value is None:
                continue

            if key not in allowed:
                allowed_keys = ", ".join(sorted(allowed)) or "none"
                raise ValueError(
                    f"{key} is not valid for {self.property_type} compounds. "
                    f"Supported compound keys: {allowed_keys}"
                )

            compound[key] = value

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

    # validation helpers
    def _validate_required_inputs(self) -> None:
        """Validate required top-level and PROPERTY input keys."""
        required_keys = set(self.metadata["required_keys"])

        missing = sorted(
            key for key in required_keys
            if (
                (key in self._top_level_keys or key in self._property_keys)
                and not hasattr(self, key)
            )
        )

        if missing:
            raise ValueError(f"{self.property_type} missing required input key(s): {', '.join(missing)}")

    # mode/config helpers
    def _mode_values(self, mode: str) -> Dict[str, Any]:
        if not self._has_mode:
            raise ValueError(f"{self.property_type} does not support mode")

        normalized_mode = mode.lower()
        if normalized_mode not in self._mode_options:
            allowed = ", ".join(self._mode_options)
            raise ValueError(f"Unsupported mode {mode!r} for {self.property_type}. Supported modes: {allowed}")

        values: Dict[str, Any] = {}
        mode_values = self._mode_input_mapping.get(normalized_mode, {})
        for _, scoped_values in mode_values.items():
            values.update(scoped_values)
        return values

    def _mode_controlled_property_keys(self) -> Set[str]:
        keys = set(self._mode_options)

        for mode in self._mode_options:
            property_values = self._mode_input_mapping.get(mode, {}).get("property", {})
            keys.update(property_values)

        return keys

    def _visible_property_keys(self) -> Set[str]:
        return self._property_keys.difference(self._mode_controlled_property_keys())

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

    def _reject_mode_keys(self, values: Dict[str, Any]) -> None:
        if not self._has_mode:
            return
        provided = set(self._mode_options).intersection(values)
        if not provided:
            return

        provided_keys = ", ".join(sorted(provided))
        allowed = ", ".join(self._mode_options)
        raise ValueError(
            f"{self.property_type} mode keys must be set with mode=, not keyword arguments: "
            f"{provided_keys}. Supported modes: {allowed}"
        )

    @property
    def method(self) -> str:
        return self._method

    @property
    def method_options(self) -> Tuple[str, ...]:
        return CRSJob.methods()

    def _set_method(self, method: str) -> None:
        object.__setattr__(self, "_method", CRSJob._normalize_method(method))

    def _configure(
        self,
        *,
        method: str = "COSMO-RS",
        mode: Optional[str] = None,
        use_defaults: bool = False,
        values: Dict[str, Any],
    ) -> "CRSInputBuilder":
        self._reject_mode_keys(values)
        self._set_method(method)

        if mode is None:
            mode = self._mode_default

        if mode is not None:
            self._set_mode(mode)

        for key, value in values.items():
            setattr(self, key, value)

        if use_defaults:
            self.apply_defaults()

        return self

    # formatting/default helpers
    @staticmethod
    def _get_default_value(metadata: Dict[str, Any]) -> Any:
        value = metadata.get("default")
        if isinstance(value, list) and len(value) == 1:
            return value[0]
        return value

    def _auto_default_value(self, key: str, scope: str, metadata: Dict[str, Any]) -> Tuple[bool, Any]:
        if key not in self._DEFAULT_INPUT_KEYS.get(scope, ()):
            return False, None
        if key == "pressure" and not getattr(self, "isobar", False):
            return False, None
        if "default" not in metadata:
            return False, None
        return True, self._get_default_value(metadata)

    def _input_value_state(
        self,
        key: str,
        scope: str,
        metadata: Dict[str, Any],
        include_values: bool,
    ) -> Tuple[Optional[str], Any]:
        if not include_values:
            return None, None

        if hasattr(self, key):
            return "set", getattr(self, key)

        has_auto_default, auto_default = self._auto_default_value(key, scope, metadata)
        if has_auto_default:
            return "auto-default", auto_default

        return "unset", None

    def _invalid_key_message(self, key: str) -> str:
        allowed_keys = set(self._flat_keys)
        allowed_keys.update(self._GLOBAL_TOP_LEVEL_KEYS)
        allowed_keys.add("method")

        if self._has_mode:
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
        type_unit = metadata.get("type", "")
        if metadata.get("unit"):
            type_unit = f"{type_unit} [{metadata['unit']}]"

        line = f"{key} [{scope}] {type_unit}: {metadata['description']}"

        if value_label is not None:
            line = f"{line} {value_label}: {value!r}"

        return line


    def _mode_description(self) -> Optional[str]:
        if not self._has_mode:
            return None

        options = []
        for mode in self._mode_options:
            label = f"{mode} (*)" if mode == self._mode else mode
            options.append(label)

        option_text = ", ".join(options)

        if self._mode_default is not None:
            return f"mode [mode]: {option_text} (default: {self._mode_default})"
        return f"mode [mode]: {option_text}"

    def _mode_hint_descriptions(self) -> Tuple[str, ...]:
        if not self._has_mode:
            return ()

        mode_hints = self.metadata["builder"].get("mode", {}).get("mode_hints", {})
        lines = []

        for mode in self._mode_options:
            hints = mode_hints.get(mode, ())
            if isinstance(hints, str):
                hints = (hints,)

            mode_label = f"{mode} (*)" if mode == self._mode else mode
            for hint in hints:
                lines.append(f"mode_hint [{mode_label}]: {hint}")

        return tuple(lines)

    def _compound_key_detail_descriptions(self) -> Tuple[str, ...]:
        compound_keys = self.metadata["input_keys"].get("compound", ())
        if not compound_keys:
            return ()

        compound_key_metadata = get_block_keys(CRS_DATA, "compound")
        return tuple(
            self._format_input_description(key, "compound", compound_key_metadata[key])
            for key in compound_keys
        )

    def _compound_keys_description(self) -> Optional[str]:
        compound_keys = self.metadata["input_keys"].get("compound", ())
        if not compound_keys:
            return None
        return f"compound_keys [compound]: {', '.join(compound_keys)}"

    def _input_hint_descriptions(self) -> Tuple[str, ...]:
        return CRSJob._property_type_input_hint_descriptions(self.metadata)

    def _note_descriptions(self) -> Tuple[str, ...]:
        return tuple(f"note: {note}" for note in self.metadata.get("notes", ()))

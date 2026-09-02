from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Callable, ClassVar, Dict, List, Literal, Mapping, Optional, Sequence, Set, Tuple, Type, TypeVar, Union

from scm.plams.core.functions import log
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.crs_definitions import CRS_METHODS, get_block_keys, get_crs_input_data


__all__ = [
    "ACTIVITYCOEFInputBuilder",
    "BINMIXCOEFInputBuilder",
    "BOILINGPOINTInputBuilder",
    "COMPOSITIONLINEInputBuilder",
    "CRSInputBuilder",
    "FLASHPOINTInputBuilder",
    "LLEInputBuilder",
    "LOGPInputBuilder",
    "PUREBOILINGPOINTInputBuilder",
    "PURESIGMAPOTENTIALInputBuilder",
    "PURESIGMAPROFILEInputBuilder",
    "PURESOLUBILITYInputBuilder",
    "PUREVAPORPRESSUREInputBuilder",
    "SIGMAPOTENTIALInputBuilder",
    "SIGMAPROFILEInputBuilder",
    "SOLUBILITYInputBuilder",
    "STABILITYInputBuilder",
    "TERNARYMIXInputBuilder",
    "VAPORPRESSUREInputBuilder",
    "input_builder",
    "methods",
]

CRSMethodName = Literal[
    "COSMO-RS",
    "COSMOSAC2013",
    "COSMOSAC2016",
    "COSMOSACDHB",
    "COSMOSACDHB-MESP",
]

_PathLike = Union[str, os.PathLike]
_SolubilityMode = Literal["solid", "liquid", "gas"]
_VLESweepMode = Literal["isotherm", "isobar", "flashpoint"]
_FloatInput = Union[float, int]
_FloatListInput = Union[_FloatInput, str, Sequence[_FloatInput]]
_IntegerListInput = Union[int, str, Sequence[int]]

_CRSInputBuilderT = TypeVar("_CRSInputBuilderT", bound="CRSInputBuilder")

_METHOD_ALIASES = {
    "COSMORS": "COSMO-RS",
    "COSMOSAC": "COSMOSAC2013",
}

_VAPOR_PRESSURE_KEYS = ("pvap", "tvap", "vp_equation", "vp_params")
_FUSION_KEYS = ("meltingpoint", "hfusion", "cpfusion")
_VLE_SWEEP_PROPERTY_KEYS = ("nfrac", "isotherm", "isobar", "flashpoint")
_SIGMA_KEYS = ("nprofile", "sigmamax")
_SIGMA_MOMENT_KEYS = (
    "sigmamomentpower",
    "sigmamomenthblevel",
    "sigmamomenthbcutoff",
    "sigmamomenthbcutoffbase",
    "sigmamomenthbcutoffstep"
)
# Technical top-level CRS inputs exposed through set_expert_options() for all properties.
# The keys and value types are validated from crs.json; only the expert-option
# grouping is kept here for now and may move to crs.json metadata later.
_EXPERT_INPUT_KEYS = (
    "usepolycombiforpolymer",
    "pdh_correction",
    "ignore_coskfatoms",
    "output_energy_components",
    "reuse_sigma_profile",
    "update_sigma_profile",
)


@dataclass(frozen=True)
class _InputRoute:
    """Route from a builder key to a PLAMS Settings location."""

    scope: str
    metadata: Mapping[str, Any]


@dataclass(frozen=True)
class _CompoundRoleConfig:
    """Compound role count limits and required keys."""

    min_count: Optional[int] = None
    max_count: Optional[int] = None
    required_keys: Tuple[str, ...] = ()


@dataclass(frozen=True)
class _ModeOptionConfig:
    """Input values enabled by one property mode."""

    input_values: Mapping[str, Any] = field(default_factory=dict)
    required_input_keys: Tuple[str, ...] = ()


@dataclass(frozen=True)
class _ModeConfig:
    """Mode-controlled input-key configuration."""

    default: Optional[str] = None
    options: Mapping[str, _ModeOptionConfig] = field(default_factory=dict)


_SOLUBILITY_MODE_CONFIG = _ModeConfig(
    default="solid",
    options={
        "solid": _ModeOptionConfig(),
        "liquid": _ModeOptionConfig(),
        "gas": _ModeOptionConfig(input_values={"isobar": True}),
    },
)

_VLE_SWEEP_MODE_CONFIG = _ModeConfig(
    default="isotherm",
    options={
        "isotherm": _ModeOptionConfig(input_values={"isotherm": True}, required_input_keys=("temperature",)),
        "isobar": _ModeOptionConfig(input_values={"isobar": True}, required_input_keys=("pressure",)),
        "flashpoint": _ModeOptionConfig(input_values={"flashpoint": True}),
    },
)


class CRSInputBuilder:
    """Base class for property-specific CRS input builders."""

    __slots__ = (
        "property_type",
        "_job_cls",
        "_values",
        "_accepted_keys",
        "_expert_keys",
        "_expert_values",
        "_mode_config",
        "_mode",
        "_active_mode_input_keys",
        "_compounds_by_role",
        "_method",
    )

    _PROPERTY_TYPE: ClassVar[Optional[str]] = None
    _DESCRIPTION: ClassVar[str] = ""

    # Explicit builder input keys. Routes are derived from crs.json so values are
    # written to either settings.input or settings.input.property. Mode config
    # controls property-specific input keys such as isobar/isotherm/flashpoint.
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ()
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ()
    _MODE_CONFIG: ClassVar[_ModeConfig] = _ModeConfig()
    # crs.json defines primitive input types. This builder-level override only
    # constrains schema float_list keys that are semantically single values here.
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ()

    # Compound kwargs accepted by add_compound/add_solvent/add_solute methods.
    # These keys can affect CRS calculations; role config defines supported roles, count limits, and required keys.
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ()
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {}


    def __init__(
        self,
        *,
        method: CRSMethodName = "COSMO-RS",
        mode: Optional[str] = None,
        job_cls: Optional[Type[Any]] = None,
        **kwargs: Any,
    ) -> None:
        if self._PROPERTY_TYPE is None:
            raise TypeError("CRSInputBuilder is a base class; use input_builder()")

        self.property_type = self._PROPERTY_TYPE
        self._job_cls = job_cls
        self._values: Dict[str, Any] = {}
        self._accepted_keys = self._build_accepted_keys()
        self._expert_keys = self._build_expert_keys()
        self._expert_values: Dict[str, Any] = {}
        self._mode_config = self._MODE_CONFIG
        self._mode: Optional[str] = None
        self._active_mode_input_keys: Set[str] = set()
        self._compounds_by_role: Dict[str, List[Settings]] = {role: [] for role in self._compound_roles()}

        self.method = method
        self._set_mode_or_default(mode)
        self.set(**kwargs)

    @property
    def method(self) -> str:
        """COSMO-RS/SAC method name."""
        return self._method

    @method.setter
    def method(self, method: CRSMethodName) -> None:
        self._method = normalize_method(method)

    @classmethod
    def _build_accepted_keys(cls) -> Dict[str, _InputRoute]:
        data = get_crs_input_data()
        property_keys = get_block_keys(data, "property")
        accepted: Dict[str, _InputRoute] = {}

        for key in cls._EXPOSED_INPUT_KEYS:
            accepted[key] = cls._input_route_from_schema(key, data, property_keys)

        return accepted

    @classmethod
    def _build_expert_keys(cls) -> Dict[str, _InputRoute]:
        data = get_crs_input_data()
        property_keys = get_block_keys(data, "property")
        expert_keys: Dict[str, _InputRoute] = {}

        for key in _EXPERT_INPUT_KEYS:
            route = cls._input_route_from_schema(key, data, property_keys)
            if route.scope != "top_level":
                raise ValueError(f"CRS expert option {key!r} must be a top-level input key")
            expert_keys[key] = route

        return expert_keys

    @staticmethod
    def _input_route_from_schema(
        key: str,
        data: Mapping[str, Any],
        property_keys: Mapping[str, Any],
    ) -> _InputRoute:
        if key in data and data[key].get("kind") == "key":
            return _InputRoute(scope="top_level", metadata=data[key])

        if key in property_keys:
            return _InputRoute(scope="property", metadata=property_keys[key])

        raise ValueError(f"{key!r} is not a CRS top-level or PROPERTY input key")

    def __dir__(self) -> List[str]:
        entries = {
            "describe",
            "expert_options",
            "get",
            "method",
            "set",
            "set_expert_options",
            "to_job",
            "to_settings",
        }
        entries.update(self._accepted_keys)
        if self._mode_config.options:
            entries.update(("mode", "mode_options"))
            entries.difference_update(self._mode_controlled_input_keys())
        for role in self._compound_roles():
            entries.add(f"add_{role}")
            entries.add(f"add_{role}_from_adfcrs_database")
        return sorted(entries)

    def get(self, key: str, default: Any = None) -> Any:
        """Return a builder input value if set."""
        if key == "method":
            return self.method
        if key == "mode":
            return self.mode
        self._validate_input_key(key)
        return self._values.get(key, default)

    def set(self: _CRSInputBuilderT, **kwargs: Any) -> _CRSInputBuilderT:
        """Set one or more schema-backed CRS input values."""
        self._reject_mode_controlled_keys(kwargs)
        for key, value in kwargs.items():
            self._set_input_value(key, value)
        return self

    def expert_options(self) -> Tuple[str, ...]:
        """Return supported expert option names."""
        return tuple(sorted(self._expert_keys))

    def set_expert_options(self: _CRSInputBuilderT, **kwargs: Any) -> _CRSInputBuilderT:
        """Set technical top-level CRS options.

        Use expert_options() to list supported names.

        Example:
            builder.set_expert_options(pdh_correction=True)
        """
        for key, value in kwargs.items():
            self._set_expert_option(key, value)
        return self

    def describe(
        self,
        *,
        include_values: bool = False,
        include_compound_details: bool = False,
        include_expert_options: bool = False,
    ) -> Tuple[str, ...]:
        """Return concise descriptions of supported inputs and compound roles."""
        lines = [f"{self.property_type}: {self._DESCRIPTION}", f"method [top_level]: {self.method}"]

        for key, route in sorted(self._accepted_keys.items()):
            if key in self._mode_controlled_input_keys():
                continue
            lines.append(self._format_key_description(key, route.scope, route.metadata, include_values=include_values))

        if self._mode_config.options:
            lines.append(f"mode: {', '.join(self._mode_config.options)}")

        compound_keys = self._CALCULATION_COMPOUND_KEYS
        if compound_keys and include_compound_details:
            compound_key_metadata = get_block_keys(get_crs_input_data(), "compound")
            for key in compound_keys:
                lines.append(self._format_key_description(key, "compound", compound_key_metadata[key]))
        elif compound_keys:
            lines.append(f"compound_keys [compound]: {', '.join(compound_keys)}")

        if include_expert_options:
            for key, route in sorted(self._expert_keys.items()):
                lines.append(
                    self._format_key_description(
                        key,
                        "expert",
                        route.metadata,
                        include_values=include_values,
                        values=self._expert_values,
                    )
                )

        return tuple(lines)

    def to_settings(self, *, include_compounds: bool = True) -> Settings:
        """Build PLAMS Settings for this CRS input."""
        self._validate_required_inputs()
        if include_compounds:
            self._validate_compound_role_counts()
            self._validate_required_compound_keys()

        settings = Settings()
        settings.input.method = self.method
        settings.input.property._h = self.property_type

        for key, route in sorted(self._accepted_keys.items()):
            if key not in self._values:
                continue
            if route.scope == "top_level":
                settings.input[key] = self._values[key]
            elif route.scope == "property":
                settings.input.property[key] = self._values[key]

        for key in sorted(self._expert_values):
            settings.input[key] = self._expert_values[key]

        if include_compounds:
            compounds = self._compound_blocks()
            if compounds:
                settings.input.compound = compounds

        return settings

    def to_job(self, name: Optional[str] = None, **kwargs: Any) -> Any:
        """Build a CRSJob from this builder."""
        if self._job_cls is None:
            raise ValueError("to_job() requires a job class; use CRSJob.input_builder()")
        if name is not None:
            kwargs["name"] = name
        return self._job_cls(settings=self.to_settings(), **kwargs)

    def _set_input_value(self, key: str, value: Any) -> None:
        self._validate_input_key(key)
        value = self._normalize_value_by_metadata(key, value, self._accepted_keys[key].metadata)
        self._values[key] = value

    def _set_expert_option(self, key: str, value: Any) -> None:
        self._validate_expert_key(key)
        route = self._expert_keys[key]
        self._expert_values[key] = self._normalize_value_by_metadata(key, value, route.metadata)

    def _normalize_value_by_metadata(self, key: str, value: Any, metadata: Mapping[str, Any]) -> Any:
        input_type = metadata.get("type")

        if key in self._SINGLE_VALUE_INPUT_KEYS:
            if input_type == "float_list":
                self._validate_float_input(key, value)
                return value
            if input_type == "integer_list":
                self._validate_int_input(key, value)
                return value

        if input_type == "float_list":
            return self._normalize_float_list_input(key, value)
        if input_type == "integer_list":
            return self._normalize_integer_list_input(key, value)
        if input_type == "float":
            self._validate_float_input(key, value)
        elif input_type == "integer":
            self._validate_int_input(key, value)
        elif input_type == "bool":
            self._validate_bool_input(key, value)
        return value

    def _validate_float_input(self, key: str, value: Any) -> None:
        if isinstance(value, bool) or not isinstance(value, (float, int)):
            raise TypeError(f"{key} must be a single float value")

    def _normalize_float_list_input(self, key: str, value: _FloatListInput) -> Union[_FloatInput, str]:
        if isinstance(value, str):
            self._validate_list_string(key, value, float, "float")
            return value
        if isinstance(value, bool):
            raise TypeError(f"{key} must be a float value or a list of float values")
        if isinstance(value, (float, int)):
            return value
        if isinstance(value, Sequence):
            if not value:
                raise ValueError(f"{key} float list must not be empty")
            if any(isinstance(item, bool) or not isinstance(item, (float, int)) for item in value):
                raise TypeError(f"{key} must be a float value or a list of float values")
            return " ".join(str(item) for item in value)

        raise TypeError(f"{key} must be a float value or a list of float values")

    def _normalize_integer_list_input(self, key: str, value: _IntegerListInput) -> Union[int, str]:
        if isinstance(value, str):
            self._validate_list_string(key, value, int, "integer")
            return value
        if isinstance(value, bool):
            raise TypeError(f"{key} must be an integer value or a list of integer values")
        if isinstance(value, int):
            return value
        if isinstance(value, Sequence):
            if not value:
                raise ValueError(f"{key} integer list must not be empty")
            if any(isinstance(item, bool) or not isinstance(item, int) for item in value):
                raise TypeError(f"{key} must be an integer value or a list of integer values")
            return " ".join(str(item) for item in value)

        raise TypeError(f"{key} must be an integer value or a list of integer values")

    def _validate_list_string(self, key: str, value: str, converter: Callable[[str], Any], value_name: str) -> None:
        fields = value.split()
        if not fields:
            raise ValueError(f"{key} {value_name} list must not be empty")

        try:
            for field in fields:
                converter(field)
        except ValueError as exc:
            raise TypeError(f"{key} string must contain only {value_name} values") from exc

    def _validate_int_input(self, key: str, value: Any) -> None:
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{key} must be an integer value")

    def _validate_bool_input(self, key: str, value: Any) -> None:
        if not isinstance(value, bool):
            raise TypeError(f"{key} must be a boolean value")

    def _validate_input_key(self, key: str) -> None:
        if key not in self._accepted_keys:
            allowed = ", ".join(sorted(self._accepted_keys))
            raise AttributeError(f"{key!r} is not valid for {self.property_type}. Allowed input keys: {allowed}")

    def _validate_expert_key(self, key: str) -> None:
        if key not in self._expert_keys:
            allowed = ", ".join(sorted(self._expert_keys))
            raise ValueError(f"{key!r} is not a CRS expert input key. Supported expert options: {allowed}")

    def _set_mode_or_default(self, mode: Optional[str]) -> None:
        selected_mode = self._mode_config.default if mode is None else mode
        if selected_mode is not None:
            self._set_mode(selected_mode)

    def _set_mode(self, mode: str) -> None:
        if not self._mode_config.options:
            raise ValueError(f"{self.property_type} does not support mode")

        normalized_mode = mode.lower()
        option = self._mode_config.options.get(normalized_mode)
        if option is None:
            allowed = ", ".join(self._mode_config.options)
            raise ValueError(f"Unsupported mode {mode!r} for {self.property_type}. Supported modes: {allowed}")

        for key in self._active_mode_input_keys:
            self._values.pop(key, None)

        values = dict(option.input_values)
        self._values.update(values)
        self._active_mode_input_keys = set(values)
        self._mode = normalized_mode

    def _mode_controlled_input_keys(self) -> Set[str]:
        keys: Set[str] = set()
        for option in self._mode_config.options.values():
            keys.update(option.input_values)
        return keys

    def _reject_mode_controlled_keys(self, values: Mapping[str, Any]) -> None:
        provided = self._mode_controlled_input_keys().intersection(values)
        if provided:
            allowed_modes = ", ".join(self._mode_config.options)
            raise ValueError(
                f"{self.property_type} mode-controlled key(s) must be set with mode=, "
                f"not directly: {', '.join(sorted(provided))}. Supported modes: {allowed_modes}"
            )

    def _add_compound_role(
        self: _CRSInputBuilderT,
        role: str,
        compound: Union[Settings, _PathLike],
        **kwargs: Any,
    ) -> _CRSInputBuilderT:
        if role not in self._compound_roles():
            allowed = ", ".join(self._compound_roles()) or "none"
            raise ValueError(f"{self.property_type} does not support {role!r} compounds. Supported roles: {allowed}")

        self._validate_compound_role_max_count(role)
        normalized = self._normalize_compound(compound)
        self._set_compound_overrides(normalized, kwargs)
        self._compounds_by_role[role].append(normalized)
        return self

    def _add_compound_role_from_adfcrs_database(
        self: _CRSInputBuilderT,
        role: str,
        name: str,
        **kwargs: Any,
    ) -> _CRSInputBuilderT:
        if self._job_cls is None:
            raise ValueError("Database compounds require a job class; use CRSJob.input_builder()")
        return self._add_compound_role(role, self._job_cls.coskf_from_adfcrs_database(name), **kwargs)

    def _normalize_compound(self, compound: Union[Settings, _PathLike]) -> Settings:
        if isinstance(compound, Settings):
            return compound.copy()

        settings = Settings()
        settings._h = os.fspath(compound)
        return settings

    def _set_compound_overrides(self, compound: Settings, overrides: Mapping[str, Any]) -> None:
        for key, value in overrides.items():
            if value is None:
                continue
            if not self._accepts_compound_key(key):
                log(
                    f"Ignoring compound key {key!r} for {self.property_type}: it does not affect this calculation.",
                    3,
                )
                continue
            compound[key] = value

    def _accepts_compound_key(self, key: str) -> bool:
        compound_keys = get_block_keys(get_crs_input_data(), "compound")
        if key not in compound_keys:
            raise ValueError(f"{key!r} is not a COMPOUND input key in crs.json")
        return key in self._CALCULATION_COMPOUND_KEYS

    def _compound_blocks(self) -> List[Settings]:
        compounds: List[Settings] = []
        for role in self._compound_roles():
            compounds.extend(self._compounds_by_role[role])
        return compounds

    def _validate_required_inputs(self) -> None:
        required_keys = list(self._REQUIRED_INPUT_KEYS)
        if self._mode is not None:
            required_keys.extend(self._mode_config.options[self._mode].required_input_keys)

        missing = [key for key in required_keys if key not in self._values]
        if missing:
            raise ValueError(f"{self.property_type} missing required input key(s): {', '.join(missing)}")

    def _validate_compound_role_max_count(self, role: str) -> None:
        config = self._COMPOUND_ROLE_CONFIG.get(role)
        if config is not None and config.max_count is not None:
            if len(self._compounds_by_role[role]) >= config.max_count:
                raise ValueError(f"{self.property_type} supports at most {config.max_count} {role} compound(s)")

    def _validate_compound_role_counts(self) -> None:
        for role in self._compound_roles():
            config = self._COMPOUND_ROLE_CONFIG.get(role, _CompoundRoleConfig())
            count = len(self._compounds_by_role[role])
            if config.min_count is not None and count < config.min_count:
                raise ValueError(
                    f"{self.property_type} requires at least {config.min_count} {role} compound(s); got {count}"
                )
            if config.max_count is not None and count > config.max_count:
                raise ValueError(
                    f"{self.property_type} supports at most {config.max_count} {role} compound(s); got {count}"
                )

    @classmethod
    def _compound_roles(cls) -> Tuple[str, ...]:
        return tuple(cls._COMPOUND_ROLE_CONFIG)

    def _validate_required_compound_keys(self) -> None:
        for role, config in self._COMPOUND_ROLE_CONFIG.items():
            required_keys = config.required_keys
            if not required_keys:
                continue
            for index, compound in enumerate(self._compounds_by_role.get(role, ()), start=1):
                missing = [key for key in required_keys if not _has_compound_key(compound, key)]
                if missing:
                    raise ValueError(
                        f"{self.property_type} {role} #{index} missing required key(s): {', '.join(missing)}"
                    )

    def _format_key_description(
        self,
        key: str,
        scope: str,
        metadata: Mapping[str, Any],
        *,
        include_values: bool = False,
        values: Optional[Mapping[str, Any]] = None,
    ) -> str:
        input_type = metadata.get("type", "")
        unit = metadata.get("unit")
        type_unit = f"{input_type} [{unit}]" if unit else input_type
        line = f"{key} [{scope}] {type_unit}: {metadata.get('description', '')}"
        if include_values:
            value_source = self._values if values is None else values
            line += f" value: {value_source.get(key)!r}"
        return line


class _SolventRoleMixin:
    __slots__ = ()

    def add_solvent(self: _CRSInputBuilderT, compound: Union[Settings, _PathLike], **kwargs: Any) -> _CRSInputBuilderT:
        """Add a COMPOUND block as a solvent."""
        return self._add_compound_role("solvent", compound, **kwargs)

    def add_solvent_from_adfcrs_database(self: _CRSInputBuilderT, name: str, **kwargs: Any) -> _CRSInputBuilderT:
        """Add a database COMPOUND block as a solvent."""
        return self._add_compound_role_from_adfcrs_database("solvent", name, **kwargs)


class _SoluteRoleMixin:
    __slots__ = ()

    def add_solute(self: _CRSInputBuilderT, compound: Union[Settings, _PathLike], **kwargs: Any) -> _CRSInputBuilderT:
        """Add a COMPOUND block as a solute."""
        return self._add_compound_role("solute", compound, **kwargs)

    def add_solute_from_adfcrs_database(self: _CRSInputBuilderT, name: str, **kwargs: Any) -> _CRSInputBuilderT:
        """Add a database COMPOUND block as a solute."""
        return self._add_compound_role_from_adfcrs_database("solute", name, **kwargs)


class _CompoundRoleMixin:
    __slots__ = ()

    def add_compound(self: _CRSInputBuilderT, compound: Union[Settings, _PathLike], **kwargs: Any) -> _CRSInputBuilderT:
        """Add a generic COMPOUND block."""
        return self._add_compound_role("compound", compound, **kwargs)

    def add_compound_from_adfcrs_database(self: _CRSInputBuilderT, name: str, **kwargs: Any) -> _CRSInputBuilderT:
        """Add a generic database COMPOUND block."""
        return self._add_compound_role_from_adfcrs_database("compound", name, **kwargs)


class _TemperatureMixin:
    __slots__ = ()

    @property
    def temperature(self) -> Optional[_FloatInput]:
        """Temperature in K."""
        return self.get("temperature")

    @temperature.setter
    def temperature(self, value: _FloatInput) -> None:
        self._set_input_value("temperature", value)


class _TemperatureListMixin:
    __slots__ = ()

    @property
    def temperature(self) -> Optional[_FloatListInput]:
        """Temperature in K, as one value or [low, high, nsteps]."""
        return self.get("temperature")

    @temperature.setter
    def temperature(self, value: _FloatListInput) -> None:
        self._set_input_value("temperature", value)


class _PressureMixin:
    __slots__ = ()

    @property
    def pressure(self) -> Optional[_FloatInput]:
        """Pressure in bar."""
        return self.get("pressure")

    @pressure.setter
    def pressure(self, value: _FloatInput) -> None:
        self._set_input_value("pressure", value)


class _PressureListMixin:
    __slots__ = ()

    @property
    def pressure(self) -> Optional[_FloatListInput]:
        """Pressure in bar, as one value or [low, high, nsteps]."""
        return self.get("pressure")

    @pressure.setter
    def pressure(self, value: _FloatListInput) -> None:
        self._set_input_value("pressure", value)


class _MassFractionMixin:
    __slots__ = ()

    @property
    def massfraction(self) -> Optional[bool]:
        """Interpret compound fractions as mass fractions instead of mole fractions."""
        return self.get("massfraction")

    @massfraction.setter
    def massfraction(self, value: bool) -> None:
        self._set_input_value("massfraction", value)


class _DensitySolventMixin:
    __slots__ = ()

    @property
    def densitysolvent(self) -> Optional[_FloatInput]:
        """Solvent density in kg/L."""
        return self.get("densitysolvent")

    @densitysolvent.setter
    def densitysolvent(self, value: _FloatInput) -> None:
        self._set_input_value("densitysolvent", value)


class _SigmaMomentMixin:
    __slots__ = ()

    @property
    def sigmamomentpower(self) -> Optional[_IntegerListInput]:
        """Powers used for sigma moment output."""
        return self.get("sigmamomentpower")

    @sigmamomentpower.setter
    def sigmamomentpower(self, value: _IntegerListInput) -> None:
        self._set_input_value("sigmamomentpower", value)

    @property
    def sigmamomenthblevel(self) -> Optional[_IntegerListInput]:
        """Hydrogen-bond sigma moment levels; level 1 uses cutoff, higher levels use cutoffbase + level * cutoffstep."""
        return self.get("sigmamomenthblevel")

    @sigmamomenthblevel.setter
    def sigmamomenthblevel(self, value: _IntegerListInput) -> None:
        self._set_input_value("sigmamomenthblevel", value)

    @property
    def sigmamomenthbcutoff(self) -> Optional[_FloatInput]:
        """Cutoff for level 1 hydrogen-bond sigma moment in e/Angstrom^2."""
        return self.get("sigmamomenthbcutoff")

    @sigmamomenthbcutoff.setter
    def sigmamomenthbcutoff(self, value: _FloatInput) -> None:
        self._set_input_value("sigmamomenthbcutoff", value)

    @property
    def sigmamomenthbcutoffbase(self) -> Optional[_FloatInput]:
        """Cutoff base for higher-level hydrogen-bond sigma moment in e/Angstrom^2."""
        return self.get("sigmamomenthbcutoffbase")

    @sigmamomenthbcutoffbase.setter
    def sigmamomenthbcutoffbase(self, value: _FloatInput) -> None:
        self._set_input_value("sigmamomenthbcutoffbase", value)

    @property
    def sigmamomenthbcutoffstep(self) -> Optional[_FloatInput]:
        """Cutoff step for higher-level hydrogen-bond sigma moment in e/Angstrom^2."""
        return self.get("sigmamomenthbcutoffstep")

    @sigmamomenthbcutoffstep.setter
    def sigmamomenthbcutoffstep(self, value: _FloatInput) -> None:
        self._set_input_value("sigmamomenthbcutoffstep", value)


class _SigmaPotentialMixin:
    __slots__ = ()

    @property
    def estpotential(self) -> Optional[bool]:
        """Estimate sigma potential where the sigma profile has no area."""
        return self.get("estpotential")

    @estpotential.setter
    def estpotential(self, value: bool) -> None:
        self._set_input_value("estpotential", value)


class _SigmaMixin:
    __slots__ = ()

    @property
    def nprofile(self) -> Optional[int]:
        """Number of sigma profile or sigma potential data points."""
        return self.get("nprofile")

    @nprofile.setter
    def nprofile(self, value: int) -> None:
        self._set_input_value("nprofile", value)

    @property
    def sigmamax(self) -> Optional[_FloatInput]:
        """Maximum sigma value for sigma profile or sigma potential output."""
        return self.get("sigmamax")

    @sigmamax.setter
    def sigmamax(self, value: _FloatInput) -> None:
        self._set_input_value("sigmamax", value)


class _SolubilityModeMixin:
    __slots__ = ()

    @property
    def mode(self) -> Optional[_SolubilityMode]:
        """Solubility mode: solid, liquid, or gas."""
        return self._mode

    @mode.setter
    def mode(self, mode: _SolubilityMode) -> None:
        self._set_mode(mode)

    @property
    def mode_options(self) -> Tuple[str, ...]:
        """Supported solubility mode names."""
        return tuple(self._mode_config.options)


class _VLESweepMixin:
    """Inputs for VLE-style sweep properties."""
    __slots__ = ()

    @property
    def mode(self) -> Optional[_VLESweepMode]:
        """VLE sweep mode: isotherm, isobar, or flashpoint."""
        return self._mode

    @mode.setter
    def mode(self, mode: _VLESweepMode) -> None:
        self._set_mode(mode)

    @property
    def mode_options(self) -> Tuple[str, ...]:
        """Supported VLE sweep mode names."""
        return tuple(self._mode_config.options)

    @property
    def nfrac(self) -> Optional[int]:
        """Mixture-fraction resolution: binary n+5, ternary (n+1)*(n+2)/2, composition line n+1."""
        return self.get("nfrac")

    @nfrac.setter
    def nfrac(self, value: int) -> None:
        self._set_input_value("nfrac", value)


class ACTIVITYCOEFInputBuilder(
    _TemperatureMixin,
    _MassFractionMixin,
    _DensitySolventMixin,
    _SolventRoleMixin,
    _SoluteRoleMixin,
    CRSInputBuilder,
):
    """Builder for ACTIVITYCOEF CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "ACTIVITYCOEF"
    _DESCRIPTION: ClassVar[str] = "Activity coefficients in a solvent mixture."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "massfraction", "densitysolvent")
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1", "density") + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "solvent": _CompoundRoleConfig(required_keys=("frac1",)),
        "solute": _CompoundRoleConfig(),
    }


class LOGPInputBuilder(
    _TemperatureMixin,
    _MassFractionMixin,
    _SolventRoleMixin,
    _SoluteRoleMixin,
    CRSInputBuilder,
):
    """Builder for LOGP CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "LOGP"
    _DESCRIPTION: ClassVar[str] = "Partition coefficients between two immiscible solvent phases."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "massfraction", "volumequotient")
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1", "frac2", "density")
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "solvent": _CompoundRoleConfig(min_count=2, required_keys=("frac1", "frac2")),
        "solute": _CompoundRoleConfig(min_count=1),
    }

    @property
    def volumequotient(self) -> Optional[_FloatInput]:
        """Molar-volume ratio of solvent 1 to solvent 2."""
        return self.get("volumequotient")

    @volumequotient.setter
    def volumequotient(self, value: _FloatInput) -> None:
        self._set_input_value("volumequotient", value)


class SOLUBILITYInputBuilder(
    _TemperatureListMixin,
    _PressureMixin,
    _MassFractionMixin,
    _DensitySolventMixin,
    _SolubilityModeMixin,
    _SolventRoleMixin,
    _SoluteRoleMixin,
    CRSInputBuilder,
):
    """Builder for SOLUBILITY CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "SOLUBILITY"
    _DESCRIPTION: ClassVar[str] = "Solubility of solutes in a solvent mixture or under gas-pressure conditions."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = (
        "temperature",
        "pressure",
        "massfraction",
        "densitysolvent",
        "isobar",
    )
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("pressure",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1", "density") + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "solvent": _CompoundRoleConfig(required_keys=("frac1",)),
        "solute": _CompoundRoleConfig(),
    }
    _MODE_CONFIG: ClassVar[_ModeConfig] = _SOLUBILITY_MODE_CONFIG

    def _validate_required_compound_keys(self) -> None:
        super()._validate_required_compound_keys()
        if self.mode != "solid":
            return

        for index, compound in enumerate(self._compounds_by_role.get("solute", ()), start=1):
            if not (_has_compound_key(compound, "meltingpoint") and _has_compound_key(compound, "hfusion")):
                raise ValueError(
                    f"{self.property_type} solute #{index} requires meltingpoint and hfusion for mode='solid'"
                )


class PURESOLUBILITYInputBuilder(
    _TemperatureListMixin,
    _PressureMixin,
    _SolubilityModeMixin,
    _SolventRoleMixin,
    _SoluteRoleMixin,
    CRSInputBuilder,
):
    """Builder for PURESOLUBILITY CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "PURESOLUBILITY"
    _DESCRIPTION: ClassVar[str] = "Solubility of a solute in pure solvents over a temperature range."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "pressure", "isobar")
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("pressure",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1", "density") + _FUSION_KEYS + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "solvent": _CompoundRoleConfig(min_count=1),
        "solute": _CompoundRoleConfig(min_count=1, max_count=1),
    }
    _MODE_CONFIG: ClassVar[_ModeConfig] = _SOLUBILITY_MODE_CONFIG

    def _validate_required_compound_keys(self) -> None:
        super()._validate_required_compound_keys()
        if self.mode != "solid":
            return

        for index, compound in enumerate(self._compounds_by_role.get("solute", ()), start=1):
            if not (_has_compound_key(compound, "meltingpoint") and _has_compound_key(compound, "hfusion")):
                raise ValueError(
                    f"{self.property_type} solute #{index} requires meltingpoint and hfusion for mode='solid'"
                )


class VAPORPRESSUREInputBuilder(
    _TemperatureMixin,
    _MassFractionMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for VAPORPRESSURE CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "VAPORPRESSURE"
    _DESCRIPTION: ClassVar[str] = "Vapor pressure of a mixture at fixed temperature."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "massfraction")
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",) + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(required_keys=("frac1",)),
    }


class PUREVAPORPRESSUREInputBuilder(
    _TemperatureListMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for PUREVAPORPRESSURE CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "PUREVAPORPRESSURE"
    _DESCRIPTION: ClassVar[str] = "Pure-compound vapor pressure over a temperature range."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(),
    }


class BOILINGPOINTInputBuilder(
    _PressureListMixin,
    _MassFractionMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for BOILINGPOINT CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "BOILINGPOINT"
    _DESCRIPTION: ClassVar[str] = "Boiling temperature of a mixture for a pressure range."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("pressure", "massfraction")
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("pressure",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",) + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(),
    }


class PUREBOILINGPOINTInputBuilder(
    _PressureListMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for PUREBOILINGPOINT CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "PUREBOILINGPOINT"
    _DESCRIPTION: ClassVar[str] = "Pure-compound boiling point over a pressure range."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("pressure",)
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("pressure",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(),
    }


class FLASHPOINTInputBuilder(
    _MassFractionMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for FLASHPOINT CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "FLASHPOINT"
    _DESCRIPTION: ClassVar[str] = "Flash point of a mixture using user-supplied pure-compound flash points."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("massfraction",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1", "flashpoint") + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(required_keys=("frac1",)),
    }


class BINMIXCOEFInputBuilder(
    _TemperatureMixin,
    _PressureMixin,
    _MassFractionMixin,
    _VLESweepMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for BINMIXCOEF CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "BINMIXCOEF"
    _DESCRIPTION: ClassVar[str] = "Binary-mixture coefficients over a composition range."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = (
        "temperature",
        "pressure",
        "massfraction",
        "isotherm",
        "isobar",
        "flashpoint",
        "nfrac",
    )
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",) + _VAPOR_PRESSURE_KEYS + ("flashpoint",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "pressure")
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(min_count=2, max_count=2),
    }
    _MODE_CONFIG: ClassVar[_ModeConfig] = _VLE_SWEEP_MODE_CONFIG


class TERNARYMIXInputBuilder(
    _TemperatureMixin,
    _PressureMixin,
    _MassFractionMixin,
    _VLESweepMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for TERNARYMIX CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "TERNARYMIX"
    _DESCRIPTION: ClassVar[str] = "Ternary mixture property sweep over composition space."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = (
        "temperature",
        "pressure",
        "massfraction",
        *_VLE_SWEEP_PROPERTY_KEYS,
    )
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "pressure")
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",) + _VAPOR_PRESSURE_KEYS + ("flashpoint",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(min_count=3, max_count=3),
    }
    _MODE_CONFIG: ClassVar[_ModeConfig] = _VLE_SWEEP_MODE_CONFIG


class COMPOSITIONLINEInputBuilder(
    _TemperatureMixin,
    _PressureMixin,
    _MassFractionMixin,
    _VLESweepMixin,
    _SolventRoleMixin,
    CRSInputBuilder,
):
    """Builder for COMPOSITIONLINE CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "COMPOSITIONLINE"
    _DESCRIPTION: ClassVar[str] = "Composition-line calculation between two endpoint phase compositions."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = (
        "temperature",
        "pressure",
        "massfraction",
        *_VLE_SWEEP_PROPERTY_KEYS,
    )
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "pressure")
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = (
        "frac1",
        "frac2",
    ) + _VAPOR_PRESSURE_KEYS + ("flashpoint",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "solvent": _CompoundRoleConfig(),
    }
    _MODE_CONFIG: ClassVar[_ModeConfig] = _VLE_SWEEP_MODE_CONFIG


class LLEInputBuilder(
    _TemperatureMixin,
    _MassFractionMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for LLE CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "LLE"
    _DESCRIPTION: ClassVar[str] = "Liquid-liquid equilibrium for a ternary mixture."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "massfraction")
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(required_keys=("frac1",)),
    }


class STABILITYInputBuilder(
    _TemperatureMixin,
    _MassFractionMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for STABILITY CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "STABILITY"
    _DESCRIPTION: ClassVar[str] = "Michelsen tangent-plane-distance stability test for a feed composition."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "massfraction")
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(required_keys=("frac1",)),
    }


class SIGMAPROFILEInputBuilder(
    _MassFractionMixin,
    _SigmaMixin,
    _SigmaMomentMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for SIGMAPROFILE CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "SIGMAPROFILE"
    _DESCRIPTION: ClassVar[str] = "Sigma profile for a solvent mixture."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("massfraction", *_SIGMA_KEYS, *_SIGMA_MOMENT_KEYS)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(required_keys=("frac1",)),
    }


class PURESIGMAPROFILEInputBuilder(
    _SigmaMixin,
    _SigmaMomentMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for PURESIGMAPROFILE CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "PURESIGMAPROFILE"
    _DESCRIPTION: ClassVar[str] = "Sigma profile for pure compounds."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = (*_SIGMA_KEYS, *_SIGMA_MOMENT_KEYS)
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(),
    }


class SIGMAPOTENTIALInputBuilder(
    _TemperatureMixin,
    _MassFractionMixin,
    _SigmaMixin,
    _SigmaPotentialMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for SIGMAPOTENTIAL CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "SIGMAPOTENTIAL"
    _DESCRIPTION: ClassVar[str] = "Sigma potential for a solvent mixture."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "massfraction", *_SIGMA_KEYS, "estpotential")
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(required_keys=("frac1",)),
    }


class PURESIGMAPOTENTIALInputBuilder(
    _TemperatureMixin,
    _SigmaMixin,
    _SigmaPotentialMixin,
    _CompoundRoleMixin,
    CRSInputBuilder,
):
    """Builder for PURESIGMAPOTENTIAL CRS input settings."""

    __slots__ = ()

    _PROPERTY_TYPE: ClassVar[str] = "PURESIGMAPOTENTIAL"
    _DESCRIPTION: ClassVar[str] = "Sigma potential for pure compounds."
    _EXPOSED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", *_SIGMA_KEYS, "estpotential")
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _SINGLE_VALUE_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(),
    }


_CRS_INPUT_BUILDER_CLASSES: Dict[str, Type[CRSInputBuilder]] = {
    "ACTIVITYCOEF": ACTIVITYCOEFInputBuilder,
    "LOGP": LOGPInputBuilder,
    "SOLUBILITY": SOLUBILITYInputBuilder,
    "PURESOLUBILITY": PURESOLUBILITYInputBuilder,
    "VAPORPRESSURE": VAPORPRESSUREInputBuilder,
    "PUREVAPORPRESSURE": PUREVAPORPRESSUREInputBuilder,
    "BOILINGPOINT": BOILINGPOINTInputBuilder,
    "PUREBOILINGPOINT": PUREBOILINGPOINTInputBuilder,
    "FLASHPOINT": FLASHPOINTInputBuilder,
    "BINMIXCOEF": BINMIXCOEFInputBuilder,
    "TERNARYMIX": TERNARYMIXInputBuilder,
    "COMPOSITIONLINE": COMPOSITIONLINEInputBuilder,
    "LLE": LLEInputBuilder,
    "STABILITY": STABILITYInputBuilder,
    "SIGMAPROFILE": SIGMAPROFILEInputBuilder,
    "PURESIGMAPROFILE": PURESIGMAPROFILEInputBuilder,
    "SIGMAPOTENTIAL": SIGMAPOTENTIALInputBuilder,
    "PURESIGMAPOTENTIAL": PURESIGMAPOTENTIALInputBuilder,
}


def input_builder(
    property_type: str,
    *,
    method: CRSMethodName = "COSMO-RS",
    mode: Optional[str] = None,
    job_cls: Optional[Type[Any]] = None,
    **kwargs: Any,
) -> CRSInputBuilder:
    """Return a property-specific CRS input builder."""
    normalized_property_type = normalize_property_type(property_type)
    builder_class = _CRS_INPUT_BUILDER_CLASSES[normalized_property_type]
    return builder_class(method=method, mode=mode, job_cls=job_cls, **kwargs)


def methods() -> Tuple[str, ...]:
    """Return supported CRS method names."""
    return CRS_METHODS


def normalize_method(method: CRSMethodName) -> str:
    """Normalize and validate a CRS method name."""
    normalized = _METHOD_ALIASES.get(method.upper(), method.upper())
    if normalized not in CRS_METHODS:
        allowed = ", ".join(CRS_METHODS)
        raise ValueError(f"Unsupported CRS method {method!r}. Supported methods: {allowed}")
    return normalized


def normalize_property_type(property_type: str) -> str:
    """Normalize and validate a supported builder property type."""
    normalized = property_type.upper()
    if normalized not in _CRS_INPUT_BUILDER_CLASSES:
        allowed = ", ".join(sorted(_CRS_INPUT_BUILDER_CLASSES))
        raise ValueError(f"Unsupported CRS input builder property {property_type!r}. Supported properties: {allowed}")
    return normalized


def _has_compound_key(compound: Settings, key: str) -> bool:
    return key in compound and compound[key] is not None

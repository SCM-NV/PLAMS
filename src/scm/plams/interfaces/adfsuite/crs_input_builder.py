from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, ClassVar, Dict, List, Literal, Mapping, Optional, Sequence, Set, Tuple, Type, TypeVar, Union

from scm.plams.core.functions import log
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.crs_definitions import CRS_METHODS, get_block_keys, get_crs_input_data


__all__ = [
    "ACTIVITYCOEFInputBuilder",
    "BINMIXCOEFInputBuilder",
    "CRSInputBuilder",
    "SOLUBILITYInputBuilder",
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

_CRSInputBuilderT = TypeVar("_CRSInputBuilderT", bound="CRSInputBuilder")

_METHOD_ALIASES = {
    "COSMORS": "COSMO-RS",
    "COSMOSAC": "COSMOSAC2013",
}

_VAPOR_PRESSURE_KEYS = ("pvap", "tvap", "vp_equation", "vp_params")
_FUSION_KEYS = ("meltingpoint", "hfusion", "cpfusion")

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

_BINMIXCOEF_MODE_CONFIG = _ModeConfig(
    default="isotherm",
    options={
        "isotherm": _ModeOptionConfig(input_values={"isotherm": True}),
        "isobar": _ModeOptionConfig(input_values={"isobar": True}),
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
    # Runtime value normalization for explicit builder inputs. These groups keep
    # constructor kwargs, set(**kwargs), and property setters consistent.
    _FLOAT_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ()
    _FLOAT_LIST_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ()
    _INT_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ()
    _BOOL_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ()

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
            "get",
            "method",
            "set",
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

    def describe(self, *, include_values: bool = False, include_compound_details: bool = False) -> Tuple[str, ...]:
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
        value = self._normalize_input_value(key, value)
        self._values[key] = value

    def _normalize_input_value(self, key: str, value: Any) -> Any:
        if key in self._FLOAT_LIST_INPUT_KEYS:
            return self._normalize_float_range_input(key, value)
        if key in self._FLOAT_INPUT_KEYS:
            self._validate_float_input(key, value)
        elif key in self._INT_INPUT_KEYS:
            self._validate_int_input(key, value)
        elif key in self._BOOL_INPUT_KEYS:
            self._validate_bool_input(key, value)
        return value

    def _validate_float_input(self, key: str, value: Any) -> None:
        if isinstance(value, bool) or not isinstance(value, (float, int)):
            raise TypeError(f"{key} must be a single float value")

    def _normalize_float_range_input(self, key: str, value: _FloatListInput) -> Union[_FloatInput, str]:
        if isinstance(value, str):
            self._validate_float_range_string(key, value)
            return value
        if isinstance(value, bool):
            raise TypeError(f"{key} must be a float value or a range [low, high, nsteps]")
        if isinstance(value, (float, int)):
            return value
        if isinstance(value, Sequence):
            if len(value) != 3:
                raise ValueError(f"{key} range must contain [low, high, nsteps]")
            low, high, nsteps = value
            if (
                isinstance(low, bool)
                or isinstance(high, bool)
                or isinstance(nsteps, bool)
                or not isinstance(low, (float, int))
                or not isinstance(high, (float, int))
                or not isinstance(nsteps, int)
            ):
                raise TypeError(f"{key} range must be [float, float, int]")
            return f"{low} {high} {nsteps}"

        raise TypeError(f"{key} must be a float value or a range [low, high, nsteps]")

    def _validate_float_range_string(self, key: str, value: str) -> None:
        fields = value.split()
        if len(fields) not in (1, 3):
            raise ValueError(f"{key} must be a single value or range: low high nsteps")

        try:
            float(fields[0])
            if len(fields) == 3:
                float(fields[1])
                int(fields[2])
        except ValueError as exc:
            raise TypeError(f"{key} range string must be: float float int") from exc

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
        missing = [key for key in self._REQUIRED_INPUT_KEYS if key not in self._values]
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
    ) -> str:
        input_type = metadata.get("type", "")
        unit = metadata.get("unit")
        type_unit = f"{input_type} [{unit}]" if unit else input_type
        line = f"{key} [{scope}] {type_unit}: {metadata.get('description', '')}"
        if include_values:
            line += f" value: {self._values.get(key)!r}"
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
        """Temperature in K as a single value."""
        return self.get("temperature")

    @temperature.setter
    def temperature(self, value: _FloatInput) -> None:
        self._set_input_value("temperature", value)


class _TemperatureListMixin:
    __slots__ = ()

    @property
    def temperature(self) -> Optional[_FloatListInput]:
        """Temperature in K as a single value or range [temperature, temperature_high, nsteps]."""
        return self.get("temperature")

    @temperature.setter
    def temperature(self, value: _FloatListInput) -> None:
        self._set_input_value("temperature", value)


class _PressureMixin:
    __slots__ = ()

    @property
    def pressure(self) -> Optional[_FloatInput]:
        """Pressure in bar as a single value."""
        return self.get("pressure")

    @pressure.setter
    def pressure(self, value: _FloatInput) -> None:
        self._set_input_value("pressure", value)


class _PressureListMixin:
    __slots__ = ()

    @property
    def pressure(self) -> Optional[_FloatListInput]:
        """Pressure in bar as a single value or range [pressure, pressure_high, nsteps]."""
        return self.get("pressure")

    @pressure.setter
    def pressure(self, value: _FloatListInput) -> None:
        self._set_input_value("pressure", value)


class _MassFractionMixin:
    __slots__ = ()

    @property
    def massfraction(self) -> Optional[bool]:
        """Use mass fractions instead of molar fractions."""
        return self.get("massfraction")

    @massfraction.setter
    def massfraction(self, value: bool) -> None:
        self._set_input_value("massfraction", value)


class _DensitySolventMixin:
    __slots__ = ()

    @property
    def densitysolvent(self) -> Optional[_FloatInput]:
        """Density of the solvent."""
        return self.get("densitysolvent")

    @densitysolvent.setter
    def densitysolvent(self, value: _FloatInput) -> None:
        self._set_input_value("densitysolvent", value)


class _NFracMixin:
    __slots__ = ()

    @property
    def nfrac(self) -> Optional[int]:
        """Number of different mixtures for mixture sweeps."""
        return self.get("nfrac")

    @nfrac.setter
    def nfrac(self, value: int) -> None:
        self._set_input_value("nfrac", value)


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
        return self._mode_config.options


class _VLESweepModeMixin:
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
        return self._mode_config.options


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
    _FLOAT_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "densitysolvent")
    _BOOL_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("massfraction",)
    _CALCULATION_COMPOUND_KEYS: ClassVar[Tuple[str, ...]] = ("frac1", "density") + _VAPOR_PRESSURE_KEYS
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "solvent": _CompoundRoleConfig(required_keys=("frac1",)),
        "solute": _CompoundRoleConfig(),
    }


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
    _FLOAT_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("pressure", "densitysolvent")
    _FLOAT_LIST_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _BOOL_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("massfraction",)
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


class BINMIXCOEFInputBuilder(
    _TemperatureMixin,
    _PressureMixin,
    _MassFractionMixin,
    _NFracMixin,
    _VLESweepModeMixin,
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
    _REQUIRED_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature",)
    _FLOAT_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("temperature", "pressure")
    _INT_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("nfrac",)
    _BOOL_INPUT_KEYS: ClassVar[Tuple[str, ...]] = ("massfraction",)
    _COMPOUND_ROLE_CONFIG: ClassVar[Mapping[str, _CompoundRoleConfig]] = {
        "compound": _CompoundRoleConfig(min_count=2, max_count=2),
    }
    _MODE_CONFIG: ClassVar[_ModeConfig] = _BINMIXCOEF_MODE_CONFIG


_CRS_INPUT_BUILDER_CLASSES: Dict[str, Type[CRSInputBuilder]] = {
    "ACTIVITYCOEF": ACTIVITYCOEFInputBuilder,
    "SOLUBILITY": SOLUBILITYInputBuilder,
    "BINMIXCOEF": BINMIXCOEFInputBuilder,
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

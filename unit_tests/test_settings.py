import pytest
from unittest.mock import patch

from scm.plams.core.settings import Settings, ConfigSettings
from scm.plams.unit_tests.test_helpers import assert_config_as_expected


class TestSettings:
    """
    Test suite for the Settings class
    """

    @pytest.fixture
    def flat_settings(self):
        """
        Simple flat settings
        """
        return self.get_flat_settings()

    def get_flat_settings(self):
        return Settings({"a": "ant", "B": "bEAR", "c ": " c at", "_d_": "dingo", "__d__": "dog", "3": "emu"})

    @pytest.fixture
    def nested_settings(self):
        """
        Nested settings with various primitive types for keys
        """
        return self.get_nested_settings()

    def get_nested_settings(self):
        return Settings(
            {
                "elements": {
                    "H": {
                        "name": "Hydrogen",
                        "num": 1,
                        "mass": 1.008,
                        "metal": False,
                        "common_isotopes": [
                            {"name": "H1", "mass": 1, "abundance": 99.99},
                            {"name": "H2", "mass": 2, "abundance": 0.01},
                        ],
                        "properties": ["gas", "dimer"],
                    },
                    "O": {
                        "name": "Oxygen",
                        "num": 8,
                        "mass": 15.999,
                        "metal": False,
                        "common_isotopes": [
                            {"name": "O16", "mass": 16, "abundance": 99.8},
                            {"name": "O18", "mass": 18, "abundance": 0.02},
                        ],
                    },
                    "Fe": {"name": "Iron", "num": 26, "mass": 55.845, "metal": True, "properties": ["common"]},
                },
                False: {"s": "string", True: "bool", 42: "int", 42.99: "float"},
                1: {1: "one", 2.0: "two", 3.0001: "three"},
            }
        )

    @pytest.fixture
    def extra_nested_settings(self):
        """
        Extra nested settings which can be patched into nested settings
        """
        return self.get_extra_nested_settings()

    def get_extra_nested_settings(self):
        return Settings(
            {
                "elements": {
                    "Fe": {
                        "gas_at_rt": False,
                        "mass": 55.846,
                        "common_isotopes": [
                            {"name": "Fe54", "mass": 54, "abundance": 5.85},
                            {"name": "Fe56", "mass": 56, "abundance": 91.8},
                            {"name": "Fe57", "mass": 57, "abundance": 2.12},
                            {"name": "Fe58", "mass": 58, "abundance": 0.28},
                        ],
                        "properties": ["forms oxides"],
                    },
                    "He": {"name": "Helium", "num": 2, "mass": 4.003, "metal": False, "gas_at_rt": True},
                },
            }
        )

    @pytest.fixture
    def nested_mutable_settings(self):
        """
        Nested settings with mutable keys
        """
        return self.get_nested_mutable_settings()

    def get_nested_mutable_settings(self):
        return Settings(
            {
                ("immutable", "key", 1): ("immutable", "value", 1),
                ("immutable", "key", 2): ["mutable", "value", 2],
                ("immutable", "key", 3): {
                    ("immutable", "key", 4): ["mutable", "value", 4],
                    ("immutable", "key", 5): (("immutable", "value", 5), ["mutable", "value", 5]),
                    ("immutable", "key", 6): [("immutable", "value", 6), ["mutable", "value", 6]],
                },
            }
        )

    def test_settings_access_as_expected(self, flat_settings, nested_settings, nested_mutable_settings):
        # Happy, has key (case-insensitive for strings), returns value
        assert flat_settings["a"] == "ant"
        assert flat_settings["B"] == "bEAR"
        assert flat_settings["c "] == " c at"
        assert flat_settings["_d_"] == "dingo"
        assert flat_settings["__d__"] == "dog"
        assert flat_settings["3"] == "emu"

        assert nested_settings["elements"]["H"]["name"] == "Hydrogen"
        assert nested_settings["eLEments"]["h"]["Num"] == 1
        assert nested_settings["ELEMENTS"]["H"]["mass"] == 1.008
        assert not nested_settings["elements"]["h"]["metal"]
        assert nested_settings["elements"]["H"]["common_isotopes"][0]["mass"] == 1
        assert nested_settings["elements"]["H"]["properties"][1] == "dimer"
        assert nested_settings[False]["S"] == "string"
        assert nested_settings[False][True] == "bool"
        assert nested_settings[False][42] == "int"
        assert nested_settings[False][42.99] == "float"
        assert nested_settings[1][1] == "one"
        assert nested_settings[1][2] == "two"
        assert nested_settings[1][2.00] == "two"
        assert nested_settings[1][3.0001] == "three"

        assert nested_mutable_settings[("immutable", "key", 1)] == ("immutable", "value", 1)
        assert nested_mutable_settings[("immutable", "key", 3)][("immutable", "key", 6)] == [
            ("immutable", "value", 6),
            ["mutable", "value", 6],
        ]

        # Unhappy, does not have key, returns empty settings object
        assert flat_settings["_"] == Settings()
        assert nested_settings["elements"]["_"]["__"] == Settings()
        assert nested_mutable_settings[("IMMUTABLE", "KEY", 1)] == Settings()

    def test_settings_dot_access_as_expected(self, flat_settings, nested_settings):
        # Happy, keys are strings without whitespace and double underscores, returns value
        assert flat_settings.a == "ant"
        assert flat_settings.B == "bEAR"
        assert flat_settings._d_ == "dingo"
        assert nested_settings.elements.H.name == "Hydrogen"
        assert nested_settings.elements.H.properties == ["gas", "dimer"]
        assert nested_settings.elements.H.common_isotopes[0].mass == 1
        assert nested_settings.ELEMENTS.Fe.MaSS == 55.845
        assert nested_settings[False].s == "string"

        # Unhappy, other keys, returns empty settings object
        assert flat_settings.c == Settings()

        # Unhappy, key does not exist, returns empty settings object
        assert flat_settings.foo == Settings()
        assert flat_settings.foo.bar == Settings()
        assert nested_settings.elements.foo.bar == Settings()

        # Unhappy, private key, gives error
        with pytest.raises(AttributeError):
            assert flat_settings.__d__ == Settings()

    def test_settings_find_case_returns_case_insensitive_key_or_original(self, flat_settings, nested_settings):
        # Happy, keys are returned with a case-insensitive match
        assert flat_settings.find_case("a") == "a"
        assert flat_settings.find_case("b") == "B"
        assert flat_settings.find_case("_D_") == "_d_"

        # Unhappy, keys are not strings/present in settings so are just parroted back
        assert flat_settings.find_case(" C At") == " C At"
        assert flat_settings.find_case("_") == "_"
        assert not nested_settings.find_case(False)
        assert nested_settings.find_case(True)
        assert nested_settings.find_case(1) == 1
        assert nested_settings.find_case(3) == 3

    def test_settings_copy_deep_copies_nested_settings_but_shallow_copies_contents(
        self, nested_settings, nested_mutable_settings
    ):
        # Copy settings
        nested_settings_copy = nested_settings.copy()
        nested_mutable_settings_copy = nested_mutable_settings.copy()

        # Settings objects are distinct objects but elements are not
        assert nested_settings.elements is not nested_settings_copy.elements
        assert nested_mutable_settings[("immutable", "key", 2)] is nested_mutable_settings_copy[("immutable", "key", 2)]
        assert (
            nested_mutable_settings[("immutable", "key", 3)]
            is not nested_mutable_settings_copy[("immutable", "key", 3)]
        )
        assert (
            nested_mutable_settings[("immutable", "key", 3)][("immutable", "key", 6)]
            is nested_mutable_settings_copy[("immutable", "key", 3)][("immutable", "key", 6)]
        )

        # Call out consequence of the above
        nested_mutable_settings_copy[("immutable", "key", 2)] = "new value"
        nested_mutable_settings_copy[("immutable", "key", 3)][("immutable", "key", 4)][1] = "new value"

        assert nested_mutable_settings[("immutable", "key", 2)] == ["mutable", "value", 2]
        assert nested_mutable_settings[("immutable", "key", 3)][("immutable", "key", 4)] == ["mutable", "new value", 4]

    def test_settings_soft_update_does_not_overwrite_existing(self, nested_settings, extra_nested_settings):
        # Soft update new settings
        nested_settings.soft_update(extra_nested_settings)

        # Existing keys should have their values unchanged, new keys should be set
        assert nested_settings.elements.Fe.mass == 55.845
        assert not nested_settings.elements.Fe.gas_at_rt
        assert nested_settings.elements.Fe.properties == ["common"]
        assert nested_settings.elements.Fe.common_isotopes[0].mass == 54
        assert nested_settings.elements.He.name == "Helium"
        assert nested_settings.elements.He.mass == 4.003
        assert nested_settings.elements.He.gas_at_rt

    def test_settings_update_does_overwrite_existing(self, nested_settings, extra_nested_settings):
        # Update new settings
        nested_settings.update(extra_nested_settings)

        # Existing keys should have their values changed, new keys should be set
        assert nested_settings.elements.Fe.mass == 55.846
        assert not nested_settings.elements.Fe.gas_at_rt
        assert nested_settings.elements.Fe.properties == ["forms oxides"]
        assert nested_settings.elements.Fe.common_isotopes[0].mass == 54
        assert nested_settings.elements.He.name == "Helium"
        assert nested_settings.elements.He.mass == 4.003
        assert nested_settings.elements.He.gas_at_rt

    def test_settings_merge_copies_and_soft_updates(self, nested_settings, extra_nested_settings):
        # Merge is a simple combination of copy and soft update
        merged = nested_settings.merge(extra_nested_settings)

        # Existing keys should have their values unchanged, new keys should be set
        assert merged.elements.Fe.mass == 55.845
        assert not merged.elements.Fe.gas_at_rt
        assert merged.elements.Fe.properties == ["common"]
        assert merged.elements.Fe.common_isotopes[0].mass == 54
        assert merged.elements.He.name == "Helium"
        assert merged.elements.He.mass == 4.003
        assert merged.elements.He.gas_at_rt

        # Settings objects should be deep copied
        assert merged.elements is not nested_settings.elements

    def test_settings_remove_deletes_matching_keys(self, nested_settings):
        nested_settings.remove(
            Settings(
                {
                    "elements": {
                        "H": {
                            "mass": 1.008,
                            "metal": False,
                            "common_isotopes": {},
                            "properties": ["gas", "dimer"],
                        },
                        "O": None,
                        "Fe": {"name": "Iron", "num": 26, "mass": 55.845, "metal": True, "properties": ["common"]},
                    },
                    False: {},
                    2: {1: "one", 2.0: "two", 3.0001: "three"},
                }
            )
        )

        assert nested_settings["elements"] == {
            "Fe": {"properties": []},
            "H": {"name": "Hydrogen", "num": 1, "properties": []},
        }
        assert False not in nested_settings

    def test_settings_difference_returns_keys_not_in_other(self, nested_settings):
        diff = nested_settings.difference(
            Settings(
                {
                    "elements": {
                        "H": {
                            "mass": 1.008,
                            "metal": False,
                            "common_isotopes": {},
                            "properties": ["gas", "dimer"],
                        },
                        "O": None,
                        "Fe": {"name": "Iron", "num": 26, "mass": 55.845, "metal": True, "properties": ["common"]},
                    },
                    False: {},
                    2: {1: "one", 2.0: "two", 3.0001: "three"},
                }
            )
        )

        assert diff["elements"] == {
            "Fe": {"properties": []},
            "H": {"name": "Hydrogen", "num": 1, "properties": []},
        }
        assert False not in diff

    def test_settings_dictionary_equivalent_methods_case_insensitive(self, nested_settings):
        # Variety of dictionary methods should behave as usual but with case-insensitivity
        assert nested_settings.get("Elements").get("FE").get("NAME") == "Iron"
        assert nested_settings.elements.Fe.pop("MaSS") == 55.845
        assert not nested_settings.elements.Fe.pop("liquid_at_rt", False)
        assert nested_settings.elements.popitem() == (
            "Fe",
            {"name": "Iron", "num": 26, "metal": True, "properties": ["common"]},
        )
        assert nested_settings.elements.o.setdefault("nAme") == "Oxygen"
        assert nested_settings.elements.o.setdefault("gas_at_rt", True)
        assert list(nested_settings.elements.keys()) == ["H", "O"]
        assert list(nested_settings.elements.H.common_isotopes[0].keys()) == ["name", "mass", "abundance"]
        assert list(nested_settings.elements.H.common_isotopes[0].values()) == ["H1", 1, 99.99]
        assert list(nested_settings.elements.H.common_isotopes[1].items()) == [
            ("name", "H2"),
            ("mass", 2),
            ("abundance", 0.01),
        ]

    def test_settings_as_dict_returns_nested_dictionary(self, nested_settings):
        # Dictionary should be accessible as expected (case-sensitive)
        nested_dict = nested_settings.as_dict()
        assert nested_dict["elements"]["H"]["name"] == "Hydrogen"
        assert nested_dict["elements"]["H"]["common_isotopes"][1]["mass"] == 2
        assert nested_dict[False]["s"] == "string"
        assert nested_dict[False][True] == "bool"
        assert nested_dict[False][42] == "int"
        assert nested_dict[False][42.99] == "float"
        assert nested_dict[1][1] == "one"
        assert nested_dict[1][2] == "two"
        assert nested_dict[1][2.00] == "two"
        assert nested_dict[1][3.0001] == "three"

    def test_settings_suppress_missing_raises_key_error(self, flat_settings, nested_settings):
        # Should raise key error if any intermediate keys missing
        with flat_settings.suppress_missing(), pytest.raises(KeyError):
            flat_settings.z.name = "Zebra"
        with nested_settings.suppress_missing(), pytest.raises(KeyError):
            nested_settings.elements.Zn.name = "Zinc"

    @pytest.mark.parametrize(
        "suppress_missing", [True, False], ids=["with_suppress_missing", "without_suppress_missing"]
    )
    def test_settings_contains_nested_as_expected(self, suppress_missing, nested_settings):
        assert nested_settings.contains_nested(("eleMENTS", "Fe", "NAME"))
        assert nested_settings.contains_nested([1, 2])
        assert nested_settings.contains_nested((False,))
        assert nested_settings.contains_nested(("elements", "H", "common_isotopes"))
        assert nested_settings.contains_nested({"elements": 1, "H": 2, "common_isotopes": 3})
        assert nested_settings.contains_nested(("elements", "H", "common_isotopes", 0, "name"))

        if suppress_missing:
            with pytest.raises(KeyError):
                nested_settings.contains_nested(("eleMENTS", "Zn", "NAME"), True)
            with pytest.raises(KeyError):
                nested_settings.contains_nested(("elements", "Fe", "num", 32), True)
        else:
            assert not nested_settings.contains_nested(("eleMENTS", "Zn", "NAME"))
            assert not nested_settings.contains_nested(("elements", "Fe", "num", 32))

        with pytest.raises(TypeError):
            nested_settings.contains_nested("elements")

    @pytest.mark.parametrize(
        "suppress_missing", [True, False], ids=["with_suppress_missing", "without_suppress_missing"]
    )
    def test_settings_get_nested_as_expected(self, suppress_missing, nested_settings):
        assert nested_settings.get_nested(("eleMENTS", "Fe", "NAME")) == "Iron"
        assert nested_settings.get_nested([1, 2]) == "two"
        assert nested_settings.get_nested((False,)) == {"s": "string", True: "bool", 42: "int", 42.99: "float"}
        assert nested_settings.get_nested(("elements", "H", "common_isotopes")) == [
            {"name": "H1", "mass": 1, "abundance": 99.99},
            {"name": "H2", "mass": 2, "abundance": 0.01},
        ]
        assert nested_settings.get_nested(("elements", "H", "common_isotopes", 0, "name")) == "H1"

        if suppress_missing:
            with pytest.raises(KeyError):
                nested_settings.get_nested(("eleMENTS", "Zn", "NAME"), True)
            with pytest.raises(KeyError):
                nested_settings.get_nested(("elements", "Fe", "num", 32), True)
        else:
            assert nested_settings.get_nested(("eleMENTS", "Zn", "NAME")) is None
            assert nested_settings.get_nested(("elements", "Fe", "num", 32), default=42) == 42

        with pytest.raises(TypeError):
            nested_settings.get_nested("elements")

    @pytest.mark.parametrize(
        "suppress_missing", [True, False], ids=["with_suppress_missing", "without_suppress_missing"]
    )
    def test_settings_set_nested_as_expected(self, suppress_missing, nested_settings):
        nested_settings.set_nested(("eleMENTS", "Fe", "NAME"), "Ferrum")
        assert nested_settings.elements.Fe.name == "Ferrum"

        nested_settings.set_nested([1, 2], "2")
        assert nested_settings.get_nested((1, 2)) == "2"

        nested_settings.set_nested((False,), {"s": "string", False: "bool", 43: "int"})
        assert nested_settings.get_nested((False,)) == {"s": "string", False: "bool", 43: "int"}

        nested_settings.set_nested(
            ("elements", "H", "common_isotopes"),
            [
                {"name": "H1", "mass": 1, "abundance": 99.999},
                {"name": "H_2", "mass": 2, "abundance": 0.001},
            ],
        )
        nested_settings.set_nested(("elements", "H", "common_isotopes", 0, "name"), "H_1")
        assert nested_settings.get_nested(("elements", "H", "common_isotopes")) == [
            {"name": "H_1", "mass": 1, "abundance": 99.999},
            {"name": "H_2", "mass": 2, "abundance": 0.001},
        ]

        if suppress_missing:
            with pytest.raises(KeyError):
                nested_settings.set_nested(("eleMENTS", "Zn", "NAME"), "Zinc", True)
        else:
            nested_settings.set_nested(("eleMENTS", "Zn", "NAME"), "Zinc")
            assert nested_settings.elements.Zn.name == "Zinc"

    @pytest.mark.parametrize(
        "suppress_missing", [True, False], ids=["with_suppress_missing", "without_suppress_missing"]
    )
    def test_settings_pop_nested_as_expected(self, suppress_missing, nested_settings):
        assert nested_settings.pop_nested(("eleMENTS", "Fe", "NAME")) == "Iron"
        assert nested_settings.pop_nested((1, 2)) == "two"
        assert nested_settings.pop_nested((False,)) == {"s": "string", True: "bool", 42: "int", 42.99: "float"}
        assert nested_settings.pop_nested(("elements", "H", "common_isotopes", 0, "name")) == "H1"
        assert nested_settings.pop_nested(("elements", "H", "common_isotopes")) == [
            {"mass": 1, "abundance": 99.99},
            {"name": "H2", "mass": 2, "abundance": 0.01},
        ]

        if suppress_missing:
            with pytest.raises(KeyError):
                nested_settings.pop_nested(("eleMENTS", "Zn", "NAME"), True)
            with pytest.raises(KeyError):
                nested_settings.pop_nested(("elements", "Fe", "num", 32), True)
        else:
            assert nested_settings.pop_nested(("eleMENTS", "Zn", "NAME")) is None
            assert nested_settings.pop_nested(("elements", "Fe", "num", 32), default=42) == 42

        with pytest.raises(TypeError):
            nested_settings.pop_nested("elements")

    def test_settings_nested_keys(self, nested_settings):
        sett = Settings()
        sett.elements = nested_settings.elements
        sett.empty  # Add empty branches
        sett.elements.Fe.empty
        sett.elements.empty_list = []
        sett.elements.half_empty = [Settings(), Settings({"k": "v", "empty": ""})]

        keys = list(sett.nested_keys(flatten_list=False))
        assert keys == [
            ("elements",),
            ("elements", "half_empty"),
            ("elements", "H"),
            ("elements", "H", "name"),
            ("elements", "H", "num"),
            ("elements", "H", "mass"),
            ("elements", "H", "common_isotopes"),
            ("elements", "H", "properties"),
            ("elements", "O"),
            ("elements", "O", "name"),
            ("elements", "O", "num"),
            ("elements", "O", "mass"),
            ("elements", "O", "common_isotopes"),
            ("elements", "Fe"),
            ("elements", "Fe", "name"),
            ("elements", "Fe", "num"),
            ("elements", "Fe", "mass"),
            ("elements", "Fe", "metal"),
            ("elements", "Fe", "properties"),
        ]
        assert all([sett.contains_nested(k) for k in keys])

        keys = list(sett.nested_keys(flatten_list=False, include_empty=True))
        assert keys == [
            ("elements",),
            ("elements", "empty_list"),
            ("elements", "half_empty"),
            ("elements", "H"),
            ("elements", "H", "name"),
            ("elements", "H", "num"),
            ("elements", "H", "mass"),
            ("elements", "H", "metal"),
            ("elements", "H", "common_isotopes"),
            ("elements", "H", "properties"),
            ("elements", "O"),
            ("elements", "O", "name"),
            ("elements", "O", "num"),
            ("elements", "O", "mass"),
            ("elements", "O", "metal"),
            ("elements", "O", "common_isotopes"),
            ("elements", "Fe"),
            ("elements", "Fe", "name"),
            ("elements", "Fe", "num"),
            ("elements", "Fe", "mass"),
            ("elements", "Fe", "metal"),
            ("elements", "Fe", "properties"),
            ("elements", "Fe", "empty"),
            ("empty",),
        ]
        assert all([sett.contains_nested(k) for k in keys])

        keys = list(sett.nested_keys())
        assert keys == [
            ("elements",),
            ("elements", "H"),
            ("elements", "H", "name"),
            ("elements", "H", "num"),
            ("elements", "H", "mass"),
            ("elements", "H", "common_isotopes"),
            ("elements", "H", "common_isotopes", 0),
            ("elements", "H", "common_isotopes", 0, "name"),
            ("elements", "H", "common_isotopes", 0, "mass"),
            ("elements", "H", "common_isotopes", 0, "abundance"),
            ("elements", "H", "common_isotopes", 1),
            ("elements", "H", "common_isotopes", 1, "name"),
            ("elements", "H", "common_isotopes", 1, "mass"),
            ("elements", "H", "common_isotopes", 1, "abundance"),
            ("elements", "H", "properties"),
            ("elements", "H", "properties", 0),
            ("elements", "H", "properties", 1),
            ("elements", "O"),
            ("elements", "O", "name"),
            ("elements", "O", "num"),
            ("elements", "O", "mass"),
            ("elements", "O", "common_isotopes"),
            ("elements", "O", "common_isotopes", 0),
            ("elements", "O", "common_isotopes", 0, "name"),
            ("elements", "O", "common_isotopes", 0, "mass"),
            ("elements", "O", "common_isotopes", 0, "abundance"),
            ("elements", "O", "common_isotopes", 1),
            ("elements", "O", "common_isotopes", 1, "name"),
            ("elements", "O", "common_isotopes", 1, "mass"),
            ("elements", "O", "common_isotopes", 1, "abundance"),
            ("elements", "Fe"),
            ("elements", "Fe", "name"),
            ("elements", "Fe", "num"),
            ("elements", "Fe", "mass"),
            ("elements", "Fe", "metal"),
            ("elements", "Fe", "properties"),
            ("elements", "Fe", "properties", 0),
            ("elements", "half_empty"),
            ("elements", "half_empty", 1),
            ("elements", "half_empty", 1, "k"),
        ]
        assert all([sett.contains_nested(k) for k in keys])

        keys = list(sett.nested_keys(include_empty=True))
        assert keys == [
            ("elements",),
            ("elements", "empty_list"),
            ("elements", "H"),
            ("elements", "H", "name"),
            ("elements", "H", "num"),
            ("elements", "H", "mass"),
            ("elements", "H", "metal"),
            ("elements", "H", "common_isotopes"),
            ("elements", "H", "common_isotopes", 0),
            ("elements", "H", "common_isotopes", 0, "name"),
            ("elements", "H", "common_isotopes", 0, "mass"),
            ("elements", "H", "common_isotopes", 0, "abundance"),
            ("elements", "H", "common_isotopes", 1),
            ("elements", "H", "common_isotopes", 1, "name"),
            ("elements", "H", "common_isotopes", 1, "mass"),
            ("elements", "H", "common_isotopes", 1, "abundance"),
            ("elements", "H", "properties"),
            ("elements", "H", "properties", 0),
            ("elements", "H", "properties", 1),
            ("elements", "O"),
            ("elements", "O", "name"),
            ("elements", "O", "num"),
            ("elements", "O", "mass"),
            ("elements", "O", "metal"),
            ("elements", "O", "common_isotopes"),
            ("elements", "O", "common_isotopes", 0),
            ("elements", "O", "common_isotopes", 0, "name"),
            ("elements", "O", "common_isotopes", 0, "mass"),
            ("elements", "O", "common_isotopes", 0, "abundance"),
            ("elements", "O", "common_isotopes", 1),
            ("elements", "O", "common_isotopes", 1, "name"),
            ("elements", "O", "common_isotopes", 1, "mass"),
            ("elements", "O", "common_isotopes", 1, "abundance"),
            ("elements", "Fe"),
            ("elements", "Fe", "name"),
            ("elements", "Fe", "num"),
            ("elements", "Fe", "mass"),
            ("elements", "Fe", "metal"),
            ("elements", "Fe", "properties"),
            ("elements", "Fe", "properties", 0),
            ("elements", "Fe", "empty"),
            ("elements", "half_empty"),
            ("elements", "half_empty", 0),
            ("elements", "half_empty", 1),
            ("elements", "half_empty", 1, "k"),
            ("elements", "half_empty", 1, "empty"),
            ("empty",),
        ]
        assert all([sett.contains_nested(k) for k in keys])

    def test_settings_block_keys(self, nested_settings):
        sett = Settings()
        sett.elements = nested_settings.elements
        sett.empty  # Add empty branches
        sett.elements.Fe.empty
        sett.elements.empty_list = []
        sett.elements.half_empty = [Settings(), Settings({"k": "v", "empty": ""})]

        keys = list(sett.block_keys(flatten_list=False))
        assert keys == [("elements",), ("elements", "H"), ("elements", "O"), ("elements", "Fe")]
        assert all([sett.contains_nested(k) for k in keys])

        keys = list(sett.block_keys(flatten_list=False, include_empty=True))
        assert keys == [
            ("elements",),
            ("elements", "H"),
            ("elements", "O"),
            ("elements", "Fe"),
            ("elements", "Fe", "empty"),
            ("empty",),
        ]
        assert all([sett.contains_nested(k) for k in keys])

        keys = list(sett.block_keys())
        assert keys == [
            ("elements",),
            ("elements", "H"),
            ("elements", "H", "common_isotopes"),
            ("elements", "H", "common_isotopes", 0),
            ("elements", "H", "common_isotopes", 1),
            ("elements", "H", "properties"),
            ("elements", "O"),
            ("elements", "O", "common_isotopes"),
            ("elements", "O", "common_isotopes", 0),
            ("elements", "O", "common_isotopes", 1),
            ("elements", "Fe"),
            ("elements", "Fe", "properties"),
            ("elements", "half_empty"),
            ("elements", "half_empty", 1),
        ]
        assert all([sett.contains_nested(k) for k in keys])

        keys = list(sett.block_keys(include_empty=True))
        assert keys == [
            ("elements",),
            ("elements", "H"),
            ("elements", "H", "common_isotopes"),
            ("elements", "H", "common_isotopes", 0),
            ("elements", "H", "common_isotopes", 1),
            ("elements", "H", "properties"),
            ("elements", "O"),
            ("elements", "O", "common_isotopes"),
            ("elements", "O", "common_isotopes", 0),
            ("elements", "O", "common_isotopes", 1),
            ("elements", "Fe"),
            ("elements", "Fe", "properties"),
            ("elements", "Fe", "empty"),
            ("elements", "half_empty"),
            ("elements", "half_empty", 0),
            ("elements", "half_empty", 1),
            ("empty",),
        ]
        assert all([sett.contains_nested(k) for k in keys])

    def test_settings_compare_added_removed_and_modified(self, nested_settings):
        no_diff = nested_settings.compare(nested_settings)

        assert no_diff == {"added": {}, "removed": {}, "modified": {}}

        other = Settings(
            {
                "H": {
                    "name": "Hydrogen",
                    "num": 1,
                    "mass": 1.00799,
                    "common_isotopes": [
                        {"name": "H1", "mass": 1, "abundance": 99.99},
                        {"name": "H2", "mass": 2, "abundance": 0.01},
                    ],
                    "properties": ["gas", "dimer", "combustible"],
                },
                "O": {
                    "name": "Oxygen",
                    "num": 8,
                    "mass": 15.9999,
                    "metal": False,
                    "colour": {"liquid": "blue", "gas": "colourless"},
                    "common_isotopes": [
                        {"name": "O16", "mass": 16, "abundance": 99.8},
                        {"name": "O18", "mass": 18, "abundance": 0.02},
                    ],
                },
                "C": {"name": "Carbon"},
            }
        )

        diff = nested_settings["elements"].compare(other)

        assert diff["added"] == {
            ("Fe", "mass"): 55.845,
            ("Fe", "metal"): True,
            ("Fe", "name"): "Iron",
            ("Fe", "num"): 26,
            ("Fe", "properties", 0): "common",
            ("H", "metal"): False,
        }

        assert diff["removed"] == {
            ("C", "name"): "Carbon",
            ("H", "properties", 2): "combustible",
            ("O", "colour", "gas"): "colourless",
            ("O", "colour", "liquid"): "blue",
        }

        assert diff["modified"] == {
            ("H", "mass"): (1.008, 1.00799),
            ("O", "mass"): (15.999, 15.9999),
        }

    def test_settings_flatten_as_expected(self, nested_settings, extra_nested_settings):
        # Flatten, with a case of not flattening lists

        settings = nested_settings.elements.merge(extra_nested_settings.elements)
        settings_flattened = settings.flatten()
        settings_flattened_not_lists = settings.flatten(False)
        assert settings_flattened[("Fe", "mass")] == 55.845
        assert settings_flattened[("Fe", "properties", 0)] == "common"
        assert settings_flattened[("Fe", "common_isotopes", 1, "name")] == "Fe56"
        assert settings_flattened_not_lists[("Fe", "metal")]
        assert settings_flattened_not_lists[("Fe", "properties")] == ["common"]
        assert settings_flattened_not_lists[("O", "common_isotopes")][0] == Settings(
            {"name": "O16", "mass": 16, "abundance": 99.8}
        )

    def test_settings_unflatten_as_expected(self):
        # Unflatten, with a case of not unflattening lists
        settings_flattened = Settings({"a": "b", ("l", 0): "m", ("l", 1): "n", ("x", "y", 0): {"xy": "z"}})
        settings = settings_flattened.unflatten()
        assert settings.a == "b"
        assert settings.l[0] == "m"
        assert settings.l[1] == "n"
        assert settings.x.y[0].xy == "z"

        settings_flattened_with_lists = Settings({"a": "b", "l": ["m", "n"], ("x", "y"): [{"xy": "z"}]})
        settings_with_lists = settings_flattened_with_lists.unflatten()
        assert settings_with_lists.a == "b"
        assert settings_with_lists.l[0] == "m"
        assert settings_with_lists.l[1] == "n"
        assert settings_with_lists.x.y[0].xy == "z"

    def test_settings_can_roundtrip_to_and_from_flattened_and_unflattened_forms(
        self, flat_settings, nested_settings, extra_nested_settings
    ):
        # Can successfully move between flatten and unflatten, also with unflatten list option
        def roundtrip_and_assert(settings):
            settings_1 = settings.copy()
            settings_flattened_1 = settings_1.flatten()
            settings_2 = settings_flattened_1.unflatten()
            settings_flattened_2 = settings_2.flatten()

            assert settings_1 == settings
            assert settings_2 == settings
            assert settings_flattened_1 == settings_flattened_2

            settings_3 = settings.copy()
            settings_flattened_3 = settings_3.flatten(False)
            settings_4 = settings_flattened_3.unflatten(False)
            settings_flattened_4 = settings_4.flatten(False)

            assert settings_3 == settings_1
            assert settings_4 == settings_1
            assert settings_flattened_3 == settings_flattened_4

        roundtrip_and_assert(flat_settings)
        roundtrip_and_assert(nested_settings.elements)
        roundtrip_and_assert(nested_settings.elements.merge(extra_nested_settings.elements))

    def test_settings_string_representation(self, nested_settings, extra_nested_settings):
        elements = nested_settings.elements.soft_update(extra_nested_settings)
        actual = str(elements)
        expected = """H: 	
  name: 	Hydrogen
  num: 	1
  mass: 	1.008
  metal: 	False
  common_isotopes: 	[name: 	H1
                   	mass: 	1
                   	abundance: 	99.99
                   	, name: 	H2
                   	mass: 	2
                   	abundance: 	0.01
                   	]
  properties: 	['gas', 'dimer']
O: 	
  name: 	Oxygen
  num: 	8
  mass: 	15.999
  metal: 	False
  common_isotopes: 	[name: 	O16
                   	mass: 	16
                   	abundance: 	99.8
                   	, name: 	O18
                   	mass: 	18
                   	abundance: 	0.02
                   	]
Fe: 	
   name: 	Iron
   num: 	26
   mass: 	55.845
   metal: 	True
   properties: 	['common']
elements: 	
         Fe: 	
            common_isotopes: 	[name: 	Fe54
                             	mass: 	54
                             	abundance: 	5.85
                             	, name: 	Fe56
                             	mass: 	56
                             	abundance: 	91.8
                             	, name: 	Fe57
                             	mass: 	57
                             	abundance: 	2.12
                             	, name: 	Fe58
                             	mass: 	58
                             	abundance: 	0.28
                             	]
            gas_at_rt: 	False
            mass: 	55.846
            properties: 	['forms oxides']
         He: 	
            gas_at_rt: 	True
            mass: 	4.003
            metal: 	False
            name: 	Helium
            num: 	2
""".replace(
            "\r\n", "\n"
        )
        assert actual == expected

    def test_settings_dir_also_returns_dynamic_attributes(self, flat_settings):
        flat_settings["f"] = "fly"
        attributes = flat_settings.__dir__()

        assert "a" in attributes
        assert "B" in attributes
        assert "c " not in attributes
        assert "_d_" in attributes
        assert "__d__" in attributes
        assert "3" not in attributes
        assert "f" in attributes


class TestConfigSettings(TestSettings):
    """
    Test suite for the config settings class.
    This should have all the same functionality of the base settings class, but with additional defined fields.
    """

    @pytest.fixture
    def flat_settings(self):
        """
        Simple flat settings
        """
        return ConfigSettings(self.get_flat_settings())

    @pytest.fixture
    def nested_settings(self):
        """
        Nested settings with various primitive types for keys
        """
        return ConfigSettings(self.get_nested_settings())

    @pytest.fixture
    def extra_nested_settings(self):
        """
        Extra nested settings which can be patched into nested settings
        """
        return ConfigSettings(self.get_extra_nested_settings())

    @pytest.fixture
    def nested_mutable_settings(self):
        """
        Nested settings with mutable keys
        """
        return ConfigSettings(self.get_nested_mutable_settings())

    @pytest.fixture(autouse=True)
    def mock_job_runner_and_manager(self):
        """
        Mock out the job runner and manager to avoid creating run directories etc.
        """
        with patch("scm.plams.core.jobrunner.JobRunner") as mock_jobrunner, patch(
            "scm.plams.core.jobmanager.JobManager"
        ) as mock_jobmanager:
            yield mock_jobrunner, mock_jobmanager

    def test_settings_string_representation(self, nested_settings, extra_nested_settings):
        elements = nested_settings.elements.soft_update(extra_nested_settings)
        actual = str(elements)
        expected = """H: 	
  name: 	Hydrogen
  num: 	1
  mass: 	1.008
  metal: 	False
  common_isotopes: 	[name: 	H1
                   	mass: 	1
                   	abundance: 	99.99
                   	, name: 	H2
                   	mass: 	2
                   	abundance: 	0.01
                   	]
  properties: 	['gas', 'dimer']
O: 	
  name: 	Oxygen
  num: 	8
  mass: 	15.999
  metal: 	False
  common_isotopes: 	[name: 	O16
                   	mass: 	16
                   	abundance: 	99.8
                   	, name: 	O18
                   	mass: 	18
                   	abundance: 	0.02
                   	]
Fe: 	
   name: 	Iron
   num: 	26
   mass: 	55.845
   metal: 	True
   properties: 	['common']
_explicit_init: 	False
daemon_threads: 	True
default_jobmanager: 	None
default_jobrunner: 	None
elements: 	
         Fe: 	
            common_isotopes: 	[name: 	Fe54
                             	mass: 	54
                             	abundance: 	5.85
                             	, name: 	Fe56
                             	mass: 	56
                             	abundance: 	91.8
                             	, name: 	Fe57
                             	mass: 	57
                             	abundance: 	2.12
                             	, name: 	Fe58
                             	mass: 	58
                             	abundance: 	0.28
                             	]
            gas_at_rt: 	False
            mass: 	55.846
            properties: 	['forms oxides']
         He: 	
            gas_at_rt: 	True
            mass: 	4.003
            metal: 	False
            name: 	Helium
            num: 	2
erase_workdir: 	False
ignore_failure: 	True
init: 	False
job: 	
    pickle: 	True
    pickle_protocol: 	-1
    keep: 	all
    save: 	all
    runscript: 	
              shebang: 	#!/bin/sh
              stdout_redirect: 	False
    link_files: 	True
jobmanager: 	
           counter_len: 	3
           hashing: 	input
           remove_empty_directories: 	True
log: 	
    file: 	5
    stdout: 	3
    csv: 	7
    time: 	True
    date: 	True
preview: 	False
saferun: 	
        repeat: 	10
        delay: 	1
sleepstep: 	5
""".replace(
            "\r\n", "\n"
        )
        assert actual == expected

    def test_new_has_correct_defaults_and_nested_types(self):
        config = ConfigSettings()
        assert_config_as_expected(config, init=False, explicit_init=False)

    def test_property_setters_and_getters_equivalent_to_dict_get_and_set_items(self):
        config = ConfigSettings()

        # Getter equal to dict get item
        def assert_equivalent_access(preview=False, stdout_redirect=False, stdout=3):
            assert config.preview == config["preview"] == preview
            assert (
                config.job.runscript.stdout_redirect == config["job"]["runscript"]["stdout_redirect"] == stdout_redirect
            )
            assert config.log.stdout == config["log"]["stdout"] == stdout

        assert_equivalent_access()

        # Modify using property setters
        config.preview = True
        config.job.runscript.stdout_redirect = True
        config.log.stdout = 1

        assert_equivalent_access(True, True, 1)

        # Modify using dict set item
        config["preview"] = False
        config["job"]["runscript"]["stdout_redirect"] = False
        config["log"]["stdout"] = 7

        assert_equivalent_access(stdout=7)

    def test_copy_has_correct_defaults_and_nested_types(self):
        config = ConfigSettings()

        copy = config.copy()
        assert_config_as_expected(copy, init=False, explicit_init=False)

    def test_get_non_existing_item_returns_empty_settings(self):
        config = ConfigSettings()
        assert config.foo == Settings()

    def test_flatten_and_unflatten_coverts_to_settings(self):
        # Round trip between flattened and unflattened forms, should give the same structure but as normal settings object
        config = ConfigSettings()
        flattened = config.flatten()
        unflattened = flattened.unflatten()

        assert_config_as_expected(unflattened, init=False, explicit_init=False, verify_derived_types=False)

    def test_new_has_lazy_job_runner_and_manager(self, mock_job_runner_and_manager):
        # Verify initially that the job runner/manager are not called
        mock_jobrunner, mock_jobmanager = mock_job_runner_and_manager
        config = ConfigSettings()

        assert not mock_jobrunner.called
        assert not mock_jobmanager.called

        # Force initialisation
        _ = config.default_jobrunner
        _ = config.default_jobmanager

        # Verify that the components are called with the default arguments
        assert mock_jobrunner.called
        assert mock_jobmanager.called
        assert mock_jobrunner.call_args.args == ()
        assert mock_jobmanager.call_args.args == (
            {"counter_len": 3, "hashing": "input", "remove_empty_directories": True},
        )

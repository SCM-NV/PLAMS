#!/usr/bin/env amspython
# coding: utf-8

# ## Task, Engine and Property Blocks

from scm.plams import Settings, AMSJob, init

# this line is not required in AMS2025+
init()


# Helper function to print the settings in the AMS input format style:


def print_input(settings):
    print(AMSJob(settings=settings).get_input())


# Tasks can be specified in the settings under the `input.AMS.Task` key:

go_settings = Settings()
go_settings.input.AMS.Task = "GeometryOptimization"
go_settings.input.AMS.GeometryOptimization.Convergence.Gradients = 1e-5
print_input(go_settings)


# Properties can be specified under the `input.AMS.Properties` key:

nm_settings = Settings()
nm_settings.input.ams.Properties.NormalModes = "Yes"
print_input(nm_settings)


# Engine settings can be specified under the `input.<Engine>` key, for the engine of interest:

lda_settings = Settings()
lda_settings.input.ADF.Basis.Type = "DZP"
print_input(lda_settings)


# ## Combining Settings

# Settings objects can also be combined for easy reuse and job setup. Settings can be merged using the `+` and `+=` operators.

settings = go_settings + lda_settings + nm_settings
print_input(settings)


# Note however that this merge is a "soft" update, so values of existing keys will not be overwritten:

pbe_settings = Settings()
pbe_settings.input.ADF.Basis.Type = "TZP"
pbe_settings.input.ADF.xc.gga = "pbe"
settings += pbe_settings  # adds "pbe" but will not set "TZP", because a Basis.Type already exists in settings
print_input(settings)


# To achieve "hard update" behaviour, the `update` method can be used, which overwrites existing keys:

settings.update(pbe_settings)  # will also set "TZP"
print_input(settings)


# ## Subtracting settings
#
# In AMS2025+, settings can also be removed using the `-` and `-=` operators:


def check_ams_version():
    try:
        from scm.plams import __version__

        return __version__ >= "2024.2"
    except ImportError:
        return False


is_ams_2025_or_higher = check_ams_version()


if is_ams_2025_or_higher:
    settings -= nm_settings
    print_input(settings)


# ## Comparison

# In AMS2025+, two settings objects can be compared to check the differences between them. The result will show the nested key and value of any added, removed and modified entries.

from pprint import pprint

if is_ams_2025_or_higher:
    settings1 = go_settings + lda_settings + nm_settings
    settings2 = Settings()
    settings2.input.AMS.Task = "SinglePoint"
    settings2.input.DFTB.Model = "GFN1-xTB"
    comparison = settings2.compare(settings1)  # dictionary with keys 'added', 'modified', and 'removed'
    pprint(comparison)


# ## Repeated blocks as lists, Hybrid engine
#
# Multiple values in a settings block can be configured using a list:

hybrid_settings = go_settings.copy()
hybrid_settings.input.Hybrid.Energy.Term = []
for i in range(5):
    factor = (-1) ** (i % 2) * 1.0
    region = "*" if i == 0 else "one" if i < 3 else "two"
    engine_id = "adf-lda" if i == 0 or factor == -1 else "adf-gga"
    term = Settings({"Factor": factor, "Region": region, "EngineID": engine_id})
    hybrid_settings.input.Hybrid.Energy.Term.append(term)
hybrid_settings.input.Hybrid.Engine = [lda_settings.input.ADF.copy(), pbe_settings.input.ADF.copy()]
hybrid_settings.input.Hybrid.Engine[0]._h = "ADF adf-lda"
hybrid_settings.input.Hybrid.Engine[1]._h = "ADF adf-gga"
print_input(hybrid_settings)


# Note also in the example below, the use of the special `_h` "header" key, which can be used to add data to the header line for a block.

# ## Nested Keys

# It can be useful to access values from a Settings object using "nested" keys, available in AMS2025+. These are tuples of keys, where each successive element of the tuple corresponds to a further layer in the settings. Lists are flattened so their elements can be accessed with the corresponding index.
#
# These options are mainly useful if you programatically need to compare or set different settings. If you manually create a Settings object it is easier to just use the normal dot notation shown in the previous example.

if is_ams_2025_or_higher:
    print(list(hybrid_settings.nested_keys()))


if is_ams_2025_or_higher:
    print(hybrid_settings.get_nested(("input", "AMS", "Task")))


if is_ams_2025_or_higher:
    if hybrid_settings.contains_nested(("input", "Hybrid", "Engine", 0)):
        hybrid_settings.set_nested(("input", "Hybrid", "Engine", 0, "NumericalQuality"), "Good")
    print_input(hybrid_settings)

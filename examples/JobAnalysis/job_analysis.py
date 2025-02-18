#!/usr/bin/env amspython
# coding: utf-8

# ## Create Example Jobs

# To begin with, create a variety of AMS jobs with different settings, engines and calculation types.

from scm.plams import from_smiles, AMSJob, PlamsError, Settings, Molecule, Atom
from scm.libbase import UnifiedChemicalSystem as ChemicalSystem
from scm.input_classes.drivers import AMS
from scm.input_classes.engines import DFTB
from scm.utils.conversions import plams_molecule_to_chemsys


def example_job_dftb(smiles, task, use_chemsys=False):
    # Generate molecule from smiles
    mol = from_smiles(smiles)
    if use_chemsys:
        mol = plams_molecule_to_chemsys(mol)

    # Set up calculation settings using PISA
    sett = Settings()
    sett.runscript.nproc = 1
    driver = AMS()
    driver.Task = task
    driver.Engine = DFTB()
    sett.input = driver
    return AMSJob(molecule=mol, settings=sett, name="dftb")


def example_job_adf(smiles, task, basis, gga=None, use_chemsys=False):
    # Generate molecule from smiles
    mol = from_smiles(smiles)
    if use_chemsys:
        mol = plams_molecule_to_chemsys(mol)

    # Set up calculation settings using standard settings
    sett = Settings()
    sett.runscript.nproc = 1
    sett.input.AMS.Task = task
    sett.input.ADF.Basis.Type = basis
    if gga:
        sett.input.ADF.XC.GGA = gga
    return AMSJob(molecule=mol, settings=sett, name="adf")


def example_job_neb(iterations, use_chemsys=False):
    # Set up molecules
    main_molecule = Molecule()
    main_molecule.add_atom(Atom(symbol="C", coords=(0, 0, 0)))
    main_molecule.add_atom(Atom(symbol="N", coords=(1.18, 0, 0)))
    main_molecule.add_atom(Atom(symbol="H", coords=(2.196, 0, 0)))
    final_molecule = main_molecule.copy()
    final_molecule.atoms[1].x = 1.163
    final_molecule.atoms[2].x = -1.078

    mol = {"": main_molecule, "final": final_molecule}

    if use_chemsys:
        mol = {k: plams_molecule_to_chemsys(v) for k, v in mol.items()}

    # Set up calculation settings
    sett = Settings()
    sett.runscript.nproc = 1
    sett.input.ams.Task = "NEB"
    sett.input.ams.NEB.Images = 9
    sett.input.ams.NEB.Iterations = iterations
    sett.input.DFTB

    return AMSJob(molecule=mol, settings=sett, name="neb")


# Now, run a selection of them.

from scm.plams import config, JobRunner

config.default_jobrunner = JobRunner(parallel=True, maxthreads=8)

smiles = ["CC", "C", "O", "CO"]
tasks = ["SinglePoint", "GeometryOptimization"]
engines = ["DFTB", "ADF"]
jobs = []
for i, s in enumerate(smiles):
    for j, t in enumerate(tasks):
        job_dftb = example_job_dftb(s, t, use_chemsys=i % 2)
        job_adf1 = example_job_adf(s, t, "DZ", use_chemsys=True)
        job_adf2 = example_job_adf(s, t, "TZP", "PBE")
        jobs += [job_dftb, job_adf1, job_adf2]

job_neb1 = example_job_neb(10)
job_neb2 = example_job_neb(100, use_chemsys=True)
jobs += [job_neb1, job_neb2]

for j in jobs:
    j.run()


# ## Job Analysis

# ### Adding and Loading Jobs
#
# Jobs can be loaded by passing job objects directly, or loading from a path.

from scm.plams import JobAnalysis


ja = JobAnalysis(jobs=jobs[:10], paths=[j.path for j in jobs[10:-2]])


# Jobs can also be added or removed after initialization.

ja.add_job(jobs[-2]).load_job(jobs[-1].path).display_table()


# ### Adding and Removing Fields

# A range of common standard fields can be added with the `add_standard_field(s)` methods.
# Custom fields can also be added with the `add_field` method, by defining a field key, value accessor and optional arguments like display name and value formatting.
# Fields can be removed by calling `remove_field` with the corresponding field key.

(
    ja.remove_field("Path")
    .add_standard_fields(["Formula", "Smiles", "CPUTime", "SysTime"])
    .add_settings_input_fields()
    .add_field("Energy", lambda j: j.results.get_energy(unit="kJ/mol"), display_name="Energy [kJ/mol]", fmt=".2f")
    .display_table(max_rows=5)
)


# In addition to the fluent syntax, both dictionary and dot syntaxes are also supported for adding and removing fields.

import numpy as np

ja["AtomType"] = lambda j: [at.symbol for at in j.results.get_main_molecule()]
ja.Charge = lambda j: j.results.get_charges()
ja.AtomCoords = lambda j: [np.array(at.coords) for at in j.results.get_main_molecule()]

del ja["Check"]
del ja.SysTime

ja.display_table(max_rows=5, max_col_width=30)


# ### Processing Data

# Once an initial analysis has been created, the data can be further processed, depending on the use case.
# For example, to inspect the difference between failed and successful jobs, jobs can be filtered down and irrelevant fields removed.

ja_neb = (
    ja.copy()
    .filter_jobs(lambda data: data["InputAmsTask"] == "NEB")
    .remove_field("AtomCoords")
    .remove_uniform_fields(ignore_empty=True)
)

ja_neb.display_table()


# Another use case may be to analyze the results from one or more jobs.
# For this, it can be useful to utilize the `expand` functionality to convert job(s) to multiple rows.
# During this process, fields selected for expansion will have their values extracted into individual rows, whilst other fields have their values duplicated.

ja_adf = (
    ja.copy()
    .filter_jobs(
        lambda data: data["InputAmsTask"] == "GeometryOptimization"
        and data["InputAdfBasisType"] is not None
        and data["Smiles"] == "O"
    )
    .expand_field("AtomType")
    .expand_field("Charge")
    .expand_field("AtomCoords")
    .remove_uniform_fields()
)

ja_adf.display_table()


# For more nested values, the depth of expansion can also be selected to further flatten the data.

(
    ja_adf.add_field("Coord", lambda j: [("x", "y", "z") for _ in j.results.get_main_molecule()], expansion_depth=2)
    .expand_field("AtomCoords", depth=2)
    .display_table()
)


# Expansion can be undone with the corresponding `collapse` method.
#
# Fields can be also further filtered, modified or reordered to customize the analysis.

(
    ja_adf.collapse_field("AtomCoords")
    .collapse_field("Coord")
    .filter_fields(lambda vals: all([not isinstance(v, list) for v in vals]))  # remove arrays
    .remove_field("Name")
    .format_field("CPUTime", ".2f")
    .format_field("Charge", ".4f")
    .rename_field("InputAdfBasisType", "Basis")
    .reorder_fields(["AtomType", "Charge", "Energy"])
    .display_table()
)


# ### Extracting Analysis Data

# Analysis data can be extracted in a variety of ways.
#
# As has been demonstrated, a visual representation of the table can be easily generated using the `to_table` method (or `display_table` in a notebook).
# The format can be selected as markdown, html or rst. This will return the data with the specified display names and formatting.

print(ja_adf.to_table(fmt="rst"))


# Alternatively, raw data can be retrieved via the `get_analysis` method, which returns a dictionary of analysis keys to values.

print(ja_adf.get_analysis())


# Data can also be easily written to a csv file using `to_csv_file`, to be exported to another program.

csv_name = "./tmp.csv"
ja_adf.to_csv_file(csv_name)

with open(csv_name) as csv:
    print(csv.read())


# Finally, for more complex data analysis, the results can be converted to a [pandas](https://pandas.pydata.org) dataframe. This is recommended for more involved data manipulations, and can be installed using amspackages i.e. using the command: `"${AMSBIN}/amspackages" install pandas`.

try:
    import pandas

    df = ja_adf.to_dataframe()
    print(df)

except ImportError:

    print(
        "Pandas not available. Please install with amspackages to run this example '${AMSBIN}/amspackages install pandas'"
    )


# ### Additional Analysis Methods

# The `JobAnalysis` class does have some additional built in methods to aid with job analysis.
#
# For example, the `get_timeline` and `display_timeline` methods show pictorially when jobs started, how long they took to run and what their status is.
#
# This can be useful for visualising the dependencies of jobs. Here you can see that the first 8 jobs started running in parallel, due to the `maxthreads` constraint, and the remaining jobs waited before starting. Also that the penultimate job failed.

ja.display_timeline(fmt="rst")

#!/usr/bin/env amspython
# coding: utf-8

import itertools
import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.testing.decorators import image_comparison
from scm.plams.core.functions import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.mol.atom import Atom
from scm.plams.mol.molecule import Molecule
from scm.plams.recipes.md.trajectoryanalysis import AMSMSDJob, AMSMSDResults
from scm.plams.tools.plot import (
    get_correlation_xy,
    linear_fit_extrapolate_to_0,
    plot_band_structure,
    plot_phonons_band_structure,
    plot_phonons_dos,
    plot_correlation,
    plot_grid_molecules,
    plot_molecule,
    plot_molecule_counts,
    plot_msd,
    plot_work_function,
    plot_energy_landscape,
)
from test_helpers import skip_if_no_scm_inputs

# Use non-interactive backend for ci env
# Results are available for inspection in the result_images directory
matplotlib.use("Agg")


def test_linear_fit_extrapolate_to_0():
    fit_x, fit_y, slope, intercept = linear_fit_extrapolate_to_0([1.0, 2.0, 3.0], [3.0, 5.0, 7.0])

    assert slope == pytest.approx(2.0)
    assert intercept == pytest.approx(1.0)
    assert fit_x.tolist() == pytest.approx([1.0, 2.0, 3.0, 0.0])
    assert fit_y.tolist() == pytest.approx([3.0, 5.0, 7.0, 1.0])

    fit_x, fit_y, slope, intercept = linear_fit_extrapolate_to_0([0.0, 1.0, 2.0], [1.0, 3.0, 5.0])

    assert slope == pytest.approx(2.0)
    assert intercept == pytest.approx(1.0)
    assert fit_x.tolist() == pytest.approx([0.0, 1.0, 2.0])
    assert fit_y.tolist() == pytest.approx([1.0, 3.0, 5.0])
    assert fit_x.tolist().count(0.0) == 1


@pytest.fixture
def run_calculations():
    run_calculations = False  # Manual toggle whether to re-run AMS calculations
    return "AMSBIN" in os.environ and run_calculations


@pytest.fixture
def rkf_tools_plot(rkf_folder):
    return Path(rkf_folder / "tools_plot")


# ----------------------------------------------------------
# Testing plot_molecule
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_molecule"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=30,
)
def test_plot_molecule():
    plt.close("all")

    glycine = from_smiles("C(C(=O)O)N")
    glycine.from_array(
        [
            [1.52585228, -0.36095771, 0.17604921],
            [0.26794267, 0.21265392, -0.21199336],
            [-0.90913905, -0.53151974, 0.26695495],
            [-0.85195920, -1.53060464, 1.04179363],
            [-2.14637457, -0.16398621, -0.26887301],
            [2.25707210, 0.02537156, -0.47070660],
            [1.79020298, -0.01286013, 1.12713073],
            [0.28233609, 0.31718008, -1.32556763],
            [0.18775611, 1.26499430, 0.16051118],
            [-2.40368941, 0.77972858, -0.49529909],
        ]
    )

    plot_molecule(glycine, rotation=("90x,45y,0z"), keep_axis=True, figsize=(3, 2))


# ----------------------------------------------------------
# Testing plot_grid_molecules
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_grid_molecules"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=25,
)
def test_plot_grid_molecules():
    plt.close("all")
    ethanol = from_smiles("CCO")
    ethanol.from_array(
        [
            [0.88165321, -0.04478118, -0.01474324],
            [-0.58159753, -0.3761451, 0.05108011],
            [-1.35694986, 0.75702494, 0.1824042],
            [1.26504668, 0.17421359, 1.01224746],
            [1.01649295, 0.87054063, -0.60898906],
            [1.48536012, -0.89981469, -0.39366191],
            [-0.78791391, -0.98893872, 0.96478502],
            [-0.83504017, -1.00512403, -0.81678084],
            [-1.08705149, 1.51302455, -0.37634174],
        ]
    )

    glycine = from_smiles("C(C(=O)O)N")
    glycine.from_array(
        [
            [1.52585228, -0.36095771, 0.17604921],
            [0.26794267, 0.21265392, -0.21199336],
            [-0.90913905, -0.53151974, 0.26695495],
            [-0.85195920, -1.53060464, 1.04179363],
            [-2.14637457, -0.16398621, -0.26887301],
            [2.25707210, 0.02537156, -0.47070660],
            [1.79020298, -0.01286013, 1.12713073],
            [0.28233609, 0.31718008, -1.32556763],
            [0.18775611, 1.26499430, 0.16051118],
            [-2.40368941, 0.77972858, -0.49529909],
        ]
    )
    molecules = [ethanol, glycine]
    _, ax = plt.subplots(figsize=(3, 2))
    ax.axis("off")
    plot_grid_molecules(molecules, ax=ax)


def test_plot_grid_molecules_options_products():
    """since of the presence of global variables, working with rdkit is difficult. I test the different options products, if no errors are raised it is already good"""

    save_svg_path = Path("mols.svg")

    def svg_plot(molecules):
        _ = plot_grid_molecules(molecules, molsPerRow=2, save_svg_path=save_svg_path)

    def ax_plot(molecules):
        _, ax = plt.subplots()
        ax = plot_grid_molecules(molecules, molsPerRow=2, ax=ax)

    def pil_notebook_plot(molecules):
        _ = plot_grid_molecules(molecules, molsPerRow=2)

    def iterate_options(molecules, options):
        for x in options:
            if x == "svg":
                svg_plot(molecules)
            if x == "ax":
                ax_plot(molecules)
            if x == "notebook":
                pil_notebook_plot(molecules)

    molecules = [from_smiles("CCO"), from_smiles("NC")]
    for options_x in itertools.product(["svg", "ax", "notebook"], repeat=2):
        iterate_options(molecules, options_x)
    save_svg_path.unlink(missing_ok=True)


# ----------------------------------------------------------
# Testing plot_band_structure
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_band_structure"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=10,
)
def test_plot_band_structure(run_calculations, rkf_tools_plot):
    plt.close("all")

    d = 2.085
    mol = Molecule()
    mol.add_atom(Atom(symbol="Ni", coords=(0, 0, 0)))
    mol.add_atom(Atom(symbol="O", coords=(d, d, d)))
    mol.lattice = [[0.0, d, d], [d, 0.0, d], [d, d, 0.0]]

    s = Settings()
    s.input.ams.task = "SinglePoint"
    s.input.band.Unrestricted = "yes"
    s.input.band.XC.GGA = "BP86"
    s.input.band.Basis.Type = "DZ"
    s.input.band.NumericalQuality = "Basic"
    s.input.band.HubbardU.Enabled = "Yes"
    s.input.band.HubbardU.UValue = "0.6 0.0"
    s.input.band.HubbardU.LValue = "2 -1"
    s.input.band.BandStructure.Enabled = "Yes"

    if run_calculations:
        job = AMSJob(settings=s, molecule=mol, name="NiO")
        job.run()
    else:
        job = AMSJob.load_external(rkf_tools_plot / "NiO")

    ax = plot_band_structure(*job.results.get_band_structure(unit="eV"), zero="vbmax")
    ax.set_ylim(-10, 10)
    ax.set_ylabel("$E - E_{VBM}$ (eV)")
    ax.set_xlabel("Path")
    ax.set_title("NiO with DFT+U")


# ----------------------------------------------------------
# Testing plot_phonons_band_structure
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_phonons_band_structure"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=10,
)
def test_plot_phonons_band_structure(run_calculations, rkf_tools_plot):
    plt.close("all")

    d = 0.44625000
    mol = Molecule()
    mol.add_atom(Atom(symbol="C", coords=(d, d, d)))
    mol.add_atom(Atom(symbol="C", coords=(-d, -d, -d)))
    mol.lattice = [[0.0, 4 * d, 4 * d], [4 * d, 0.0, 4 * d], [4 * d, 4 * d, 0.0]]

    s = Settings()
    s.input.ams.task = "SinglePoint"
    s.input.ams.Properties.Phonons = "Yes"
    s.input.ams.Phonons.Method = "Numerical"
    s.input.ams.NumericalPhonons.SuperCell._1 = "2 0 0"
    s.input.ams.NumericalPhonons.SuperCell._2 = "0 2 0"
    s.input.ams.NumericalPhonons.SuperCell._3 = "0 0 2"
    s.input.ams.NumericalPhonons.AutomaticBZPath = "Yes"
    s.input.DFTB = Settings()

    if run_calculations:
        job = AMSJob(settings=s, molecule=mol, name="Diamond")
        job.run()
    else:
        job = AMSJob.load_external(rkf_tools_plot / "Diamond")

    ax = plot_phonons_band_structure(*job.results.get_phonons_band_structure(unit="cm^-1"))
    ax.set_ylabel("$E$ (cm$^{-1}$)")
    ax.set_xlabel("Path")
    ax.set_title("Diamond Phonons with DFTB")


# ----------------------------------------------------------
# Testing plot_phonons_dos
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_phonons_dos"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=12,
)
def test_plot_phonons_dos(run_calculations, rkf_tools_plot):
    plt.close("all")

    d = 0.44625000
    mol = Molecule()
    mol.add_atom(Atom(symbol="C", coords=(d, d, d)))
    mol.add_atom(Atom(symbol="C", coords=(-d, -d, -d)))
    mol.lattice = [[0.0, 4 * d, 4 * d], [4 * d, 0.0, 4 * d], [4 * d, 4 * d, 0.0]]

    s = Settings()
    s.input.ams.task = "SinglePoint"
    s.input.ams.Properties.Phonons = "Yes"
    s.input.ams.Phonons.Method = "Numerical"
    s.input.ams.NumericalPhonons.SuperCell._1 = "2 0 0"
    s.input.ams.NumericalPhonons.SuperCell._2 = "0 2 0"
    s.input.ams.NumericalPhonons.SuperCell._3 = "0 0 2"
    s.input.ams.NumericalPhonons.AutomaticBZPath = "Yes"
    s.input.DFTB = Settings()

    if run_calculations:
        job = AMSJob(settings=s, molecule=mol, name="Diamond")
        job.run()
    else:
        job = AMSJob.load_external(rkf_tools_plot / "Diamond")

    ax = plot_phonons_dos(*job.results.get_phonons_dos(unit="cm^-1"))
    ax.set_ylabel("DOS (1/cm$^{-1}$)")
    ax.set_xlabel("Energy (cm$^{-1}$)")
    ax.set_title("Diamond Phonons DOS with DFTB")


# ----------------------------------------------------------
# Testing plot_correlation & get_correlation_xy
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_correlation"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=12,
)
def test_plot_correlation(run_calculations, rkf_tools_plot):
    plt.close("all")

    glycine = from_smiles("C(C(=O)O)N")

    e1 = Settings()
    e1.input.GFNFF

    e2 = Settings()
    e2.input.DFTB.Model = "GFN1-xTB"

    sp = Settings()
    sp.input.ams.Task = "SinglePoint"
    sp.input.ams.Properties.Gradients = "Yes"

    if run_calculations:
        job1 = AMSJob(settings=sp + e1, name="glycine-engine1", molecule=glycine)
        job2 = AMSJob(settings=sp + e2, name="glycine-engine2", molecule=glycine)
        job1.run()
        job2.run()
    else:
        job1 = AMSJob.load_external(rkf_tools_plot / "glycine-engine1")
        job2 = AMSJob.load_external(rkf_tools_plot / "glycine-engine2")

    plot_correlation(job1, job2, section="AMSResults", variable="Gradients", file="engine")

    x, y = get_correlation_xy(job1, job2, section="AMSResults", variable="Gradients", file="engine")

    x0 = [
        -0.02618337,
        -0.02185398,
        -0.00569558,
        -0.00728149,
        -0.02743453,
        0.00760765,
        0.02840358,
        0.05185629,
        -0.02249499,
        -0.00930743,
        -0.04960230,
        0.03370857,
        0.00630104,
        -0.00342449,
        -0.00515473,
        0.01226365,
        0.01222411,
        -0.01690943,
        0.00674109,
        0.01271709,
        0.01819203,
        0.00787463,
        0.00516595,
        -0.01108063,
        -0.00154521,
        0.01105898,
        0.00377830,
        -0.01726648,
        0.00929288,
        -0.00195118,
    ]

    y0 = [
        -0.03408318,
        -0.01360583,
        -0.00908411,
        -0.00699156,
        -0.03623220,
        0.01994845,
        0.03482779,
        0.08117168,
        -0.04550755,
        0.00640144,
        -0.08007524,
        0.06008221,
        -0.01476529,
        -0.01611956,
        -0.01090110,
        0.01657241,
        0.01123652,
        -0.01511969,
        0.00813057,
        0.01138332,
        0.01821275,
        0.00715580,
        0.00602574,
        -0.01219420,
        0.00199320,
        0.00620757,
        0.00153412,
        -0.01924117,
        0.03000799,
        -0.00697089,
    ]

    assert np.allclose(x, x0, 1e-5)
    assert np.allclose(y, y0, 1e-5)


# ----------------------------------------------------------
# Testing plot_molecule_counts
# ----------------------------------------------------------
def test_plot_molecule_counts(rkf_tools_plot):
    plt.close("all")

    job = AMSJob.load_external(rkf_tools_plot / "md")
    frames, time_fs, counts = job.results.get_molecule_count_history(species=["H2O"])
    assert frames[:2].tolist() == [1, 2]
    assert time_fs[:2].tolist() == pytest.approx([0.0, 0.5])
    assert counts["H2O"][:2].tolist() == [16, 16]

    ax = plot_molecule_counts(job, species=["H2O"], time_unit="ps")

    assert ax.get_xlabel() == "Time (ps)"
    assert ax.get_ylabel() == "Molecule count"
    assert ax.lines[0].get_label() == "H2O"
    assert ax.lines[0].get_xdata().tolist()[:2] == pytest.approx([0.0, 0.0005])
    assert ax.lines[0].get_ydata().tolist()[:2] == [16, 16]

    ax = plot_molecule_counts(job, x_axis="frame")
    assert ax.get_xlabel() == "Frame"
    assert ax.lines[0].get_xdata().tolist()[:2] == [1, 2]


# ----------------------------------------------------------
# Testing plot_msd
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_msd"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=20,
)
def test_plot_msd(run_calculations, rkf_tools_plot, xyz_folder):
    plt.close("all")

    mol = Molecule(xyz_folder / "water_box.xyz")
    s = Settings()
    s.input.ams.Task = "MolecularDynamics"
    s.input.ReaxFF.ForceField = "Water2017.ff"
    s.input.ams.RNGSeed = "1 2 3 4 5 6 7 8 9"
    s.input.ams.MolecularDynamics.CalcPressure = "Yes"
    s.input.ams.MolecularDynamics.InitialVelocities.Temperature = 300
    s.input.ams.MolecularDynamics.Trajectory.SamplingFreq = 1
    s.input.ams.MolecularDynamics.TimeStep = 0.5
    s.input.ams.MolecularDynamics.NSteps = 200

    if run_calculations:
        s.runscript.nproc = 1
        os.environ["OMP_NUM_THREADS"] = "1"
        job = AMSJob(settings=s, molecule=mol, name="md")
        job.run()
        md_job = AMSMSDJob(job)
        md_job.run()
    else:
        # Cannot load the AMSMSDJob directly, so simulate running the job by loading the kf into the results
        job = AMSJob.load_external(rkf_tools_plot / "md/ams.rkf", settings=s)
        md_job = AMSMSDJob(job, name="msd")
        md_job.prerun()
        md_job.path = Path(rkf_tools_plot / "md/msd")
        results = AMSMSDResults(md_job)
        results.finished.set()
        results.done.set()
        results.collect()
        md_job.results = results

    plot_msd(md_job)


@image_comparison(
    baseline_images=["plot_msd_scm_inputs"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=10,
)
def test_plot_msd_with_scm_inputs(run_calculations, rkf_tools_plot, xyz_folder):
    skip_if_no_scm_inputs()

    from scm.inputs import Analysis

    plt.close("all")

    mol = Molecule(xyz_folder / "water_box.xyz")
    s = Settings()
    s.input.ams.Task = "MolecularDynamics"
    s.input.ReaxFF.ForceField = "Water2017.ff"
    s.input.ams.RNGSeed = "1 2 3 4 5 6 7 8 9"
    s.input.ams.MolecularDynamics.CalcPressure = "Yes"
    s.input.ams.MolecularDynamics.InitialVelocities.Temperature = 300
    s.input.ams.MolecularDynamics.Trajectory.SamplingFreq = 1
    s.input.ams.MolecularDynamics.TimeStep = 0.5
    s.input.ams.MolecularDynamics.NSteps = 200

    sets = Analysis()
    msd = Analysis.MeanSquareDisplacementBlock()
    msd.UseAllValues = True
    msd.Property = "Coords"
    msd.Atoms.Element = ["O"]
    sets.MeanSquareDisplacement = [msd]
    s_msd = Settings()
    s_msd.input = sets

    if run_calculations:
        s.runscript.nproc = 1
        os.environ["OMP_NUM_THREADS"] = "1"
        job = AMSJob(settings=s, molecule=mol, name="md")
        job.run()
        msd_job = AMSMSDJob(job, settings=s_msd)
        msd_job.run()
        msd_scm_inputs_job = AMSMSDJob(job, settings=sets)
        msd_scm_inputs_job.run()
    else:
        # Cannot load the AMSMSDJob directly, so simulate running the job by loading the kf into the results
        job = AMSJob.load_external(rkf_tools_plot / "md/ams.rkf", settings=s)
        msd_job = AMSMSDJob(job, name="msd", settings=s_msd)
        msd_job.prerun()
        msd_job.path = Path(rkf_tools_plot / "md/msd")
        results = AMSMSDResults(msd_job)
        results.finished.set()
        results.done.set()
        results.collect()
        msd_job.results = results

        msd_scm_inputs_job = AMSMSDJob(job, name="msd", settings=sets)
        msd_scm_inputs_job.prerun()
        msd_scm_inputs_job.path = Path(rkf_tools_plot / "md/msd")
        results = AMSMSDResults(msd_scm_inputs_job)
        results.finished.set()
        results.done.set()
        results.collect()
        msd_scm_inputs_job.results = results

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    plot_msd(msd_job, ax=ax1)
    plot_msd(msd_scm_inputs_job, ax=ax2)


# ----------------------------------------------------------
# Testing plot_work_function
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_work_function"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=20,
)
def test_plot_work_function(run_calculations, rkf_tools_plot):
    plt.close("all")

    mol = Molecule()
    mol.add_atom(Atom(symbol="Al", coords=(0.0, 0.0, 12.2)))
    mol.add_atom(Atom(symbol="Al", coords=(1.4, 1.4, 10.1)))
    mol.add_atom(Atom(symbol="Al", coords=(0.0, 0.0, 8.1)))
    mol.add_atom(Atom(symbol="Al", coords=(1.4, 1.4, 6.1)))
    mol.lattice = [[2.9, 0.0, 0.0], [0.0, 2.9, 0.0], [0.0, 0.0, 18.3]]

    s = Settings()
    s.input.ams.task = "SinglePoint"
    s.input.QuantumEspresso.K_Points._h = "automatic"
    s.input.QuantumEspresso.K_Points._1 = "3 3 1 0 0 0"
    s.input.QuantumEspresso.System.ecutwfc = 15.0
    s.input.QuantumEspresso.System.occupations = "smearing"
    s.input.QuantumEspresso.System.degauss = 0.05
    s.input.QuantumEspresso.Properties.WorkFunction = "Yes"
    s.input.QuantumEspresso.WorkFunction.idir = 3
    s.input.QuantumEspresso.WorkFunction.awin = 3.85

    if run_calculations:
        job = AMSJob(settings=s, molecule=mol, name="Al_surface")
        job.run()
    else:
        job = AMSJob.load_external(rkf_tools_plot / "Al_surface")

    wf_results = job.results.get_work_function_results("eV", "Angstrom")
    ax = plot_work_function(*wf_results)

    ax.set_title("Electrostatic Potential Profile", fontsize=14)
    ax.set_xlabel("Length (Angstroms)", fontsize=13)
    ax.set_ylabel("Energy (eV)", fontsize=13)


def force_remove_text(ax):
    # Workaround: Matplotlib's @image_comparison misses some text in subplot grids.
    # We manually hide the remaining text here to ensure visual regression tests
    # pass consistently across different operating systems.
    for axis in ax.flat:
        axis.set_xlabel("")
        axis.set_ylabel("")
        for txt in axis.texts:
            txt.set_visible(False)


# ----------------------------------------------------------
# Testing plot_energy_landscape layouts for molecules
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_energy_landscape_molecules_layouts"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=20,
)
def test_plot_energy_landscape_molecules_layouts(run_calculations, rkf_tools_plot, xyz_folder):
    plt.close("all")

    mol = Molecule()
    mol.add_atom(Atom(symbol="H", coords=(0.26799604, 1.56164318, 0.80172174)))
    mol.add_atom(Atom(symbol="O", coords=(0.70317302, 1.15034441, 0.02980438)))
    mol.add_atom(Atom(symbol="N", coords=(0.07385403, -1.26509181, -0.02723174)))
    mol.add_atom(Atom(symbol="C", coords=(0.33895691, -0.13354579, 0.04083563)))

    sett = Settings()
    sett.input.ams.UseSymmetry = "No"
    sett.input.ams.Task = "PESExploration"
    sett.input.ams.PESExploration.RandomSeed = 1
    sett.input.ams.PESExploration.Job = "ProcessSearch"
    sett.input.ams.PESExploration.NumExpeditions = 500
    sett.input.ams.PESExploration.NumExplorers = 4
    sett.input.ams.PESExploration.SaddleSearch.MaxEnergy = 6.0
    sett.input.ams.PESExploration.SaddleSearch.MinEnergyBarrier = 0.1
    sett.input.ams.PESExploration.StructureComparison.UseCovalent = "Yes"
    sett.input.ams.UseSymmetry = "No"
    sett.input.MOPAC.Model = "AM1"

    if run_calculations:
        job = AMSJob(name="HCNO", molecule=mol, settings=sett)
        job.run()
    else:
        job = AMSJob.load_external(rkf_tools_plot / "HCNO")

    energy_landscape = job.results.get_energy_landscape()

    _, ax = plt.subplots(3, 2, figsize=(20, 10))
    plot_energy_landscape(
        energy_landscape, ax=ax[0, 0], layout="auto", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[0, 1], layout="dfs", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[1, 0], layout="bfs", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[1, 1], layout="longest_path", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[2, 0], layout="force", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[2, 1], layout="crossings", show_molecules=True, molecule_plot_backend="plot_molecule"
    )

    force_remove_text(ax)


# ----------------------------------------------------------
# Testing plot_energy_landscape layouts for surfaces
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_energy_landscape_surfaces_layouts"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=20,
)
def test_plot_energy_landscape_surfaces_layouts(run_calculations, rkf_tools_plot, xyz_folder):
    plt.close("all")

    mol = Molecule(xyz_folder / "MetOH_Cu111.xyz")

    sett = Settings()
    sett.runscript.preamble_lines = ["export OMP_NUM_THREADS=1"]
    sett.input.ams.Task = "PESExploration"
    sett.input.ams.PESExploration.RandomSeed = 10
    sett.input.ams.PESExploration.Job = "ProcessSearch"
    sett.input.ams.PESExploration.NumExpeditions = 50
    sett.input.ams.PESExploration.NumExplorers = 4
    sett.input.ams.PESExploration.SaddleSearch.MaxEnergy = 3.0
    sett.input.ams.PESExploration.SaddleSearch.MinEnergyBarrier = 0.1
    sett.input.ams.PESExploration.SaddleSearch.DisplaceAlongNormalModesWeight = 0.7
    sett.input.ams.PESExploration.StructureComparison.DistanceDifference = 0.5
    sett.input.ams.PESExploration.StructureComparison.EnergyDifference = 0.5
    sett.input.ReaxFF.ForceField = "CuCHO.ff"
    sett.input.ReaxFF.Charges.Solver = "Direct"
    sett.input.ams.Constraints.FixedRegion = "surface"

    if run_calculations:
        job = AMSJob(name="MetOH_Cu111", molecule=mol, settings=sett)
        job.run()
    else:
        job = AMSJob.load_external(rkf_tools_plot / "MetOH_Cu111")

    energy_landscape = job.results.get_energy_landscape()

    _, ax = plt.subplots(3, 2, figsize=(40, 20))
    plot_energy_landscape(
        energy_landscape, ax=ax[0, 0], layout="auto", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[0, 1], layout="dfs", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[1, 0], layout="bfs", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[1, 1], layout="longest_path", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[2, 0], layout="force", show_molecules=True, molecule_plot_backend="plot_molecule"
    )
    plot_energy_landscape(
        energy_landscape, ax=ax[2, 1], layout="crossings", show_molecules=True, molecule_plot_backend="plot_molecule"
    )

    force_remove_text(ax)


# ----------------------------------------------------------
# Testing plot_energy_landscape functions
# ----------------------------------------------------------
@image_comparison(
    baseline_images=["plot_energy_landscape_functions"],
    remove_text=True,
    extensions=["png"],
    style="mpl20",
    tol=20,
)
def test_plot_energy_landscape_functions(run_calculations, rkf_tools_plot, xyz_folder):
    plt.close("all")

    mol = Molecule()
    mol.add_atom(Atom(symbol="H", coords=(0.26799604, 1.56164318, 0.80172174)))
    mol.add_atom(Atom(symbol="O", coords=(0.70317302, 1.15034441, 0.02980438)))
    mol.add_atom(Atom(symbol="N", coords=(0.07385403, -1.26509181, -0.02723174)))
    mol.add_atom(Atom(symbol="C", coords=(0.33895691, -0.13354579, 0.04083563)))

    sett = Settings()
    sett.input.ams.UseSymmetry = "No"
    sett.input.ams.Task = "PESExploration"
    sett.input.ams.PESExploration.RandomSeed = 1
    sett.input.ams.PESExploration.Job = "ProcessSearch"
    sett.input.ams.PESExploration.NumExpeditions = 500
    sett.input.ams.PESExploration.NumExplorers = 4
    sett.input.ams.PESExploration.SaddleSearch.MaxEnergy = 6.0
    sett.input.ams.PESExploration.SaddleSearch.MinEnergyBarrier = 0.1
    sett.input.ams.PESExploration.StructureComparison.UseCovalent = "Yes"
    sett.input.ams.UseSymmetry = "No"
    sett.input.MOPAC.Model = "AM1"

    if run_calculations:
        job = AMSJob(name="HCNO", molecule=mol, settings=sett)
        job.run()
    else:
        job = AMSJob.load_external(rkf_tools_plot / "HCNO")

    energy_landscape = job.results.get_energy_landscape()
    energy_landscape_select_states = energy_landscape.select_states([1, 6, 5, 7, 3], keep_original_ids=True)
    energy_landscape_accessible_states = energy_landscape.accessible_states(3, 3.5, unit="eV", keep_original_ids=True)

    _, ax = plt.subplots(2, 1, figsize=(10, 5))
    plot_energy_landscape(
        energy_landscape_select_states,
        ax=ax[0],
        layout="auto",
        show_molecules=True,
        molecule_y_offset=0.2,
        molecule_scale=0.4,
        molecule_plot_backend="plot_molecule",
        molecule_plot_kwargs={"rotation": "0x,0y,0z"},
        molecule_plot_kwargs_by_state={
            7: {"rotation": "90x,0y,0z"},
            3: {"rotation": "0x,90y,0z"},
        },
    )

    plot_energy_landscape(
        energy_landscape_accessible_states,
        ax=ax[1],
        layout="auto",
        show_molecules=True,
        molecule_y_offset=0.2,
        molecule_scale=0.3,
        molecule_plot_backend="plot_molecule",
        molecule_plot_kwargs={"rotation": "0x,0y,0z"},
        molecule_plot_kwargs_by_state={
            7: {"rotation": "90x,0y,0z"},
            3: {"rotation": "0x,90y,0z"},
        },
    )

    force_remove_text(ax)

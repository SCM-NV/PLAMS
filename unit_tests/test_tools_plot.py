#!/usr/bin/env amspython
# coding: utf-8

import os
import matplotlib.pyplot as plt
import pytest
from matplotlib.testing.decorators import image_comparison
from pathlib import Path
import numpy as np

from scm.plams.core.functions import Settings
from scm.plams.mol.molecule import Molecule
from scm.plams.mol.atom import Atom
from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.recipes.md.trajectoryanalysis import AMSMSDJob, AMSMSDResults
from scm.plams.tools.plot import (
    plot_molecule,
    plot_correlation,
    plot_band_structure,
    plot_work_function,
    get_correlation_xy,
    plot_msd,
)


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
@image_comparison(baseline_images=["plot_molecule"], remove_text=True, extensions=["png"], style="mpl20", tol=0.1)
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

    plt.pause(2)


# ----------------------------------------------------------
# Testing plot_band_structure
# ----------------------------------------------------------
@image_comparison(baseline_images=["plot_band_structure"], remove_text=True, extensions=["png"], style="mpl20", tol=3)
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

    plt.pause(2)


# ----------------------------------------------------------
# Testing plot_correlation & get_correlation_xy
# ----------------------------------------------------------
@image_comparison(baseline_images=["plot_correlation"], remove_text=True, extensions=["png"], style="mpl20", tol=1.0)
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

    plt.pause(2)

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
# Testing plot_msd
# ----------------------------------------------------------
@image_comparison(baseline_images=["plot_msd"], remove_text=True, extensions=["png"], style="mpl20", tol=4)
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

    plt.pause(2)


# ----------------------------------------------------------
# Testing plot_work_function
# ----------------------------------------------------------
@image_comparison(baseline_images=["plot_work_function"], remove_text=True, extensions=["png"], style="mpl20", tol=11)
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

    plt.pause(2)

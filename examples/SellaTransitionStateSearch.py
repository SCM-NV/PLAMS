#!/usr/bin/env amspython
import os
import time

from ase.io import read
from scm.plams import *
from scm.plams.interfaces.adfsuite.ase_calculator import AMSCalculator
from sella import Sella

"""

This examples runs a transition state search with the Sella python package and
compares results and timings to the builtin transition state search in AMS (using a DFTB hessian to initialize the DFT transition state search).

This example requires an ADF license.

"""


def adf_settings():
    s = Settings()
    s.input.adf.basis.type = "DZP"
    s.input.adf.basis.core = "None"
    s.input.adf.xc.gga = "pbe"
    s.input.adf.symmetry = "nosym"
    return s


def run_sella(molecule):
    """Run a TS search with Sella"""
    s = adf_settings()
    s.input.ams.task = "SinglePoint"
    s.input.ams.properties.gradients = "Yes"

    atoms = toASE(molecule)

    if os.path.exists("sella.log"):
        os.remove("sella.log")
    if os.path.exists("sella.traj"):
        os.remove("sella.traj")

    with AMSCalculator(settings=s, amsworker=True) as calc:
        atoms.calc = calc
        opt = Sella(atoms, internal=True, logfile="sella.log", trajectory="sella.traj")
        opt.run(fmax=0.027, steps=50)  # 0.027 eV/ang roughly matches AMS default criterion 0.001 Ha/ang

    traj = read("sella.traj", ":")
    print(f"Sella converged after {len(traj)} single-point calculations")
    print(f"Final energy: {atoms.get_potential_energy()} eV")

    optimized_mol = fromASE(atoms)
    return optimized_mol


def run_ams(molecule):
    """Run DFT transition state search but use initial hessian calculated with DFTB"""
    # this line is not required in AMS2025+
    init()

    s = adf_settings()
    s.input.ams.task = "TransitionStateSearch"
    s.input.ams.geometryoptimization.initialhessian.type = "CalculateWithFastEngine"
    job = AMSJob(settings=s, molecule=molecule, name="ams_ts")
    job.run()
    print(f"AMS finished after {len(job.results.get_history_property('Energy'))} iterations")
    print(f"AMS optimized energy: {job.results.get_energy(unit='eV')} eV")


def get_molecule():
    return AMSJob.from_input(
        """
        System
          atoms
             C         0.049484    0.042994    0.000000
             H        -0.068980    0.638928   -0.915972
             H        -0.068980    0.638928    0.915972
             H        -0.841513   -0.626342    0.000000
             H         0.555494   -1.148227    0.000000
             Hg        2.303289   -0.007233    0.000000
             Cl        4.429752    0.776056    0.000000
             Cl        1.342057   -2.676083    0.000000
          end
        End

    """
    ).molecule[""]


def main():
    # this line is not required in AMS2025+
    init()

    os.environ["OMP_NUM_THREADS"] = "1"
    mol = get_molecule()

    start = time.time()
    run_sella(mol)
    print(f"Sella finished in {time.time()-start} seconds")

    start = time.time()
    run_ams(mol)
    print(f"AMS finished in {time.time()-start} seconds")


if __name__ == "__main__":
    main()

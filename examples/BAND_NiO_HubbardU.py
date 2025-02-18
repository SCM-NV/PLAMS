#!/usr/bin/env amspython
from scm.plams import *

""" Run as: $AMSBIN/amspython BAND_NiO_HubbardU.py """


def main():
    # this line is not required in AMS2025+
    init()

    d = 2.085
    mol = Molecule()
    mol.add_atom(Atom(symbol="Ni", coords=(0, 0, 0)))
    mol.add_atom(Atom(symbol="O", coords=(d, d, d)))
    mol.lattice = [[0.0, d, d], [d, 0.0, d], [d, d, 0.0]]

    s = Settings()
    s.input.ams.task = "SinglePoint"
    s.input.band.Unrestricted = "yes"
    s.input.band.XC.GGA = "BP86"
    s.input.band.Basis.Type = "DZP"
    s.input.band.KSpace.Quality = "Basic"
    s.input.band.NumericalQuality = "Normal"
    s.input.band.DOS.CalcPDOS = "Yes"
    s.input.band.HubbardU.Enabled = "Yes"
    s.input.band.HubbardU.UValue = "0.6 0.0"
    s.input.band.HubbardU.LValue = "2 -1"

    job = AMSJob(settings=s, molecule=mol, name="NiO")
    job.run()

    toeV = Units.convert(1.0, "hartree", "eV")
    topvb = job.results.readrkf("BandStructure", "TopValenceBand", file="engine") * toeV
    bottomcb = job.results.readrkf("BandStructure", "BottomConductionBand", file="engine") * toeV
    gap = bottomcb - topvb

    log("Results:")
    log(f"Top of valence band:       {topvb:.2f} eV")
    log(f"Bottom of conduction band: {bottomcb:.2f} eV")
    log(f"Band gap:                  {gap:.2f} eV")


if __name__ == "__main__":
    main()

#!/usr/bin/env amspython
import scm.plams as plams

""" Run as: $AMSBIN/amspython BAND_NiO_HubbardU.py """


def main():

    # this line is not required in AMS2025+
    plams.init()

    d = 2.085
    mol = plams.Molecule()
    mol.add_atom(plams.Atom(symbol="Ni", coords=(0, 0, 0)))
    mol.add_atom(plams.Atom(symbol="O", coords=(d, d, d)))
    mol.lattice = [[0.0, d, d], [d, 0.0, d], [d, d, 0.0]]

    ams_settings = plams.Settings()
    ams_settings.input.ams.task = "SinglePoint"
    ams_settings.input.band.Unrestricted = "yes"
    ams_settings.input.band.XC.GGA = "BP86"
    ams_settings.input.band.Basis.Type = "TZP"
    ams_settings.input.band.KSpace.Type = "Symmetric"
    ams_settings.input.band.KSpace.Symmetric = "KInteg=3"
    ams_settings.input.band.NumericalQuality = "Normal"
    ams_settings.input.band.HubbardU.Enabled = "Yes"
    ams_settings.input.band.HubbardU.PrintOccupations = "Yes"
    ams_settings.input.band.HubbardU.Atom = [plams.Settings()]
    ams_settings.input.band.HubbardU.Atom[0].Element = "Ni"
    ams_settings.input.band.HubbardU.Atom[0].UValue = 0.3
    ams_settings.input.band.HubbardU.Atom[0].LValue = "d"
    ams_settings.input.band.Print = "AtomicChargesDetails"
    ams_settings.input.band.DOS.CalcPDOS = "Yes"

    job = plams.AMSJob(settings=ams_settings, molecule=mol, name="NiO")
    job.run()

    toeV = plams.Units.convert(1.0, "hartree", "eV")
    topvb = job.results.readrkf("BandStructure", "TopValenceBand", file="engine") * toeV
    bottomcb = job.results.readrkf("BandStructure", "BottomConductionBand", file="engine") * toeV
    gap = bottomcb - topvb

    plams.log("Results:")
    plams.log(f"Top of valence band:       {topvb:7.2f} eV")
    plams.log(f"Bottom of conduction band: {bottomcb:7.2f} eV")
    plams.log(f"Band gap:                  {gap:7.2f} eV")


if __name__ == "__main__":
    main()

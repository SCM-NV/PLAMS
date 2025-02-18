from scm.plams import init, Molecule, Settings, AMSJob
from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob

init(folder="PLAMS_CCSD_N2_DZ")

mol = Molecule("N2.xyz")
adfsett = Settings()
adfsett.input.ams.task = "SinglePoint"
adfsett.input.adf.IntegralsToFile = "SERENITY"
adfsett.input.adf.TotalEnergy = ""
adfsett.input.adf.NumericalQuality = "VeryGood"
adfsett.input.adf.basis.core = "None"
adfsett.input.adf.basis.type = "DZ"
adfsett.input.adf.basis.CreateOutput = "yes"
adfsett.input.adf.relativity = "Level=None"
adfsett.input.adf.symmetry = "NoSym"
adfsett.input.adf.xc.hartreefock = ""

adfjob = AMSJob(molecule=mol, settings=adfsett, name="ADF_N2_HF")
adfjob.run()
bonding_energy = adfjob.results.get_energy()
print(f"ADF bonding energy: {bonding_energy} hartree")


sersett = SerenitySettings()
sersett.input.system.N2.geometry = "N2.xyz"

sersett.input.task.CC.system = "N2"
sersett.input.task.CC.level = "CCSD"

serjob = SerenityJob(settings=sersett, name="Serenity_N2_CCSD")
serjob.run()
ccsd_correction = serjob.results.get_ccsd_energy_correction()
print(f"Serenity CCSD energy correction: {ccsd_correction} hartree")

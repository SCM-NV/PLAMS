from scm.plams import init, Molecule, Settings, AMSJob
from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob

init(folder="PLAMS_CCSD_T_HF_TZP")

mol = Molecule("HF.xyz")
adfsett = Settings()
adfsett.input.ams.task = "SinglePoint"
adfsett.input.adf.IntegralsToFile = "SERENITY"
adfsett.input.adf.TotalEnergy = ""
adfsett.input.adf.NumericalQuality = "VeryGood"
adfsett.input.adf.basis.core = "None"
adfsett.input.adf.basis.type = "TZP"
adfsett.input.adf.basis.CreateOutput = "yes"
adfsett.input.adf.relativity = "Level=None"
adfsett.input.adf.symmetry = "NoSym"
adfsett.input.adf.xc.hartreefock = ""

adfjob = AMSJob(molecule=mol, settings=adfsett, name="ADF_HF_HF")
adfjob.run()
bonding_energy = adfjob.results.get_energy()
print(f"ADF bonding energy: {bonding_energy} hartree")


sersett = SerenitySettings()
sersett.input.system.HF.geometry = "HF.xyz"

sersett.input.task.CC.system = "HF"
sersett.input.task.CC.level = "CCSD(T)"

serjob = SerenityJob(settings=sersett, name="Serenity_HF_CCSD_T")
serjob.run()
ccsd_correction = serjob.results.get_ccsd_energy_correction()
print(f"Serenity CCSD energy correction: {ccsd_correction} hartree")

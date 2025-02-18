from scm.plams import init
from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob

init(folder="test_localization")

sersett = SerenitySettings()
# The path to the geometry file
sersett.input.system.Riboflavin.geometry = "riboflavin.xyz"
sersett.input.system.Riboflavin.method = "DFT"
sersett.input.system.Riboflavin.dft.functional = "PBE"
sersett.input.system.Riboflavin.basis.label = "Def2-SVP"
# Performs a SCF calculation for the given active (act) system
sersett.input.task.SCF.act = "Riboflavin"
# Localize the orbtials
sersett.input.task.loc.act = "Riboflavin"
sersett.input.task.loc.locType = "IBO"
# Print the orbitals to file
sersett.input.task.cube.act = "Riboflavin"
sersett.input.task.cube.occOrbitals = "true"

serjob = SerenityJob(settings=sersett, name="Serenity_localization")
serjob.run()

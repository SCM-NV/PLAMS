from scm.plams import init
from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob

init(folder="test_mp2")

sersett = SerenitySettings()
# The path to the geometry file
sersett.input.system.gly_gly_gly.geometry = "gly-gly-gly.xyz"
sersett.input.system.gly_gly_gly.method = "HF"
sersett.input.system.gly_gly_gly.charge = "0"
sersett.input.system.gly_gly_gly.spin = "0"
sersett.input.system.gly_gly_gly.basis.label = "Def2-SVP"
# Performs a SCF calculation for the given active (act) system
sersett.input.task.SCF.act = "gly_gly_gly"
# Localize the orbtials
sersett.input.task.MP2.act = "gly_gly_gly"

serjob = SerenityJob(settings=sersett, name="Serenity_mp2")
serjob.run()

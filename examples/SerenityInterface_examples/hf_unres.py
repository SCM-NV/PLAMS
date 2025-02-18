from scm.plams import init
from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob

init(folder="test_methylradical")

sersett = SerenitySettings()
# The path to the geometry file
sersett.input.system.MethylRadical.geometry = "methylRadical.xyz"
# You can define the excess number of alpha electrons here
sersett.input.system.MethylRadical.spin = "1"
# This means that we allow singly occupied orbitals which leads to different
# spatial-orbitals for alpha and beta electrons.
# Serenity will choose this automatically as soon as spin != 0.
sersett.input.system.MethylRadical.scfMode = "UNRESTRICTED"
sersett.input.system.MethylRadical.basis.label = "Def2-SVP"
# Performs a SCF calculation for the given active (act) system
sersett.input.task.SCF.act = "MethylRadical"

serjob = SerenityJob(settings=sersett, name="Serenity_methylradical")
serjob.run()

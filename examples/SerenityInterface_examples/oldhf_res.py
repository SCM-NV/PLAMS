from scm.plams import init
from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob

init(folder="test_gly-gly-gly")

# instead of giving the geometry key you can also do:
# mol = Molecule('gly-gly-gly.xyz')
# serjob = SerenityJob(molecule=mol, settings=sersett, name='Serenity_gly-gly-gly')
# the name of the molecule does not have to match the name of the system it will be given as mol.xyz not depending on the given name


sersett = SerenitySettings()
# The path to the geometry file
sersett.input.system.gly_gly_gly.geometry = "gly-gly-gly.xyz"
# You can define the charge and the excess number of alpha electrons here
sersett.input.system.gly_gly_gly.charge = "0"
sersett.input.system.gly_gly_gly.spin = "0"
# The electronic structure method
sersett.input.system.gly_gly_gly.method = "HF"
sersett.input.system.gly_gly_gly.basis.label = "6-31GS"
# Performs a SCF calculation for the given active (act) system
sersett.input.task.SCF.act = "gly_gly_gly"

serjob = SerenityJob(settings=sersett, name="Serenity_gly-gly-gly")
serjob.run()

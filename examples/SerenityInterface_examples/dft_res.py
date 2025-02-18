from scm.plams import init
from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob

init(folder="test_riboflavin")

# instead of giving the geometry key and the file as value you can also do:
# Riboflavin = Molecule()
# Riboflavin.add_atom(Atom(symbol='C', coords=(-2.335,1.498,5.390)))
# Riboflavin.add_atom(Atom(symbol='C', coords=(-2.573,0.058,4.841)))
# Riboflavin.add_atom(Atom(symbol='C', coords=(-3.154,0.011,3.393)))
# ...
# name of the molecule does not necessarily have to match the name of the system. You have to give the molecule to the job
# serjob = SerenityJob(molecule=Riboflavin, settings=sersett, name='Serenity_Riboflavin')


sersett = SerenitySettings()
# The path to the geometry file
sersett.input.system.Riboflavin.geometry = "riboflavin.xyz"
# The electronic structure method
sersett.input.system.Riboflavin.method = "DFT"
sersett.input.system.Riboflavin.basis.label = "6-31GS"
# Settings used for the DFT calculation.
sersett.input.system.Riboflavin.dft.functional = "PBE"
# Performs a SCF calculation for the given active (act) system
sersett.input.task.SCF.act = "Riboflavin"
# Every task in Serenity has a print-level which can be tuned to
# adjust the amount of output generated.
sersett.input.task.SCF.printLevel = "normal"

serjob = SerenityJob(settings=sersett, name="Serenity_Riboflavin")
serjob.run()

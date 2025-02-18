from scm.plams import init, Molecule, Atom
from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob

init(folder="test")

mol1 = Molecule()
mol1.add_atom(Atom(symbol="O", coords=(0, 0, 0)))
mol1.add_atom(Atom(symbol="H", coords=(0.758602, 0, 0.504284)))
mol1.add_atom(Atom(symbol="H", coords=(0.260455, 0, -0.872893)))

mol2 = Molecule()
mol2.add_atom(Atom(symbol="O", coords=(3, 0.5, 0)))
mol2.add_atom(Atom(symbol="H", coords=(3.758602, 0.5, 0.504284)))
mol2.add_atom(Atom(symbol="H", coords=(3.260455, 0.5, -0.872893)))

sersett = SerenitySettings()
sersett.input.system.water1.method = "dft"
sersett.input.system.water1.dft.functional = "pw91"

sersett.input.system.water2.method = "dft"
sersett.input.system.water2.dft.functional = "pw91"

sersett.input.task.fde.act = "water1"
sersett.input.task.fde.env = "water2"
sersett.input.task.fde.emb.naddxcfunc = "pw91"
sersett.input.task.fde.emb.naddkinfunc = "pw91k"
# when using a dictionary and the molecule class to provide geometries you have no make sure that the given names match the systems that should use the respective geometry
serjob = SerenityJob(molecule={"water1": mol1, "water2": mol2}, settings=sersett, name="water_dimer")
serjob.run()

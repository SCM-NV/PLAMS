from scm.plams import Settings
from scm.plams import fromASE
from scm.plams.recipes.bandfragment import BANDFragmentJob

# build the surface
from ase.build import fcc111, add_adsorbate

mol = fcc111("Au", size=(2, 2, 3))
add_adsorbate(mol, "H", 1.5, "ontop")
mol.center(vacuum=10.0, axis=2)

# separate the fragments
surface = mol.copy()
symbols = surface.get_chemical_symbols()
del surface[[i for i in range(len(symbols)) if symbols[i] != "Au"]]
adsorbate = mol.copy()
del adsorbate[[i for i in range(len(symbols)) if symbols[i] == "Au"]]

# optional: load the optimized molecules
from ase import io

surface_opt = io.read("surface_opt.xyz")
adsorbate_opt = io.read("adsorbate_opt.xyz")
assert len(surface_opt) == len(surface)
assert len(adsorbate_opt) == len(adsorbate)

# settings for job
base_settings = Settings()
base_settings.input.ams.task = "SinglePoint"
base_settings.input.band.basis.type = "TZP"
base_settings.input.band.basis.core = "Medium"
base_settings.input.band.dos.calcdos = "No"
base_settings.input.band.kspace.regular.numberofpoints = "5 5 1"
base_settings.input.band.beckegrid.quality = "Good"
base_settings.input.band.zlmfit.quality = "Good"
base_settings.input.band.usesymmetry = False
base_settings.input.band.xc.gga = "PBE"
base_settings.input.band.xc.dispersion = "Grimme4"

eda_settings = Settings()
eda_settings.input.band.peda = ""

eda_job = BANDFragmentJob(
    fragment1=fromASE(surface),
    fragment2=fromASE(adsorbate),
    settings=base_settings,
    full_settings=eda_settings,
    fragment1_opt=fromASE(surface_opt),
    fragment2_opt=fromASE(adsorbate_opt),
)

results = eda_job.run()
eda_res = eda_job.results.get_energy_decomposition()
print("{:<20} {:>10}".format("Term", "Energy [kJ/mol]"))
for key, value in eda_res.items():
    print("{:<20} {:>10.4f}".format(key, value))

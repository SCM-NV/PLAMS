from scm.plams import Settings, Molecule, init
from scm.plams.recipes.adffragment import ADFFragmentJob

# this line is not required in AMS2025+
init()

common = Settings()  # common settings for all 3 jobs
common.input.ams.Task = "SinglePoint"
common.input.adf.basis.type = "DZP"
common.input.adf.xc.gga = "PBE"
common.input.adf.symmetry = "NOSYM"

full = Settings()  # additional settings for full system calculation
full.input.adf.etsnocv  # empty block
full.input.adf.print = "etslowdin"

mol1 = Molecule("ethene.xyz")
mol2 = Molecule("butadiene.xyz")

j = ADFFragmentJob(fragment1=mol1, fragment2=mol2, settings=common, full_settings=full)
r = j.run()

print("Energy decomposition:")
decom = r.get_energy_decomposition()
for k in decom:
    print(k, decom[k])
print(j.full.results.readrkf("NOCV", "NOCV_eigenvalues_restricted", "engine"))

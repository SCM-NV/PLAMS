#!/usr/bin/env amspython
# coding: utf-8

# ## Initialization

#!/usr/bin/env amspython
from scm.plams import Molecule, Settings, ReorganizationEnergyJob, init, AMSJob

# this line is not required in AMS2025+
init()


# ## Define molecule
#
# Normally you would read it from an xyz file, but here is for convenience explicit code

# molecule = Molecule("pyrrole.xyz")


def get_molecule(input_string):
    job = AMSJob.from_input(input_string)
    return job.molecule[""]


molecule = get_molecule("""
System
    Atoms
        C      -1.12843000       0.00000000      -0.35463200
        C      -0.71293000       0.00000000       0.96463800
        C       0.71293000       0.00000000       0.96463800
        C       1.12843000       0.00000000      -0.35463200
        N       0.00000000       0.00000000      -1.14563200
        H       0.00000000       0.00000000      -2.15713200
        H      -2.12074000       0.00000000      -0.79100200
        H      -1.36515000       0.00000000       1.83237800
        H       1.36515000       0.00000000       1.83237800
        H       2.12074000       0.00000000      -0.79100200
    End
End
""")

molecule.properties.name = "pyrrole"  # normally the name of the xyz file


# ## Setup and run job

# Generic settings of the calculation
# (for quantitatively better results, use better settings)
common_settings = Settings()
common_settings.input.adf.Basis.Type = "DZ"

# Specific settings for the neutral calculation.
# Nothing special needs to be done for the neutral calculation,
# so we just use an empty settings.
neutral_settings = Settings()

# Specific settings for the anion calculation:
anion_settings = Settings()
anion_settings.input.ams.System.Charge = -1
anion_settings.input.adf.Unrestricted = "Yes"
anion_settings.input.adf.SpinPolarization = 1

# Create and run the ReorganizationEnergyJob:
job = ReorganizationEnergyJob(
    molecule, common_settings, neutral_settings, anion_settings, name=molecule.properties.name
)
job.run()


# ## Fetch and print the results:

energy_unit = "eV"
energies = job.results.get_all_energies(energy_unit)
reorganization_energy = job.results.reorganization_energy(energy_unit)

print("")
print("== Results ==")
print("")
print(f"Molecule: {molecule.properties.name}")
print("State A: neutral")
print("State B: anion")
print("")
print(f"Reorganization energy: {reorganization_energy:.6f} [{energy_unit}]")
print("")
print(f"|   State   | Optim Geo | Energy [{energy_unit}]")
print(f'|     A     |     A     | {energies["state A geo A"]:.6f}')
print(f'|     A     |     B     | {energies["state A geo B"]:.6f}')
print(f'|     B     |     A     | {energies["state B geo A"]:.6f}')
print(f'|     B     |     B     | {energies["state B geo B"]:.6f}')
print("")

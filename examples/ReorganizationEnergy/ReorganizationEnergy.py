#!/usr/bin/env amspython
from scm.plams import Molecule, Settings, ReorganizationEnergyJob, init

# this line is not required in AMS2025+
init()

# Compute the neutral-anion reorganization energy of pyrrole
# using ADF as computational engine

molecule = Molecule("pyrrole.xyz")

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

# Fetch and print the results:

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

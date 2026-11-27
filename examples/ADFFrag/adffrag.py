#!/usr/bin/env amspython
# coding: utf-8

# ## Initialization

from scm.plams import Settings, Molecule, init, AMSJob, Units
from scm.plams.recipes.adffragment import ADFFragmentJob

# this line is not required in AMS2025+
init()


# ## Define the molecules
# For convenience we define here two molecules, normally you would read them from xyz files


def get_molecule(input_string):
    job = AMSJob.from_input(input_string)
    return job.molecule[""]


mol1 = get_molecule("""
System
    Atoms
        C      -0.75086900       1.37782400      -2.43303700
        C      -0.05392100       2.51281000      -2.41769100
        H      -1.78964800       1.33942600      -2.09651100
        H      -0.30849400       0.43896500      -2.76734700
        H      -0.49177100       3.45043100      -2.06789100
        H       0.98633900       2.54913500      -2.74329400
    End
End
""")


mol2 = get_molecule("""
System
    Atoms
        C       0.14667300      -0.21503500       0.40053800
        C       1.45297400      -0.07836900       0.12424400
        C       2.23119700       1.15868100       0.12912100
        C       1.78331500       2.39701500       0.38779700
        H      -0.48348000       0.63110600       0.67664100
        H      -0.33261900      -1.19332100       0.35411600
        H       2.01546300      -0.97840100      -0.14506700
        H       3.29046200       1.03872500      -0.12139700
        H       2.45728900       3.25301000       0.35150400
        H       0.74193400       2.60120700       0.64028800
    End
End
""")


# ## Setup and run the job

common = Settings()  # common settings for all 3 jobs
common.input.ams.Task = "SinglePoint"
common.input.adf.basis.type = "DZP"
common.input.adf.xc.gga = "PBE"
common.input.adf.symmetry = "NOSYM"

full = Settings()  # additional settings for full system calculation
full.input.adf.etsnocv  # empty block
full.input.adf.print = "etslowdin"

# normally you would read here the two molecules from xyz files.
# mol1 = Molecule("ethene.xyz")
# mol2 = Molecule("butadiene.xyz")

j = ADFFragmentJob(fragment1=mol1, fragment2=mol2, settings=common, full_settings=full)
r = j.run()


# ## Print the results


def print_eterm(energy_term, energy):
    print(
        f'{energy_term:>30s} {energy:16.4f} {Units.convert(energy, "au", "eV"):16.3f} {Units.convert(energy, "au", "kcal/mol"):16.2f} {Units.convert(energy, "au", "kJ/mol"):16.2f}'
    )


def print_bonding_energy_terms(r):
    print("Energy terms contributing to the bond energy (with respect to the fragments):")

    bond_energy = r.get_energy()
    decom = r.get_energy_decomposition()
    print(f'\n{"term":>30s} {"Hartree":>16s} {"eV":>16s} {"kcal/mol":>16s} {"kJ/mol":>16s}')
    for energy_term, energy in decom.items():
        print_eterm(energy_term, energy)

    print_eterm("total bond energy", bond_energy)
    print("")


def print_eda_terms(job):
    bond_energy = job.full.results.readrkf("Energy", "Bond Energy", "adf")
    steric_interaction = job.full.results.readrkf("Energy", "Steric Total", "adf")
    orbital_interaction = job.full.results.readrkf("Energy", "Orb.Int. Total", "adf")
    print("\nFragment based energy decomposition analysis of the bond energy:")
    print(f'\n{"term":>30s} {"Hartree":>16s} {"eV":>16s} {"kcal/mol":>16s} {"kJ/mol":>16s}')
    print_eterm("Steric interaction", steric_interaction)
    print_eterm("Orbital interaction", orbital_interaction)
    print_eterm("total bond energy", bond_energy)
    print("")


def print_nocv_decomposition():
    print("NOCV decomposition of the orbital interaction term\n")

    print("The NOCV eigenvalues are occupation numbers, they should come in pairs,")
    print("with one negative value mirrored by a positive value.")
    print("The orbital interaction energy contribution is calculated for each NOCV pair.")
    print("")

    nocv_eigenvalues = j.full.results.readrkf("NOCV", "NOCV_eigenvalues_restricted", "engine")
    nocv_orbitalinteraction = j.full.results.readrkf("NOCV", "NOCV_oi_restricted", "engine")

    n_pairs = int(len(nocv_eigenvalues) / 2)
    threshold = 0.001

    print(f'{"index":>9s} {"neg":>9s} {"pos":>9s} {"kcal/mol":>10s}')
    for index in range(n_pairs):
        pop1 = nocv_eigenvalues[index]
        pop2 = nocv_eigenvalues[len(nocv_eigenvalues) - index - 1]

        if (abs(pop1) + abs(pop2)) < threshold:
            continue

        orbitalinteraction = (
            nocv_orbitalinteraction[index] + nocv_orbitalinteraction[len(nocv_orbitalinteraction) - index - 1]
        )
        print(f"{index:9d} {pop1:9.3f} {pop2:9.3f} {orbitalinteraction:10.2f}")


print_bonding_energy_terms(r)

print_eda_terms(j)

print_nocv_decomposition()

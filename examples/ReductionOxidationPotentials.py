#!/usr/bin/env amspython
from scm.plams import from_smiles, CRSJob, Settings, init
from scm.plams.recipes.redox import AMSRedoxScreeningJob, AMSRedoxDirectJob, AMSRedoxThermodynamicCycleJob


def main():
    # this line is not required in AMS2025+
    init()

    mol = from_smiles("C1=CC(=O)C=CC1=O", forcefield="uff")  # benzoquinone
    solvent_name = "Water"  # Solvent for AMSRedoxDirectJob and AMSRedoxThermodynamicCycleJob. See the ADF Solvation documentation for which solvents are available.
    solvent_coskf = CRSJob.database() + "/Water.coskf"  # .coskf file or AMSRedoxScreeningJob

    vibrations = True  # set to False to turn off vibrations for AMSRedoxDirectJob and AMSRedoxThermodynamicCycleJob
    reduction = True
    oxidation = True  # not really relevant for benzoquinone, set to False to not calculate

    # DFT settings for AMSRedoxDirectJob and AMSRedoxThermodynamicCycleJob
    # Note: These are for demonstration purposes only. Use a better functional/settings for production purposes.
    s = Settings()
    s.input.adf.basis.type = "DZP"
    s.input.adf.xc.gga = "PBE"
    s.input.ams.GeometryOptimization.Convergence.Quality = "Basic"

    jobs = [
        AMSRedoxScreeningJob(
            name="quick_screening", molecule=mol, reduction=reduction, oxidation=oxidation, solvent_coskf=solvent_coskf
        ),
        AMSRedoxDirectJob(
            name="direct_best_method",
            molecule=mol,
            settings=s,
            reduction=reduction,
            oxidation=oxidation,
            vibrations=vibrations,
            solvent=solvent_name,
        ),
        AMSRedoxThermodynamicCycleJob(
            name="thermodynamic_cycle",
            molecule=mol,
            settings=s,
            reduction=reduction,
            oxidation=oxidation,
            vibrations=vibrations,
            solvent=solvent_name,
        ),
    ]

    for job in jobs:
        job.run()
        print_results([job])

    print("Final summary:")
    print_results(jobs)


def print_results(jobs):
    SHE = 4.42  # standard hydrogen electrode in eV on absolute scale

    print("The experimental reduction potential of benzoquinone is +0.10 V vs. SHE")
    print(
        "{:24s} {:24s} {:24s} {:24s} {:24s}".format(
            "Jobname", "Eox(vib,rel-to-SHE)[V]", "Ered(vib,rel-to-SHE)[V]", "Eox(rel-to-SHE)[V]", "Ered(rel-to-SHE)[V]"
        )
    )

    for job in jobs:
        s = f"{job.name:24s}"
        for vibrations in [True, False]:
            try:
                Eox = job.results.get_oxidation_potential(vibrations=vibrations) - SHE
                Eox = f"{Eox:.2f}"
            except:
                Eox = "N/A"
            s += f" {Eox:24s}"

            try:
                Ered = job.results.get_reduction_potential(vibrations=vibrations) - SHE
                Ered = f"{Ered:.2f}"
            except:
                Ered = "N/A"
            s += f" {Ered:24s}"

        print(s)


if __name__ == "__main__":
    main()

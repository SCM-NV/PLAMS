""" Deprecated, do not use """

import argparse
import os

from scm.plams import (
    AMSJob,
    CRSJob,
    KFFile,
    Molecule,
    Results,
    Settings,
    Units,
    finish,
    init,
)

# import redox_defaults as defaults


__all__ = ["redox_potential"]


### ==== SETTINGS ==== ###


### ==== DFT SETTINGS ==== ###
def DFT_defaults():
    s = Settings()
    s.input.ams.task = "GeometryOptimization"
    s.input.adf.basis.type = "TZ2P"
    s.input.adf.basis.core = "None"
    s.input.adf.xc.hybrid = "B3LYP"
    s.input.adf.xc.Dispersion = "GRIMME3 BJDAMP"
    s.input.adf.Relativity.Level = "None"
    s.input.adf.NumericalQuality = "Good"
    s.input.adf.Symmetry = "NOSYM"
    s.input.ams.UseSymmetry = "No"
    s.input.adf.Unrestricted = "No"
    s.input.adf.SpinPolarization = 0
    s.input.ams.System.Charge = 0
    return s

    ### ==== DFTB SETTINGS ==== ###


def DFTB_defaults():
    s = Settings()
    s.input.ams.task = "GeometryOptimization"
    s.input.DFTB
    s.input.DFTB.Model = "GFN1-xTB"
    s.input.ams.System.Charge = 0
    return s

    ### ==== FREQ SETTINGS ==== ###


def frequencies_defaults():
    s = Settings()
    s.input.ams.properties.NormalModes = "Yes"
    s.input.ams.Properties.PESPointCharacter = "No"
    s.input.ams.NormalModes.ReScanFreqRange = "-1000 0"
    s.input.ams.PESPointCharacter.NegativeEigenvalueTolerance = -0.001
    return s


def get_settings(
    molecule: Molecule = None,
    task: str = "GeometryOptimization",
    state: str = "neutral",
    phase: str = "vacuum",
    frequencies: bool = True,
    use_dftb: bool = False,
    use_COSMORS: bool = False,
    init_charge: int = 0,
) -> Settings:
    """
    Method that generates settings for jobs based on settings provided
    """

    if use_COSMORS:
        sett = DFT_defaults()

        solvation_block = {
            "surf": "Delley",
            "solv": "name=CRS cav0=0.0 cav1=0.0",
            "charged": "method=Conj corr",
            "c-mat": "Exact",
            "scf": "Var All",
            "radii": {
                "H": 1.30,
                "C": 2.00,
                "N": 1.83,
                "O": 1.72,
                "F": 1.72,
                "Si": 2.48,
                "P": 2.13,
                "S": 2.16,
                "Cl": 2.05,
                "Br": 2.16,
                "I": 2.32,
            },
        }

        sett.input.adf.solvation = solvation_block
    else:
        if use_dftb:
            sett = DFTB_defaults()
        else:
            sett = DFT_defaults()

        # load cosmo solvent
        if phase != "vacuum" and not use_COSMORS:
            sett.input.adf.Solvation.Solv = f"name={phase}"

    sett.input.ams.task = task

    if frequencies:
        sett.soft_update(frequencies_defaults())

    # set the charge
    if state == "oxidation":
        charge = init_charge + 1
    elif state == "reduction":
        charge = init_charge - 1
    else:
        charge = init_charge

    if not use_dftb:
        # if we are doing DFT we have to change the spinpolarization and whether to use unrestricted DFT
        num_electrons = sum(atom.atnum for atom in molecule) - charge
        spinpol = num_electrons % 2
        if spinpol == 0:
            unrestricted = "No"
        else:
            unrestricted = "Yes"
        sett.input.adf.SpinPolarization = spinpol
        sett.input.adf.Unrestricted = unrestricted

    sett.input.ams.System.Charge = charge

    return sett


### ==== UTILITY FUNCTIONS ==== ###


def check_termination_succes(result: Results) -> bool:
    term = result.readrkf("General", "termination status", "ams")
    if term == "NORMAL TERMINATION":
        return True
    elif "NORMAL TERMINATION" in term:
        return "WARNING"
    return


### ==== CALCULATIONS ==== ###


def COSMORS_property(solvent_path: str, solute_path: str, name: str, temperature: float = 298.15) -> float:
    """This method runs a COSMORS property job to obtain the activity coefficient
    which will also calculate G solute which we need to calculate the redox
    potential.
    """

    sett = Settings()
    sett.input.property._h = "ACTIVITYCOEF"
    compounds = [Settings(), Settings()]
    compounds[0]._h = solvent_path
    compounds[1]._h = solute_path
    compounds[0].frac1 = 1
    compounds[1].frac1 = 0

    sett.input.temperature = str(temperature)
    sett.input.compound = compounds

    res = CRSJob(settings=sett, name=name).run().get_results()
    if res:
        # convert the Gibbs energy to hartree (COSMORS gives in kcal/mol)
        return float(Units.convert(res["G solute"][1], "kcal/mol", "hartree"))
    else:
        return False


def calculation_step(
    molecule: Molecule,
    task: str = "GeometryOptimization",
    state: str = "neutral",
    phase: str = "vacuum",
    frequencies: bool = False,
    use_dftb: bool = False,
    use_COSMORS: bool = False,
    name: str = None,
    solvent_path: str = None,
    init_charge: int = 0,
) -> dict:
    """Method used to optimize the geometry of molecule using DFT (by default B3LYP)
    Other settings may be supplied using settings which will be soft-updated using DFT_defaults
    State specifies whether the molecule is neutral, oxidised or reduced
    Phase specifies whether the system is in vacuum or solvated
    if use_dftb, the system will be optimised using DFTB (by default GFN1-xTB) instead of DFT
    """

    if use_COSMORS:
        phase = "solvent"

    settings = get_settings(
        molecule=molecule,
        task=task,
        use_dftb=use_dftb,
        use_COSMORS=use_COSMORS,
        state=state,
        phase=phase,
        frequencies=frequencies,
        init_charge=init_charge,
    )

    # summarize job in one string
    task_abbrev = {"GeometryOptimization": "GO", "SinglePoint": "SP"}[task]
    job_desc = f"{task_abbrev}_{state}_{phase}"
    if use_COSMORS:
        job_desc += "_COSMO-RS"
    if use_dftb:
        job_desc += "_DFTB"

    print(f'\nStarting calculation {name + "_" + job_desc}')
    print(f"\ttask                 = {task}")
    print(f"\tuse_dftb             = {use_dftb}")
    print(f"\tuse_COSMORS          = {use_COSMORS}")
    print(f"\tfrequencies          = {frequencies}")
    print(f"\tstate                = {state}")
    print(f"\tphase                = {phase}")
    print(f"\tcharge               = {settings.input.ams.System.Charge}")
    if not use_dftb:
        print(f"\tunrestricted         = {settings.input.adf.Unrestricted}")
        print(f"\tspin polarization    = {settings.input.adf.SpinPolarization}")

    # run the job
    job = AMSJob(molecule=molecule, settings=settings, name=name + "_" + job_desc)
    res = job.run()

    result_dict = {}
    # pull out results
    if check_termination_succes(res):
        print(f"\tSuccessfull          = {check_termination_succes(res)}")  # True or WARNING
        # set some default values
        bond_energy = None
        gibbs_energy = None

        # If we are doing COSMO calculations then we need to run an additional job to obtain the activity coefficient
        # when calculating the activity coefficient, the G solute is also calculated.
        if use_COSMORS:
            resfile = KFFile(res["adf.rkf"])
            cosmo_data = resfile.read_section("COSMO")
            coskf = KFFile(os.path.join(job.path, "solute.coskf"))
            for k, v in cosmo_data.items():
                coskf.write("COSMO", k, v)
            res.collect()
            bond_energy = res.readrkf("AMSResults", "Energy", "adf")
            solute_path = os.path.join(job.path, "solute.coskf")
            gibbs_energy = COSMORS_property(solvent_path, solute_path, job.name + "_ACTIVITYCOEF")
        # if we dont use COSMO-RS we can just extract the Gibbs and bonding energies from the regular job
        else:
            if use_dftb:
                bond_energy = res.readrkf("AMSResults", "Energy", "dftb")
                if frequencies:
                    gibbs_energy = res.readrkf("Thermodynamics", "Gibbs free Energy", "dftb")
            else:
                bond_energy = res.readrkf("Energy", "Bond Energy", "adf")
                if frequencies:
                    gibbs_energy = res.readrkf("Thermodynamics", "Gibbs free Energy", "adf")

        print("\tResults:")
        if not bond_energy is None:
            result_dict["bond_energy"] = Units.convert(bond_energy, "hartree", "eV")
            print(f'\t\tBond Energy  = {result_dict["bond_energy"]:.4f} eV')
        if not gibbs_energy is None:
            result_dict["gibbs_energy"] = Units.convert(gibbs_energy, "hartree", "eV")
            print(f'\t\tGibbs Energy = {result_dict["gibbs_energy"]:.4f} eV')

        # extract also optimised molecule
        if task == "GeometryOptimization":
            result_dict["geometry"] = res.get_main_molecule()

        # and if the phase is solvent we also need the solvation gibbs energy change
        if phase != "vacuum":
            dG_solvation = res.readrkf("Energy", "Solvation Energy (el)", "adf") + res.readrkf(
                "Energy", "Solvation Energy (cd)", "adf"
            )
            result_dict["dG_solvation"] = Units.convert(dG_solvation, "hartree", "eV")
            print(f'\t\tdG_solvation = {result_dict["dG_solvation"]:.4f} eV')

    else:
        print("\tSuccessfull          = False")

    return result_dict


### ==== MAIN FUNCTION ==== ###


def redox_potential(
    molecule: Molecule,
    mode: str,
    method: str = "screening",
    name: str = None,
    COSMORS_solvent_path: str = None,
    solvent: str = "Dichloromethane",
    init_charge: int = 0,
) -> float:

    assert mode in ["oxidation", "reduction"]
    assert method in [
        "DC",
        "TC-COSMO",
        "TC-COSMO-RS",
        "screening",
    ], 'Argument "method" must be "DC", "TC-COSMO", "TC-COSMO-RS" or "screening"'
    assert os.path.exists(COSMORS_solvent_path), f"Could not find coskf file {COSMORS_solvent_path}"

    # set name
    if name is None:
        name = molecule.properties.name

    print("========================================================================")
    print(f"Starting redox potential calculation for molecule {name}:\n")

    print("\nInitial coordinates:")
    print(molecule)
    print("Settings:")
    print(f"\tName:              {name}")
    print(f"\tMode:              {mode}")
    print(f"\tMethod:            {method}")
    print(f"\tSolvent:           {solvent}")
    print(f"\tInitial Charge:    {init_charge}")

    if mode == "oxidation":
        Gelectron = -0.0375
    elif mode == "reduction":
        Gelectron = 0.0375

    # some settings that every calculation step requires
    general_settings = {
        "name": name,
        "init_charge": init_charge,
    }

    # get on with the actual calculations
    if method == "DC":
        GO_os = calculation_step(molecule, state=mode, phase=solvent, frequencies=True, **general_settings)
        GO_ns = calculation_step(molecule, phase=solvent, frequencies=True, **general_settings)

        redoxpot = GO_os["gibbs_energy"] - GO_ns["gibbs_energy"] + Gelectron

    elif method == "TC-COSMO":
        GO_nv = calculation_step(molecule, frequencies=True, **general_settings)
        SP_nv_ns = calculation_step(GO_nv["geometry"], task="SinglePoint", phase=solvent, **general_settings)
        GO_ns = calculation_step(molecule, phase=solvent, **general_settings)
        SP_nv_nv = calculation_step(GO_ns["geometry"], task="SinglePoint", **general_settings)

        GO_ov = calculation_step(molecule, state=mode, frequencies=True, **general_settings)
        SP_ov_os = calculation_step(
            GO_ov["geometry"], task="SinglePoint", state=mode, phase=solvent, **general_settings
        )
        GO_os = calculation_step(molecule, state=mode, phase=solvent, **general_settings)
        SP_ov_ov = calculation_step(GO_os["geometry"], task="SinglePoint", state=mode, **general_settings)

        redox_part = GO_ov["gibbs_energy"] + SP_ov_os["dG_solvation"] + (SP_ov_ov["bond_energy"] - GO_ov["bond_energy"])
        neutral_part = (
            GO_nv["gibbs_energy"] + SP_nv_ns["dG_solvation"] + (SP_nv_nv["bond_energy"] - GO_nv["bond_energy"])
        )
        redoxpot = redox_part - neutral_part + Gelectron

    elif method == "TC-COSMO-RS":
        assert os.path.exists(COSMORS_solvent_path), f"Solvent database {COSMORS_solvent_path} does not exist"
        GO_nv = calculation_step(molecule, frequencies=True, **general_settings)
        SP_nv_ns = calculation_step(
            GO_nv["geometry"],
            task="SinglePoint",
            use_COSMORS=True,
            solvent_path=COSMORS_solvent_path,
            **general_settings,
        )
        GO_ns = calculation_step(molecule, use_COSMORS=True, solvent_path=COSMORS_solvent_path, **general_settings)
        SP_ns_nv = calculation_step(GO_ns["geometry"], task="SinglePoint", **general_settings)

        GO_ov = calculation_step(molecule, state=mode, frequencies=True, **general_settings)
        SP_ov_os = calculation_step(
            GO_ov["geometry"],
            task="SinglePoint",
            state=mode,
            use_COSMORS=True,
            solvent_path=COSMORS_solvent_path,
            **general_settings,
        )
        GO_os = calculation_step(
            molecule, state=mode, use_COSMORS=True, solvent_path=COSMORS_solvent_path, **general_settings
        )
        SP_os_ov = calculation_step(GO_os["geometry"], task="SinglePoint", state=mode, **general_settings)

        redox_part = GO_ov["gibbs_energy"] + SP_ov_os["dG_solvation"] + (SP_os_ov["bond_energy"] - GO_ov["bond_energy"])
        neutral_part = (
            GO_nv["gibbs_energy"] + SP_nv_ns["dG_solvation"] + (SP_ns_nv["bond_energy"] - GO_nv["bond_energy"])
        )
        redoxpot = redox_part - neutral_part + Gelectron

    elif method == "screening":
        assert os.path.exists(COSMORS_solvent_path), f"Solvent database {COSMORS_solvent_path} does not exist"
        GO_nv = calculation_step(molecule, use_dftb=True, **general_settings)
        COSMO_nv = calculation_step(
            GO_nv["geometry"],
            task="SinglePoint",
            use_COSMORS=True,
            solvent_path=COSMORS_solvent_path,
            **general_settings,
        )

        GO_ov = calculation_step(molecule, state=mode, use_dftb=True, **general_settings)
        COSMO_ov = calculation_step(
            GO_ov["geometry"],
            task="SinglePoint",
            state=mode,
            use_COSMORS=True,
            solvent_path=COSMORS_solvent_path,
            **general_settings,
        )

        redoxpot = COSMO_ov["gibbs_energy"] - COSMO_nv["gibbs_energy"] + Gelectron

    print(f"\nOxidation potential: {redoxpot:.4f} eV")
    return redoxpot


def main(mode):
    job_dir = "./Test"
    solvent = "DMSO"
    mol_files = [
        "./molecules/20410-46-2.xyz",
        "./molecules/109-75-1.xyz",
        "./molecules/167951-80-6.xyz",
        "./molecules/92-52-4.xyz",
    ]
    mol_names = [
        "Phosphorodifluoridate (8CI,9CI)",
        "3-Butenenitrile",
        "3,3,3-Trifluoropropylene carbonate",
        "1,1-Biphenyl",
    ]
    init_charges = [-1, 0, 0, 0]

    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    COSMORS_solvent_path = os.path.abspath(f"coskf/{solvent}.coskf")

    results = {}
    for file, name, charge in zip(mol_files, mol_names, init_charges):
        results[name] = {}
        for method in ["screening", "TC-COSMO", "TC-COSMO-RS", "DC"]:

            job_name = f"{name}_{method}_{mode}"

            # calculation part
            init(path=job_dir, folder=job_name)

            mol = Molecule(file)
            redoxpot = redox_potential(
                mol,
                mode,
                name=name,
                method=method,
                COSMORS_solvent_path=COSMORS_solvent_path,
                solvent=solvent,
                init_charge=charge,
            )
            results[name][method] = redoxpot

            finish()

    ### PRINT POTENTIALS
    print(f"{mode.capitalize()} Potentials:")
    name_len = max(len("System"), max(len(n) for n in results))
    methods = list(results[list(results.keys())[0]])
    method_lens = [max(9, len(m)) for m in methods]
    print(f'{"System".ljust(name_len)} | {" | ".join([m.ljust(l) for m, l in zip(methods, method_lens)])}')
    for name, res in results.items():
        s = f'{name.ljust(name_len)} | {" | ".join([(str(round(res[m],3)) + " V").rjust(l) for m, l in zip(methods, method_lens)])}'
        print(s)

    print("\nCalculations coomplete!\a")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Redox potential workflow.")
    ap.add_argument("molecule", type=str)
    main("oxidation")

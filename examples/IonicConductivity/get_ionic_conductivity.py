import sys
import os
from scm.plams import AMSJob, Molecule, Settings, Units, AMSAnalysisJob
from scm.plams.recipes.assign_charged_ions import assign_charged_ions


def main(filename):
    """
    The main body of the script
    """
    # Get the molecule
    job = AMSJob.load_external(filename)
    mol = job.molecule

    # Assign the charges to the ions
    mol = assign_charged_ions(mol, use_system_charges=False)

    # Read the temperature from KF
    job = AMSJob.load_external(filename)
    T = job.results.readrkf("MDResults", "MeanTemperature")
    kBT = Units.constants["Boltzmann"] * T
    print(f"Average temperature {T:f} K")

    # Get the molecular system and extract the ions
    iontypes, ioncharges, nions, formulas = get_ions(mol)
    print(f"{'Ion':>8s} {'N':>8s} {'Charge':>10s}")
    for k, indices in iontypes.items():
        print(f"{formulas[k]:>8s} {len(indices):8d} {ioncharges[k]:10.5f}")

    # Compute diffusion coefficient for each ion
    diffusion_coeffs = {}
    for label, atoms in iontypes.items():
        s = Settings()
        s.input.Task = "MeanSquareDisplacement"
        s.input.TrajectoryInfo.Trajectory.KFFilename = filename
        atsettings = [iat + 1 for iat in atoms]
        s.input.MeanSquareDisplacement.Atoms.Atom = atsettings

        job = AMSAnalysisJob(settings=s)
        results = job.run()
        D = results._kf.read("Slope(1)", "Final")
        diffusion_coeffs[label] = D
        units = results._kf.read("Slope(1)", "Final(units)")
        print(formulas[label], D, units)

    # Compute the number density for each ion
    rho = {}
    for label, ni in nions.items():
        rho[label] = ni / mol.unit_cell_volume(unit="m")

    # Compute the ionic conductivity
    sigma = 0.0
    for label, D in diffusion_coeffs.items():
        q = ioncharges[label] * Units.constants["electron_charge"]
        s = q**2 * rho[label] * D / kBT
        sigma += s
    return sigma


def get_ions(mol):
    """
    Extract the ions and their charges from the atom properties
    """
    ions = {}
    for i, at in enumerate(mol.atoms):
        regions = []
        if "region" in at.properties:
            regions = at.properties.region
            if not isinstance(regions, set):
                regions = set([regions])
            regions = [s for s in regions if s[:3] == "ion"]
        if len(regions) == 0:
            continue
        name = regions[0]
        if not name in ions:
            ions[name] = []
        ions[name].append(i)

    # Extract the ion info
    iontypes = {}
    ioncharges = {}
    nions = {}
    formulas = {}
    for name, atoms in ions.items():
        charge = 0
        for iat in atoms:
            if "analysis" in mol.atoms[iat].properties:
                if "charge" in mol.atoms[iat].properties.analysis:
                    charge += mol.atoms[iat].properties.analysis.charge
                else:
                    raise PlamsError("Not all charges present in system file")
            elif "forcefield" in mol.atoms[iat].properties:
                if "charge" in mol.atoms[iat].properties.forcefield:
                    charge += mol.atoms[iat].properties.forcefield.charge
                else:
                    raise PlamsError("Not all charges present in system file")
        typename = name.split("_")[0]
        if not typename in iontypes:
            ion = mol.get_fragment(atoms)
            formulas[typename] = ion.get_formula()
            ioncharges[typename] = charge
            iontypes[typename] = []
            nions[typename] = 0
        iontypes[typename] += atoms
        nions[typename] += 1

    return iontypes, ioncharges, nions, formulas


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: amspython get_ionic_conductivity.py path/to/ams.rkf")
        sys.exit(0)
    filename = os.path.abspath(sys.argv[1])
    sigma = main(filename)
    print("Ionic conductivity: %20.10e Siemens/m" % (sigma))

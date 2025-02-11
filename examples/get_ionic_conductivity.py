import sys
from scm.plams import AMSJob, Settings, Units, AMSAnalysisJob, init


def main(filename, chargelines):
    """
    The main body of the script
    """
    # this line is not required in AMS2025+
    init()

    # Get the molecular system
    job = AMSJob.load_external(filename)
    mol = job.molecule
    elements = [at.symbol for at in mol.atoms]

    # Read the temperature from KF
    T = job.results.readrkf("MDResults", "MeanTemperature")
    kBT = Units.constants["Boltzmann"] * T
    print("Average temperaturs %f K" % (T))

    # Read the charges from the input file
    ions = {}
    ioncharges = {}
    nions = {}
    formulas = {}
    for line in chargelines:
        words = line.split()
        if len(words) == 0:
            continue
        q = float(line.split(":")[0])
        atoms = [int(iat) for iat in line.split(":")[-1].split()]
        submol = mol.get_fragment(atoms)
        label = submol.label()
        if not label in ions.keys():
            ions[label] = []
            ioncharges[label] = q * Units.constants["electron_charge"]
            nions[label] = 0
            formulas[label] = submol.separate()[0].get_formula()
        ions[label] += atoms
        nions[label] += len(atoms)

    # Compute diffusion coefficient for each ion
    diffusion_coeffs = {}
    for label, atoms in ions.items():
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
        s = ioncharges[label] ** 2 * rho[label] * D / kBT
        sigma += s
    return sigma


if __name__ == "__main__":
    if len(sys.argv) == 1 or sys.stdin.isatty():
        print("Usage: amspython get_ionic_conductivity.py path/to/ams.rkf < charges.in")
        sys.exit(0)
    chargelines = sys.stdin.readlines()
    filename = sys.argv[1]
    sigma = main(filename, chargelines)
    print("Ionic conductivity: %20.10e Siemens/m" % (sigma))

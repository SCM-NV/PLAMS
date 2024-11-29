import numpy as np

from scm.plams.core.functions import config, delete_job, finish, init
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.amsworker import AMSWorker
from scm.plams.mol.molecule import Molecule
from typing import Optional

__all__ = ["preoptimize", "refine_density", "refine_lattice", "shakemd"]


def preoptimize(
    molecule: Molecule, model: str = "UFF", settings: Settings = None, nproc: int = 1, maxiterations: int = 100
):
    """
    Returns an optimized Molecule (or list of optimized molecules)

    molecule: Molecule or list of Molecules
        Molecule to optimize
    model: str
        Shorthand for some model, e.g. 'UFF'
    settings: Settings
        Custom engine settings (overrides ``model``)
    nproc: int
        Number of processes
    maxiterations: int
        Maximum number of iterations for the geometry optimization.
    """
    my_settings = _get_quick_settings(model, settings, nproc)
    single_molecule = isinstance(molecule, Molecule)

    input_molecules = [molecule] if single_molecule else molecule
    output_molecules = []
    with AMSWorker(my_settings) as worker:
        for i, mol in enumerate(input_molecules):
            results = worker.GeometryOptimization(
                name=f"preoptimization_{i}", molecule=mol, maxiterations=maxiterations, pretendconverged=True
            )
            output_molecules.append(results.get_main_molecule())

    if single_molecule:
        return output_molecules[0]
    else:
        return output_molecules


def shakemd(
    molecule: Molecule,
    density: Optional[float] = None,  # kg/m^3
    nsteps: int = 4000,
    model: str = "UFF",
    settings: Optional[Settings] = None,
    nproc: int = 1,
) -> Molecule:
    """
    Performs SHAKE MD simulation with all bonds fixed. If ``density`` is specified, a lattice deformation will also be applied.

    The temperature is 100 K and the time step 0.5 fs.

    Returns: a Molecule corresponding to the final MD step.

    molecule: Molecule
        The initial structure
    density: float, optional
        If specified, the molecule must have a 3D lattice. Target density in kg/m^3   (1000x the density in g/cm^3)
    nsteps: int
        Number of MD steps.
    model: str
        e.g. 'UFF'
    settings: Settings
        Engine settings (overrides ``model``)
    """
    called_plams_init = _ensure_init()
    s = _get_quick_settings(model, settings, nproc)
    s.input.ams.Task = "MolecularDynamics"
    s.input.ams.MolecularDynamics.NSteps = nsteps
    s.input.ams.MolecularDynamics.Thermostat.Type = "Berendsen"
    s.input.ams.MolecularDynamics.Thermostat.Tau = 100
    s.input.ams.MolecularDynamics.Thermostat.Temperature = 100
    s.input.ams.MolecularDynamics.TimeStep = 0.5
    s.input.ams.MolecularDynamics.InitialVelocities.Temperature = 100
    s.input.ams.MolecularDynamics.Shake.All = "bonds * *"
    if density is not None:
        temp_mol = molecule.copy()
        temp_mol.set_density(density)
        l = temp_mol.lattice
        target_lattice_str = f"""
            {l[0][0]} {l[0][1]} {l[0][2]}
            {l[1][0]} {l[1][1]} {l[1][2]}
            {l[2][0]} {l[2][1]} {l[2][2]}
        """
        s.input.ams.MolecularDynamics.Deformation.StartStep = max(nsteps - 3000, 1)
        s.input.ams.MolecularDynamics.Deformation.StopStep = nsteps
        s.input.ams.MolecularDynamics.Deformation.TargetLattice._1 = target_lattice_str

    job = AMSJob(settings=s, molecule=molecule, name="shakemd")
    job.run()

    out = job.results.get_main_molecule()

    delete_job(job)

    if called_plams_init:
        finish()

    return out


def refine_density(
    molecule: Molecule,
    density: float,
    step_size=50,
    model: str = "UFF",
    settings: Settings = None,
    nproc: int = 1,
    maxiterations: int = 100,
):
    """

    Performs a series of geometry optimizations with densities approaching
    ``density``. This can be useful if you want to compress a system to a
    given density, but cannot just use apply_strain() (because
    apply_strain() also scales the bond lengths).

    This function can be useful if for example packmol does not succeed to
    pack molecules with the desired density. Packmol can then generate a
    structure with a lower density, and this function can be used to
    increase the density to the desired value.

    Returns: a Molecule with the requested density.

    molecule: Molecule
        The molecule must have a 3D lattice
    density: float
        Target density in kg/m^3   (1000x the density in g/cm^3)
    step_size: float
        Step size for the density (in kg/m^3). Set step_size to a large number to only use 1 step.
    model: str
        e.g. 'UFF'
    settings: Settings
        Engine settings (overrides ``model``)
    maxiterations: int
        maximum number of iterations for the geometry optimization.


    """

    assert len(molecule.lattice) == 3
    tolerance = 1e-3
    my_settings = _get_quick_settings(model, settings, nproc)
    current_density = molecule.get_density()  # kg/m^3
    output_molecule = molecule.copy()
    counter = 0
    with AMSWorker(my_settings) as worker:
        while current_density < density - tolerance or current_density > density + tolerance:
            counter += 1
            if current_density < density - step_size:
                new_density = current_density + step_size
            elif current_density > density + step_size:
                new_density = current_density - step_size
            else:
                new_density = density

            output_molecule.set_density(new_density)

            results = worker.GeometryOptimization(
                name=f"preoptimization_{counter}",
                molecule=output_molecule,
                maxiterations=maxiterations,
                pretendconverged=True,
            )

            output_molecule = results.get_main_molecule()
            current_density = output_molecule.get_density()  # kg/m^3

    return output_molecule


def refine_lattice(
    molecule: Molecule,
    lattice,
    n_points=None,
    max_strain=0.15,
    model: str = "UFF",
    settings: Settings = None,
    nproc: int = 1,
    maxiterations: int = 10,
):
    """

    Returns a ``Molecule`` for which the lattice of the ``molecule`` is
    transformed to ``lattice``, by performing short geometry optimizations
    (each for at most ``maxiterations``) on gradually distorted lattices
    (linearly interpolating from the original lattice to the new lattice
    using ``n_points`` points).

    This can be useful for transforming an orthorhombic box of a liquid
    into a non-orthorhombic box of a liquid, where the gradual transformation
    of the lattice ensures that the molecules do not become too distorted.

    If init() has been called before calling this function, the job will be run in
    the current PLAMS working directory and will be deleted when the job finishes.

    Returns: a Molecule with the requested lattice. If the refinement
    fails, ``None`` is returned.

    molecule: Molecule
        The initial molecule
    lattice: list of list of float
        List with 1, 2, or 3 elements. Each element is a list of float with 3 elements each. For example, ``lattice=[[10, 0, 0],[-5, 5, 0],[0, 0, 12]]``.
    n_points: None or int >=2
        Number of points used for the linear interpolation. If None, n_points will be chosen such that the maximum strain for any step is at most ``max_strain`` compared to the original lattice vector lengths.
    max_strain: float
        Only if ``n_points=None``, use this value to determine the maximum allowed strain from one step to the next (as a fraction of the length of the original lattice vectors).
    model: str
        e.g. 'UFF'
    settings: Settings
        Engine settings (overrides ``model``)
    nproc: int
        Number of processes used by the job
    maxiterations: int
        maximum number of iterations for the geometry optimizations.

    """
    assert len(lattice) == len(
        molecule.lattice
    ), f"Different number of lattice vectors: len(molecule.lattice) = {len(molecule.lattice)}, len(lattice) = {len(lattice)}"
    assert all(len(x) == 3 for x in lattice), f"Lattice vectors must have three components. Lattice: {lattice}"
    assert all(len(x) == 3 for x in molecule.lattice), f"Lattice vectors must have three components. Lattice: {lattice}"
    assert (
        len(lattice) >= 1 and len(lattice) <= 3
    ), f"{len(lattice)} lattice vectors given but must be between 1 and 3. Lattice: {lattice}"

    def lattice2str(latt):
        return "\n".join(" ".join(str(j) for j in i) for i in latt)

    if n_points is None:
        lat1 = np.array(molecule.lattice)
        lat2 = np.array(lattice)
        strain = np.linalg.norm(lat2, axis=1) / np.linalg.norm(lat1, axis=1)
        n_points = np.ceil(np.max(np.abs(strain - 1) / max_strain))
        n_points = int(max(n_points, 2))

    assert n_points >= 2, f"Expected n_points >=2, but received {n_points}"

    called_plams_init = _ensure_init()

    s = _get_quick_settings(model, settings, nproc)
    s.input.ams.task = "PESScan"
    s.input.ams.PESScan.ScanCoordinate.FromLattice._1 = lattice2str(molecule.lattice)
    s.input.ams.PESScan.ScanCoordinate.ToLattice._1 = lattice2str(lattice)
    s.input.ams.PESScan.ScanCoordinate.NPoints = n_points
    s.input.ams.GeometryOptimization.MaxIterations = maxiterations
    s.input.ams.GeometryOptimization.PretendConverged = "Yes"

    job = AMSJob(settings=s, molecule=molecule, name="refine_density")
    job.run()

    final_molecule = None
    if job.ok():
        results = job.results.get_pesscan_results(molecules=True)
        final_molecule = results["Molecules"][-1]

    delete_job(job)

    if called_plams_init:
        finish()

    return final_molecule


def _ensure_init():
    if config.init:
        called_plams_init = False
    else:
        s = Settings()
        s.erase_workdir = True
        s.log.stdout = 0
        init(config_settings=s)
        called_plams_init = True

    return called_plams_init


def model_to_settings(model: str):
    """
    Returns Settings
    """
    settings = Settings()
    model = model.lower()
    if model == "uff":
        settings.input.ForceField.Type = "UFF"
    elif model == "gaff":
        settings.input.ForceField.Type = "GAFF"
        settings.input.ForceField.AnteChamberIntegration = "Yes"
    elif model == "gfnff":
        settings.input.GFNFF
    elif model == "ani-2x" or model == "ani2x":
        settings.input.MLPotential.Model = "ANI-2x"
    elif model == "gfn1xtb" or model == "gfn1-xtb":
        settings.input.DFTB.Model = "GFN1-xTB"
    elif model == "m3gnet-up-2022" or model == "m3gnetup2022":
        settings.input.MLPotential.Model = "M3GNet-UP-2022"
    else:
        raise ValueError("Unknown model: {}".format(model))

    return settings


def _get_quick_settings(model, settings, nproc):
    if settings is None:
        my_settings = model_to_settings(model)
    else:
        my_settings = settings.copy()

    if nproc:
        my_settings.runscript.nproc = nproc

    return my_settings

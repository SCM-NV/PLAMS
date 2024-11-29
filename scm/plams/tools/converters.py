import os
import re
import tempfile
from typing import Optional

from scm.plams.interfaces.molecule.ase import toASE
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.kftools import KFFile
from scm.plams.tools.units import Units
from scm.plams.core.functions import requires_optional_package
from scm.plams.trajectories.rkffile import RKFTrajectoryFile
from scm.plams.trajectories.rkfhistoryfile import RKFHistoryFile

__all__ = [
    "traj_to_rkf",
    "vasp_output_to_ams",
    "qe_output_to_ams",
    "gaussian_output_to_ams",
    "rkf_to_ase_traj",
    "rkf_to_ase_atoms",
    "file_to_traj",
]


@requires_optional_package("ase")
def traj_to_rkf(trajfile, rkftrajectoryfile, task=None, timestep: float = 0.25):
    """
    Convert ase .traj file to .rkf file. NOTE: The order of atoms (or the number of atoms) cannot change between frames!

    trajfile : str
        path to a .traj file
    rkftrajectoryfile : str
        path to the output .rkf file (will be created)
    task : str
        Which task to write. If None it is auto-determined.

    timestep: float
        Which timestep to write when task == 'moleculardynamics'


    Returns : 2-tuple (coords, cell)
        The final coordinates and cell in angstrom
    """
    import warnings

    from ase.io import Trajectory

    warnings.filterwarnings("ignore", "Creating an ndarray from ragged nested sequences")
    traj = Trajectory(trajfile)
    rkfout = RKFTrajectoryFile(rkftrajectoryfile, mode="wb")
    rkfout.store_historydata()
    rkfout.store_mddata()

    # Store the units properly
    mdunits = {}
    mdunits["PotentialEnergy"] = "hartree"
    mdunits["KineticEnergy"] = "hartree"
    mdunits["TotalEnergy"] = "hartree"
    rkfout._set_mdunits(mdunits)

    energy_converter = Units.convert(1.0, "eV", "hartree")
    gradients_converter = Units.convert(1.0, "eV/angstrom", "hartree/bohr")
    stress_converter = Units.convert(1.0, "eV/angstrom^3", "hartree/bohr^3")

    coords, cell = None, None
    try:
        for i, atoms in enumerate(traj):
            if i == 0:
                rkfout.set_elements(atoms.get_chemical_symbols())
            coords = atoms.get_positions()  # angstrom
            cell = atoms.get_cell()
            try:
                gradients = -atoms.get_forces() * gradients_converter
            except:
                gradients = None
            try:
                stresstensor = atoms.get_stress(voigt=False) * stress_converter
            except:
                stresstensor = None
            mddata = {}
            try:
                energy = atoms.get_potential_energy() * energy_converter
                mddata["PotentialEnergy"] = energy
            except:
                energy = None

            try:
                kinetic_energy = atoms.get_kinetic_energy() * energy_converter
                mddata["KineticEnergy"] = kinetic_energy
            except:
                kinetic_energy = None

            if "PotentialEnergy" in mddata and "KineticEnergy" in mddata:
                mddata["TotalEnergy"] = mddata["PotentialEnergy"] + mddata["KineticEnergy"]

            if str(task).lower() == "moleculardynamics" or len(mddata) > 0:
                mddata["Time"] = timestep * i

            final_mddata = None if len(mddata) == 0 else mddata

            # Create a historydata dictionary, to go into the History section
            historydata = {}
            if gradients is not None:
                historydata["Gradients"] = gradients
                historydata["EngineGradients"] = gradients
            if stresstensor is not None:
                historydata["StressTensor"] = stresstensor
            if len(historydata) == 0:
                historydata = {}

            rkfout.write_next(coords=coords, cell=cell, historydata=historydata, mddata=final_mddata)

    finally:
        rkfout.close()

    # the below is needed to be able to load the .rkf file with AMSJob.load_external()
    kf = KFFile(rkftrajectoryfile)
    kf["EngineResults%nEntries"] = 0
    kf["General%program"] = "ams"
    if task is None:
        if len(traj) == 1:
            task = "singlepoint"
        elif kinetic_energy is None or kinetic_energy == 0:
            task = "geometryoptimization"
        else:
            task = "moleculardynamics"
    kf["General%task"] = task
    kf["General%user input"] = "\xFF".join([f"Task {task}", "Engine External", "EndEngine"])

    return coords, cell


@requires_optional_package("ase")
def file_to_traj(outfile, trajfile):
    """
    outfile : str
        path to existing file (OUTCAR, qe.out, etc.)
    trajfile : str
        will be created
    """
    from ase.io import read, write

    if os.path.exists(trajfile):
        os.remove(trajfile)

    atoms = read(outfile, ":")
    write(trajfile, atoms)

    if not os.path.exists(trajfile):
        raise RuntimeError("Couldn't write {}".format(trajfile))

    return trajfile


def _remove_or_raise(file, overwrite):
    if os.path.exists(file):
        if overwrite:
            os.remove(file)
        else:
            raise RuntimeError("{} already exists, specify overwrite=True to overwrite".format(file))


def _write_engine_rkf(kffile, enginefile):
    kf = KFFile(kffile)
    enginerkf = KFFile(enginefile)
    # write engine.rkf
    # copy General, Molecule, InputMolecule from ams.rkf
    for sec in ["General", "Molecule", "InputMolecule"]:
        secdict = kf.read_section(sec)
        for k, v in secdict.items():
            enginerkf[sec + "%" + k] = v
    enginerkf["General%program"] = "plams"
    nEntries = kf["History%nEntries"]
    suffix = "({})".format(nEntries)
    if ("History", "Energy" + suffix) in kf:
        enginerkf["AMSResults%Energy"] = kf["History%Energy" + suffix]
    if ("History", "Gradients" + suffix) in kf:
        enginerkf["AMSResults%Gradients"] = kf["History%Gradients" + suffix]
    if ("History", "StressTensor" + suffix) in kf:
        enginerkf["AMSResults%StressTensor"] = kf["History%StressTensor" + suffix]


def _postprocess_vasp_amsrkf(kffile, outcar):
    # add extra info to the kffile
    kf = KFFile(kffile, autosave=False)
    try:
        kf["EngineResults%nEntries"] = 1
        kf["EngineResults%Title(1)"] = "vasp"
        kf["EngineResults%Description(1)"] = "Standalone VASP run. Data from {}".format(os.path.abspath(outcar))
        kf["EngineResults%Files(1)"] = "vasp.rkf"
        kf["General%user input"] = "!VASP"

        # read the INCAR
        incarfile = os.path.join(os.path.dirname(outcar), "INCAR")
        userinput = ["!VASP", "Engine External", "  Input", "  !INCAR"]
        if os.path.exists(incarfile):
            with open(incarfile) as incar:
                for line in incar:
                    line = line.split("!")[0]
                    line = line.split("#")[0]
                    line = line.strip()
                    if line.lower().startswith("end"):  # "End" is reserved to end the block
                        line = "!" + line
                    if len(line) > 0:
                        userinput.append("    " + line)
            userinput.append("  !EndINCAR")
        userinput.append("  EndInput")  # end of the Free block
        userinput.append("EndEngine")
        userinput.append("Task {}".format(kf["General%task"]))
        kf["General%user input"] = "\xFF".join(userinput)

    finally:
        kf.save()


def vasp_output_to_ams(
    vasp_folder: str,
    wdir: Optional[str] = None,
    overwrite: bool = False,
    write_engine_rkf: bool = True,
    task: Optional[str] = None,
    timestep: float = 0.25,
):
    """
    Converts VASP output (OUTCAR, ...) to AMS output (ams.rkf, vasp.rkf)

    Returns: a string containing the directory where ams.rkf was written

    vasp_folder : str
        path to a directory with an OUTCAR, INCAR, POTCAR etc. files

    wdir : str or None
        directory in which to write the ams.rkf and vasp.rkf files
        If None, a subdirectory "AMSJob" of vasp_folder will be created

    overwrite : bool
        if False, first check if wdir already contains ams.rkf and vasp.rkf, in which case do nothing
        if True, overwrite if exists

    write_engine_rkf : bool
        If True, also write vasp.rkf alongside ams.rkf. The vasp.rkf file will only contain an AMSResults section (energy, gradients, stress tensor). It will not contain the DOS or the band structure.

    task : str
        Which task to write to ams.rkf. If None it is auto-determined (probably set to 'geometryoptimization')

    timestep : float
        If task='moleculardynamics', which timestep (in fs) between frames to write
    """
    if not os.path.isdir(vasp_folder):
        raise ValueError("Directory {} does not exist".format(vasp_folder))

    outcar = os.path.join(vasp_folder, "OUTCAR")
    if not os.path.exists(outcar):
        if os.path.exists(os.path.join(vasp_folder, "XDATCAR")):
            outcar = os.path.join(vasp_folder, "XDATCAR")
        else:
            raise ValueError("File {} does not exist, should be an OUTCAR file.".format(outcar))

    if wdir is None:
        wdir = os.path.join(os.path.dirname(outcar), "AMSJob")
        os.makedirs(wdir, exist_ok=True)

    # exit early if ams.rkf already exists
    if os.path.exists(os.path.join(wdir, "ams.rkf")) and not overwrite:
        return wdir

    # convert OUTCAR to a .traj file inside wdir
    trajfile = file_to_traj(outcar, os.path.join(wdir, "vasp.traj"))

    # remove the target files first if overwrite
    kffile = os.path.join(wdir, "ams.rkf")
    enginefile = os.path.join(wdir, "vasp.rkf")
    _remove_or_raise(kffile, overwrite)
    _remove_or_raise(enginefile, overwrite)

    # convert the .traj file to ams.rkf
    traj_to_rkf(trajfile, kffile, task=task, timestep=timestep)

    _postprocess_vasp_amsrkf(kffile, outcar)
    if write_engine_rkf:
        _write_engine_rkf(kffile, enginefile)

    if os.path.exists(trajfile):
        os.remove(trajfile)

    return wdir


def _postprocess_qe_amsrkf(kffile, qe_outfile):
    # add extra info to the kffile
    kf = KFFile(kffile, autosave=False)
    try:
        kf["EngineResults%nEntries"] = 1
        kf["EngineResults%Title(1)"] = "qe"
        kf["EngineResults%Description(1)"] = "Standalone Quantum ESPRESSO run. Data from {}".format(
            os.path.abspath(qe_outfile)
        )
        kf["EngineResults%Files(1)"] = "qe.rkf"

        userinput = [
            "!QuantumESPRESSO",
            "Engine External",
            "  Input",
            "    Unknown Quantum ESPRESSO input",
            "  EndInput",
            "EndEngine",
        ]
        kf["General%user input"] = "\xFF".join(userinput)

    finally:
        kf.save()


def _postprocess_gaussian_amsrkf(kffile, gaussian_outfile):
    # add extra info to the kffile
    kf = KFFile(kffile, autosave=False)
    try:
        kf["EngineResults%nEntries"] = 1
        kf["EngineResults%Title(1)"] = "gaussian"
        kf["EngineResults%Description(1)"] = "Standalone Gaussian. Data from {}".format(
            os.path.abspath(gaussian_outfile)
        )
        kf["EngineResults%Files(1)"] = "gaussian.rkf"

        userinput = ["!Gaussian", "Engine External", "  Input", "    Unknown Gaussian input", "  EndInput", "EndEngine"]
        kf["General%user input"] = "\xFF".join(userinput)

    finally:
        kf.save()


def text_out_file_to_ams(qe_outfile, wdir=None, overwrite=False, write_engine_rkf=True, enginename="qe"):
    """
    Converts a qe .out or gaussian .out file to ams.rkf and qe.rkf/gaussian.rkf

    Do not use this function directly, instaead call qe_output_to_ams or gaussian_output_to_ams
    """
    if not os.path.exists(qe_outfile) or os.path.isdir(qe_outfile):
        raise FileNotFoundError(qe_outfile)

    basename = os.path.basename(qe_outfile)
    basename_no_suffix = basename
    if basename.endswith(".out"):
        basename_no_suffix = re.sub(".out$", "", basename)
    dirname = os.path.abspath(os.path.dirname(qe_outfile))

    if wdir is None:
        if os.path.isdir(os.path.join(dirname, basename_no_suffix + ".results")):
            wdir = os.path.join(dirname, basename_no_suffix + ".results", "AMSJob")
        else:
            wdir = os.path.join(dirname, basename_no_suffix + ".AMSJob")

    if os.path.exists(os.path.join(wdir, "ams.rkf")) and not overwrite:
        return wdir

    os.makedirs(wdir, exist_ok=True)

    # convert to a .traj file inside wdir
    # first trim the qe_outfile to the first occurrence of "JOB DONE". This is needed because
    # running standalone QE via the AMS GUI will print multiple jobs into the same output file
    # i.e. both the geo opt and band structure calculation into the same file, which causes
    # the ASE qe.out parser to crash
    with tempfile.NamedTemporaryFile(delete=False) as tf:
        with open(qe_outfile, "r") as instream:
            for line in instream:
                tf.write(line.encode())
                if "JOB DONE" in line:
                    break
        tf.flush()
        tf.file.close()
        trajfile = file_to_traj(tf.name, os.path.join(wdir, "out.traj"))
        os.remove(tf.name)

    # remove the target files first if overwrite
    kffile = os.path.join(wdir, "ams.rkf")
    enginefile = os.path.join(wdir, f"{enginename}.rkf")
    _remove_or_raise(kffile, overwrite)
    _remove_or_raise(enginefile, overwrite)

    # convert the .traj file to ams.rkf
    traj_to_rkf(trajfile, kffile)

    if write_engine_rkf:
        # here one could also run $AMSBIN/tokf to get things like the DOS
        # $AMSBIN/tokf qe qe_outfile enginefile
        # But that would need to reread the entire trajectory
        # and one would need to postprocess it with the AMSResults section
        _write_engine_rkf(kffile, enginefile)

    if os.path.exists(trajfile):
        os.remove(trajfile)

    return wdir


def qe_output_to_ams(qe_outfile, wdir=None, overwrite=False, write_engine_rkf=True):
    """
    Converts a qe .out file to ams.rkf and qe.rkf.

    Returns: a string containing the directory where ams.rkf was written

    If the filename ends in .out, check if a .results directory exists. In that case, place
    the AMSJob subdirectory in the .results directory.

    Otherwise, create a new directory called filename.AMSJob

    qe_outfile : str
        path to the qe output file

    """
    wdir = text_out_file_to_ams(
        qe_outfile, wdir, overwrite=overwrite, write_engine_rkf=write_engine_rkf, enginename="qe"
    )
    _postprocess_qe_amsrkf(os.path.join(wdir, "ams.rkf"), qe_outfile)
    return wdir


def gaussian_output_to_ams(outfile, wdir=None, overwrite=False, write_engine_rkf=True):
    """
    Converts a Gaussian .out file to ams.rkf and gaussian.rkf.

    Returns: a string containing the directory where ams.rkf was written

    If the filename ends in .out, check if a .results directory exists. In that case, place
    the AMSJob subdirectory in the .results directory.

    Otherwise, create a new directory called filename.AMSJob

    outfile : str
        path to the gaussian output file

    """
    wdir = text_out_file_to_ams(
        outfile, wdir, overwrite=overwrite, write_engine_rkf=write_engine_rkf, enginename="gaussian"
    )
    _postprocess_gaussian_amsrkf(os.path.join(wdir, "ams.rkf"), outfile)
    return wdir


@requires_optional_package("ase")
def rkf_to_ase_atoms(rkf_file, get_results=True):
    """
    Convert an ams.rkf trajectory to a list of ASE atoms

    rkf_file: str
        Path to an ams.rkf file

    get_results: bool
        Whether to include results like energy, forces, and stress in the trajectory.

    Returns: a list of all the ASE Atoms objects.
    """
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator
    import numpy as np

    bohr2angstrom = Units.convert(1.0, "bohr", "angstrom")
    hartree2eV = Units.convert(1.0, "hartree", "eV")

    def get_ase_atoms(elements, crd, cell, energy, gradients, stress):

        pbc = None
        if cell is not None:
            cell = np.array(cell).reshape(-1, 3)
            pbc = ["T"] * len(cell) + ["F"] * (3 - len(cell))
        atoms = Atoms(symbols=elements, positions=np.array(crd).reshape(-1, 3), cell=cell, pbc=pbc)
        if get_results:
            calculator = SinglePointCalculator(atoms)
            atoms.set_calculator(calculator)

            if energy:
                atoms.calc.results["energy"] = energy * hartree2eV
            if gradients:
                forces = -np.array(gradients).reshape(-1, 3) * hartree2eV / bohr2angstrom
                atoms.calc.results["forces"] = forces
            if stress:
                n = len(stress)
                if n == 9:
                    stress = np.array(stress).reshape(3, 3) * hartree2eV / bohr2angstrom**3
                    atoms.calc.results["stress"] = np.array(
                        [stress[0][0], stress[1][1], stress[2][2], stress[1][2], stress[0][2], stress[0][1]]
                    )

        return atoms

    rkf_filename = rkf_file
    kf = KFFile(rkf_filename)
    if "History" in kf.keys():
        if "ChemicalSystem(1)" in kf.keys():
            rkf = RKFHistoryFile(rkf_filename)
        else:
            rkf = RKFTrajectoryFile(rkf_filename)

        rkf.store_historydata()
        all_atoms = []
        for crd, cell in rkf:
            energy, stress = None, None
            if get_results:
                energy = rkf.historydata.get("EngineEnergy", None)
                if energy is None:
                    energy = rkf.historydata.get("Energy", None)
                gradients = rkf.historydata.get("EngineGradients", None)
                if gradients is None:
                    gradients = rkf.historydata.get("Gradients", None)
                stress = rkf.historydata.get("StressTensor", None)
            atoms = get_ase_atoms(rkf.elements, crd, cell, energy, gradients, stress)
            all_atoms.append(atoms)
    else:
        atoms = toASE(Molecule(rkf_filename))
        all_atoms = [atoms]

    return all_atoms


@requires_optional_package("ase")
def rkf_to_ase_traj(rkf_file, out_file, get_results=True):
    """
    Convert an ams.rkf trajectory to a different trajectory format (.xyz, .traj, anything supported by ASE)

    rkf_file: str
        Path to an ams.rkf file

    out_file: str
        Path to the .traj or .xyz file that will be created. If the file exists it will be overwritten. If a .xyz file is specified it will use the normal ASE format (not the AMS format).

    get_results: bool
        Whether to include results like energy, forces, and stress in the trajectory.

    """

    from ase.io import write

    all_atoms = rkf_to_ase_atoms(rkf_file, get_results=get_results)
    write(out_file, all_atoms)
    return all_atoms

#!/usr/bin/env amspython
from scm.plams import AMSJob, Atom, Molecule, Settings, finish, init


def main():
    init()

    mol = get_system()
    s = get_settings()

    # Note: do NOT give the job the name "quantumespresso", as that will cause file name clashes
    job = AMSJob(molecule=mol, settings=s, name="afm-FeO")
    print(job.get_input())

    job.run()

    print("Final energy (hartree)")
    print(job.results.get_energy())

    finish()


def get_settings() -> Settings:
    s = Settings()

    # s.runscript.preamble_lines = ['export SCM_QESTARTCMD="custom_start_command"']

    s.input.ams.Task = "SinglePoint"

    # Defining the Pseudopotentials
    s.input.QuantumEspresso.Pseudopotentials.Family = "SSSP-Efficiency"
    s.input.QuantumEspresso.Pseudopotentials.Functional = "PBE"

    # Use ._h and ._1 for unstructured blocks like K_Points
    s.input.QuantumEspresso.K_Points._h = "automatic"
    s.input.QuantumEspresso.K_Points._1 = "6 6 6 0 0 0"

    # Hubbard U specified in the new QE 7.1 format (different from QE 7.0)
    s.input.QuantumEspresso.Hubbard._h = "ortho-atomic"
    s.input.QuantumEspresso.Hubbard._1 = """
        U Fe1-3d 4.6
        U Fe2-3d 4.6
    """

    # When initializing many keys inside the same block
    # it can be convenient to group them like this
    s.input.QuantumEspresso.System = Settings(
        ecutwfc=40.0,
        ecutrho=240.0,
        occupations="smearing",
        smearing="gaussian",
        degauss=0.02,
        nspin=2,
        starting_magnetization=[Settings(Label="Fe1", Value=1.0), Settings(Label="Fe2", Value=-1.0)],
    )

    # You may also just use the normal PLAMS Settings dot notation
    s.input.QuantumEspresso.Electrons.conv_thr = 1.0e-8
    s.input.QuantumEspresso.Electrons.mixing_beta = 0.3

    return s


def get_system() -> Molecule:
    """

    Returns a PLAMS Molecule for FeO with the QE.Label properties set to 'Fe1'
    and 'Fe2' for the two Fe atoms

    """

    d = 4.33
    mol = Molecule()
    mol.add_atom(Atom(symbol="Fe", coords=(0, 0, 0)))
    mol.add_atom(Atom(symbol="Fe", coords=(d / 2, d / 2, 0)))
    mol.add_atom(Atom(symbol="O", coords=(0, d / 2, 0)))
    mol.add_atom(Atom(symbol="O", coords=(d / 2, 0, 0)))

    mol.lattice = [[d, 0, 0], [0, d, 0], [d / 2, 0, d / 2]]

    mol[1].properties.QE.Label = "Fe1"
    mol[2].properties.QE.Label = "Fe2"

    return mol


if __name__ == "__main__":
    main()

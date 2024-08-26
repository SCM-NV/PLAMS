try:
    import dill as pickle  # noqa: F401
except ImportError:
    pass

from scm.plams import AMSJob, Settings
from scm.plams.unit_tests.test_helpers import skip_if_no_ams_installation


def test_hybrid_engine_input():
    """test :meth:`AMSJob.get_input` for writing sub engines in a hybrid engine block of ams."""
    AMSinput = """system
  Atoms
              C      -1.6447506665       1.4391568332       0.0000000000 
              C      -0.5773632247       0.3290989412      -0.0074515501 
              H      -1.2586274675       2.3064141938       0.5149500419 
              H      -2.5299121499       1.0840525749       0.5067446240 
              H      -1.8942698087       1.7054858887      -1.0164689035 
              H       0.3077982588       0.6842031995      -0.5141961741 
              H      -0.9634864236      -0.5381584194      -0.5224015919 
              H      -0.3278440824       0.0627698856       1.0090173534 
  End
  BondOrders
     2 1 1.0
     3 1 1.0
     4 1 1.0
     5 1 1.0
     6 2 1.0
     7 2 1.0
     8 2 1.0
  End
End

task SinglePoint

Engine hybrid
  energy
    term
      engineid ADF-1
      factor 0.5
      region *
    End
    term
      engineid DFTB-1
      factor 0.5
      region *
    End
  End
  engine ADF ADF-1
  EndEngine
  engine DFTB DFTB-1
  EndEngine

EndEngine

"""
    skip_if_no_ams_installation()
    job = AMSJob.from_input(AMSinput)
    assert job.get_input() == AMSinput


def test_list_block_input():
    """test :meth:`AMSJob.get_input` for writing list-blocks"""
    AMSinput = """system
  Atoms
              C      -1.6447506665       1.4391568332       0.0000000000 
  End
End

task SinglePoint

Engine hybrid
  energy
    term
      engineid ADF-1
      factor 0.5
      region *
    End
    term
      engineid ADF-2
      factor 0.5
      region *
    End
    term
      engineid DFTB-1
      factor 0.5
      region *
    End
  End
  engine ADF ADF-1
     Save TAPE10
     Save TAPE11
     Basis
     PerAtomType
     Symbol H
     End
     PerAtomType
     Symbol O
     End
     End
  EndEngine
  engine ADF ADF-2
  EndEngine
  engine DFTB DFTB-1
  EndEngine

EndEngine

"""
    skip_if_no_ams_installation()
    job = AMSJob.from_input(AMSinput)
    assert job.get_input() == AMSinput


def test_settings_to_molecule():
    s = Settings()
    s.input.ams.system.Atoms._1 = "         H       0.0000000000       0.0000000000       0.0000000000 "
    s.input.ams.system.Atoms._2 = "         O       1.0000000000       0.0000000000       0.0000000000 "

    t = s.copy()
    mol = AMSJob.settings_to_mol(t)[""]

    assert AMSJob(settings=s).get_input() == AMSJob(settings=t, molecule=mol).get_input()

    s = Settings()
    s.input.ams.system.Atoms._1 = [
        "         H       0.0000000000       0.0000000000       0.0000000000 ",
        "         O       1.0000000000       0.0000000000       0.0000000000 ",
    ]

    t = s.copy()
    mol = AMSJob.settings_to_mol(t)[""]

    assert AMSJob(settings=s).get_input() == AMSJob(settings=t, molecule=mol).get_input()

    # test for headers:
    s.input.ams.system.Atoms._h = "[A]"
    s.input.ams.system.Atoms._1 = [
        "         O      -0.0509000000      -0.2754000000       0.6371000000 ForceField.Charge=-0.834 ForceField.Type=OW",
        "         H       0.0157000000       0.5063000000       0.0531000000 ForceField.Charge=0.417 ForceField.Type=HW",
        "         H      -0.0055000000      -1.0411000000       0.0658000000 ForceField.Charge=0.417 ForceField.Type=HW",
        "         O       0.0981000000       1.7960000000      -1.2550000000 ForceField.Charge=-0.834 ForceField.Type=OW",
        "         H      -0.6686000000       2.2908000000      -1.5343000000 ForceField.Charge=0.417 ForceField.Type=HW",
        "         H       0.8128000000       2.3488000000      -1.5619000000 ForceField.Charge=0.417 ForceField.Type=HW",
    ]
    t = s.copy()
    # currently plams completely ignores the [A] unit specifier, it might also not print it if it does
    # get implemented and instead just convert the values in the molecule.
    del s.input.ams.system.Atoms._h

    mol = AMSJob.settings_to_mol(t)[""]
    assert AMSJob(settings=s).get_input() == AMSJob(settings=t, molecule=mol).get_input()

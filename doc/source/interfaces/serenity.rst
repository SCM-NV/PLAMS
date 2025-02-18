Serenity
-------------------------

.. currentmodule:: scm.plams.interfaces.thirdparty.serenity

Serenity is a quantum chemistry code originally developed in the group of Johannes Neugebauer at the University of Münster with a strong focus on quantum chemical subsystem/embedding methods. See the `serenity page <https://github.com/qcserenity/serenity>`_ for more details.

PLAMS offers a simple Serenity interface and is capable of running Serenity calculations.
The relevant classes are |SerenityJob|, |SerenitySettings|, and |SerenityResults|.

It is also possible to run post-Hartree-Fock (HF) Serenity calculations starting from a ground-state Self-Consistent-Field calculation perfromed with ADF.


Serenity calculation
~~~~~~~~~~~~~~~~~~~~~~~~~

Preparing an instance of a |SerenityJob| follows general principles for |SingleJob|.
Here is an example of a simple HF calculation::
  
    from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob
    
    sersett = SerenitySettings()
    # The path to the geometry file
    sersett.input.system.gly_gly_gly.geometry = "mygeometry.xyz"
    # Charge and spin
    sersett.input.system.gly_gly_gly.charge = "0"
    sersett.input.system.gly_gly_gly.spin = "0"
    # The electronic structure method
    sersett.input.system.gly_gly_gly.method = "HF"
    sersett.input.system.gly_gly_gly.basis.label = "6-31G"
    # Performs a SCF calculation 
    sersett.input.task.SCF.act = "mymolecule"
    
    serjob = SerenityJob(settings=sersett, name="Serenitycalc")
    serjob.run()

This input performs an HF/6-31G calculation on the chosen molecule.


ADF interface
~~~~~~~~~~~~~~~~~~~~~~~~~

Interfacing ADF with Serenity is done through a special file, saved by ADF, which contains information about the Hamiltonian matrix which can be used by Serenity for post-HF calculations.
Here is a simple example::


    from scm.plams import init, Molecule, Settings, AMSJob
    from scm.plams.interfaces.thirdparty.serenity import SerenitySettings, SerenityJob
    
    init(folder="PLAMS_CCSD_N2_DZ")
    
    mol = Molecule("N2.xyz")
    adfsett = Settings()
    adfsett.input.ams.task = "SinglePoint"
    adfsett.input.adf.IntegralsToFile = "SERENITY"
    adfsett.input.adf.TotalEnergy = ""
    adfsett.input.adf.basis.core = "None"
    adfsett.input.adf.basis.type = "DZ"
    adfsett.input.adf.basis.CreateOutput = "yes"
    adfsett.input.adf.relativity = "Level=None"
    adfsett.input.adf.symmetry = "NoSym"
    adfsett.input.adf.xc.hartreefock = ""
    
    adfjob = AMSJob(molecule=mol, settings=adfsett, name="ADF_N2_HF")
    adfjob.run()
    bonding_energy = adfjob.results.get_energy()
    print(f"ADF bonding energy: {bonding_energy} hartree")
    
    
    sersett = SerenitySettings()
    sersett.input.system.N2.geometry = "N2.xyz"
    
    sersett.input.task.CC.system = "N2"
    sersett.input.task.CC.level = "CCSD"
    
    serjob = SerenityJob(settings=sersett, name="Serenity_N2_CCSD")
    serjob.run()
    ccsd_correction = serjob.results.get_ccsd_energy_correction()
    print(f"Serenity CCSD energy correction: {ccsd_correction} hartree")


This script produces the CCSD/DZ energy of the nitrogen molecules using integrals computed by ADF with a Slater-type orbital basis.

API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: SerenityJob
    :exclude-members: _result_type, _get_ready
.. autoclass:: SerenityResults

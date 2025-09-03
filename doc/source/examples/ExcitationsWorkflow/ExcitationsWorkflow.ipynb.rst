Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

   from scm.plams import AMSResults, Units, Settings, read_molecules, AMSJob, init

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/examples/ExcitationsWorkflow/plams_workdir

Helper Functions
~~~~~~~~~~~~~~~~

Set up a couple of useful functions for extracting results.

.. code:: ipython3

   def get_excitations(results):
       """Returns excitation energies (in eV) and oscillator strengths (in Debye)."""
       if results.job.ok():
           exci_energies_au = results.readrkf("Excitations SS A", "excenergies", file="engine")
           oscillator_str_au = results.readrkf("Excitations SS A", "oscillator strengths", file="engine")
           # The results are stored in atomic units. Convert them to more convenient units:
           exci_energies = Units.convert(exci_energies_au, "au", "eV")
           oscillator_str = Units.convert(oscillator_str_au, "au", "Debye")
           return exci_energies, oscillator_str
       else:
           return [], []

.. code:: ipython3

   def has_good_excitations(results, min_energy, max_energy, oscillator_str_threshold=1e-4):
       """Returns True if there is at least one excitation with non-vanishing oscillator strength
       in the energy range [min_energy, max_energy]. Unit for min_energy and max energy: eV."""
       exci_energies, oscillator_str = get_excitations(results)
       for e, o in zip(exci_energies, oscillator_str):
           if min_energy < e < max_energy and o > oscillator_str_threshold:
               return True
       return False

Calculation settings
~~~~~~~~~~~~~~~~~~~~

Configure the settings for the various jobs.

.. code:: ipython3

   # Settings for geometry optimization with the AMS driver:
   go_sett = Settings()
   go_sett.input.ams.Task = "GeometryOptimization"
   go_sett.input.ams.GeometryOptimization.Convergence.Gradients = 1.0e-4

.. code:: ipython3

   # Settings for single point calculation with the AMS driver
   sp_sett = Settings()
   sp_sett.input.ams.Task = "SinglePoint"

.. code:: ipython3

   # Settings for the DFTB engine (including excitations)
   dftb_sett = Settings()
   dftb_sett.input.dftb.Model = "SCC-DFTB"
   dftb_sett.input.dftb.ResourcesDir = "QUASINANO2015"
   dftb_sett.input.dftb.Properties.Excitations.TDDFTB.calc = "singlet"
   dftb_sett.input.dftb.Properties.Excitations.TDDFTB.lowest = 10
   dftb_sett.input.dftb.Occupation.Temperature = 5.0

.. code:: ipython3

   # Settings for the geometry optimization with the ADF engine
   adf_sett = Settings()
   adf_sett.input.adf.Basis.Type = "DZP"
   adf_sett.input.adf.NumericalQuality = "Basic"

.. code:: ipython3

   # Settings for the excitation calculation using the ADF engine
   adf_exci_sett = Settings()
   adf_exci_sett.input.adf.Basis.Type = "TZP"
   adf_exci_sett.input.adf.XC.GGA = "PBE"
   adf_exci_sett.input.adf.NumericalQuality = "Basic"
   adf_exci_sett.input.adf.Symmetry = "NoSym"
   adf_exci_sett.input.adf.Excitations.lowest = 10
   adf_exci_sett.input.adf.Excitations.OnlySing = ""

Load Molecules
~~~~~~~~~~~~~~

Import all xyz files in the folder ‘molecules’.

.. code:: ipython3

   molecules = read_molecules("molecules")

DFTB Prescreen
~~~~~~~~~~~~~~

Perform an initial prescreen of all molecules with DFTB.

.. code:: ipython3

   promising_molecules = {}

.. code:: ipython3

   for name, mol in molecules.items():
       dftb_job = AMSJob(name="DFTB_" + name, molecule=mol, settings=go_sett + dftb_sett)
       dftb_results = dftb_job.run()

       if has_good_excitations(dftb_results, 1, 6):
           promising_molecules[name] = dftb_results.get_main_molecule()

::

   [02.09|15:36:29] JOB DFTB_H2O STARTED
   [02.09|15:36:29] JOB DFTB_H2O RUNNING
   [02.09|15:36:30] JOB DFTB_H2O FINISHED
   [02.09|15:36:30] JOB DFTB_H2O SUCCESSFUL
   [02.09|15:36:30] JOB DFTB_NH3 STARTED
   [02.09|15:36:30] JOB DFTB_NH3 RUNNING
   [02.09|15:36:31] JOB DFTB_NH3 FINISHED
   [02.09|15:36:31] JOB DFTB_NH3 SUCCESSFUL
   [02.09|15:36:31] JOB DFTB_S2Cl2 STARTED
   [02.09|15:36:31] JOB DFTB_S2Cl2 RUNNING
   ... (PLAMS log lines truncated) ...

.. code:: ipython3

   print(f"Found {len(promising_molecules)} promising molecules with DFTB")

::

   Found 2 promising molecules with DFTB

Optimization and excitations calculation with ADF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each of the molecules identified in the prescreen, run a further calculation with ADF.

.. code:: ipython3

   for name, mol in promising_molecules.items():
       adf_go_job = AMSJob(name="ADF_GO_" + name, molecule=mol, settings=go_sett + adf_sett)
       adf_go_job.run()

       optimized_mol = adf_go_job.results.get_main_molecule()

       adf_exci_job = AMSJob(name="ADF_exci_" + name, molecule=optimized_mol, settings=sp_sett + adf_exci_sett)
       adf_exci_results = adf_exci_job.run()

       if has_good_excitations(adf_exci_results, 2, 4):
           print(f"Molecule {name} has excitation(s) satysfying our criteria!")
           print(optimized_mol)
           exci_energies, oscillator_str = get_excitations(adf_exci_results)
           print("Excitation energy [eV], oscillator strength:")
           for e, o in zip(exci_energies, oscillator_str):
               print(f"{e:8.4f}, {o:8.4f}")

::

   [02.09|15:36:34] JOB ADF_GO_S2Cl2 STARTED
   [02.09|15:36:34] JOB ADF_GO_S2Cl2 RUNNING
   [02.09|15:36:41] JOB ADF_GO_S2Cl2 FINISHED
   [02.09|15:36:41] JOB ADF_GO_S2Cl2 SUCCESSFUL
   [02.09|15:36:41] JOB ADF_exci_S2Cl2 STARTED
   [02.09|15:36:41] JOB ADF_exci_S2Cl2 RUNNING
   [02.09|15:36:47] JOB ADF_exci_S2Cl2 FINISHED
   [02.09|15:36:47] JOB ADF_exci_S2Cl2 SUCCESSFUL
   Molecule S2Cl2 has excitation(s) satysfying our criteria!
     Atoms: 
       1         S      -0.658306      -0.316643       0.909151
       2         S      -0.658306       0.316643      -0.909151
       3        Cl       0.758306       0.752857       2.053019
       4        Cl       0.758306      -0.752857      -2.053019

   Excitation energy [eV], oscillator strength:
     3.4107,   0.0114
     3.5386,   0.0160
     3.5400,   0.0011
     3.9864,   0.1105
     4.3225,   0.0049
     4.3513,   0.2551
     4.7544,   0.0011
     4.9414,   0.0105
     5.3188,   0.0036
     5.3272,   0.0721
   [02.09|15:36:47] JOB ADF_GO_CSCl2 STARTED
   [02.09|15:36:47] JOB ADF_GO_CSCl2 RUNNING
   ... (PLAMS log lines truncated) ...

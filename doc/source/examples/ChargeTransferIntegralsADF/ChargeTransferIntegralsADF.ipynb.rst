Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

   from scm.plams import Settings, ADFFragmentResults, Molecule, log, ADFFragmentJob, add_to_class, init

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/examples/ChargeTransferIntegralsADF/plams_workdir

Helper Functions
~~~~~~~~~~~~~~~~

Add helper results extraction method.

.. code:: ipython3

   @add_to_class(ADFFragmentResults)
   def get_transfer_integrals(self):
       return self.job.full.results.read_rkf_section("TransferIntegrals", file="adf")

Configure Settings
~~~~~~~~~~~~~~~~~~

Set up common settings for all 3 jobs, the specific settings for full system job. Note that in the interest of computational speed we use a minimal basis set. For more quantitatively meaningful results, you should use a larger basis.

.. code:: ipython3

   common = Settings()
   common.input.ams.Task = "SinglePoint"
   common.input.adf.Basis.Type = "SZ"
   common.input.adf.Basis.Core = "None"
   common.input.adf.Symmetry = "NoSym"

.. code:: ipython3

   full = Settings()
   full.input.adf.transferintegrals = True

Load Molecule
~~~~~~~~~~~~~

Load benzene dimer from xyz file and separate it into 2 fragments. Alternatively one could simply load fragments from separate xyz files:

::

   mol1 = Molecule('fragment1.xyz')
   mol2 = Molecule('fragment2.xyz')

.. code:: ipython3

   mol = Molecule("BenzeneDimer.xyz")
   mol.guess_bonds()
   fragments = mol.separate()
   if len(fragments) != 2:
       log("ERROR: Molecule {} was split into {} fragments".format(mol.properties.name, len(fragments)))
       import sys

       sys.exit(1)
   else:
       mol1, mol2 = fragments

Run Job and Get Results
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   job = ADFFragmentJob(name="ADFTI", fragment1=mol1, fragment2=mol2, settings=common, full_settings=full)
   results = job.run()

::

   [10.02|15:07:14] JOB ADFTI STARTED
   [10.02|15:07:14] JOB ADFTI RUNNING
   [10.02|15:07:14] JOB ADFTI/frag1 STARTED
   [10.02|15:07:14] JOB ADFTI/frag1 RUNNING
   [10.02|15:07:20] JOB ADFTI/frag1 FINISHED
   [10.02|15:07:20] JOB ADFTI/frag1 SUCCESSFUL
   [10.02|15:07:20] JOB ADFTI/frag2 STARTED
   [10.02|15:07:20] JOB ADFTI/frag2 RUNNING
   [10.02|15:07:25] JOB ADFTI/frag2 FINISHED
   [10.02|15:07:25] JOB ADFTI/frag2 SUCCESSFUL
   ... (PLAMS log lines truncated) ...

.. code:: ipython3

   # TI is a dictionary with the whole TransferIntegrals section from adf.rkf
   print("== Results ==")
   TI = results.get_transfer_integrals()
   for key, value in sorted(TI.items()):
       print("{:<28}: {:>12.6f}".format(key, value))

::

   == Results ==
   J(charge recombination 12)  :     0.010744
   J(charge recombination 21)  :     0.010744
   J(electron)                 :     0.012050
   J(hole)                     :    -0.034988
   S(charge recombination 12)  :    -0.016229
   S(charge recombination 21)  :    -0.016229
   S(electron)                 :    -0.018755
   S(hole)                     :     0.049852
   V(charge recombination 12)  :     0.009867
   V(charge recombination 21)  :     0.009867
   V(electron)                 :     0.013129
   V(hole)                     :    -0.026790
   Vtot(charge recombination 12):     0.013193
   Vtot(charge recombination 21):     0.013193
   Vtot(electron)              :     0.021464
   Vtot(hole)                  :     0.034178
   e1(electron)                :     0.057279
   e1(hole)                    :    -0.165788
   e2(electron)                :     0.057280
   e2(hole)                    :    -0.165790

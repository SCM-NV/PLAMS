.. _ExcitationsWorkflowExample:

Workflow: filtering molecules based on excitation energies
==========================================================

Identifying systems with certain physical properties out of a large database of molecules is a typical task that can be easily automatized.
In this example we will show how this can be achieved with a simple PLAMS script.


**Introducing the case study**

In our toy study, we want to scan a database of structures, looking for all molecules that absorb light within a certain energy range. For example, between 2 and 4 eV.

A simple approach would be to calculate the excitation energies (and the corresponding oscillator strengths) using ADF's `time-dependent DFT <../../../ADF/Input/Excitation_energies.html>`__ (TD-DFT) for all the molecules in out database of structures.
Since TD-DFT is an expensive method, this procedure can be computationally demanding for large numbers of molecules.

A faster approach would be to pre-screen the large database of molecules by first using a less accurate but more efficient method (e.g. time-dependent `DFTB <../../../DFTB/index.html>`__ (TD-DFTB)) and then run the expensive TD-DFT calculations with ADF only for the systems exhibiting the most promising absorption energies in the faster-but-less-accurate calculations.

The basic workflow may look as follows:

**Step 1**

* Optimize the structure of all molecules in the database with DFTB and calculate excitation energies and oscillator strengths with TD-DFTB
* Select the molecules with electronic excitations of energies between 1 and 6 eV and non-zero oscillator strenght (since TD-DFTB is less accurate than TD-DFT, we opt for a larger energy range in this step). This should significantly reduce the number of molecules to be considered in the following steps.

**Step 2**

* Optimize the structures of the promising molecules with ADF.
* Compute the electronic excitations energies using ADF's TD-DFT and select the molecules with excitations of energies between 2 and 4 eV (and non-zero oscillator strength).

Note: in this example we focus on the **scripting** aspects rather than studying **physical phenomena**. The computational methods used here might not be entirely appropriate to describe the physical properties of interest.

.. include:: ExcitationsWorkflow.common_header.rst
Then download :download:`molecules.tar <../../../../examples/ExcitationsWorkflow/molecules.tar>` and extract it.

.. include:: ExcitationsWorkflow.ipynb.rst
.. include:: ExcitationsWorkflow.common_footer.rst
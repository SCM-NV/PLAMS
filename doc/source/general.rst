.. _General:

General
============


What's new in PLAMS for AMS2023?
--------------------------------------

* The :ref:`AMSCalculator` class for running any AMS engine with ASE (see: :ref:`AMSCalculatorExample`)

* Classes for calculating :ref:`reduction and oxidation potentials  <RedoxExample>` with ADF and optionally COSMO-RS

* The :ref:`ADFCOSMORSCompoundJob <ADFCOSMORSCompound>` class for running jobs equivalent to "Task COSMO-RS Compound" in the AMS GUI. Such a job generates a .coskf file for use with COSMO-RS.

* Classes for running and restarting :ref:`molecular dynamics (MD) jobs with AMS <AMSMDJob>`

* A class for generating and analyzing :ref:`conformers <conformers_interface>`

* :ref:`Quick jobs <Quickjobs>`, like for example the ``preoptimize()`` function let you quickly optimize a Molecule

* :ref:`Packmol interface <PackmolInterface>` for generating liquid and gas mixtures, solid-liquid interfaces, and microsolvation spheres

* :ref:`FileFormatConversionTools` for converting VASP, Gaussian, or Quantum ESPRESSO output to ams.rkf and engine.rkf files that can be opened with the AMS GUI

* :ref:`PlottingTools` for plotting a molecule or ASE Atoms inside a Jupyter notebook

* :ref:`PlottingTools` for plotting the :ref:`electronic band structure <BandStructureExample>`

* Additions to |AMSResults|: get_homo_energies(), get_lumo_energies, get_smallest_homo_lumo_gap()

* Additions to |Molecule|: guess_atomic_charges(), set_density(), get_unique_bonds(), get_unique_angles()

* Many new :ref:`examples`

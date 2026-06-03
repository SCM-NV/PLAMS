.. _PremadeRecipes:

Recipes
-------

These recipe workflows are shipped with PLAMS and can be imported directly.
Each recipe includes details and API documentation.

COSMO-RS Compound
~~~~~~~~~~~~~~~~~

Generate ``.coskf`` files equivalent to AMS GUI task "COSMO-RS Compound".

Use this recipe when you want an easy Python interface for generating COSMO-RS compound files
(``.coskf``) from one or many input structures.

.. currentmodule:: scm.plams.recipes.adfcosmorscompound

.. autoclass:: ADFCOSMORSCompoundJob
   :noindex:
   :no-special-members:

.. autoclass:: ADFCOSMORSCompoundResults
   :noindex:
   :no-special-members:

COSMO-RS Conformers
~~~~~~~~~~~~~~~~~~~

Build customizable conformer workflows that end in ADF/COSMO calculations for COSMO-RS.

This recipe builds a configurable conformer-generation pipeline and produces ``.coskf`` files for
the selected conformers.

.. currentmodule:: scm.plams.recipes.adfcosmorsconformers

.. autoclass:: ADFCOSMORSConfJob
   :noindex:
   :no-special-members:

.. autoclass:: ADFCOSMORSConfFilter
   :noindex:
   :no-special-members:

.. autoclass:: ADFCOSMORSConfResults
   :noindex:
   :no-special-members:

MD Jobs
~~~~~~~

Convenience job classes for MD simulation setup, restart/spawning, and common trajectory analyses.

These classes wrap common AMS MolecularDynamics patterns (NVE/NVT/NPT), restarts and spawned
ensembles, plus trajectory analysis jobs.

.. currentmodule:: scm.plams.recipes.md.amsmdjob

.. autoclass:: AMSMDJob
   :noindex:
   :no-private-members:

.. autoclass:: AMSNVEJob
   :noindex:
   :no-private-members:

.. autoclass:: AMSNVTJob
   :noindex:
   :no-private-members:

.. autoclass:: AMSNPTJob
   :noindex:
   :no-private-members:

.. autoclass:: AMSNPTResults
   :noindex:
   :no-private-members:

.. currentmodule:: scm.plams.recipes.md.nvespawner

.. autoclass:: AMSNVESpawnerJob
   :noindex:
   :no-private-members:

.. currentmodule:: scm.plams.recipes.md.scandensity

.. autoclass:: AMSMDScanDensityJob
   :noindex:
   :no-private-members:

.. currentmodule:: scm.plams.recipes.md.trajectoryanalysis

.. autoclass:: AMSMSDJob
   :noindex:
   :no-private-members:

.. autoclass:: AMSRDFJob
   :noindex:
   :no-private-members:

.. autoclass:: AMSVACFJob
   :noindex:
   :no-private-members:

.. autoclass:: AMSViscosityFromBinLogJob
   :noindex:
   :no-private-members:

ADF Fragment
~~~~~~~~~~~~

Run fragment-based ADF analysis as a dedicated multi-job workflow.

This recipe orchestrates two fragment calculations and one full-system ADF calculation and exposes
energy-decomposition and related properties through a dedicated results class.

.. currentmodule:: scm.plams.recipes.adffragment

.. autoclass:: ADFFragmentJob
   :noindex:
   :no-special-members:

.. autoclass:: ADFFragmentResults
   :noindex:
   :no-special-members:

BAND Fragment
~~~~~~~~~~~~~

Run BAND energy decomposition workflows for periodic systems with optional NOCV analysis.

This recipe extends the fragment approach to periodic BAND systems and supports NOCV workflows via
``NOCVBandFragmentJob``.

.. currentmodule:: scm.plams.recipes.bandfragment

.. autoclass:: BANDFragmentJob
   :noindex:
   :no-special-members:

.. autoclass:: BANDFragmentResults
   :noindex:
   :no-special-members:

.. autoclass:: NOCVBandFragmentJob
   :noindex:
   :no-special-members:

Reorganization Energy
~~~~~~~~~~~~~~~~~~~~~

Compute Marcus reorganization energies using four coordinated AMS calculations.

The workflow evaluates energies at optimized geometries of two electronic states and combines them
into the reorganization energy used in Marcus theory.

.. currentmodule:: scm.plams.recipes.reorganization_energy

.. autoclass:: ReorganizationEnergyJob
   :noindex:
   :no-special-members:

.. autoclass:: ReorganizationEnergyResults
   :noindex:
   :no-special-members:

ADFNBO
~~~~~~

Extend ADF runs with NBO6/gennbo execution through an AMS job wrapper.

``ADFNBOJob`` extends the regular AMS/ADF runscript to execute ``adfnbo`` and ``gennbo6`` while
ensuring required ADF input is present.

.. currentmodule:: scm.plams.recipes.adfnbo

.. autoclass:: ADFNBOJob
   :noindex:
   :no-special-members:


pyAHFCDOS
~~~~~~~~~

Compute vibronic density of states from two vibronic spectra at AH-FC level.

The ``FCFDOS`` class combines absorption and emission vibronic data to build a DOS profile with
configurable broadening.

.. currentmodule:: scm.plams.recipes.fcf_dos

.. autoclass:: FCFDOS
   :noindex:
   :no-special-members:

ADFVibronicDOS
~~~~~~~~~~~~~~

Practical ADF workflow that uses AH-FC vibronic spectra to build a DOS profile.

This workflow runs the vibronic calculations and then applies ``FCFDOS`` to produce and analyze the
vibronic density of states.

.. currentmodule:: scm.plams.recipes.fcf_dos

.. autoclass:: FCFDOS
   :noindex:
   :no-special-members:

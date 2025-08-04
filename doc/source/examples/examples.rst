.. _examples:

Examples
========

The PLAMS examples use the Amsterdam Modeling Suite. Get a license or free trial from www.scm.com.

.. dropdown:: Note about $AMSBIN for external Python environments

   If you run PLAMS with ``amspython`` (recommended), then you most likely
   have already set the $AMSBIN environment correctly and can ignore the
   rest of this message.

   If you do not use ``amspython`` but a custom Python environment, make sure
   that the ``AMSBIN`` environment variable is properly set. 

   You can test this by typing ``$AMSBIN/amspython -h`` in a terminal: this
   should print the help message. 

   If this is not the case, see the `AMS Installation
   documentation <../../Installation/index.html>`__ for details.

Getting Started
----------------

.. toctree::
   :maxdepth: 1
 
   WaterOptimization/WaterOptimization
   AMSSettingsSystem/AMSSettingsSystem
   AMSSettingsInput/AMSSettingsInput
   He2DissociationCurve/He2DissociationCurve
   ManyJobsInParallel/ManyJobsInParallel
   Logging/Logging
   JobAnalysis/JobAnalysis

Molecule analysis
---------------------

.. toctree::
   :maxdepth: 1

   MoleculeFormats/MoleculeFormats
   MoleculesFromRKFTrajectory
   MoleculesTable
   MoleculeSubstitution/MoleculeSubstitution
   ConvertToAMSRKFTrajectory
   PlotCorrelation/PlotCorrelation
   MapMoleculesAndConvertToDCD
   HydrogenBondsFromMD
   PlotReaction2D/PlotReaction2D

MD trajectory analysis
----------------------

.. toctree::
    :maxdepth: 1

    IonicConductivityFromMD.rst
    ViscosityGreenKubo/ViscosityGreenKubo
    IRSpectrumFromMD/IRSpectrumFromMD.rst
    IRSpectrumFromMDH2ODimer.rst

Benchmarks
-----------------

.. toctree::
   :maxdepth: 1

   BasisSetBenchmark/BasisSetBenchmark
   ReactionEnergyBenchmark/ReactionEnergyBenchmark

Workflows
------------------

.. toctree::
   :maxdepth: 1

   RedoxPotential
   ExcitationsWorkflow/ExcitationsWorkflow
   AMSTSWorkflow/AMSTSWorkflow
   ChargeTransferIntegralsADF/ChargeTransferIntegralsADF
   TuningRangeSeparation/TuningRangeSeparation
   ConformersGeneration/ConformersGeneration
   ConformersMultipleMolecules/ConformersMultipleMolecules
   ReactionsDiscovery/ReactionsDiscovery

COSMO-RS and property prediction
-----------------------------------

For more examples, see the `COSMO-RS documentation <../../COSMO-RS/Python_Examples.html>`__.

.. toctree::
   :maxdepth: 1

   PropertyPrediction/PropertyPrediction
   ams_crs

Packmol and AMS-ASE interfaces
-------------------------------

.. toctree::
   :maxdepth: 1

   PackMolExample/PackMolExample
   CustomASECalculator
   AMSCalculator/ASECalculator
   ChargedAMSCalculator/ChargedAMSCalculator
   i-PI-AMS
   SellaTransitionStateSearch

.. N.B. the ASE calculator example is linked to via an external site, so the path should not be updated (or a redirect needs to be added)

ParAMS, Simple Active Learning, and pyZacros
---------------------------------------------------

See the respective documentation pages:

* `ParAMS <../../params/index.html>`__ 
* `Simple Active Learning <../../Workflows/SimpleActiveLearning/SimpleActiveLearning.html>`__
* `pyZacros <../../pyzacros/index.html>`__ 

Other AMS calculations
------------------------

.. toctree::
   :maxdepth: 1

   BAND_NiO_HubbardU
   BandStructure/BandStructure
   AMSPlumedMD/AMSPlumedMD
   QE_AMS_AFM_HubbardU
   BasicMDAnalysis/BasicMDAnalysis
   UseLowestEnergy/UseLowestEnergy
   M3GNet/M3GNet
   ConstrainedGOAMSWorker/ConstrainedGOAMSWorker

Pymatgen
----------------

.. toctree::
    :maxdepth: 1
    
    XRD/XRD

.. _recipes:

Pre-made recipes
----------------

The examples presented in here are simple job types built using basic PLAMS elements.
They are shipped with PLAMS in the ``recipes`` subpackage and can be directly used in your scripts.
In other words, the code presented there is already included in PLAMS and (unlike examples from two other sections) does not need to be copied to your script.
The source code of ``recipes`` modules is presented here to demonstrate how easy it is to build on top of existing PLAMS elements and create your own fully customized job types.

.. toctree::
    :maxdepth: 2

    COSMORSCompound/COSMORSCompound
    COSMORSConformers/COSMORSConformers
    MDJobs
    ADFFrag/ADFFrag
    BandFrag/BandFrag
    ReorganizationEnergy/ReorganizationEnergy
    ADFNBO/ADFNBO
    NumGrad/NumGrad
    NumHess/NumHess
    pyAHFCDOS
    ADFVibronicDOS/ADFVibronicDOS.rst
    ReuseForceFieldParams/ReuseForceFieldParams


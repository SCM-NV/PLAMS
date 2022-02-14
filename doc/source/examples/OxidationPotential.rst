Oxidation Potential
===================

A recent paper by Belic, J. et al. (Phys. Chem. Chem. Phys., **24 (1)**, 197–210 (2022)) has investigated ways to calculate oxidation potentials using AMS and has shown that they are very accurate to reference. The paper describes four ways to calculate the oxidation potential which are described in the tabs below. In general, :math:`g_p^i` denotes the optimised geometry of the molecule in phase :math:`p` (solvent (sol) or gaseous (gas)) and in state :math:`i` (oxidised (+) or neutral (0)), :math:`G_P^I(g_p^i)` denotes the Gibbs free energy of the geometry :math:`g_p^i` calculated in phase :math:`P` and state :math:`I`, similarly :math:`E_P(g_p^i)` denotes the bond energy of the geometry :math:`g_p^i` calculated in phase :math:`P`. The electron Gibbs free energy :math:`G_{gas}(e^-)=0.0375 \text{eV}` is also used in these calculations.

.. tabs::
 
  .. tab:: Direct (DC)

    Direct method of calculating oxidation potential using the COSMO solvation model.
    
    .. math::
    	\Delta G^{DC}_{COSMO} = G_{sol}^+(g_{sol}^+) + G_{gas}(e^-) - G_{sol}^0(g_{sol}^0)


  .. tab:: Thermodynamic Cycle (TC)

  	Thermodynamic cycle using either the COSMO or the COSMO-RS solvation model. 

  	.. math::
  		\Delta G^{TC}_{COSMO/COSMO-RS} = G_{sol}^+(g_{gas}^+) + G_{gas}(e^-) - G_{sol}^0(g_{gas}^0)

  	where we have 

  	.. math::
  		G_{sol}^i = G_{gas}^i(g_{gas}^i) + \Delta G_{sol}^i(g_{gas}^i) + [E_{gas}(g_{sol}^i) - E_{gas}(g_{gas}^i)].


  .. tab:: Screening

  	This method also implements the thermodynamic cycle, but uses lower theory methods (DFTB) to greatly increase speed of calculation. This method uses COSMO-RS as the solvation method.

  	.. math::
  		\Delta G_{COSMO-RS}^{screening} = G_{sol, CRS}^+(g_{gas}^+) + G_{gas}(e^-) - G_{sol, CRS}^0(g_{gas}^0)

  	where we have

  	.. math::
  		G_{sol, CRS}^i = E_{sol}^i(g_{gas}^i) + \Delta G_{sol, CRS}^i(g_{gas}^i])

This recipe implements these methods for your use. The class ``OxidationPotentialCalculator`` implements all the required methods for calculating the oxidation potential for any molecule using one of the methods described above.

To perform a calculation, create an ``OxidationPotentialCalculator`` object and then call it. You must pass an ``scm.plams.Molecule`` object, optionally you can specify the method, which must be one of ``['DC', 'TS-COSMO', 'TC-COSMO-RS', 'screening']`` (by default uses ``'screening'``). You can also specify the job name and job directory for PLAMS. Once the ``OxidationPotentialCalculator`` object has been created you can change the default settings by accessing the ``scm.plams.Settings`` attributes set in the ``set_default_settings`` method. By default the calculator uses B3LYP-D3(BJ)/TZ2P level of theory for the DFT calculations and GFN1-xTB for the DFTB calculations. The default COSMO solvent is dichloromethane, to change the COSMO-RS solvent you will need to perform a COSMO-RS job before the oxidation potential calculation and provide the ``.coskf`` file (`COSMO-RS tutorials <../../Tutorials/COSMO-RS/index.html>`__).

.. literalinclude:: ../../../recipes/ox_potential.py
	:language: python

**Example usage:** (:download:`Download input files <../../../examples/OxidationPotentialFiles.zip>`)

.. literalinclude:: ../../../examples/OxidationPotential.py
	:language: python

**Example output:**

.. literalinclude:: ./OxidationPotential.out

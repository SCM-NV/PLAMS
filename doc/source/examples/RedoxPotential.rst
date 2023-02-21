.. _RedoxExample:

Reduction and oxidation potentials
====================================

**Note**: This example requires AMS2023 or later.

Definitions and introduction
------------------------------

**Reduction potential**: A + e\ :sup:`–` → A\ :sup:`–` ; ΔG⁰ = –nFE⁰

**Oxidation potential**: A → A\ :sup:`+` + e\ :sup:`–` ; ΔG⁰ = nFE⁰

There are three PLAMS recipes for calculating one-electron reduction or oxidation potentials in implicit solvent:

* ``AMSRedoxDirectJob``: The best method. Geometry optimizations (and optionally frequencies) are calculated for both neutral and reduced/oxidized species in implicit solvent. The solvent must be supported by ADF. Requires an ADF license.

* ``AMSRedoxThermodynamicCycleJob``: Only useful if you include the vibrations (frequencies) and the molecule is large (in which case it is faster but less accurate than AMSRedoxDirectJob). The frequencies are only calculated for the gasphase molecule. A thermodynamic cycle gives the reduction or oxidation potential. The solvent must be supported by ADF. Requires an ADF license.

* ``AMSRedoxScreeningJob``: The fastest (and least accurate) method. Geometry optimizations are performed at the GFN1-xTB level of theory. The solvation free energy is evaluated by COSMO-RS. Vibrational effects are always implicitly accounted for. A ``.coskf`` file for the solvent is required - this can either be obtained from the ADFCRS-2018 database or generated with an ``ADFCOSMORSCompoundJob``. Requires DFTB, ADF, and COSMO-RS licenses.

The ``AMSRedoxScreeningJob`` workflow was developed by `Belic et al., Phys. Chem. Chem. Phys. 24, 197–210 (2022) <https://pubs.rsc.org/en/content/articlehtml/2021/cp/d1cp04218a>`__ for calculating oxidation potentials. The paper also describes the other workflows in detail. 

.. note::

    The ``AMSRedoxScreeningJob`` uses the recommended ADF settings for generating .coskf files. Belic et al. used a different density functional. Using the class will give slightly different results compared to Belic et al.

The **free energy of the electron** is set to be –0.0375 eV (see Belic et al.).

In this example the above three workflows are used to evaluate the **reduction potential of benzoquinone in water**.
The experimental value is E⁰ = +0.10 V relative to SHE (**standard hydrogen electrode**).

`Ho et al. <https://doi.org/10.1201/b19122>`__  (`preprint pdf <https://comp.chem.umn.edu/truhlar/docs/C86preprint.pdf>`__)
evaluated E⁰ = -0.40 V or -0.28 V for benzoquinone relative to SHE  using different computational methods.

The reduction and oxidation potentials calculated by the PLAMS classes are given on an **absolute scale**. On this scale,
the SHE is at +4.42 or +4.28 V (see Ho et al.). To get potentials relative to SHE, we therefore subtract 4.42 V from the calculated values.

.. note::

    Energy differences in eV correspond directly to potentials in V.

.. tip::

    If you compare many calculated reduction potentials, use one of them as the reference state. See for example `Huyhn et al., J. Am. Chem. Soc. 2016, 138, 49, 15903–15910 <https://doi.org/10.1021/jacs.6b05797>`__ in the case of substituted quinones.

.. important::

    These classes can only be used for 
    
    * molecules that do not undergo big conformational changes (or dissociation) upon oxidation/reduction.

    * one-electron reduction/oxidation. For multiple electrons, run the same script with different initial charges for the molecule

    * molecules with 0 or 1 unpaired electrons

Redox potential results
-------------------------

**Results**: The below script outputs the following table:

.. code-block:: none

    The experimental reduction potential of benzoquinone is +0.10 V vs. SHE
    Jobname                  Eox(vib,rel-to-SHE)[V]   Ered(vib,rel-to-SHE)[V]  Eox(rel-to-SHE)[V]       Ered(rel-to-SHE)[V]
    quick_screening          2.95                     0.44                     2.95                     0.44
    direct_best_method       2.62                     -0.08                    2.72                     -0.09
    thermodynamic_cycle      2.64                     -0.07                    2.71                     -0.10

The first two columns give the oxidation and reduction potentials incorporating
vibrational effects in the calculation. The last two columns ignore the
vibrational effects. Here, we are only interested in the two **Ered** columns
(reduction potentials). We see:

* Compared to the experimental E⁰ = +0.10 V, the screening method gives a higher value (+0.44 V) and the more accurate methods a lower value (-0.08 V or -0.07 V).

* For this molecule, vibrational effects on the reduction potential are very small (e.g. the vibrations cause E⁰ to go from -0.09 V to -0.08 V for the direct method). Because calculating the vibrations is computationally expensive, they could have been turned off by setting ``vibrations = False`` in the below script.

* The oxidation potentials (Eox) are given for demonstration purposes only. Turn them off by setting ``oxidation = False`` in the below script.

Code example
--------------------

**Download** :download:`ReductionOxidationPotentials.py <../../../examples/ReductionOxidationPotentials.py>`

.. literalinclude:: ../../../examples/ReductionOxidationPotentials.py

The PLAMS recipes/classes
---------------------------

This section contains some details about the implementation.

The oxidation and reduction potentials are calculated as follows:

.. literalinclude:: ../../../recipes/redox.py
    :pyobject: AMSRedoxDirectResults 

.. literalinclude:: ../../../recipes/redox.py
    :pyobject: AMSRedoxThermodynamicCycleResults  

where 

* ``go_0_vacuum`` is a **g**\ eometry **o**\ ptimization for the neutral (**0**) molecule in vacuum, 
* ``go_0_vaccum_sp_solvated`` is a **s**\ ingle **p**\ oint in implicit solvent on the previous structure,
* etc.

.. literalinclude:: ../../../recipes/redox.py
    :pyobject: AMSRedoxScreeningResults

The complete class definitions:

.. literalinclude:: ../../../recipes/redox.py
  :language: python

Details about individual calculations
----------------------------------------

The paper by Belic et al. describes four ways to calculate the oxidation potential which are described in the tabs below, however we have extended them to allow for calculation of the reduction potential. In general, :math:`g_p^i` denotes the optimised geometry of the molecule in phase :math:`p` (solvent [sol] or gaseous [gas]) and in state :math:`i` (oxidised [+]/reduced [-] or neutral [0]), :math:`G_P^I(g_p^i)` denotes the Gibbs free energy of the geometry :math:`g_p^i` calculated in phase :math:`P` and state :math:`I`, similarly :math:`E_P(g_p^i)` denotes the bond energy of the geometry :math:`g_p^i` calculated in phase :math:`P`. The electron Gibbs free energy :math:`G_{gas}(e^-)=0.0375 \text{eV}` is also used in these calculations.

.. tabs::

  .. tab:: Direct (DC)

    Direct method of calculating oxidation potential using the COSMO solvation model.
    
    .. math::
      \Delta G^{DC}_{COSMO} = G_{sol}^\pm(g_{sol}^\pm) \pm G_{gas}(e^-) - G_{sol}^0(g_{sol}^0)

    The following steps are calculated

    .. list-table::
      :header-rows: 1

      * - Step 
        - Task
        - Structure
        - Method
        - I
        - P
        - Frequencies

      * - 1
        - GO
        - 
        - DFT
        - 0
        - sol
        - Yes

      * - 2
        - GO
        - 
        - DFT
        - :math:`\pm`
        - sol
        - Yes



  .. tab:: Thermodynamic Cycle (TC)

    Thermodynamic cycle using either the COSMO or the COSMO-RS solvation model. 

    .. math::
      \Delta G^{TC}_{COSMO/COSMO-RS} = G_{sol}^\pm(g_{gas}^\pm) \pm G_{gas}(e^-) - G_{sol}^0(g_{gas}^0)

    where we have 

    .. math::
      G_{sol}^i = G_{gas}^i(g_{gas}^i) + \Delta G_{sol}^i(g_{gas}^i) + [E_{gas}(g_{sol}^i) - E_{gas}(g_{gas}^i)].

    The following steps are calculated

    .. list-table::
      :header-rows: 1

      * - Step 
        - Task
        - Structure
        - Method
        - I
        - P
        - Frequencies

      * - 1
        - GO
        - 
        - DFT
        - 0
        - gas
        - Yes

      * - 2
        - SP
        - :math:`g^0_{gas}`, step 1
        - DFT
        - 0
        - sol
        - No

      * - 3
        - GO
        - 
        - DFT
        - 0
        - sol
        - No

      * - 4
        - SP
        - :math:`g^0_{sol}`, step 3
        - DFT
        - 0
        - gas
        - No

      * - 5
        - GO
        - 
        - DFT
        - :math:`\pm`
        - gas
        - Yes

      * - 6
        - SP
        - :math:`g^\pm_{gas}`, step 5
        - DFT
        - :math:`\pm`
        - sol
        - No

      * - 7
        - GO
        - 
        - DFT
        - :math:`\pm`
        - sol
        - No

      * - 8
        - SP
        - :math:`g^\pm_{sol}`, step 7
        - DFT
        - :math:`\pm`
        - gas
        - No


  .. tab:: Screening

    This method also implements the thermodynamic cycle, but uses lower theory methods (DFTB) to greatly increase speed of calculation. This method uses COSMO-RS as the solvation method.

    .. math::
      \Delta G_{COSMO-RS}^{screening} = G_{sol, CRS}^\pm(g_{gas}^\pm) \pm G_{gas}(e^-) - G_{sol, CRS}^0(g_{gas}^0)

    where we have

    .. math::
      G_{sol, CRS}^i = E_{sol}^i(g_{gas}^i) + \Delta G_{sol, CRS}^i(g_{gas}^i])

    .. list-table::
      :header-rows: 1

      * - Step 
        - Task
        - Structure
        - Method
        - I
        - P
        - Frequencies

      * - 1
        - GO
        - 
        - DFTB
        - 0
        - gas
        - No

      * - 2
        - SP
        - :math:`g^0_{gas}`, step 1
        - DFT
        - 0
        - sol
        - No

      * - 3
        - GO
        - 
        - DFTB
        - :math:`\pm`
        - gas
        - No

      * - 4
        - SP
        - :math:`g^\pm_{gas}`, step 3
        - DFT
        - :math:`\pm`
        - sol
        - No


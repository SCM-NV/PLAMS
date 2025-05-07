.. _ams_crs_workflow:

ADF and COSMO-RS workflow
-------------------------

**Note**: This example requires AMS2023 or later.

This example uses ADF to generate the .coskf file for benzene. You can also
modify it to instead use the Benzene.coskf from the ADFCRS-2018 database.
Note that you first need to install the ADFCRS-2018 database.

The script then plots the sigma profile. This is not a necessary step but is done for illustration purposes only.

Then a solubility calculation is performed for benzene in water between 0 and
10 degrees C. The melting point and enthalpy of fusion can either be estimated
using the property prediction tool, or the experimental numbers can be given
(recommended).

**Example usage**: AMS2023+ required (:download:`ams_crs.py <../../../examples/ams_crs.py>`)

This script is built on the PLAMS. For more detailed useage of COSMO-RS with PLAMS, refer to the `documentation <../../COSMO-RS/PLAMS_COSMO-RS_scripting.html#python-scripting-for-cosmo-rs-with-plams>`__.

.. raw:: html

   <details>
   <summary style="color:#008000;cursor:pointer">[show/hide code]</summary>
.. literalinclude:: ../../../examples/ams_crs.py
   :language: python

.. raw:: html

    </details>


**pyCRS example usage**: AMS2024+ required (:download:`ams_pyCRS.py <../../../examples/ams_pyCRS.py>`)

This script utilizes pyCRS, a convenient python wrapper built on PLAMS and SQL database for various thermodynamic calculations. For more detatils, refer to the `pyCRS Overview <../../COSMO-RS/pyCRS_Overview.html#overview>`__.

.. raw:: html

   <details>
   <summary style="color:#008000;cursor:pointer">[show/hide code]</summary>
.. literalinclude:: ../../../examples/ams_pyCRS.py
   :language: python

.. raw:: html

    </details>


**Results**

Both scripts above will generate the following results.

.. image:: /_static/ams_crs_benzene_sigmaprofile.png

.. parsed-literal::

   273.15 1.6797
   274.15 1.7258
   275.15 1.7731
   276.15 1.8215
   277.15 1.8711
   278.15 1.9220
   279.15 1.9603
   280.15 1.9823
   281.15 2.0046
   282.15 2.0273
   283.15 2.0503

.. image:: /_static/ams_crs_benzene_solubility.png

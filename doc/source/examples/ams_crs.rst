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

**Example usage**: (:download:`ams_crs.py <../../../examples/ams_crs.py>`)

.. literalinclude:: ../../../examples/ams_crs.py



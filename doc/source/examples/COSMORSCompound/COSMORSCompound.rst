.. _ADFCOSMORSCompound:

ADF: Task COSMO-RS Compound
===========================

The |ADFCOSMORSCompoundJob| class generates results identical to the "Task COSMO-RS Compound" in the AMS ADF graphical user interface.  This python interface allows users to easily generate the .coskf files for one or many structures.  A possible usage is given in :ref:`ams_crs_workflow`.

.. include:: COSMORSCompound.common_header.rst
.. include:: cosmors_compound.ipynb.rst
.. include:: COSMORSCompound.common_footer.rst

Source code for ``ADFCOSMORSCompound``
--------------------------------------

.. dropdown::

    .. literalinclude:: ../../../../src/scm/plams/recipes/adfcosmorscompound.py

Brief API Documentation
-----------------------

.. automodule:: scm.plams.recipes.adfcosmorscompound
    :no-special-members:
    :exclude-members: _result_type, __init__, new_children, postrun, _get_radii, adf_settings
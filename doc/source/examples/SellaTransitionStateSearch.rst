.. _SellaExample:

Sella transition state search with AMS
================================================

**Note**: This example requires AMS2023 or later.

This example shows how to couple the `Sella <https://doi.org/10.1021/acs.jctc.2c00395>`__ code to AMS using the 
AMS :ref:`ASE calculator <ASECalculatorExample>`.

Sella implements special algorithms and coordinate systems for performing transition state searches.

AMS also internally implements several methods. The example compares

* Sella with AMS/ADF providing the energies and forces
* AMS-only, where the initial hessian is calculated with DFTB and the TS search is done with ADF

For more information about Sella, refer to the Sella documentation and examples.

.. important::

    Sella is not included with AMS and not supported by SCM. To install it into the AMS python environment, run

    .. code-block:: none
        
        $AMSBIN/amspython -m pip install sella

**Example usage:** (:download:`Download SellaTransitionStateSearch.py <../../../examples/SellaTransitionStateSearch.py>`)

.. literalinclude:: ../../../examples/SellaTransitionStateSearch.py
    :language: python


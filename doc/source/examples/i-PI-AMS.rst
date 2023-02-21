.. _iPIExample:

i-PI path integral MD with AMS
================================================

**Note**: This example requires AMS2023 or later.

This example shows how to couple the `i-PI <http://ipi-code.org>`__ code to AMS using the 
AMS :ref:`ASE calculator <ASECalculatorExample>`.

i-PI can be used to run for example path integral molecular dynamics. The
example shows how to run thermostatted ring polymer molecular dynamics for a
water molecule together with the UFF force field.

For more information about i-PI, refer to the i-PI documentation and examples.

The below example only runs on Linux as it uses a unix socket.

.. important::

    i-PI is not included with AMS and is not supported by SCM. 

**Example usage:** (:download:`Download run-ase.py <../../../examples/i-PI-AMS/run-ase.py>` and auxiliary files :download:`input.xml <../../../examples/i-PI-AMS/input.xml>`, :download:`firstframe.xyz <../../../examples/i-PI-AMS/firstframe.xyz>`, :download:`run-server.sh <../../../examples/i-PI-AMS/run-server.sh>`)

.. literalinclude:: ../../../examples/i-PI-AMS/run-ase.py
    :language: python


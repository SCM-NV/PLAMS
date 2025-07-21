.. _ReorganizationEnergyExample:

Reorganization Energy
=====================

One of the ingredients for computing hopping rates in Marcus theory is the  reorganization energy :math:`\lambda`, defined as

.. math::

   \lambda = & (E^\text{state B}_\text{at optimal geometry of state A} - E^\text{state A}_\text{at optimal geometry of state A}) + 

   & (E^\text{state A}_\text{at optimal geometry of state B} - E^\text{state B}_\text{at optimal geometry of state B})

where states A and B are two electronic configurations, e.g. neutral and anion (**see the example usage below**).

In this recipe we build a job class ``ReorganizationEnergyJob`` by extending |MultiJob|. Our job will perform four |AMSJob| calculations: two geometry optimizations for states A anb B, followed by two single point calculations (state A at the optimal geometry of state B and state B at the optimal geometry of state A).

In ``ReorganizationEnergyResults``, the reorganization energy is computed by fetching and combining the results from the children jobs.

.. literalinclude:: ../../../../src/scm/recipes/reorganization_energy.py


.. include:: ReorganizationEnergy.common_header.rst
.. include:: ReorganizationEnergy.ipynb.rst
.. include:: ReorganizationEnergy.common_footer.rst



.. _TuningRangeSeparationExample:

Tuning the range separation
===========================

In this example we optimize the value of *gamma* parameter for long-range corrected XC functional (in our case: LCY-PBE) in ADF.
Long-range corrected XC functionals can be used in ADF with XCfun (see `ADF manual  <../../../ADF/Input/Density_Functional.html#range-separated-hybrids>`_).

The optimal range separation parameter *gamma* yields the HOMO energy equal to the ionization potential (IP).
Given a molecular system, we simultaneously minimize the difference between HOMO and IP for that system (N) and its anion (A) (system with one more electron). We define the J function as:

.. math::

   J = \sqrt{N^2+A^2}

and find the value of *gamma* (within a certain range) which minimizes J. See also `this article by Kronik and coworkers <https://doi.org/10.1063/1.4807325>`__.

We first define a new job type ``GammaJob`` by extending |MultiJob|.
The goal of ``GammaJob`` is to calculate the J function for one fixed value of *gamma*
To do that we need to perform 3 different single point calculations: 1 for the given system (let's call it S0), 1 for the system with one more electron (S-) and 1 for the system with one less electron (S+).
S+ calculation is needed to find the ionization potential of S0.

The constructor (``__init__``) of ``GammaJob`` accepts several new arguments and simply stores them.
These new arguments define: the value of *gamma*, the |Molecule| together with its initial charge, and the values of spin for S-, S0 and S+ (as a tuple of length 3).
Then the |prerun| method is used to create three children jobs with different values of total charge and spin multiplicity.
A dedicated |Results| subclass features a simple method for extracting the value of J based on results on three children jobs.

We can then treat our newly defined ``GammaJob`` as a blackbox with simple interface: input *gamma* -> run -> extract J.
The next step is to create multiple instances of ``GammaJob`` for a range of different *gammas*.
That task can be conveniently wrapped in a simple function ``gamma_scan``.

.. include:: TuningRangeSeparation.common_header.rst
.. include:: TuningRangeSeparation.ipynb.rst
.. include:: TuningRangeSeparation.common_footer.rst
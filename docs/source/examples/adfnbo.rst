NBO with ADF
------------

NBO 6.0 is a tool for Natural Bond order analysis that uses the results of an ADF calculation.
More information about NBO can be found in the corresponding section of the `ADF manual  <../../ADF/Input/Advanced_analysis.html#adfnbo-gennbo-nbo-analysis>`_.

NBO analysis is performed based on a prior ADF calculation (with some special keywords) by using two separate binary executables: ``adfnbo`` and ``gennbo6``.
In this case no special job type is created for these binaries.
Instead of that we extend the |AMSJob| class in such a way that calls of ``adfnbo`` and ``gennbo6`` are appended to the usual ADF runscript.
We also make sure that all the required ADF keywords are present in the initial ADF input file.
Input keywords for ``adfnbo`` are taken from ``myjob.settings.adfnbo``.
All this work happens in |prerun|.
No specialized |Results| subclass is defined for ``ADFNBOJob``.


The source code of the whole module:

.. literalinclude:: ../../../recipes/adfnbo.py

An example usage:
(:download:`methane.xyz <../../../examples/ADFNBO/methane.xyz>`,
:download:`adfnbo_test.py <../../../examples/ADFNBO/adfnbo.py>`)


.. literalinclude:: ../../../examples/ADFNBO/methane.xyz
.. literalinclude:: ../../../examples/ADFNBO/adfnbo.py

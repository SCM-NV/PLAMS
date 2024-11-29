Reaction energies with many different engines
==============================================

.. literalinclude:: ../../../examples/ReactionEnergyBenchmark.py
   :language: python

.. note::
   The Amsterdam Modeling Suite requires the installation of additional Python
   packages to run the machine learning potential backends.
   You can use the following command to install all available machine learning potentials:

   .. code-block:: sh

      $ "$AMSBIN"/amspackages install mlpotentials

.. note::
    To execute this PLAMS script:
    
    * :download:`Download ReactionEnergyBenchmark.py <../../../examples/ReactionEnergyBenchmark.py>` and :download:`download and unzip molecules_ReactionEnergyBenchmark.zip <../../../examples/molecules_ReactionEnergyBenchmark.zip>`
    * ``$AMSBIN/amspython ReactionEnergyBenchmark.py``

    This PLAMS script might take a while to run – around 4 hours on a modern laptop. Feel free to grab a coffee or take a break while it executes!


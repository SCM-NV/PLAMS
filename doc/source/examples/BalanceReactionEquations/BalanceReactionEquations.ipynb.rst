Worked Example
--------------

Test 1. Balance simple reactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we present eleven different sets of reactants and products, not necessarily based on real reactions. The molecules can be passed as formulas, or as PLAMS Molecule objects.

.. code:: ipython3

   from scm.plams import ReactionEquation
   from scm.plams import from_smiles

Aspirin in water
~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["C9H8O4", from_smiles("O")]
   products = ["C2H4O2", "C7H6O3"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   1 C9H8O4 + 1 H2O => 1 C2H4O2 + 1 C7H6O3 | Charge = 0.00
   Coefficients:  [1 1 1 1]

By default, a native method is used to compute the nullspace. Optionally, this can be done with sympy. The result should be the same

.. code:: ipython3

   reaction.method = 'sympy'

.. code:: ipython3

   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)

::

   1 C9H8O4 + 1 H2O => 1 C2H4O2 + 1 C7H6O3 | Charge = 0.00

Aspirin in water with additional reactants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["C9H8O4", "H2O", "HO", "HO"]
   products = ["C2H4O2", "C7H6O3", "H2O2"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   1 C9H8O4 + 1 H2O => 1 C2H4O2 + 1 C7H6O3 | Charge = 0.00
   Coefficients:  [1 1 0 0 1 1 0]

Aspirin in water with different products
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["C9H8O4", "H2O"]
   products = ["CH2O", "C7H6O3"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   1 C9H8O4 + 1 H2O => 2 CH2O + 1 C7H6O3 | Charge = 0.00
   Coefficients:  [1 1 2 1]

Oxydation of ethane
~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["C2H6", "O2"]
   products = ["CO2", "H2O"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   2 C2H6 + 7 O2 => 4 CO2 + 6 H2O | Charge = 0.00
   Coefficients:  [2 7 4 6]

Exampe that is not easily balanced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example SMILES string for the reactant is O=CC(=O)c1cccc(c1)C(=O)O

.. code:: ipython3

   reactants = ["C9O4H6", "OH"]
   products = ["C2O2H3", "C7O3H2"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   9 C9O4H6 + 25 OH => 23 C2O2H3 + 5 C7O3H2 | Charge = 0.00
   Coefficients:  [ 9 25 23  5]

Carbonmonoxide with carbondioxide and hydrogen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["CO", "CO2", "H2"]
   products = ["CH4", "H2O"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   1 CO + 3 H2 => 1 CH4 + 1 H2O | Charge = 0.00
   Coefficients:  [1 0 3 1 1]

Example of a reaction that cannot be balanced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["FeS2", "HNO3"]
   products = ["Fe2S3O12", "NO", "H2SO4"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)
   print (reaction.message)

::

   Inconsistent system: No solution
   Coefficients:  None
   Empty nullspace

Potassiumnitrate and methane
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["KNO3", "C"]
   products = ["K2CO3", "CO", "N2"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   2 KNO3 + 4 C => 1 K2CO3 + 3 CO + 1 N2 | Charge = 0.00
   Coefficients:  [2 4 1 3 1]

Pyrite and nitric acid
~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["FeS2", "HNO3"]
   products = ["Fe2S4O12", "N2H2"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   2 FeS2 + 4 HNO3 => 1 Fe2S4O12 + 2 N2H2 | Charge = 0.00
   Coefficients:  [2 4 1 2]

Reacting iron ligands
~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["FeS2O6N5H3"]
   products = ["Fe2S4O12N10H6"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   2 FeS2O6N5H3 => 1 Fe2S4O12N10H6 | Charge = 0.00
   Coefficients:  [2 1]

Example of a reaction with non-matching elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["KNO3", "C"]
   products = ["Fe"]
   reaction = ReactionEquation(reactants, products)
   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)
   print (reaction.message)

::

   Inconsistent system: No solution
   Coefficients:  None
   Empty nullspace

Test 2. Balance charged reactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Aspirin and hydroxide ion
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["C9H8O4", "OH"]
   products = ["C2H3O2", "C7H6O3"]
   reaction = ReactionEquation(reactants, products)
   rcharges = [0, -1]
   pcharges = [-1, 0]
   reaction.set_charges(rcharges, pcharges)

   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   1 C9H8O4 + 1 OH => 1 C2H3O2 + 1 C7H6O3 | Charge = 0.00
   Coefficients:  [1 1 1 1]

If the charges provided cannot result in a neutral reaction equation, no solution will be found

.. code:: ipython3

   rcharges = [0, -1]
   pcharges = [0, 0]
   reaction.set_charges(rcharges, pcharges)

   coeffs = reaction.balance()
   print(reaction)
   print("Coefficients: ", coeffs)

::

   Inconsistent system: No solution
   Coefficients:  None

Test 3. Reaction with many possible products
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Aspirin can react with water to form a very wide range of products, but most likely not all in the same reaction. Here, we supply many possible products at once, and then balance the equation towards each product in turn.

.. code:: ipython3

   import numpy
   reactant_smiles = ["O", "CC(=O)Oc1ccccc1C(=O)O"]
   product_smiles = [
               "CC(=O)O",
               "O=C(O)c1ccccc1O", 
               "OO",
               "O=C(O)O",
               "O=C=O",
               "CO",
               "C",
               "O=C(O)C1=CC(O)C=CC1=O", 
               "O=C(O)C1=CC(O)C(O)C=C1O",
               "O=C(O)C(C=CO)=C(O)C=CO",
               "O=C(CO)Oc1ccccc1",
               "Oc1ccccc1", 
               "O=C1C=CC=CC1", 
               "OC1=CC(O)C=CC1",
               "C=CC=CC=C=O",
               "C1=C=CC=CC=1",
               "O=C=CO",
               "O=CO",
               "O=C=C(O)O",
               "O=C=C=O",
               "O=C1OC1=O",
               "O=COc1ccccc1",
               "O=C1CC=CC(O)C1",
               "C=CC(O)C(O)C=C=O",
               "OC1=CC=CC(O)C1=C(O)O",
               "O=C(O)C12C(=O)C1C=CC2O",
               "O=C(O)CO",
               "O=CC(=O)O",
               "O=C1C=CC(O)=C(C1)C(=O)O",
               "C=CC=CC(=O)O",
           ]

.. code:: ipython3

   reactants = [from_smiles(smiles) for smiles in reactant_smiles]
   products = [from_smiles(smiles) for smiles in product_smiles]

   # Create the Reaction object with all the molecules
   reaction = ReactionEquation(reactants, products)

   print("Starting loop over products..")
   nmols = len(reactants) + len(products)
   nreactants = len(reactants)
   for iprod, _ in enumerate(products):
       print(f"{iprod:8d} {product_smiles[iprod]:30s}: ", end="")
       min_coeffs = numpy.zeros(nmols)
       min_coeffs[nreactants + iprod] = 1
       coeffs = reaction.balance(min_coeffs)
       print(reaction)

::

   Starting loop over products..
          0 CC(=O)O                       : 1 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 C7H6O3 | Charge = 0.00
          1 O=C(O)c1ccccc1O               : 1 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 C7H6O3 | Charge = 0.00
          2 OO                            : 2 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 H2O2 + 1 C7H6O2 | Charge = 0.00
          3 O=C(O)O                       : 1 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 CH2O3 + 1 C6H4 | Charge = 0.00
          4 O=C=O                         : 1 C9H8O4 => 1 C2H4O2 + 1 CO2 + 1 C6H4 | Charge = 0.00
          5 CO                            : 1 C9H8O4 => 1 CH4O + 1 C6H4 + 1 C2O3 | Charge = 0.00
          6 C                             : 1 H2O + 1 C9H8O4 => 1 C7H6O3 + 1 CO2 + 1 CH4 | Charge = 0.00
          7 O=C(O)C1=CC(O)C=CC1=O         : 2 H2O + 1 C9H8O4 => 1 CH4 + 1 C7H6O4 + 1 CH2O2 | Charge = 0.00
          8 O=C(O)C1=CC(O)C(O)C=C1O       : 3 H2O + 1 C9H8O4 => 1 CH4 + 1 C7H8O5 + 1 CH2O2 | Charge = 0.00
          9 O=C(O)C(C=CO)=C(O)C=CO        : 3 H2O + 1 C9H8O4 => 1 CH4 + 1 C7H8O5 + 1 CH2O2 | Charge = 0.00
         10 O=C(CO)Oc1ccccc1              : 1 H2O + 1 C9H8O4 => 1 C8H8O3 + 1 CH2O2 | Charge = 0.00
         11 Oc1ccccc1                     : 1 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 CO2 + 1 C6H6O | Charge = 0.00
         12 O=C1C=CC=CC1                  : 1 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 CO2 + 1 C6H6O | Charge = 0.00
         13 OC1=CC(O)C=CC1                : 2 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 CO2 + 1 C6H8O2 | Charge = 0.00
         14 C=CC=CC=C=O                   : 1 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 CO2 + 1 C6H6O | Charge = 0.00
         15 C1=C=CC=CC=1                  : 1 C9H8O4 => 1 C2H4O2 + 1 CO2 + 1 C6H4 | Charge = 0.00
         16 O=C=CO                        : 1 C9H8O4 => 1 C2H2O2 + 1 C7H6O2 | Charge = 0.00
         17 O=CO                          : 1 H2O + 1 C9H8O4 => 1 C8H8O3 + 1 CH2O2 | Charge = 0.00
         18 O=C=C(O)O                     : 1 H2O + 2 C9H8O4 => 2 C8H8O3 + 1 C2H2O3 | Charge = 0.00
         19 O=C=C=O                       : 2 C9H8O4 => 2 C8H8O3 + 1 C2O2 | Charge = 0.00
         20 O=C1OC1=O                     : 1 C9H8O4 => 1 CH4O + 1 C6H4 + 1 C2O3 | Charge = 0.00
         21 O=COc1ccccc1                  : 1 C9H8O4 => 1 C2H2O2 + 1 C7H6O2 | Charge = 0.00
         22 O=C1CC=CC(O)C1                : 2 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 CO2 + 1 C6H8O2 | Charge = 0.00
         23 C=CC(O)C(O)C=C=O              : 2 H2O + 1 C9H8O4 => 1 CH4O + 1 C2O2 + 1 C6H8O3 | Charge = 0.00
         24 OC1=CC=CC(O)C1=C(O)O          : 2 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 C7H8O4 | Charge = 0.00
         25 O=C(O)C12C(=O)C1C=CC2O        : 2 H2O + 1 C9H8O4 => 1 CH4 + 1 CH2O2 + 1 C7H6O4 | Charge = 0.00
         26 O=C(O)CO                      : 1 H2O + 1 C9H8O4 => 1 C7H6O2 + 1 C2H4O3 | Charge = 0.00
         27 O=CC(=O)O                     : 1 H2O + 2 C9H8O4 => 2 C8H8O3 + 1 C2H2O3 | Charge = 0.00
         28 O=C1C=CC(O)=C(C1)C(=O)O       : 2 H2O + 1 C9H8O4 => 1 CH4 + 1 CH2O2 + 1 C7H6O4 | Charge = 0.00
         29 C=CC=CC(=O)O                  : 2 H2O + 1 C9H8O4 => 1 C2H4O2 + 1 C2H2O2 + 1 C5H6O2 | Charge = 0.00

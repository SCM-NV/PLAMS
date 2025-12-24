Worked Example
--------------

Balance simple reactions
~~~~~~~~~~~~~~~~~~~~~~~~

Here we present eleven different sets of reactants and products, not necessarily based on real reactions. The molecules can be passed as formulas, or as PLAMS Molecule objects.

.. code:: ipython3

   from scm.plams.tools.reaction import ReactionEquation, balance

Example: Aspirin in water
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["C9H8O4", "H2O"]
   products = ["C2H4O2", "C7H6O3"]

   reaction = balance(reactants, products)
   print(reaction)

::

   1 C9H8O4 + 1 H2O => 1 C2H4O2 + 1 C7H6O3

``reaction`` is of type ``ReactionEquation`` and has the following attributes:

.. code:: ipython3

   print(f"{reaction.coeffs=}")
   print(f"{reaction.message=}")
   print(f"{reaction.reactants=}")
   print(f"{reaction.products=}")

::

   reaction.coeffs=array([1, 1, 1, 1])
   reaction.message='Success'
   reaction.reactants=['C9H8O4', 'H2O']
   reaction.products=['C2H4O2', 'C7H6O3']

Specify reactants/products as Species
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An alternative to the string representation of reactants and products is to use Species, which also allows you to set the molecular charges and minimum coefficients (see later examples):

.. code:: ipython3

   from scm.plams.tools.reaction import Species as S

   reactants = [S("C9H8O4", charge=0), S("H2O", charge=0)]
   products = [S("C2H4O2", charge=0), S("C7H6O3", charge=0)]
   reaction = balance(reactants, products)
   print(reaction)

::

   1 C9H8O4 + 1 H2O => 1 C2H4O2 + 1 C7H6O3

Initialize species from SMILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   from scm.plams.tools.reaction import Species as S

   reactants = [S.from_smiles("C"), S.from_smiles("O=O")]
   products = [S.from_smiles("O=C=O"), S.from_smiles("O")]
   reaction = balance(reactants, products)
   print(reaction)

::

   1 CH4 + 2 O2 => 1 CO2 + 2 H2O

You can also format the reaction with the SMILES strings (requires that all reactants and products have the ``smiles`` attribute set):

.. code:: ipython3

   print(f"{reaction:smiles}")

::

   1 C + 2 O=O => 1 O=C=O + 2 O

Superfluous reactants and products (multiple possible solutions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When there are multiple possible solutions, a single one is chosen.

.. code:: ipython3

   reactants = ["C9H8O4", "H2O", "HO", "HO2"]
   products = ["C2H4O2", "C7H6O3", "H2O2"]
   reaction = balance(reactants, products)
   print(reaction)
   print("Some coefficients are zero:")
   print(f"{reaction.coeffs=}")

::

   1 C9H8O4 + 1 H2O => 1 C2H4O2 + 1 C7H6O3
   Some coefficients are zero:
   reaction.coeffs=array([1, 1, 0, 0, 1, 1, 0])

Exampe that is not easily balanced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["C9O4H6", "OH"]
   products = ["C2O2H3", "C7O3H2"]
   reaction = balance(reactants, products)
   print(reaction)

::

   9 C9O4H6 + 25 OH => 23 C2O2H3 + 5 C7O3H2

Reaction that cannot be balanced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   reactants = ["FeS2", "HNO3"]
   products = ["Fe2S3O12", "NO", "H2SO4"]
   reaction = balance(reactants, products)
   print(reaction)
   print(f"{reaction.coeffs=}")
   print(f"{reaction.message=}")

::

   Inconsistent system: No solution
   reaction.coeffs=None
   reaction.message='Empty nullspace'

Sympy method
~~~~~~~~~~~~

By default, a native method is used to compute the nullspace. Optionally, this can be done with sympy. The result should be the same.

.. code:: ipython3

   reactants = ["C9H8O4", "H2O"]
   products = ["C2H4O2", "C7H6O3"]

   reaction = balance(reactants, products, method="sympy")
   print(reaction)

::

   1 C9H8O4 + 1 H2O => 1 C2H4O2 + 1 C7H6O3

Balance reactions with charged species
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Charges are printed within square brackets:

.. code:: ipython3

   from scm.plams.tools.reaction import Species as S

   reactants = [S("OH", charge=-1), S("H3O", charge=+1)]
   products = [S("H2O", charge=0)]
   reaction = balance(reactants, products)
   print(reaction)

::

   1 OH[-1] + 1 H3O[+1] => 2 H2O

If the charges cannot be balanced, no solution will be found:

.. code:: ipython3

   from scm.plams.tools.reaction import Species as S

   reactants = [S("OH", charge=-1), S("H3O", charge=0)]  # "neutral" H3O!
   products = [S("H2O", charge=0)]
   reaction = balance(reactants, products)
   print(reaction)

::

   Inconsistent system: No solution

The charges are also inferred from the SMILES strings:

.. code:: ipython3

   from scm.plams.tools.reaction import Species as S

   reactants = [S.from_smiles("[H+]"), S.from_smiles("[OH-]")]
   products = [S.from_smiles("O")]
   reaction = balance(reactants, products)
   print(reaction)

::

   1 H[+1] + 1 HO[-1] => 1 H2O

.. code:: ipython3

   print(f"{reaction:smiles}")

::

   1 [H+] + 1 [OH-] => 1 O

Reaction with many possible products: Set minimum coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can set minimum values for the coefficients with the ``min_coeff`` attribute for a ``Species``. This is useful if you have a list of many possible products, and want to find a balanced reaction that includes one or more specific products:

.. code:: ipython3

   from scm.plams.tools.reaction import Species as S

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
   ]

   reactants = [S.from_smiles(x, min_coeff=1) for x in reactant_smiles]  # both reactants must appear
   products = [S.from_smiles(x) for x in product_smiles]
   products[1].min_coeff = 1  # the second product must appear
   products[2].min_coeff = 3  # the third product must have a coefficient of at least 3
   reaction = balance(reactants, products)
   print(reaction)

::

   6 H2O + 1 C9H8O4 => 1 C7H6O3 + 3 H2O2 + 1 CH4O + 1 CH4

.. code:: ipython3

   print(f"{reaction:smiles}")

::

   6 O + 1 CC(=O)Oc1ccccc1C(=O)O => 1 O=C(O)c1ccccc1O + 3 OO + 1 CO + 1 C

Direct usage of the ReactionEquation class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also directly use the ReactionEquation class. This can be useful and save some time if you want to loop over many possible values of ``min_coeffs``, for example.

In this case, do not use ``Species`` but provide formulas or PLAMS Molecules directly, use the ``set_charges`` method to set charges, and set the ``min_coeffs`` array:

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

   from scm.plams import from_smiles

   reactants = [from_smiles(smiles) for smiles in reactant_smiles]  # PLAMS Molecules
   products = [from_smiles(smiles) for smiles in product_smiles]  # PLAMS Molecules

   # Create the Reaction object with all the molecules
   reaction = ReactionEquation(reactants, products)

   print("Starting loop over products..")
   nmols = len(reactants) + len(products)
   nreactants = len(reactants)
   for iprod, _ in enumerate(products):
       print(f"{product_smiles[iprod]:>25s}: ", end="")
       min_coeffs = numpy.zeros(nmols)
       min_coeffs[nreactants + iprod] = 1
       reaction.balance(min_coeffs=min_coeffs)
       print(f"{reaction:smiles}")

::

   Starting loop over products..
                     CC(=O)O: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C(O)C1=CC=CC=C1O
             O=C(O)c1ccccc1O: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C(O)C1=CC=CC=C1O
                          OO: 2 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 OO + 1 O=COC1C=CC=CC=1
                     O=C(O)O: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C(O)O + 1 C1=C=CC=CC=1
                       O=C=O: 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C=O + 1 C1=C=CC=CC=1
                          CO: 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CO + 1 C1=C=CC=CC=1 + 1 O=C1OC1=O
                           C: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 O=C(O)C1=CC=CC=C1O + 1 O=C=O + 1 C
       O=C(O)C1=CC(O)C=CC1=O: 2 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 C + 1 O=C(O)C1=CC(O)C=CC1=O + 1 O=CO
     O=C(O)C1=CC(O)C(O)C=C1O: 3 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 C + 1 O=C(O)C1=CC(O)C(O)C=C1O + 1 O=CO
      O=C(O)C(C=CO)=C(O)C=CO: 3 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 C + 1 O=C(O)C(C=CO)=C(O)C=CO + 1 O=CO
            O=C(CO)Oc1ccccc1: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 O=C(CO)OC1C=CC=CC=1 + 1 O=CO
                   Oc1ccccc1: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C=O + 1 OC1C=CC=CC=1
                O=C1C=CC=CC1: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C=O + 1 O=C1C=CC=CC1
              OC1=CC(O)C=CC1: 2 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C=O + 1 OC1=CC(O)C=CC1
                 C=CC=CC=C=O: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C=O + 1 C=CC=CC=C=O
                C1=C=CC=CC=1: 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C=O + 1 C1=C=CC=CC=1
                      O=C=CO: 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 O=C=CO + 1 O=COC1C=CC=CC=1
                        O=CO: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 O=C(CO)OC1C=CC=CC=1 + 1 O=CO
                   O=C=C(O)O: 1 O + 2 CC(=O)OC1=CC=CC=C1C(=O)O => 2 O=C(CO)OC1C=CC=CC=1 + 1 O=C=C(O)O
                     O=C=C=O: 2 CC(=O)OC1=CC=CC=C1C(=O)O => 2 O=C(CO)OC1C=CC=CC=1 + 1 O=C=C=O
                   O=C1OC1=O: 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CO + 1 C1=C=CC=CC=1 + 1 O=C1OC1=O
                O=COc1ccccc1: 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 O=C=CO + 1 O=COC1C=CC=CC=1
              O=C1CC=CC(O)C1: 2 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C=O + 1 O=C1CC=CC(O)C1
            C=CC(O)C(O)C=C=O: 2 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CO + 1 O=C=C=O + 1 C=CC(O)C(O)C=C=O
        OC1=CC=CC(O)C1=C(O)O: 2 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 OC1=CC=CC(O)C1=C(O)O
      O=C(O)C12C(=O)C1C=CC2O: 2 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 C + 1 O=CO + 1 O=C(O)C12C(=O)C1C=CC2O
                    O=C(O)CO: 1 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 O=COC1C=CC=CC=1 + 1 O=C(O)CO
                   O=CC(=O)O: 1 O + 2 CC(=O)OC1=CC=CC=C1C(=O)O => 2 O=C(CO)OC1C=CC=CC=1 + 1 O=CC(=O)O
     O=C1C=CC(O)=C(C1)C(=O)O: 2 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 C + 1 O=CO + 1 O=C1C=CC(O)=C(C1)C(=O)O
                C=CC=CC(=O)O: 2 O + 1 CC(=O)OC1=CC=CC=C1C(=O)O => 1 CC(=O)O + 1 O=C=CO + 1 C=CC=CC(=O)O

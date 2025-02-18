Worked Example
--------------

Task, Engine and Property Blocks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    from scm.plams import Settings, AMSJob

To start with, in order to see the input that will be passed to AMS from
the Settings, make a function to create an AMSJob and print the input:

.. code:: ipython3

    def print_input(settings):
        print(AMSJob(settings=settings).get_input())

Tasks can be specified in the settings under the ``input.AMS.Task`` key:

.. code:: ipython3

    go_settings = Settings()
    go_settings.input.AMS.Task = "GeometryOptimization"
    go_settings.input.AMS.GeometryOptimization.Convergence.Gradients = 1e-5
    print_input(go_settings)


.. parsed-literal::

    GeometryOptimization
      Convergence
        Gradients 1e-05
      End
    End
    
    Task GeometryOptimization
    
    


Properties can be specified under the ``input.AMS.Properties`` key:

.. code:: ipython3

    nm_settings = Settings()
    nm_settings.input.ams.Properties.NormalModes = "Yes"
    print_input(nm_settings)


.. parsed-literal::

    Properties
      NormalModes Yes
    End
    
    


Engine settings can be specified under the ``input.AMS.<Engine>`` key,
for the engine of interest:

.. code:: ipython3

    lda_settings = Settings()
    lda_settings.input.ADF.Basis.Type = "DZP"
    print_input(lda_settings)


.. parsed-literal::

    
    Engine ADF
      Basis
        Type DZP
      End
    EndEngine
    
    


Combining Settings
~~~~~~~~~~~~~~~~~~

Settings objects can also be combined for easy reuse and job setup.
Settings can be merged using the ``+`` and ``+=`` operators.

.. code:: ipython3

    settings = go_settings + lda_settings + nm_settings
    print_input(settings)


.. parsed-literal::

    GeometryOptimization
      Convergence
        Gradients 1e-05
      End
    End
    
    Properties
      NormalModes Yes
    End
    
    Task GeometryOptimization
    
    
    Engine ADF
      Basis
        Type DZP
      End
    EndEngine
    
    


Note however that this merge is a “soft” update, so values of existing
keys will not be overwritten:

.. code:: ipython3

    pbe_settings = Settings()
    pbe_settings.input.ADF.Basis.Type = "TZP"
    pbe_settings.input.ADF.xc.gga = "pbe"
    settings += pbe_settings
    print_input(settings)


.. parsed-literal::

    GeometryOptimization
      Convergence
        Gradients 1e-05
      End
    End
    
    Properties
      NormalModes Yes
    End
    
    Task GeometryOptimization
    
    
    Engine ADF
      Basis
        Type DZP
      End
      xc
        gga pbe
      End
    EndEngine
    
    


To achieve “hard update” behaviour, the ``update`` method can be used,
which overwrites existing keys:

.. code:: ipython3

    settings.update(pbe_settings)
    print_input(settings)


.. parsed-literal::

    GeometryOptimization
      Convergence
        Gradients 1e-05
      End
    End
    
    Properties
      NormalModes Yes
    End
    
    Task GeometryOptimization
    
    
    Engine ADF
      Basis
        Type TZP
      End
      xc
        gga pbe
      End
    EndEngine
    
    


Settings can also be removed using the ``-`` and ``-=`` operators:

.. code:: ipython3

    settings -= nm_settings
    print_input(settings)


.. parsed-literal::

    GeometryOptimization
      Convergence
        Gradients 1e-05
      End
    End
    
    Properties
    End
    
    Task GeometryOptimization
    
    
    Engine ADF
      Basis
        Type TZP
      End
      xc
        gga pbe
      End
    EndEngine
    
    


Multiple values in a settings block can be configured using a list:

.. code:: ipython3

    hybrid_settings = go_settings.copy()
    hybrid_settings.input.AMS.Hybrid.Energy.Term = []
    for i in range(5):
        factor = (-1) ** (i % 2) * 1.0
        region = "*" if i == 0 else "one" if i < 3 else "two"
        engine_id = "adf-lda" if i == 0 or factor == -1 else "adf-gga"
        term = Settings({"Factor": factor, "Region": region, "EngineID": engine_id})
        hybrid_settings.input.AMS.Hybrid.Energy.Term.append(term)
    hybrid_settings.input.AMS.Hybrid.Engine = [lda_settings.input.ADF.copy(), pbe_settings.input.ADF.copy()]
    hybrid_settings.input.AMS.Hybrid.Engine[0]._h = "ADF adf-lda"
    hybrid_settings.input.AMS.Hybrid.Engine[1]._h = "ADF adf-gga"
    print_input(hybrid_settings)


.. parsed-literal::

    GeometryOptimization
      Convergence
        Gradients 1e-05
      End
    End
    
    Hybrid
      Energy
        Term
          EngineID adf-lda
          Factor 1.0
          Region *
        End
        Term
          EngineID adf-lda
          Factor -1.0
          Region one
        End
        Term
          EngineID adf-gga
          Factor 1.0
          Region one
        End
        Term
          EngineID adf-lda
          Factor -1.0
          Region two
        End
        Term
          EngineID adf-gga
          Factor 1.0
          Region two
        End
      End
      Engine ADF adf-lda
        Basis
          Type DZP
        End
      EndEngine
      Engine ADF adf-gga
        Basis
          Type TZP
        End
        xc
          gga pbe
        End
      EndEngine
    
    End
    
    Task GeometryOptimization
    
    


Note also in the example below, the use of the special ``_h`` “header”
key, which can be used to add data to the header line for a block.

Nested Keys
~~~~~~~~~~~

It can be useful to access values from a Settings object using “nested”
keys. These are tuples of keys, where each successive element of the
tuple corresponds to a further layer in the settings. Lists are
flattened so their elements can be accessed with the corresponding
index.

.. code:: ipython3

    list(hybrid_settings.nested_keys())




.. parsed-literal::

    [('input',),
     ('input', 'AMS'),
     ('input', 'AMS', 'Task'),
     ('input', 'AMS', 'GeometryOptimization'),
     ('input', 'AMS', 'GeometryOptimization', 'Convergence'),
     ('input', 'AMS', 'GeometryOptimization', 'Convergence', 'Gradients'),
     ('input', 'AMS', 'Hybrid'),
     ('input', 'AMS', 'Hybrid', 'Energy'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 0),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 0, 'Factor'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 0, 'Region'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 0, 'EngineID'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 1),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 1, 'Factor'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 1, 'Region'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 1, 'EngineID'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 2),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 2, 'Factor'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 2, 'Region'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 2, 'EngineID'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 3),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 3, 'Factor'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 3, 'Region'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 3, 'EngineID'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 4),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 4, 'Factor'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 4, 'Region'),
     ('input', 'AMS', 'Hybrid', 'Energy', 'Term', 4, 'EngineID'),
     ('input', 'AMS', 'Hybrid', 'Engine'),
     ('input', 'AMS', 'Hybrid', 'Engine', 0),
     ('input', 'AMS', 'Hybrid', 'Engine', 0, '_h'),
     ('input', 'AMS', 'Hybrid', 'Engine', 0, 'Basis'),
     ('input', 'AMS', 'Hybrid', 'Engine', 0, 'Basis', 'Type'),
     ('input', 'AMS', 'Hybrid', 'Engine', 1),
     ('input', 'AMS', 'Hybrid', 'Engine', 1, '_h'),
     ('input', 'AMS', 'Hybrid', 'Engine', 1, 'Basis'),
     ('input', 'AMS', 'Hybrid', 'Engine', 1, 'Basis', 'Type'),
     ('input', 'AMS', 'Hybrid', 'Engine', 1, 'xc'),
     ('input', 'AMS', 'Hybrid', 'Engine', 1, 'xc', 'gga')]



.. code:: ipython3

    hybrid_settings.get_nested(("input", "AMS", "Task"))




.. parsed-literal::

    'GeometryOptimization'



.. code:: ipython3

    if hybrid_settings.contains_nested(("input", "AMS", "Hybrid", "Engine", 0)):
        hybrid_settings.set_nested(("input", "AMS", "Hybrid", "Engine", 0, "Basis", "Type"), "TZP")
    print(hybrid_settings.get_nested(("input", "AMS", "Hybrid", "Engine", 0, "Basis")))


.. parsed-literal::

    Type: 	TZP
    


Comparison
~~~~~~~~~~

Two settings objects can be compared to check the differences between
them. The result will show the nested key and value of any added,
removed and modified entries.

.. code:: ipython3

    import os
    
    settings1 = go_settings + lda_settings + nm_settings
    settings2 = go_settings.copy()
    settings2.input.AMS.Task = "SinglePoint"
    settings2.input.DFTB.Model = "GFN1-xTB"
    comparison = settings2.compare(settings1)
    print(
        f"Items in settings2 not in settings1:{os.linesep}{os.linesep.join(f'  - {k}: {v}' for k, v in comparison['added'].items())}"
    )
    print(
        f"Items in settings1 not in settings2:{os.linesep}{os.linesep.join(f'  - {k}: {v}' for k, v in comparison['removed'].items())}"
    )
    print(
        f"Items modified from settings1 to settings2:{os.linesep}{os.linesep.join(f'  - {k}: {v[1]} -> {v[0]}' for k, v in comparison['modified'].items())}"
    )


.. parsed-literal::

    Items in settings2 not in settings1:
      - ('input', 'DFTB', 'Model'): GFN1-xTB
    Items in settings1 not in settings2:
      - ('input', 'ADF', 'Basis', 'Type'): DZP
      - ('input', 'AMS', 'Properties', 'NormalModes'): Yes
    Items modified from settings1 to settings2:
      - ('input', 'AMS', 'Task'): GeometryOptimization -> SinglePoint



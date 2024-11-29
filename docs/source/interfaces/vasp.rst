VASP
-------------------------

(*contributed by*  `Patrick Melix <https://www.researchgate.net/profile/Patrick_Melix>`_\)

.. currentmodule:: scm.plams.interfaces.thirdparty.vasp

Results
~~~~~~~~~~~~~~~~~~~~

It is highly recommended to use the `ASE <https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html>`_ features for accessing the results of a |VASPJob| through the `vasprun.xml`.
Use something like the snippet below to create a dummy ASE calculator and retrieve the results you need.
Remember that ASE needs the path to the POTCARs in an environment variable. ::

    ase_calc = Vasp(directory=str(job.path), xc='PBE') # set xc to anything, just needed for the automatisms of ASE
    ase_calc.read()
    print(ase_calc.get_potential_energy(force_consistent=False))
    forces = ase_calc.get_forces()
    ...


API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: VASPJob()
.. autoclass:: VASPResults()

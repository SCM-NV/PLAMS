CP2K
-------------------------

(*contributed by* `Felipe Zapata <https://www.researchgate.net/profile/Felipe_Zapata>`_\, extended by `Patrick Melix <https://www.researchgate.net/profile/Patrick_Melix>`_\)

.. currentmodule:: scm.plams.interfaces.thirdparty.cp2k

Settings
~~~~~~~~~~~~~~

The cp2k_ input is rather complex one compared to other computational codes, but its input is structured as a set of nested block and sub-blocks that can be easily represented by the |Settings| class.
Like the other interfaces, the *CP2K* input file is generated using the *input* branch of the job settings.
For instance, a single point calculation for pentacene::

    penta = Settings()

    penta.input.force_eval.dft.basis_set_file_name = "/path/to/basis"
    penta.input.force_eval.dft.potential_set_file_name = "/path/to/potential"
    penta.input.force_eval.dft.mgrid.cutoff = 400
    penta.input.force_eval.dft.mgrid.ngrids = 4
    penta.input.force_eval.dft.poisson = Settings()
    penta.input.force_eval.dft.localize._h = "T"


    penta.input.force_eval.dft.qs.method = "GPW"
    penta.input.force_eval.dft.scf.eps_scf = 1e-6
    penta.input.force_eval.dft.scf.max_scf = 200
    penta.input.force_eval.dft.xc.xc_functional._h = "PBE"

    penta.input.force_eval.subsys.cell.A = '16.11886919 0.07814137 -0.697284243'
    penta.input.force_eval.subsys.cell.B = '-0.215317662 4.389405268 1.408951791'
    penta.input.force_eval.subsys.cell.C = '-0.216126961 1.732808365 9.74896108'
    penta.input.force_eval.subsys.cell.periodic = 'XYZ'
    penta.input.force_eval.subsys.kind.C.basis_set = "DZVP-MOLOPT-SR-GTH-Q4"
    penta.input.force_eval.subsys.kind.C.potential = "GTH-PBE-Q4"
    penta.input.force_eval.subsys.kind.H.basis_set = "DZVP-MOLOPT-SR-GTH-Q1"
    penta.input.force_eval.subsys.kind.H.potential = "GTH-PBE-q1"
    penta.input.force_eval.subsys.topology.coord_file_name = "./penta.xyz"
    penta.input.force_eval.subsys.topology.coordinate = "xyz"

    penta.input['global'].print_level = "low"
    penta.input['global'].project  = "example"
    penta.input['global'].run_type = "energy_force"

The input generated during the execution of the cp2k_ job is similar to: ::

    &FORCE_EVAL
      &DFT
        BASIS_SET_FILE_NAME  /path/to/basis
        &MGRID
          CUTOFF  400
          NGRIDS  4
        &END
        &POISSON
        &END
        &LOCALIZE T
        &END
        POTENTIAL_FILE_NAME  /path/to/potential
        &QS
          METHOD  GPW
        &END
        &SCF
          EPS_SCF  1e-06
          MAX_SCF  200
        &END
        &XC
          &XC_FUNCTIONAL PBE
          &END
        &END
      &END
      &SUBSYS
        &CELL
          A  16.11886919 0.07814137 -0.697284243
          B  -0.215317662 4.389405268 1.408951791
          C  -0.216126961 1.732808365 9.748961085
          PERIODIC  XYZ
        &END
        &KIND  C
          BASIS_SET  DZVP-MOLOPT-SR-GTH-q4
          POTENTIAL  GTH-PBE-q4
        &END
        &KIND  H
          BASIS_SET  DZVP-MOLOPT-SR-GTH-q1
          POTENTIAL  GTH-PBE-q1
        &END
        &TOPOLOGY
          COORD_FILE_NAME  ./geometry.xyz
          COORDINATE  XYZ
        &END
      &END
    &END

    &GLOBAL
      PRINT_LEVEL  LOW
      PROJECT  example
      RUN_TYPE  ENERGY_FORCE
    &END

PLAMS automatically creates the indented structure of the previous example together with the special character *&* at the beginning and end of each section, and finally the keyword *END* at the end of each section.

Notice that *CP2K* requires the explicit declaration of the basis set together with the charge and the name of the potential used for each one of the atoms.
In the previous example the basis for the carbon is *DZVP-MOLOPT-SR-GTH*, while the potential is *GTH-PBE* and the charge *q4*.

The input parser also allows for header values like in this example *&LOCALIZE T  &END* by using the ``<some_section>._h`` (for header) notation in the |Settings| instance.

Inclusion of external files by using the *@INCLUDE* notation of Cp2k is supported. Also the *@SET* and *@IF* keys can be used, just replace the *@* sign by ``at_`` in the definition of your |Settings|.
If you need some files to be copied to the actual execution directory, pass them to the constructor using the ``copy=`` option. See the API below.


Molecule parsing
~~~~~~~~~~~~~~~~~~~~

Molecules can be parsed into the input automatically by passing them to the |Cp2kJob| constructor.
Leaving the molecule unset is supported, make sure you then set your molecule by hand (e.g. using the *TOPOLOGY* section).
Parsing of the molecule can be avoided by setting ``job.settings.ignore_molecule = True`` in the |Settings| of the job.

Also, the simulation cell can be specified using the x, y, z vectors like in the previous example.
A cubic box can be easily specified by::

  penta.input.force_eval.subsys.cell.ABC = "[angstrom] 50 50 50"

That results in a simulation cube of 50 cubic angstroms.
For a more detailed description of the *cp2k* input see manual_.

.. _cp2k: https://www.cp2k.org/

.. _manual: https://manual.cp2k.org/#gsc.tab=0

Loading jobs
~~~~~~~~~~~~~~~~~~~

Calculations done without PLAMS can be loaded using the |load_external| functionality. The |Cp2kResults| class supports reading input files into |Settings| objects.

Just do ``Cp2kJob.load_external(path)`` to get the settings from the file.

Molecule loading
~~~~~~~~~~~~~~~~~~~~

A |Molecule| can be recreated from a |Settings| instance using the |Cp2kSettings2Mol| function.


API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: Cp2kJob()

.. autoclass:: Cp2kResults()

.. autofunction:: Cp2kSettings2Mol()

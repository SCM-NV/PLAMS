Worked Example
--------------

Create Example Jobs
~~~~~~~~~~~~~~~~~~~

To begin with, create a variety of AMS jobs with different settings, engines and calculation types.

.. code:: ipython3

   from scm.plams import from_smiles, AMSJob, PlamsError, Settings, Molecule, Atom
   from scm.libbase import ChemicalSystem
   from scm.input_classes.drivers import AMS
   from scm.input_classes.engines import DFTB
   from scm.utils.conversions import plams_molecule_to_chemsys


   def example_job_dftb(smiles, task, use_chemsys=False):
       # Generate molecule from smiles
       mol = from_smiles(smiles)
       if use_chemsys:
           mol = plams_molecule_to_chemsys(mol)

       # Set up calculation settings using PISA
       sett = Settings()
       sett.runscript.nproc = 1
       driver = AMS()
       driver.Task = task
       driver.Engine = DFTB()
       sett.input = driver
       return AMSJob(molecule=mol, settings=sett, name="dftb")


   def example_job_adf(smiles, task, basis, gga=None, use_chemsys=False):
       # Generate molecule from smiles
       mol = from_smiles(smiles)
       if use_chemsys:
           mol = plams_molecule_to_chemsys(mol)

       # Set up calculation settings using standard settings
       sett = Settings()
       sett.runscript.nproc = 1
       sett.input.AMS.Task = task
       sett.input.ADF.Basis.Type = basis
       if gga:
           sett.input.ADF.XC.GGA = gga
       return AMSJob(molecule=mol, settings=sett, name="adf")


   def example_job_neb(iterations, use_chemsys=False):
       # Set up molecules
       main_molecule = Molecule()
       main_molecule.add_atom(Atom(symbol="C", coords=(0, 0, 0)))
       main_molecule.add_atom(Atom(symbol="N", coords=(1.18, 0, 0)))
       main_molecule.add_atom(Atom(symbol="H", coords=(2.196, 0, 0)))
       final_molecule = main_molecule.copy()
       final_molecule.atoms[1].x = 1.163
       final_molecule.atoms[2].x = -1.078

       mol = {"": main_molecule, "final": final_molecule}

       if use_chemsys:
           mol = {k: plams_molecule_to_chemsys(v) for k, v in mol.items()}

       # Set up calculation settings
       sett = Settings()
       sett.runscript.nproc = 1
       sett.input.ams.Task = "NEB"
       sett.input.ams.NEB.Images = 9
       sett.input.ams.NEB.Iterations = iterations
       sett.input.DFTB

       return AMSJob(molecule=mol, settings=sett, name="neb")

Now, run a selection of them.

.. code:: ipython3

   from scm.plams import config, JobRunner

   config.default_jobrunner = JobRunner(parallel=True, maxthreads=8)

   smiles = ["CC", "C", "O", "CO"]
   tasks = ["SinglePoint", "GeometryOptimization"]
   engines = ["DFTB", "ADF"]
   jobs = []
   for i, s in enumerate(smiles):
       for j, t in enumerate(tasks):
           job_dftb = example_job_dftb(s, t, use_chemsys=i % 2)
           job_adf1 = example_job_adf(s, t, "DZ", use_chemsys=True)
           job_adf2 = example_job_adf(s, t, "TZP", "PBE")
           jobs += [job_dftb, job_adf1, job_adf2]

   job_neb1 = example_job_neb(10)
   job_neb2 = example_job_neb(100, use_chemsys=True)
   jobs += [job_neb1, job_neb2]

   for j in jobs:
       j.run()

::

   [15.05|09:44:38] JOB dftb STARTED
   [15.05|09:44:38] JOB adf STARTED
   [15.05|09:44:38] JOB adf STARTED
   [15.05|09:44:38] JOB dftb STARTED
   [15.05|09:44:38] Renaming job adf to adf.002
   [15.05|09:44:38] JOB adf STARTED
   [15.05|09:44:38] JOB adf STARTED
   [15.05|09:44:38] JOB dftb STARTED
   [15.05|09:44:38] JOB adf STARTED
   [15.05|09:44:38] JOB dftb RUNNING
   ... (PLAMS log lines truncated) ...

Job Analysis
~~~~~~~~~~~~

Adding and Loading Jobs
~~~~~~~~~~~~~~~~~~~~~~~

Jobs can be loaded by passing job objects directly, or loading from a path.

.. code:: ipython3

   from scm.plams import JobAnalysis

.. code:: ipython3

   ja = JobAnalysis(jobs=jobs[:10], paths=[j.path for j in jobs[10:-2]])

::

   [15.05|09:44:48] Waiting for job adf.003 to finish
   [15.05|09:44:49] JOB neb FINISHED
   [15.05|09:44:49] Job neb reported errors. Please check the output
   [15.05|09:44:49] JOB neb FAILED
   [15.05|09:44:49] Job neb reported errors. Please check the output
   [15.05|09:44:49] Error message for job neb was:
       NEB optimization did NOT converge
   [15.05|09:44:49] Job neb reported errors. Please check the output
   [15.05|09:44:49] Job neb reported errors. Please check the output
   [15.05|09:44:49] JOB neb.002 FINISHED
   [15.05|09:44:49] JOB neb.002 SUCCESSFUL
   [15.05|09:44:52] JOB adf.014 FINISHED
   [15.05|09:44:52] JOB adf.014 SUCCESSFUL
   [15.05|09:44:53] JOB adf.008 FINISHED
   [15.05|09:44:53] JOB adf.008 SUCCESSFUL
   [15.05|09:44:57] JOB adf.015 FINISHED
   [15.05|09:44:57] JOB adf.015 SUCCESSFUL
   ... (PLAMS log lines truncated) ...
   [15.05|09:45:00] Waiting for job adf.004 to finish

Jobs can also be added or removed after initialization.

.. code:: ipython3

   ja = ja.add_job(jobs[-2]).load_job(jobs[-1].path)
   ja.display_table()

======================================================= ======== ===== ===== =================================
Path                                                    Name     OK    Check ErrorMsg
======================================================= ======== ===== ===== =================================
/path/plams/examples/JobAnalysis/plams_workdir/dftb     dftb     True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf      adf      True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.002  adf.002  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/dftb.002 dftb.002 True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.003  adf.003  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.004  adf.004  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/dftb.003 dftb.003 True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.005  adf.005  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.006  adf.006  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/dftb.004 dftb.004 True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.007  adf.007  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.008  adf.008  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/dftb.005 dftb.005 True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.009  adf.009  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.010  adf.010  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/dftb.006 dftb.006 True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.011  adf.011  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.012  adf.012  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/dftb.007 dftb.007 True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.013  adf.013  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.014  adf.014  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/dftb.008 dftb.008 True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.015  adf.015  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/adf.016  adf.016  True  True  None
/path/plams/examples/JobAnalysis/plams_workdir/neb      neb      False False NEB optimization did NOT converge
/path/plams/examples/JobAnalysis/plams_workdir/neb.002  neb.002  True  True  None
======================================================= ======== ===== ===== =================================

Adding and Removing Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~

A range of common standard fields can be added with the ``add_standard_field(s)`` methods. In addition, fields deriving from the job settings can be added with the ``add_settings_input_fields`` method, and fields from the output rkfs with the ``add_rkf_field`` method. Custom fields can also be added with the ``add_field`` method, by defining a field key, value accessor and optional arguments like display name and value formatting.

Fields can be removed by calling ``remove_field`` with the corresponding field key.

.. code:: ipython3

   ja = (
       ja.remove_field("Path")
       .add_standard_fields(["Formula", "Smiles", "CPUTime", "SysTime"])
       .add_rkf_field("General", "engine")
       .add_settings_input_fields()
       .add_field("Energy", lambda j: j.results.get_energy(unit="kJ/mol"), display_name="Energy [kJ/mol]", fmt=".2f")
   )
   ja.display_table(max_rows=5)

======= ===== ===== ================================= ================= ================= ========= ======== ================ ==================== ================= ============= ================= ===================== ===============
Name    OK    Check ErrorMsg                          Formula           Smiles            CPUTime   SysTime  AmsGeneralEngine InputAmsTask         InputAdfBasisType InputAdfXcGga InputAmsNebImages InputAmsNebIterations Energy [kJ/mol]
======= ===== ===== ================================= ================= ================= ========= ======== ================ ==================== ================= ============= ================= ===================== ===============
dftb    True  True  None                              C2H6              CC                0.189760  0.047422 dftb             SinglePoint          None              None          None              None                  -19594.01
adf     True  True  None                              C2H6              CC                3.702566  0.122056 adf              SinglePoint          DZ                None          None              None                  -3973.29
…       …     …     …                                 …                 …                 …         …        …                …                    …                 …             …                 …                     …
adf.016 True  True  None                              CH4O              CO                18.829835 0.833251 adf              GeometryOptimization TZP               PBE           None              None                  -2900.38
neb     False False NEB optimization did NOT converge : CHN, final: CHN : C=N, final: C#N 0.471924  0.050081 dftb             NEB                  None              None          9                 10                    None
neb.002 True  True  None                              : CHN, final: CHN : C=N, final: C#N 0.730961  0.090290 dftb             NEB                  None              None          9                 100                   -14936.53
======= ===== ===== ================================= ================= ================= ========= ======== ================ ==================== ================= ============= ================= ===================== ===============

In addition to the fluent syntax, both dictionary and dot syntaxes are also supported for adding and removing fields.

.. code:: ipython3

   import numpy as np

   ja["AtomType"] = lambda j: [at.symbol for at in j.results.get_main_molecule()]
   ja.Charge = lambda j: j.results.get_charges()
   ja.AtomCoords = lambda j: [np.array(at.coords) for at in j.results.get_main_molecule()]

   del ja["Check"]
   del ja.SysTime

   ja.display_table(max_rows=5, max_col_width=30)

======= ===== =============================== ================= ================= ========= ================ ==================== ================= ============= ================= ===================== =============== =============================== =============================== ===============================
Name    OK    ErrorMsg                        Formula           Smiles            CPUTime   AmsGeneralEngine InputAmsTask         InputAdfBasisType InputAdfXcGga InputAmsNebImages InputAmsNebIterations Energy [kJ/mol] AtomType                        Charge                          AtomCoords
======= ===== =============================== ================= ================= ========= ================ ==================== ================= ============= ================= ===================== =============== =============================== =============================== ===============================
dftb    True  None                            C2H6              CC                0.189760  dftb             SinglePoint          None              None          None              None                  -19594.01       [‘C’, ‘C’, ‘H’, ‘H’, ‘H’, ‘H’,… [-0.07293185 -0.07372966 0.02…  [array([-0.74763668, 0.041837…
adf     True  None                            C2H6              CC                3.702566  adf              SinglePoint          DZ                None          None              None                  -3973.29        [‘C’, ‘C’, ‘H’, ‘H’, ‘H’, ‘H’,… [-0.83243445 -0.83187828 0.27…  [array([-0.74763668, 0.041837…
…       …     …                               …                 …                 …         …                …                    …                 …             …                 …                     …               …                               …                               …
adf.016 True  None                            CH4O              CO                18.829835 adf              GeometryOptimization TZP               PBE           None              None                  -2900.38        [‘C’, ‘O’, ‘H’, ‘H’, ‘H’, ‘H’]  [ 0.58673094 -0.60299606 -0.10… [array([-0.36298962, -0.021487…
neb     False NEB optimization did NOT conve… : CHN, final: CHN : C=N, final: C#N 0.471924  dftb             NEB                  None              None          9                 10                    None            [‘C’, ‘N’, ‘H’]                 None                            [array([0.46884763, 0.20209473…
neb.002 True  None                            : CHN, final: CHN : C=N, final: C#N 0.730961  dftb             NEB                  None              None          9                 100                   -14936.53       [‘C’, ‘N’, ‘H’]                 [-0.00732595 -0.21157426 0.21…  [array([0.56218708, 0.20551051…
======= ===== =============================== ================= ================= ========= ================ ==================== ================= ============= ================= ===================== =============== =============================== =============================== ===============================

Processing Data
~~~~~~~~~~~~~~~

Once an initial analysis has been created, the data can be further processed, depending on the use case. For example, to inspect the difference between failed and successful jobs, jobs can be filtered down and irrelevant fields removed.

.. code:: ipython3

   ja_neb = (
       ja.filter_jobs(lambda data: data["InputAmsTask"] == "NEB")
       .remove_field("AtomCoords")
       .remove_uniform_fields(ignore_empty=True)
   )

   ja_neb.display_table()

======= ===== ======== =====================
Name    OK    CPUTime  InputAmsNebIterations
======= ===== ======== =====================
neb     False 0.471924 10
neb.002 True  0.730961 100
======= ===== ======== =====================

Another use case may be to analyze the results from one or more jobs. For this, it can be useful to utilize the ``expand`` functionality to convert job(s) to multiple rows. During this process, fields selected for expansion will have their values extracted into individual rows, whilst other fields have their values duplicated.

.. code:: ipython3

   ja_adf_expanded = (
       ja.filter_jobs(
           lambda data: data["InputAmsTask"] == "GeometryOptimization"
           and data["InputAdfBasisType"] is not None
           and data["Smiles"] == "O"
       )
       .expand_field("AtomType")
       .expand_field("Charge")
       .expand_field("AtomCoords")
       .remove_uniform_fields()
   )

   ja_adf_expanded.display_table()

======= ======== ================= ============= =============== ======== =================== ===============================================
Name    CPUTime  InputAdfBasisType InputAdfXcGga Energy [kJ/mol] AtomType Charge              AtomCoords
======= ======== ================= ============= =============== ======== =================== ===============================================
adf.011 2.697854 DZ                None          -1316.30        O        -0.8416865250737331 [-2.17062120e-04 3.82347777e-01 0.00000000e+00]
adf.011 2.697854 DZ                None          -1316.30        H        0.42084716070260286 [-0.81250923 -0.19167629 0. ]
adf.011 2.697854 DZ                None          -1316.30        H        0.4208393643711281  [ 0.8127263 -0.19067148 0. ]
adf.012 4.089876 TZP               PBE           -1363.77        O        -0.6739805275850443 [-2.46726007e-04 4.01580956e-01 0.00000000e+00]
adf.012 4.089876 TZP               PBE           -1363.77        H        0.33698188085180536 [-0.76455997 -0.2012764 0. ]
adf.012 4.089876 TZP               PBE           -1363.77        H        0.33699864673323343 [ 0.76480669 -0.20030455 0. ]
======= ======== ================= ============= =============== ======== =================== ===============================================

For more nested values, the depth of expansion can also be selected to further flatten the data.

.. code:: ipython3

   ja_adf_expanded2 = ja_adf_expanded.add_field(
       "Coord", lambda j: [("x", "y", "z") for _ in j.results.get_main_molecule()], expansion_depth=2
   ).expand_field("AtomCoords", depth=2)

   ja_adf_expanded2.display_table()

======= ======== ================= ============= =============== ======== =================== ======================= =====
Name    CPUTime  InputAdfBasisType InputAdfXcGga Energy [kJ/mol] AtomType Charge              AtomCoords              Coord
======= ======== ================= ============= =============== ======== =================== ======================= =====
adf.011 2.697854 DZ                None          -1316.30        O        -0.8416865250737331 -0.00021706211955194217 x
adf.011 2.697854 DZ                None          -1316.30        O        -0.8416865250737331 0.38234777653349844     y
adf.011 2.697854 DZ                None          -1316.30        O        -0.8416865250737331 0.0                     z
adf.011 2.697854 DZ                None          -1316.30        H        0.42084716070260286 -0.8125092343354401     x
adf.011 2.697854 DZ                None          -1316.30        H        0.42084716070260286 -0.19167629390344054    y
adf.011 2.697854 DZ                None          -1316.30        H        0.42084716070260286 0.0                     z
adf.011 2.697854 DZ                None          -1316.30        H        0.4208393643711281  0.8127262964549918      x
adf.011 2.697854 DZ                None          -1316.30        H        0.4208393643711281  -0.19067148263005784    y
adf.011 2.697854 DZ                None          -1316.30        H        0.4208393643711281  0.0                     z
adf.012 4.089876 TZP               PBE           -1363.77        O        -0.6739805275850443 -0.00024672600727009935 x
adf.012 4.089876 TZP               PBE           -1363.77        O        -0.6739805275850443 0.40158095623473306     y
adf.012 4.089876 TZP               PBE           -1363.77        O        -0.6739805275850443 0.0                     z
adf.012 4.089876 TZP               PBE           -1363.77        H        0.33698188085180536 -0.7645599672263915     x
adf.012 4.089876 TZP               PBE           -1363.77        H        0.33698188085180536 -0.2012764045590436     y
adf.012 4.089876 TZP               PBE           -1363.77        H        0.33698188085180536 0.0                     z
adf.012 4.089876 TZP               PBE           -1363.77        H        0.33699864673323343 0.7648066932336616      x
adf.012 4.089876 TZP               PBE           -1363.77        H        0.33699864673323343 -0.20030455167568945    y
adf.012 4.089876 TZP               PBE           -1363.77        H        0.33699864673323343 0.0                     z
======= ======== ================= ============= =============== ======== =================== ======================= =====

Expansion can be undone with the corresponding ``collapse`` method.

Fields can be also further filtered, modified or reordered to customize the analysis.

.. code:: ipython3

   ja_adf = (
       ja_adf_expanded2.collapse_field("AtomCoords")
       .collapse_field("Coord")
       .filter_fields(lambda vals: all([not isinstance(v, list) for v in vals]))  # remove arrays
       .remove_field("Name")
       .format_field("CPUTime", ".2f")
       .format_field("Charge", ".4f")
       .rename_field("InputAdfBasisType", "Basis")
       .reorder_fields(["AtomType", "Charge", "Energy"])
   )
   ja_adf.display_table()

======== ======= =============== ======= ===== =============
AtomType Charge  Energy [kJ/mol] CPUTime Basis InputAdfXcGga
======== ======= =============== ======= ===== =============
O        -0.8417 -1316.30        2.70    DZ    None
H        0.4208  -1316.30        2.70    DZ    None
H        0.4208  -1316.30        2.70    DZ    None
O        -0.6740 -1363.77        4.09    TZP   PBE
H        0.3370  -1363.77        4.09    TZP   PBE
H        0.3370  -1363.77        4.09    TZP   PBE
======== ======= =============== ======= ===== =============

Extracting Analysis Data
~~~~~~~~~~~~~~~~~~~~~~~~

Analysis data can be extracted in a variety of ways.

As has been demonstrated, a visual representation of the table can be easily generated using the ``to_table`` method (or ``display_table`` in a notebook). The format can be selected as markdown, html or rst. This will return the data with the specified display names and formatting.

.. code:: ipython3

   print(ja_adf.to_table(fmt="rst"))

::

   +----------+---------+-----------------+---------+-------+---------------+
   | AtomType | Charge  | Energy [kJ/mol] | CPUTime | Basis | InputAdfXcGga |
   +==========+=========+=================+=========+=======+===============+
   | O        | -0.8417 | -1316.30        | 2.70    | DZ    | None          |
   +----------+---------+-----------------+---------+-------+---------------+
   | H        | 0.4208  | -1316.30        | 2.70    | DZ    | None          |
   +----------+---------+-----------------+---------+-------+---------------+
   | H        | 0.4208  | -1316.30        | 2.70    | DZ    | None          |
   +----------+---------+-----------------+---------+-------+---------------+
   | O        | -0.6740 | -1363.77        | 4.09    | TZP   | PBE           |
   +----------+---------+-----------------+---------+-------+---------------+
   | H        | 0.3370  | -1363.77        | 4.09    | TZP   | PBE           |
   +----------+---------+-----------------+---------+-------+---------------+
   | H        | 0.3370  | -1363.77        | 4.09    | TZP   | PBE           |
   +----------+---------+-----------------+---------+-------+---------------+

Alternatively, raw data can be retrieved via the ``get_analysis`` method, which returns a dictionary of analysis keys to values.

.. code:: ipython3

   print(ja_adf.get_analysis())

::

   {'AtomType': ['O', 'H', 'H', 'O', 'H', 'H'], 'Charge': [-0.8416865250737331, 0.42084716070260286, 0.4208393643711281, -0.6739805275850443, 0.33698188085180536, 0.33699864673323343], 'Energy': [-1316.2997406426532, -1316.2997406426532, -1316.2997406426532, -1363.766294275197, -1363.766294275197, -1363.766294275197], 'CPUTime': [2.697854, 2.697854, 2.697854, 4.089876, 4.089876, 4.089876], 'InputAdfBasisType': ['DZ', 'DZ', 'DZ', 'TZP', 'TZP', 'TZP'], 'InputAdfXcGga': [None, None, None, 'PBE', 'PBE', 'PBE']}

Data can also be easily written to a csv file using ``to_csv_file``, to be exported to another program.

.. code:: ipython3

   csv_name = "./tmp.csv"
   ja_adf.to_csv_file(csv_name)

   with open(csv_name) as csv:
       print(csv.read())

::

   AtomType,Charge,Energy,CPUTime,InputAdfBasisType,InputAdfXcGga
   O,-0.8416865250737331,-1316.2997406426532,2.697854,DZ,
   H,0.42084716070260286,-1316.2997406426532,2.697854,DZ,
   H,0.4208393643711281,-1316.2997406426532,2.697854,DZ,
   O,-0.6739805275850443,-1363.766294275197,4.089876,TZP,PBE
   H,0.33698188085180536,-1363.766294275197,4.089876,TZP,PBE
   H,0.33699864673323343,-1363.766294275197,4.089876,TZP,PBE

Finally, for more complex data analysis, the results can be converted to a `pandas <https://pandas.pydata.org>`__ dataframe. This is recommended for more involved data manipulations, and can be installed using amspackages i.e. using the command: ``"${AMSBIN}/amspackages" install pandas``.

.. code:: ipython3

   try:
       import pandas

       df = ja_adf.to_dataframe()
       print(df)

   except ImportError:

       print(
           "Pandas not available. Please install with amspackages to run this example '${AMSBIN}/amspackages install pandas'"
       )

::

     AtomType    Charge       Energy   CPUTime InputAdfBasisType InputAdfXcGga
   0        O -0.841687 -1316.299741  2.697854                DZ          None
   1        H  0.420847 -1316.299741  2.697854                DZ          None
   2        H  0.420839 -1316.299741  2.697854                DZ          None
   3        O -0.673981 -1363.766294  4.089876               TZP           PBE
   4        H  0.336982 -1363.766294  4.089876               TZP           PBE
   5        H  0.336999 -1363.766294  4.089876               TZP           PBE

Additional Analysis Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``JobAnalysis`` class does have some additional built in methods to aid with job analysis.

For example, the ``get_timeline`` and ``display_timeline`` methods show pictorially when jobs started, how long they took to run and what their status is.

This can be useful for visualizing the dependencies of jobs. Here you can see that the first 8 jobs started running in parallel, due to the ``maxthreads`` constraint, and the remaining jobs waited before starting. Also that the penultimate job failed.

.. code:: ipython3

   ja.display_timeline(fmt="rst")

::

   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | JobName  | ↓2025-05-15 09:44:37 | ↓2025-05-15 09:44:49 | ↓2025-05-15 09:45:00 | ↓2025-05-15 09:45:12 | ↓2025-05-15 09:45:23 | WaitDuration | RunDuration | TotalDuration |
   +==========+======================+======================+======================+======================+======================+==============+=============+===============+
   | dftb     | ==>                  |                      |                      |                      |                      | 0s           | 0s          | 1s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf      | ========>            |                      |                      |                      |                      | 0s           | 4s          | 4s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.002  | ================>    |                      |                      |                      |                      | 0s           | 9s          | 9s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | dftb.002 | ==>                  |                      |                      |                      |                      | 0s           | 0s          | 1s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.003  | ==================== | ===================> |                      |                      |                      | 0s           | 22s         | 22s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.004  | ==================== | ==================== | ==================== | ===================* | >                    | 0s           | 45s         | 45s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | dftb.003 | ===>                 |                      |                      |                      |                      | 0s           | 1s          | 1s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.005  | ======>              |                      |                      |                      |                      | 0s           | 3s          | 3s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.006  | --=======>           |                      |                      |                      |                      | 0s           | 5s          | 5s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | dftb.004 | ..=>                 |                      |                      |                      |                      | 1s           | 0s          | 1s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.007  | ..-===========>      |                      |                      |                      |                      | 1s           | 6s          | 8s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.008  | ...+================ | =======>             |                      |                      |                      | 1s           | 13s         | 15s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | dftb.005 | ...---=>             |                      |                      |                      |                      | 2s           | 1s          | 3s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.009  | ......-===>          |                      |                      |                      |                      | 3s           | 2s          | 5s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.010  | .......-====>        |                      |                      |                      |                      | 4s           | 2s          | 6s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | dftb.006 | ........-=>          |                      |                      |                      |                      | 4s           | 1s          | 5s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.011  | .........-=====>     |                      |                      |                      |                      | 5s           | 3s          | 8s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.012  | ..........========>  |                      |                      |                      |                      | 5s           | 4s          | 10s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | dftb.007 | ..........-->        |                      |                      |                      |                      | 5s           | 1s          | 6s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.013  | ............-=====>  |                      |                      |                      |                      | 6s           | 3s          | 10s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.014  | .............-====== | ======>              |                      |                      |                      | 7s           | 7s          | 14s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | dftb.008 | ..............-=>    |                      |                      |                      |                      | 7s           | 1s          | 9s            |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.015  | ...............-==== | ===============>     |                      |                      |                      | 8s           | 11s         | 19s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | adf.016  | ................==== | ==================== | =============>       |                      |                      | 9s           | 20s         | 29s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | neb      | ................--== | X                    |                      |                      |                      | 9s           | 1s          | 10s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+
   | neb.002  | ..................== | >                    |                      |                      |                      | 10s          | 1s          | 11s           |
   +----------+----------------------+----------------------+----------------------+----------------------+----------------------+--------------+-------------+---------------+

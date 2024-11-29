import os

from scm.plams.core.basejob import JobError
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob, AMSResults
from scm.plams.interfaces.adfsuite.crs import CRSJob, CRSResults
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.kftools import KFFile

__all__ = ["run_crs_ams"]


def run_crs_ams(settings_ams, settings_crs, solvents, solutes=None, return_amsresults=False, **kwargs):
    """A workflow for running COSMO-RS calculations with AMS (*i.e.* DFT) COSMO surface charges.

    The workflow consists of four distinct steps:

    1. Perform gas-phase |AMSJob| calculations on the solvents and solutes (see *settings_ams*).
    2. Perform COSMO |AMSJob| calculations using the .rkf file from step 1
    3. Extract the relevant COSMO information from the results and create .coskf files.
    4. Perform a COSMO-RS calculations with the .coskf files produced in steps 2 and 3
       (see *settings_crs*).
       This calculation is conducted for all possible solvent/solute pairs,
       assuming solutes have been specified by the user.

    The ams solvation block (*ams_settings.input.solvation*) is soft updated with suitable
    settings for constructing COSMO-RS compatible surface charges
    (see :func:`.add_solvation_block`).
    No user-specified values are overwritten during this process.

    .. admonition:: Examples

        An example value for *settings_ams*:

        .. code:: python

            >>> from scm.plams import Settings

            >>> settings_ams = Settings()
            >>> settings_ams.input.adf.basis.type = 'TZ2P'
            >>> settings_ams.input.adf.xc.gga = 'BP86'
            >>> settings_ams.input.adf.scf.converge = '1.0e-06'

        An example value for *settings_crs* (`activity coefficient`_ calculation):

        .. code:: python

            >>> settings_crs = Settings()
            >>> settings_crs.input.temperature = 298.15
            >>> settings_crs.input.property._h = 'activitycoef'

        And finally the actual calculation with methanol, ethanol and propanol as solvents and
        acetic acid as solute:

        .. code:: python

            >>> solvents = [Molecule('methanol.xyz'), Molecule('ethanol.xyz'), Molecule('propanol.xyz')]
            >>> solutes = Molecule('acetic_acid.xyz')

            >>> crs_dict = run_crs_ams(settings_ams, settings_crs, solvents, solutes)
            >>> for key in crs_dict:
            >>>     print("\\nSystem:", key)
            >>>     res = crs_dict[key].get_results(settings_crs.input.property._h.upper())
            >>>     print("activity coefficients:")
            >>>     print(res["gamma"])

            System: CRSJob.methanol.acetic_acid
            activity coefficients:
            [[1.        ]
             [0.48067852]]

            System: CRSJob.ethanol.acetic_acid
            activity coefficients:
            [[1.        ]
             [0.43502696]]

            System: CRSJob.propanol.acetic_acid
            activity coefficients:
            [[1.        ]
             [0.58399355]]

    :type settings_ams: :class:`.Settings`
    :parameter settings_ams:
        A Settings instance with settings for :class:`.AMSJob` (see Examples).

    :type settings_crs: :class:`.Settings`
    :parameter settings_crs:
        A Settings instance with settings for :class:`.CRSJob` (see Examples).

    :type solvents: |Molecule| or :class:`list` [|Molecule|]
    :parameter solvents: A Molecule or list of one or more Molecules representing solvents.

    :type solutes: |Molecule| or :class:`list` [|Molecule|], optional
    :parameter solutes: An optional Molecule or list of one or more Molecules representing solutes.

    :type return_amsresults: :class:`bool`
    :parameter return_amsresults: If ``True``, return both the solvent and solute AMS results in addition to the final COSMO-RS.

    :parameter \**kwargs, optional:
        Optional keyword arguments that will be passed to all calls of :meth:`.Job.run`.
        For example, one could consider passing a custom :ref:`job_runners` or :ref:`job_manager`.

    :returns: A dictionary with the resulting COSMO-RS output.
        The `name` of each :class:`.CRSResults` instance is used as key.
        If ``return_amsresults=True``, return both the COSMO-RS and AMS solvent and solute results.
    :rtype: :class:`dict` or :class:`tuple` [:class:`dict`, :class:`list`, :class:`list`]

    .. _`activity coefficient`: ../../COSMO-RS/Properties.html#activity-coefficients-solvent-and-solute
    """

    solvents = [solvents] if isinstance(solvents, Molecule) else solvents
    solutes = [solutes] if isinstance(solutes, Molecule) else solutes

    # soft update to defaults
    _settings_ams = update_adf_defaults(settings_ams)
    _settings_ams = add_solvation_block(_settings_ams)
    # Validate arguments
    _validate_settings_ams(_settings_ams)

    # Decapitalize the "solvation" and "allpoints" keys
    if "solvation" in _settings_ams.settings.input.adf:
        _settings_ams.input.adf.solvation = _settings_ams.settings.input.adf.pop("solvation")

    # Create the COSMO surfaces for the solute
    solvent_list = [run_amsjob(mol, _settings_ams, **kwargs) for mol in solvents]
    solute_list = [run_amsjob(mol, _settings_ams, **kwargs) for mol in solutes]
    if not solutes:
        solute_list = [None]

    # Start the and return the COSMO-RS job
    crs = {}
    for solvent in solvent_list:
        for solute in solute_list:
            results = run_crsjob(solvent, settings_crs, solute=solute, **kwargs)
            crs[results.job.name] = results

    if not return_amsresults:
        return crs
    else:
        return crs, solvent_list, solute_list


def run_amsjob(mol: Molecule, s: Settings, **kwargs) -> AMSResults:
    """Run an :class:`.AMSJob` on *mol* using the settings provided in *s*."""
    name = "AMSJob." + mol.properties.name if "name" in mol.properties else "AMSJob.mol"

    # Do a geometry optimization in the gas phase
    _s = s.copy()
    _s.input.ams.soft_update({"Task": "GeometryOptimization"})
    job1 = AMSJob(molecule=mol, settings=_s, name=name + ".gas")
    del job1.settings.input.adf.solvation
    results1 = job1.run(**kwargs)
    ams_out = results1["ams.rkf"]
    adf_out = results1["adf.rkf"]

    # Construct the AMS COSMO surface
    _s.input.ams.Task = "SinglePoint"
    _s.input.ams.EngineRestart = adf_out
    _s.input.ams.LoadSystem.File = ams_out
    job2 = AMSJob(molecule=None, settings=_s, depend=[job1], name=name)
    results2 = job2.run(**kwargs)

    coskf_name = job2.name.split(".")[1]
    convert_to_coskf(results2, coskf_name)
    return results2


def convert_to_coskf(res: AMSResults, name: str):

    f = KFFile(res["adf.rkf"])
    cosmo = f.read_section("COSMO")
    # print (os.path.dirname(res.rkfpath(file='adf')), name)
    coskf = KFFile(os.path.join(os.path.dirname(res.rkfpath(file="adf")), name + ".coskf"))
    for k, v in cosmo.items():
        coskf.write("COSMO", k, v)
    res.collect()


def run_crsjob(solvent: AMSResults, s: Settings, solute: AMSResults = None, **kwargs) -> CRSResults:
    """Run an :class:`.CRSJob` on with *solvent* and, optionally, *solute* using the settings provided in *s*."""
    name = "CRSJob." + solvent.job.name.split(".")[1]

    solv_coskf = solvent.job.name.split(".")[1] + ".coskf"
    if solute is not None:
        name += "." + solute.job.name.split(".")[1]
        solute_coskf = solute.job.name.split(".")[1] + ".coskf"
        job = CRSJob(settings=s, depend=[solvent.job, solute.job], name=name)
        set_header(job.settings, solvent[solv_coskf], solute[solute_coskf])
    else:
        job = CRSJob(settings=s, depend=[solvent.job], name=name)
        set_header(job.settings, solvent[solv_coskf])

    return job.run(**kwargs)


def set_header(s: Settings, *values: str) -> None:
    """Assign *value* to the ``["_h"]`` key in *s.input.compound*."""
    s.input.compound = []
    for item in values:
        s.input.compound.append(Settings({"_h": item}))
    s.input.compound[0].frac1 = 1.0  # The first item in *values should be the solvent


def add_solvation_block(ams_settings: Settings) -> None:
    """Add the solvation block to *ams_settings*, returning a copy of the new settings.

    The solvation block (*ams_settings.input.solvation*) is soft updated
    with the following Settings:

    .. code::

        c-mat:  Exact
        charged:        method=Conj
        radii:
            Br:       2.16
            C:        2.0
            Cl:       2.05
            F:        1.72
            H:        1.3
            I:        2.32
            N:        1.83
            O:        1.72
            P:        2.13
            S:        2.16
            Si:       2.48
        scf:    Var All
        solv:   name=CRS cav0=0.0 cav1=0.0
        surf:   Delley

    """

    # Find the solvation key
    solvation = ams_settings.input.adf.find_case("solvation")
    solvation_block = ams_settings.input.adf[solvation]

    # Find all keys for within the solvation block
    keys = ("surf", "solv", "charged", "c-mat", "scf", "radii")
    surf, solv, charged, cmat, scf, radii = [solvation_block.find_case(item) for item in keys]

    # Construct the default solvation block
    solvation_block_new = {
        surf: "Delley",
        solv: "name=CRS cav0=0.0 cav1=0.0",
        charged: "method=Conj",
        cmat: "Exact",
        scf: "Var All",
        radii: {
            "H": 1.30,
            "C": 2.00,
            "N": 1.83,
            "O": 1.72,
            "F": 1.72,
            "Si": 2.48,
            "P": 2.13,
            "S": 2.16,
            "Cl": 2.05,
            "Br": 2.16,
            "I": 2.32,
        },
    }

    # Copy ams_settings and perform a soft update
    ret = ams_settings.copy()
    ret.input.adf[solvation].soft_update(solvation_block_new)
    return ret


def update_adf_defaults(ams_settings: Settings) -> None:
    """Add the COSMO-RS compound defaults to *ams_settings*, returning a copy of the new settings.

    The engine block (*ams_settings.input.adf*) is soft updated
    with the following Settings:

    .. code::

        Engine ADF
            Basis
                Type TZP
                Core Small
            End
            XC
                GGA BP86
            End
            Relativity
                Level Scalar
            End
            BeckeGrid
                Quality Good
            End
        EndEngine

    """

    # Find the solvation key
    # solvation = ams_settings.input.adf.find_case('solvation')
    # solvation_block = ams_settings.input.adf[solvation]

    adf = ams_settings.input.find_case("adf")
    adf_block = ams_settings.input[adf]

    # Find all keys for within the adf block
    keys = ("basis", "xc", "relativity", "BeckeGrid")
    basis, xc, relativity, BeckeGrid = [adf_block.find_case(item) for item in keys]

    # Construct the default solvation block
    adf_defaults_block = Settings(
        {
            basis: {"Type": "TZP", "Core": "Small"},
            xc: {"GGA": "BP86"},
            relativity: {"Level": "Scalar"},
            BeckeGrid: {"Quality": "Good"},
        }
    )
    # Copy ams_settings and perform a soft update
    ret = ams_settings.copy()
    ret.input.adf.soft_update(adf_defaults_block)
    return ret


def _validate_settings_ams(s: Settings) -> None:
    """Validate the *settings_ams* argument in :func:`run_crs_ams`."""
    solvation = s.input.adf.find_case("solvation")
    solv = s.input.adf[solvation].find_case("solv")

    if solvation not in s.input.adf:
        raise JobError("run_crs_ams: The 'solvation' key is absent" " from settings_ams.input'")

    if solv not in s.input.adf[solvation]:
        raise JobError("run_crs_ams: The 'solv' key is absent" " from settings_ams.input.solvation")

    if "name=crs" not in s.input.adf[solvation][solv].lower():
        raise JobError("run_crs_ams: The 'name=CRS' value is absent" " from settings_ams.input.solvation.solv")

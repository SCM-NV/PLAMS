ORCA
-------------------------

.. currentmodule:: scm.plams.interfaces.thirdparty.orca

Settings
~~~~~~~~~~~~~~
The ``input.main`` branch of the |Settings| class corresponds to the lines starting with the special character ``!`` in the ORCA input.
See the orcalibrary_ and the manual_ for details.
The ``input.molecule`` branch can be used to manually set the coordinate input line (e.g. to use an external file), it overwrites
any existing molecule of the |ORCAJob|. If not set, the molecule of the |ORCAJob| will be parsed using the ``xyz`` option of ORCA.

An *end* keyword is mandatory for only a subset of sections. For instance the following orca input shows the keywords *methods*
and *basis* use of end::

            job = ORCAJob(molecule=Molecule(<Path/to/molecule>), copy="input.xyz")
            job.settings.input.main = "UKS B3LYP/G SV(P) SV/J TightSCF Direct Grid3 FinalGrid4"
            job.settings.input.maxcore = "2000"
            job.settings.input.method.SpecialGridAtoms = 26
            job.settings.input.method.SpecialGridIntAcc = 7

            job.settings.input.basis.NewGTO._end = '26 "CP(PPP)"'
            job.settings.input.basis.NewAuxGTO = '26 "TZV/J" end'

            job.settings.input.molecule = "xyzfile +2 1 input.xyz" #overwrites the molecule given above

The input generated during the execution of the ORCA job is similar to: ::

            ! UKS B3LYP/G SV(P) SV/J TightSCF Direct Grid3 FinalGrid4
            %maxcore 2000
            %method SpecialGridAtoms 26
                    SpecialGridIntAcc 7
                    end
            %basis NewGTO 26 "CP(PPP)" end
                   NewAuxGTO 26 "TZV/J" end
                   end
            * xyzfile +2 1 input.xyz

Additional input files (molecules, restart files, ...) can be automatically copied to the jobs rundir by passing them to the
|ORCAJob| initilization under ``copy`` (string or list). To reduce disk usage, symlinks can be used if the filesystem permits.

.. _orcalibrary: https://sites.google.com/site/orcainputlibrary
.. _manual: https://orcaforum.kofo.mpg.de/app.php/portal

Loading jobs
~~~~~~~~~~~~~~~~~~~

Calculations done without PLAMS can be loaded using the |load_external| functionality. The |ORCAResults| class does not
support reading input files into |Settings| objects.

API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: ORCAJob(copy=None, copy_symlink=False)

.. autoclass:: ORCAResults()

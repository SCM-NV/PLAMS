Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

    import multiprocessing
    import numpy as np
    from scm.plams import Settings, Results, MultiJob, JobRunner, config, Molecule, Atom, AMSJob

Helper Classes
~~~~~~~~~~~~~~

.. code:: ipython3

    class GammaResults(Results):
    
        @staticmethod
        def get_difference(job, jobplus):
            """Calculate the difference between HOMO and IP.
            *jobplus* should be the counterpart of *job* with one less electron."""
            homo = job.results.readrkf("Properties", "HOMO", file="engine")
            IP = jobplus.results.get_energy() - job.results.get_energy()
            return IP + homo
    
        def get_J(self):
            N = GammaResults.get_difference(self.job.children[1], self.job.children[2])
            A = GammaResults.get_difference(self.job.children[0], self.job.children[1])
            return (N**2 + A**2) ** 0.5

.. code:: ipython3

    class GammaJob(MultiJob):
        _result_type = GammaResults
    
        def __init__(self, molecule, gamma, charge, spins, **kwargs):
            MultiJob.__init__(self, **kwargs)
            self.molecule = molecule
            self.charge = charge
            self.spins = spins
            self.gamma = gamma
    
        def prerun(self):
            charges = [self.charge - 1, self.charge, self.charge + 1]
            for charge, spin in zip(charges, self.spins):
                name = "{}_charge_{}".format(self.name, charge)
                name = name.replace("-", "minus")
                newjob = AMSJob(name=name, molecule=self.molecule, settings=self.settings)
                newjob.molecule.properties.charge = charge
                newjob.settings.input.adf.xc.rangesep = "gamma={:f}".format(self.gamma)
                if spin != 0:
                    newjob.settings.input.adf.unrestricted = True
                    newjob.settings.input.adf.SpinPolarization = spin
    
                self.children.append(newjob)

.. code:: ipython3

    def gamma_scan(gammas, settings, molecule, name="scan", charge=0, spins=(1, 0, 1)):
        """Calculate values of J function for given range of gammas.
    
        Arguments:
        gammas   - list of gamma values to calculate the J function for
        settings - Settings object for an ADF calculation
        molecule - Molecule object with the system of interest
        name     - base name of all the jobs
        charge   - base charge of the system of interest. The J function is going to be
                   calculated based on two systems: with charge, and charge-1
        spins    - values of spin polarization for jobs with, respectively, charge-1,
                   charge and charge +1
    
        In other words, if charge=X and spins=(a,b,c) the three resulting jobs
        are going to have the following values for charge and spin:
    
        Charge=X-1  SpinPolarization=a
        Charge=X    SpinPolarization=b
        Charge=X+1  SpinPolarization=c
    
        Returns a list of pairs (gamma, J) of the same length as the parameter *gammas*
        """
        jobs = [
            GammaJob(
                molecule=molecule, settings=settings, gamma=g, charge=charge, spins=spins, name=name + "_gamma_" + str(g)
            )
            for g in gammas
        ]
        results = [j.run() for j in jobs]
        js = [r.get_J() for r in results]
        return list(zip(gammas, js))

Configure Parallel JobRunner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set up the default jobrunner to run in parallel, with as many jobs as
there are cores.

.. code:: ipython3

    config.default_jobrunner = JobRunner(parallel=True, maxjobs=multiprocessing.cpu_count())

Settings of the ADF calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure settings object for the calculation.

.. code:: ipython3

    s = Settings()
    s.input.ams.task = "SinglePoint"
    s.input.adf.basis.type = "DZP"
    s.input.adf.basis.core = "None"
    s.input.adf.xc.gga = "PBE"
    s.input.adf.xc.xcfun = True
    s.runscript.nproc = 1

Set Up Molecule
~~~~~~~~~~~~~~~

Create a toy hydrogen dimer.

.. code:: ipython3

    mol = Molecule()
    mol.add_atom(Atom(symbol="H", coords=(0, 0, -0.3540)))
    mol.add_atom(Atom(symbol="H", coords=(0, 0, 0.3540)))

Calculate Gamma Values
~~~~~~~~~~~~~~~~~~~~~~

Perform a scan of a few values for gamma. In practice, you want to scan
a wider range and smaller step.

.. code:: ipython3

    gammas = np.around(np.arange(1.2, 1.9, 0.2), decimals=3)

.. code:: ipython3

    results = gamma_scan(gammas, s, mol)


.. parsed-literal::

    [13.08|17:22:19] JOB scan_gamma_1.2 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.4 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.6 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.8 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.2 RUNNING
    [13.08|17:22:19] Waiting for job scan_gamma_1.2 to finish
    [13.08|17:22:19] JOB scan_gamma_1.4 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.2/scan_gamma_1.2_charge_minus1 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.6 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.4/scan_gamma_1.4_charge_minus1 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.2/scan_gamma_1.2_charge_0 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.8 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.6/scan_gamma_1.6_charge_minus1 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.4/scan_gamma_1.4_charge_0 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.2/scan_gamma_1.2_charge_1 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.8/scan_gamma_1.8_charge_minus1 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.6/scan_gamma_1.6_charge_0 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.4/scan_gamma_1.4_charge_1 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.8/scan_gamma_1.8_charge_0 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.2/scan_gamma_1.2_charge_minus1 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.6/scan_gamma_1.6_charge_1 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.4/scan_gamma_1.4_charge_minus1 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.8/scan_gamma_1.8_charge_1 STARTED
    [13.08|17:22:19] JOB scan_gamma_1.2/scan_gamma_1.2_charge_0 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.6/scan_gamma_1.6_charge_minus1 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.4/scan_gamma_1.4_charge_0 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.2/scan_gamma_1.2_charge_1 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.8/scan_gamma_1.8_charge_minus1 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.6/scan_gamma_1.6_charge_0 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.4/scan_gamma_1.4_charge_1 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.8/scan_gamma_1.8_charge_0 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.6/scan_gamma_1.6_charge_1 RUNNING
    [13.08|17:22:19] JOB scan_gamma_1.8/scan_gamma_1.8_charge_1 RUNNING


.. code:: ipython3

    print("== Results ==")
    print("gamma \t J")
    for g, j in results:
        print("{:.4f} \t {:.8f}".format(g, j))
    print("Optimal gamma value: {:.4f}".format(min(results, key=lambda x: x[1])[0]))

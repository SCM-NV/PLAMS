import os, sys
from scm.plams import *

__all__ = ['RedOxPotentialCalculator']


class RedOxPotentialCalculator:
    def __init__(self, 
                 logfile = sys.stdout) -> None:

        self.logfile = logfile #file to send prints to, default None disables logging (print still works)
        self.set_default_settings() #after this function has been called it is possible to overwrite settings


    def __call__(self, *args, **kwargs) -> float:
        return self.redox_potential(*args, **kwargs)


    def set_default_settings(self) -> None:
        #### DEFAULT SETTINGS, these can be changed after creating a new RedOxPotentialCalculator object
        #DFTB GO settings
        self.DFTB_defaults = Settings()
        self.DFTB_defaults.input.ams.task   = 'GeometryOptimization'
        self.DFTB_defaults.input.DFTB
        self.DFTB_defaults.input.DFTB.Model = "GFN1-xTB" 

        #default settings for optimization and singlepoint using DFT
        self.DFT_defaults = Settings()  
        self.DFT_defaults.input.ams.task             = 'GeometryOptimization'
        self.DFT_defaults.input.adf.basis.type       = 'TZ2P'
        self.DFT_defaults.input.adf.basis.core       = 'None'
        self.DFT_defaults.input.adf.xc.hybrid        = 'B3LYP'
        self.DFT_defaults.input.adf.xc.Dispersion    = 'GRIMME3 BJDAMP'
        self.DFT_defaults.input.adf.Relativity.Level = 'None'
        self.DFT_defaults.input.adf.NumericalQuality = 'Good'
        self.DFT_defaults.input.adf.Symmetry         = 'NOSYM'
        self.DFT_defaults.input.ams.UseSymmetry      = 'No'

        #frequency calculation settings
        self.frequencies_defaults = Settings()
        self.frequencies_defaults.input.ams.properties.NormalModes   = 'Yes'
        self.frequencies_defaults.input.ams.Properties.PESPointCharacter     = 'No'
        self.frequencies_defaults.input.ams.NormalModes.ReScanFreqRange      = '-1000 0'
        self.frequencies_defaults.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = -20

        #default solvent settings for optimization
        self.COSMO_defaults = Settings()
        self.COSMO_defaults.input.adf.Solvation.Solv = "Name=Dichloromethane"

        #default settings for oxidized molecules
        self.DFT_oxidized_defaults = Settings()
        self.DFT_oxidized_defaults.input.adf.Unrestricted     = 'Yes'
        self.DFT_oxidized_defaults.input.adf.SpinPolarization = '1.0'
        self.DFT_oxidized_defaults.input.ams.System.Charge    = '1.0'
        #default settings for DFTB for oxidized molecules
        self.DFTB_oxidized_defaults = Settings()
        self.DFTB_oxidized_defaults.input.ams.System.Charge   = '1.0'

        #default settings for reduced molecules
        self.DFT_reduced_defaults = Settings()
        self.DFT_reduced_defaults.input.adf.Unrestricted      = 'Yes'
        self.DFT_reduced_defaults.input.adf.SpinPolarization  = '1.0'
        self.DFT_reduced_defaults.input.ams.System.Charge     = '-1.0'
        #default settings for DFTB for reduced molecules
        self.DFTB_reduced_defaults = Settings()
        self.DFTB_reduced_defaults.input.ams.System.Charge    = '-1.0'


    def check_termination_succes(self,
                                 result :Results) -> bool:
        term = result.readrkf('General', 'termination status', 'ams')
        if term == 'NORMAL TERMINATION':
            return True
        elif 'NORMAL TERMINATION' in term:
            return 'WARNING'
        return 


    def redox_potential(self, 
                            molecule            :Molecule, 
                            mode                :str,
                            method              :str        = 'screening',
                            name                :str        = None,
                            COSMORS_solvent_path:str        = None,
                                ) -> float:

        assert mode in ['oxidation', 'reduction']
        assert method in ['DC', 'TC-COSMO', 'TC-COSMO-RS', 'screening'], 'Argument "method" must be "DC", "TC-COSMO", "TC-COSMO-RS" or "screening"'

        #set name
        if name is None:
            name = molecule.properties.name

        self.log(f'========================================================================')
        self.log(f'Starting redox potential calculation for molecule {name}:\n')

        self.log('\nInitial coordinates:')
        self.log(molecule)
        self.log('Settings:')
        self.log(f'\tName:              {name}')
        self.log(f'\tMode:              {mode}')
        self.log(f'\tMethod:            {method}')

        if mode == 'oxidation':
            Gelectron = -0.0375
        elif mode == 'reduction':
            Gelectron = 0.0375
        #get on with the actual calculations
        if method == 'DC':
            GO_os    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='solvent')
            GO_ns    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral',  phase='solvent')

            redoxpot = GO_os['gibbs_energy'] - GO_ns['gibbs_energy'] + Gelectron
            
        elif method == 'TC-COSMO':
            GO_nv    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum')
            SP_nv_ns = self._calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', phase='solvent', frequencies=False)
            GO_ns    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='solvent', frequencies=False)
            SP_nv_nv = self._calculation_step(GO_ns['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum',  frequencies=False)

            GO_ov    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='vacuum')
            SP_ov_os = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state=mode, phase='solvent', frequencies=False)
            GO_os    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='solvent', frequencies=False)
            SP_ov_ov = self._calculation_step(GO_os['geometry'], task='SinglePoint', name=name, state=mode, phase='vacuum',  frequencies=False)

            redox_part = GO_ov['gibbs_energy'] + SP_ov_os['dG_solvation'] + (SP_ov_ov['bond_energy'] - GO_ov['bond_energy'])
            neutral_part  = GO_nv['gibbs_energy'] + SP_nv_ns['dG_solvation'] + (SP_nv_nv['bond_energy'] - GO_nv['bond_energy'])
            redoxpot = redox_part - neutral_part + Gelectron

        elif method == 'TC-COSMO-RS':
            self.COSMORS_solvent_path = COSMORS_solvent_path
            assert os.path.exists(self.COSMORS_solvent_path), f'Solvent database {self.COSMORS_solvent_path} does not exist'
            GO_nv    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum',   frequencies=False)
            SP_nv_ns = self._calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', use_COSMORS=True, frequencies=False)
            GO_ns    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', use_COSMORS=True, frequencies=False)
            SP_ns_nv = self._calculation_step(GO_ns['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum',   frequencies=False)

            GO_ov    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='vacuum', frequencies=False)
            SP_ov_os = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state=mode, use_COSMORS=True, frequencies=False)
            GO_os    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, use_COSMORS=True, frequencies=False)
            SP_os_ov = self._calculation_step(GO_os['geometry'], task='SinglePoint', name=name, state=mode, phase='vacuum',   frequencies=False)

            redox_part = GO_ov['gibbs_energy'] + SP_ov_os['dG_solvation'] + (SP_os_ov['bond_energy'] - GO_ov['bond_energy'])
            neutral_part  = GO_nv['gibbs_energy'] + SP_nv_ns['dG_solvation'] + (SP_ns_nv['bond_energy'] - GO_nv['bond_energy'])
            redoxpot = redox_part - neutral_part + Gelectron

        elif method == 'screening':
            self.COSMORS_solvent_path = COSMORS_solvent_path
            assert os.path.exists(self.COSMORS_solvent_path), f'Solvent database {self.COSMORS_solvent_path} does not exist'
            GO_nv    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum', use_dftb=True, frequencies=False)
            COSMO_nv = self._calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', use_COSMORS=True, frequencies=False)

            GO_ov    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='vacuum', use_dftb=True, frequencies=False)
            COSMO_ov = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state=mode, use_COSMORS=True, frequencies=False)

            redoxpot = COSMO_ov['gibbs_energy'] - COSMO_nv['gibbs_energy'] + Gelectron

        self.log(f"\nOxidation potential: {redoxpot:.4f} eV")
        return redoxpot

    def _get_settings(self,
                task                :str        = 'GeometryOptimization',
                state               :str        = 'neutral',
                phase               :str        = 'vacuum',
                frequencies         :bool       = True,
                use_dftb            :bool       = False,
                use_COSMORS         :bool       = False,
                    ) -> Settings:
        '''
        Method that generates settings for jobs based on 
        '''

        assert state in ['neutral', 'oxidation', 'reduction'], 'argument "state" must be "neutral", "oxidation" or "reduction"'
        assert phase in ['vacuum', 'solvent'], 'argument "phase" must be "vacuum" or "solvent"'
        assert task  in ['GeometryOptimization', 'SinglePoint'], 'argument "task" must be "GeometryOptimization", "SinglePoint" or "COSMO"'

        if use_COSMORS:
            defaults = self.DFT_defaults.copy()

            solvation_block = {
                'surf': 'Delley',
                'solv': 'name=CRS cav0=0.0 cav1=0.0',
                'charged': 'method=Conj corr',
                'c-mat': 'Exact',
                'scf': 'Var All',
                'radii': {
                    'H': 1.30,
                    'C': 2.00,
                    'N': 1.83,
                    'O': 1.72,
                    'F': 1.72,
                    'Si': 2.48,
                    'P': 2.13,
                    'S': 2.16,
                    'Cl': 2.05,
                    'Br': 2.16,
                    'I': 2.32
                }}

            defaults.input.adf.solvation = solvation_block
        else:
            if use_dftb:
                defaults = self.DFTB_defaults.copy()
            else:
                defaults = self.DFT_defaults.copy()

            #load cosmo solvent
            if phase == 'solvent' and not use_COSMORS:
                defaults.soft_update(self.COSMO_defaults)

        defaults.input.ams.task = task

        if frequencies:
            defaults.soft_update(self.frequencies_defaults)

        #handle state, if neutral the settings are already correct
        if state == 'oxidation':
            if use_dftb:
                defaults.soft_update(self.DFTB_oxidized_defaults)
            else:
                defaults.soft_update(self.DFT_oxidized_defaults)
        if state == 'reduction':
            if use_dftb:
                defaults.soft_update(self.DFTB_reduced_defaults)
            else:
                defaults.soft_update(self.DFT_reduced_defaults)

        return defaults


    def _COSMORS_property(self, 
                solvent_path        :str, 
                solute_path         :str,
                name                :str,
                temperature         :float      = 298.15
                    ) -> float:
        """This method runs a COSMORS property job to obtain the activity coefficient
        which will also calculate G solute which we need to calculate the oxidation 
        potential
        """

        defaults = Settings()
        defaults.input.property._h = 'ACTIVITYCOEF'
        compounds = [Settings(), Settings()]
        compounds[0]._h = solvent_path
        compounds[1]._h = solute_path
        compounds[0].frac1 = 1
        compounds[1].frac1 = 0

        defaults.input.temperature = str(temperature)
        defaults.input.compound = compounds

        res = CRSJob(settings=defaults, name=name).run().get_results()
        if res:
            return float(Units.convert(res["G solute"][1], 'kcal/mol', 'hartree'))
        else:
            return False


    def _calculation_step(self, 
                molecule            :Molecule,
                task                :str        = 'GeometryOptimization',
                state               :str        = 'neutral',
                phase               :str        = 'vacuum',
                frequencies         :bool       = True,
                use_dftb            :bool       = False,
                use_COSMORS         :bool       = False,
                name                :str        = None,
                    ) -> dict:
        """ Method used to optimize the geometry of molecule using DFT (by default B3LYP)
        Other settings may be supplied using settings which will be soft-updated using self.DFT_defaults
        State specifies whether the molecule is neutral, oxidised or reduced
        Phase specifies whether the system is in vacuum or solvated
        if use_dftb, the system will be optimised using DFTB (by default GFN1-xTB) instead of DFT
        """

        assert state in ['neutral', 'oxidation', 'reduction'], 'argument "state" must be "neutral", "oxidation", "reduction"'
        assert phase in ['vacuum', 'solvent'], 'argument "phase" must be "vacuum" or "solvent"'
        assert task  in ['GeometryOptimization', 'SinglePoint'], 'argument "task" must be "GeometryOptimization" or "SinglePoint"'

        if use_COSMORS: phase = 'solvent'

        settings = self._get_settings(task=task, 
                                      use_dftb=use_dftb, 
                                      use_COSMORS=use_COSMORS,
                                      state=state, 
                                      phase=phase, 
                                      frequencies=frequencies)

        #summarize job in one string
        task_abbrev = {"GeometryOptimization":"GO", "SinglePoint":"SP"}[task]
        job_desc = f'{task_abbrev}_{state}_{phase}'
        if use_COSMORS:
            job_desc += '_COSMO-RS'
        if use_dftb:
            job_desc += '_DFTB'

        self.log(f'\nStarting calculation {name + "_" + job_desc}')
        self.log(f'\ttask                 = {task}')
        self.log(f'\tuse_dftb             = {use_dftb}')
        self.log(f'\tuse_COSMORS          = {use_COSMORS}')
        self.log(f'\tfrequencies          = {frequencies}')
        self.log(f'\tstate                = {state}')
        self.log(f'\tphase                = {phase}')

        #run the job
        job = AMSJob(molecule   = molecule,
                     settings   = settings,
                     name       = name + '_' + job_desc)
        res = job.run()

        result_dict = {}
        #pull out results
        if self.check_termination_succes(res):
            self.log(f'\tSuccessfull          = {self.check_termination_succes(res)}') #True or WARNING
            #set some default values
            bond_energy = None 
            gibbs_energy = None

            #If we are doing COSMO calculations then we need to run an additional job to obtain the activity coefficient
            #when calculating the activity coefficient, the G solute is also calculated.
            if use_COSMORS:
                resfile = KFFile(res['adf.rkf'])
                cosmo_data = resfile.read_section('COSMO')
                coskf = KFFile(os.path.join(job.path, job.name + '.coskf'))
                for k, v in cosmo_data.items():
                    coskf.write('COSMO', k, v)
                res.collect()
                bond_energy = res.readrkf('AMSResults', 'Energy', 'adf')
                gibbs_energy = self._COSMORS_property(self.COSMORS_solvent_path, os.path.join(job.path, job.name + '.coskf'), job.name + '_ACTIVITYCOEF')
            #if we dont use COSMO-RS we can just extract the Gibbs and bonding energies from the regular job
            else:
                if use_dftb:
                    bond_energy = res.readrkf('AMSResults', 'Energy', 'dftb')
                    if frequencies: 
                        gibbs_energy = res.readrkf('Thermodynamics', 'Gibbs free Energy', 'dftb')
                else:
                    bond_energy = res.readrkf('Energy', 'Bond Energy', 'adf')
                    if frequencies: 
                        gibbs_energy = res.readrkf('Thermodynamics', 'Gibbs free Energy', 'adf')
            
            self.log(f'\tResults:')
            if not bond_energy is None: 
                result_dict['bond_energy'] = Units.convert(bond_energy, 'hartree', 'eV')
                self.log(f'\t\tBond Energy  = {result_dict["bond_energy"]:.4f} eV')
            if not gibbs_energy is None: 
                result_dict['gibbs_energy'] = Units.convert(gibbs_energy, 'hartree', 'eV')
                self.log(f'\t\tGibbs Energy = {result_dict["gibbs_energy"]:.4f} eV')

            #extract also optimised molecule
            if task == 'GeometryOptimization':
                result_dict['geometry'] = res.get_main_molecule()

            #and if the phase is solvent we also need the solvation gibbs energy change
            if phase == 'solvent':
                dG_solvation = res.readrkf('Energy','Solvation Energy (el)','adf') + res.readrkf('Energy','Solvation Energy (cd)','adf')
                result_dict['dG_solvation'] = Units.convert(dG_solvation, 'hartree', 'eV')
                self.log(f'\t\tdG_solvation = {result_dict["dG_solvation"]:.4f} eV')

        else:
            self.log(f'\tSuccessfull          = False')

        return result_dict


    def log(self, line, newline=True):
        if not self.logfile is None: 
            self.logfile.write(str(line) + '\n'*newline)
            self.logfile.flush()
        print(line)



if __name__ == '__main__':
    mode = 'oxidation'
    job_dir = './Test'
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    COSMORS_solvent_path = os.path.abspath('Dichloromethane.coskf')

    results = {}
    for mol_file in ['./molecules/NDI44.xyz']:
        mol_name = os.path.basename(mol_file).split('.')[0]
        results[mol_name] = {}
        for method in ['screening', 'TC-COSMO', 'TC-COSMO-RS', 'DC']:
            job_name = None #if set to None, a name will be generated
            
            if job_name is None:
                job_name = mol_name + '_' + method

            #calculation part
            init(path=job_dir, folder=job_name)

            workdir = config.default_jobmanager.workdir
            logfile = open(f'{workdir}/{job_name}_python.log', 'w')

            RedOxPotCalc = RedOxPotentialCalculator(logfile=logfile)
            mol = Molecule(mol_file)
            redoxpot = RedOxPotCalc(mol, mode, method=method, COSMORS_solvent_path=COSMORS_solvent_path)
            results[mol_name][method] = redoxpot
            finish()

    print('RedOx Potentials:')
    name_len = max(len('System'), max(len(n) for n in results))
    methods = list(results[list(results.keys())[0]])
    method_lens = [max(9, len(m)) for m in methods]
    print(f'{"System".ljust(name_len)} | {" | ".join([m.ljust(l) for m, l in zip(methods, method_lens)])}')
    for name, res in results.items():
        s = f'{name.ljust(name_len)} | {" | ".join([(str(round(res[m],3)) + " eV").rjust(l) for m, l in zip(methods, method_lens)])}'
        print(s)

    print('\nCalculations coomplete!\a')
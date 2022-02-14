import os, sys
from scm.plams import *

__all__ = ['OxidationPotentialCalculator']

class OxidationPotentialCalculator:
    def __init__(self, 
                 logfile = sys.stdout) -> None:

        self.logfile = logfile #file to send prints to, default None disables logging (print still works)
        self.set_default_settings() #after this function has been called it is possible to overwrite settings


    def __call__(self, *args, **kwargs) -> float:
        return self.oxidation_potential(*args, **kwargs)


    def set_default_settings(self) -> None:
        #### DEFAULT SETTINGS these can be changed after creating a new OxidationPotentialCalculator object
        #default settings for pre-optimization using DFTB
        self.pre_optimize_defaults = Settings()
        self.pre_optimize_defaults.input.ams.task = 'GeometryOptimization' 
        self.pre_optimize_defaults.input.ams.Properties.NormalModes           = 'Yes'
        self.pre_optimize_defaults.input.ams.Properties.PESPointCharacter     = 'Yes'
        self.pre_optimize_defaults.input.ams.NormalModes.ReScanFreqRange      = '-1000 0'
        self.pre_optimize_defaults.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = -20
        self.pre_optimize_defaults.input.DFTB
        self.pre_optimize_defaults.input.DFTB.Model           = "DFTB3" 
        self.pre_optimize_defaults.input.DFTB.ResourcesDir    = 'DFTB.org/3ob-3-1'

        #DFTB GO settings
        self.DFTB_defaults = Settings()
        self.DFTB_defaults.input.ams.task       = 'GeometryOptimization'
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
        self.DFT_defaults.input.adf.NumericalQuality         = 'Good'
        self.DFT_defaults.input.adf.RIHartreeFock.UseMe      = 'Yes'
        self.DFT_defaults.input.adf.RIHartreeFock.Quality    = 'Normal'  
        self.DFT_defaults.input.adf.Symmetry                 = 'NOSYM'
        self.DFT_defaults.input.ams.UseSymmetry              = 'No'

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
        self.DFTB_oxidized_defaults.input.ams.System.Charge    = '1.0'


    def set_paths(self, 
                  job_dir:str = os.getcwd()
                  ) -> None:

        #default paths
        self.job_dir = job_dir
        if not os.path.exists(job_dir):                 os.makedirs(job_dir)
        self.log(f'Job directory:                             {job_dir}')

        self.geo_preopt_dir     = os.path.join(job_dir,'Geometries', 'Preoptimization', 'Preopt_Geo')
        if not os.path.exists(self.geo_preopt_dir):     os.makedirs(self.geo_preopt_dir)
        self.log(f'Preoptimization geometry directory:        {self.geo_preopt_dir}')

        self.not_geo_preopt_dir     = os.path.join(job_dir,'Geometries', 'Preoptimization', 'Not_Preopt_Geo')
        if not os.path.exists(self.not_geo_preopt_dir): os.makedirs(self.not_geo_preopt_dir)
        self.log(f'Failed preoptimization geometry directory: {self.not_geo_preopt_dir}')

        self.geo_opt_dir    = os.path.join(job_dir,'Geometries', "Optimized_Geo")
        if not os.path.exists(self.geo_opt_dir):        os.makedirs(self.geo_opt_dir)
        self.log(f'DFT optimized geometry directory:          {self.geo_opt_dir}')

        self.not_geo_opt_dir = os.path.join(job_dir,'Geometries', "Not_Optimized_Geo")
        if not os.path.exists(self.not_geo_opt_dir):    os.makedirs(self.not_geo_opt_dir)
        self.log(f'Failed DFT optimized geometry directory:   {self.not_geo_opt_dir}')


    def check_termination_succes(self,
                                 result :Results) -> bool:
        term = result.readrkf('General', 'termination status', 'ams')
        if term == 'NORMAL TERMINATION':
            return True
        elif 'NORMAL TERMINATION' in term:
            return 'WARNING'
        return 

    
    def oxidation_potential(self, 
                            molecule            :Molecule, 
                            method              :str        = 'screening',
                            name                :str        = None,
                            pre_optimize        :bool       = True,
                            job_dir             :str        = os.getcwd(),
                            COSMORS_solvent_path:str        = None,
                                ) -> float:
        

        assert method in ['DC', 'TC-COSMO', 'TC-COSMO-RS', 'screening'], 'Argument "method" must be "DC", "TC-COSMO", "TC-COSMO-RS" or "screening"'

        #set name
        if name is None:
            name = molecule.properties.name

        self.log(f'========================================================================')
        self.log(f'Starting oxidation potential calculation for molecule {name}:\n')

        #set paths for the jobs:
        self.log('Setting paths ...')
        self.set_paths(job_dir)

        self.log('\nInitial coordinates:')
        self.log(molecule)
        self.log('Settings:')
        self.log(f'\tName:              {name}')
        self.log(f'\tMethod:            {method}')
        self.log(f'\tPreoptimize:       {pre_optimize}')

        if pre_optimize:
            self.log('\nPre-optimizing the molecule ...', newline=False)
            molecule = self._pre_optimize(molecule, name=name)
            self.log(' done!')

        #get on with the actual calculations
        if method == 'DC':
            GO_os    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='oxidized', phase='solvent')
            GO_ns    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral',  phase='solvent')

            oxpot = GO_os['gibbs_energy'] - GO_ns['gibbs_energy'] + 0.0375
            
        elif method == 'TC-COSMO':
            GO_nv    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum')
            SP_ns_nv = self._calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', phase='solvent')
            GO_ns    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='solvent')
            SP_nv_ns = self._calculation_step(GO_ns['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum')

            GO_ov    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='oxidized', phase='vacuum')
            SP_os_ov = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state='oxidized', phase='solvent')
            SP_nv_ov = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum')
            GO_os    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='oxidized', phase='solvent')
            SP_nv_os = self._calculation_step(GO_os['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum')

            oxidized_part = GO_nv['gibbs_energy'] + SP_ns_nv['gibbs_energy'] + (GO_nv['bond_energy'] - SP_nv_ns['bond_energy'])
            neutral_part  = GO_ov['gibbs_energy'] + SP_os_ov['gibbs_energy'] + (SP_nv_ov['bond_energy'] - SP_nv_os['bond_energy'])
            oxpot = oxidized_part - neutral_part + 0.0375

        elif method == 'TC-COSMO-RS':
            GO_nv    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum')
            SP_ns_nv = self._calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', phase='solvent')
            GO_ns    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='solvent')
            SP_nv_ns = self._calculation_step(GO_ns['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum')

            GO_ov    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='oxidized', phase='vacuum')
            SP_os_ov = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state='oxidized', phase='solvent')
            SP_nv_ov = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum')
            GO_os    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='oxidized', phase='solvent')
            SP_nv_os = self._calculation_step(GO_os['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum')

            oxidized_part = GO_nv['gibbs_energy'] + SP_ns_nv['gibbs_energy'] + (GO_nv['bond_energy'] - SP_nv_os['bond_energy'])
            neutral_part  = GO_ov['gibbs_energy'] + SP_os_ov['gibbs_energy'] + (SP_nv_ov['bond_energy'] - SP_nv_os['bond_energy'])
            oxpot = oxidized_part - neutral_part + 0.0375

        elif method == 'screening':
            self.COSMORS_solvent_path = COSMORS_solvent_path
            assert os.path.exists(self.COSMORS_solvent_path), f'Solvent database {self.COSMORS_solvent_path} does not exist'
            GO_nv    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum', use_dftb=True)
            SP_nv_nv = self._calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum', frequencies=False)
            COSMO_nv = self._calculation_step(GO_nv['geometry'], task='COSMO', name=name, state='neutral')

            GO_ov    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='oxidized', phase='vacuum', use_dftb=True)
            SP_ov_ov = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state='oxidized', phase='vacuum', frequencies=False)
            COSMO_ov = self._calculation_step(GO_ov['geometry'], task='COSMO', name=name, state='oxidized')

            oxpot = COSMO_ov['gibbs_energy'] - COSMO_nv['gibbs_energy'] + 0.0375

        self.log(f"\nOxidation potential: {oxpot:.4f} eV")
        return oxpot


    def _pre_optimize(self,
                  molecule  :Molecule,
                  name      :str        = None,
                     ) -> Molecule:
        """Method used to pre-optimize the molecule using DFTB3
        """
        job = AMSJob(molecule   = molecule,
                     settings   = self.pre_optimize_defaults,
                     name       = name + '_preoptimization')
        res = job.run()

        #if the optimization is succesful we return the new molecule
        if self.check_termination_succes(res):
            molecule = res.get_main_molecule()
            molecule.write(os.path.join(self.geo_preopt_dir, name + '.xyz'))
        else:
            molecule.write(os.path.join(self.not_geo_preopt_dir, name + '.xyz'))
        
        return molecule


    def _get_settings(self,
                task                :str        = 'GeometryOptimization',
                use_dftb            :bool       = False,
                state               :str        = 'neutral',
                phase               :str        = 'vacuum',
                frequencies         :bool        = True,
                    ) -> Settings:
        '''
        Method that generates settings for jobs based on 
        '''

        assert state in ['neutral', 'oxidized'], 'argument "state" must be "neutral" or "oxidized"'
        assert phase in ['vacuum', 'solvent'], 'argument "phase" must be "vacuum" or "solvent"'
        assert task  in ['GeometryOptimization', 'SinglePoint', 'COSMO'], 'argument "task" must be "GeometryOptimization", "SinglePoint" or "COSMO"'

        if task == 'COSMO':
            defaults = self.DFT_defaults.copy()
            defaults.input.ams.task = 'SinglePoint'

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

            defaults.input.ams.task = task
    
            #load cosmo solvent
            if phase == 'solvent':
                defaults.soft_update(self.COSMO_defaults)

            if frequencies:
                defaults.soft_update(self.frequencies_defaults)

        #handle state, if neutral the settings are already correct
        if state == 'oxidized':
            if use_dftb:
                defaults.soft_update(self.DFTB_oxidized_defaults)
            else:
                defaults.soft_update(self.DFT_oxidized_defaults)

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
                use_dftb            :bool       = False,
                state               :str        = 'neutral',
                phase               :str        = 'vacuum',
                frequencies         :bool       = True,
                name                :str        = None,
                    ) -> dict:
        """ Method used to optimize the geometry of molecule using DFT (by default B3LYP)
        Other settings may be supplied using settings which will be soft-updated using self.DFT_defaults
        State specifies whether the molecule is neutral or oxidised
        Phase specifies whether the system is in vacuum or solvated
        if use_dftb, the system will be optimised using DFTB (by default GFN1-xTB) instead of DFT
        """

        assert state in ['neutral', 'oxidized'], 'argument "state" must be "neutral" or "oxidized"'
        assert phase in ['vacuum', 'solvent'], 'argument "phase" must be "vacuum" or "solvent"'
        assert task  in ['GeometryOptimization', 'SinglePoint', 'COSMO'], 'argument "task" must be "GeometryOptimization", "SinglePoint" or "COSMO"'

        if task == 'COSMO': phase = 'solvent'

        settings = self._get_settings(task=task, 
                                      use_dftb=use_dftb, 
                                      state=state, 
                                      phase=phase, 
                                      frequencies=frequencies)

        #summarize job in one string
        task_abbrev = {"GeometryOptimization":"GO", "SinglePoint":"SP", "COSMO":"COSMO"}[task]
        job_desc = f'{task_abbrev}_{state}_{phase}'
        if use_dftb:
            job_desc += '_DFTB'

        self.log(f'\nStarting calculation {name + "_" + job_desc}')
        self.log(f'\ttask                 = {task}')
        self.log(f'\tuse_dftb             = {use_dftb}')
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
            if task == 'COSMO':
                resfile = KFFile(res['adf.rkf'])
                cosmo_data = resfile.read_section('COSMO')
                coskf = KFFile(os.path.join(job.path, job.name + '.coskf'))
                for k, v in cosmo_data.items():
                    coskf.write('COSMO', k, v)
                res.collect()
                bond_energy = res.readrkf('AMSResults', 'Energy', 'adf')
                gibbs_energy = self._COSMORS_property(self.COSMORS_solvent_path, os.path.join(job.path, job.name + '.coskf'), job.name + '_ACTIVITYCOEF')
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

            if task == 'GeometryOptimization':
                opt_mol = res.get_main_molecule()
                opt_mol.write(os.path.join(self.geo_opt_dir, name + '_' + job_desc + '.xyz'))
                result_dict['geometry'] = res.get_main_molecule()

        else:
            self.log(f'\tSuccessfull          = False')
            if task == 'GeometryOptimization':
                molecule.write(os.path.join(self.not_geo_opt_dir, name + '_' + job_desc + '.xyz'))

        return result_dict


    def log(self, line, newline=True):
        if not self.logfile is None: 
            self.logfile.write(str(line) + '\n'*newline)
            self.logfile.flush()
        print(line)



if __name__ == '__main__':
    job_dir = './Screening'
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    COSMORS_solvent_path = os.path.abspath('Dichloromethane.coskf')

    results = {}
    for mol_file in ["./molecules/NDI.xyz", "./molecules/NDI44.xyz", "./molecules/NDI55.xyz", "./molecules/NDI54.xyz", "./molecules/PDI.xyz"]:
        for method in ['screening']:
            job_name = None #if set to None, a name will be generated
        
            if job_name is None:
                job_name = os.path.basename(mol_file).split('.')[0] + '_' + method

            #calculation part
            init(path=job_dir, folder=job_name)

            workdir = config.default_jobmanager.workdir
            logfile = open(f'{workdir}/{job_name}_python.log', 'w')

            ox_potential_calc = OxidationPotentialCalculator(logfile=logfile)
            mol = Molecule(mol_file)
            oxpot = ox_potential_calc.oxidation_potential(mol, job_dir=workdir, method=method, COSMORS_solvent_path=COSMORS_solvent_path)
            results[job_name] = oxpot
            finish()

    name_len = max(len(n) for n in results)
    print(f'\n{"System".ljust(name_len)} | Oxidation potential')
    for n, o in results.items():
        print(f'{n:{name_len}} | {o:.3f} eV')

    print('\nCalculations coomplete!\a')

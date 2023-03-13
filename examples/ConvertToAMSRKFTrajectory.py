#!/usr/bin/env amspython
from scm.plams import *
import os

def main():
    # the first argument needs to be a file readable by ASE
    convert_to_ams_rkf_with_bond_guessing('OUTCAR', 'converted-to-ams.rkf')
    #convert_to_ams_rkf_with_bond_guessing('somefile.xyz', 'converted-to-ams.rkf')

def convert_to_ams_rkf_with_bond_guessing(filename, outfile='out.rkf', task='moleculardynamics', timestep=0.5):
    temp_traj = 'out.traj'
    file_to_traj(filename, temp_traj)
    traj_to_rkf(temp_traj, outfile, task=task, timestep=timestep)

    init()
    #config.log.stdout = 0
    #config.erase_workdir = True   # to remove workdir, only use this if you're not already inside another PLAMS workflow

    s = Settings()
    s.input.ams.task = 'replay'
    s.input.lennardjones
    s.input.ams.replay.file = os.path.abspath(outfile)
    s.input.ams.properties.molecules = 'yes'
    s.input.ams.properties.bondorders = 'yes'
    s.runscript.nproc = 1
    job = AMSJob(settings=s, name='rep')
    job.run()
    job.results.wait()
    cpkf = os.path.expandvars('$AMSBIN/cpkf')
    os.system(f'sh "{cpkf}" "{job.results.rkfpath()}" "{outfile}" History Molecules')
    delete_job(job)
    finish()

    os.remove(temp_traj)

if __name__ == '__main__':
    main()


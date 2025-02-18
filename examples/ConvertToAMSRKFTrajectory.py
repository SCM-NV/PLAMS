#!/usr/bin/env amspython
import os

from scm.plams import *


def main():
    # the first argument needs to be a file readable by ASE
    original_trajectory = os.path.expandvars(
        "$AMSHOME/scripting/scm/params/tests/fixtures/results_importer/vasp/OUTCAR"
    )
    convert_to_ams_rkf_with_bond_guessing(original_trajectory, "converted-to-ams.rkf")
    # convert_to_ams_rkf_with_bond_guessing('somefile.xyz', 'converted-to-ams.rkf')


def convert_to_ams_rkf_with_bond_guessing(filename, outfile="out.rkf", task="moleculardynamics", timestep=0.5):
    # this line is not required in AMS2025+
    init()

    temp_traj = "out.traj"
    file_to_traj(filename, temp_traj)
    traj_to_rkf(temp_traj, outfile, task=task, timestep=timestep)

    # config.log.stdout = 0
    # config.erase_workdir = True   # to remove workdir, only use this if you're not already inside another PLAMS workflow

    s = Settings()
    s.input.ams.task = "replay"
    s.input.lennardjones
    s.input.ams.replay.file = os.path.abspath(outfile)
    s.input.ams.properties.molecules = "yes"
    s.input.ams.properties.bondorders = "yes"
    s.runscript.nproc = 1
    job = AMSJob(settings=s, name="rep")
    job.run()
    job.results.wait()
    cpkf = os.path.expandvars("$AMSBIN/cpkf")
    os.system(f'sh "{cpkf}" "{job.results.rkfpath()}" "{outfile}" History Molecules')
    delete_job(job)

    os.remove(temp_traj)


if __name__ == "__main__":
    main()

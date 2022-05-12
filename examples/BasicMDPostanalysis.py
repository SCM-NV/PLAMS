#!/usr/bin/env plams
from scm.plams import *
import numpy as np
import os
import matplotlib.pyplot as plt

def run_md():
    mol = packmol_liquid(from_smiles('O'), n_molecules=16, density=1.0)
    s = Settings()
    s.input.ams.Task = 'MolecularDynamics'
    s.input.ReaxFF.ForceField = 'Water2017.ff'
    s.input.ams.MolecularDynamics.CalcPressure = 'Yes'
    s.input.ams.MolecularDynamics.InitialVelocities.Temperature = 300
    s.input.ams.MolecularDynamics.Trajectory.SamplingFreq = 1
    s.input.ams.MolecularDynamics.TimeStep = 0.5
    s.input.ams.MolecularDynamics.NSteps = 2000
    s.runscript.nproc = 1
    os.environ['OMP_NUM_THREADS'] = '1'
    job = AMSJob(settings=s, molecule=mol, name='md')
    job.run()
    return job

def plot_results(results):
    plt.clf()
    times, vacf = results.get_velocity_acf(start_fs=0, max_dt_fs=250, normalize=False)
    normalized_vacf = vacf / vacf[0]
    plt.plot(times, normalized_vacf)
    plt.xlabel("Time (fs)")
    plt.ylabel("Velocity autocorrelation function")
    plt.title("Velocity autocorrelation function")
    plt.savefig("plams_vacf.png")
    A = np.stack((times, normalized_vacf), axis=1)
    np.savetxt("plams_vacf.txt", A, header="Time(fs) VACF")

    plt.clf()
    t_D, D = results.get_diffusion_coefficient_from_velocity_acf(times, vacf)
    plt.plot(t_D, D)
    plt.xlabel("Time (fs)")
    plt.ylabel("D (m²s⁻¹)")
    plt.title("Diffusion coefficient")
    plt.savefig("plams_vacf_D.png")
    A = np.stack((t_D, D), axis=1)
    np.savetxt("plams_vacf_D.txt", A, header="time(fs) D(m^2*s^-1)")

    plt.clf()
    freq, intensities = results.get_power_spectrum(times, vacf, number_of_points=1000)
    plt.plot(freq, intensities)
    plt.xlabel("Frequency (cm⁻¹)")
    plt.ylabel("Power spectrum (arbitrary units)")
    plt.title("Power spectrum")
    plt.savefig("plams_power_spectrum.png")
    A = np.stack((freq, intensities), axis=1)
    np.savetxt("plams_power_spectrum.txt", A, header="Frequency(cm^-1) PowerSpectrum")

    plt.clf()
    t, viscosity = results.get_green_kubo_viscosity(start_fs=0, max_dt_fs=250) # do not do this for NPT simulations
    plt.plot(t, viscosity)
    plt.xlabel("Time (fs)")
    plt.ylabel("Viscosity (mPa*s)")
    plt.title("Viscosity")
    plt.savefig("plams_green_kubo_viscosity.png")
    A = np.stack((t, viscosity), axis=1)
    np.savetxt("plams_green_kubo_viscosity.txt", A, header="Time(fs) Viscosity(mPa*s)")

def main():
    job = run_md()
    # alternatively:
    # job = AMSJob.load_external('/path/to/ams.rkf')
    plot_results(job.results)

if __name__ == '__main__':
    main()


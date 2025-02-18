from scm.plams import Molecule, Settings, AMSNVTJob, AMSNVEJob, AMSAnalysisJob, init

text = """6

O          0.0000000000        0.0000000000        0.0000000000
H          0.5580000000        0.5980000000       -0.5130000000
H          0.0280000000       -0.8150000000       -0.5120000000
O         -0.0660000000        0.0230000000        2.7240000000
H          0.0000000000        0.0000000000        1.7490000000
H         -0.9560000000        0.3570000000        2.8380000000
"""


def main():
    """
    The main script
    """
    # this line is not required in AMS2025+
    init()

    # Write the XYZ file
    xyzfile = open("2h2o.xyz", "w")
    xyzfile.write(text)
    xyzfile.close()

    # Set up and run the equilibration at 300K
    mol = Molecule("2h2o.xyz")
    settings = Settings()
    settings.input.DFTB.Model = "GFN1-xTB"
    job = AMSNVTJob(molecule=mol, settings=settings, nsteps=10000, timestep=0.5, thermostat="Berendsen")
    results = job.run()

    # Run the production run, without thermostat
    job = AMSNVEJob.restart_from(job, binlog_dipolemoment=True)
    results = job.run()

    # Set up and run the analysis
    ansett = Settings()
    ansett.input.TrajectoryInfo.Trajectory.KFFilename = results.rkfpath()
    ansett.input.Task = "AutoCorrelation"
    ansett.input.AutoCorrelation.Property = "DipoleMomentFromBinLog"
    ansett.input.AutoCorrelation.WritePropertyToKF = "Yes"
    ansett.input.AutoCorrelation.UseTimeDerivative.Enabled = "Yes"
    ansett.input.AutoCorrelation.WritePropertyToKF = "Yes"

    anjob = AMSAnalysisJob(settings=ansett)
    anresults = anjob.run()

    # Write the analysis plots
    plots = anresults.get_all_plots()
    for xy in plots:
        xy.write("%s" % (xy.section + ".txt"))


if __name__ == "__main__":
    main()

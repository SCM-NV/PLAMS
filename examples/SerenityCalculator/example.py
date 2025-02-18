import os
from scm.plams import from_smiles, Settings, AMSJob, log

molecule = from_smiles("[HH]")

# first do a pesscan with pure AMS
s = Settings()
s.input.ams.Task = "PESScan"
s.input.ams.PESScan.ScanCoordinate.nPoints = 5
s.input.ams.PESScan.ScanCoordinate.Distance = "1 2 0.65 0.8"
s.input.ForceField.Type = "UFF"
psjob = AMSJob(settings=s, molecule=molecule, name="ams_pesscan")
psjob.run()

# then run a replay with the serenity calculator

s = Settings()
s.input.ams.Task = "Replay"
s.input.ams.Replay.File = psjob.results.rkfpath()
s.input.ASE.File = os.path.abspath("my_calculator.py")

job = AMSJob(settings=s, molecule=molecule, name="serenity_replay")
job.run()

errormsg = job.get_errormsg()
if errormsg:
    log(errormsg)

replayresults = job.results.get_pesscan_results()
log(replayresults)

""" Deprecated, do not use """

from scm.plams import Settings


### ==== DFT SETTINGS ==== ###
def DFT():
    s = Settings()
    s.input.ams.task = "GeometryOptimization"
    s.input.adf.basis.type = "TZ2P"
    s.input.adf.basis.core = "None"
    s.input.adf.xc.hybrid = "B3LYP"
    s.input.adf.xc.Dispersion = "GRIMME3 BJDAMP"
    s.input.adf.Relativity.Level = "None"
    s.input.adf.NumericalQuality = "Good"
    s.input.adf.Symmetry = "NOSYM"
    s.input.ams.UseSymmetry = "No"
    s.input.adf.Unrestricted = "No"
    s.input.adf.SpinPolarization = 0
    s.input.ams.System.Charge = 0
    return s


### ==== DFTB SETTINGS ==== ###
def DFTB():
    s = Settings()
    s.input.ams.task = "GeometryOptimization"
    s.input.DFTB
    s.input.DFTB.Model = "GFN1-xTB"
    s.input.ams.System.Charge = 0
    return s


### ==== FREQ SETTINGS ==== ###
def frequencies():
    s = Settings()
    s.input.ams.properties.NormalModes = "Yes"
    s.input.ams.Properties.PESPointCharacter = "No"
    s.input.ams.NormalModes.ReScanFreqRange = "-1000 0"
    s.input.ams.PESPointCharacter.NegativeEigenvalueTolerance = -0.001
    return s

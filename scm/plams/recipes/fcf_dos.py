# from .scmjob import SCMJob, SCMResults
import math

import numpy
from scm.plams.tools.kftools import KFFile

__all__ = ["FCFDOS"]


class FCFDOS:
    """
    A class for calculating the convolution of two FCF spectra
    """

    # Conversion factor to cm-1 == ( me^2 * c * alpha^2 ) / ( 100 * 2pi * amu * hbar )
    G2F = 120.399494933
    # J2Ha = 1.0 / (scipy.constants.m_e * scipy.constants.speed_of_light**2 * 0.0072973525693**2)
    J2Ha = 2.293712278400752e17

    def __init__(self, absrkf, emirkf, absele=0.0, emiele=0.0):
        """
        absrkf : KFFile from an FCF absorption calculation for the acceptor, as name, path, or KFFile object
        emirkf : KFFile from an FCF emission calculation for the donor, as name, path, or KFFile object
        absele : Acceptor absorption electronic energy cm-1
        emiele : Donor emission electronic energy in cm-1
        """
        if isinstance(absrkf, str):
            self.absrkf = KFFile(absrkf)
        else:
            self.absrkf = absrkf
        if isinstance(emirkf, str):
            self.emirkf = KFFile(emirkf)
        else:
            self.emirkf = emirkf
        self.absele = absele
        self.emiele = emiele
        self.spc = None

    def _getstick(self, source):
        spclen = source.read("Fcf", "nspectrum")
        rawspc = source.read("Fcf", "spectrum")
        stick = numpy.reshape(numpy.array(rawspc), (2, spclen)).transpose()
        # Reorder spectrum if decreasing
        if stick[0, 0] > stick[-1, 0]:
            for ix in range(numpy.size(stick, 0) // 2):
                XX = numpy.copy(stick[ix, :])
                stick[ix, :] = stick[-ix - 1, :]
                stick[-ix - 1, :] = XX
        # Get frequencies in cm-1
        frq1 = numpy.array(source.read("Fcf", "gamma1")) * self.G2F
        frq2 = numpy.array(source.read("Fcf", "gamma2")) * self.G2F
        # Calculate energy in Joules
        # factor = scipy.constants.pi * scipy.constants.hbar * scipy.constants.speed_of_light * 100.0
        factor = 9.932229285744643e-24
        viben1 = sum(frq1) * factor
        viben2 = sum(frq2) * factor
        # Convert to Hartree
        viben1 = viben1 * self.J2Ha
        viben2 = viben2 * self.J2Ha
        # Add ZPE and electronic energy
        stick[:, 0] = stick[:, 0] + viben2 - viben1 + self.absele + self.emiele
        return stick

    def _spcinit(self):
        # Get absorption and emission spectra
        absspc = self._getstick(self.absrkf)
        emispc = self._getstick(self.emirkf)
        # Find spectrum bounds
        absmin = numpy.amin(absspc[:, 0])
        absmax = numpy.amax(absspc[:, 0])
        emimin = numpy.amin(emispc[:, 0])
        emimax = numpy.amax(emispc[:, 0])
        spcmin = min(absmin, emimin)
        spcmax = max(absmax, emimax)
        # Find spectrum length and grain
        absdlt = abs(absspc[0, 0] - absspc[1, 0])
        emidlt = abs(absspc[0, 0] - absspc[1, 0])
        spcgrn = min(absdlt, emidlt)
        spclen = math.floor((spcmax - spcmin) / spcgrn) + 1
        # Initialize spectrum
        self.spc = self.newspc(spcmin, spcmax, spclen)
        return None

    def dos(self, lineshape="GAU", HWHM=100.0):
        """
        Calculate density of states by computing the overlap of the two FCF spectra
        The two spectra are broadened using the chosen lineshape and Half-Width at Half-Maximum in cm-1
        """
        # Initialize spectrum
        self._spcinit()
        # Get stick spectra
        absstick = self._getstick(self.absrkf)
        emistick = self._getstick(self.emirkf)
        # Convolute with gaussians
        absspc = numpy.copy(self.spc)
        emispc = numpy.copy(self.spc)
        absspc = self.convolute(absspc, absstick, lineshape, HWHM)
        emispc = self.convolute(emispc, emistick, lineshape, HWHM)
        # Integrate
        # Calculate DOS
        self.spc[:, 1] = absspc[:, 1] * emispc[:, 1]
        dos = self.trapezoid(self.spc)
        return dos

    def newspc(self, spcmin, spcmax, spclen=1000):
        # Dimensions of the spectrum
        delta = (spcmax - spcmin) / (spclen - 1)
        spc = numpy.zeros((spclen, 2), dtype=float)
        for ix in range(spclen):
            spc[ix, 0] = spcmin + delta * ix
        return spc

    def convolute(self, spc, stick, lineshape=None, HWHM=None):
        """
        Convolute stick spectrum with the chosen width and lineshape
        lineshape : Can be Gaussian or Lorentzian
        HWHM      : expressed in cm-1
        """
        if HWHM is None:
            raise ValueError("HWHM not defined")
        if lineshape is None:
            raise ValueError("Lineshape not defined")
        # Data for the convolution
        delta = spc[1, 0] - spc[0, 0]
        if lineshape[0:3].upper() == "GAU":
            # Gaussian lineshape
            idline = 1
            # This includes the Gaussian prefactor and the factor to account for the reduced lineshape width
            # factor = 1. / scipy.special.erf(2*math.sqrt(math.log(2)))
            factor = 1.0188815852036244
            factA = math.sqrt(numpy.log(2.0) / math.pi) * factor / HWHM
            factB = -numpy.log(2.0) / HWHM**2
            # We only convolute between -2HWHM and +2HWHM which accounts for 98.1% of the area
            ishft = math.floor(2 * HWHM / delta)
        elif lineshape[0:3].upper() == "LOR":
            # Lorentzian lineshape
            idline = 2
            # This includes the Lorentzian prefactor and the factor to account for the reduced lineshape width
            factA = (math.pi / math.atan(12) / 2) * HWHM / math.pi
            factB = 0.0
            # We only convolute between -12HWHM and +12HWHM which accounts for 94.7% of the area
            ishft = math.floor(12 * HWHM / delta)
        else:
            raise ValueError("invalid lineshape")
        # Loop over peaks in the stick spectrum
        spclen = numpy.size(spc, 0)
        for ix in range(numpy.size(stick[:, 0])):
            # Find peak position in the convoluted spectrum
            peakpos = 1 + math.floor((stick[ix, 0] - spc[0, 0]) / delta)
            # Convolution interval, limited to save computational time
            i1 = max(peakpos - ishft, 1)
            i2 = min(peakpos + ishft, spclen)
            factor = factA * stick[ix, 1]
            if idline == 1:  # Gaussian
                for i in range(i1, i2):
                    spc[i, 1] = spc[i, 1] + factor * math.exp(factB * (spc[i, 0] - stick[ix, 0]) ** 2)
            elif idline == 2:  # Lorentzian
                for i in range(i1, i2):
                    spc[i, 1] = spc[i, 1] + factor / (HWHM + (spc[i, 0] - stick[ix, 0]) ** 2)
        return spc

    def trapezoid(self, spc):
        """
        Integrate spectrum using the trapezoid rule
        """
        value = (spc[0, 1] + spc[-1, 1]) / 2
        value = value + sum(spc[1:-1, 1])
        # The abscissas must be equally spaced for this to work
        value = value * (spc[1, 0] - spc[0, 0])
        return value

    def __str__(self):
        string = f"Absorption from {self.absrkf.path}\nEmission from {self.emirkf.path}\nAcceptor absorption electronic energy = {self.absele} cm-1\nDonor emission electronic energy = {self.emiele} cm-1"
        return string

"""
Module providing interfaces to atmospheric forward model generators.
"""
from . import _chimera

from astropy.modeling import Fittable1DModel, Parameter
import numpy as np
from scipy.optimize import least_squares


def _instrument_non_uniform_tophat(wlgrid, wno, Fp):
    szmod = wlgrid.shape[0]

    delta = np.zeros(szmod)
    Fint = np.zeros(szmod)
    delta[0:-1] = wlgrid[1:] - wlgrid[:-1]
    delta[szmod - 1] = delta[szmod - 2]
    # pdb.set_trace()
    for i in range(szmod - 1):
        i = i + 1
        loc = np.where((1E4 / wno >= wlgrid[i] - 0.5 * delta[i - 1]) & (
        1E4 / wno < wlgrid[i] + 0.5 * delta[i]))
        Fint[i] = np.mean(Fp[loc])

    loc = np.where((1E4 / wno > wlgrid[0] - 0.5 * delta[0]) & (
    1E4 / wno < wlgrid[0] + 0.5 * delta[0]))
    Fint[0] = np.mean(Fp[loc])

    return Fint

class LineForwardModel(Fittable1DModel):
    inputs = ('wavelength',)
    outputs = ('transmission',)

    Rp = Parameter(default=1.359)
    Rstar = Parameter(default=1.155)
    Mp = Parameter(default=0.69)
    Tirr = Parameter(default=1200.)
    logKir = Parameter(default=-1.5)
    logg1 = Parameter(default=-1)
    logMet = Parameter(default=0)
    logCtoO = Parameter(default=-0.26)
    logPQC = Parameter(default=-5)
    logPQN = Parameter(default=-5)
    logRayAmp = Parameter(default=0)
    RaySlp = Parameter(default=4)
    logPc = Parameter(default=1.5)

    def __init__(self, abscoeff_dir, cea_path,
                 Rp=Rp.default, Rstar=Rstar.default, Mp=Mp.default,
                 Tirr=Tirr.default, logKir=logKir.default,
                 logg1=logg1.default, logMet=logMet.default,
                 logCtoO=logCtoO.default, logPQC=logPQC.default, logPQN=logPQN.default,
                 logRayAmp=logRayAmp.default, RaySlp=RaySlp.default, logPc=logPc.default, abund_path=None):

        self.abscoeff_dir = abscoeff_dir
        self.cea_path = cea_path
        self.abund_path = abund_path
        self.Pgrid, self.Tgrid, self.wno, self.gord, self.wts, self.xsecarr = _chimera.fm.xsects(self.abscoeff_dir)

        super(AstropyForwardModel, self).__init__(Rp, Rstar, Mp, Tirr, logKir,
                                                  logg1, logMet, logCtoO, logPQC,
                                                  logPQN, logRayAmp, RaySlp, logPc)

    def evaluate(self, wavelength, Rp, Rstar, Mp, Tirr, logKir, logg1, logMet, logCtoO,
                 logPQC, logPQN, logRayAmp, RaySlp, logPc):
        # seting up input state vector. Must be in this order as indicies are hard wired in fx inside fm
        #  0    1        2       3     4      5              6           7     8    9     10     11     12

        # for some reason astropy models pass parameters as 1D arrays, which breaks CEA code
        x = np.array([Tirr[0], logKir[0], logg1[0], logMet[0], logCtoO[0], logPQC[0],
                      logPQN[0], Rp[0], Rstar[0], Mp[0], logRayAmp[0], RaySlp[0], logPc[0]])
        # calling forward model
        # thermochemical gas profile scaling factors
        # 0   1    2    3   4    5    6     7    8    9   10    11   12   13    14   15   16   17   18
        # H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He
        gas_scale = np.array(
        [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
         1.])  # can be made free params if desired (won't affect mmw)

        # Return model spectrum, wavenumber grid, and vertical abundance
        # profiles from chemistry
        F, wnocrop, atm = _chimera.fm.fx(x, gas_scale, self.Pgrid,
                                                    self.Tgrid, self.wno,
                                                    self.gord, self.wts,
                                                    self.xsecarr,
                                                    cea_path=self.cea_path,
                                                    abund_path=self.abund_path)

        transmission = _instrument_non_uniform_tophat(wavelength, wnocrop, F)

        return transmission

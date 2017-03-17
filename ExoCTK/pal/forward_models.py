"""
Module providing interfaces to atmospheric forward model generators.
"""
from . import _chimera

import numpy as np
from scipy.optimize import least_squares

class LineForwardModel(object):
    """
    Forward model as calculated by Mike Line's CHIMERA code.
    """

    def __init__(self, abscoeff_dir, cea_path=None, abund_path=None):
        """
        Model initializer.  User's must supply either cea_path or abund_path for
        calculating atmospheric abundances.

        Parameters
        ----------
        abscoeff_dir: str
            path to cross section tables
        cea_path: str
            path to a Chemical Equilibrium with Applications (CEA) executable
        abund_path: str
            path to an abundance grid for interpolation
        """
        self.abscoeff_dir = abscoeff_dir
        self.cea_path = cea_path
        self.abund_path = abund_path
        self.Pgrid, self.Tgrid, self.wno, self.gord, self.wts, self.xsecarr = _chimera.fm.xsects(self.abscoeff_dir)


    def __call__(self, Rp=1.359, Rstar=1.155, Mp=0.690, Tirr=1200., logKir=-1.5,
                 logg1=-1, logMet=0., logCtoO=-0.26, logPQC=-5, logPQN=-5,
                 logRayAmp=0., RaySlp=4., logPc=1.5):
        """
        Calculate a transmission spectrum using Mike Line's CHIMERA code.

        Parameters
        ----------
        Rp: float
            Radius of the planet in units of Jupiter Radius
        Rstar: float
            The stellar radius in units of Solar Radius
        Mp: float
            The planet mass in units of Jupiter Mass
        Tirr: float
            The isothermal temperature at the terminator in Kelvin. If full
            redistribution, this is equilibrium temperature.
        logKir: float
            The TP profile IR opacity, the "vertical" location of the gradient
        logg1: float
            The single channel Vis/IR opacity. Controls the delta T between
            deep T and TOA T
        logMet: float
            The metallicity relative to solar log--solar is 0, 10x=1,
            0.1x = -1 used -1.01*log10(M)+0.6
        logCtoO: float
            The log(C/O) ratio. Log solar is -0.26
        logPQC: float
            CH4, CO, H2O Qunech pressure--forces CH4, CO, and H2O to
            constant value at quench pressure value
        logPQN: float
            N2, NH3 Quench pressure--forces N2 and NH3 to
            --ad hoc for chemical kinetics-- reasonable assumption
        RayAmp: float
            log Rayleigh Haze Amplitude (relative to H2)
        RaySlp: float
            haze slope--4 is Rayeigh, 0 is "gray" or flat
        logPc: float
            log of the hard gray cloud top pressure

        Returns
        -------
        wavelength: np.ndarray
            Wavelength grid in micron,
        transmission: np.ndarray
            Transmission spectrum (R/R
        """
        # seting up input state vector. Must be in this order as indicies are hard wired in fx inside fm
        #  0    1        2       3     4      5              6           7     8    9     10     11     12

        x = np.array([Tirr, logKir, logg1, logMet, logCtoO, logPQC,
                      logPQN, Rp, Rstar, Mp, logRayAmp, RaySlp, logPc])
        # calling forward model
        # thermochemical gas profile scaling factors
        # 0   1    2    3   4    5    6     7    8    9   10    11   12   13    14   15   16   17   18
        # H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He
        gas_scale = np.array(
        [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
         1.])  # can be made free params if desired (won't affect mmw)

        # Return model spectrum, wavenumber grid, and vertical abundance
        # profiles from chemistry
        transmission, wnocrop, atm = _chimera.fm.fx(x, gas_scale, self.Pgrid,
                                                    self.Tgrid, self.wno,
                                                    self.gord, self.wts,
                                                    self.xsecarr,
                                                    cea_path=self.cea_path,
                                                    abund_path=self.abund_path)
        return 1e4/wnocrop, transmission

    def residuals(self, p, y_meas, err):
        """
        the residuals between a forward model and some observed data.

        Parameters
        ----------
        p: np.ndarray
            the forward model parameters
        y_meas: np.ndarray
            The measured values
        err: np.ndarray
            the error on the measured value

        Returns
        -------
        loglikelihood: float

        """
        _, y_mod = self.__call__(*p)
        res = (y_meas - y_mod) / err
        return res

    def fit(self, y_meas, err, p_init=None):
        """
        Calculate the Maximum Likelihood Estimate model for an observed spectrum
        using the Levenberg–Marquardt algorithm.

        Parameters
        ----------
        y_obs: np.ndarray
            the observed tranmission spectrum
        p_init: np.ndarray
            inial parameters to being the optimization.  If None the default
            model parameters are used.

        Returns
        -------
        y_mle: np.ndarray

        """
        if p_init is None:
            p_init = [1.359, 1.155, 0.690, 1200., -1.5, -1, 0., -0.26, -5,
                  -5, 0., 4., 1.5]

        # Minimizing so need the negative log-likelihood
        result = least_squares(self.residuals, p_init, args=(y_meas, err))
        return result

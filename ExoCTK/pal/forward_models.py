from . import chimera

import os
import logging

from astropy.modeling.models import custom_model
import numpy as np
import numba
from scipy.interpolate import griddata


# logger = logging.getLogger()
# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# logger.handlers[0].setFormatter(formatter)
# logger.setLevel(logging.DEBUG)

k_b = 1.380658E-23
pi = 3.14159265358979
amu = 1.6605402E-27
hplanck = 6.6260755E-34
clight = 2.99792458E+08


def read_kempton_CIA_cross_section(fname):
    """
    Read a Collision Induced Absorption cross section table provided with
    Exo_Transmit.

    Parameters
    ----------
    fname: str
        Filename of table to be read

    Returns
    -------
    T: np.ndarray
        The temperature grid points.
    P: np.ndarray
        The pressure grid points.
    wavelength: np.ndarray
        The wavelength grid points.
    sigma: np.ndarray
        The grid of cross section values.
    """
    npressure = 13
    ntemperature = 30
    nlambda = 4616
    T = np.zeros(ntemperature)
    P = np.zeros(npressure)
    wavelength = np.zeros(nlambda)
    data = np.zeros((nlambda, 15, ntemperature))
    with open(fname, 'r') as f:
        contents = f.read().split()
        for k in range(ntemperature):
            T[k] = float(contents.pop(0))
        for k in range(ntemperature):
            # T[k] = float(contents.pop(0))
            data[:, :, k] = np.array(contents[(4616 * 15 + 1) * k + 1: (4616 * 15 + 1) * (k + 1)], dtype=float).reshape((4616, 15))

    absorption = {
        'H2_H2' : data[:, 1, :],
        'H2_He' : data[:, 2, :],
        'H2_H' : data[:, 3, :],
        'H2_CH4' : data[:, 4, :],
        # 'CH4_Ar' : data[:, 5, :],
        'CH4_CH4' : data[:, 6, :],
        'CO2_CO2' : data[:, 7, :],
        'He_H' : data[:, 8, :],
        'N2_CH4' : data[:, 9, :],
        'N2_H2' : data[:, 10, :],
        'N2_N2' : data[:, 11, :],
        'O2_CO2' : data[:, 12, :],
        'O2_N2' : data[:, 13, :],
        'O2_O2' :data[:, 14, :]
    }

    return absorption


def read_kempton_cross_section(fname):
    """
    Read a cross section table provided with Exo_Transmit.

    Parameters
    ----------
    fname: str
        Filename of table to be read.

    Returns
    -------
    T: np.ndarray
        The temperature grid points.
    P: np.ndarray
        The pressure grid points.
    wavelength: np.ndarray
        The wavelength grid points.
    sigma: np.ndarray
        The grid of cross section values.
    """

    npressure = 13
    ntemperature = 30
    nlambda = 4616
    T = np.zeros(ntemperature)
    P = np.zeros(npressure)
    wavelength = np.zeros(nlambda)
    sigma = np.zeros((nlambda, npressure, ntemperature))
    with open(fname, 'r') as f:
        contents = f.read().split()
        for k in range(ntemperature):
            T[k] = float(contents.pop(0))

        for j in range(npressure):
            P[j] = float(contents.pop(0))

        for i in range(nlambda):
            size = (ntemperature + 1) * npressure + 1
            wavelength[i] = contents[i * size]
            sigma[i] = np.array(contents[(1 + i * size):(i + 1) * size], dtype=np.float64).reshape(13, 31)[:, 1:]

    return T, P, wavelength, sigma


def read_kempton_abundance(fname):
    """
    Read a chemistry table provided with Exo_Transmit

    Parameters
    ----------
    fname: str
        Filename of table to be read.

    Returns
    -------
    abundances: dict
        A dictionary of abundance grids for each species
    """
    NPRESSURE=13
    NTEMP=30
    P = np.zeros(NPRESSURE)
    T = np.zeros(NTEMP)
    total = np.zeros((NPRESSURE, NTEMP))
    abundances = {}
    abundances["C"] = np.zeros((NPRESSURE, NTEMP))
    abundances["CH4"] = np.zeros((NPRESSURE, NTEMP))
    abundances["CO"] = np.zeros((NPRESSURE, NTEMP))
    abundances["CO2"] = np.zeros((NPRESSURE, NTEMP))
    abundances["C2H2"] = np.zeros((NPRESSURE, NTEMP))
    abundances["C2H4"] = np.zeros((NPRESSURE, NTEMP))
    abundances["C2H6"] = np.zeros((NPRESSURE, NTEMP))
    abundances["H"] = np.zeros((NPRESSURE, NTEMP))
    abundances["HCN"] = np.zeros((NPRESSURE, NTEMP))
    abundances["HCl"] = np.zeros((NPRESSURE, NTEMP))
    abundances["HF"] = np.zeros((NPRESSURE, NTEMP))
    abundances["H2"] = np.zeros((NPRESSURE, NTEMP))
    abundances["H2CO"] = np.zeros((NPRESSURE, NTEMP))
    abundances["H2O"] = np.zeros((NPRESSURE, NTEMP))
    abundances["H2S"] = np.zeros((NPRESSURE, NTEMP))
    abundances["He"] = np.zeros((NPRESSURE, NTEMP))
    abundances["K"] = np.zeros((NPRESSURE, NTEMP))
    abundances["MgH"] = np.zeros((NPRESSURE, NTEMP))
    abundances["N"] = np.zeros((NPRESSURE, NTEMP))
    abundances["N2"] = np.zeros((NPRESSURE, NTEMP))
    abundances["NO2"] = np.zeros((NPRESSURE, NTEMP))
    abundances["NH3"] = np.zeros((NPRESSURE, NTEMP))
    abundances["NO"] = np.zeros((NPRESSURE, NTEMP))
    abundances["Na"] = np.zeros((NPRESSURE, NTEMP))
    abundances["O"] = np.zeros((NPRESSURE, NTEMP))
    abundances["O2"] = np.zeros((NPRESSURE, NTEMP))
    abundances["O3"] = np.zeros((NPRESSURE, NTEMP))
    abundances["OCS"] = np.zeros((NPRESSURE, NTEMP))
    abundances["OH"] = np.zeros((NPRESSURE, NTEMP))
    abundances["PH3"] = np.zeros((NPRESSURE, NTEMP))
    abundances["SH"] = np.zeros((NPRESSURE, NTEMP))
    abundances["SO2"] = np.zeros((NPRESSURE, NTEMP))
    abundances["SiH"] = np.zeros((NPRESSURE, NTEMP))
    abundances["SiO"] = np.zeros((NPRESSURE, NTEMP))
    abundances["TiO"] = np.zeros((NPRESSURE, NTEMP))
    abundances["VO"] = np.zeros((NPRESSURE, NTEMP))
    with open(fname, 'r') as f:
        contents = f.read().split()
        for i in range(38):
            species = contents.pop(0)

        for i in range(NPRESSURE - 1, -1, -1):
            P[i] = float(contents.pop(0))
            for j in range(NTEMP - 1, -1, -1):
                T[j] = float(contents.pop(0))
                total[i][j] = float(contents.pop(0))
                abundances["C"][i][j] = float(contents.pop(0))
                abundances["CH4"][i][j] = float(contents.pop(0))
                abundances["CO"][i][j] = float(contents.pop(0))
                abundances["OCS"][i][j] = float(contents.pop(0))
                abundances["CO2"][i][j] = float(contents.pop(0))
                abundances["C2H2"][i][j] = float(contents.pop(0))
                abundances["C2H4"][i][j] = float(contents.pop(0))
                abundances["C2H6"][i][j] = float(contents.pop(0))
                abundances["H"][i][j] = float(contents.pop(0))
                abundances["HCN"][i][j] = float(contents.pop(0))
                abundances["HCl"][i][j] = float(contents.pop(0))
                abundances["HF"][i][j] = float(contents.pop(0))
                abundances["H2"][i][j] = float(contents.pop(0))
                abundances["H2CO"][i][j] = float(contents.pop(0))
                abundances["H2O"][i][j] = float(contents.pop(0))
                abundances["H2S"][i][j] = float(contents.pop(0))
                abundances["He"][i][j] = float(contents.pop(0))
                abundances["K"][i][j] = float(contents.pop(0))
                abundances["MgH"][i][j] = float(contents.pop(0))
                abundances["N"][i][j] = float(contents.pop(0))
                abundances["N2"][i][j] = float(contents.pop(0))
                abundances["NO2"][i][j] = float(contents.pop(0))
                abundances["NH3"][i][j] = float(contents.pop(0))
                abundances["NO"][i][j] = float(contents.pop(0))
                abundances["Na"][i][j] = float(contents.pop(0))
                abundances["O"][i][j] = float(contents.pop(0))
                abundances["O2"][i][j] = float(contents.pop(0))
                abundances["O3"][i][j] = float(contents.pop(0))
                abundances["OH"][i][j] = float(contents.pop(0))
                abundances["PH3"][i][j] = float(contents.pop(0))
                abundances["SH"][i][j] = float(contents.pop(0))
                abundances["SO2"][i][j] = float(contents.pop(0))
                abundances["SiH"][i][j] = float(contents.pop(0))
                abundances["SiO"][i][j] = float(contents.pop(0))
                abundances["TiO"][i][j] = float(contents.pop(0))
                abundances["VO"][i][j] = float(contents.pop(0))

    return abundances


molar_masses = {
    'H2': 2.0158,
    'H' : 1.0079,
    'He' : 4.002602,
    'H2O' : 18.0152,
    'CH4' : 16.0423,
    'CO' : 28.010,
    'CO2' : 44.010,
    'O' : 15.9994,
    'C' : 12.0107,
    'N' : 14.0067,
    'NH3' : 17.031,
    'N2' : 28.0134,
    'O2' : 31.9988,
    'O3' : 47.9982,
    'C2H2' : 26.0373,
    'C2H4' : 28.0532,
    'C2H6' : 30.0690,
    'H2CO' : 30.0260,
    'H2S' : 34.0809,
    'HCl' : 36.4609,
    'HCN' : 27.0253,
    'HF' : 20.0063,
    'MgH' : 25.3129,
    'NO' : 30.0061,
    'NO2' : 46.0055,
    'OCS' : 60.0751,
    'OH' : 17.0073,
    'PH3' : 33.9976,
    'SH' : 33.0729,
    'SiH' : 29.0934,
    'SiO' : 44.0849,
    'SO2' : 64.0638,
    'TiO' : 63.8664,
    'VO' : 66.9409,
    'Na' : 22.988977,
    'K' : 39.0983
}


class LineForwardModel(object):
    """
    Forward model as calculated by Mike Line's CHIMERA code.
    """

    def __init__(self, abscoeff_dir):
        self.abscoeff_dir = abscoeff_dir


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

        # Return model spectrum, wavenumber grid, and vertical abundance profiles from chemistry
        transmission, wnocrop, atm = chimera.fm.fx(x, gas_scale, self.abscoeff_dir)
        return 1e4/wnocrop, transmission


class KemptonForwardModel(object):
    """
    Forward model as calculated by Eliza Kempton's Exo_Transmit code.
    """
    def __init__(self, opac_dir):
        """
        Intialize the model by loading all the cross sections.

        Parameters
        ----------
        opac_dir: str
            Directory containing Freedman opacity tables as provided in Exo_Transmit

        """
        self.cross_sections = {}
        for f in os.listdir(opac_dir):
            path = os.path.join(opac_dir, f)
            species = f[4:-4]
            if species == 'CIA':
                logger.info('Reading {} cross section table from {}'.format(species, path))
                self.CIA_absorption = read_kempton_CIA_cross_section(path)
            else:
                logger.info('Reading {} cross section table from {}'.format(species, path))
                self.T, self.P, self.wavelength, sigma = read_kempton_cross_section(path)
                self.cross_sections[species] = sigma

    def __call__(self, eos_file, t_p_file, chem_selection, g=9.8,
                 R_planet=6.40e+06, R_star=7.00e+08, threshold=None,
                 rayleigh=1.0, CIA=True, scattering=True):
        """

        Parameters
        ----------
        eos_file: str
            Exo_Transmit Equation of State file.
        t_p_file: str
            File containing the T-P profile
        chem_selection: list
            List of species to include in calculation.
        g: float
            The surface gravity of the planet in m / s^2
        R_planet: float
            Radius of the Planet in meters
        R_star: float
            Radius of the Star in meters
        threshold: float
            Maximum pressure of the atmosphere.  If None the whole
            T-P profile will be used
        rayleigh: float
            Multiplicative factor to augment Rayleigh scattering.
        CIA: boolean
            If True include collision-induced absorption in opacity.
        scattering: boolean
            If True include scattering in opacity

        Returns
        -------
        wavelength: np.ndarray
            The wavelength values for the spectrum
        transmission: np.ndarray
            The transmission spectrum
        """

        total_opacity = self._total_opacity(eos_file, chem_selection, rayleigh,
                                            CIA, scattering)
        ntau, P_atmos, T_atmos, mu_atmos = self._make_atmosphere(t_p_file,
                                                                 threshold)
        wavelength, transmission = self._transmission(ntau, total_opacity,
                                                      mu_atmos, T_atmos,
                                                      P_atmos,
                                                      g, R_planet, R_star)

        return wavelength, transmission

    def _transmission(self, ntau, total_opacity, mu_atmos, T_atmos, P_atmos, g,
                      R_planet, R_star):
        """
        Calculate the transmission through the atmosphere.
        Parameters
        ----------
        ntau: int
            Number of atmosphere points
        total_opacity: np.ndarray
            Opacity grid
        mu_atmos: np.ndarray
            Molar mass at each point in the atmosphere
        T_atmos: np.ndarray
            Temperature at each point in the atmosphere
        P_atmos: np.ndarray
            Pressure at each point in the atmosphere
        g: float
            The surface gravity of the planet in m / s^2
        R_planet: float
            Radius of the Planet in meters
        R_star: float
            Radius of the Star in meters

        Returns
        -------
        wavelength: np.ndarray
            The wavelength values for the spectrum
        transmission: np.ndarray
            The transmission spectrum
        """
        R = 0
        nlambda = total_opacity.shape[0]
        kappa_nu = np.zeros((nlambda, ntau))
        ds = np.zeros(ntau)
        ds[0] = P_atmos[0] * ((k_b * T_atmos[0])
                              / (mu_atmos[0] * amu * P_atmos[0] * g))
        ds[1:] = np.diff(P_atmos) * ((k_b * T_atmos[1:])
                                     / (mu_atmos[1:] * amu * P_atmos[1:] * g))

        # interpolate total opacity onto atmosphere T-P points
        T_grid, P_grid = np.meshgrid(self.T, self.P)
        for i in range(nlambda):
            kappa_nu[i] = griddata((T_grid.flatten(), P_grid.flatten()),
                                   total_opacity[i].flatten(),
                                   (T_atmos, P_atmos))

        tau_tr = _tau_los(kappa_nu, ds, R_planet, nlambda)
        theta, dtheta = _angles(ds, R_planet, R_star)

        R -= R_planet + np.sum(ds)
        logger.info("R {:f}\n".format(1.0 - R ** 2 / R_star ** 2))

        t_star = 6000.0

        flux_pl = np.zeros_like(self.wavelength)
        flux_st = np.zeros_like(self.wavelength)
        flux_tr = np.zeros_like(self.wavelength)
        intensity = np.zeros_like(self.wavelength)
        for i in range(nlambda):
            flux_pl[i] = 0.0
            intensity[i] = Planck(t_star, self.wavelength[i])

            for j in range(ntau):
                flux_pl[i] += (intensity[i] * np.exp(-tau_tr[i][j]) *
                               np.cos(theta[j]) * np.sin(theta[j]) * dtheta[j])

            flux_pl[i] *= 2 * pi

            flux_st[i] = pi * intensity[i]
            flux_tr[i] = ((1.0 - (R_planet + np.sum(ds)) ** 2 / R_star ** 2) *
                          flux_st[i] + flux_pl[i])

        return self.wavelength, 100.0 * (1.0 - flux_tr / flux_st)

    def _make_atmosphere(self, t_p_file, threshold):
        """
        Get the pressure, temperature and molar mass of atmosphere
        layers.

        Parameters
        ----------
        t_p_file: str
            File containing the T-P profile
        threshold: float
            Maximum pressure of the atmosphere.  If None the whole T-P
            profile will be used

        Returns
        -------
        ntau: int
            Number of layers in the atmosphere
        P_atmos: np.ndarray
            Pressure at each atmosphere layer
        T_atmos: np.ndarray
            Temperature at each atmosphere layer
        mu_atmos: np.ndarray
            Molar mass at each atmosphere layer
        """

        _, P_atmos, T_atmos = np.loadtxt(t_p_file, unpack=True, skiprows=1)

        if threshold:
            ntau = (P_atmos < threshold).sum() + 1
            T_atmos = T_atmos[:ntau]
            P_atmos = P_atmos[:ntau]
        else:
            ntau = P_atmos.shape[0]

        T_grid, P_grid = np.meshgrid(self.T, self.P)
        mu_atmos = griddata((T_grid.flatten(), P_grid.flatten()),
                            self.mu.flatten(), (T_atmos, P_atmos))


        if threshold:
            proportion = ((np.log10(threshold) - np.log10(P_atmos[-2])) /
                          (np.log10(P_atmos[-1]) - np.log10(P_atmos[-2])))
            P_atmos[-1] = threshold
            T_atmos[-1] = proportion*(T_atmos[-1] - T_atmos[-2]) + T_atmos[-2]
            mu_atmos[-1] = proportion*(mu_atmos[-1] -
                                       mu_atmos[-2]) + mu_atmos[-2]

        return ntau, P_atmos, T_atmos, mu_atmos

    def _total_opacity(self, eos_file, chem_selection, rayleigh, CIA,
                       scattering):
        """
        Calculate the total opacity grid.

        Parameters
        ----------
        eos_file: str
            Exo_Transmit Equation of State file.
        chem_selection: list
            List of species to include in calculation.
        rayleigh: float
            Multiplicative factor to augment Rayleigh scattering.
        CIA: boolean
            If True include collision-induced absorption in opacity.
        scattering: boolean
            If True include scattering in opacity

        Returns
        -------
        total_opacity: np.ndarray
            The total opacity grid
        """
        abundances = read_kempton_abundance(eos_file)
        total_opacity = np.zeros_like(self.cross_sections['H2O'])
        self.mu = np.zeros_like(self.cross_sections['H2O'][0])
        for species, abundance in abundances.items():
            self.mu += molar_masses[species] * abundance

        for species in chem_selection:
            opac = self.cross_sections[species] \
                   * (abundances[species] *
                      np.divide.outer(self.P, k_b * self.T))[np.newaxis, :, :]
            total_opacity += opac

        if CIA:
            for pair, absorption in self.CIA_absorption.items():
                s1, s2 = pair.split('_')
                logger.info('Including collision-induced absorption between {} and {}'.format(s1, s2))
                total_opacity += absorption[:, np.newaxis, :] * (abundances[s1] * abundances[s2] * (np.divide.outer(P, (k_b * T)))**2)[np.newaxis, :, :]

        if scattering:
            logger.info('Including scattering')
            scat_opacity = np.zeros_like(total_opacity)
            scat_opacity += ((8.0*pi/3.0) * 0.80e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['H2'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 0.21e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['He'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 1.74e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['N2'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 1.45e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['H2O'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 1.95e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['CO'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 2.91e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['CO2'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 2.26e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['NH3'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 2.59e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['CH4'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 1.58e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['O2'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 3.21e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['O3'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 3.33e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['C2H2'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 4.25e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['C2H4'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 4.47e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['C2H6'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 2.59e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['HCN'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 2.63e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['HCl'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 0.80e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['HF'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 1.70e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['NO'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 3.02e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['NO2'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity += ((8.0*pi/3.0) * 4.84e-30**2 * (2.0*pi/ self.wavelength)**4)[:, np.newaxis, np.newaxis] * \
            (abundances['PH3'] * np.divide.outer(self.P, (k_b * self.T)))[np.newaxis, :, :]

            scat_opacity *= rayleigh
            total_opacity += scat_opacity

        return total_opacity


@numba.jit
def _tau_los(kappa_nu, ds, Rpl, NLam):
    R = Rpl + np.sum(ds)

    tau_tr = np.zeros((NLam, len(ds)))
    dl = np.zeros_like(ds)
    a = R
    for j in range(len(ds)):
        a -= ds[j]
        b = R

        for k in range(j + 1):
            dl[k] = 2.0 * np.sqrt(b ** 2 - a ** 2)
            b -= ds[k]

        for k in range(j + 1):
            if (k != j):
                dl[k] -= dl[k + 1]

        for k in range(j + 1):
            for i in range(NLam):
                tau_tr[i][j] += kappa_nu[i][k] * dl[k]

    return tau_tr

def _angles(ds, Rpl, Rst):

    R = Rpl + np.sum(ds)
    h = R
    theta = np.zeros_like(ds)
    dtheta = np.zeros_like(ds)
    for j in range(len(ds)):
        h -= ds[j];
        theta[j] = np.arcsin(h/Rst);

        if j == 0:
            dtheta[j] = np.arcsin(R/Rst) - theta[j]
        else:
            dtheta[j] = theta[j-1] - theta[j]

    return theta, dtheta

def Planck(T, lam):
    MAX_EXPONENT = 400.0

    hc_Tkla = (hplanck * clight) / (T * k_b * lam)
    twohnu3_c2 = (2.0 * hplanck * clight) / lam**3

    if (hc_Tkla <= MAX_EXPONENT):
        Bnu = twohnu3_c2 / (np.exp(hc_Tkla) - 1.0)
    else:
        Bnu = 0.0

    return Bnu
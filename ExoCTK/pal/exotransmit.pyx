"""
This module provides a wrapper for the Exo_Transmit C code which 
generates transmission spectra to study exoplanet atmospheres.
"""

import os

import numpy as np
cimport numpy as np
from libc.stdio cimport printf
from libc.stdlib cimport malloc, free


DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

cdef extern from "vars.h":
    ctypedef struct vars:
        int NTAU

        int NTEMP
        double TLOW
        double THIGH

        int NPRESSURE
        double PLOW
        double PHIGH
        double THRESHOLD
        double RAYLEIGH

        int NLAMBDA
        double G
        double R_PLANET
        double R_STAR
        double T_STAR

        char* tpfname
        char* eosfname

        char** fileArray

        int chemselection[32]

cdef extern from "getVars.c":
    vars getVars()


cdef extern from "atmos.h":
    pass

cdef extern from "opac.h":
    pass

cdef extern from "nrutil.h":
    pass

cdef extern from "nrutil.c":
    pass

cdef extern from "constant.h":
    pass

cdef extern from "prototypes.h":
    pass

cdef extern from "readchemtable.c":
    pass

cdef extern from "read_t_p.c":
    pass

cdef extern from "rt_transmission.c":
    pass

cdef extern from "geometry.c":
    pass

cdef extern from "utils.c":
    pass

cdef extern from "planck.c":
    pass

cdef extern from "totalopac.c":
    pass 

cdef extern from "getChemSelection.c":
    pass

cdef extern from "getFileArray.c":
    char** getFileArray()

cdef extern from "getNTau.c":
    pass

cdef extern from "interpol.c":
    pass

cdef extern from "stdio.h":
    pass

cdef extern from "readopactable.c":
    pass

cdef extern from "main_transmission.c":
    void transmission(vars variables, double* wavelength, double* flux)

def _make_chem_selection(chem):
    """
    Create the a boolean chemSelection array in the order Exo_Transmit expects.

    Parameters
    ----------
    chem: list
        A list of strings specifying species/phenomena to be included in opacity calcuation

    Returns
    -------
    array
        int array specifying chemistry selection in the correct order 
    """
    cdef np.ndarray[long, ndim=1] chemselection = np.zeros(32, dtype=long)

    chemselection[0]  = "CH4" in chem  
    chemselection[1]  = "CO2" in chem 
    chemselection[2]  = "CO" in chem 
    chemselection[3]  = "H2O" in chem 
    chemselection[4]  = "NH3" in chem 
    chemselection[5]  = "O2" in chem 
    chemselection[6]  = "O3" in chem 
    chemselection[7]  = "C2H2" in chem 
    chemselection[8]  = "C2H4" in chem 
    chemselection[9]  = "C2H6" in chem 
    chemselection[10] = "H2CO" in chem 
    chemselection[11] = "H2S" in chem 
    chemselection[12] = "HCl" in chem 
    chemselection[13] = "HCN" in chem 
    chemselection[14] = "HF" in chem 
    chemselection[15] = "MgH" in chem 
    chemselection[16] = "N2" in chem 
    chemselection[17] = "NO" in chem 
    chemselection[18] = "NO2" in chem 
    chemselection[19] = "OCS" in chem 
    chemselection[20] = "OH" in chem 
    chemselection[21] = "PH3" in chem 
    chemselection[22] = "SH" in chem 
    chemselection[23] = "SiH" in chem 
    chemselection[24] = "SiO" in chem 
    chemselection[25] = "SO2" in chem 
    chemselection[26] = "TiO" in chem 
    chemselection[27] = "VO" in chem 
    chemselection[28] = "Na" in chem 
    chemselection[29] = "K" in chem
    chemselection[30] = "Scattering" in chem
    chemselection[31] = "Collision Induced Absorption" in chem
    return chemselection 

cdef char ** to_cstring_array(list_str):
    """
    Convert a Python list of strings to a C char** array points.

    Parameters
    ----------
    list_str: list
        list containing strings

    Returns
    -------
    char**
        pointer to a 2d char array
    """
    cdef char **ret = <char **>malloc(len(list_str) * sizeof(char *))
    b_list = [s.encode() for s in list_str]
    for i in range(len(list_str)):
        b_entry = b_list[i]
        ret[i] = b_entry

    return ret

def _get_spectrum(n_tau=334, n_temp=30, T_low=100.0, T_high=3000.0, n_pressure=13, P_low=1.0e-4, 
    P_high=1.0e8, threshold=0.0, rayleigh=1.0, n_lambda=7454, g=9.8, R_planet=6.40e+6, R_star=7.00e+8, 
    tpfname='T_P/t_p_800K.dat', eosfname='EOS/eos_1Xsolar_cond.dat', file_array=None,
    chemistry=["CH4", "CO2", "CO", "H2O", "NH3", "O2", "O3", "C2H2", "C2H4", "C2H6", "H2CO", "H2S", "HCl", 
    "HCN", "HF", "MgH", "N2", "NO", "NO2", "OCS", "OH", "PH3", "SH", "SiH", "SiO", "SO2", "TiO", "VO", 
    "Na", "K", "Scattering", "Collision Induced Absorption"]):
    """
    Run ExoTransmit main to produce a transmission spectrum.
    """
    cdef vars variables
    variables.NTAU = n_tau

    variables.NTEMP = n_temp
    variables.TLOW = T_low
    variables.THIGH = T_high

    variables.NPRESSURE = n_pressure
    variables.PLOW = P_low
    variables.PHIGH = P_high
    variables.THRESHOLD = threshold
    variables.RAYLEIGH = rayleigh

    variables.NLAMBDA = n_lambda
    variables.G = g
    variables.R_PLANET = R_planet
    variables.R_STAR = R_star

    chemselection = _make_chem_selection(chemistry).data
    variables.chemselection = chemselection 

    cdef bytes b_tpfname = tpfname.encode()
    variables.tpfname = b_tpfname

    cdef bytes b_eosfname = eosfname.encode()
    variables.eosfname = b_eosfname

    cdef char** string_buf = to_cstring_array(file_array)

    variables.fileArray = string_buf

    cdef np.ndarray[DTYPE_t, ndim=1] wavelength = np.empty(variables.NLAMBDA, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] flux = np.empty(variables.NLAMBDA, dtype=DTYPE)

    transmission(variables, <double*> wavelength.data, <double*> flux.data)
    free(string_buf)
    return wavelength, flux

class ExoTransmit(object):
    """
    Creates an object for running ExoTransmit.

    Attributes
    ----------

    tpfname: str
        Temperature-pressure data file path
    eosfname: str
        Chemistry (gas abundance) data file path
    g: float
        Planet surface gravity [m / s^2]
    R_planet: float
        Planet radius [m]
    R_star: float
        Steller radius [m]
    chemistry: list
        List of species/phenomenon to take into account when calculating opacity.
    rayleigh: float
        Rayleigh scattering augmentation factor
    threshold: float
        Pressure [Pa] of the top of an optically thick cloud deck.
    """

    def __init__(self):
        """
        Initializes the ExoTransmit object with some default parameters.

        threshold = 0.0

        rayleigh = 1.0
        
        g = 9.8
        
        R_planet = 6.40e6
        
        R_star = 7.00e8
        
        tpfname = ExoCTK.pal.__file__/data/T_P/t_p_800K.dat
        
        tpfname = ExoCTK.pal.__file__/data/EOS/eos_0p1XSolar_cond.dat
        
        chemistry = ["CH4", "CO2", "CO", "H2O", "NH3", "O2", "O3", "C2H2", "C2H4", "C2H6", "H2CO", "H2S", "HCl", 
        "HCN", "HF", "MgH", "N2", "NO", "NO2", "OCS", "OH", "PH3", "SH", "SiH", "SiO", "SO2", "TiO", "VO", "Na", "K", 
        "Scattering", "Collision Induced Absorption"]
        """
        self.threshold = 0.0
        self.rayleigh = 1.0
        self.g = 9.8
        self.R_planet = 6.40e6
        self.R_star = 7.00e8
        self.tpfname = os.path.join(os.path.dirname(__file__), "data/T_P/t_p_800K.dat")
        self.eosfname = os.path.join(os.path.dirname(__file__), "data/EOS/eos_0p1Xsolar_cond.dat")
        self.chemistry = ["CH4", "CO2", "CO", "H2O", "NH3", "O2", "O3", "C2H2", "C2H4", "C2H6", "H2CO", "H2S", "HCl", 
        "HCN", "HF", "MgH", "N2", "NO", "NO2", "OCS", "OH", "PH3", "SH", "SiH", "SiO", "SO2", "TiO", "VO", "Na", "K", 
        "Scattering", "Collision Induced Absorption"]

    def __str__(self):
        parts = []
        parts.append("tpfname = {}".format(self.tpfname))
        parts.append("eosfname = {}".format(self.eosfname))
        parts.append("g = {}".format(self.g))
        parts.append("R_planet = {}".format(self.R_planet))
        parts.append("R_star = {}".format(self.R_star))
        parts.append("chemistry = {}".format(self.chemistry))
        parts.append("threshold = {}".format(self.threshold))
        parts.append("rayleigh = {}".format(self.rayleigh))

        return '\n'.join(parts)

    def __call__(self):
        """
        Generate a transmission spectrum.

        Parameters should be set by modifying the ExoTransmit object attributes.

        Returns
        -------
        wavelength: array
            The spectrum wavelengths [m]
        transmission: array
            The transmission [percent]
        """

        # fileArray expected by Exo_Transmit including opacity tables. ORDER MATTERS!
        fileArray = [os.path.join(os.path.dirname(__file__), "data/T_P/t_p_800K.dat"), 
            os.path.join(os.path.dirname(__file__), "data/EOS/eos_0p1Xsolar_cond.dat"),
            os.path.join(os.path.dirname(__file__), "data/Spectra/test3.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacCH4.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacC2H2.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacC2H4.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacC2H6.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacCO.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacCO2.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacH2CO.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacH2O.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacH2S.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacHCN.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacHCl.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacHF.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacMgH.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacN2.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacNH3.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacNO.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacNO2.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacO2.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacO3.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacOCS.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacOH.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacPH3.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacSH.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacSO2.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacSiH.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacSiO.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacTiO.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacVO.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacNa.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacK.dat"),
            os.path.join(os.path.dirname(__file__), "data/Opac/opacCIA.dat"),
            os.path.join(os.path.dirname(__file__), "data")]
        
        # Would be specified in OtherInput.in not meant to be specified by user
        n_tau = 334
        n_temp = 30
        T_low = 100.0
        T_high = 3000.0
        n_pressure = 13
        P_low = 1.0e-4
        P_high = 1.0e8
        n_lambda = 7454
        
        # Call the C code
        spec = _get_spectrum(n_tau=n_tau, n_temp=n_temp, T_low=T_low, T_high=T_high, 
            n_pressure=n_pressure, P_low=P_low, P_high=P_high, threshold=self.threshold, 
            rayleigh=self.rayleigh, n_lambda=n_lambda, g=self.g, R_planet=self.R_planet, R_star=self.R_star, 
            tpfname=self.tpfname, eosfname=self.eosfname, file_array=fileArray, chemistry=self.chemistry)

        return spec



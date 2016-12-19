import numpy as np
cimport numpy as np
from libc.stdio cimport printf

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
    pass

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

chem = ["CO2", "CO", "H2O", "NH3", "O2", "O3", "C2H2", "C2H4", "C2H6", "H2CO", 
"H2S", "HCl", "HCN", "HF", "MgH", "N2", "NO", "NO2", "OCS", "OH", "PH3", "SH", "SiH", 
"SiO", "SO2", "TiO", "VO", "Na", "K", "Scattering", "Collision Induced Absorption"]

def make_chem_selection(chem):
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

def get_spectrum(n_tau=334, n_temp=30, T_low=100.0, T_high=3000.0, n_pressure=13, P_low=1.0e-4, 
    P_high=1.0e8, threshold=0.0, rayleigh=1.0, n_lambda=7454, g=9.8, R_planet=6.40e+6, R_star=7.00e+8, 
    tpfname='T_P/t_p_800K.dat', eosfname='EOS/eos_1Xsolar_cond.dat', 
    chemistry=["CH4", "CO2", "CO", "H2O", "NH3", "O2", "O3", "C2H2", "C2H4", "C2H6", "H2CO", "H2S", "HCl", 
    "HCN", "HF", "MgH", "N2", "NO", "NO2", "OCS", "OH", "PH3", "SH", "SiH", "SiO", "SO2", "TiO", "VO", 
    "Na", "K", "Scattering", "Collision Induced Absorption"]):
    """
    Run ExoTransmit main to produce a transmission spectrum
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

    chemselection = make_chem_selection(chemistry).data
    variables.chemselection = chemselection 

    cdef bytes b_tpfname = tpfname.encode()
    variables.tpfname = b_tpfname

    cdef bytes b_eosfname = eosfname.encode()
    variables.eosfname = b_eosfname

    cdef np.ndarray[DTYPE_t, ndim=1] wavelength = np.empty(variables.NLAMBDA, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] flux = np.empty(variables.NLAMBDA, dtype=DTYPE)

    transmission(variables, <double*> wavelength.data, <double*> flux.data)
    return wavelength, flux

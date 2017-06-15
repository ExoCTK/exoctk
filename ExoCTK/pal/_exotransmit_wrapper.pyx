"""
This module provides a wrapper for the Exo_Transmit C code which 
generates transmission spectra to study exoplanet atmospheres.
"""

cdef extern from "include/main_transmission.c":
    void main()

def exotransmit():
    main()

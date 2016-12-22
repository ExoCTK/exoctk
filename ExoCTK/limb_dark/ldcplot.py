#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module of plotting tools for the limb darkening subpackage.
"""
import matplotlib
import matplotlib.pyplot as plt
import os
from astropy.io import fits

def ld_v_mu(data):
    """
    Plot the limb darkening versus mu
    
    Parameters
    ----------
    data: array-like, str
        The list of (coefficient, mu) values or the path
        to a FITS file with COEFFS and MU extensions
    
    """
    # Create the plot
    plt.figure()
    
    # Get the data
    if isinstance(data, (list,tuple,np.ndarray)):
        coeffs, mu = data[:2]
    
    else:
        HDU = fits.open()
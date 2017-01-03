#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module of plotting tools for the limb darkening subpackage.
"""
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import os
from astropy.io import fits

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#params = [np.unique(model_grid.data[p]) for p in ['Teff','logg','FeH']]

def ld_v_mu(data, params):
    """
    Plot the limb darkening versus mu
    
    Parameters
    ----------
    data: array-like, str
        The list of (coefficient, mu) values or the path
        to a FITS file with COEFFS and MU extensions
    
    """    
    try:
    
        # Get the data directly...
        if isinstance(data, (list,tuple,np.ndarray)):
            coeffs, mu = data[:2]

        # ...or from a FITS file
        else:
            coeffs = fits.getdata(path, extname='COEFFS')
            mu = fits.getdata(path, extname='MU')
        
        # Create the plot
        plt.figure()
        plt.xlabel(r'$\mu$m')
        plt.ylabel=('Limb Darkening')
        
        # Get the param ranges
        
        for c,m in zip(coeffs, mu):
            plt.scatter(coeffs, mu)
    
    except:
        print('Sorry, that input is not plottable!')
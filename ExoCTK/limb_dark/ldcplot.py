#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module of plotting tools for the limb darkening subpackage.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, cm
import os
import copy
from astropy.io import fits
try:
    from ExoCTK import core
    from ExoCTK.limb_dark import ldcfit
except ImportError:
    from ExoCTK.ExoCTK import core
    from ExoCTK.ExoCTK.limb_dark import ldcfit

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def ld_v_mu(model_grid, compare, profiles=('quadratic','nonlinear'), **kwargs):
    """
    Produce a plot of mu vs. limb darkening values for a range of 
    model grid parameters compared to some interpolated value
    
    Parameters
    ----------
    model_grid: core.ModelGrid
        The model grid to plot. The ModelGrid.customize() method can be run
        before hand or arguments can be passed via **kwargs argument
    compare: list
        The list of off-grid parameters to compare to the
        model grid plots
    profiles: list
        The limb darkening profiles to include in the comparison
        
    Example
    -------
    >>> model_comparison(model_grid, (2550, 5.22, 0), **{'Teff_rng':(2500,2600), 
                         'logg_rng':(5,5.5), 'FeH_rng':(0,0)})
    
    """
    # Make a copy of the ModelGrid so it doesn't change the 
    # input object in the Python session
    grid = copy.copy(model_grid)
    grid.customize(**kwargs)
    
    # Plot the grids
    for p,ls in zip(profiles,['-','--',':','-.']):
        _ = ldcfit.ldc_grid(grid, p, plot=plt.gcf(), **{'ls':ls})
    
    # Plot the interpolated comparison
    for p,ls in zip(profiles,['-','--',':','-.']):
        _ = ldcfit.ldc(*compare, grid, p, plot=plt.gcf(), **{'ls':ls, 'c':'k', 'lw':3})
    
    # Delete the copy
    del grid
    
    # Plot labels
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$I(\mu)/I(\mu =0)$')

def ldc_v_wavelength(model_grid, wave_ranges, profile, **kwargs):
    """
    Plot limb darkening coefficients vs. wavelength as resolution of intensity spectra.
    Overplot mean value for specified filters and/or wavelength intervals
    """
    # Draw plot
    fig = plt.figure()
    plt.xlabel('Wavelength')
    plt.ylabel('Coefficients')
    
    means, wavelength, unc = [], [], []
    
    # Calculate the grid in the given wavelength ranges
    for wr in wave_ranges:
        
        # Get wave center and range
        w = (wr[0]+wr[1])/2.
        w_unc = wr[1]-w
        
        # Make a copy of the ModelGrid so it doesn't change the 
        # input object in the Python session
        grid = copy.copy(model_grid)
        
        # Get wavelength segment
        grid.customize(wave_rng=wr)
        
        # Calculate the coefficients
        c_grid, m_grid, r_grid = ldcfit.ldc_grid(grid, profile)
        
        # Calculate mean and store
        #m = calculate_coeff_means_somehow
        means.append(m)
        wavelength.append(w)
        unc.append(w_unc)
        
        # Plot the values
        plt.plot(wavelength, coeffs)
        
        del grid
        
    # Plot the mean
    plt.errorbar(wavelength, means, xerr=unc, c='k')
    

#def ld_v_mu(data, params):
#    """
#    Plot the limb darkening versus mu
#    
#    Parameters
#    ----------
#    data: array-like, str
#        The list of (coefficient, mu) values or the path
#        to a FITS file with COEFFS and MU extensions
#    
#    """    
#    try:
#    
#        # Get the data directly...
#        if isinstance(data, (list,tuple,np.ndarray)):
#            coeffs, mu = data[:2]
#
#        # ...or from a FITS file
#        else:
#            coeffs = fits.getdata(path, extname='COEFFS')
#            mu = fits.getdata(path, extname='MU')
#        
#        # Create the plot
#        plt.figure()
#        plt.xlabel(r'$\mu$m')
#        plt.ylabel=('Limb Darkening')
#        
#        # Get the param ranges
#        
#        for c,m in zip(coeffs, mu):
#            plt.scatter(coeffs, mu)
#    
#    except:
#        print('Sorry, that input is not plottable!')
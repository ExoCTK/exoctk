#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module of plotting tools for the limb darkening subpackage.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import copy
import inspect
from bokeh.plotting import figure
from matplotlib import rc, cm
from astropy.io import fits
try:
    from .. import core
except:
    from ExoCTK import core
try:
    from . import ldcfit
except:
    from ExoCTK.ldc import ldcfit

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def ld_plot(coeffs, ldfunc, fig=None, params=[], **kwargs):
    """
    Make a LD plot in Bokeh
    
    Parameters
    ----------
    
    """
    # Evaluate the limb darkening profile
    mu_vals = np.linspace(0, 1, 1000)
    ld_vals = ldfunc(mu_vals, *coeffs)
    
    # Plot stuff
    label = ', '.join(map(str,params))
    
    if not fig or fig==True:
        
        fig = plt.figure(figsize=(10,4))
        
    # Is it a matplotlib plot?
    if isinstance(fig, matplotlib.figure.Figure):
        
        # Make axes
        ax = fig.add_subplot(111)
        
        # Draw the curve
        p = ax.plot(mu_vals, ld_vals, label=label, **kwargs)
        color = p[0].get_color()
        
        # Plot the mu_values
        # try:
        #     ld_raw = np.mean(flux, axis=1)/np.mean(flux, axis=1)[np.where(mu_raw==1)]
        #     ld_err = np.std(flux, axis=1)/np.std(flux, axis=1)[np.where(mu_raw==1)]
        #     ax.errorbar(mu_raw, ld_raw, yerr=ld_err, c=color, ls='None', marker='o',
        #              markeredgecolor=color, markerfacecolor='None')
        #
        # except:
        #     pass
        
        plt.xlim(0,1)
        plt.ylim(0,1)
    
    # Otherwise it mush be bokeh!
    else:
        
        fig.line(mu_vals, ld_vals, legend=label)
        

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
        t, g, m = compare
        _ = ldcfit.ldc(t, g, m, grid, p, plot=plt.gcf(), **{'ls':ls, 'c':'k', 'lw':3})
    
    # Delete the copy
    del grid
    
    # Plot labels
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$I(\mu)/I(\mu =0)$')

def ldc_v_wavelength(model_grid, wave_ranges, profile, **kwargs):
    """
    Plot limb darkening coefficients vs. wavelength as resolution of intensity spectra.
    Overplot mean value for specified filters and/or wavelength intervals
    
    Parameters
    ----------
    model_grid: core.ModelGrid
        The model grid to use for the calculation
    wave_ranges: list
        The list of wavelength range tuples in [um],
        e.g. [(0.6,0.9), (0.9,1.2), (1.2,1.5)]
    profile: str
        The limb darkening profile to use
    
    """
    # Get the number of coefficients for the limb darkening profile
    rows = len(inspect.getargspec(ldcfit.ld_profile(profile)).args)-1
    cols = len(wave_ranges)
    
    # Initialize limb darkening coefficient, mu, and effecive radius grids
    T = model_grid.Teff_vals
    G = model_grid.logg_vals
    M = model_grid.FeH_vals
    coeff_grid = np.zeros((cols,rows,len(T),len(G),len(M)))
    mu_grid = np.zeros((cols,len(T),len(G),len(M)))
    r_grid = np.zeros((cols,len(T),len(G),len(M)))
    w_mean, w_unc = [], []
    
    # Calculate the grid in the given wavelength ranges
    for n,wr in enumerate(wave_ranges):
        
        # Get wave center and range
        w = (wr[0]+wr[1])/2.
        w_mean.append(w)
        w_unc.append(wr[1]-w)
        
        # Make a copy of the ModelGrid so it doesn't change the 
        # input object in the Python session
        grid = copy.copy(model_grid)
        
        # Apply wavelength segment to model grid
        grid.customize(wave_rng=wr)
        
        # Calculate the coefficient, mu, and radius grids
        cg, mg, rg = ldcfit.ldc_grid(grid, profile)
        
        # Add them to the arrays
        coeff_grid[n] = cg
        mu_grid[n] = mg
        r_grid[n] = rg
        
        del grid
    
    # Draw plot
    C = ['c{}'.format(i+1) for i in range(rows)]
    fig = core.multiplot(rows, 1, xlabel='Wavelength', sharex=True, sharey=False, ylabel=C, 
                         title=profile.title())
    
    # Reshape the coeff grid so that the coefficients are the first dimension
    coeff_grid = coeff_grid.reshape(rows,len(T),len(G),len(M),cols)
    
    # For each coefficient, make a plot
    for n,coeffs in enumerate(coeff_grid):
        for nt,t in enumerate(T):
            for ng,g in enumerate(G):
                for nm,m in enumerate(M):
                    fig[n+1].plot(w_mean, coeff_grid[n,nt,ng,nm], label=[t,g,m], marker='o')
    
    # Plot a legend
    fig[-1].legend(loc=0, frameon=False, fontsize=15)
    diff = np.diff(w_mean)[0]/4.
    plt.xlim(min(w_mean)-diff,max(w_mean)+diff)

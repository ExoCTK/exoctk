#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module to calculate limb darkening coefficients from a grid of model spectra
"""
import numpy as np
        
def calculate_ldc(model_grid, orders, mu_min=0.02):
    """
    Calculates the limb darkening coefficients for a given 
    grid of synthetic spectra
    
    Reference for limb-darkening laws:
        http://www.astro.ex.ac.uk/people/sing/David_Sing/Limb_Darkening.html
    
    Parameters
    ----------
    model_grid: core.ModelGrid object
        The grid of synthetic spectra from which the coefficients will
        be calculated 
    orders: int
        The polynomial order, i.e. the number of coefficients
    mu_min: float
        The minimum mu value to consider
    
    Returns
    -------
    list
        The list of limb darkening coefficients, mu values, and effective 
        radii calculated from the input core.ModelGrid
    
    """
    # Initialize limb darkening coefficient, mu, and effecive radius grids
    T, G, M = [np.unique(model_grid.data[p]) for p in ['Teff','logg','FeH']]
    ldc = np.zeros((orders,len(T),len(G),len(M)))
    mu0 = np.zeros((len(T),len(G),len(M)))
    r_eff = np.zeros((len(T),len(G),len(M)))
    
    # Iterate through spectra files and populate grids
    for f in model_grid.data:
        
        # Get the physical parameters for this model
        t, g, m = [f[p] for p in ['Teff','logg','FeH']]
        
        # Locate the grid position for this model
        t_idx, g_idx, m_idx = [np.where(A==a)[0][0] for A,a in 
                               zip([T,G,M],[t,g,m])]
        
        # Retrieve the wavelength, flux, mu, and effective radius
        wave, flux, mu, radius = model_grid.get(t, g, m)
        
        # Put the effective radius value in the radius grid
        r_eff[t_idx,g_idx,m_idx] = radius
        
        # Calculate mean intensity vs. mu
        mean_i = np.mean(flux, axis=1)
        
        # Calculate limb darkening, I[mu]/I[1] vs. mu
        ld = mean_i/mean_i[np.where(mu==1)]
        
        # Rescale mu values. Spherical Phoenix models extend beyond limb
        muz = np.interp(0.01, ld, mu)
        mu = (mu-muz)/(1-muz)
        mu0[t_idx,g_idx,m_idx] = muz
        
        # Trim to useful mu range
        imu = np.where(mu>mu_min)
        mu, ld = mu[imu], ld[imu]
        
        # Fit limb darkening to get limb darkening coefficients (LDCs)
        err = 1.
        ldc0 = np.ones(orders)
        if orders==2:
            ldc[:,t_idx,g_idx,m_idx] = mpfitfun('ld2func',mu,ld,err,ldc0)
        elif orders==4:
            ldc[:,t_idx,g_idx,m_idx] = mpfitfun('ld4func',mu,ld,err,ldc0)
        else:
           print('Order number must be 2 or 4.')
           return
        
    return ldc, mu0, r_eff
    
def ldfunc(mu, coeffs, order=2):
    """
    Define the 2nd or 4th order function to fit to the limb darkening profile
    
    Parameters
    ----------
    mu: float
        The mu value for the limb darkening function
    coeffs: array-like
        The polynomial coefficients for the function
    order: int
        The order of the polynomial, must be 2 or 4
    
    Returns
    -------
    float
        The value of the evaluated function with the given mu value and
        coefficients
    """
    if order==2:    
        return 1. - coeffs[0]*(1.-mu)\
                  - coeffs[1]*(1.-mu)^2
    
    elif order==4:
        return 1. - coeffs[0]*(1.-mu^0.5)\
                  - coeffs[1]*(1.-mu)    \
                  - coeffs[2]*(1.-mu^1.5)\
                  - coeffs[3]*(1.-mu^2)
                  
    else:
        return
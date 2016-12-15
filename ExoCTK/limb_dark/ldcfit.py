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
    #Initialize coefficient grid
    T, G, M = [np.unique(model_grid.data[p]) for p in ['Teff','logg','FeH']]
    ldc = np.zeros((orders,len(T),len(G),len(M)))
    mu0 = np.zeros((len(T),len(G),len(M)))
    r_eff = np.zeros((len(T),len(G),len(M)))
    
    # Iterate through spectra files and populate coefficient grid
    for f in model_grid.data:
        
        # Get the physical parameters for this model
        t, g, m = [f[p] for p in ['Teff','logg','FeH']]
        
        # Locate the grid position for this model
        t_idx, g_idx, m_idx = [np.where(A==a)[0][0] for A,a in zip([T,G,M],[t,g,m])]
        
        # Retrieve the wavelength, flux, mu, and effective radius
        wave, flux, mu, radius = model_grid.get(t, g, m)
        
        # Put the effective radius value in the radius grid
        r_eff[t_idx,g_idx,m_idx] = radius
        
        # Calculate mean intensity vs. mu
        mean_i = np.mean(flux, axis=0)
        
        # Calculate limb darkening, I[mu]/I[1] vs. mu
        center_i, = np.where(mu0==1)
        ld = mean_i/mean_i[center_i]
        
        # Rescale mu values. Spherical Phoenix models extend beyond limb
        muz = np.interp(mu, ld, 0.01)
        mu = (mu-muz)/(1-muz)
        mu0[t_idx,g_idx,m_idx] = muz
        
        # Trim to useful mu range
        imu = np.where(mu>mu_min)
        mu, ld = mu[imu], ld[imu]
        
        # Fit limb darkening to get limb darkening coefficients (LDCs)
        # ===================================================================
        # What the hell is ldc0?
        # 
        # Need to write mpfitfun in Python!
        # ===================================================================
        err = 1.
        #ldc0 = replicate(1d0,npar) # ?????????????????????
        if orders==2:
            ldc[:,t_idx,g_idx,m_idx] = mpfitfun('ld2func',mu,ld,err,ldc0)
        elif orders==4:
            ldc[:,t_idx,g_idx,m_idx] = mpfitfun('ld4func',mu,ld,err,ldc0)
        else:
            print('Order number must be 2 or 4.')
            return
        
    return ldc, mu0, r_eff
    
    
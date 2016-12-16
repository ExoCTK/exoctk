#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module to calculate limb darkening coefficients from a grid of model spectra
"""
import numpy as np
from scipy.optimize import curve_fit
        
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
        
        # Define the fitting function given the number of orders
        if orders==2:
            def ldfunc(m, c1, c2):
                return 1. - c1*(1.-m) - c2*(1.-m)**2
        elif orders==4:
            def ldfunc(m, c1, c2, c3, c4):
                return 1. - c1*(1.-m**0.5) - c2*(1.-m) \
                          - c3*(1.-m**1.5) - c4*(1.-m**2)
        else:
           print('Order number must be 2 or 4.')
           return
        
        # Fit limb darkening to get limb darkening coefficients (LDCs)
        ldc[:,t_idx,g_idx,m_idx] = curve_fit(ldfunc, mu, ld, method='lm')[0]
        
    return ldc, mu0, r_eff
    
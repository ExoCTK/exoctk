#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module to calculate limb darkening coefficients from a grid of model spectra
"""
import numpy as np
import inspect
from ExoCTK.ExoCTK import core
from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator

def ld_profile(name='quadratic'):
    """
    Define the function to fit the limb darkening profile
    
    Reference:
        https://www.cfa.harvard.edu/~lkreidberg/batman/tutorial.html#limb-darkening-options
        
    Parameters
    ----------
    name: str
        The name of the limb darkening profile function to use, 
        including 'uniform', 'linear', 'quadratic', 'square-root', 
        'logarithmic', 'exponential', and 'nonlinear'
        
    Returns
    -------
    function
        The corresponding function for the given profile
        
    """
    # Supported profiles a la BATMAN
    names = ['uniform','linear','quadratic','square-root',
             'logarithmic','exponential','nonlinear']
    
    # Check that the profile is supported
    if name in names:

        # Uniform
        if name=='uniform':
            def profile(m):
                return 1.
            
        # Linear
        if name=='linear':
            def profile(m, c1):
                return 1. - c1*(1.-m)
        
        # Quadratic
        if name=='quadratic':
            def profile(m, c1, c2):
                return 1. - c1*(1.-m) - c2*(1.-m)**2
            
        # Square-root
        if name=='square-root':
            def profile(m, c1, c2):
                return 1. - c1*(1.-m) - c2*m*log(m)
            
        # Exponential
        if name=='exponential':
            def profile(m, c1, c2):
                return 1. - c1*(1.-m) - c2/(1.-e**m)
        
        # Nonlinear
        if name=='nonlinear':
            def profile(m, c1, c2, c3, c4):
                return 1. - c1*(1.-m**0.5) - c2*(1.-m) \
                          - c3*(1.-m**1.5) - c4*(1.-m**2)
        
        return profile
        
    else:
        print(name,'is not a supported profile. Try',names)
        return
        

def ldc(teff, logg, FeH, model_grid, profile, mu_min=0.02):
    """
    Calculates the limb darkening coefficients for a given synthetic spectrum.
    If the model grid does not contain a spectrum of the given parameters, the
    grid is interpolated to those parameters.
    
    Reference for limb-darkening laws:
        http://www.astro.ex.ac.uk/people/sing/David_Sing/Limb_Darkening.html
    
    Parameters
    ----------
    teff: int
        The effective temperature of the model
    logg: float
        The logarithm of the surface gravity
    FeH: float
        The logarithm of the metallicity
    model_grid: core.ModelGrid object
        The grid of synthetic spectra from which the coefficients will
        be calculated 
    profile: str
        The name of the limb darkening profile function to use, 
        including 'uniform', 'linear', 'quadratic', 'square-root', 
        'logarithmic', 'exponential', and 'nonlinear'
    mu_min: float
        The minimum mu value to consider
    
    Returns
    -------
    np.ndarray
        The list of limb darkening coefficients, mu values, and effective 
        radius calculated from the model of the given parameters from the
        input core.ModelGrid 
    
    """
    # Define the limb darkening profile function
    ldfunc = ld_profile(profile)
    
    if not ldfunc:
        return
        
    else:
    
        # See if the model with the desired parameters is on the grid
        in_grid = model_grid.data[[(model_grid.data['Teff']==teff)&
                                   (model_grid.data['logg']==logg)&
                                   (model_grid.data['FeH']==FeH)]]\
                                   in model_grid.data

        # If a model with the given parameters exists, calculate it
        if in_grid:

            # Retrieve the wavelength, flux, mu, and effective radius
            spec_dict = model_grid.get(teff, logg, FeH)
            wave = spec_dict.get('wave')
            flux = spec_dict.get('flux')
            mu = spec_dict.get('mu')
            radius = spec_dict.get('r_eff')

            # Calculate mean intensity vs. mu
            mean_i = np.mean(flux, axis=1)

            # Calculate limb darkening, I[mu]/I[1] vs. mu
            ld = mean_i/mean_i[np.where(mu==1)]

            # Rescale mu values. Spherical Phoenix models extend beyond limb
            muz = np.interp(0.01, ld, mu)
            mu = (mu-muz)/(1-muz)

            # Trim to useful mu range
            imu = np.where(mu>mu_min)
            mu, ld = mu[imu], ld[imu]

            # Fit limb darkening to get limb darkening coefficients (LDCs)
            coeffs = curve_fit(ldfunc, mu, ld, method='lm')[0]

        # If a model with the given parameters doesn't exist, 
        # calculate ALL grid values and interpolate
        else:

            # Print that it has to calculate
            print('Teff:', teff, ' logg:', logg, ' FeH:', FeH, 
                  ' model not in grid. Calculating...')

            # Get values for the entire model grid
            coeff_grid, mu_grid, r_grid = ldc_grid(model_grid, orders, 
                                                   mu_min=mu_min)

            # Cretae a grid of the parameter values to interpolate over
            params = [np.array(np.unique(model_grid.data[p])) 
                      for p in ['Teff','logg','FeH']]

            # Interpolate mu value
            interp_muz = RegularGridInterpolator(params, mu_grid)
            muz, = interp_muz(np.array([teff,logg,FeH]))

            # Interpolate effective radius value
            interp_r = RegularGridInterpolator(params, r_grid)
            radius, = interp_r(np.array([teff,logg,FeH]))

            # Interpolate coefficients
            coeffs = []
            for c_grid in coeff_grid:
                interp_coeff = RegularGridInterpolator(params, c_grid)
                coeffs.append(interp_coeff(np.array([teff,logg,FeH])))
            coeffs = np.array(coeffs).flatten()

        return coeffs, muz, radius
    
        
def ldc_grid(model_grid, profile, write_to='', mu_min=0.02):
    """
    Calculates the limb darkening coefficients for a given 
    grid of synthetic spectra
    
    Parameters
    ----------
    model_grid: core.ModelGrid object
        The grid of synthetic spectra from which the coefficients will
        be calculated 
    profile: str
        The name of the limb darkening profile function to use, 
        including 'uniform', 'linear', 'quadratic', 'square-root', 
        'logarithmic', 'exponential', and 'nonlinear'
    write_to: str
        The path and filename to write the results to
    mu_min: float
        The minimum mu value to consider
    
    Returns
    -------
    list
        The list of limb darkening coefficients, mu values, and effective 
        radii calculated from the input core.ModelGrid
    
    """
    # Get the arguments for the limb darkening profile
    C = inspect.getargspec(ld_profile(profile)).args
    
    # Initialize limb darkening coefficient, mu, and effecive radius grids
    T, G, M = [np.unique(model_grid.data[p]) for p in ['Teff','logg','FeH']]
    coeff_grid = np.zeros((len(C)-1,len(T),len(G),len(M)))
    mu_grid = np.zeros((len(T),len(G),len(M)))
    r_grid = np.zeros((len(T),len(G),len(M)))
    
    # Iterate through spectra files and populate grids
    for f in model_grid.data:
        
        # Get the physical parameters for this model
        t, g, m = [f[p] for p in ['Teff','logg','FeH']]
        
        # Locate the grid position for this model
        t_idx, g_idx, m_idx = [np.where(A==a)[0][0] for A,a in 
                               zip([T,G,M],[t,g,m])]
                               
        # Fit limb darkening to get limb darkening coefficients (LDCs)
        coeffs, muz, radius = ldc(t, g, m, model_grid, profile, mu_min)
        
        # Add the coefficients, mu values and effective radius to grids
        coeff_grid[:,t_idx,g_idx,m_idx] = coeffs
        mu_grid[t_idx,g_idx,m_idx] = muz
        r_grid[t_idx,g_idx,m_idx] = radius
    
    # Write the results to file
    if write_to:
        
        # FITS file format
        if write_to.endswith('.fits'):
            extensions = {k:v for k,v in zip(['COEFFS','MU','RADII'],
                                             [coeff_grid, mu_grid,r_grid])}
            core.writeFITS(write_to, extensions, headers=())
    
        # ASCII? Numpy? JSON?
        else:
            pass
    
    # Or return them
    else:
        return coeff_grid, mu_grid, r_grid
    
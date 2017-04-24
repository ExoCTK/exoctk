#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module to calculate limb darkening coefficients from a grid of model spectra
"""
import numpy as np
import inspect
import datetime
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.optimize import curve_fit
from . import ldcplot as lp
from .. import core

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
   
def ld_profile(name='quadratic', latex=False):
    """
    Define the function to fit the limb darkening profile
    
    Reference:
        https://www.cfa.harvard.edu/~lkreidberg/batman/tutorial.html#limb-darkening-options
        
    Parameters
    ----------
    name: str
        The name of the limb darkening profile function to use, 
        including 'uniform', 'linear', 'quadratic', 'square-root', 
        'logarithmic', 'exponential', '3-parameter', and '4-parameter'
    latex: bool
        Return the function as a LaTeX formatted string
        
    Returns
    -------
    function, str
        The corresponding function for the given profile
        
    """
    # Supported profiles a la BATMAN
    names = ['uniform','linear','quadratic','square-root',
             'logarithmic','exponential','3-parameter','4-parameter']
    
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
                return 1. - c1*(1.-m) - c2*(1.-np.sqrt(m))
        
        # Logarithmic
        if name=='logarithmic':
            def profile(m, c1, c2):
                return 1. - c1*(1.-m) - c2*m*(1.-np.log(m))
            
        # Exponential
        if name=='exponential':
            def profile(m, c1, c2):
                return 1. - c1*(1.-m) - c2/(1.-np.e**m)
        
        # 3-parameter
        if name=='3-parameter':
            def profile(m, c1, c2, c3):
                return 1. -  c1*(1.-m) - c2*(1.-m**1.5) - c3*(1.-m**2)
        
        # 4-parameter
        if name=='4-parameter':
            def profile(m, c1, c2, c3, c4):
                return 1. - c1*(1.-m**0.5) - c2*(1.-m) \
                          - c3*(1.-m**1.5) - c4*(1.-m**2)
        
        if latex:
            profile = inspect.getsource(profile).replace('\n','')
            profile = profile.replace('\\','').split('return ')[1]
            
            for i,j in [('**','^'),('m','\mu'),(' ',''),('np.','\\'),
                        ('0.5','{0.5}'),('1.5','{1.5}')]:
                profile = profile.replace(i,j)
        
        return profile
        
    else:
        print("'{}' is not a supported profile. Try".format(name),names)
        return
        

def ldc(Teff, logg, FeH, model_grid, profiles, mu_min=0.05, ld_min=0.001, 
        bandpass='', plot=False, **kwargs):
    """
    Calculates the limb darkening coefficients for a given synthetic spectrum.
    If the model grid does not contain a spectrum of the given parameters, the
    grid is interpolated to those parameters.
    
    Reference for limb-darkening laws:
        http://www.astro.ex.ac.uk/people/sing/David_Sing/Limb_Darkening.html
    
    Parameters
    ----------
    Teff: int
        The effective temperature of the model
    logg: float
        The logarithm of the surface gravity
    FeH: float
        The logarithm of the metallicity
    model_grid: core.ModelGrid object
        The grid of synthetic spectra from which the coefficients will
        be calculated 
    profiles: str, list
        The name(s) of the limb darkening profile function to use, 
        including 'uniform', 'linear', 'quadratic', 'square-root', 
        'logarithmic', 'exponential', and '4-parameter'
    mu_min: float
        The minimum mu value to consider
    ld_min: float
        The minimum limb darkening value to consider
    bandpass: core.Filter() (optional)
        The photometric filter through which the limb darkening
        is to be calculated
    plot: bool, matplotlib.figure.Figure, bokeh.plotting.figure.Figure
        Plot mu vs. limb darkening for this model in an existing
        figure or in a new figure
    
    Returns
    -------
    np.ndarray
        The list of limb darkening coefficients, mu values, and effective 
        radius calculated from the model of the given parameters from the
        input core.ModelGrid 
    
    """              
    # Get the model, interpolating if necessary
    grid_point = model_grid.get(Teff, logg, FeH)
    
    # If the model exists, continue
    if grid_point:
            
        # Retrieve the wavelength, flux, mu, and effective radius
        wave = grid_point.get('wave')
        flux = grid_point.get('flux')
        mu = grid_point.get('mu').squeeze()
        radius = grid_point.get('r_eff')
        
        # Apply the filter if any
        if isinstance(bandpass, core.Filter):
            flux = bandpass.apply([wave,flux])
            wave = bandpass.rsr[0]
            
            if bandpass.WavelengthMin/10000<model_grid.wavelength[0]\
            or bandpass.WavelengthMax/10000>model_grid.wavelength[-1]:
                print('Bandpass {} not covered by'.format(bandpass.filterID))
                print('model grid of wavelength range',model_grid.wave_rng)
                
                return
            
        # Calculate mean intensity vs. mu
        mean_i = np.mean(flux, axis=1)
        
        # Calculate limb darkening, I[mu]/I[1] vs. mu
        ld = mean_i/mean_i[np.where(mu==1)]
        
        # Rescale mu values to make f(mu=0)=ld_min
        # for the case where spherical models extend beyond limb
        muz = np.interp(ld_min, ld, mu) if any(ld<ld_min) else 0
        mu = (mu-muz)/(1-muz)
        
        # Trim to useful mu range
        mu_raw = mu.copy()
        imu = np.where(mu>mu_min)
        mu, ld = mu[imu], ld[imu]
        
        # Iterate through the requested profiles
        if isinstance(profiles, str):
            profiles = [profiles]
        for profile in profiles:
                        
            # Define the limb darkening profile function
            ldfunc = ld_profile(profile)
            
            if not ldfunc:
                return
                
            else:
                
                # Make dict for profile
                grid_point[profile] = {}
                
                # Fit limb darkening to get limb darkening coefficients
                coeffs, cov = curve_fit(ldfunc, mu, ld, method='lm')
                
                # Add coeffs and errors to the dictionary
                grid_point[profile]['coeffs'] = coeffs
                grid_point[profile]['err'] = np.sqrt(np.diag(cov))
                
        # Add raw data and inputs
        grid_point['flux'] = flux
        grid_point['wave'] = wave
        grid_point['mu_min'] = mu_min
        grid_point['r_eff'] = radius
        grid_point['profiles'] = profiles
        grid_point['bandpass'] = bandpass
        
        if plot:
            
            # Make a list of LD functions and plot them
            ldfuncs = [ld_profile(profile) for profile in profiles]
            lp.ld_plot(ldfuncs, grid_point, fig=plot, **kwargs)
        
        return grid_point
            
    # If the desired params are not within the grid bounds, return
    else:
        # Print that it cannot calculate
        print('Teff:', Teff, ' logg:', logg, ' FeH:', FeH, 
              ' model not within grid bounds', model_grid.Teff_rng,
              model_grid.logg_rng, model_grid.FeH_rng)
              
        return

def ldc_grid(model_grid, profile, write_to='', mu_min=0.05, plot=False, **kwargs):
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
        'logarithmic', 'exponential', and '4-parameter'
    write_to: str
        The path and filename to write the results to
    mu_min: float
        The minimum mu value to consider
    plot: bool, matplotlib.figure.Figure
        Plot mu vs. limb darkening for this model in an existing
        figure or in a new figure
        
    Returns
    -------
    list
        The list of limb darkening coefficients, mu values, and effective 
        radii calculated from the input core.ModelGrid
    
    """
    # Get the arguments for the limb darkening profile
    C = inspect.getargspec(ld_profile(profile)).args
    
    # Initialize limb darkening coefficient, mu, and effecive radius grids
    T = model_grid.Teff_vals
    G = model_grid.logg_vals
    M = model_grid.FeH_vals
    coeff_grid = np.zeros((len(C)-1,len(T),len(G),len(M)))
    mu_grid = np.zeros((len(T),len(G),len(M)))
    r_grid = np.zeros((len(T),len(G),len(M)))
    
    if plot:
        
        # If a figure is not passed, make one
        if not isinstance(plot, plt.Figure):
            fig = plt.figure()
            plt.xlabel(r'$\mu$')
            plt.ylabel(r'$I(\mu)/I(\mu =0)$')
        
        # If a figure is passed, proceed
        else:
            fig = plot
    
    else:
        
        # No figures for me, thank you!
        fig = None
    
    # Iterate through spectra files and populate grids
    for f in model_grid.data:
        
        # Get the physical parameters for this model
        t, g, m = [f[p] for p in ['Teff','logg','FeH']]
        
        # Locate the grid position for this model
        t_idx, g_idx, m_idx = [np.where(A==a)[0][0] for A,a in 
                               zip([T,G,M],[t,g,m])]
                               
        # Fit limb darkening to get limb darkening coefficients (LDCs)
        coeffs, muz, radius = ldc(t, g, m, model_grid, profile, 
                                  mu_min, plot=fig, **kwargs)
        
        # Add the coefficients, mu values and effective radius to grids
        coeff_grid[:,t_idx,g_idx,m_idx] = coeffs
        mu_grid[t_idx,g_idx,m_idx] = muz
        r_grid[t_idx,g_idx,m_idx] = radius 
        
    # Write legend
    if plot and not isinstance(plot, plt.Figure):
        plt.legend(loc=0, frameon=False)
    
    # Write the results to file
    if write_to:
        
        # Collect keys for the header
        hdr = []
        date = str(datetime.datetime.now())

        # From this calculation
        hdr.append(('PROFILE', profile, 'The limb darkening profile used'))
        hdr.append(('DATE', date, 'The data the file was generated'))

        # ...and from the ModelGrid() object
        for k,v in model_grid.__dict__.items():
            if isinstance(v, (list,str,int,float,tuple)):
                if isinstance(v, (list,tuple)):
                    v = repr(v)
                hdr.append((k.upper()[:8], v, 'core.ModelGrid() attribute'))
        
        # FITS file format
        if write_to.endswith('.fits'):
            
            # Create the extensions
            extensions = {k:v for k,v in zip(['COEFFS','MU','RADII'],
                                             [coeff_grid, mu_grid,r_grid])}
            
            # Write the FITS file
            core.writeFITS(write_to, extensions, headers=hdr)
            
        # ASCII? Numpy? JSON?
        else:
            pass
            
    # Or return them
    else:
        return coeff_grid, mu_grid, r_grid
    

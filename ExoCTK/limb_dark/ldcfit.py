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
from scipy.interpolate import RegularGridInterpolator
try:
    from ExoCTK import core
except ImportError:
    from ExoCTK.ExoCTK import core

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
    
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
        

def ldc(teff, logg, FeH, model_grid, profile, mu_min=0.05, bandpass='', 
        plot=False, **kwargs):
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
    bandpass: core.Filter() (optional)
        The photometric filter through which the limb darkening
        is to be calculated
    plot: bool, matplotlib.figure.Figure
        Plot mu vs. limb darkening for this model in an existing
        figure or in a new figure
    
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
    
        # See if the model with the desired parameters is witin the grid
        in_grid = all([(teff>=model_grid.Teff_rng[0])&
                       (teff<=model_grid.Teff_rng[1])&
                       (logg>=model_grid.logg_rng[0])&
                       (logg<=model_grid.logg_rng[1])&
                       (FeH>=model_grid.FeH_rng[0])&
                       (FeH<=model_grid.FeH_rng[1])])
                       
        # Caluclate if the parameters are within the grid
        if in_grid:
            
            # See if the model with the desired parameters is a true grid point
            on_grid = model_grid.data[[(model_grid.data['Teff']==teff)&
                                       (model_grid.data['logg']==logg)&
                                       (model_grid.data['FeH']==FeH)]]\
                                       in model_grid.data
                                       
            # If a model is a true grid point, just calculate it
            if on_grid:
                
                # Retrieve the wavelength, flux, mu, and effective radius
                spec_dict = model_grid.get(teff, logg, FeH)
                wave = spec_dict.get('wave')
                flux = spec_dict.get('flux')
                mu = spec_dict.get('mu')
                radius = spec_dict.get('r_eff')
                
                # Apply the filter if any
                if isinstance(bandpass, core.Filter):
                    flux = bandpass.convolve([wave,flux])
                    wave = bandpass.rsr[0]

                # Calculate mean intensity vs. mu
                mean_i = np.mean(flux, axis=1)
                
                # Calculate limb darkening, I[mu]/I[1] vs. mu
                ld = mean_i/mean_i[np.where(mu==1)]
                
                # Rescale mu values. Spherical Phoenix models extend beyond limb
                muz = np.interp(0.01, ld, mu)
                mu = (mu-muz)/(1-muz)
                mu_raw = mu.copy()
                
                # Trim to useful mu range
                imu = np.where(mu>mu_min)
                mu, ld = mu[imu], ld[imu]
                
                # Fit limb darkening to get limb darkening coefficients (LDCs)
                coeffs = curve_fit(ldfunc, mu, ld, method='lm')[0]
                
            # If a model with the given parameters is not a true grid point 
            # but is within the grid range, calculate ALL grid values 
            # and interpolate
            else:
                
                # Print that it has to calculate
                print('Teff:', teff, ' logg:', logg, ' FeH:', FeH, 
                      ' model not in grid. Calculating...')
                
                # Get values for the entire model grid
                coeff_grid, mu_grid, r_grid = ldc_grid(model_grid, profile, mu_min=mu_min)
                                                       
                # Create a grid of the parameter values to interpolate over,
                # eliminating parameters that can't be interpolated
                params, values = [], []
                for p,v in zip([model_grid.Teff_vals, 
                                model_grid.logg_vals,
                                model_grid.FeH_vals],
                               [teff, logg, FeH]):
                    if len(p)>1:
                        params.append(p)
                        values.append(v)
                          
                # Interpolate mu value
                interp_muz = RegularGridInterpolator(params, mu_grid)
                muz, = interp_muz(np.array(values))
                
                # Interpolate effective radius value
                interp_r = RegularGridInterpolator(params, r_grid)
                radius, = interp_r(np.array(values))
                
                # Interpolate coefficients
                coeffs = []
                for c_grid in coeff_grid:
                    interp_coeff = RegularGridInterpolator(params, c_grid)
                    coeffs.append(interp_coeff(np.array(values)))
                coeffs = np.array(coeffs).flatten()
                
            if plot:
                
                # Make a new figure if necessary
                if not isinstance(plot, plt.Figure):
                    fig = plt.figure(figsize=(10,4))
                    ax = plt.subplot()
                    plt.xlabel(r'$\mu$')
                    plt.ylabel(r'$I(\mu)/I(\mu =0)$')
                    
                    # Print the calculation parameters
                    info = r'${}-{}$ $\mu m$'.format(*model_grid.wave_rng)
                    if isinstance(bandpass, core.Filter):
                        info = bandpass.params['filterID']
                    info += '\n'
                    info += ', '.join([r'$c_{}={:.2f}$'.format(n+1,c) 
                                       for n,c in enumerate(coeffs)])
                    ax.text(0.95, 0.1, info, ha='right', transform=ax.transAxes,
                            linespacing=2)
                    
                    
                # Evaluate the limb darkening profile
                mu_vals = np.linspace(0, 1, 1000)
                ld_vals = ldfunc(mu_vals, *coeffs)
                
                # Draw the curve
                label = '{}, {}, {}'.format(teff, logg, FeH)
                p = plt.plot(mu_vals, ld_vals, label=label, **kwargs)
                color = p[0].get_color()
                
                # Plot the mu_values
                try:
                    ld_raw = np.mean(flux, axis=1)/np.mean(flux, axis=1)[np.where(mu_raw==1)]
                    plt.plot(mu_raw, ld_raw, c=color, ls='None', marker='o',
                             markeredgecolor=color, markerfacecolor='None')
                except:
                    pass
                
                # New figure stuff
                if not isinstance(plot, plt.Figure):
                    
                    # Plot the minimum mu cutoff and legend
                    plt.axvline(x=mu_min, c='0.5', ls=':', label=r'$\mu$ cutoff')
                    plt.legend(loc=2, frameon=False)
                
                plt.xlim(0,1)
                plt.ylim(0,1)
                    
            return coeffs, muz, radius
            
        # If the desired params are not within the grid bounds, return
        else:
            # Print that it cannot calculate
            print('Teff:', teff, ' logg:', logg, ' FeH:', FeH, 
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
        'logarithmic', 'exponential', and 'nonlinear'
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
    
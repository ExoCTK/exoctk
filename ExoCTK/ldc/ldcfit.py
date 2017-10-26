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
import astropy.table as at
import astropy.units as q
from matplotlib import rc
from scipy.optimize import curve_fit
from . import ldcplot as lp
from .. import core
from .. import svo

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
            def profile(m, c1):
                return c1
            
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
        

def ldc(Teff, logg, FeH, model_grid, profiles, mu_min=0.05, ld_min=1E-6, 
        bandpass='', grid_point='', plot=False, save=False, **kwargs):
    """
    Calculates the limb darkening coefficients for a given synthetic spectrum.
    If the model grid does not contain a spectrum of the given parameters, the
    grid is interpolated to those parameters.
    
    Reference for limb-darkening laws:
        http://www.astro.ex.ac.uk/people/sing/David_Sing/Limb_Darkening.html
    
    Parameters
    ----------
    Teff: int, sequence
        The effective temperature of the model
    logg: float, sequence
        The logarithm of the surface gravity
    FeH: float, sequence
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
    bandpass: svo.Filter() (optional)
        The photometric filter through which the limb darkening
        is to be calculated
    grid_point: dict (optional)
        A previously computed model grid point, rather
        than providing Teff, logg, and FeH
    plot: bool, matplotlib.figure.Figure, bokeh.plotting.figure.Figure
        Plot mu vs. limb darkening for this model in an existing
        figure or in a new figure
    save: str
        Save the plot and the table of coefficients to file
    
    Returns
    -------
    np.ndarray
        The list of limb darkening coefficients, mu values, and effective 
        radius calculated from the model of the given parameters from the
        input core.ModelGrid 
    
    """              
    # Get the model, interpolating if necessary
    if not grid_point:
        grid_point = model_grid.get(Teff, logg, FeH)
    
    # If the model exists, continue
    if grid_point:
            
        # Retrieve the wavelength, flux, mu, and effective radius
        wave = grid_point.get('wave')
        flux = grid_point.get('flux')
        mu = grid_point.get('mu').squeeze()
        radius = grid_point.get('r_eff')
        
        # Check if a bandpass is provided
        if isinstance(bandpass, svo.Filter):
            
            # Make sure the bandpass has coverage
            if bandpass.WavelengthMin*q.Unit(bandpass.WavelengthUnit)\
                <model_grid.wave_rng[0]*model_grid.wl_units\
            or bandpass.WavelengthMax*q.Unit(bandpass.WavelengthUnit)\
                >model_grid.wave_rng[-1]*model_grid.wl_units:
                print('Bandpass {} not covered by'.format(bandpass.filterID))
                print('model grid of wavelength range',model_grid.wave_rng)
                
                return
            
            else:
                # Apply the filter
                flux = bandpass.apply([wave,flux])
                
                # Make rsr curve 3 dimensions if there is only one
                # wavelength bin, then get wavelength only
                bp = bandpass.rsr
                if len(bp.shape)==2:
                    bp = bp[None,:]
                wave = bp[:,0,:]
            
        # Calculate mean intensity vs. mu
        wave = wave[None,:] if len(wave.shape)==1 else wave
        flux = flux[None,:] if len(flux.shape)==2 else flux
        mean_i = np.nanmean(flux, axis=-1)
        mean_i[mean_i==0] = np.nan
        
        # Calculate limb darkening, I[mu]/I[1] vs. mu
        ld = mean_i/mean_i[:,np.where(mu==max(mu))].squeeze(axis=-1)
        
        # Rescale mu values to make f(mu=0)=ld_min
        # for the case where spherical models extend beyond limb
        ld_avg = np.nanmean(ld, axis=0)
        muz = np.interp(ld_min, ld_avg, mu) if any(ld_avg<ld_min) else 0
        mu = (mu-muz)/(1-muz)
        grid_point['scaled_mu'] = mu
        grid_point['ld_raw'] = ld
        
        # Trim to useful mu range
        mu_raw = mu.copy()
        imu, = np.where(mu>mu_min)
        mu, ld = mu[imu], ld[:,imu]
        
        # Add raw data and inputs
        grid_point['flux'] = flux
        grid_point['wave'] = wave
        grid_point['mu_min'] = mu_min
        grid_point['r_eff'] = radius
        grid_point['bandpass'] = bandpass
        
        if isinstance(profiles, str):
            profiles = [profiles]
        grid_point['profiles'] = profiles
        
        if isinstance(bandpass, svo.Filter):
            grid_point['n_bins'] = bandpass.n_bins
            grid_point['pixels_per_bin'] = bandpass.pixels_per_bin
            grid_point['centers'] = bandpass.centers.round(5)
        else:
            grid_point['n_bins'] = 1
            grid_point['pixels_per_bin'] = wave.shape[-1]
            grid_point['centers'] = np.array([(wave[-1]+wave[0])/2.]).round(5)
        
        # Iterate through the requested profiles
        for profile in profiles:
                        
            # Define the limb darkening profile function
            ldfunc = ld_profile(profile)
            
            if not ldfunc:
                return
                
            else:
                
                # Make dict for profile
                grid_point[profile] = {}
                
                # Fit limb darkening to get limb darkening
                # coefficients for each wavelength bin
                all_coeffs, all_errs = [], []
                cen = grid_point['centers'][0]
                c = range(len(inspect.signature(ldfunc).parameters)-1)
                for w,l in zip(cen,ld):
                    coeffs, cov = curve_fit(ldfunc, mu, l, method='lm')
                    err = np.sqrt(np.diag(cov))
                    all_coeffs.append([w]+list(coeffs))
                    all_errs.append(list(err))
                    
                # Make a table of coefficients
                all_coeffs = list(zip(*all_coeffs))
                c_cols = ['wavelength']+['c{}'.format(n+1) for n in c]
                c_table = at.Table(all_coeffs, names=c_cols)
                
                # Make a table of errors
                all_errs = list(zip(*all_errs))
                e_cols = ['e{}'.format(n+1) for n in c]
                e_table = at.Table(all_errs, names=e_cols)
                
                # Combine, format, and store tables
                cols = ['wavelength']+','.join(['c{0},e{0}'.format(n+1) for n in c]).split(',')
                grid_point[profile]['coeffs'] = at.hstack([c_table,e_table])[cols]
                for k in c_cols[1:]+e_cols:
                    grid_point[profile]['coeffs'][k].format = '%.3f'
                
        # Make a table for each profile then stack them so that
        # the columns are ['Profile','c0','e0',...,'cn','en']
        for p in grid_point['profiles']:
            print(p,':')
            grid_point[p]['coeffs'].pprint(max_width=-1)
            print('\r')
            
            # Write the table to file
            if save:
                with open(save, 'a') as f:
                    f.write('Profile: '+p+'\n')
                    grid_point[p]['coeffs'].write(f, format='ascii.ipac')
                    f.write('\r')

        if plot:

            # Make a list of LD functions and plot them
            ldfuncs = [ld_profile(profile) for profile in profiles]
            lp.ld_plot(ldfuncs, grid_point, fig=plot, **kwargs)
        
        return grid_point
            
    # If the desired params are not within the grid bounds, return
    else:
        # Print that it cannot calculate
        print('Teff:', Teff, ' logg:', logg, ' FeH:', FeH, 
              ' model not within grid bounds', (min(model_grid.Teff_vals),
              max(model_grid.Teff_vals)), (min(model_grid.logg_vals),
              max(model_grid.logg_vals)), (min(model_grid.FeH_vals),
              max(model_grid.FeH_vals)))
              
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
    

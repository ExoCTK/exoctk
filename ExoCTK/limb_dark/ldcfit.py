# Reference for model spectra: Husser, T.-O. et al. (2013)
#   http://adsabs.harvard.edu/abs/2013A%26A...553A...6H
#   http://phoenix.astro.physik.uni-goettingen.de/?page_id=73
# 


from glob import glob
from astropy.io import fits
import astropy.table as at
import numpy as np
import os

class ModelGrid(object):
    """
    Creates a ModelGrid object which contains a multi-parameter
    grid of model spectra and its references
    """
    
    def __init__(self, spec_files, bibcode='2013A&A...553A...6H'):
        """
        Initializes the model grid by creating a table with a column
        for each parameter and ingests the spectra
    
        Parameters
        ----------
        spec_files: str
            The path plus wildcard filename for the FITS files of spectra
        bibcode: str, array-like
            The bibcode or list of bibcodes for this data set
        """
        self.path = os.path.dirname(spec_files)+'/'
        self.refs = bibcode
        self.wavelength_range = (0,1E10)
        
        # Get list of spectral intensity files
        files = glob(spec_files)
        if not files:
            print('No files match',spec_files,'.')
            return
    
        # Parse the filenames and make a table for the grid
        t = [int(f.split('/')[-1][3:8]) for f in files]
        g = [float(f.split('/')[-1][9:13]) for f in files]
        m = [float(f.split('/')[-1][13:17]) for f in files]
        f = [f.split('/')[-1] for f in files]
        self.data = at.Table([t,g,m,f], 
                             names=['Teff','logg','FeH','filename'])
        
    def trim(self, teff_range=(2300,2800), logg_range=(4.5,6), 
             FeH_range=(-0.5,0.5), wavelength_range=(1.1,1.7)):
        """
        Trims the model grid by the given ranges in effective temperature,
        surface gravity, metallicity, and wavelength
    
        Parameters
        ----------
        teff_range: array-like
            The lower and upper inclusive bounds for the effective
            temperature (K)
        logg_range: array-like
            The lower and upper inclusive bounds for the logarithm of the
            surface gravity (dex)
        FeH_range: array-like
            The lower and upper inclusive bounds for the logarithm of the
            ratio of the metallicity and solar metallicity (dex)
        wavelength_range: array-like
            The lower and upper inclusive bounds for the wavelength (microns)
    
        Example
        -------
        ModelGrid.trim([2300,2800], [4.5,6], [-0.5,0.5], [1.1,1.7])
        """
        # Make a copy of the grid
        grid = self.data.copy()
        self.wavelength_range = wavelength_range
        
        # Filter grid by given parameters
        self.data = grid[[(grid['Teff']>=teff_range[0])&
                          (grid['Teff']<=teff_range[1])&
                          (grid['logg']>=logg_range[0])&
                          (grid['logg']<=logg_range[1])&
                          (grid['FeH']>=FeH_range[0])&
                          (grid['FeH']<=FeH_range[1])]]
    
        # Print a summary of the returned grid
        print('{}/{}'.format(len(self.data),len(grid)),
              'spectra in parameter range',
              'Teff:', teff_range, ', logg:',logg_range,
              ', FeH:', FeH_range, ', wavelength:', wavelength_range)
        
        # Clear the grid copy from memory
        del grid
              
    def get(self, teff, logg, FeH, verbose=False):
        """
        Retrieve the wavelength, flux, and effective radius 
        for the spectrum of the given parameters
        
        Parameters
        ----------
        teff: int
            The effective temperature (K)
        logg: float
            The logarithm of the surface gravity (dex)
        FeH: float
            The logarithm of the ratio of the metallicity 
            and solar metallicity (dex)
        verbose: bool
            Print some information about the extracted spectrum
        
        Returns
        -------
        spectrum: np.ndarray
            An array of the wavelength and flux
        """
        # Get the filepath
        filepath = self.path+str(self.data[[(self.data['Teff']==2300)&
                                (self.data['logg']==5.0)&
                                (self.data['FeH']==0.5)]]
                                ['filename'][0])
                                
        # Open the FITS file
        HDU = fits.open(filepath)
        if verbose:
            HDU.info()
            
        # Get the flux and mu arrays
        raw_flux = HDU[0].data
        mu = HDU[1].data
        
        # Extract starting wavelength, wavelength increment, and
        # effective radius [cm] from FITS header
        w0 = HDU[0].header.get('CRVAL1')
        dw = HDU[0].header.get('CDELT1')
        radius = HDU[0].header.get('PHXREFF')
        
        # Close the FITS file
        HDU.close()
        
        # Construct full wavelength scale and convert
        # from Angstroms to microns
        raw_wave = (w0 + dw*np.arange(len(raw_flux[0])))/1E4
        
        # Trim the wavelength and flux arrays
        idx = np.where(np.logical_and(raw_wave>=self.wavelength_range[0],
                                      raw_wave<=self.wavelength_range[1]))
        flux = raw_flux[:,idx]
        wave = raw_wave[idx]
        
        return wave, flux, mu, radius
        
def calculate_ldc(model_grid, orders, mu_min=0.02):
    """
    This function calculates the limb darkening coefficients for a given 
    grid of synthetic spectra
    
    Reference for limb-darkening laws:
        http://www.astro.ex.ac.uk/people/sing/David_Sing/Limb_Darkening.html
    
    Parameters
    ----------
    model_grid: ModelGrid object
        The grid of synthetic spectra from which the coefficients will
        be calculated
    orders: int
        The polynomial order, i.e. the number of coefficients
    mu_min: float
        The minimum mu value to consider
    
    Returns
    -------
    
    
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
        muz = np.interp(mu, ld, 0,01)
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
        ldc0 = replicate(1d0,npar) # ?????????????????????
        if orders==2:
            ldc[:,t_idx,g_idx,m_idx] = mpfitfun('ld2func',mu,ld,err,ldc0)
        elif orders==4:
            ldc[:,t_idx,g_idx,m_idx] = mpfitfun('ld4func',mu,ld,err,ldc0)
        else:
            print('Order number must be 2 or 4.')
            return
        
    return ldc, mu0, r_eff
    
    
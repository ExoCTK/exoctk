# Reference for model spectra: Husser, T.-O. et al. (2013)
#   http://adsabs.harvard.edu/abs/2013A%26A...553A...6H
#   http://phoenix.astro.physik.uni-goettingen.de/?page_id=73
# 
# Reference for limb-darkening laws.
#   http://www.astro.ex.ac.uk/people/sing/David_Sing/Limb_Darkening.html
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
    
    def __init__(self, spec_files, bibcode):
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
                             names=['Teff','logg','Fe/H','filename'])
        
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
        
        # Filter grid by given parameters
        self.data = grid[[(grid['Teff']>=teff_range[0])&
                          (grid['Teff']<=teff_range[1])&
                          (grid['logg']>=logg_range[0])&
                          (grid['logg']<=logg_range[1])&
                          (grid['Fe/H']>=FeH_range[0])&
                          (grid['Fe/H']<=FeH_range[1])]]
    
        # Print a summary of the returned grid
        print('{}/{}'.format(len(self.data),len(grid)),
              'spectra in parameter range',
              'Teff:', teff_range, ', logg:',logg_range,
              ', Fe/H:', FeH_range, ', wavelength:', wavelength_range)
        
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
                                (self.data['Fe/H']==0.5)]]
                                ['filename'][0])
                                
        # Open the FITS file
        HDU = fits.open(filepath)
        if verbose:
            HDU.info()
            
        # Get the flux
        flux = HDU[0].data
        
        # Extract starting wavelength, wavelength increment, and
        # effective radius [cm] from FITS header
        w0 = HDU[0].header.get('CRVAL1')
        dw = HDU[0].header.get('CDELT1')
        radius = HDU[0].header.get('PHXREFF')
        
        # Construct full wavelength scale and convert
        # from Angstroms to microns
        wave = (w0 + dw*np.arange(len(flux)))/1E4
        
        # Close the FITS files
        HDU.close()
        
        return wave, flux, radius
        
def ldc_grid(n_coeffs, minmu=0.02, *args):
    """
    This function calculates the limb darkening coefficients for a given 
    grid of synthetic spectra
    
    Parameters
    ----------
    n_coeffs: int
        The number of 
    teff_range: array-like
        The lower and upper inclusive bounds for the effective temperature (K)
    logg_range: array-like
        The lower and upper inclusive bounds for the logarithm of the 
        surface gravity (dex)
    FeH_range: array-like
        The lower and upper inclusive bounds for the logarithm of the
        ratio of the metallicity divided by solar metallicity (dex)
    wavelength_range: array-like
        The lower and upper inclusive bounds for the wavelength (microns)
    minmu: float
    
    Returns
    -------
    
    
    Example
    -------
    ldcgrid(4, [2300,2800], [4.5,6], [-0.5,0.5], [1.1,1.7])
    """
    pass

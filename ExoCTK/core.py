#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for classes and functions used across all ExoCTK subpackages
"""
from glob import glob
from astropy.io import fits
import astropy.table as at
import numpy as np
import os

class ModelGrid(object):
    """
    Creates a ModelGrid object which contains a multi-parameter
    grid of model spectra and its references
    
    Attributes
    ----------
    path: str
        The path to the directory of FITS files used to create the ModelGrid
    refs: list, str
        The references for the data contained in the ModelGrid
    wavelength_range: array-like
        The lower and upper inclusive bounds of the ModelGrid wavelength
        in microns
    data: astropy.table.Table
        The table of parameters for the ModelGrid
    
    """
    
    def __init__(self, spec_files, bibcode='2013A&A...553A...6H'):
        """
        Initializes the model grid by creating a table with a column
        for each parameter and ingests the spectra
    
        Parameters
        ----------
        spec_files: str
            The path plus wildcard filename for the FITS files of spectra
        bibcode: str, array-like (optional)
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
            Print some information about the spectrum
        
        Returns
        -------
        list
            A list of arrays of the wavelength, flux, and 
            mu values and the effective radius for the given model
        
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
        idx, = np.where(np.logical_and(raw_wave>=self.wavelength_range[0],
                                      raw_wave<=self.wavelength_range[1]))
        flux = raw_flux[:,idx]
        wave = raw_wave[idx]
        
        return wave, flux, mu, radius

class References(object):
    """
    Creates and manages a References object to track references 
    within an ExoCTK user session
    
    Attributes
    ----------
    bibfile: str
        The path to the bibtex file from which the references will be read
    bibtex: dict
        The dictionary of bibtex entries
    refs: list
        The list of bibcodes saved during the user session
    """
    def __init__(self, bibfile):
        """
        Initializes an empty References object which points to a
        .bib file
    
        Parameters
        ----------
        bibfile: str
          The path to the bibtex file from which the references will be read
    
        """
        # Attributes for the filepath, file contents, and references
        self.bibfile = bibfile
        self.bibtex = pickle.load(open(bibfile, 'rb'))
        self.refs = []

    def add(self, bibcode):
        """
        Adds a bibcode to the References object
        
        Parameters
        ----------
        bibcode: str
            The unique compact identifier for the reference to be added
        
        """
        # Check that the bibcode is i the bibtex file
        if bibcode in self.bibtex:
            self.refs += bibcode
            print(bibcode,'added to list of references.')
        
        # Suggest adding it to the bibfile
        else:
            print(bibcode,'not in bibfile at',self.bibfile)
            print('Add the bibtex entry to',self.bibfile,'and try agin.')
            
    def remove(self, bibcode):
        """
        Removes a bibcode from the References object
        
        Parameters
        ----------
        bibcode: str
            The unique compact identifier for the reference to be removed
        
        """
        # Check that the bibcode is i the bibtex file
        if bibcode in self.bibtex:
            self.refs = [r for r in self.refs if r!=bibcode]
            print(bibcode,'removed from list of references.')
        
        # Suggest adding it to the bibfile
        else:
            print(bibcode,'not in bibfile at',self.bibfile)
                
    def write(self, filepath=''):
        """
        Write the .bib file
        
        Parameters
        ----------
        filepath: str
            The directory to which the .bib file should be written.
            Can be an existing .bib file or just a directory
        
        """
        # Use existing .bib file or create new one
        bibfile = filepath if filepath.endswith('.bib') else filepath+'biblio.bib'
        
        # Iterate through the references and write the relevant bibtex to file
        for bibcode in self.refs:
            with open(bibfile, 'a') as bib:
                bib.write(self.bibtex[bibcode])

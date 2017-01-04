#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for classes and functions used across all ExoCTK subpackages
"""
from glob import glob
from astropy.io import fits
import astropy.table as at
import astropy.io.votable as vo
import numpy as np
import urllib
import os

class Filter(object):
    """
    Creates a Filters object to store photometric filter profiles
    
    Attributes
    ----------
    
    """
    def __init__(self, band, filter_directory='./ExoCTK/filters/'):
        """
        Creates a list of all the filters in the filter_directory
        """        
        # Get list of filters
        filters = glob(filter_directory+'*')
        filepath = filter_directory+band
        
        # If the filter is missing, ask what to do
        if filepath not in filters:
        
            print('No filters match',filter_directory+band)
            dl = input('Would you like me to download it? [y/n] ')
            
            if dl.lower()=='y':
                
                # Prompt for new filter
                print('\nA full list of available filters from the\n'\
                      'SVO Filter Profile Service can be found at\n'\
                      'http://svo2.cab.inta-csic.es/theory/fps3/\n')
                band = input('Enter the band name to retrieve (e.g. 2MASS/2MASS.J): ')
                
                # Download the XML (VOTable) file
                baseURL = 'http://svo2.cab.inta-csic.es/svo/theory/fps/fps.php?ID='
                filepath = filter_directory+os.path.basename(band)
                _ = urllib.request.urlretrieve(baseURL+band, filepath)
                
                # Print the new filepath
                print('Band stored as',filepath)
            
        # Try to read filter info
        try:

            # Parse the XML file
            vot = vo.parse_single_table(filepath)
            self.rsr = np.array([list(i) for i in vot.array]).T
            
            # Convert to microns
            self.rsr *= np.array([[0.0001],[1.]])
            
            # Store PARAM dict as attribute
            self.params = {}
            params = [str(p).split() for p in vot.params]
            
            for p in params:
                
                # Extract the key/value pairs
                key = p[1].split('"')[1]
                val = p[-1].split('"')[1]
                
                # Do some formatting
                if p[2].split('"')[1]=='float'\
                or p[3].split('"')[1]=='float':
                    val = float(val)
                
                else:
                    val = val.replace('b&apos;','')\
                             .replace('&apos','')\
                             .replace('&amp;','&')\
                             .strip(';')
                
                # Garbage
                if key!='Description':
                    self.params[key] = val
            
            # Create some attributes
            self.path = filepath
            self.refs = []

        # If empty, delete XML file
        except:

            print('No filter named',band)
            if os.path.isfile(filepath):
                os.remove(filepath)

            return
        
                 
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
    teff_rng: tuple
        The range of effective temperatures [K]
    logg_rng: tuple
        The range of surface gravities [dex]
    FeH_rng: tuple
        The range of metalicities [dex]
    wave_rng: array-like
        The wavelength range of the models [um]
    n_bins: int
        The number of bins for the ModelGrid wavelength array
    data: astropy.table.Table
        The table of parameters for the ModelGrid
    
    """
    def __init__(self, model_directory, bibcode='2013A&A...553A...6H',
                 names={'Teff':'PHXTEFF', 'logg':'PHXLOGG',
                       'FeH':'PHXM_H', 'mass':'PHXMASS',
                       'r_eff':'PHXREFF', 'Lbol':'PHXLUM'}):
        """
        Initializes the model grid by creating a table with a column
        for each parameter and ingests the spectra
    
        Parameters
        ----------
        model_directory: str
            The path to the directory of FITS files of spectra,
            which may include a filename with a wildcard caharacter
        bibcode: str, array-like (optional)
            The bibcode or list of bibcodes for this data set
        names: dict (optional)
            A dictionary to rename the table columns. The Phoenix
            model keywords are given as an example
        """
        # Make sure we can use glob if a directory 
        # is given without a wildcard
        if '*' not in model_directory:
            model_directory += '*'
        
        # Create some attributes
        self.path = os.path.dirname(model_directory)+'/'
        self.refs = bibcode
        self.wave_rng = (0,40)
        self.n_bins = 1E10
        
        # Get list of spectral intensity files
        files = glob(model_directory)
        filenames = []
        if not files:
            print('No files match',model_directory,'.')
            return
        
        # Parse the FITS headers
        vals, dtypes = [], []
        for f in files:
            if f.endswith('.fits'):
                try:
                    header = fits.getheader(f)
                    keys = np.array(header.cards).T[0]
                    dtypes = [type(i[1]) for i in header.cards]
                    v = [header.get(k) for k in keys]
                    vals.append(list(v))
                    filenames.append(f.split('/')[-1])
                except:
                    print(f,'could not be read into the model grid.')
            
        # Fix data types and make the table
        dtypes = [str if d==bool else d for d in dtypes]
        table = at.Table(np.array(vals), names=keys, dtype=dtypes)
        
        # Add the filenames as a column
        table['filename'] = filenames
        
        # Rename any columns
        for new,old in names.items():
            try:
                table.rename_column(old, new)
            except:
                print('No column named',old)
        
        # Remove columns where the values are all the same
        # and store value as attribute instead
        for n in table.colnames:
            val = table[n][0]
            if list(table[n]).count(val) == len(table[n]):
                setattr(self, n, val)
                table.remove_column(n)
        
        # Store the table in the data attribute
        self.data = table
        
        # Store the parameter ranges
        self.Teff_rng = (min(table['Teff']),max(table['Teff']))
        self.logg_rng = (min(table['logg']),max(table['logg']))
        self.FeH_rng = (min(table['FeH']),max(table['FeH']))
        
    def get(self, Teff, logg, FeH, verbose=False):
        """
        Retrieve the wavelength, flux, and effective radius 
        for the spectrum of the given parameters
        
        Parameters
        ----------
        Teff: int
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
        # Get the row index and filepath
        try: 
            row, = np.where((self.data['Teff']==Teff)
                          & (self.data['logg']==logg)
                          & (self.data['FeH']==FeH))[0]
        except ValueError:
            print('Teff:', Teff, ' logg:', logg, ' FeH:', FeH, 
                  ' model not in grid.')
            return
        filepath = self.path+str(self.data[row]['filename'])
        
        # Get the flux, mu, and abundance arrays
        raw_flux = fits.getdata(filepath, 0)
        mu = fits.getdata(filepath, 1)
        #abund = fits.getdata(filepath, 2)
        
        # Construct full wavelength scale and convert to microns
        raw_wave = (self.CRVAL1+self.CDELT1*np.arange(len(raw_flux[0])))/1E4
        
        # Trim the wavelength and flux arrays
        idx, = np.where(np.logical_and(raw_wave>=self.wave_rng[0],
                                      raw_wave<=self.wave_rng[1]))
        flux = raw_flux[:,idx]
        wave = raw_wave[idx]
        
        # Make a dictionary of parameters
        # This should really be a core.Spectrum() object!
        spec_dict = dict(zip(self.data.colnames, self.data[row].as_void()))
        spec_dict['wave'] = wave
        
        # Bin the spectrum if necessary
        if self.n_bins>0 and self.n_bins<len(wave):
            pass
            
        spec_dict['flux'] = flux
        spec_dict['mu'] = mu
        #spec_dict['abund'] = abund
        
        return spec_dict
        
    def customize(self, Teff_rng=(0,1E4), logg_rng=(0,6), 
                  FeH_rng=(-3,3), wave_rng=(0,40), 
                  n_bins=''):
        """
        Trims the model grid by the given ranges in effective temperature,
        surface gravity, and metallicity. Also sets the wavelength range
        and number of bins for retrieved model spectra.
        
        Parameters
        ----------
        Teff_rng: array-like
            The lower and upper inclusive bounds for the effective
            temperature (K)
        logg_rng: array-like
            The lower and upper inclusive bounds for the logarithm of the
            surface gravity (dex)
        FeH_rng: array-like
            The lower and upper inclusive bounds for the logarithm of the
            ratio of the metallicity and solar metallicity (dex)
        wave_rng: array-like
            The lower and upper inclusive bounds for the wavelength (microns)
        n_bins: int
            The number of bins for the wavelength axis
        
        """
        # Make a copy of the grid
        grid = self.data.copy()
        self.wave_rng = wave_rng
        self.n_bins = n_bins or self.n_bins
        
        # Filter grid by given parameters
        self.data = grid[[(grid['Teff']>=Teff_rng[0])
                         & (grid['Teff']<=Teff_rng[1])
                         & (grid['logg']>=logg_rng[0])
                         & (grid['logg']<=logg_rng[1])
                         & (grid['FeH']>=FeH_rng[0])
                         & (grid['FeH']<=FeH_rng[1])]]
        
        # Print a summary of the returned grid
        print('{}/{}'.format(len(self.data),len(grid)),
              'spectra in parameter range',
              'Teff:', Teff_rng, ', logg:',logg_rng,
              ', FeH:', FeH_rng, ', wavelength:', wave_rng)
        
        # Clear the grid copy from memory
        del grid
              
              
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
        # Check that the bibcode is in the bibtex file
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
        # Check that the bibcode is in the bibtex file
        if bibcode in self.bibtex:
            self.refs = [r for r in self.refs if r!=bibcode]
            print(bibcode,'removed from list of references.')
        
        # Nothing to remove!
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

def writeFITS(filename, extensions, headers=()):
    '''
    Write some data to a new FITS file
    
    Parameters
    ----------
    filename: str
        The filename of the output FITS file
    extensions: dict
        The extension name and associated data to include
        in the file
    headers: array-like
        The (keyword,value,comment) groups for the PRIMARY
        header extension
        
    '''
    # Write the arrays to a FITS file
    prihdu = fits.PrimaryHDU()
    prihdu.name = 'PRIMARY'
    hdulist = fits.HDUList([prihdu])
    
    # Write the header to the PRIMARY HDU
    hdulist['PRIMARY'].header.extend(headers, end=True)
    
    # Write the data to the HDU
    for k,v in extensions.items():
        hdulist.append(fits.ImageHDU(data=v, name=k))

    # Write the file
    hdulist.writeto(filename, clobber=True)
    hdulist.close()
    
    # Insert END card to prevent header error
    #hdulist[0].header.tofile(filename, endcard=True, clobber=True)
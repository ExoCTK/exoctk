#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for classes and functions used across all ExoCTK subpackages
"""
from __future__ import print_function

from glob import glob
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from scipy.interpolate import splmake, spleval
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import zoom
from functools import partial
import multiprocessing
import bibtexparser as bt
import astropy.table as at
import astropy.io.votable as vo
import astropy.io.ascii as ii
import astropy.units as q
import matplotlib.pyplot as plt
import pkg_resources
import pickle
import warnings
import numpy as np
import urllib
import os
import time
import h5py

warnings.simplefilter('ignore', category=AstropyWarning)

def interp_flux(mu, flux, params, values):
    """
    Interpolate a cube of synthetic spectra for a
    given index of mu
    
    Parameters
    ----------
    mu: int
        The index of the (Teff, logg, FeH, *mu*, wavelength)
        data cube to interpolate
    flux: np.ndarray
        The 5D data array
    params: list
        A list of each free parameter range
    values: list
        A list of each free parameter values
    
    Returns
    -------
    tu
        The array of new flux values
    """
    # Iterate over each wavelength (-1 index of flux array)
    l = flux.shape[-1]
    flx = np.zeros(l)
    generators = []
    for lam in range(l):
        interp_f = RegularGridInterpolator(params, flux[:,:,:,mu,lam])
        f, = interp_f(values)
        
        flx[lam] = f
        generators.append(interp_f)
    
    return flx, generators

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
    inv_file: str
        An inventory file to more quickly load the database
    
    """
    def __init__(self, model_directory, bibcode='2013A&A...553A...6H',
                 names={'Teff':'PHXTEFF', 'logg':'PHXLOGG',
                       'FeH':'PHXM_H', 'mass':'PHXMASS',
                       'r_eff':'PHXREFF', 'Lbol':'PHXLUM'}, 
                 resolution='', wl_units=q.um, **kwargs):
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
        resolution: int (optional)
            The desired wavelength resolution (lambda/d_lambda) 
            of the grid spectra
        wl_units: astropy.units.quantity
        """
        # Make sure we can use glob if a directory 
        # is given without a wildcard
        if '*' not in model_directory:
            model_directory += '*'
            
        # Check for a precomputed pickle of this ModelGrid
        model_grid = ''
        if model_directory.endswith('/*'):
            # Location of model_grid pickle
            file = model_directory.replace('*','model_grid.p')
            
            try:
                model_grid = pickle.load(open(file, 'rb'))
            except:
                pass
        
        # Instantiate the precomputed model grid
        if model_grid:
            
            for k,v in vars(model_grid).items():
                setattr(self, k, v)
            
            self.flux_file = self.path+'model_grid_flux.hdf5'
            self.flux = ''
            self.wavelength = ''
            self.r_eff = ''
            self.mu = ''
            
            del model_grid
            
        # Or compute it from scratch
        else:
            
            # Print update...
            if model_directory.endswith('/*'):
                print("Indexing models. Loading this model grid will be MUCH faster next time!")
            
            # Create some attributes
            self.path = os.path.dirname(model_directory)+'/'
            self.refs = ''
            self.wave_rng = (0,40)
            self.flux_file = self.path+'model_grid_flux.hdf5'
            self.flux = ''
            self.wavelength = ''
            self.r_eff = ''
            self.mu = ''
            
            # Save the refs to a References() object
            if bibcode:
                if isinstance(bibcode, (list,tuple)):
                    pass
                elif bibcode and isinstance(bibcode, str):
                    bibcode = [bibcode]
                else:
                    pass
                    
                self.refs = bibcode
                # _check_for_ref_object()
            
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
                        vals.append([header.get(k) for k in keys])
                        filenames.append(f.split('/')[-1])
                    except:
                        print(f,'could not be read into the model grid.')
                        
            # Fix data types, trim extraneous values, and make the table
            dtypes = [str if d==bool else d for d in dtypes]
            vals = [v[:len(dtypes)] for v in vals]
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
                if list(table[n]).count(val) == len(table[n])\
                and n not in ['Teff','logg','FeH']:
                    setattr(self, n, val)
                    table.remove_column(n)
                    
            # Store the table in the data attribute
            self.data = table
            
            # Store the parameter ranges
            self.Teff_vals = np.asarray(np.unique(table['Teff']))
            self.logg_vals = np.asarray(np.unique(table['logg']))
            self.FeH_vals = np.asarray(np.unique(table['FeH']))
            
            # Write an inventory file to this directory for future table loads
            if model_directory.endswith('/*'):
                self.file = file
                try:
                    pickle.dump(self, open(self.file, 'wb'))
                except IOError:
                    print('Could not write model grid to',self.file)
                    
        # Print something
        print(len(self.data),'models loaded from',self.path)
        
        # In case no filter is used
        self.n_bins = 1
        
        # Set the wavelength_units
        self.wl_units = q.AA
        if wl_units:
            self.set_units(wl_units)
        else:
            self.const = 1
            
        # Save the desired resolution
        self.resolution = resolution
        
        # Customize from the get-go
        if kwargs:
            self.customize(**kwargs)
            
    def get(self, Teff, logg, FeH, resolution='', interp=True):
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
        resolution: int (optional)
            The desired wavelength resolution (lambda/d_lambda) 
        interp: bool
            Interpolate the model if possible
        
        Returns
        -------
        dict
            A dictionary of arrays of the wavelength, flux, and 
            mu values and the effective radius for the given model
        
        """
        # See if the model with the desired parameters is witin the grid
        in_grid = all([(Teff>=min(self.Teff_vals))&
                       (Teff<=max(self.Teff_vals))&
                       (logg>=min(self.logg_vals))&
                       (logg<=max(self.logg_vals))&
                       (FeH>=min(self.FeH_vals))&
                       (FeH<=max(self.FeH_vals))])
                       
        if in_grid:
            
            # See if the model with the desired parameters is a true grid point
            on_grid = self.data[[(self.data['Teff']==Teff)&
                                 (self.data['logg']==logg)&
                                 (self.data['FeH']==FeH)]]\
                                 in self.data
            
            # Grab the data if the point is on the grid
            if on_grid:
                
                # Get the row index and filepath
                row, = np.where((self.data['Teff']==Teff)
                              & (self.data['logg']==logg)
                              & (self.data['FeH']==FeH))[0]
                              
                filepath = self.path+str(self.data[row]['filename'])
                
                # Get the flux, mu, and abundance arrays
                raw_flux = fits.getdata(filepath, 0)
                mu = fits.getdata(filepath, 1)
                #abund = fits.getdata(filepath, 2)
                
                # Construct full wavelength scale and convert to microns
                if self.CRVAL1=='-':
                    # Try to get data from WAVELENGTH extension...
                    raw_wave = np.array(fits.getdata(filepath, ext=-1)).squeeze()
                else:
                    # ...or try to generate it
                    l = len(raw_flux[0])
                    raw_wave = np.array(self.CRVAL1+self.CDELT1*np.arange(l)).squeeze()
                    
                # Convert from A to desired units
                raw_wave *= self.const
                
                # Trim the wavelength and flux arrays
                idx, = np.where(np.logical_and(raw_wave>=self.wave_rng[0],
                                              raw_wave<=self.wave_rng[1]))
                flux = raw_flux[:,idx]
                wave = raw_wave[idx]
                
                # Bin the spectrum if necessary
                if resolution or self.resolution:
                    
                    # Calculate zoom
                    z = _calc_zoom(resolution or self.resolution, wave)
                    wave = zoom(wave, z)
                    flux = zoom(flux, (1, z))
                    
                # Make a dictionary of parameters
                # This should really be a core.Spectrum() object!
                spec_dict = dict(zip(self.data.colnames, self.data[row].as_void()))
                spec_dict['wave'] = wave
                spec_dict['flux'] = flux
                spec_dict['mu'] = mu
                spec_dict['r_eff'] = ''
                #spec_dict['abund'] = abund
                
            # If not on the grid, interpolate to it
            else:
                # Call grid_interp method
                if interp:
                    spec_dict = self.grid_interp(Teff, logg, FeH)
                else:
                    return
                    
            return spec_dict
            
        else:
            print('Teff:', Teff, ' logg:', logg, ' FeH:', FeH, 
                  ' model not in grid.')
            return
    
    def grid_interp(self, Teff, logg, FeH, plot=False):
        """
        Interpolate the grid to the desired parameters
        
        Parameters
        ----------
        Teff: int
            The effective temperature (K)
        logg: float
            The logarithm of the surface gravity (dex)
        FeH: float
            The logarithm of the ratio of the metallicity 
            and solar metallicity (dex)
        plot: bool
            Plot the interpolated spectrum along
            with the 8 neighboring grid spectra
        
        Returns
        -------
        dict
            A dictionary of arrays of the wavelength, flux, and 
            mu values and the effective radius for the given model
        """
        # Load the fluxes
        if isinstance(self.flux,str):
            self.load_flux()
            
        # Get the flux array
        flux = self.flux.copy()
        
        # Get the interpolable parameters
        params, values = [], []
        for p,v in zip([self.Teff_vals, self.logg_vals, self.FeH_vals],
                       [Teff, logg, FeH]):
            if len(p)>1:
                params.append(p)
                values.append(v)
        values = np.asarray(values)
        label = '{}/{}/{}'.format(Teff,logg,FeH)
        
        try:
            # Interpolate flux values at each wavelength
            # using a pool for multiple processes
            print('Interpolating grid point [{}]...'.format(label))
            processes = 4
            mu_index = range(flux.shape[-2])
            start = time.time()
            pool = multiprocessing.Pool(processes)
            func = partial(interp_flux, flux=flux, params=params, values=values)
            new_flux, generators = zip(*pool.map(func, mu_index))
            pool.close()
            pool.join()
            
            # Clean up and time of execution
            new_flux = np.asarray(new_flux)
            generators = np.asarray(generators)
            print('Run time in seconds: ', time.time()-start)
            
            # Interpolate mu value
            interp_mu = RegularGridInterpolator(params, self.mu)
            mu = interp_mu(np.array(values)).squeeze()
            
            # Interpolate r_eff value
            interp_r = RegularGridInterpolator(params, self.r_eff)
            r_eff = interp_r(np.array(values)).squeeze()
            
            # Make a dictionary to return
            grid_point = {'Teff':Teff, 'logg':logg, 'FeH':FeH,
                          'mu': mu, 'r_eff': r_eff,
                          'flux':new_flux, 'wave':self.wavelength,
                          'generators':generators}
                          
            return grid_point
            
        except IOError:
            print('Grid too sparse. Could not interpolate.')
            return
            
    def load_flux(self, reset=False):
        """
        Retrieve the flux arrays for all models 
        and load into the ModelGrid.array attribute
        with shape (Teff, logg, FeH, mu, wavelength)
        """
        if reset:
            
            # Delete the old file and clear the flux attribute
            if os.path.isfile(self.flux_file):
                os.remove(self.flux_file)
            self.flux = ''
            
        if isinstance(self.flux,str):
            
            print('Loading flux into table...')
            
            if os.path.isfile(self.flux_file):
                
                # Load the flux from the HDF5 file
                f = h5py.File(self.flux_file, "r")
                self.flux = f['flux'][:]
                f.close()
                
            else:
                
                # Get array dimensions
                T, G, M = self.Teff_vals, self.logg_vals, self.FeH_vals
                shp = [len(T),len(G),len(M)]
                n, N = 1, np.prod(shp)
                
                # Iterate through rows
                for nt,teff in enumerate(T):
                    for ng,logg in enumerate(G):
                        for nm,feh in enumerate(M):
                            
                            try:
                                
                                # Retrieve flux using the `get()` method
                                d = self.get(teff, logg, feh, interp=False)
                                
                                if d:
                                    
                                    # Make sure arrays exist
                                    if isinstance(self.flux,str):
                                        self.flux = np.zeros(shp+list(d['flux'].shape))
                                    if isinstance(self.r_eff,str):
                                        self.r_eff = np.zeros(shp)
                                    if isinstance(self.mu,str):
                                        self.mu = np.zeros(shp+list(d['mu'].shape))
                                        
                                    # Add data to respective arrays
                                    self.flux[nt,ng,nm] = d['flux']
                                    self.r_eff[nt,ng,nm] = d['r_eff'] or np.nan
                                    self.mu[nt,ng,nm] = d['mu'].squeeze()
                                    
                                    # Get the wavelength array
                                    if isinstance(self.wavelength,str):
                                        self.wavelength = d['wave']
                                        
                                    # Garbage collection
                                    del d
                                    
                                    # Print update
                                    n += 1
                                    print("{:.2f} percent complete.".format(n*100./N), end='\r')
                                    
                            except IOError:
                                # No model computed so reduce total
                                N -= 1
                                
                # Load the flux into an HDF5 file
                f = h5py.File(self.flux_file, "w")
                dset = f.create_dataset('flux', data=self.flux)
                f.close()
                del dset
                print("100.00 percent complete!", end='\n')
                
        else:
            print('Data already loaded.')
            
    def customize(self, Teff_rng=(2300,8000), logg_rng=(0,6), 
                  FeH_rng=(-2,1), wave_rng=(0,40), n_bins=''):
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
        
        # Do nothing if he cut leaves the grid empty
        if len(self.data)==0:
            self.data = grid
            print('The given parameter ranges would leave 0 models in the grid.')
            print('The model grid has not been updated. Please try again.')
            return
            
        # Update the wavelength and flux attributes
        if isinstance(self.wavelength,np.ndarray):
            w = self.wavelength
            W_idx, = np.where((w>=wave_rng[0])&(w<=wave_rng[1]))
            T_idx, = np.where((self.Teff_vals>=Teff_rng[0])&(self.Teff_vals<=Teff_rng[1]))
            G_idx, = np.where((self.logg_vals>=logg_rng[0])&(self.logg_vals<=logg_rng[1]))
            M_idx, = np.where((self.FeH_vals>=FeH_rng[0])&(self.FeH_vals<=FeH_rng[1]))
            
            # Trim arrays
            self.wavelength = w[W_idx]
            self.flux = self.flux[T_idx[0]:T_idx[-1]+1,G_idx[0]:G_idx[-1]+1,M_idx[0]:M_idx[-1]+1,:,W_idx[0]:W_idx[-1]+1]
            self.mu = self.mu[T_idx[0]:T_idx[-1]+1,G_idx[0]:G_idx[-1]+1,M_idx[0]:M_idx[-1]+1]
            self.r_eff = self.r_eff[T_idx[0]:T_idx[-1]+1,G_idx[0]:G_idx[-1]+1,M_idx[0]:M_idx[-1]+1]
        
        # Update the parameter attributes
        self.Teff_vals = np.unique(self.data['Teff'])
        self.logg_vals = np.unique(self.data['logg'])
        self.FeH_vals = np.unique(self.data['FeH'])
        
        # Reload the flux array with the new grid parameters
        self.load_flux(reset=True)
        
        # Clear the grid copy from memory
        del grid
        
    def info(self):
        """
        Print a table of info about the current ModelGrid
        """
        # Get the info from the class 
        tp = (int, bytes, bool, str, float, tuple, list, np.ndarray)
        info = [[k,str(v)] for k,v in vars(self).items() if isinstance(v, tp)]

        # Make the table
        table = at.Table(np.asarray(info).reshape(len(info),2),
                 names=['Attributes','Values'])
        
        # Sort and print
        table.sort('Attributes')
        table.pprint(max_width=-1, align=['>','<'])
        
    def reset(self):
        """
        Reset the current grid to the original state
        """
        try:
            os.remove(self.path+'model_grid_flux.hdf5')
        except:
            pass
        self.__init__(self.path)
        
    def set_units(self, wl_units=q.um):
        """
        Set the wavelength and flux units
        
        Parameters
        ----------
        wl_units: str, astropy.units.core.PrefixUnit, astropy.units.core.CompositeUnit
            The wavelength units
        """
        # Set wavelength units
        old_unit = self.wl_units
        self.wl_units = q.Unit(wl_units)
            
        # Update the wavelength
        self.const = (old_unit/self.wl_units).decompose()._scale
        
def _calc_zoom(R_f, arr):
    """
    Calculate the zoom factor required to make the given
    array into the given resolution
    
    Parameters
    ----------
    R_f: int
        The desired final resolution of the wavelength array
    arr: array-like
        The array to zoom
    """
    # Get initial resolution
    lam = arr[-1]-arr[0]
    d_lam_i = np.nanmean(np.diff(arr))
    R_i = lam/d_lam_i
    
    # Calculate zoom
    d_lam_f = lam/R_f
    z = d_lam_i/d_lam_f
    
    return z
    

def rebin_spec(spec, wavnew, oversamp=100, plot=False):
    """
    Rebin a spectrum to a new wavelength array while preserving 
    the total flux
    
    Parameters
    ----------
    spec: array-like
        The wavelength and flux to be binned
    wavenew: array-like
        The new wavelength array
        
    Returns
    -------
    np.ndarray
        The rebinned flux
    
    """
    wave, flux = spec
    nlam = len(wave)
    x0 = np.arange(nlam, dtype=float)
    x0int = np.arange((nlam-1.)*oversamp + 1., dtype=float)/oversamp
    w0int = np.interp(x0int, x0, wave)
    spec0int = np.interp(w0int, wave, flux)/oversamp

    # Set up the bin edges for down-binning
    maxdiffw1 = np.diff(wavnew).max()
    w1bins = np.concatenate(([wavnew[0]-maxdiffw1], .5*(wavnew[1::]+wavnew[0:-1]), [wavnew[-1]+maxdiffw1]))
    
    # Bin down the interpolated spectrum:
    w1bins = np.sort(w1bins)
    nbins = len(w1bins)-1
    specnew = np.zeros(nbins)
    inds2 = [[w0int.searchsorted(w1bins[ii], side='left'), w0int.searchsorted(w1bins[ii+1], side='left')] for ii in range(nbins)]

    for ii in range(nbins):
        specnew[ii] = np.sum(spec0int[inds2[ii][0]:inds2[ii][1]])
    
    if plot:
        plt.figure()
        plt.loglog(wave, flux, c='b')    
        plt.loglog(wavnew, specnew, c='r')
        
    return specnew

class References(object):
    """
    Creates and manages a References object to track references 
    within an ExoCTK user session
    
    Attributes
    ----------
    bibfile: str
        The path to the bibtex file from which the references will be read
    refs: list
        The list of bibcodes saved during the user session
    database: bibtexparser.bibdatabase.BibDatabase object
        The database of parsed bibtex entries
    bibcodes: list
        The list of all bibcodes in the database
        
    """
    def __init__(self, bibfile=''):
        """
        Initializes an empty References object which points to a
        .bib file
        
        Parameters
        ----------
        bibfile: str
          The path to the bibtex file from which the references will be read
        
        """
        bibfile = bibfile or \
            pkg_resources.resource_filename('ExoCTK', 'data/core/bibtex.bib')
        
        # Attributes for the filepath and references
        self.bibfile = bibfile
        self.refs = []
        
        # Load the bibtex into a database
        bf = open(bibfile)
        self.database = bt.load(bf)
        bf.close()
        
        # The list of all bibcodes in the bibfile
        self.bibcodes = [i['ID'] for i in self.database.entries]
        
    def add(self, bibcode):
        """
        Adds a bibcode to the References object
        
        Parameters
        ----------
        bibcode: str
            The unique compact identifier for the reference to be added
        
        """
        # Check that the bibcode is in the bibtex file
        if bibcode in self.bibcodes:
            self.refs += [bibcode]
            print(bibcode,'added to list of references.')
        
        # Suggest adding it to the bibfile
        else:
            print(bibcode,'not in bibfile at',self.bibfile)
            print('Add the bibtex entry to the file and try agin.')
            
    def remove(self, bibcode):
        """
        Removes a bibcode from the References object
        
        Parameters
        ----------
        bibcode: str
            The unique compact identifier for the reference to be removed
        
        """
        # Check that the bibcode is in the bibtex file
        if bibcode in self.bibcodes:
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

        # Create a new database instance
        final = bt.bibdatabase.BibDatabase()
        
        # Get the relevant bibtex entries
        final.entries = [d for d in self.database.entries 
                         if d['ID'] in list(set(self.refs))]
        
        # Write the bibtex to file
        with open(bibfile, 'w') as out:
            out.write(bt.bwriter.BibTexWriter().write(final))

def multiplot(rows, columns, ylabel='', xlabel='', sharey=True, sharex=True, 
              fontsize=20, figsize=(15, 7), title='', **kwargs):
    """
    Creates subplots with given number or *rows* and *columns*.
    
    Parameters
    ----------
    rows: int
        The number of rows in the figure
    columns: int
        The number of columns in the figure
    ylabel: str, list
        The shared y-label or list of y-labels for each column
    xlabel: str, list
        The shared y-label or list of y-labels for each column
    sharey: bool
        Same y-axis limits
    sharex: bool
        Same x-axis limits
    fontsize: int
        The fontsize to use throughout the figure
    figsize: tuple, list
        The (x,y) dimenstions of the figure
    title: str
        The title of the figure
    Returns
    -------
    list
        A list of the figure and axes objects for the current figure
    
    Example
    -------
    >>> fig, (ax11, ax12, ax13), (ax21, ax22, ax23) = multiplot(2, 3)
    >>> ax11.plot(x, y, label='Row 1, Col 1 Plot')
    """
    # Initialize the plot
    fig, axes = plt.subplots(rows, columns, sharey=sharey, sharex=sharex, figsize=figsize)
    plt.rc('text', usetex=True)
    plt.rc('font', size=fontsize)
    
    # Set the y-label(s)
    if ylabel:
        if isinstance(ylabel, str):
            fig.text(0.04, 0.54, ylabel, ha='center', va='center', rotation='vertical', **kwargs)
        else:
            if columns > 1:
                for a, l in zip(axes, ylabel):
                    a[0].set_ylabel(l, fontsize=fontsize, labelpad=fontsize)
            else:
                for a, l in zip(axes, ylabel):
                    a.set_ylabel(l, fontsize=fontsize, labelpad=fontsize)
    
    # Set the x-label(s)
    if xlabel:
        if isinstance(xlabel, str):
            fig.text(0.54, 0.04, xlabel, ha='center', va='center', **kwargs)
        else:
            if rows > 1:
                for a, l in zip(axes, xlabel):
                    a[0].set_xlabel(l, fontsize=fontsize, labelpad=fontsize)
            else:
                for a, l in zip(axes, xlabel):
                    a.set_xlabel(l, fontsize=fontsize, labelpad=fontsize)
    
    # Plot formatting
    plt.suptitle(title)
    plt.subplots_adjust(right=0.96, top=0.93 if title else 0.96, bottom=0.15, left=0.12, hspace=0, wspace=0)
    fig.canvas.draw()
    
    return [fig] + list(axes)

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

def smooth(x,window_len=10,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string

	Source: http://www.scipy.org/Cookbook/SignalSmooth		2009-03-13 
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    #s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    s=np.r_[2*np.median(x[0:window_len/5])-x[window_len:1:-1],x,2*np.median(x[-window_len/5:])-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def medfilt(x, window_len):
    """
    Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.
    
    Parameters
    ----------
    x: np.array
        The 1D array to smooth
    window_len: int
        The size of the smoothing window
    
    Returns
    -------
    np.ndarray
        The smoothed 1D array
    """
    # assert x.ndim == 1, "Input must be one-dimensional."
    if window_len % 2 == 0:
        print("Median filter length ("+str(window_len)+") must be odd. Adding 1.")
        window_len += 1
    window_len = int(window_len)
    k2 = int((window_len - 1)//2)
    s = np.r_[2*np.median(x[0:int(window_len/5)])-x[window_len:1:-1],x,2*np.median(x[int(-window_len/5):])-x[-1:-window_len:-1]]
    y = np.zeros((len(s), window_len), dtype=s.dtype)
    
    y[:,k2] = s
    for i in range (k2):
        j = k2 - i
        y[j:,i] = s[:-j]
        y[:j,i] = s[0]
        y[:-j,-(i+1)] = s[j:]
        y[-j:,-(i+1)] = s[-1]
    return np.median(y[window_len-1:-window_len+1], axis=1)

def find_closest(axes, points, n=1, values=False):
    """
    Find the n-neighboring elements of a given value in an array
        
    Parameters
    ----------
    axes: list, np.array
        The array(s) to search
    points: array-like, float
        The point(s) to search for
    n: int
        The number of values to the left and right of the points
    Returns
    -------
    np.ndarray
        The n-values to the left and right of 'points' in 'axes'
    """
    results = []
    if not isinstance(axes,list):
        axes = [axes]
        points = [points]
        
    for i,(axis,point) in enumerate(zip(axes,points)):
        if point>=min(axis) and point<=max(axis):
            axis = np.asarray(axis)
            idx = np.clip(axis.searchsorted(point), 1, len(axis)-1)
        
            if values:
                result = axis[max(0,idx-n):min(idx+n,len(axis))]
            else:
                result = np.arange(0,len(axis))[max(0,idx-n):min(idx+n,len(axis))].astype(int)
                
            results.append(result)
        else:
            print('Point {} outside grid.'.format(point))
            return

    return results

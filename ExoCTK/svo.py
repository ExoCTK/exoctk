#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A Python wrapper for the SVO Filter Profile Service
"""
from astropy.utils.exceptions import AstropyWarning
from glob import glob
import astropy.table as at
import astropy.io.votable as vo
import astropy.units as q
import astropy.constants as ac
import matplotlib.pyplot as plt
import warnings
import pickle
import inspect
import pkg_resources
import numpy as np
import urllib
import os
from . import utils

warnings.simplefilter('ignore', category=AstropyWarning)
WL_KEYS = ['FWHM', 'WavelengthCen', 'WavelengthEff', 'WavelengthMax',
           'WavelengthMean', 'WavelengthMin', 'WavelengthPeak', 
           'WavelengthPhot', 'WavelengthPivot', 'WidthEff',
           'wl_min', 'wl_max']

class Filter(object):
    """
    Creates a Filter object to store a photometric filter profile
    and metadata
    
    Attributes
    ----------
    path: str
        The absolute filepath for the bandpass data, an ASCII file with
        a wavelength column in Angstroms and a response column of values 
        ranging from 0 to 1
    refs: list, str
        The references for the bandpass data
    rsr: np.ndarray
        The wavelength and relative spectral response (RSR) arrays
    Band: str
        The band name
    CalibrationReference: str
        The paper detailing the calibration
    FWHM:float
        The FWHM for the filter
    Facility: str
        The telescope facility
    FilterProfileService: str
        The SVO source
    MagSys: str
        The magnitude system
    PhotCalID: str
        The calibration standard
    PhotSystem: str
        The photometric system
    ProfileReference: str
        The SVO reference
    WavelengthCen: float
        The center wavelength
    WavelengthEff: float
        The effective wavelength
    WavelengthMax: float
        The maximum wavelength
    WavelengthMean: float
        The mean wavelength
    WavelengthMin: float
        The minimum wavelength
    WavelengthPeak: float
        The peak wavelength
    WavelengthPhot: float
        The photon distribution based effective wavelength
    WavelengthPivot: float
        The wavelength pivot
    WavelengthUCD: str
        The SVO wavelength unit
    WavelengthUnit: str
        The wavelength unit
    WidthEff: float
        The effective width
    ZeroPoint: float
        The value of the zero point flux
    ZeroPointType: str
        The system of the zero point
    ZeroPointUnit: str
        The units of the zero point
    filterID: str
        The SVO filter ID
    
    """
    def __init__(self, band, filter_directory=pkg_resources.resource_filename('ExoCTK', 'data/filters/'), 
                 wl_units=q.um, zp_units=q.erg/q.s/q.cm**2/q.AA, DELETE=False, **kwargs):
        """
        Loads the bandpass data into the Filter object
        
        Parameters
        ----------
        band: str
            The bandpass filename (e.g. 2MASS.J)
        filter_directory: str
            The directory containing the filter files
        wl_units: str, astropy.units.core.PrefixUnit  (optional)
            The wavelength units
        zp_units: str, astropy.units.core.PrefixUnit  (optional)
            The zeropoint flux units
        DELETE: bool
            Delete the given filter
        """
        # If the filter_directory is an array, load it as the filter
        if isinstance(filter_directory, (list,np.ndarray)):
            
            self.raw = np.array([filter_directory[0].data,filter_directory[1]])
            self.WavelengthUnit = str(q.um)
            self.ZeroPointUnit = str(q.erg/q.s/q.cm**2/q.AA)
            x, f = self.raw
            
            # Get a spectrum of Vega
            vega = np.genfromtxt(pkg_resources.resource_filename('ExoCTK', 'data/core/vega.txt'), unpack=True)[:2]
            vega = utils.rebin_spec(vega, x)*q.erg/q.s/q.cm**2/q.AA
            self.ZeroPoint = (np.trapz((vega[1]*f).to(q.erg/q.s/q.cm**2/q.AA), x=x)/np.trapz(f, x=x)).to(q.erg/q.s/q.cm**2/q.AA).value
            
            # Calculate the filter's properties
            self.filterID = band
            self.WavelengthPeak = np.max(self.raw[0])
            self.WavelengthMin = np.interp(max(f)/100.,f[:np.where(np.diff(f)>0)[0][-1]],x[:np.where(np.diff(f)>0)[0][-1]])
            self.WavelengthMax = np.interp(max(f)/100.,f[::-1][:np.where(np.diff(f[::-1])>0)[0][-1]],x[::-1][:np.where(np.diff(f[::-1])>0)[0][-1]])
            self.WavelengthEff = np.trapz(f*x*vega, x=x)/np.trapz(f*vega, x=x)
            self.WavelengthMean = np.trapz(f*x, x=x)/np.trapz(f, x=x)
            self.WidthEff = np.trapz(f*x, x=x)
            self.WavelengthPivot = np.sqrt(np.trapz(f*x, x=x)/np.trapz(f/x, x=x))
            self.WavelengthPhot = np.trapz(f*vega*x**2, x=x)/np.trapz(f*vega*x, x=x)

            # Fix these two:
            self.WavelengthCen = self.WavelengthMean
            self.FWHM = self.WidthEff
            
            # Add missing attributes
            self.rsr = self.raw.copy()
            self.path = ''
            self.pixels_per_bin = self.rsr.shape[-1]
            self.n_pixels = self.rsr.shape[-1]
            self.n_bins = 1
            self.wl_min = self.WavelengthMin
            self.wl_max = self.WavelengthMax
            
        else:
            
            # Check if TopHat
            if band.lower().replace('-','').replace(' ','')=='tophat':
            
                # check kwargs for limits
                wl_min = kwargs.get('wl_min')
                wl_max = kwargs.get('wl_max')
                filepath = ''
            
                if not wl_min and not wl_max:
                    print("Please provide **{'wl_min','wl_max'} to create top hat filter.")
                    return
                else:
                    # Load the filter
                    self.load_TopHat(wl_min, wl_max, kwargs.get('n_pixels',100))
                    
            else:
                
                # Get list of filters
                files = glob(filter_directory+'*')
                bands = [os.path.basename(b) for b in files]
                filepath = filter_directory+band
                
                # If the filter is missing, ask what to do
                if filepath not in files:
                    
                    print('Current filters:',
                          ', '.join(bands),
                          '\n')
                      
                    print('No filters match',filepath)
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
                    
                    else:
                        return
                    
                # Try to read filter info
                try:
                    
                    # Parse the XML file
                    vot = vo.parse_single_table(filepath)
                    self.rsr = np.array([list(i) for i in vot.array]).T
                    
                    # Parse the filter metadata
                    for p in [str(p).split() for p in vot.params]:
                        
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
                                 
                        # Set the attribute
                        if key!='Description':
                            setattr(self, key, val)
                            
                    # Create some attributes
                    self.path = filepath
                    self.pixels_per_bin = self.rsr.shape[-1]
                    self.n_pixels = self.rsr.shape[-1]
                    self.n_bins = 1
                    self.raw = self.rsr.copy()
                    self.wl_min = self.WavelengthMin
                    self.wl_max = self.WavelengthMax
                    
                # If empty, delete XML file
                except IOError:
                    
                    print('No filter named',band)
                    # if os.path.isfile(filepath):
                    #     os.remove(filepath)
                    
                    return
                    
        # Get the bin centers
        w_cen = np.nanmean(self.rsr[0])
        f_cen = np.nanmean(self.rsr[1])
        self.centers = np.asarray([[w_cen],[f_cen]])
        
        # Set the wavelength units
        if wl_units:
            self.set_wl_units(wl_units)
            
        # Set zeropoint flux units
        if zp_units!=self.ZeroPointUnit:
            self.set_zp_units(zp_units)
            
        # Get references
        try:
            self.refs = [self.CalibrationReference.split('=')[-1]]
        except:
            self.refs = []
            
        # Bin
        if kwargs:
            bwargs = {k:v for k,v in kwargs.items() if k in \
                      inspect.signature(self.bin).parameters.keys()}
            self.bin(**bwargs)
                
    def delete(self):
        """
        Delete the current filter
        """
        do_it = input('Are you sure you want to delete {} ? [y/n] '.format(self.filterID))
        if do_it.lower()=='y':
            os.remove(self.path)
            
            filters(update=True)
    
    def apply(self, spectrum, plot=False):
        """
        Apply the filter to the given spectrum
        
        Parameters
        ----------
        spectrum: array-like
            The wavelength [um] and flux of the spectrum
            to apply the filter to
        plot: bool
            Plot the original and filtered spectrum
        
        Returns
        -------
        np.ndarray
            The filtered spectrum
            
        """
        # Make into iterable arrays
        wav, flx = [np.asarray(i) for i in spectrum]
        
        # Make flux 2D
        if len(flx.shape)==1:
            flx = np.expand_dims(flx, axis=0)
        
        # Make throughput 3D
        rsr = np.copy(self.rsr)
        if len(rsr.shape)==2:
            rsr = np.expand_dims(rsr, axis=0)
        
        # Make empty filtered array
        filtered = np.zeros((rsr.shape[0],flx.shape[0],rsr.shape[2]))
        
        # Rebin the input spectra to the filter wavelength array
        # and apply the RSR curve to the spectrum
        for i,bn in enumerate(rsr):
            for j,f in enumerate(flx):
                filtered[i][j] = np.interp(bn[0], wav, f)*bn[1]
                
        if plot:
            plt.loglog(wav, flx[0])
            for n,bn in enumerate(rsr):
                plt.loglog(bn[0], filtered[n][0])
        
        del rsr, wav, flx
        
        return filtered.squeeze()
        
    def bin(self, n_bins=1, pixels_per_bin='', bin_throughput='',
            wl_min='', wl_max=''):
        """
        Break the filter up into bins and apply a throughput to each bin,
        useful for G141, G102, and other grisms
        
        Parameters
        ----------
        n_bins: int
            The number of bins to dice the throughput curve into
        n_cahnnels: int (optional)
            The number of channels per bin, which will be used 
            to calculate n_bins
        bin_throughput: array-like (optional)
            The throughput for each bin (top hat by default)
            must be of length pixels_per_bin
        wl_min: astropy.units.quantity (optional)
            The minimum wavelength to use
        wl_max: astropy.units.quantity (optional)
            The maximum wavelength to use
        """
        # Set n_bins and pixels_per_bin
        self.n_bins = 1
        self.pixels_per_bin = self.raw.shape[-1]
        self.n_pixels = self.raw.shape[-1]
        
        # Get wavelength limits
        unit = q.Unit(self.WavelengthUnit)
        if not wl_min:
            wl_min = self.wl_min*unit
        if not wl_max:
            wl_max = self.wl_max*unit
            
        # Apply wavelength unit
        wl_min = wl_min.to(unit)
        wl_max = wl_max.to(unit)
        r = self.raw
        
        # Trim the rsr by the given min and max
        self.rsr = r[:,np.logical_and(r[0]*unit>=wl_min,r[0]*unit<=wl_max)]
        print('Bandpass trimmed to',
              '{} - {}'.format(wl_min,wl_max))
              
        # Calculate the number of bins and channels
        rsr = len(self.rsr[0])
        if pixels_per_bin and isinstance(pixels_per_bin,int):
            self.pixels_per_bin = int(pixels_per_bin)
            self.n_bins = int(rsr/self.pixels_per_bin)
        elif n_bins and isinstance(n_bins,int):
            self.n_bins = int(n_bins)
            self.pixels_per_bin = int(rsr/self.n_bins)
        elif not n_bins and not pixels_per_bin \
        and isinstance(bin_throughput, (list,tuple,np.ndarray)):
            pass
        else:
            print('Please specify n_bins or pixels_per_bin as integers.')
            return
            
        print('{} bins of {} pixels each.'.format(self.n_bins,self.pixels_per_bin))
        
        # Trim throughput edges so that there are an integer number of bins
        new_len = self.n_bins*self.pixels_per_bin
        start = (rsr-new_len)//2
        self.rsr = np.copy(self.rsr[:,start:new_len+start])
        
        # Reshape the throughput array
        self.rsr = self.rsr.reshape(2,self.n_bins,self.pixels_per_bin)
        self.rsr = self.rsr.swapaxes(0,1)
        
        # Get the bin centers
        w_cen = np.nanmean(self.rsr[:,0,:], axis=1)
        f_cen = np.nanmean(self.rsr[:,1,:], axis=1)
        self.centers = np.asarray([w_cen,f_cen])
        
        # Get the bin throughput function
        if not isinstance(bin_throughput, (list,tuple,np.ndarray)):
            bin_throughput = np.ones(self.pixels_per_bin)
            
        # Make sure the shape is right
        if len(bin_throughput)==self.pixels_per_bin:
            
            # Save the attribute
            self.bin_throughput = np.asarray(bin_throughput)
            
            # Apply the bin throughput
            self.rsr[:,1] *= self.bin_throughput
            
        else:
            print('bin_throughput must be an array of length',self.pixels_per_bin)
            print('Using top hat throughput for each bin.')
                
    def plot(self):
        """
        Plot the filter
        """
        # If the filter is binned, plot each with bin centers
        try:
            for x,y in self.rsr:
                plt.plot(x, y)
            plt.plot(*self.centers, ls='None', marker='.', c='k')
            plt.plot(self.raw[0], self.raw[1], lw=6, alpha=0.1, zorder=0)
            
        # Otherwise just plot curve
        except:
            plt.plot(*self.rsr)
            
        plt.xlabel('Wavelength [{}]'.format(str(self.WavelengthUnit)))
        plt.ylabel('Throughput')
        
    def info(self, fetch=False):
        """
        Print a table of info about the current filter
        """
        # Get the info from the class 
        tp = (int, bytes, bool, str, float, tuple, list, np.ndarray)
        exclude = ['rsr', 'bin_throughput', 'raw', 'centers']
        info = [[k,str(v)] for k,v in vars(self).items() if isinstance(v, tp)
                and k not in exclude]
                
        # Make the table
        table = at.Table(np.asarray(info).reshape(len(info),2),
                 names=['Attributes','Values'])
        
        # Sort and print
        table.sort('Attributes')
        
        if fetch:
            return table
        else:
            table.pprint(max_width=-1, max_lines=-1, align=['>','<'])
        
    def set_wl_units(self, wl_units):
        """
        Set the wavelength and flux units
        
        Parameters
        ----------
        wl_units: str, astropy.units.core.PrefixUnit
            The wavelength units
        """
        # Set wavelength units
        old_unit = q.Unit(self.WavelengthUnit)
        new_unit = q.Unit(wl_units)
        for key in WL_KEYS:
            old_val = getattr(self, key)*old_unit
            setattr(self, key, round(old_val.to(new_unit).value, 5))
            
        # Update the rsr curve
        const = (old_unit/new_unit).decompose()._scale 
        self.raw[0] *= const
        
        if len(self.rsr.shape)==2:
            self.rsr[0] *= const
        else:
            self.rsr[:,0] *= const
        
        self.centers[0] *= const
        self.WavelengthUnit = str(new_unit)
        
    def set_zp_units(self, zp_units):
        """
        Set the wavelength and flux units
        
        Parameters
        ----------
        zp_units: str, astropy.units.core.PrefixUnit
            The units of the zeropoint flux density
        """
        # Set zeropoint flux units
        old_unit = q.Unit(self.ZeroPointUnit)
        new_unit = q.Unit(zp_units)
        
        f_nu = self.ZeroPoint*old_unit
        lam = self.WavelengthEff*q.Unit(self.WavelengthUnit)
        f_lam = (f_nu*ac.c/lam**2).to(new_unit)
        
        # Update the attributes curve
        self.ZeroPoint = f_lam.value
        self.ZeroPointUnit = str(new_unit)
        
    def load_TopHat(self, wl_min, wl_max, n_pixels=100):
        """
        Loads a top hat filter given wavelength min and max values
        
        Parameters
        ----------
        wl_min: astropy.units.quantity (optional)
            The minimum wavelength to use
        wl_max: astropy.units.quantity (optional)
            The maximum wavelength to use
        n_pixels: int
            The number of pixels for the filter
        """
        if not isinstance(wl_min, q.quantity.Quantity) \
        and not isinstance(wl_max, q.quantity.Quantity):
            print('Please provide an astropy.units.quantity.Quantity for wl_min and wl_max.')
            return
            
        # Get min, max, effective wavelengths and width
        self.n_pixels = n_pixels
        self.n_bins = 1
        self.WavelengthUnit = str(wl_min.unit)
        self.wl_min = wl_min.value
        self.wl_max = wl_max.value
        wl_eff = (self.wl_min+self.wl_max)/2.
        width = self.wl_max-self.wl_min
        
        # Create the RSR curve
        wave = np.linspace(self.wl_min, self.wl_max, n_pixels)
        rsr = np.ones(n_pixels)
        self.raw = np.array([wave,rsr])
        self.rsr = self.raw
        
        # Add the attributes
        self.WavelengthMin = self.wl_min
        self.WavelengthMax = self.wl_max
        self.path = ''
        self.refs = ''
        self.Band = 'Top Hat'
        self.CalibrationReference = ''
        self.FWHM = width
        self.Facility = '-'
        self.FilterProfileService = '-'
        self.MagSys = '-'
        self.PhotCalID = ''
        self.PhotSystem = ''
        self.ProfileReference = ''
        self.WavelengthCen = wl_eff
        self.WavelengthEff = wl_eff
        self.WavelengthMean = wl_eff
        self.WavelengthPeak = wl_eff
        self.WavelengthPhot = wl_eff
        self.WavelengthPivot = wl_eff
        self.WavelengthUCD = ''
        self.WidthEff = width
        self.ZeroPoint = 0
        self.ZeroPointType = ''
        self.ZeroPointUnit = 'Jy'
        self.filterID = 'Top Hat'

def filters(filter_directory=pkg_resources.resource_filename('ExoCTK', 'data/filters/'), update=False, fmt='table', **kwargs):
    """
    Get a list of the available filters
    
    Parameters
    ----------
    filter_directory: str
        The directory containing the filter relative spectral response curves
    update: bool
        Check the filter directory for new filters and generate pickle of table
    fmt: str
        The format for the returned table
        
    Returns
    -------
    list
        The list of band names
    """
    # Get the pickle path and make sure file exists
    p_path = filter_directory.split('/filters/')[0]+'/filter_list.p'
    updated = False
    if not os.path.isfile(p_path):
        os.system('touch {}'.format(p_path))
        
    if update:
        
        print('Loading filters into table...')
        
        # Get all the filters
        files = glob(filter_directory+'*')
        bands = [os.path.basename(b) for b in files]
        tables = []
        
        for band in bands:
            
            # Load the filter
            filt = Filter(band, **kwargs)
            filt.Band = band
            
            # Put metadata into table with correct dtypes
            info = filt.info(True)
            vals = [float(i) if i.replace('.','').replace('-','')\
                .replace('+','').isnumeric() else i for i in info['Values']]
            dtypes = np.array([type(i) for i in vals])
            table = at.Table(np.array([vals]), names=info['Attributes'], dtype=dtypes)
            
            tables.append(table)
            
            del filt, info, table
            
        # Write to the pickle
        with open(p_path, 'wb') as file:
            pickle.dump(at.vstack(tables), file)
        
    # Load the saved pickle
    try:
        with open(p_path, 'rb') as file:
            data = pickle.load(file)
    except:
        data = {}
    
    # Return the data 
    if data:
        
        if fmt=='dict':
            data = {r[0]:{k:r[k].value if hasattr(r[k],'unit') else r[k] \
                    for k in data.keys()[1:]} for r in data}
                    
        # Add Band as index
        data.add_index('Band')
        
        return data
        
    # Or try to generate it once
    else:
        if not updated:
            updated = True
            filters(update=True)
        else:
            print('No filters found in',filter_directory)
            
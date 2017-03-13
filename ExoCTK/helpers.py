#!/usr/bin/python
# -*- coding: latin-1 -*-
from astropy.io import fits, ascii
from shutil import copyfile
from glob import glob
from pkg_resources import resource_filename
import numpy as np
import os

def external_files():
    """
    A snippet to propagate the external files directory
    to the submodules
    """
    try:
        from ConfigParser import ConfigParser
    except ImportError:
        from configparser import ConfigParser

    conf = ConfigParser()
    conf.read(['../setup.cfg'])
    metadata = dict(conf.items('metadata'))

    return metadata.get('external_files')

def convert_ATLAS9(filepath, destination='', template=resource_filename('ExoCTK', 'data/core/ModelGrid_tmp.fits')):
    """
    Split ATLAS9 FITS files into separate files containing one Teff, log(g), and Fe/H
    
    ACES models are in [erg/s/cm2/cm] whereas ATLAS9 models are in [erg/cm2/s/hz/ster]
    
    Parameters
    ----------
    filepath: str
        The path to the ATLAS9 FITS file to convert
    destination: str
        The destination for the split files
    template: str
        The path to the FITS template file to use
    """        
    # Get all the data
    L = open(filepath).readlines()
    
    # Get the indexes of each log(g) chunk
    start = []
    for idx,l in enumerate(L):
        if l.startswith('EFF'):
            start.append(idx)
    
    # Break up into chunks
    for n,idx in enumerate(start):
        
        # Get the parameters
        h = L[idx].strip().split()
        teff = int(h[1].split('.')[0])
        logg = float(h[3][:3])
        vturb = float(h[8])
        xlen = float(h[11])
        feh = 0.
        
        # Parse the data
        try:
            end = start[n+1]
        except:
            end = -1
        data = L[idx+3:end-4]
        
        # Fix column spacing
        for n,l in enumerate(data):
            data[n] = l[:18]+' '+' '.join([l[idx:idx+6] for idx in np.arange(18,len(l),6)])
            
        cols = ['wl']+L[idx+2].strip().split()
        data = ascii.read(data, names=cols)

        # Put intensity array for each mu value in a cube
        data_cube = np.array([data[cols][n] for n in cols[1:]])[::-1]

        # mu values
        mu = list(map(float,cols[1:]))[::-1]

        # Get the wavelength and convert from nm to A
        wave = np.array(data['wl'])*10.

        # Copy the old HDU list
        logg_txt = str(abs(int(logg*10.))).zfill(2)
        feh_txt = '{}{}'.format('m' if feh<0 else 'p', str(abs(int(feh*10.))).zfill(2))
        new_file = destination+'ATLAS9_{}_{}_{}.fits'.format(teff,logg_txt,feh_txt)
        HDU = fits.open(template)
        
        # Write the new data
        HDU[0].data = data_cube
        HDU[1].data = mu
        
        # Write the new key/values
        HDU[0].header['PHXTEFF'] = teff
        HDU[0].header['PHXLOGG'] = logg
        HDU[0].header['PHXM_H'] = feh
        HDU[0].header['PHXXI_L'] = vturb
        HDU[0].header['PHXXI_M'] = vturb
        HDU[0].header['PHXXI_N'] = vturb        
        HDU[0].header['PHXEOS'] = 'ATLAS9'
        HDU[0].header['PHXMXLEN'] = xlen
        HDU[0].header['PHXREFF'] = '-'
        HDU[0].header['PHXBUILD'] = '-'
        HDU[0].header['PHXVER'] = '-'
        HDU[0].header['DATE'] = '-'
        HDU[0].header['PHXMASS'] = '-'
        HDU[0].header['PHXLUM'] = '-'
        HDU[0].header['CRVAL1'] = '-'
        HDU[0].header['CDELT1'] = '-'

        # Create a WAVELENGTH extension
        ext = fits.ImageHDU(wave)
        ext.update_ext_name('WAVELENGTH')
        HDU.append(ext)
        
        # Write the new file
        fits.HDUList(HDU).writeto(new_file, clobber=True)
        
        HDU.close()

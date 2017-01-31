#!/usr/bin/python
# -*- coding: latin-1 -*-
from astropy.io import fits

def convert_ATLAS9(filepath, destination=''):
    """
    Split ATLAS9 FITS files into separate files containing one Teff, log(g), and Fe/H
    
    Parameters
    ----------
    filepath: str
        The path to the ATLAS9 FITS file to convert
    destination: str
        The destination for the split files
    """
    # Open the fits file
    HDU = fits.open('/Users/jfilippazzo/Desktop/ckp00_3500.fits')
    D = HDU[1].data
    primary = HDU[0].header
    hdr = HDU[1].header
    
    # Get the common data/parameters
    wave = D['WAVELENGTH']
    teff = primary['TEFF']
    feh = primary['LOG_Z']
    
    # Slice data by log(g) value
    for g in ['TTYPE{}'.format(n) for n in range(3,13)]:
        
        # Get the log(g) value and data
        logg = int(hdr[g][1:])/10.
        vals = D[hdr[g]]
        
        # Make the new HDU
        new_hdu = HDU.copy()
        print(wave,vals)
        new_hdu[1].data = np.ndarray([wave,vals])
        
        # Make the new header
        new_hdr = hdr.copy()
        new_hdr['LOGG'] = logg
        new_hdu[0].header = new_hdr
        
        # Write the new file
        filename = 'ATLAS9_{}_{}_{}.fits'.format(teff,int(logg*10),int(feh*10))
        new_hdu.writeto(destination+filename, overwrite=True)
    
    # Close the file
    HDU.close()
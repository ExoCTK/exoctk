#!/usr/bin/python

"""A module for helpful code snippets.

Authors
-------

    - Joe Filippazzo

Use
---

    This module is intended to be imported and used by other python
    modules/scripts, e.g.

    ::

        from exoctk import helpers
        helpers.convert_ATLAS9()

Dependencies
------------

    - ``astropy``
    - ``numpy``
"""

from glob import glob
import os
from pkg_resources import resource_filename
from shutil import copyfile

import astropy.constants as ac
from astropy.io import fits, ascii
import astropy.units as q
import numpy as np


def external_files():
    """A snippet to propagate the external files directory to the
    to the submodules


    Returns
    -------
    ext_files : str
        The external files directory
    """
    try:
        from ConfigParser import ConfigParser
    except ImportError:
        from configparser import ConfigParser

    conf = ConfigParser()
    conf.read(['../setup.cfg'])
    metadata = dict(conf.items('metadata'))
    ext_files = metadata.get('external_files')

    return ext_files


def convert_ATLAS9(filepath, destination='', template=resource_filename('ExoCTK', 'data/core/ModelGrid_tmp.fits')):
    """Split ``ATLAS9`` FITS files into separate files containing one
    ``Teff``, ``log(g)``, and ``Fe/H``.

    ``ACES`` models are in ``[erg/s/cm2/cm]`` whereas ``ATLAS9`` models
    are in ``[erg/cm2/s/hz/ster]``.

    Parameters
    ----------
    filepath : str
        The path to the ``ATLAS9`` FITS file to convert
    destination : str
        The destination for the split files
    template : str
        The path to the FITS template file to use
    """

    # Get all the data
    L = open(filepath, encoding='utf-8').readlines()

    # Get the indexes of each log(g) chunk
    start = []
    for idx,l in enumerate(L):
        if l.startswith('TEFF'):
            start.append(idx)

    # Break up into chunks
    for n,idx in enumerate(start):

        try:

            # Get the parameters
            h = L[idx].strip().split()
            teff = int(h[1].split('.')[0])
            logg = float(h[3][:3])
            vturb = float(h[8])
            xlen = float(h[11])
            feh = float(h[6].replace('[','').replace(']',''))

            # Parse the data
            try:
                end = start[n+1]
            except:
                end = -1
            data = L[idx+3:end-4]

            # Fix column spacing
            for n,l in enumerate(data):
                data[n] = l[:19]+' '+' '.join([l[idx:idx+6] for idx in np.arange(19,len(l),6)])

            # Get cols and data
            cols = ['wl']+L[idx+2].strip().split()
            data = ascii.read(data, names=cols)

            # Put intensity array for increasing mu values in a cube
            data_cube = np.array([data[cols][n] for n in cols[1:]])[::-1]

            # Scale the flux values by the flux(mu=1) value
            data_cube[:-1] *= data_cube[-1]
            data_cube[-1] *= 1E5

            # Apply units
            data_cube = data_cube*q.erg/q.cm**2/q.s/q.steradian/q.Hz

            # mu values
            mu = list(map(float,cols[1:]))[::-1]

            # Get the wavelength and convert from nm to A
            wave = np.array(data['wl'])*q.nm.to(q.AA)

            # Convert the flux from [erg/cm2/s/hz/ster] to [erg*m/cm**2/Hz/s**2/sr]
            # by multiplying by c/lambda**2
            data_cube = data_cube*ac.c/(wave**2)

            # Convert [m/sr] to [cm-1]
            data_cube = data_cube*q.steradian/q.cm**2

            # Convert to [erg/s/cm2/cm]
            data_cube = data_cube.to(q.erg/q.s/q.cm**3)*1E16

            # Copy the old HDU list
            logg_txt = str(abs(int(logg*10.))).zfill(2)
            feh_txt = '{}{}'.format('m' if feh<0 else 'p', str(abs(int(feh*10.))).zfill(2))
            new_file = destination+'ATLAS9_{}_{}_{}.fits'.format(teff,logg_txt,feh_txt)
            HDU = fits.open(template)

            # Write the new data
            HDU[0].data = data_cube.value
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

        except:
            pass

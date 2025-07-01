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

from pkg_resources import resource_filename
import importlib
import os

import astropy.constants as ac
from astropy.io import fits, ascii
import astropy.units as q
from bokeh.plotting import figure, show
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import stpsf

from exoctk.modelgrid import ACES
from exoctk.utils import medfilt


def generate_DHS_traces(teff_values, blocking_filter='F150W2', pupil_mask='DHS_01', outdir=None, plot=False):
    """
    Generate a library of NIRCam DHS traces in the given Teff parameter space

    Example
    -------
    # This generates and saves the FITS files to the user's EXOCTK_DATA directory
    import numpy as np
    from exoctk.helpers import generate_DHS_traces
    teffs = np.arange(2400, 8000, 200)
    generate_DHS_traces(teffs)

    Parameters
    ----------
    teff_values: list[int]
        The Teff values for each trace
    blocking_filter: str
        The filter to use, ['F070W', 'F090W', 'F115W', 'F150W2', 'F200W']
    pupil_mask: str
        The pupil to use, ['DHS_01', 'DHS_02', ..., 'DHS_10']
    """
    # Set up the PSF generator
    nrc = stpsf.NIRCam()
    nrc.filter = blocking_filter
    nrc.pupil_mask = pupil_mask  # Different for 1-10?
    nrc.detector = 'NRCA1'

    # Make a PSF cube
    wave_min = 1.01
    wave_max = 2.25
    waves = np.linspace(wave_min, wave_max, 100) * 1E-6
    psf_cube = nrc.calc_datacube(waves, oversample=1)[0].data
    psf_cube = psf_cube[:, 50:-50, 50:-50]
    # psf_cube = np.rot90(psf_cube, k=1, axes=(1, 2))

    # Get wavelengths for each detector column given the dispersion scale of NIRCam
    disp_scale = 0.290 * 0.001  # um/pixel
    pixel_wavelengths = np.arange(wave_min, wave_max, disp_scale)

    # Get throughputs and interpolate to pixel_wavelengths
    with importlib.resources.files('exoctk.data.throughputs').joinpath(f'NIRCam.{blocking_filter}.CLEARP.SW.txt').open('r') as f:
        data = np.loadtxt(f)
    throughputs = np.interp(pixel_wavelengths, data[:, 0], data[:, 1])

    # 1. Define the original grid points
    x_original = np.arange(psf_cube.shape[0])
    y_original = np.arange(psf_cube.shape[1])
    z_original = np.arange(psf_cube.shape[2])

    # 2. Define the interpolation grid
    # Create a meshgrid of the new interpolation points
    X_interp, Y_interp, Z_interp = np.meshgrid(pixel_wavelengths, y_original, z_original, indexing='ij')

    # Reshape the meshgrid to be a list of points (each point is [x, y, z])
    interpolation_points = np.vstack([X_interp.ravel(), Y_interp.ravel(), Z_interp.ravel()]).T

    # 3. Create the interpolator object
    interp = RegularGridInterpolator((x_original, y_original, z_original), psf_cube)

    # 4. Perform the interpolation
    interpolated_values = interp(interpolation_points)

    # Reshape the interpolated values back to the desired output shape
    interpolated_data_cube = interpolated_values.reshape((len(pixel_wavelengths), 60, 60))

    # Load model grid once
    mg = ACES()

    # Iterate over list of Teff values
    for teff in teff_values:

        try:

            # Fetch the (smoothed) SED and interpolate to the pixel wavelengths
            model = mg.get(teff, 4.5, 0)
            flx = medfilt(model['flux'][0], 11)
            interpolated_sed = np.interp(pixel_wavelengths, model['wave'], flx)

            # Multiply by the filter throughputs and SED
            final_psf_cube = interpolated_data_cube * throughputs[:, None, None] * interpolated_sed[:, None, None]

            # Determine the dimensions of the output 2D array
            num_slices, slice_height, slice_width = final_psf_cube.shape
            output_width = slice_width + (num_slices - 1)  # Width with offset
            output_height = slice_height

            # Create the result array
            trace = np.zeros((output_height, output_width), dtype=final_psf_cube.dtype)

            # Iterate and place slices
            for i in range(num_slices):
                offset = i  # Offset to the right by i columns
                trace[:, offset:offset + slice_width] += final_psf_cube[i, :, :]

            # Save to file
            primary = fits.PrimaryHDU()
            tracehdu = fits.ImageHDU(data=trace, name='TRACE')
            wavhdu = fits.ImageHDU(data=pixel_wavelengths, name='WAV')
            thruhdu = fits.ImageHDU(data=throughputs, name='THRU')
            sedhdu = fits.ImageHDU(data=interpolated_sed, name='SED')
            hdulist = fits.HDUList([primary, tracehdu, wavhdu, thruhdu, sedhdu])

            # Write the file
            filename = f'NRCA5_DHS_{blocking_filter}_{teff}.fits'
            if outdir is None:
                outdir = f"{os.environ['EXOCTK_DATA']}exoctk_contam/traces/NRCA5_DHS_{blocking_filter}"
            filepath = os.path.join(outdir, filename)
            hdulist.writeto(filepath, overwrite=True)
            hdulist.close()

            if plot:
                trace_data = fits.getdata(filepath)
                tr = figure(width=900, height=100)
                tr.image([trace_data], x=0, y=0, dw=trace_data.shape[1], dh=trace_data.shape[0])
                show(tr)
            print(f"Created trace file at {filepath}")

        except:
            print(f"Could not create trace file for Teff={teff}")


def convert_ATLAS9(filepath, destination='', template=resource_filename('exoctk', 'data/core/ModelGrid_tmp.fits')):
    """
    Split ATLAS9 FITS files into separate files containing one Teff,
    log(g), and Fe/H

    ACES models are in [erg/s/cm2/cm] whereas ATLAS9 models are in
    [erg/cm2/s/hz/ster]

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
    for idx, l in enumerate(L):
        if l.startswith('TEFF'):
            start.append(idx)

    # Break up into chunks
    for n, idx in enumerate(start):

        try:

            # Get the parameters
            h = L[idx].strip().split()
            teff = int(h[1].split('.')[0])
            logg = float(h[3][:3])
            vturb = float(h[8])
            xlen = float(h[11])
            feh = float(h[6].replace('[', '').replace(']', ''))

            # Parse the data
            try:
                end = start[n + 1]
            except:
                end = -1
            data = L[idx + 3:end - 4]

            # Fix column spacing
            for n, l in enumerate(data):
                data[n] = l[:19] + ' ' + ' '.join([l[idx:idx + 6] for idx in np.arange(19, len(l), 6)])

            # Get cols and data
            cols = ['wl'] + L[idx + 2].strip().split()
            data = ascii.read(data, names=cols)

            # Put intensity array for increasing mu values in a cube
            data_cube = np.array([data[cols][n] for n in cols[1:]])[::-1]

            # Scale the flux values by the flux(mu=1) value
            data_cube[:-1] *= data_cube[-1]
            data_cube[-1] *= 1E5

            # Apply units
            data_cube = data_cube * q.erg / q.cm ** 2 / q.s / q.steradian / q.Hz

            # mu values
            mu = list(map(float, cols[1:]))[::-1]

            # Get the wavelength and convert from nm to A
            wave = np.array(data['wl']) * q.nm.to(q.AA)

            # Convert the flux from [erg/cm2/s/hz/ster] to [erg*m/cm**2/Hz/s**2/sr]
            # by multiplying by c/lambda**2
            data_cube = data_cube * ac.c / (wave ** 2)

            # Convert [m/sr] to [cm-1]
            data_cube = data_cube * q.steradian / q.cm ** 2

            # Convert to [erg/s/cm2/cm]
            data_cube = data_cube.to(q.erg / q.s / q.cm ** 3) * 1E16

            # Copy the old HDU list
            logg_txt = str(abs(int(logg * 10.))).zfill(2)
            feh_txt = '{}{}'.format('m' if feh < 0 else 'p', str(abs(int(feh * 10.))).zfill(2))
            new_file = destination + 'ATLAS9_{}_{}_{}.fits'.format(teff, logg_txt, feh_txt)
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

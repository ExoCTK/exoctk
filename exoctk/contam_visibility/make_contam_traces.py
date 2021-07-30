import numpy as np
import os
from astropy.io import fits
from bokeh.plotting import figure, show
from bokeh.io import output_notebook
from hotsoss.plotting import plot_frame
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation


def generate_jwst_traces(min_teff=2800, max_teff=6000, increment=100, norm_mag=10., outdir=None):
    """
    Generate the precomputed traces for a range of Teff values
    to be used by the contamination tool

    Parameters
    ----------
    min_teff: int
        The minimum Teff to calculate
    max_teff: int
        The maximum Teff to calculate
    increment: int
        The increment in Teff space to use
    norm_mag: float
        The magnitude to normalize to
    outdir: str
        The path for the generated files

    """
    modes = {'NIS_SUBSTRIP256': {'inst': 'niriss', 'mode': 'soss', 'subarray': 'substrip256'},
             'NIS_SUBSTRIP96': {'inst': 'niriss', 'mode': 'soss', 'subarray': 'substrip96'},
             'MIRIM_SLITLESSPRISM': {'inst': 'miri', 'mode': 'lrsslitless'},
             'NRCA5_GRISM256_F322W2': {'inst': 'nircam', 'mode': 'wfgrism', 'filter': 'f322w2'},
             'NRCA5_GRISM256_F444W': {'inst': 'nircam', 'mode': 'wfgrism', 'filter': 'f444w'}}

    for name, mode in modes.items():

        # Configure the instrument
        configuration = build_default_calc("jwst", mode['inst'], mode['mode'])
        configuration['configuration']['instrument']['filter'] = mode.get('filter')
        subarray = mode.get('subarray', configuration['configuration']['detector']['subarray'])
        configuration['configuration']['detector']['subarray'] = subarray

        # Set the scene
        scene = {}
        scene['position'] = {'x_offset': 0., 'y_offset': 0., 'orientation': 0., 'position_parameters': ['x_offset', 'y_offset', 'orientation']}
        scene['shape'] = {'geometry': 'point'}
        scene['spectrum'] = {'name': 'Phoenix Spectrum', 'spectrum_parameters': ['sed', 'normalization']}
        scene['spectrum']['normalization'] = {'type': 'jwst', 'bandpass': 'niriss,imaging,f480m', 'norm_flux': norm_mag, 'norm_fluxunit': 'vegamag'}

        # Perform the calculation
        for teff in np.arange(min_teff, max_teff + increment, increment):

            print("Generating {} {}...".format(name, teff))

            # Set the temperature
            scene['spectrum']['sed'] = {'sed_type': 'phoenix', 'teff': teff, 'log_g': 5.0, 'metallicity': 0.0}

            # Set the scene
            configuration['scene'][0] = scene

            # Perform calculation
            report = perform_calculation(configuration, webapp=False)
            trace = np.rot90(report['2d']['detector'])

            # Set directory for output
            if outdir is None:
                outdir = os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/traces')

            fullpath = os.path.join(outdir, name)
            if not os.path.exists(fullpath):
                os.system('mkdir {}'.format(fullpath))

            # Save the file
            hdu0 = fits.PrimaryHDU()
            hdu1 = fits.ImageHDU([trace])
            hdulist = fits.HDUList([hdu0, hdu1])
            hdulist.writeto(os.path.join(fullpath, '{}_{}.fits'.format(name, int(teff))), overwrite=True)
            print("Saved {}".format(fullpath))
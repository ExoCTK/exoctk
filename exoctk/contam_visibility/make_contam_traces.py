from copy import copy
import numpy as np
import os
from exoctk.utils import add_array_at_position
from bokeh.plotting import figure, show
from scipy.interpolate import interp1d
import stpsf
from hotsoss.plotting import plot_frame
from jwst.extract_1d.soss_extract.pastasoss import get_soss_traces
from astropy.io import fits
from svo_filters import Filter

from exoctk.contam_visibility.field_simulator import APERTURES


def make_SOSS_trace_template():
    """
    Generate a NIRISS SOSS trace template file that can be scaled for each source
    """

    # Get params from APERTURES dict
    xdim = 2048
    ydim_256 = 256

    # Initialize the STPSF instance
    niriss = stpsf.NIRISS()
    niriss.filter = 'CLEAR'
    niriss.detector = 'NIS'
    niriss.pupil_mask = 'GR700XD'

    # Make a cube of PSFs to interpolate
    nwave = 100
    wavelengths_um = np.linspace(0.6, 2.9, nwave)
    fov_pixels = 65
    oversample = 1
    cube = np.zeros((nwave, fov_pixels, fov_pixels))
    for i, wave_um in enumerate(wavelengths_um):
        hdul = niriss.calc_psf(monochromatic=wave_um * 1e-6, fov_pixels=fov_pixels, oversample=oversample)
        psf = hdul[0].data
        psf /= psf.sum()
        psf = np.rot90(psf, k=-1)
        cube[i] = psf

    psf_interp = interp1d(wavelengths_um, cube, axis=0, kind='linear', bounds_error=False, fill_value='extrapolate')

    # Make SUBSTRIP256 traces
    substrip256_traces = np.zeros((3, ydim_256 + 1, xdim))
    spectrace256_file = os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/wavecal/jwst_niriss_spectrace_0023.fits')
    for order in [1, 2, 3]:
        frame = np.zeros((ydim_256, xdim))
        _, x, y, w = get_soss_traces(245.76, order)
        thru256 = fits.getdata(spectrace256_file, ext=order)
        thru = np.interp(w, thru256['WAVELENGTH'], thru256['THROUGHPUT'])

        # Orders are not the same dispersion so normalize PSF by bin
        dw = np.abs(np.gradient(w))

        # Trace is constructed left to right, in reverse dispersion direction
        # so start with PSF for the longest possible wavelength in the cube
        w_prev = cube[-1]
        for i, (xv, yv, wv, dv, tv) in enumerate(zip(x, y, w, dw, thru)):
            try:
                psf = psf_interp(wv) * tv * dv
                w_prev = psf
            except:
                print(i, wave, 'Using previous PSF')
                psf = w_prev * tv
            frame = add_array_at_position(frame, psf, int(xv), round(yv), centered=True)

        substrip256_traces[order - 1, 1:, :] = frame

        # Save wavelengths to first row. field_simulator.get_trace pulls it out later for scaling
        substrip256_traces[order - 1, 0, 4:4 + len(w)] = w

    np.save(os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/traces/NIS_SUBSTRIP256.npy'), substrip256_traces)


def make_DHS_trace_template(aperture='NRCA5_41STRIPE1_DHS_F322W2'):
    """
    Generate a NIRCam DHS mode SOSS trace template file that can be scaled for each source

    Parameters
    ----------
    aperture: str
        The DHS short wavelength channel aperture name, ['NRCA5_41STRIPE1_DHS_F322W2', 'NRCA5_41STRIPE1_DHS_F444W']
    """

    # Get aperture params
    wavecal_file = os.path.join(os.environ['EXOCTK_DATA'], f'exoctk_contam/wavecal/{aperture}_wavecal.npy')
    all_traces = np.load(wavecal_file)[:, :, 11:-11]
    xdim, ydim = 4257, 4257
    y0, y1 = 1512, 2744

    # Container for traces
    dhs_traces = np.zeros((10, ydim, xdim))

    # Get F150W2 throughput
    f150w2 = Filter('JWST/NIRCam.F150W2')
    thru_w, thru_a = f150w2.rsr[0]

    # Make PSF cube
    nircam = stpsf.NIRCam()
    nircam.filter = 'F150W2'
    nircam.detector = 'NRCA3'
    pupils = [f'DHS_{str(n).zfill(2)}' for n in [5, 4, 3, 2, 1, 6, 7, 8, 9, 10]]

    # Make a cube of PSFs to interpolate
    nwave = 100
    wavelengths_um = np.linspace(0.9, 2.3, nwave)
    fov_pixels = 65
    oversample = 1
    for order, pupil in enumerate(pupils):
        cube = np.zeros((nwave, fov_pixels, fov_pixels))
        for i, wave_um in enumerate(wavelengths_um):
            hdul = nircam.calc_psf(monochromatic=wave_um * 1e-6, fov_pixels=fov_pixels, oversample=oversample)
            psf = hdul[0].data
            psf /= psf.sum()
            psf = np.rot90(psf, k=-1)
            cube[i] = psf

        psf_interp = interp1d(wavelengths_um, cube, axis=0, kind='linear', bounds_error=False, fill_value='extrapolate')

        # Add PSFs to frame
        x, y, w = all_traces[order]
        nircam.pupil_mask = pupils[order]
        thru = np.interp(w, thru_w, thru_a)
        frame = np.zeros((xdim, ydim))

        for i, (xv, yv, wv, tv) in enumerate(zip(x, y, w, thru)):
            try:
                psf = psf_interp(wv) * tv
                w_prev = psf
            except:
                print(i, wave, 'Using previous PSF')
                psf = w_prev * tv
            frame = add_array_at_position(frame, psf, int(xv), round(yv), centered=True)

        # Add frame to cube
        dhs_traces[order, :, :] = frame

        # Add the wavelength values for the trace to the bottom row for easy access
        dhs_traces[order, y0-1, :len(w)] = w

    # Trim detector with no signal
    dhs_traces = dhs_traces[:, y0-1:y1, :]

    # Save the traces to file
    np.save(os.path.join(os.environ['EXOCTK_DATA'], f'exoctk_contam/traces/{aperture}.npy'), dhs_traces)


def generate_pandeia_traces(min_teff=2800, max_teff=6000, increment=100, norm_mag=10., outdir=None):
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
    from pandeia.engine.calc_utils import build_default_calc
    from pandeia.engine.perform_calculation import perform_calculation

    modes = {'NIS_SUBSTRIP256': {'inst': 'niriss', 'mode': 'soss', 'subarray': 'substrip256'},
             'NIS_SUBSTRIP96': {'inst': 'niriss', 'mode': 'soss', 'subarray': 'substrip96'},
             'MIRIM_SLITLESSPRISM': {'inst': 'miri', 'mode': 'lrsslitless'},
             'NRS_S1600A1_SLIT_PRISM_CLEAR': {'inst': 'nirspec', 'mode': 'bots', 'disperser': 'prism', 'filter': 'clear'},
             'NRS_S1600A1_SLIT_G140M_F070LP': {'inst': 'nirspec', 'mode': 'bots', 'disperser': 'g140m', 'filter': 'f070lp'},
             'NRS_S1600A1_SLIT_G140M_F100LP': {'inst': 'nirspec', 'mode': 'bots', 'disperser': 'g140m', 'filter': 'f100lp'},
             'NRS_S1600A1_SLIT_G235M_F170LP': {'inst': 'nirspec', 'mode': 'bots', 'disperser': 'g235m', 'filter': 'f170lp'},
             'NRS_S1600A1_SLIT_G395M_F290LP': {'inst': 'nirspec', 'mode': 'bots', 'disperser': 'g395m', 'filter': 'f290lp'},
             'NRS_S1600A1_SLIT_G140H_F070LP': {'inst': 'nirspec', 'mode': 'bots', 'disperser': 'g140h', 'filter': 'f070lp'},
             'NRS_S1600A1_SLIT_G140H_F100LP': {'inst': 'nirspec', 'mode': 'bots', 'disperser': 'g140h', 'filter': 'f100lp'},
             'NRS_S1600A1_SLIT_G235H_F170LP': {'inst': 'nirspec', 'mode': 'bots', 'disperser': 'g235h', 'filter': 'f170lp'},
             'NRS_S1600A1_SLIT_G395H_F290LP': {'inst': 'nirspec', 'mode': 'bots', 'disperser': 'g395h', 'filter': 'f290lp'},
             'NRCA5_GRISM256_F322W2': {'inst': 'nircam', 'mode': 'wfgrism', 'filter': 'f322w2'},
             'NRCA5_GRISM256_F356W': {'inst': 'nircam', 'mode': 'wfgrism', 'filter': 'f356w'},
             'NRCA5_GRISM256_F277W': {'inst': 'nircam', 'mode': 'wfgrism', 'filter': 'f277w'},
             'NRCA5_GRISM256_F444W': {'inst': 'nircam', 'mode': 'wfgrism', 'filter': 'f444w'}}

    for name, mode in modes.items():

        # Configure the instrument
        configuration = build_default_calc("jwst", mode['inst'], mode['mode'])
        if 'filter' in mode:
            configuration['configuration']['instrument']['filter'] = mode.get('filter')
        if 'disperser' in mode:
            configuration['configuration']['instrument']['disperser'] = mode.get('disperser')
        if 'subarray' in mode:
            configuration['configuration']['detector']['subarray'] = mode.get('subarray', configuration['configuration']['detector']['subarray'])

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

            # Set directory for output
            if outdir is None:
                outdir = os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/traces')

            fullpath = os.path.join(outdir, name)
            if not os.path.exists(fullpath):
                os.system('mkdir {}'.format(fullpath))

            # Save the file
            hdu0 = fits.PrimaryHDU()
            hdu1 = fits.ImageHDU([report['2d']['detector']])
            hdulist = fits.HDUList([hdu0, hdu1])
            hdulist.writeto(os.path.join(fullpath, '{}_{}.fits'.format(name, int(teff))), overwrite=True)
            print("Saved {}".format(fullpath))
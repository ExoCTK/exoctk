#! /usr/bin/env python

"""Generate Pandeia trace assets for MIRI/LRS contamination calculations."""

import argparse
import copy
from contextlib import contextmanager
import json
import os
import re
import subprocess
from unittest.mock import patch

from astropy.io import fits
import numpy as np

from . import miri_lrs


def pandeia_provenance(require_miri=True):
    """Validate and report the active Pandeia data provenance.

    Parameters
    ----------
    require_miri : bool, optional
        Require the Pandeia engine and reference data versions used to build
        the MIRI/IP trace assets.

    Returns
    -------
    dict
        Engine, source-revision, reference-data, and PSF-library versions.

    Raises
    ------
    EnvironmentError
        If Pandeia reference data or its PSF library is not configured.
    RuntimeError
        If ``require_miri`` is true and the active engine or reference data
        does not match the required MIRI/IP generation version.
    """

    import pandeia.engine
    from pandeia.engine import config

    module_path = os.path.realpath(pandeia.engine.__file__)
    refdata = config.default_refdata()
    psfs = config.default_psfs()
    if not refdata or not psfs:
        raise EnvironmentError('Set pandeia_refdata and PSF_DIR before generation')
    with open(os.path.join(refdata, 'VERSION_DATA')) as handle:
        data_version = handle.readline().strip()
    with open(os.path.join(psfs, 'VERSION_PSF')) as handle:
        psf_version = handle.readline().strip()

    # Editable/source injection over an older environment can leave stale
    # importlib metadata.  The setup.py beside the imported source is the
    # authoritative RC declaration and both values are reported.
    setup_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(module_path))), 'setup.py')
    declared_version = None
    source_revision = None
    if os.path.isfile(setup_path):
        with open(setup_path) as handle:
            match = re.search(r'version=["\']([^"\']+)', handle.read())
        declared_version = match.group(1) if match else None
        try:
            source_revision = subprocess.check_output(
                ['git', '-C', os.path.dirname(setup_path), 'describe',
                 '--tags', '--always', '--dirty'], text=True).strip()
        except (OSError, subprocess.CalledProcessError):
            pass

    provenance = {
        'engine_declared_version': declared_version or pandeia.engine.__version__,
        'engine_metadata_version': pandeia.engine.__version__,
        'engine_source_revision': source_revision,
        'reference_data_version': data_version,
        'psf_version': psf_version,
    }
    print(f'Pandeia module:  {module_path}')
    print(f'Engine declared: {provenance["engine_declared_version"]}')
    print(f'Engine metadata: {provenance["engine_metadata_version"]}')
    print(f'Engine revision: {source_revision or "unavailable"}')
    print(f'Reference data: {os.path.realpath(refdata)} ({data_version})')
    print(f'PSF library:    {os.path.realpath(psfs)} ({psf_version})')

    if require_miri and provenance['engine_declared_version'] != '2026.7':
        raise RuntimeError('MIRI/IP assets require Pandeia engine 2026.7')
    if require_miri and data_version != '2026.7':
        raise RuntimeError('MIRI/IP assets require Pandeia refdata 2026.7')
    return provenance


def validate_miri_ip_calculation():
    """Run a small MIRI/IP calculation to verify Pandeia compatibility.

    Returns
    -------
    provenance : dict
        Version information returned by :func:`pandeia_provenance`.
    report : dict
        Pandeia calculation report for the validation scene.

    Raises
    ------
    RuntimeError
        If Pandeia reports compatibility warnings or its required versions do
        not match the asset-generation environment.
    """

    from pandeia.engine.calc_utils import build_default_calc
    from pandeia.engine.perform_calculation import perform_calculation

    provenance = pandeia_provenance()
    calculation = build_default_calc('jwst', 'miri', 'lrsslitless')
    calculation['configuration']['detector']['subarray'] = 'slitlessprism_ip'
    calculation['strategy']['background_subtraction'] = False
    report = perform_calculation(calculation, webapp=False)
    if report.get('warnings'):
        raise RuntimeError(f'MIRI/IP compatibility warnings: {report["warnings"]}')
    shape = np.asarray(report['2d']['detector']).shape
    wave = np.asarray(report['1d']['wave_pix'])
    print('MIRI configuration: lrsslitless / slitlessprism_ip available')
    print(f'MIRI calculation:   compatible; grid={shape}, wavelengths={wave.size}')
    print(f'Extraction aperture: {report["scalar"]["aperture_size"]} arcsec; '
          f'{report["scalar"]["extraction_area"]} pixels')
    return provenance, report


def _scene(teff, norm_mag):
    """Build a point-source scene normalized in Gaia G.

    Parameters
    ----------
    teff : int or float
        Phoenix effective temperature in kelvin.
    norm_mag : float
        Vega magnitude in the Gaia G band.

    Returns
    -------
    dict
        Pandeia scene-source configuration.
    """
    return {
        'position': {'x_offset': 0., 'y_offset': 0., 'orientation': 0.},
        'shape': {'geometry': 'point'},
        'spectrum': {
            'name': 'Phoenix Spectrum',
            'normalization': {
                'type': 'photsys', 'bandpass': 'gaia,g',
                'norm_flux': norm_mag, 'norm_fluxunit': 'vegamag'},
            'sed': {'sed_type': 'phoenix', 'teff': int(teff), 'log_g': 5.0,
                    'metallicity': 0.0},
        },
    }


def _miri_product(report, source_detector, aperture_size,
                  calibrated_wave=None):
    """Map Pandeia's expanded slitless grid into the 384x68 IP SCI frame.

    ``report['2d']['detector']`` contains a random noise realization and must
    not be used as a deterministic trace template.  ``source_detector`` is
    Pandeia's noiseless source-rate plane, already flipped into its displayed
    P750L orientation.  Register its central IP crop with the wavelength grid
    and retain the expanded rows on both sides so shifted-source PSF wings are
    cropped only after translation.

    Parameters
    ----------
    report : dict
        Pandeia report containing its one-dimensional wavelength grid and
        scalar extraction metadata.
    source_detector : numpy.ndarray
        Expanded, noiseless Pandeia source-rate plane in displayed P750L
        orientation.
    aperture_size : float
        Extraction-aperture size in arcseconds.
    calibrated_wave : array-like, optional
        Wavelength samples to register on the native IP subarray.  By default,
        all projected Pandeia wavelength samples are used.

    Returns
    -------
    trace : numpy.ndarray
        Source-rate image on the native MIRI/IP subarray.
    wavelength : numpy.ndarray
        Wavelength associated with each native detector row.
    valid : numpy.ndarray
        Boolean mask selecting calibrated wavelength rows.
    mask : numpy.ndarray
        Target extraction mask.
    blue_extension, red_extension : numpy.ndarray
        Source-rate rows lying beyond the native detector boundaries.
    registration : dict
        Metadata describing wavelength and detector registration.

    Raises
    ------
    RuntimeError
        If the Pandeia product is invalid or cannot be registered on the
        native MIRI/IP subarray.
    """

    # Validate the deterministic source-rate plane before using it as a
    # reusable trace template.  Tiny negative values from numerical noise are
    # harmless, but appreciably negative source rates indicate a bad product.
    detector = np.asarray(source_detector, dtype=float)
    if detector.ndim != 2 or not np.all(np.isfinite(detector)):
        raise RuntimeError('Invalid Pandeia noiseless source-rate plane')
    negative_tolerance = 1.e-12 * max(float(detector.max()), 1.)
    if detector.min() < -negative_tolerance:
        raise RuntimeError(
            'Pandeia noiseless source-rate plane is negative: '
            f'min={detector.min():.6g}')
    detector = np.maximum(detector, 0.)
    # The Pandeia grid may include extrapolated red-end samples, while
    # calibrated_wave can select a shorter calibrated subset.  Keep both
    # lengths because they determine different parts of the registration.
    projected_wave = np.asarray(report['1d']['wave_pix'], dtype=float)
    wave = (projected_wave if calibrated_wave is None else
            np.asarray(calibrated_wave, dtype=float))
    n_projected = len(projected_wave)
    n_wave = len(wave)
    if n_wave > n_projected:
        raise RuntimeError('Calibrated wavelength grid exceeds Pandeia grid')
    # Align the wavelength-bearing rows with the native IP subarray, then
    # align the narrow Pandeia detector plane with the SIAF reference column.
    pandeia_row0 = (detector.shape[0] - n_projected) // 2 + 1
    output_row0 = (miri_lrs.SHAPE[0] - n_wave) // 2
    detector_row0 = pandeia_row0 - output_row0
    reference_x = 34  # zero-based form of SIAF XSciRef=34.5
    pandeia_x = detector.shape[1] // 2
    output_x0 = reference_x - pandeia_x

    # Split the expanded Pandeia plane into the on-detector template and the
    # two off-detector extensions.  The extensions can later enter the native
    # frame when a neighboring source is translated along the dispersion axis.
    detector_row1 = detector_row0 + miri_lrs.SHAPE[0]
    if detector_row0 < 0 or detector_row1 > detector.shape[0]:
        raise RuntimeError('Pandeia detector scene cannot cover the MIRI/IP subarray')
    trace = np.zeros(miri_lrs.SHAPE, dtype=float)
    trace[:, output_x0:output_x0 + detector.shape[1]] = detector[
        detector_row0:detector_row1]
    blue_extension = np.zeros((detector_row0, miri_lrs.SHAPE[1]), dtype=float)
    blue_extension[:, output_x0:output_x0 + detector.shape[1]] = detector[
        :detector_row0]
    red_extension = np.zeros(
        (detector.shape[0] - detector_row1, miri_lrs.SHAPE[1]), dtype=float)
    red_extension[:, output_x0:output_x0 + detector.shape[1]] = detector[
        detector_row1:]
    # Register wavelengths only on rows covered by the selected grid; NaNs
    # distinguish the remaining detector rows from calibrated spectral data.
    wavelength = np.full(miri_lrs.SHAPE[0], np.nan)
    # Pandeia reports wave_pix from red to blue for this mode, while detector
    # row index increases from the blue end to the red end.  Its extracted 1-D
    # flux therefore matches the detector rows only when wave_pix is reversed.
    wavelength[output_row0:output_row0 + n_wave] = wave[::-1]
    valid = np.isfinite(wavelength)

    # Pandeia SpecApPhot reports 12 pixels for its 1.32 arcsec strip.  Preserve
    # that definition and center it on the SIAF target reference column.
    width = int(round(float(report['scalar']['extraction_area'])))
    if not np.isclose(width, report['scalar']['extraction_area']):
        raise RuntimeError('Non-integral MIRI extraction aperture needs weighted-mask support')
    mask = np.zeros(miri_lrs.SHAPE, dtype=float)
    x0 = reference_x - width // 2 + 1
    mask[valid, x0:x0 + width] = 1.0
    # Store enough registration metadata to audit how the Pandeia product was
    # cropped and to place its off-detector extensions during later rendering.
    registration = {
        'pandeia_wavelength_row0': pandeia_row0,
        'pandeia_projected_wavelength_rows': int(n_projected),
        'pandeia_calibrated_wavelength_rows': int(n_wave),
        'pandeia_detector_row0': detector_row0,
        'pandeia_detector_product': 'noiseless_source_rate',
        'pandeia_expanded_detector_rows': int(detector.shape[0]),
        'ip_wavelength_row0': output_row0,
        'ip_cross_dispersion_col0': output_x0,
        'blue_extension_y0': int(-detector_row0),
        'blue_extension_rows': int(blue_extension.shape[0]),
        'red_extension_y0': int(miri_lrs.SHAPE[0]),
        'red_extension_rows': int(red_extension.shape[0]),
        'extension_method': 'pandeia_noiseless_source_rate',
        'extraction_width_pixels': width,
        'extraction_aperture_arcsec': aperture_size,
        'registration_status': 'experimental_requires_manual_validation',
    }
    return (trace, wavelength, valid, mask, blue_extension, red_extension,
            registration)


def _extended_miri_wave_pix(wave_pix, zero_wavelength, fit_rows=32):
    """Extrapolate the P750L pixel grid through its terminal response zero.

    Pandeia's MIRI/IP dispersion grid currently stops at about 13.86 microns,
    although the P750L response remains nonzero for several more detector
    pixels.  Continue the local wavelength solution by detector pixel so the
    engine can project those wavelengths instead of leaving only the PSF wing
    of its final calibrated sample.

    Parameters
    ----------
    wave_pix : array-like
        Native, monotonically increasing wavelength-per-pixel grid in microns.
    zero_wavelength : float
        Wavelength of the first terminal zero in the P750L response table.
    fit_rows : int, optional
        Number of red-end native pixels used for the local quadratic fit.

    Returns
    -------
    numpy.ndarray
        Native wavelength grid extended exactly to ``zero_wavelength``.

    Raises
    ------
    RuntimeError
        If the input grid is invalid, extrapolation ceases to increase, or the
        response zero cannot be reached within the safety limit.
    """

    wave = np.asarray(wave_pix, dtype=float)
    if (wave.ndim != 1 or len(wave) < fit_rows or
            not np.all(np.isfinite(wave)) or not np.all(np.diff(wave) > 0.)):
        raise RuntimeError('Invalid P750L wavelength/pixel grid')
    if zero_wavelength <= wave[-1]:
        raise RuntimeError('P750L response zero is not redward of pixel grid')

    pixels = np.arange(len(wave), dtype=float)
    coefficients = np.polyfit(
        pixels[-fit_rows:], wave[-fit_rows:], 2)
    extended = wave.tolist()
    for pixel in range(len(wave), len(wave) + 128):
        value = float(np.polyval(coefficients, pixel))
        if value <= extended[-1]:
            raise RuntimeError('P750L wavelength extrapolation is not increasing')
        extended.append(min(value, float(zero_wavelength)))
        if extended[-1] == float(zero_wavelength):
            break
    else:
        raise RuntimeError('Could not reach P750L response zero in 128 pixels')
    return np.asarray(extended)


@contextmanager
def _extended_miri_pandeia_grid(zero_wavelength):
    """Temporarily extend the MIRI/IP pixel grid to its response-table zero.

    The engine limits MIRI/IP to the end of its calibrated dispersion table,
    while the P750L response table remains nonzero slightly farther redward.
    Extend only this mode's wavelength-to-pixel relation and use its local
    pixel spacing as the added dispersion. Pandeia continues to evaluate its
    unmodified response table, SED, remaining throughput terms, quantum yield,
    and PSF projection at the added wavelengths.

    Parameters
    ----------
    zero_wavelength : float
        Wavelength of the terminal P750L response zero, in microns.

    Yields
    ------
    None
        Control while MIRI LRS wavelength-grid methods are temporarily
        patched.  Original methods are restored on context exit.
    """

    from pandeia.engine.jwst import MIRI

    original_wave_pix = MIRI.get_wave_pix
    original_wave_range = MIRI.get_wave_range
    original_dispersion = MIRI.get_dispersion

    def get_dispersion(instance, wave):
        """Evaluate native or extrapolated dispersion at input wavelengths.

        Parameters
        ----------
        instance : pandeia.engine.jwst.MIRI
            Instrument instance being evaluated.
        wave : array-like
            Wavelengths in microns.

        Returns
        -------
        numpy.ndarray
            Dispersion in microns per detector pixel.
        """
        if instance.mode != 'lrsslitless':
            return original_dispersion(instance, wave)
        native = np.asarray(original_wave_pix(instance), dtype=float)
        values = np.asarray(wave, dtype=float)
        red = values > native[-1]
        # Do not ask the native dispersion interpolator to evaluate outside
        # its calibrated wavelength grid. Its value there is replaced below
        # by the local pixel spacing of the extrapolated wavelength solution.
        dispersion = np.empty_like(values, dtype=float)
        if np.any(~red):
            dispersion[~red] = original_dispersion(instance, values[~red])
        if np.any(red):
            extended = _extended_miri_wave_pix(native, zero_wavelength)
            pixel_width = np.gradient(extended)
            dispersion[red] = np.interp(
                values[red], extended, pixel_width)
        return dispersion

    with (patch.object(
              MIRI, 'get_wave_pix',
              lambda instance: _extended_miri_wave_pix(
                  original_wave_pix(instance), zero_wavelength)
              if instance.mode == 'lrsslitless'
              else original_wave_pix(instance)),
          patch.object(
              MIRI, 'get_wave_range',
              lambda instance: dict(
                  original_wave_range(instance),
                  wmax=float(zero_wavelength))
              if instance.mode == 'lrsslitless'
              else original_wave_range(instance)),
          patch.object(MIRI, 'get_dispersion', get_dispersion)):
        yield


def _miri_red_response_limit(calculation):
    """Find the P750L terminal response zero and native grid limit.

    Parameters
    ----------
    calculation : dict
        Pandeia calculation configuration for MIRI LRS slitless mode.

    Returns
    -------
    dict
        Terminal zero wavelength, native red wavelength limit, and response
        filename.

    Raises
    ------
    RuntimeError
        If the response lacks a terminal zero or that zero is not redward of
        the native wavelength grid.
    """

    from pandeia.engine.instrument_factory import InstrumentFactory

    instrument = InstrumentFactory(config=calculation['configuration'])
    key = '{}_{}'.format(
        instrument.instrument['aperture'],
        instrument.instrument['disperser'])
    path = os.path.join(instrument.ref_dir, instrument.paths[key])
    response = fits.getdata(path, 1)
    wave = np.asarray(response['WAVELENGTH'], dtype=float)
    throughput = np.asarray(response['THROUGHPUT'], dtype=float)
    positive = np.flatnonzero(np.isfinite(throughput) & (throughput > 0.))
    if not len(positive) or positive[-1] + 1 >= len(wave):
        raise RuntimeError('P750L response has no terminal zero-throughput sample')
    reference_zero_index = positive[-1] + 1
    if throughput[reference_zero_index] != 0.:
        raise RuntimeError('P750L terminal response sample is not zero')

    native_red_limit = float(instrument.get_wave_range()['wmax'])
    terminal_zero = float(wave[reference_zero_index])
    if terminal_zero <= native_red_limit:
        raise RuntimeError(
            'P750L response zero is not redward of the native pixel grid')
    return {
        'zero': terminal_zero,
        'native_red_limit': native_red_limit,
        'response_file': os.path.basename(path),
    }


def generate_miri_ip_traces(min_teff=2800, max_teff=6000, increment=100,
                            norm_mag=20., outdir=None):
    """Generate ``MIRIM_SLITLESSPRISM_IP`` trace assets with Pandeia.

    Parameters
    ----------
    min_teff, max_teff : int, optional
        Inclusive effective-temperature range in kelvin.
    increment : int, optional
        Temperature spacing in kelvin.
    norm_mag : float, optional
        Vega magnitude used to normalize every Phoenix scene in Gaia G.
    outdir : str or path-like, optional
        Output directory.  By default, assets are written to
        :func:`miri_lrs.trace_directory`.

    Returns
    -------
    None

    Raises
    ------
    EnvironmentError
        If Pandeia reference data or the ExoCTK data root is unavailable.
    RuntimeError
        If the Pandeia environment is incompatible, a calculation emits
        warnings, or a detector product cannot be registered safely.
    """

    from pandeia.engine.calc_utils import build_default_calc
    from pandeia.engine.perform_calculation import perform_calculation

    provenance = pandeia_provenance()
    output = outdir or miri_lrs.trace_directory()
    os.makedirs(output, exist_ok=True)

    base = build_default_calc('jwst', 'miri', 'lrsslitless')
    base['configuration']['detector']['subarray'] = 'slitlessprism_ip'
    base['strategy']['background_subtraction'] = False
    aperture_size = base['strategy']['aperture_size']
    red_response = _miri_red_response_limit(base)
    red_zero = red_response['zero']
    native_red_limit = red_response['native_red_limit']
    with _extended_miri_pandeia_grid(red_zero):
        for teff in np.arange(min_teff, max_teff + increment, increment):
            calculation = copy.deepcopy(base)
            calculation['scene'] = [_scene(teff, norm_mag)]
            report_object = perform_calculation(
                calculation, dict_report=False, webapp=False)
            report = report_object.as_dict()
            if report.get('warnings'):
                raise RuntimeError(
                    f'Pandeia warnings for MIRI/IP: {report["warnings"]}')
            # Match Pandeia's public P750L display orientation without using
            # its stochastic detector report, which adds a random noise
            # realization. Keep the native calibrated target wavelengths;
            # the added red rows are solely an off-detector source extension.
            source_detector = np.asarray(
                report_object.signal.rate, dtype=float)[::-1]
            projected_wave = np.asarray(
                report['1d']['wave_pix'], dtype=float)
            calibrated_wave = projected_wave[
                projected_wave <= native_red_limit]
            (trace, wavelength, valid, mask, blue_extension, red_extension,
             registration) = _miri_product(
                report, source_detector, aperture_size,
                calibrated_wave=calibrated_wave)
            metadata = dict(provenance, **registration)
            metadata.update({
                'source_teff': int(teff),
                'normalization_magnitude': float(norm_mag),
                'normalization_bandpass': 'gaia,g',
                'shape': list(miri_lrs.SHAPE), 'dispersion_axis': 0,
                'cross_dispersion_axis': 1,
                'wavelength_orientation': 'increasing_toward_negative_y_sci',
                'short_wavelength_foldover': 'none_in_pandeia_wave_pix',
                'p750l_response_file': red_response['response_file'],
                'p750l_native_red_limit_micron': float(native_red_limit),
                'p750l_terminal_zero_wavelength_micron': red_zero,
                'red_extension_method':
                    'pandeia_extended_pixel_grid_reference_response',
            })
            primary = fits.PrimaryHDU()
            primary.header['APERTURE'] = miri_lrs.APERTURE
            primary.header['REFYPIX'] = 290
            primary.header['REFXPIX'] = 34
            primary.header['METAJSON'] = json.dumps(
                metadata, separators=(',', ':'))
            hdul = fits.HDUList([
                primary,
                fits.ImageHDU(trace.astype('float32'), name='TRACE'),
                fits.ImageHDU(
                    wavelength.astype('float64'), name='WAVELENGTH'),
                fits.ImageHDU(
                    mask.astype('float32'), name='EXTRACTION_MASK'),
                fits.ImageHDU(
                    valid.astype('uint8'), name='VALID_WAVELENGTH'),
                fits.ImageHDU(blue_extension.astype('float32'),
                              name='BLUE_EXTENSION'),
                fits.ImageHDU(red_extension.astype('float32'),
                              name='RED_EXTENSION'),
            ])
            path = os.path.join(
                output, f'{miri_lrs.APERTURE}_{int(teff)}.fits')
            hdul.writeto(path, overwrite=True)
            print(f'Saved {path}')


def main():
    """Run the MIRI/IP trace-generation command-line interface.

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validate-only', action='store_true')
    parser.add_argument('--min-teff', type=int, default=2800)
    parser.add_argument('--max-teff', type=int, default=6000)
    parser.add_argument('--increment', type=int, default=100)
    parser.add_argument('--norm-mag', type=float, default=20.)
    parser.add_argument('--outdir')
    args = parser.parse_args()
    if args.validate_only:
        validate_miri_ip_calculation()
    else:
        generate_miri_ip_traces(args.min_teff, args.max_teff, args.increment,
                                args.norm_mag, args.outdir)


if __name__ == '__main__':
    main()

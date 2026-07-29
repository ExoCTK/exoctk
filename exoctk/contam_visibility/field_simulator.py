#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module to calculate the contamination and visibility of a target on a JWST detector
"""

from copy import copy
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache, partial
import glob
import logging
import os
from pathlib import Path
import pickle
import re
import sys
import time
import logging
import requests
from urllib.parse import quote_plus

import astropy.coordinates as crd
from astropy.io import fits
import astropy.units as u
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroquery.xmatch import XMatch
from bokeh.plotting import figure, show
from bokeh.embed import json_item
from bokeh.layouts import gridplot, column
from bokeh.models import Range1d, LinearColorMapper, LogColorMapper, Label, ColorBar, ColumnDataSource, HoverTool, Slider, CustomJS, VArea, CrosshairTool, TapTool, OpenURL, Span, Legend
from bokeh.palettes import PuBu, Spectral6
from bokeh.transform import linear_cmap
from scipy.ndimage import rotate
import h5py
import numpy as np
import pysiaf
import regions

from ..utils import (
    add_array_at_position, add_scaled_array_inplace, check_for_data,
    get_canonical_name, get_env_variables, get_target_data, replace_NaNs)
from ..pkgdata import resource_filename
from ..modelgrid import ACES
from .new_vis_plot import build_visibility_plot, get_exoplanet_positions
from . import contamination_figure as cf
from . import miri_lrs
from .gaia_tap import GaiaFailoverTAP
from .precompute import save_exoplanet_data

log_file = 'contam_tool.log'
logging.basicConfig(
    filename=log_file,
    filemode='w',   # <-- THIS forces overwrite on every run
    level=logging.INFO,
    format='%(asctime)s %(message)s',
    force=True
)

_last_time = datetime.now()


@dataclass(frozen=True)
class DHSContaminationResult(Sequence):
    """Compact fractional-contamination result for NIRCam DHS.

    ``order_fractions`` contains one float64 array per DHS order. Each array
    has axes ``(V3 position angle, detector wavelength column)`` and shape
    ``(360, n_columns)``. ``position_angles`` maps the first axis to integer
    V3 position angles. At unobservable angles, wavelength columns covered by
    the extraction mask contain zero contamination and uncovered columns
    remain NaN, matching reduction of an all-zero detector plane.

    The sequence interface preserves callers that iterate over the former
    list of order arrays while providing an explicit, inspectable schema.
    Target-order detector images remain the first value returned by
    :func:`field_simulation`; this object represents contamination only.
    """

    order_fractions: tuple
    position_angles: np.ndarray

    def __post_init__(self):
        order_fractions = tuple(
            np.asarray(order, dtype=float) for order in self.order_fractions)
        position_angles = np.asarray(self.position_angles, dtype=int)
        if position_angles.ndim != 1:
            raise ValueError("position_angles must be one-dimensional")
        for order in order_fractions:
            if order.ndim != 2 or order.shape[0] != len(position_angles):
                raise ValueError(
                    "Each DHS order must have axes "
                    "(position angle, wavelength column)")
            order.setflags(write=False)
        position_angles.setflags(write=False)
        object.__setattr__(self, "order_fractions", order_fractions)
        object.__setattr__(self, "position_angles", position_angles)

    def __getitem__(self, index):
        return self.order_fractions[index]

    def __len__(self):
        return len(self.order_fractions)


def _cached_contam_result_available(group, compact_dhs=False):
    """Return whether an HDF5 target group has a complete compatible result."""
    if compact_dhs:
        return (
            "goodPA_list" in group.attrs
            and "target_trace" in group
            and "dhs_order_fractions" in group
            and "dhs_position_angles" in group)
    return "goodPA_list" in group.attrs


def log_checkpoint(message):
    global _last_time
    now = datetime.now()
    elapsed = (now - _last_time).total_seconds()
    logging.info(f'{message} (Elapsed: {elapsed:.2f} seconds)')
    _last_time = now

def parse_log():
    timestamps = []
    messages = []
    with open(log_file, 'r') as f:
        for line in f:
            parts = line.strip().split(' ', 1)
            if len(parts) == 2:
                timestamps.append(parts[0])
                messages.append(parts[1])

    for ts, msg in zip(timestamps, messages):
        print(f'{ts}: {msg}')


APERTURES = {'NIS_SOSSFULL': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.066, 'rad': 2.5, 'lam': [0.8, 2.8],
                              'c0x0': 905, 'c0y0': 1467, 'c1x0': -0.013, 'c1y0': -0.1, 'c1y1': 0.12, 'c1x1': -0.03, 'c2y1': -0.011,
                              'subarr_x': [0, 2048, 2048, 0], 'subarr_y':[0, 0, 2048, 2048], 'trim': [127, 126, 252, 1],
                              'lft': 700, 'rgt': 3022, 'top': 2050, 'bot': 1400, 'blue_ext': -150, 'red_ext': 200,
                              'xord0to1': -2886, 'yord0to1': 68, 'empirical_scale': [1, 1.5, 1.5, 1.5],
                              'tracex_offset': 0, 'tracey_offset': 0,
                              'cutoffs': [2048, 1820, 1130], 'trace_names': ['Order 1', 'Order 2', 'Order 3'],
                              'coeffs': [[1.68975801e-11, -4.60822060e-08, 4.94623886e-05, -5.93935390e-02, 8.67263818e+01],
                                         [3.95721278e-11, -7.40683643e-08, 6.88340922e-05, -3.68009540e-02, 1.06704335e+02],
                                         [1.06699517e-11, 3.36931077e-08, 1.45570667e-05, 1.69277607e-02, 1.45254339e+02]]},
             'NIS_SUBSTRIP96': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.066, 'rad': 2.5, 'lam': [0.8, 2.8],
                                'c0x0': 905, 'c0y0': 1467, 'c1x0': -0.013, 'c1y0': -0.1, 'c1y1': 0.12, 'c1x1': -0.03, 'c2y1': -0.011,
                                'subarr_x': [0, 2048, 2048, 0], 'subarr_y':[1792, 1792, 1888, 1888], 'trim': [47, 46, 0, 1],
                                'lft': 700, 'rgt': 3022, 'top': 2050, 'bot': 1400, 'blue_ext': -150, 'red_ext': 200,
                                'xord0to1': -2886, 'yord0to1': 68, 'empirical_scale': [0.001, 1, 1, 1],
                                'tracex_offset': 0, 'tracey_offset': 0,
                                'cutoffs': [2048, 1820, 1130], 'trace_names': ['Order 1', 'Order 2', 'Order 3'],
                                'coeffs': [[1.68975801e-11, -4.60822060e-08, 4.94623886e-05, -5.93935390e-02, 8.67263818e+01],
                                           [3.95721278e-11, -7.40683643e-08, 6.88340922e-05, -3.68009540e-02, 1.06704335e+02],
                                           [1.06699517e-11, 3.36931077e-08, 1.45570667e-05, 1.69277607e-02, 1.45254339e+02]]},
             'NIS_SUBSTRIP256': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.066, 'rad': 2.5, 'lam': [0.8, 2.8],
                                 'c0x0': 905, 'c0y0': 1467, 'c1x0': -0.013, 'c1y0': -0.1, 'c1y1': 0.12, 'c1x1': -0.03, 'c2y1': -0.011,
                                 'subarr_x': [0, 2048, 2048, 0], 'subarr_y':[1792, 1792, 2048, 2048], 'trim': [127, 126, 0, 1],
                                 'lft': 700, 'rgt': 3022, 'top': 2050, 'bot': 1400, 'blue_ext': -150, 'red_ext': 200,
                                 'xord0to1': -2886, 'yord0to1': 68, 'empirical_scale': [0.001, 1, 1, 1],
                                 'tracex_offset': 0, 'tracey_offset': 0,
                                 'cutoffs': [2048, 1820, 1130], 'trace_names': ['Order 1', 'Order 2', 'Order 3'],
                                 'coeffs': [[1.68975801e-11, -4.60822060e-08, 4.94623886e-05, -5.93935390e-02, 8.67263818e+01],
                                            [3.95721278e-11, -7.40683643e-08, 6.88340922e-05, -3.68009540e-02, 1.06704335e+02],
                                            [1.06699517e-11, 3.36931077e-08, 1.45570667e-05, 1.69277607e-02, 1.45254339e+02]]},
             'NRCA5_41STRIPE1_DHS_F322W2': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.031, 'rad': 2.5, 'lam': [1.07, 2.01],
                                            'subarr_x': [0, 4257, 4257, 0], 'subarr_y': [1512, 1512, 2744, 2744], 'trim': [0, 1, 0, 1],
                                            'c0x0': 1800, 'c0y0': 2116, 'c1x0': 0, 'c1y0': 0, 'c1y1': 0, 'c1x1': 0, 'c2y1': 0,
                                            'lft': 0, 'rgt': 4300, 'top': 4000, 'bot': 0, 'blue_ext': 0, 'red_ext': 0,
                                            'xord0to1': -2448, 'yord0to1': -2116, 'empirical_scale': [1] * 11,
                                            'tracex_offset': 0, 'tracey_offset': 0,
                                            'cutoffs': [3324]*10, 'trace_names': ['DHS5', 'DHS4', 'DHS3', 'DHS2', 'DHS1', 'DHS6', 'DHS7', 'DHS8', 'DHS9', 'DHS10'],
                                            'coeffs': [[ 2.81442298e-06,  1.48386268e-04,  1.14039356e+03],
                                                       [ 2.74877495e-06, -5.87120866e-04,  1.02983316e+03],
                                                       [ 2.91085676e-06, -3.23811050e-03,  0.91454889e+03],
                                                       [ 2.79902255e-06, -4.32186046e-03,  0.78638977e+03],
                                                       [ 2.59964469e-06, -4.05419001e-03,  0.66387590e+03],
                                                       [ 2.40026684e-06, -3.78651956e-03,  0.54136203e+03],
                                                       [ 2.20088898e-06, -3.51884910e-03,  0.41884815e+03],
                                                       [ 1.99460079e-06, -3.82768919e-03,  0.30647599e+03],
                                                       [ 2.13733247e-06, -6.65892407e-03,  0.18528819e+03],
                                                       [ 2.22766890e-06, -8.81816227e-03,  0.05636400e+03]]},
             'NRCA5_41STRIPE1_DHS_F444W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.031, 'rad': 2.5, 'lam': [1.07, 2.01],
                                           'subarr_x': [0, 4257, 4257, 0], 'subarr_y': [1512, 1512, 2744, 2744], 'trim': [0, 1, 0, 1],
                                           'c0x0': 900, 'c0y0': 2116, 'c1x0': 0, 'c1y0': 0, 'c1y1': 0, 'c1x1': 0, 'c2y1': 0,
                                           'lft': 0, 'rgt': 4300, 'top': 4000, 'bot': 0, 'blue_ext': 0, 'red_ext': 0,
                                           'xord0to1': -2000, 'yord0to1': -2116, 'empirical_scale': [1] * 11,
                                           'tracex_offset': 0, 'tracey_offset': 0,
                                           'cutoffs': [3324]*10, 'trace_names': ['DHS5', 'DHS4', 'DHS3', 'DHS2', 'DHS1', 'DHS6', 'DHS7', 'DHS8', 'DHS9', 'DHS10'],
                                           'coeffs': [[3.52279496e-06, -1.22488543e-03,  1.15263106e+03],
                                                     [ 2.93109211e-06, -1.22659432e-03,  1.04145870e+03],
                                                     [ 2.93877311e-06, -2.71105945e-03,  0.92296188e+03],
                                                     [ 2.59336772e-06, -2.78664360e-03,  0.79192880e+03],
                                                     [ 2.76081057e-06, -3.64650141e-03,  0.66955923e+03],
                                                     [ 2.92825343e-06, -4.50635923e-03,  0.54718966e+03],
                                                     [ 3.09569628e-06, -5.36621705e-03,  0.42482010e+03],
                                                     [ 2.12912185e-06, -3.98286955e-03,  0.31058917e+03],
                                                     [ 2.12853718e-06, -6.00260391e-03,  0.18598335e+03],
                                                     [ 2.47832831e-06, -8.20469374e-03,  0.05372619e+03]]},
             'NRCA5_GRISM256_F277W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [2.395, 3.179], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F322W2': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [2.413, 4.083], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F356W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [3.100, 4.041], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F444W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [3.835, 5.084], 'trim': [0, 1, 1250, 1]},
             'MIRIM_SLITLESSPRISM': {'inst': 'MIRI', 'full': 'MIRIM_FULL', 'scale': 0.11, 'rad': 2.0, 'lam': [5, 12], 'trim': [6, 5, 0, 1]},
             'MIRIM_SLITLESSPRISM_IP': {
                 'inst': 'MIRI', 'full': 'MIRIM_FULL', 'scale': 0.11,
                 'source_search_radius_arcsec':
                     miri_lrs.SOURCE_SEARCH_RADIUS_ARCSEC,
                 'shape': miri_lrs.SHAPE,
                 'dispersion_axis': miri_lrs.DISPERSION_AXIS,
                 'cross_dispersion_axis': miri_lrs.CROSS_DISPERSION_AXIS,
                 'trace_path': 'exoctk_contam/traces/MIRIM_SLITLESSPRISM_IP',
                 'trace_reference_pixel': (290, 34),
                 'science_bounds': ((0, 384), (0, 68)),
                 'wavelength_source': 'trace_asset:WAVELENGTH',
                 'extraction_mask_source': 'trace_asset:EXTRACTION_MASK',
             }}

WEB_CONTAMINATION_APERTURES = frozenset({
    'NIS_SUBSTRIP96',
    'NIS_SUBSTRIP256',
    'NRCA5_41STRIPE1_DHS_F322W2',
    'NRCA5_41STRIPE1_DHS_F444W',
    'MIRIM_SLITLESSPRISM_IP',
})

# Require a majority probability before rendering a source as extended.
GAIA_DSC_GALAXY_THRESHOLD = 0.5


def contamination_supported(aperture):
    """Return whether the contamination web interface supports an aperture."""

    return aperture in WEB_CONTAMINATION_APERTURES


def source_query_width(aperture):
    """Return the Gaia query width needed to cover an aperture's sources.

    Parameters
    ----------
    aperture : str
        Supported contamination aperture name.

    Returns
    -------
    astropy.units.Quantity
        Side length of the square Gaia query region.  Its circumscribed circle
        covers the mode's source-search radius.
    """

    if aperture == miri_lrs.APERTURE:
        radius = APERTURES[aperture]['source_search_radius_arcsec']
        # GaiaFailoverTAP circumscribes the requested square, so a square with
        # side sqrt(2)*radius produces exactly the desired circular radius.
        return np.sqrt(2.) * radius * u.arcsec
    return 7.5 * u.arcmin


def calculation_v3pas(aperture, observable_pas):
    """Choose the V3 position angles at which contamination is calculated.

    Parameters
    ----------
    aperture : str
        Supported contamination aperture name.
    observable_pas : array-like
        V3 position angles observable in the requested epoch.

    Returns
    -------
    numpy.ndarray
        All integer angles from 0 through 359 for MIRI LRS, or the supplied
        observable angles for other modes.
    """

    if aperture == miri_lrs.APERTURE:
        return np.arange(360)
    return np.asarray(observable_pas)


def unobservable_v3pas(pa_results):
    """Derive year-specific inaccessible V3 position angles.

    Parameters
    ----------
    pa_results : sequence of int or miri_lrs.PositionAngleResults
        Successfully calculated angles, optionally with explicit inaccessible
        angle metadata.

    Returns
    -------
    list of int
        V3 position angles to shade as not observable.
    """

    if isinstance(pa_results, miri_lrs.PositionAngleResults):
        return list(pa_results.inaccessible)
    return [pa for pa in np.arange(360) if pa not in pa_results]


# Gaia color-Teff relation
GAIA_TEFFS = np.asarray(np.genfromtxt(resource_filename('exoctk', 'data/contam_visibility/predicted_gaia_colour.txt'), unpack=True))

# Gaia TAP instance
GAIA_TAP = GaiaFailoverTAP()


@lru_cache(maxsize=1)
def _get_aces_grid():
    """Load the ACES grid only when a trace actually needs scaling."""

    return ACES()


def NIRCam_detector_gap():
    """
    Binary mask for 161px wide NIRCam detector gap
    """
    det = np.zeros((4257, 4257))
    det[2048:2209, :] = 1
    det[:, 2048:2209] = 1

    return det[1512:2744, :]


def get_trace_mask(aperture, radius=20, plot=False):
    """
    Construct a trace mask for the given aperture

    Parameters
    ----------
    radius: int
        The radius in pixels of the extraction aperture

    Returns
    -------
    np.ndarray
        The trace masks
    """
    ydim = APERTURES[aperture]['subarr_y'][2] - APERTURES[aperture]['subarr_y'][1]
    xdim = APERTURES[aperture]['subarr_x'][1] - APERTURES[aperture]['subarr_x'][0]
    coeffs = APERTURES[aperture]['coeffs']

    x = np.arange(xdim)
    traces = np.array([np.polyval(coeff, x) for coeff in coeffs]).astype(int)

    y_grid = np.arange(ydim)[None, :, None]
    trace_grid = traces[:, None, :]

    masks = np.abs(y_grid - trace_grid) < radius

    if plot:

        from bokeh.models import LinearColorMapper
        from bokeh.palettes import Greys256
        color_mapper = LinearColorMapper(palette=Greys256, low=0, high=1)
        f = figure()
        for trace, mask in zip(traces, masks):
            f.image([mask.astype(np.uint8)], y=APERTURES[aperture]['subarr_y'][1], dh=ydim, x=0, dw=xdim, color_mapper=color_mapper, alpha=0.1)
            f.line(x, trace, color='red', line_width=2)
        show(f)

    return masks


def _target_source_index(stars, target_coordinate, input_epoch=2000):
    """Identify the Gaia source matching an input-epoch target coordinate.

    Gaia DR3 coordinates are reported at each source's ``ref_epoch`` (normally
    2016), while target coordinates supplied to ExoCTK are conventionally
    J2000.  Matching the unpropagated catalog positions can select a nearby
    field source instead of a high-proper-motion target.
    """

    if not len(stars):
        raise ValueError("Gaia returned no sources near the target")

    match_ra = np.asarray(stars['ra'], dtype=float).copy()
    match_dec = np.asarray(stars['dec'], dtype=float).copy()
    for index, row in enumerate(stars):
        motion = (row['pmra'], row['pmdec'], row['ref_epoch'])
        if any(np.ma.is_masked(value) or not np.isfinite(value)
               for value in motion):
            continue
        match_ra[index], match_dec[index] = calculate_current_coordinates(
            row['ra'], row['dec'], row['pmra'], row['pmdec'],
            row['ref_epoch'], target_date=input_epoch)

    catalog_at_input_epoch = crd.SkyCoord(
        ra=match_ra * u.deg, dec=match_dec * u.deg)
    return int(np.argmin(target_coordinate.separation(catalog_at_input_epoch)))


def observable_v3pa_ranges(position_table, max_gap=7):
    """Return integer V3PAs spanning GTVT visibility windows.

    GTVT's instrument-specific columns are aperture position angles: they
    already include the aperture's ``V3IdlYAngle``. The contamination renderer
    instead accepts V3 position angles and applies the SIAF conversion itself,
    so the V3PA columns must be used here to avoid applying that offset twice.
    """

    columns = [
        'V3PA_min_pa_angle',
        'V3PA_nominal_angle',
        'V3PA_max_pa_angle',
    ]
    missing = [name for name in columns if name not in position_table]
    if missing:
        raise KeyError(f'Missing GTVT V3PA columns: {missing}')

    values = []
    for name in columns:
        column_values = np.asarray(position_table[name], dtype=float)
        values.extend(column_values[np.isfinite(column_values)])
    if not values:
        raise ValueError('GTVT returned no observable V3 position angles')

    sampled = np.sort(np.unique(np.asarray(values).astype(int) % 360))
    groups = [[sampled[0]]]
    for value in sampled[1:]:
        if value - groups[-1][-1] <= max_gap:
            groups[-1].append(value)
        else:
            groups.append([value])

    bounds = [(int(group[0]), int(group[-1])) for group in groups]
    position_angles = np.concatenate([
        np.arange(lower, upper + 1, dtype=int)
        for lower, upper in bounds
    ])
    return position_angles, bounds, sampled


def aperture_pa_from_v3pa(v3pa, aperture):
    """Convert a V3 position angle to an aperture position angle."""

    return (v3pa + aperture.V3IdlYAngle) % 360


def _normalize_sdss_type(value):
    """Map an explicit SDSS class to a contamination rendering class."""

    if np.ma.is_masked(value):
        return None

    value = str(value).strip().upper()
    if value in ('STAR', 'QSO'):
        return 'STAR'
    if value == 'GALAXY':
        return 'GALAXY'
    return None


def _sdss_type_lookup(xmatch_result):
    """Return unanimous usable SDSS classes keyed by Gaia source ID."""

    classifications = {}
    for source_id, value in zip(
            xmatch_result['source_id'], xmatch_result['spCl']):
        normalized = _normalize_sdss_type(value)
        if normalized is not None:
            classifications.setdefault(source_id, set()).add(normalized)

    return {
        source_id: next(iter(values))
        for source_id, values in classifications.items() if len(values) == 1
    }


def _finite_positive_scalar(value):
    """Return a finite positive float, or ``None`` for unusable values."""

    if np.ma.is_masked(value):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) and value > 0 else None


def _filter_valid_flux_sources(stars, column, label):
    """Retain sources with usable fluxes, requiring a usable target flux."""

    if column not in stars.colnames:
        raise ValueError(f'Source table is missing the {label} column.')

    valid_flux = np.array([
        _finite_positive_scalar(value) is not None for value in stars[column]
    ])
    if not len(valid_flux) or not valid_flux[0]:
        raise ValueError(f'The identified target does not have a valid {label}.')

    rejected = np.count_nonzero(~valid_flux)
    if rejected:
        logging.info('Ignoring %d sources without valid %s.', rejected, label)
    return stars[valid_flux]


def find_sources(ra=None, dec=None, target=None, width=5*u.arcmin, target_date=Time.now(), pm_corr=True, plot=False):
    """
    Find all the stars in the vicinity and estimate temperatures

    Parameters
    ----------
    ra : float
        The RA of the target in decimal degrees
    dec : float
        The Dec of the target in decimal degrees
    target: str
        The name of the target to resolve in ExoMAST
    width: astropy.units.quantity
        The width of the square search box
    target_date: Time, int, str
        The target epoch year of the observation, e.g. '2025'
    pm_corr: bool
        Correct source coordinates based on their proper motion

    Returns
    -------
    astropy.table.Table
        The table of stars
    """
    # Resolve target in ExoMAST if possible
    if target is not None:
        targname = get_canonical_name(target)
        data, _ = get_target_data(targname)
        ra = data.get('RA')
        dec = data.get('DEC')
        logging.info(f"Resolved '{targname}' (RA={ra}, Dec={dec}) from '{target}' in ExoMAST.")

    # Converting to degrees and query for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=u.deg if isinstance(ra, float) and isinstance(dec, float) else (u.hour, u.deg))

    # Search Gaia for stars
    logging.info('Searching Gaia DR3 to find all stars within {} of RA={}, Dec={}...'.format(width, ra, dec))

    # Query Gaia from several potential endpoints
    stars = GAIA_TAP.query_region(targetcrd, width=width, height=width)

    # Preserve the intended target as row zero throughout flux normalization,
    # proper-motion correction, and detector rendering. Gaia query order alone
    # is unsafe for targets that have moved since the input-coordinate epoch.
    target_index = _target_source_index(stars, targetcrd)
    order = np.concatenate(([target_index], np.delete(
        np.arange(len(stars)), target_index)))
    stars = stars[order]

    # Sources without measured Gaia G fluxes cannot be normalized. In
    # particular, masked values must not reach the detector arithmetic, where
    # their underlying data could be treated as a valid flux scale.
    stars = _filter_valid_flux_sources(
        stars, 'phot_g_mean_flux', 'Gaia G-band flux')

    try:
        # Perform XMatch between Gaia and SDSS DR16
        xmatch_result = XMatch.query(cat1=stars, cat2='vizier:V/154/sdss16', max_distance=2 * u.arcsec, colRA1='ra', colDec1='dec', colRA2='RA_ICRS', colDec2='DE_ICRS')

        # XMatch can return multiple counterparts per Gaia source. Associate
        # unanimous explicit classes by ID without reordering or expanding.
        sdss_types = _sdss_type_lookup(xmatch_result)
        stars['type'] = [
            sdss_types.get(source_id, '') for source_id in stars['source_id']]

    except Exception as e:
        logging.info(f"Could not perform SDSS crossmatch: {e}")

    # Or infer galaxy from parallax
    stars['type'] = [classify_source(row) for row in stars]

    # Derived from K. Volk
    stars['Teff'] = [GAIA_TEFFS[0][(np.abs(GAIA_TEFFS[1] - row['bp_rp'])).argmin()] for row in stars]

    # # Calculate relative flux in Jband (Hold off on this)
    # GtoJ_coeffs = [-1.22568340e-11, 2.41639448e-07, -1.70092031e-03, 3.15459542e+00] # Measured from synthetic colors
    # log10_gj = np.polyval(GtoJ_coeffs, stars['Teff'])
    # gj_factor_interp = 10 ** log10_gj
    # stars['estimated_j_flux'] = stars['phot_g_mean_flux'] * gj_factor_interp
    # stars['fluxscale'] = stars['estimated_j_flux'] / stars['estimated_j_flux'][0]

    # Calculate relative flux in Gband
    stars['fluxscale'] = stars['phot_g_mean_flux'] / stars['phot_g_mean_flux'][0]

    # Star names
    stars['name'] = [str(i) for i in stars['source_id']]

    # Catalog name
    cat = 'I/355/gaiadr3'

    # Add URL (before PM correcting coordinates)
    search_radius = 1
    urls = ['https://vizier.u-strasbg.fr/viz-bin/VizieR-5?-ref=VIZ62fa613b20f3fc&-out.add=.&-source={}&-c={}&eq=ICRS&rs={}&-out.orig=o'.format(cat, quote_plus(f"{ra_deg} {dec_deg}"), search_radius) for ra_deg, dec_deg in zip(stars['ra'], stars['dec'])]
    stars.add_column(urls, name='url')

    # Cope coordinates to new column
    stars.add_column(stars['ra'], name='obs_ra')
    stars.add_column(stars['dec'], name='obs_dec')

    # Set target data
    if target_date is None:
        target_date = Time.now()

    # Update RA and Dec using proper motion data
    if pm_corr:
        for row in stars:
            new_ra, new_dec = calculate_current_coordinates(row['ra'], row['dec'], row['pmra'], row['pmdec'], row['ref_epoch'], target_date=target_date)

            if not hasattr(new_ra, 'mask'):
                row['ra'] = new_ra
            if not hasattr(new_dec, 'mask'):
                row['dec'] = new_dec

    # Find spherical distance from the identified target to each source.
    coordinates = crd.SkyCoord(
        ra=np.asarray(stars['ra'], dtype=float) * u.deg,
        dec=np.asarray(stars['dec'], dtype=float) * u.deg)
    distances = coordinates[0].separation(coordinates).to_value(u.arcsec)
    stars.add_column(distances, name='distance')
    stars.sort('distance')

    # Add detector location to the table
    stars.add_columns(np.zeros((10, len(stars))), names=['xtel', 'ytel', 'xdet', 'ydet', 'xsci', 'ysci', 'xord0', 'yord0', 'xord1', 'yord1'])

    if plot:
        ra0 = float(targetcrd.ra.deg)
        dec0 = float(targetcrd.dec.deg)

        ra_vals = np.array(stars['ra'])
        dec_vals = np.array(stars['dec'])

        # Convert to relative arcsec offsets
        cos_dec = np.cos(np.deg2rad(dec0))
        dra = (ra_vals - ra0 + 180) % 360 - 180

        x_arcsec = dra * cos_dec * 3600.0
        y_arcsec = (dec_vals - dec0) * 3600.0

        # Marker size scaled by flux (robust scaling)
        flux = np.array(stars['fluxscale'])
        size = 5 + 15 * (flux / np.max(flux))**0.5

        # Make the plot
        source = ColumnDataSource(data=dict(x=x_arcsec, y=y_arcsec, name=stars['name'], size=size, flux=flux))
        p = figure(title="Source Field", x_axis_label="?RA (arcsec)", y_axis_label="?Dec (arcsec)", width=500,
                   height=500, match_aspect=True, tools="pan,wheel_zoom,box_zoom,reset,hover")
        p.scatter('x', 'y', size='size', source=source, alpha=0.6)
        p.scatter(0, 0, size=15, marker='cross')

        show(p)

    return stars


def calculate_current_coordinates(ra, dec, pm_ra, pm_dec, epoch, target_date=None):
    """
    Get the proper motion corrected coordinates of a source

    Parameters
    ----------
    ra: float
        The RA of the source
    dec: float
        The Dec of the source
    pm_ra: float
        The RA proper motion
    pm_dec: float
        The Dec proper motion
    epoch: float
        The epoch of the observation
    target_date: float
        The target epoch

    Returns
    -------
    new_ra, new_dec
        The corrected RA and Dec for the source
    """
    motion = (pm_ra, pm_dec, epoch)
    if any(np.ma.is_masked(value) or not np.isfinite(value)
           for value in motion):
        return np.ma.masked, np.ma.masked

    # Set target data
    if target_date is None:
        target_date = Time.now()

    # Convert observation year to Time object
    if isinstance(target_date, (int, float, str)):
        target_date = Time('{}-01-01'.format(int(target_date)))

    # Convert observation year to Time object
    date_obs = Time('{}-01-01'.format(int(epoch)))

    # Calculate time difference in years
    dt_years = (target_date - date_obs).to(u.year).value

    # Convert proper motion from mas/yr to degrees/yr
    pm_RA_deg_per_year = pm_ra / (3600 * 1000)  # 1 degree = 3600 * 1000 mas
    pm_Dec_deg_per_year = pm_dec / (3600 * 1000)

    # Calculate new coordinates
    new_ra = ra + (pm_RA_deg_per_year * dt_years / np.cos(np.deg2rad(dec)))
    new_dec = dec + (pm_Dec_deg_per_year * dt_years)

    return new_ra, new_dec


def add_source(startable, name, ra, dec, teff=None, fluxscale=None, delta_mag=None, dist=None, pa=None, type='STAR'):
    """
    Add a star to the star table

    Parameters
    ----------
    startable: astropy.table.Table
        The table of stars to add to
    name: str
        An identifier for the star
    ra: float
        The RA in decimal degrees
    dec: float
        The Dec in decimal degrees
    teff: float
        The effective temperature of the star
    fluxscale: float
        The star's flux relative to the target flux
    delta_mag: float
        The star's magnitude relative to the target magnitude
    dist: float
        The distance of the new star from the given RA and Dec in arcseconds
    pa: float
        The position angle of the new star relative to the given RA and Dec in degrees
    type: str
        The source type, ['STAR', 'GALAXY']

    Returns
    -------
    astropy.table.Table
        The updated table of stars
    """
    # Default
    fluxscale = fluxscale or 1

    # Convert mag to flux if necessary
    if delta_mag is not None:
        fluxscale = 10.0 ** (-0.4 * delta_mag)

    # Apply offset and position angle
    if dist is not None and pa is not None:
        coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
        newcoord = coord.directional_offset_by(pa * u.deg, dist * u.arcsec)
        ra = newcoord.ra.degree
        dec = newcoord.dec.degree

    # Add the row to the table
    startable.add_row({'name': name, 'designation': name, 'ra': ra, 'dec': dec, 'obs_ra': ra, 'obs_dec': dec, 'Teff': teff, 'fluxscale': fluxscale, 'type': type, 'distance': dist})
    startable.sort('distance')

    return startable


def _finite_float(row, column_name):
    """Return a finite scalar table value, or NaN when it is unusable."""

    if (column_name not in row.colnames
            or np.ma.is_masked(row[column_name])):
        return np.nan

    try:
        value = float(row[column_name])
    except (TypeError, ValueError):
        return np.nan

    return value if np.isfinite(value) else np.nan


def classify_source(row):
    """
    Classify a Gaia DR3 source as STAR or GALAXY for contamination rendering.

    This binary rendering choice is not a complete astrophysical taxonomy:
    quasars are point-like here and therefore render as ``STAR``.

    Parameters
    ----------
    row : astropy.table.Row
        A single row from an Astropy Table with Gaia DR3 columns.

    Returns
    -------
    str
        ``GALAXY`` for a confidently extended source; otherwise ``STAR``.
    """
    # Endpoint-specific Gaia queries must retain and normalize these exact
    # gaiadr3.gaia_source column names.
    probabilities = [
        _finite_float(row, 'classprob_dsc_combmod_star'),
        _finite_float(row, 'classprob_dsc_combmod_galaxy'),
        _finite_float(row, 'classprob_dsc_combmod_quasar'),
    ]
    if np.all(np.isfinite(probabilities)):
        star_probability, galaxy_probability, quasar_probability = (
            probabilities)
        stype = 'GALAXY' if (
            galaxy_probability >= GAIA_DSC_GALAXY_THRESHOLD
            and galaxy_probability > star_probability
            and galaxy_probability > quasar_probability
        ) else 'STAR'
        logging.info('Type=%s, Gaia DSC probabilities=%s', stype,
                     probabilities)
        return stype

    upstream_type = row['type'] if 'type' in row.colnames else None
    if (not np.ma.is_masked(upstream_type)
            and upstream_type in ('STAR', 'GALAXY')):
        return upstream_type

    parallax = _finite_float(row, 'parallax')
    excess_noise = _finite_float(row, 'astrometric_excess_noise')
    excess_factor = _finite_float(row, 'phot_bp_rp_excess_factor')
    color = _finite_float(row, 'bp_rp')

    # Retain the existing excess-noise and BP/RP-excess thresholds only as
    # legacy extension proxies. Parallax alone is not galaxy evidence.
    if (np.isfinite(excess_noise) and excess_noise > 1.0) or (
            np.isfinite(excess_factor) and np.isfinite(color)
            and excess_factor > (1.0 + 0.015 * color**2)):
        stype = 'GALAXY'
    else:
        stype = 'STAR'

    msg = f"Type={stype}, parallax={parallax}, excess noise="
    msg += f"{excess_noise}, excess_factor={excess_factor}"
    logging.info(msg)

    return stype


def fraction_contaminated(aperture, targframes, starcube, trace_masks=None,
                          collapse_axis=None):
    """
    Calculate the fraction of contamination per spectral trace

    Parameters
    ----------
    aperture: str
        The name of the aperture
    targframes: np.ndarray
        The target simulation frame(s)
    starcube: np.ndarray
        The simulation data cube
    trace_masks : sequence of numpy.ndarray or numpy.ndarray, optional
        Extraction mask for each target trace.  When omitted, masks are loaded
        for ``aperture``.
    collapse_axis : int, optional
        Detector-image axis to average over.  The mode-specific default is the
        cross-dispersion axis.

    Returns
    -------
    list
        The list of fractional contamination arrays
    """
    # Check for 2D
    if starcube.ndim == 2:
        starcube = starcube[None, :, :]

    # Get the trace masks
    if trace_masks is None:
        if aperture == miri_lrs.APERTURE:
            trace_masks = miri_lrs.load_reference_trace().extraction_mask
        else:
            trace_masks = get_trace_mask(aperture, radius=20, plot=False)
    if np.asarray(trace_masks).ndim == 2:
        trace_masks = [trace_masks]
    if collapse_axis is None:
        collapse_axis = (miri_lrs.CROSS_DISPERSION_AXIS
                         if aperture == miri_lrs.APERTURE else 0)

    # Average only pixels inside each trace mask. Multiplying off-mask pixels
    # by zero would leave them finite and incorrectly count them in the mean.
    pctlines = []
    for tframe, mask in zip(targframes, trace_masks):
        # Process one spectral trace at a time.  DHS has ten detector-sized
        # target traces, and retaining all sums and fractions simultaneously
        # multiplies the peak by the number of traces.
        simframe = tframe + starcube
        pframe = np.divide(
            starcube, simframe,
            out=np.full_like(starcube, np.nan),
            where=(simframe != 0) & ~np.isnan(simframe))
        masked = np.where(mask > 0, pframe * mask, np.nan)
        # pframe always has a leading PA/cube axis, so detector axis N is
        # cube axis N+1. Existing SOSS/DHS detector Y remains cube axis 1.
        # Explicitly retain NaN for completely empty spectral channels without
        # emitting NumPy's ``Mean of empty slice`` warning for every PA.
        reduction_axis = collapse_axis + 1
        valid = ~np.isnan(masked)
        numerator = np.nansum(masked, axis=reduction_axis)
        denominator = np.sum(valid, axis=reduction_axis)
        mean_line = np.divide(
            numerator, denominator,
            out=np.full(np.shape(numerator), np.nan, dtype=float),
            where=denominator > 0,
        )
        pctlines.append(mean_line)
        del simframe, pframe, masked, valid, numerator, denominator
    return pctlines


def _project_sources_to_detector(attitude, stars, aperture, aper):
    """Project all non-target sources for one PA using pySIAF array inputs.

    ``JwstAperture.tel_to_det`` constructs the distortion models it needs for
    every call.  Calling it once for an array of sources therefore avoids
    rebuilding identical immutable models for every Gaia row.  The placement
    equations intentionally retain the scalar implementation's independent
    integer conversions, because those conversions define the rendered pixel
    positions.
    """
    if len(stars) <= 1:
        return

    source_slice = slice(1, None)
    ra = np.asarray(stars['ra'][source_slice], dtype=float)
    dec = np.asarray(stars['dec'][source_slice], dtype=float)

    v2, v3 = pysiaf.utils.rotations.sky_to_tel(attitude, ra, dec)
    xtel = np.asarray(v2.to_value(u.arcsec), dtype=float)
    ytel = np.asarray(v3.to_value(u.arcsec), dtype=float)
    xdet, ydet = aperture.tel_to_det(xtel, ytel)
    xsci, ysci = aperture.det_to_sci(xdet, ydet)
    xsci = np.asarray(xsci, dtype=float)
    ysci = np.asarray(ysci, dtype=float)

    stars['xtel'][source_slice] = xtel
    stars['ytel'][source_slice] = ytel
    stars['xdet'][source_slice] = xdet
    stars['ydet'][source_slice] = ydet
    stars['xsci'][source_slice] = xsci
    stars['ysci'][source_slice] = ysci

    target_xsci = stars['xsci'][0]
    target_ysci = stars['ysci'][0]
    target_xord0 = stars['xord0'][0]
    target_yord0 = stars['yord0'][0]

    xord0 = np.asarray(
        xsci + aper['c0x0'] + aper['c1x0'] * (target_xsci - xsci),
        dtype=int)
    yord0 = np.asarray(
        ysci + aper['c0y0'] + aper['c1y0'] * (target_ysci - ysci),
        dtype=int)
    x_shift = np.asarray(aper['c1x1'] * (target_xord0 - xord0), dtype=int)
    y_shift = np.add(
        np.asarray(aper['c1y1'] * (target_yord0 - yord0), dtype=int),
        np.asarray(aper['c2y1'] * (target_xord0 - xord0), dtype=int))

    stars['xord0'][source_slice] = xord0
    stars['yord0'][source_slice] = yord0
    stars['xord1'][source_slice] = xord0 + aper['xord0to1'] + x_shift
    stars['yord1'][source_slice] = yord0 + aper['yord0to1'] + y_shift


def calc_v3pa(V3PA, stars, aperture, data=None, tilt=0, plot=False, POM=False,
              include_target=True):
    """
    Calculate the V3 position angle for each target at the given PA

    Parameters
    ----------
    V3PA: float
        The PA in V3
    stars: astropy.table.Table
        The table of stars in the target vicinity
    aperture: pysiaf.aperture.JwstAperture, str
        The aperture object for the given mode
    data: sequence (optional)
        The data to use instead of making a simulation (to check accuracy or ID sources)
    plot: bool
        Plot the full frame and subarray bounds with all traces

    Returns
    -------
    targframe, starframe
        The frame containing the target trace and a frame containing all contaminating star traces
    """
    logging.info("Checking PA={} with {} stars in the vicinity".format(V3PA, len(stars['ra'])))

    if isinstance(aperture, str):

        logging.info("Getting aperture info from pysiaf...")

        # Aperture names
        if aperture not in APERTURES:
            raise ValueError("Aperture '{}' not supported. Try {}".format(aperture, list(APERTURES.keys())))

        # Instantiate a pySIAF object
        inst = APERTURES[aperture]
        siaf = pysiaf.Siaf(inst['inst'])

        # Get the full and subarray apertures
        full = siaf.apertures[inst['full']]
        aperture = siaf.apertures[aperture]

        # Full frame pixel positions
        rows, cols = full.corners('det')
        aperture.minrow, aperture.maxrow = rows.min(), rows.max()
        aperture.mincol, aperture.maxcol = cols.min(), cols.max()

    if aperture.AperName == miri_lrs.APERTURE and plot:
        result = miri_lrs.calc_v3pa(V3PA, stars, aperture)
        return result, cf.miri_single_pa_plot(result)
    if aperture.AperName == miri_lrs.APERTURE:
        return miri_lrs.calc_v3pa(V3PA, stars, aperture)

    # Convert the renderer's V3PA to the aperture PA required by SIAF.
    APA = aperture_pa_from_v3pa(V3PA, aperture)

    # Aperture info
    aper = APERTURES[aperture.AperName]
    subY, subX = aper['subarr_y'][2] - aper['subarr_y'][1], aper['subarr_x'][1] - aper['subarr_x'][0]

    # Calculate corresponding V2/V3 (TEL) coordinates for Sweetspot
    stars['xdet'][0], stars['ydet'][0] = aperture.reference_point('det')
    stars['xtel'][0], stars['ytel'][0] = aperture.det_to_tel(stars['xdet'][0], stars['ydet'][0])
    stars['xsci'][0], stars['ysci'][0] = aperture.det_to_sci(stars['xdet'][0], stars['ydet'][0])

    # Order 0 location of target relative to pysiaf SCI coordinates
    stars['xord0'][0] = int(stars['xsci'][0] + aper['c0x0'])
    stars['yord0'][0] = int(stars['ysci'][0] + aper['c0y0'])

    # Order 1 location of target relative to order 0
    stars['xord1'][0] = stars['xord0'][0] + aper['xord0to1']
    stars['yord1'][0] = stars['yord0'][0] + aper['yord0to1']

    # Get target's attitude matrix for each` Position Angle
    attitude = pysiaf.utils.rotations.attitude_matrix(stars['xtel'][0], stars['ytel'][0], stars['ra'][0], stars['dec'][0], APA)

    # Get relative coordinates of the stars based on target attitude
    logging.info("Getting star locations for {} stars at PA={} from pysiaf...".format(len(stars), APA))

    _project_sources_to_detector(attitude, stars, aperture, aper)

    logging.info(f'Calculated target and {len(stars)-1} source sci coordinates.')

    # Just sources in FOV (Should always have at least 1, the target)
    lft, rgt, top, bot = aper['lft'], aper['rgt'], aper['top'], aper['bot']
    FOVstars = stars[(lft < stars['xord0']) & (stars['xord0'] < rgt) & (bot < stars['yord0']) & (stars['yord0'] < top)]

    # ``find_sources`` filters Gaia results, but callers may pass a custom
    # table straight to this renderer. Apply the same validation at this
    # boundary before multiplying detector traces by ``fluxscale``.
    FOVstars = _filter_valid_flux_sources(FOVstars, 'fluxscale', 'flux scale')

    # Remove Teff value for GALAXY type
    FOVstars['Teff'] = [np.nan if t == 'GALAXY' else i for i, t in zip(FOVstars['Teff'], FOVstars['type'])]
    FOVstars['Teff_str'] = ['---' if t == 'GALAXY' else str(int(i)) for i, t in zip(FOVstars['Teff'], FOVstars['type'])]

    logging.info("Calculating contamination from {} other sources in the FOV".format(len(FOVstars) - 1))

    # Make frame for the target and a frame for all the other stars
    dhs_mode = 'DHS' in aperture.AperName
    if dhs_mode:
        n_traces = len(_get_trace_template_cached(aperture.AperName)[0])
    else:
        # Existing modes use small precomputed detector images. Keep their
        # established rendering path while DHS uses bounded working buffers.
        star_traces = [
            get_trace(aperture.AperName, temp, typ)
            for temp, typ in zip(FOVstars['Teff'], FOVstars['type'])
        ]
        FOVstars['traces'] = star_traces
        n_traces = len(star_traces[0])
    targframes = ([np.zeros((subY, subX)) for _ in range(n_traces)]
                  if include_target else None)
    starframe = np.zeros((subY, subX))

    # Iterate over all stars in the FOV and add their scaled traces to the correct frame
    for idx, star in enumerate(FOVstars):
        fluxscale = float(star['fluxscale'])
        if dhs_mode:
            traces = _iter_scaled_dhs_traces(
                aperture.AperName, star['Teff'], star['type'])
        else:
            traces = ((trace, None) for trace in star['traces'])

        # Add each target trace to it's own frame
        if idx == 0:
            if include_target:
                for n, (trace, spectral_scale) in enumerate(traces):
                    add_scaled_array_inplace(
                        targframes[n], trace, 0, 0,
                        fluxscale=fluxscale,
                        spectral_scale=spectral_scale)

        # Add all orders to the same frame (if it is a STAR)
        else:

            # Get correct order 0
            order0 = get_order0(aperture.AperName, star['Teff'], stype=star['type'])  # Scaling factor based on observations

            # Scale the order 0 image and add it to the starframe
            scale0 = copy(order0) * fluxscale * aper['empirical_scale'][0]
            starframe = add_array_at_position(starframe, scale0, int(star['xord0'] - aper['subarr_x'][0]), int(star['yord0'] - aper['subarr_y'][1]), centered=True)

            # NOTE: Take this conditional out if you want to see galaxy traces!
            if star['type'] == 'STAR':
                for trace, spectral_scale in traces:
                    add_scaled_array_inplace(
                        starframe, trace,
                        int(star['xord1'] - stars['xord1'][0]),
                        int(star['yord1'] - stars['yord1'][0]),
                        fluxscale=fluxscale,
                        spectral_scale=spectral_scale)

    logging.info(f'Added {len(FOVstars)} sources to the simulated frames.')

    # Make results dict
    result = {
        'pa': V3PA,
        'target': (np.sum(targframes, axis=0)
                   if include_target else None),
        'target_traces': targframes,
        'contaminants': starframe,
    }
    logging.info('Compiled final results.')

    if plot:
        # The interactive single-PA plot includes one contamination line per
        # spectral order. Full simulations reduce these later, but this path
        # must calculate them before constructing the Bokeh data source.
        pctlines = [
            line[0] for line in fraction_contaminated(
                aperture.AperName, targframes, starframe)]

        # Make the plot
        tools = ['pan', 'reset', 'box_zoom', 'wheel_zoom', 'save']
        tips = [('Name', '@name'), ('Type', '@type'), ('RA', '@ra'), ('DEC', '@dec'), ('trace', '@trace'),
                ('scale', '@fluxscale'), ('Teff [K]', '@Teff_str'), ('distance [mas]', '@distance')]
        fig = figure(title='Generated FOV from Gaia EDR3', width=900, height=450 if inst['inst'] == 'NIRCam' else max(subY, 120), match_aspect=True, tools=tools)
        fig.title = '({}, {}) at PA={} in {}'.format(stars[0]['ra'], stars[0]['dec'], V3PA, aperture.AperName)

        # Plot config
        scale = 'log'
        color_map = 'Viridis256'

        # Plot the real or simulated frame
        simframe = data if data is not None else np.sum(targframes, axis=0) + starframe

        # Replace negatives
        simframe[simframe < 0] = 0

        # Rotate for PWCPOS
        # simframe = rotate(simframe, tilt)

        # Plot the image data or simulation
        positive = simframe[simframe > 0]
        vmin = np.percentile(positive, 20)
        vmax = np.percentile(positive, 99.9)
        mapper = LogColorMapper(palette=color_map, low=vmin, high=vmax) if scale == 'log' else LinearColorMapper(palette=color_map, low=vmin, high=vmax)
        imgsource = ColumnDataSource(data={'sim': [simframe]})
        fig.image(image='sim', x=aper['subarr_x'][0], dw=subX, y=aper['subarr_y'][1], dh=subY, source=imgsource, name="image", color_mapper=mapper)

        # Plot the detector gaps and reference pixels for visual inspection
        refframe = NIRCam_detector_gap() if 'NRCA5' in aperture.AperName else np.zeros((subY, subX))
        refsource = ColumnDataSource(data={'ref': [refframe]})
        fig.image(image='ref', x=aper['subarr_x'][0], dw=subX, y=aper['subarr_y'][1], dh=subY, source=refsource, name="ref", color_mapper=LinearColorMapper(palette=["white", "black"], low=0, high=1), alpha=0.1)

        # Plot order 0 locations of stars
        FOVstars_only = FOVstars[FOVstars['type'] == 'STAR']
        source0_stars = ColumnDataSource(data={'Teff_str': FOVstars_only['Teff_str'], 'distance': FOVstars_only['distance'], 'xord0': FOVstars_only['xord0'],
                                         'yord0': FOVstars_only['yord0'], 'ra': FOVstars_only['ra'], 'dec': FOVstars_only['dec'], 'name': FOVstars_only['name'],
                                         'type': FOVstars_only['type'], 'url': FOVstars_only['url'], 'fluxscale': FOVstars_only['fluxscale'],
                                         'trace': ['Order 0'] * len(FOVstars_only)})
        order0_stars = fig.scatter('xord0', 'yord0', color='red', size=20, line_width=3, fill_color=None, name='order0', source=source0_stars)

        # Plot the POM footprint
        if POM:
            fig.varea(x=[lft, rgt], y1=[bot, bot], y2=[top, top], fill_color='blue', fill_alpha=0.1)

        # Plot order 0 locations of galaxies
        FOVstars_gal = FOVstars[FOVstars['type'] == 'GALAXY']
        order0_gal = None
        if len(FOVstars_gal) > 0:
            source0_gal = ColumnDataSource(
                data={'Teff_str': FOVstars_gal['Teff_str'], 'distance': FOVstars_gal['distance'], 'xord0': FOVstars_gal['xord0'],
                      'yord0': FOVstars_gal['yord0'], 'ra': FOVstars_gal['ra'], 'dec': FOVstars_gal['dec'],
                      'name': FOVstars_gal['name'], 'type': FOVstars_gal['type'],
                      'url': FOVstars_gal['url'], 'fluxscale': FOVstars_gal['fluxscale'],
                      'trace': ['Order 0'] * len(FOVstars_gal)})
            order0_gal = fig.scatter('xord0', 'yord0', color='pink', size=20, line_width=3, fill_color=None, name='order0',
                                      source=source0_gal)

        # Plot the target order 0
        fig.scatter([stars[0]['xord0']], [stars[0]['yord0']], size=8, line_width=3, fill_color=None, line_color='black')

        # Get the nominal x and y values for the trace centroids
        n_pts = 1000
        x_ranges = [np.linspace(inst['blue_ext'], inst['cutoffs'][n] + inst['red_ext'], n_pts) for n in np.arange(n_traces)]
        y_ranges = [np.polyval(inst['coeffs'][n], xr) for n, xr in enumerate(x_ranges)]

        lines = []
        for idx, star in enumerate(FOVstars):

            # Trace locations relative to order 0
            for trx in np.arange(n_traces):

                # Make the base dict for this star
                data_dict = {'name': ['Target' if idx == 0 else star['designation']] * n_pts,
                             'type': [star['type']] * n_pts,
                             'ra': [star['ra']] * n_pts,
                             'dec': [star['dec']] * n_pts,
                             'fluxscale': [star['fluxscale']] * n_pts,
                             'Teff_str': [star['Teff_str']] * n_pts,
                             'distance': [star['distance']] * n_pts,
                             'trace': [aper['trace_names'][trx]] * n_pts,
                             'url': [star['url']] * n_pts}

                # Add the offset traces for this star
                for n in np.arange(n_traces):
                    data_dict[f'x{n}'] = x_ranges[n] + star['xord1']
                    data_dict[f'y{n}'] = y_ranges[n] + star['yord1']

                source = ColumnDataSource(data=data_dict)
                line = fig.line('x{}'.format(trx), 'y{}'.format(trx), source=copy(source), color='pink' if star['type'] == 'GALAXY' else 'red', name='traces', line_dash='solid' if idx == 0 else 'dashed', width=3 if idx == 0 else 1)
                lines.append(line)

        # Add order 0 hover and taptool
        if order0_gal is not None:
            fig.add_tools(HoverTool(renderers=[order0_stars, order0_gal], tooltips=tips, name='order0', mode='mouse'))

        fig.add_tools(TapTool(behavior='select', name='order0', callback=OpenURL(url="@url")))

        # Add traces hover and taptool
        fig.add_tools(HoverTool(renderers=lines, tooltips=tips, name='traces', mode='mouse'))
        fig.add_tools(TapTool(behavior='select', name='traces', callback=OpenURL(url="@url")))

        # Show the figure
        pad = 20
        fig.x_range = Range1d(aper['subarr_x'][0] - pad, aper['subarr_x'][1] + pad)
        fig.y_range = Range1d(aper['subarr_y'][1] - pad, aper['subarr_y'][2] + pad)

        # Source for ratio plot
        data = {f'pct_{n}': pct for n, pct in enumerate(pctlines)}
        data.update({'x': np.arange(subX), 'zeros': np.zeros(subX)})
        rsource = ColumnDataSource(data=data)

        # Make plot
        rfig = figure(title='Target Contamination', width=900, height=200, match_aspect=True, tools=tools, x_range=fig.x_range)
        colors = ['blue', 'red', 'green', 'cyan', 'dodgerblue', 'purple', 'orange', 'lime', 'yellow', 'magenta']
        trace_names = inst['trace_names']
        for n in np.arange(n_traces):
            rfig.line('x', f'pct_{n}', color=colors[n], legend_label=trace_names[n], source=copy(rsource))
            glyph = VArea(x='x', y1='zeros', y2=f'pct_{n}', fill_color=colors[n], fill_alpha=0.3)
            rfig.add_glyph(copy(rsource), glyph)
        rfig.y_range = Range1d(0, 1) #min(1, max(pctline_o1.max(), pctline_o2.max(), pctline_o3.max())))
        rfig.yaxis.axis_label = 'Contam / Total Counts'
        rfig.xaxis.axis_label = 'Detector Column'

        # Color bar
        # color_bar = ColorBar(color_mapper=mapper['transform'], width=10, location=(0, 0), title="Teff")
        # fig.add_layout(color_bar, 'right')

        # Plot grid
        gp = gridplot([[fig], [rfig]])

        return result, gp

    return result


def update_task(task, new_state):
    if task is not None:
        task.update_state(state=new_state)


def _compact_dhs_results(aperture, pa_results):
    """Stream DHS detector results into order-by-PA contamination fractions.

    Parameters
    ----------
    aperture : str
        NIRCam DHS aperture name.
    pa_results : iterable of dict
        Results from :func:`calc_v3pa`. The iterable may be a generator, so at
        most one PA's detector images need to exist at a time. Its first item
        must include target traces; later items may omit them.

    Returns
    -------
    tuple
        Target-order detector images and a :class:`DHSContaminationResult`.
    """

    targframes = None
    compact_pctlines = None
    for result in pa_results:
        if targframes is None:
            if result['target_traces'] is None:
                raise ValueError(
                    "The first DHS PA result must include target traces")
            targframes = [
                np.asarray(trace) for trace in result['target_traces']]
            empty_lines = fraction_contaminated(
                aperture, targframes, np.zeros_like(targframes[0]))
            compact_pctlines = [
                np.repeat(line, 360, axis=0) for line in empty_lines]

        pa_lines = fraction_contaminated(
            aperture, targframes, result['contaminants'])
        for order_lines, line in zip(compact_pctlines, pa_lines):
            order_lines[int(result['pa'])] = line[0]

    if targframes is None:
        raise ValueError("At least one DHS PA result is required")

    return targframes, DHSContaminationResult(
        order_fractions=tuple(compact_pctlines),
        position_angles=np.arange(360, dtype=int))


def field_simulation(ra=None, dec=None, aperture=None, targname=None,
                     binComp=None, target_date=None, plot=False, task=None,
                     title='My Target', target_db=None, slider=False):
    """Produce a contamination field simulation at the given sky coordinates

    Parameters
    ----------
    ra : float
        The RA of the target
    dec : float
        The Dec of the target
    aperture: str
        The aperture to use, ['NIS_SUBSTRIP96', 'NIS_SUBSTRIP256', 'NRCA5_GRISM256_F444W', 'NRCA5_GRISM256_F322W2']
    targname: str
        The name of the target to look up in ExoMAST
    binComp : dict
        A dictionary of parameters for a binary companion with keys {'name', 'ra', 'dec', 'fluxscale', 'teff'}
    target_date: Time, int, str
        The target epoch year of the observation, e.g. '2025'
    plot: bool
        Return a plot
    title: str
        The plot title to use
    target_db: str
        The path to the precomputed .h5 database of results
    slider: bool
        Make the PA slider plot instead of the legacy wavelength vs. PA plots
    Returns
    -------
    targframes : list of numpy.ndarray
        One target detector image per spectral order.
    contamination : numpy.ndarray or DHSContaminationResult
        Existing modes return the established dense detector cube. NIRCam DHS
        returns compact fractional contamination by order, with axes described
        by :class:`DHSContaminationResult`.
    position_angles : object
        Observable position-angle metadata, or a plot when ``plot=True``.

    Example
    -------
    from exoctk.contam_visibility import field_simulator as fs
    ra, dec = 91.872242, -25.594934
    targframe, starcube, results = fs.field_simulation(ra, dec, 'NIS_SUBSTRIP256')
    """
    # Aperture names
    if aperture not in APERTURES:
        raise ValueError("Aperture '{}' not supported. Try {}".format(aperture, list(APERTURES.keys())))

    # Check for contam tool data
    check_for_data('exoctk_contam')

    # Bookkeeping
    logging.info("Setting up simulation")
    start = time.time()

    # Resolve target in ExoMAST if possible
    if targname is not None:
        targname = get_canonical_name(targname)
        data, _ = get_target_data(targname)
        ra = data.get('RA')
        dec = data.get('DEC')
        logging.info(f"Resolved '{targname}' (RA={ra}, Dec={dec}) in ExoMAST.")

    # Check to see if there is a precomputed DB for this aperture in the user's
    # environment variables if they don't explicitly supply one as 'target_db'
    if target_db is None:
        db_path = os.environ.get('EXOCTK_CONTAM_CACHE', None)
        if db_path is not None:
            aperture_dbs = glob.glob(os.path.join(db_path, f'{aperture}*.h5'))
            if len(aperture_dbs) > 0:
                target_db = aperture_dbs[0]
            else:
                target_db = os.path.join(db_path, f"{aperture}_db.h5")

    # Check to see if the planet is in the DB
    # Require target_db and targname
    # Require None for binComp and target_date, since these change the results
    precomputed = False
    bounded_dhs = 'DHS' in aperture
    if target_db is not None:
        logging.info(f"Found target DB {target_db}")
        if targname is not None:
            logging.info(f"Target name is {targname}")
            if binComp is None:
                logging.info(f"No binary companion included")
                if target_date is None or str(target_date) == datetime.now().strftime("%Y"):
                    logging.info(f"Looking for target in database")
                    grp_name = get_canonical_name(targname).strip().replace("/", "_")
                    with h5py.File(target_db, "a") as f:
                        if grp_name in f:
                            precomputed = _cached_contam_result_available(
                                f[grp_name], compact_dhs=bounded_dhs)
                else:
                    logging.info("Can't precompute with non-current epoch")
            else:
                logging.info("Can't precompute with binary companion")
        else:
            logging.info(f"Target name {targname} is None")
    else:
        logging.warning(f"Precomputed database {target_db} not found")

    # Grab data from DB if precomputed
    if precomputed:
        targframes, starcube, attrs = fetch_contam_results(targname, target_db)
        if aperture == miri_lrs.APERTURE:
            goodPA_list = miri_lrs.PositionAngleResults(
                successful=attrs['goodPA_list'],
                inaccessible=attrs.get('inaccessiblePA_list', []))
        else:
            goodPA_list = attrs["goodPA_list"]
        logging.info(f'Using precomputed data for {targname} from {target_db}.')

    # Otherwise calculate it now
    else:
        logging.info('Target has not been precomputed. Computing now...')

        # Instantiate a pySIAF object
        logging.info(f'Getting info from pysiaf for {aperture} aperture...')

        targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=u.deg)
        inst = APERTURES[aperture]
        siaf = pysiaf.Siaf(inst['inst'])

        # Get the full and subarray apertures
        full = siaf.apertures[inst['full']]
        aper = siaf.apertures[aperture]
        subX, subY = aper.XSciSize, aper.YSciSize

        # Full frame pixel positions
        rows, cols = full.corners('det')
        aper.minrow, aper.maxrow = rows.min(), rows.max()
        aper.minrow, aper.maxrow = rows.min(), rows.max()
        aper.mincol, aper.maxcol = cols.min(), cols.max()

        # Find stars in the vicinity
        update_task(task, "RUNNING SOURCE QUERY")
        if target_date is None:
            target_date = Time.now()
        stars = find_sources(
            ra, dec, width=source_query_width(aperture),
            target_date=target_date)

        # Add stars manually
        if isinstance(binComp, dict):
            stars = add_source(stars, **binComp)

        # Get full list from ephemeris
        ra_hms, dec_dms = re.sub('[a-z]', ':', targetcrd.to_string('hmsdms')).split(' ')
        goodPAs = get_exoplanet_positions(ra_hms, dec_dms, in_FOR=True)

        # Contamination geometry is evaluated in V3PA. Instrument columns in
        # GTVT are already rotated aperture PAs and must not be converted twice.
        observable_pa_list, good_group_bounds, goodPA_ints = (
            observable_v3pa_ranges(goodPAs))

        # MIRI/LRS is inexpensive enough to evaluate every orientation. Keep
        # the year-specific visibility list separately for plot shading.
        calculation_pa_list = calculation_v3pas(
            aperture, observable_pa_list)

        log_checkpoint(
            f'Found {len(goodPA_ints)}/360 visible position angles')

        # Time it
        logging.info(
            'Calculating target contamination from %d neighboring sources '
            'at %d position angles; observable ranges are %s...',
            len(stars), len(calculation_pa_list), good_group_bounds)

        # Calculate contamination of all stars at each PA
        # -----------------------------------------------
        # To multiprocess, or not to multiprocess. That is the question.
        # Whether 'tis nobler in the code to suffer
        # The slings and arrows of outrageous list comprehensions,
        # Or to take arms against a sea of troubles,
        # And by multiprocessing end them?
        def calculate_pa_results():
            for i, pa in enumerate(calculation_pa_list):
                update_task(
                    task,
                    f"CALCULATING PA {i + 1} OF {len(calculation_pa_list)}")
                logging.info(
                    f"Calculating PA {i + 1} of {len(calculation_pa_list)}")
                yield calc_v3pa(
                    pa, stars=stars, aperture=aper, plot=False,
                    include_target=(not bounded_dhs or i == 0))

        results = []
        if bounded_dhs:
            targframes, starcube = _compact_dhs_results(
                aperture, calculate_pa_results())
        else:
            results = list(calculate_pa_results())

        if aperture == miri_lrs.APERTURE:
            observable = set(map(int, observable_pa_list))
            successful = [int(result['pa']) for result in results]
            goodPA_list = miri_lrs.PositionAngleResults(
                successful=successful,
                inaccessible=(pa for pa in range(360)
                              if pa not in observable))
        else:
            goodPA_list = observable_pa_list

        if not bounded_dhs:
            # We only need one target frame frames
            targframes = [
                np.asarray(trace)
                for trace in results[0]['target_traces']]

            # Make sure starcube is of shape (PA, rows, cols)
            starcube = np.zeros(
                (360, targframes[0].shape[0], targframes[0].shape[1]))

            # Copy good PA results into completed starcube
            for result in results:
                starcube[result['pa'], :, :] = result['contaminants']

        should_cache = all((targname is not None, target_db is not None))
        if should_cache:
            logging.info(f"Saving {targname} to cache {target_db}")
            save_exoplanet_data(
                target_db, targname, aperture, ra, dec, targframes, starcube,
                goodPA_list=goodPA_list)

        # We don't need this anymore
        del results

    logging.info('Contamination calculation complete: {} {}'.format(round(time.time() - start, 3), 's'))

    # Make contam plot
    if plot:

        logging.info("Doing a plot here")

        # MIRI calculates inaccessible orientations too, so its visibility
        # shading is explicit rather than inferred from missing cube planes.
        badPAs = unobservable_v3pas(goodPA_list)

        if aperture == miri_lrs.APERTURE:
            pctlines = fraction_contaminated(
                aperture, targframes, starcube)
            asset = miri_lrs.load_reference_trace()
            contam_plot = cf.contam_slider_plot(
                pctlines, badPAs, wavelength=asset.wavelength,
                trace_names=['MIRI LRS'],
                contamination_labels=['Spectrum'], y_max=0.1)

        # Make slider contam plot
        elif slider or bounded_dhs:
            pctlines = (starcube if bounded_dhs else
                        fraction_contaminated(
                            aperture, targframes, starcube))
            contam_plot = cf.contam_slider_plot(
                pctlines, badPAs, instrument=aperture)

        # Make old contam plot
        else:
            starcube_targ = np.zeros((362, targframes[0].shape[1], targframes[0].shape[0]))
            starcube_targ[0, :, :] = (targframes[0]).T[::-1, ::-1]
            starcube_targ[1, :, :] = (targframes[1]).T[::-1, ::-1]
            starcube_targ[2:, :, :] = starcube.swapaxes(1, 2)[:, ::-1, ::-1]
            contam_plot = cf.contam(starcube_targ, aperture, targetName=title, badPAs=badPAs)

        return targframes, starcube, contam_plot

    return targframes, starcube, goodPA_list


def fetch_contam_results(exoplanet_name, db_filename):
    """
    Load target trace and reconstruct full contamination cube.

    Returns
    -------
    target_trace : ndarray (n_traces, nrows, ncols)
    contamination : ndarray or DHSContaminationResult
        A dense detector cube for legacy modes, or compact order fractions
        for NIRCam DHS.
    attrs : dict (metadata)
    """
    name = get_canonical_name(exoplanet_name)
    grp_name = name.strip().replace("/", "_")

    with h5py.File(db_filename, "r") as f:
        if grp_name not in f:
            raise KeyError(f"Exoplanet '{exoplanet_name}' (canonical: '{grp_name}') not found in {filename}. Available: {list(f.keys())}")

        grp = f[grp_name]
        target_trace = grp["target_trace"][:]
        if ("dhs_order_fractions" in grp
                and "dhs_position_angles" in grp):
            fractions = grp["dhs_order_fractions"]
            contamination = DHSContaminationResult(
                order_fractions=tuple(
                    fractions[order] for order in range(fractions.shape[0])),
                position_angles=grp["dhs_position_angles"][:])
        else:
            stored = grp["contamination"][:]
            plane_index = grp["plane_index"][:]

            # Reconstruct contamination cube
            contamination = np.zeros(
                (360, target_trace.shape[1], target_trace.shape[2]),
                dtype=stored.dtype)
            if len(plane_index) > 0:
                contamination[plane_index] = stored

        attrs = dict(grp.attrs)
        for key in ('wavelength', 'valid_wavelength', 'extraction_mask'):
            if key in grp:
                attrs[key] = grp[key][:]

    return target_trace, contamination, attrs


@lru_cache(maxsize=128)
def _get_order0_cached(aperture, teff, stype):
    """Load an immutable order-zero template for repeated rendering.

    Parameters
    ----------
    aperture: str
        The aperture to use

    Returns
    -------
    np.ndarray
        The 2D order 0 image
    """
    if 'DHS' in aperture:
        trace = np.zeros((50,50))
    else:

        if stype == 'STAR':

            # Get the path to the trace files
            traces_path = os.path.join(os.environ['EXOCTK_DATA'], f'exoctk_contam/order0/NIS_order0_*.npy')

            # Glob the file names
            trace_files = glob.glob(traces_path)

            # Get closest Teff
            teffs = np.array([int(os.path.basename(file).split('_')[-1][:-4]) for file in trace_files])
            trace_file = trace_files[np.argmin((teffs - teff)**2)]
            logging.info(f'Fetching {aperture} {teffs[np.argmin((teffs - teff)**2)]}K trace from {trace_file}')

            # Make frame
            trace = np.load(trace_file)

        else:

            # Get stand-in for galaxy order 0
            gal_path = os.path.join(os.environ['EXOCTK_DATA'], f'exoctk_contam/order0/NIS_gal_order0.npy')
            logging.info('Fetching {} galaxy trace from {}'.format(aperture, gal_path))
            trace = np.load(gal_path)

    trace.setflags(write=False)
    return trace


def get_order0(aperture, teff, stype='STAR'):
    """Get an order-zero image, reusing the cached source template."""

    cache_teff = int(round(float(teff))) if stype == 'STAR' else 0
    return _get_order0_cached(aperture, cache_teff, stype)


def _trace_cache_temperature(teff, stype):
    """Return a stable cache key for a source template temperature."""

    return int(round(float(teff))) if stype == 'STAR' else 0


def get_trace(aperture, teff, stype, plot=False):
    """Get prepared traces for a source.

    DHS detector templates are materialized here for API compatibility.
    Production rendering uses :func:`_iter_scaled_dhs_traces` so the cache
    retains only one base template and compact wavelength scaling vectors.
    """

    if 'DHS' in aperture:
        traces = []
        for trace, spectral_scale in _iter_scaled_dhs_traces(
                aperture, teff, stype):
            prepared = np.array(trace, copy=True)
            if spectral_scale is not None:
                prepared *= spectral_scale[np.newaxis, :]
            prepared.setflags(write=False)
            traces.append(prepared)
    else:
        traces = list(_get_trace_cached(
            aperture, _trace_cache_temperature(teff, stype), stype))

    if plot:
        f = figure(width=900, height=450)
        final = np.sum(traces, axis=0)
        f.image([final], x=APERTURES[aperture]['subarr_x'][0],
                y=APERTURES[aperture]['subarr_y'][1],
                dw=final.shape[1], dh=final.shape[0])
        show(f)

    return traces


@lru_cache(maxsize=4)
def _get_trace_template_cached(aperture):
    """Load and prepare one immutable trace template per aperture."""
    trace_file = os.path.join(
        os.environ['EXOCTK_DATA'],
        f'exoctk_contam/traces/{aperture}.npy')
    data = np.load(trace_file, mmap_mode='r')
    waves = data[:, 0, :]
    traces = data[:, 1:, :]
    # Current DHS templates contain no NaNs. Preserve the mmap views in that
    # normal case; retain the established neighbor-fill behavior for older or
    # user-generated templates that do contain them.
    if np.isnan(traces).any():
        traces = replace_NaNs(traces)
    waves.setflags(write=False)
    traces.setflags(write=False)
    return waves, traces


@lru_cache(maxsize=128)
def _get_dhs_spectral_scales_cached(aperture, teff):
    """Return compact per-column stellar scaling for every DHS trace."""
    waves, _ = _get_trace_template_cached(aperture)
    model = _get_aces_grid().get(
        teff, 5.5, 0, mu1=True, interp=False)
    model_w = np.asarray(model['wave'])
    model_f = np.array(model['flux'], copy=True)
    model_f /= np.trapezoid(model_f, model_w)
    scales = []
    for wave in waves:
        scaled_f = np.interp(wave, model_w, model_f)
        scaled_f /= np.nansum(scaled_f)
        scaled_f.setflags(write=False)
        scales.append(scaled_f)
    return tuple(scales)


def _iter_scaled_dhs_traces(aperture, teff, stype):
    """Yield base DHS traces and compact column scales for one source."""
    if stype == 'GALAXY':
        for trace in get_trace_mask(aperture):
            yield trace, None
        return

    _, traces = _get_trace_template_cached(aperture)
    scales = _get_dhs_spectral_scales_cached(
        aperture, _trace_cache_temperature(teff, stype))
    yield from zip(traces, scales)


@lru_cache(maxsize=128)
def _get_trace_cached(aperture, teff, stype):
    """Get the trace for the given aperture at the given temperature

    Parameters
    ----------
    aperture: str
        The aperture to use
    teff: int
        The temperature [K]
    stype: str
        The source type, ['STAR', 'GALAXY']
    plot: bool
        Plot the trace

    Returns
    -------
    np.ndarray
        The 2D trace
    """

    if stype == 'GALAXY':
        traces = get_trace_mask(aperture)

    else:
        # Get the template trace file, which has the wavelength-dependent throughput encoded
        trace_file = os.path.join(os.environ['EXOCTK_DATA'], f'exoctk_contam/traces/{aperture}.npy')

        # Load the template traces (ntraces, ydim+1, xdim) and replace the NaN values
        data = np.load(trace_file)
        waves = data[:, 0, :] # Wavelengths in the first rows: shape=(ntraces, xdim)
        traces = data[:, 1:, :] # Traces in the rest: shape=(ntraces, ydim, xdim)
        traces = replace_NaNs(traces)

        # Get the normalized stellar model for this source
        model = _get_aces_grid().get(teff, 5.5, 0, mu1=True, interp=False)
        model_w, model_f = model['wave'], model['flux']
        model_f /= np.trapezoid(model_f, model_w)

        # Multiply each template trace by the interpolated stellar model (assumes isowavelength columns)
        for idx, (wave, trace) in enumerate(zip(waves, traces)):
            valid = wave > 0
            scaled_f = np.zeros_like(wave)
            scaled_f[valid] = np.interp(wave[valid], model_w, model_f)
            traces[idx] *= scaled_f[np.newaxis, :]

    immutable_traces = tuple(traces)
    for trace in immutable_traces:
        trace.setflags(write=False)

    return immutable_traces


def old_plot_contamination(targframe_o1, targframe_o2, targframe_o3, starcube, wlims, badPAs=[], title=''):
    """
    Plot the contamination

    Parameters
    ----------
    targframe: np.ndarray
        The frame of target data
    starcube: np.ndarray
        The cube of star data at each PA
    wlims: tuple
        The wavelength min and max
    badPAs: list
        The list of position angles with no visibility

    Returns
    -------
    bokeh.layouts.gridplot
        The contamination figure
    """
    # Data dimensions
    PAs, rows, cols = starcube.shape

    for targframe in [targframe_o1, targframe_o2, targframe_o3]:


        # Remove background values < 1 as it can blow up contamination
        targframe = np.where(targframe < 1, 0, targframe)

        # The width of the target trace
        peak = targframe.max()
        low_lim_col = np.where(targframe > 0.0001 * peak)[1].min()
        high_lim_col = np.where(targframe > 0.0001 * peak)[1].max()

        # Using the starcube of shape (PAs, rows, wave), make a frame of (wave, pa)
        contam = np.zeros([rows, PAs])
        for row in np.arange(rows):

            # Get the
            peakX = np.argmax(targframe[row, :])
            left = peakX - low_lim_col
            right = peakX + high_lim_col

            # Calculate weights
            tr = targframe[row, left:right]
            wt = tr / np.sum(tr**2)
            ww = np.tile(wt, PAs).reshape([PAs, tr.size])

            # Add to contam figure
            contam[row, :] = np.sum(starcube[:, row, left:right] * ww, axis=1, where=~np.isnan(starcube[:, row, left:right] * ww))

        # Log plot contamination, clipping small values
        contam = np.log10(np.clip(contam, 1.e-10, 1.))

    # Hover tool
    hover = HoverTool(tooltips=[("Wavelength", "$x"), ("PA", "$y"), ('Value', '@data')], name='contam')
    tools = ['pan', 'box_zoom', 'crosshair', 'reset', hover]
    trplot = figure(tools=tools, width=600, height=500, title=title, x_range=Range1d(*wlims), y_range=Range1d(0, PAs))

    # Colors
    color_mapper = LinearColorMapper(palette=PuBu[8][::-1][2:], low=-4, high=1)
    color_mapper.low_color = 'white'
    color_mapper.high_color = 'black'

    # Make the trace plot
    source = dict(data=[contam])
    trplot.image(source=source, image='data', x=wlims[0], y=0, dw=wlims[1] - wlims[0], dh=PAs, color_mapper=color_mapper, name='contam')
    trplot.xaxis.axis_label = 'Wavelength (um)'
    trplot.yaxis.axis_label = 'Aperture Position Angle (degrees)'
    color_bar = ColorBar(color_mapper=color_mapper, orientation="horizontal", location=(0, 0))
    trplot.add_layout(color_bar, 'below')

    # Shade bad position angles on the trace plot
    nbadPA = len(badPAs)
    if nbadPA > 0:
        tops = [np.max(badPA) for badPA in badPAs]
        bottoms = [np.min(badPA) for badPA in badPAs]
        left = [wlims[0]] * nbadPA
        right = [wlims[1]] * nbadPA
        trplot.quad(top=tops, bottom=bottoms, left=left, right=right, color='#555555', alpha=0.6)

    # # Make a figure summing the contamination at a given PA
    # sumplot = figure(tools=tools, width=150, height=500, x_range=Range1d(0, 100), y_range=trplot.y_range, title=None)
    # sumplot.line(100 * np.sum(contam >= 0.001, axis=1) / rows, np.arange(PAs) - 0.5, line_color='blue', legend_label='> 0.001')
    # sumplot.line(100 * np.sum(contam >= 0.01, axis=1) / rows, np.arange(PAs) - 0.5, line_color='green', legend_label='> 0.01')
    # sumplot.xaxis.axis_label = '% channels contam.'
    # sumplot.yaxis.major_label_text_font_size = '0pt'

    return trplot#gridplot(children=[[trplot, sumplot]])


if __name__ == '__main__':
    ra, dec = "04 25 29.0162", "-30 36 01.603"  # Wasp 79
    field_simulation(ra, dec, 'NIS_SUBSTRIP256')

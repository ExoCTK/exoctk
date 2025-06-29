#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module to calculate the contamination and visibility of a target on a JWST detector
"""

from copy import copy
from functools import partial
import glob
from multiprocessing import pool, cpu_count, set_start_method
import os
import re
import json
import time
from pkg_resources import resource_filename
import logging
import datetime
from urllib.parse import quote_plus

import astropy.coordinates as crd
from astropy.io import fits
from astropy.table import join
import astropy.units as u
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroquery.irsa import Irsa
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
from astroquery.gaia import Gaia
from bokeh.plotting import figure, show
from bokeh.embed import json_item
from bokeh.layouts import gridplot, column
from bokeh.models import Range1d, LinearColorMapper, LogColorMapper, Label, ColorBar, ColumnDataSource, HoverTool, Slider, CustomJS, VArea, CrosshairTool, TapTool, OpenURL, Span, Legend
from bokeh.palettes import PuBu, Spectral6
from bokeh.transform import linear_cmap
from scipy.ndimage.interpolation import rotate
import numpy as np
import pysiaf
import regions

from ..utils import get_env_variables, check_for_data, add_array_at_position, replace_NaNs
from .new_vis_plot import build_visibility_plot, get_exoplanet_positions
from . import contamination_figure as cf
from exoctk import utils

try:
    set_start_method('spawn')
except RuntimeError:
    pass

import warnings
warnings.filterwarnings("ignore", message="Mean of empty slice")

log_file = 'contam_tool.log'
logging.basicConfig(
    filename=log_file,
    filemode='w',   # <-- THIS forces overwrite on every run
    level=logging.INFO,
    format='%(asctime)s %(message)s',
    force=True
)

_last_time = datetime.datetime.now()
def log_checkpoint(message):
    global _last_time
    now = datetime.datetime.now()
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

Vizier.columns = ["**", "+_r"]
Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # DR2 is default catalog
Gaia.ROW_LIMIT = 100

APERTURES = {'NIS_SOSSFULL': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.8, 2.8],
                              'c0x0': 905, 'c0y0': 1467, 'c1x0': -0.013, 'c1y0': -0.1, 'c1y1': 0.12, 'c1x1': -0.03, 'c2y1': -0.011,
                              'subarr_x': [0, 2048, 2048, 0], 'subarr_y':[0, 0, 2048, 2048], 'trim': [127, 126, 252, 1],
                              'lft': 700, 'rgt': 3022, 'top': 2050, 'bot': 1400, 'blue_ext': -150, 'red_ext': 200,
                              'xord0to1': -2886, 'yord0to1': 68, 'empirical_scale': [1, 1.5, 1.5, 1.5],
                              'cutoffs': [2048, 1820, 1130], 'trace_names': ['Order 1', 'Order 2', 'Order 3'],
                              'coeffs': [[1.68975801e-11, -4.60822060e-08, 4.94623886e-05, -5.93935390e-02, 8.67263818e+01],
                                         [3.95721278e-11, -7.40683643e-08, 6.88340922e-05, -3.68009540e-02, 1.06704335e+02],
                                         [1.06699517e-11, 3.36931077e-08, 1.45570667e-05, 1.69277607e-02, 1.45254339e+02]]},
             'NIS_SUBSTRIP96': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.8, 2.8],
                                'c0x0': 905, 'c0y0': 1467, 'c1x0': -0.013, 'c1y0': -0.1, 'c1y1': 0.12, 'c1x1': -0.03, 'c2y1': -0.011,
                                'subarr_x': [0, 2048, 2048, 0], 'subarr_y':[1792, 1792, 1888, 1888], 'trim': [47, 46, 0, 1],
                                'lft': 700, 'rgt': 3022, 'top': 2050, 'bot': 1400, 'blue_ext': -150, 'red_ext': 200,
                                'xord0to1': -2886, 'yord0to1': 68, 'empirical_scale': [1, 1.5, 1.5, 1.5],
                                'cutoffs': [2048, 1820, 1130], 'trace_names': ['Order 1', 'Order 2', 'Order 3'],
                                'coeffs': [[1.68975801e-11, -4.60822060e-08, 4.94623886e-05, -5.93935390e-02, 8.67263818e+01],
                                           [3.95721278e-11, -7.40683643e-08, 6.88340922e-05, -3.68009540e-02, 1.06704335e+02],
                                           [1.06699517e-11, 3.36931077e-08, 1.45570667e-05, 1.69277607e-02, 1.45254339e+02]]},
             'NIS_SUBSTRIP256': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.8, 2.8],
                                 'c0x0': 905, 'c0y0': 1467, 'c1x0': -0.013, 'c1y0': -0.1, 'c1y1': 0.12, 'c1x1': -0.03, 'c2y1': -0.011,
                                 'subarr_x': [0, 2048, 2048, 0], 'subarr_y':[1792, 1792, 2048, 2048], 'trim': [127, 126, 0, 1],
                                 'lft': 700, 'rgt': 3022, 'top': 2050, 'bot': 1400, 'blue_ext': -150, 'red_ext': 200,
                                 'xord0to1': -2886, 'yord0to1': 68, 'empirical_scale': [1, 1.5, 1.5, 1.5],
                                 'cutoffs': [2048, 1820, 1130], 'trace_names': ['Order 1', 'Order 2', 'Order 3'],
                                 'coeffs': [[1.68975801e-11, -4.60822060e-08, 4.94623886e-05, -5.93935390e-02, 8.67263818e+01],
                                            [3.95721278e-11, -7.40683643e-08, 6.88340922e-05, -3.68009540e-02, 1.06704335e+02],
                                            [1.06699517e-11, 3.36931077e-08, 1.45570667e-05, 1.69277607e-02, 1.45254339e+02]]},
             'NRCA5_40STRIPE1_DHS_F322W2': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.031, 'rad': 2.5, 'lam': [0.8, 2.8],
                                            'subarr_x': [0, 4257, 4257, 0], 'subarr_y': [1064, 1064, 3192, 3192], 'trim': [0, 1, 0, 1],
                                            'c0x0': 1800, 'c0y0': 2116, 'c1x0': 0, 'c1y0': 0, 'c1y1': 0.12, 'c1x1': -0.03, 'c2y1': -0.011,
                                            'lft': 0, 'rgt': 4300, 'top': 4000, 'bot': 0, 'blue_ext': 0, 'red_ext': 0,
                                            'xord0to1': -2300, 'yord0to1': -2116, 'empirical_scale': [1.] * 11,
                                            'cutoffs': [3324]*10, 'trace_names': ['DHS5', 'DHS4', 'DHS3', 'DHS2', 'DHS1', 'DHS6', 'DHS7', 'DHS8', 'DHS9', 'DHS10'],
                                            'coeffs': [[5.31914894e-03, 2.65331915e+03], [5.31914894e-03, 2.54031915e+03],
                                                       [5.31914894e-03, 2.42031915e+03], [5.31914894e-03, 2.28931915e+03],
                                                       [5.31914894e-03, 2.15831915e+03], [5.31914894e-03, 2.03531915e+03],
                                                       [3.54609929e-03, 1.92521277e+03], [3.54609929e-03, 1.81121277e+03],
                                                       [4.43262411e-03, 1.68126596e+03], [5.31914894e-03, 1.54531915e+03]]},
             'NRCA5_40STRIPE1_DHS_F444W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.031, 'rad': 2.5, 'lam': [0.8, 2.8],
                                           'subarr_x': [0, 4257, 4257, 0], 'subarr_y':[1064, 1064, 3192, 3192], 'trim': [0, 1, 0, 1],
                                           'c0x0': 900, 'c0y0': 2116, 'c1x0': 0, 'c1y0': 0, 'c1y1': 0.12, 'c1x1': -0.03, 'c2y1': -0.011,
                                           'lft': 0, 'rgt': 4300, 'top': 4000, 'bot': 0, 'blue_ext': 0, 'red_ext': 0,
                                           'xord0to1': -2000, 'yord0to1': -2116, 'empirical_scale': [1.] * 11,
                                           'cutoffs': [3324]*10, 'trace_names': ['DHS5', 'DHS4', 'DHS3', 'DHS2', 'DHS1', 'DHS6', 'DHS7', 'DHS8', 'DHS9', 'DHS10'],
                                           'coeffs': [[5.31914894e-03, 2.65331915e+03], [5.31914894e-03, 2.54031915e+03],
                                                      [5.31914894e-03, 2.42031915e+03], [5.31914894e-03, 2.28931915e+03],
                                                      [5.31914894e-03, 2.15831915e+03], [5.31914894e-03, 2.03531915e+03],
                                                      [3.54609929e-03, 1.92521277e+03], [3.54609929e-03, 1.81121277e+03],
                                                      [4.43262411e-03, 1.68126596e+03], [5.31914894e-03, 1.54531915e+03]]},
             'NRCA5_GRISM256_F277W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [2.395, 3.179], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F322W2': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [2.413, 4.083], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F356W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [3.100, 4.041], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F444W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [3.835, 5.084], 'trim': [0, 1, 1250, 1]},
             'MIRIM_SLITLESSPRISM': {'inst': 'MIRI', 'full': 'MIRIM_FULL', 'scale': 0.11, 'rad': 2.0, 'lam': [5, 12], 'trim': [6, 5, 0, 1]}}

DHS_STRIPES = {'NRCA5_40STRIPE1_DHS_F322W2': {'DHS5': {'x0': 2196, 'x1': 3324, 'y0': 2665, 'y1': 2671},
                                              'DHS4': {'x0': 2196, 'x1': 3324, 'y0': 2552, 'y1': 2558},
                                              'DHS3': {'x0': 2196, 'x1': 3324, 'y0': 2432, 'y1': 2438},
                                              'DHS2': {'x0': 2196, 'x1': 3324, 'y0': 2301, 'y1': 2307},
                                              'DHS1': {'x0': 2196, 'x1': 3324, 'y0': 2170, 'y1': 2176},
                                              'DHS6': {'x0': 2196, 'x1': 3324, 'y0': 2047, 'y1': 2053},
                                              'DHS7': {'x0': 2196, 'x1': 3324, 'y0': 1933, 'y1': 1937},
                                              'DHS8': {'x0': 2196, 'x1': 3324, 'y0': 1819, 'y1': 1823},
                                              'DHS9': {'x0': 2196, 'x1': 3324, 'y0': 1691, 'y1': 1696},
                                              'DHS10': {'x0': 2196, 'x1': 3324, 'y0': 1557, 'y1': 1563}},
                'NRCA5_40STRIPE1_DHS_F444W': {'DHS5': {'x0': 2196, 'x1': 3324, 'y0': 2665, 'y1': 2671},
                                              'DHS4': {'x0': 2196, 'x1': 3324, 'y0': 2552, 'y1': 2558},
                                              'DHS3': {'x0': 2196, 'x1': 3324, 'y0': 2432, 'y1': 2438},
                                              'DHS2': {'x0': 2196, 'x1': 3324, 'y0': 2301, 'y1': 2307},
                                              'DHS1': {'x0': 2196, 'x1': 3324, 'y0': 2170, 'y1': 2176},
                                              'DHS6': {'x0': 2196, 'x1': 3324, 'y0': 2047, 'y1': 2053},
                                              'DHS7': {'x0': 2196, 'x1': 3324, 'y0': 1933, 'y1': 1937},
                                              'DHS8': {'x0': 2196, 'x1': 3324, 'y0': 1819, 'y1': 1823},
                                              'DHS9': {'x0': 2196, 'x1': 3324, 'y0': 1691, 'y1': 1696},
                                              'DHS10': {'x0': 2196, 'x1': 3324, 'y0': 1557, 'y1': 1563}},
               }

# Gaia color-Teff relation
GAIA_TEFFS = np.asarray(np.genfromtxt(resource_filename('exoctk', 'data/contam_visibility/predicted_gaia_colour.txt'), unpack=True))

def NIRCam_DHS_trace_mask(aperture, gap_value=0, ref_value=0, substripe_value=1, det_value=0, combined=False, plot=False):
    """
    Construct a trace mask for NIRCam DHS mode

    Parameters
    ----------
    aperture: str
        The sperture to use, [
    gap_value: float
        The value in the detector gap
    ref_value: float
        The value of the reference pixels
    substripe_value: float
        The value in the mask area
    combined: bool
        Return a single flattened image
    plot: bool
        Plot the final array

    Returns
    -------
    list, array
        The final array or list of arrays for each trace
    """
    # The DHS uses four detectors
    starting = 0  # pixel
    ref_pixels = 8  # pixel
    det_size_full = 2048  # pixel
    det_size = det_size_full - ref_pixels  # pixel

    # DHS has 2 field points configuration and each field points can be paired to 6 filter position
    # Depending on the field point/filter combination, the spectra will not fall on the same part of the 4 detectors nor have the same length
    stripes = DHS_STRIPES[aperture]
    pixel_scale = 0.031  # arcsec/pixel (on sky)
    gap_fov = 5  # arcsec
    gap = int(gap_fov / pixel_scale)  # pixel

    # full = np.ones((det_size_full + det_size_full + gap, det_size_full + det_size_full + gap))
    full = np.ones((det_size_full + det_size_full + gap, det_size_full + det_size_full + gap)) * det_value

    # now let's populate the full array with the corresponding gap, reference or substripe locations
    full[det_size_full:det_size_full + gap, :] = gap_value
    full[:, det_size_full:det_size_full + gap] = gap_value

    # Add reference pixels
    full[starting:det_size_full, starting:ref_pixels] = ref_value
    full[det_size_full + gap:, starting:ref_pixels] = ref_value
    full[starting:det_size_full, det_size:det_size_full] = ref_value
    full[det_size_full + gap:, det_size:det_size_full] = ref_value
    full[starting:det_size_full, det_size_full + gap:det_size_full + gap + ref_pixels] = ref_value
    full[det_size_full + gap:, det_size_full + gap:det_size_full + gap + ref_pixels] = ref_value
    full[starting:det_size_full, det_size_full + gap + det_size:] = ref_value
    full[det_size_full + gap:, det_size_full + gap + det_size:] = ref_value

    full[starting:ref_pixels, starting:det_size_full] = ref_value
    full[starting:ref_pixels, det_size_full + gap:] = ref_value
    full[det_size:det_size_full, starting:det_size_full] = ref_value
    full[det_size:det_size_full, det_size_full + gap:] = ref_value
    full[det_size_full + gap:det_size_full + gap + ref_pixels, starting:det_size_full] = ref_value
    full[det_size_full + gap:det_size_full + gap + ref_pixels, det_size_full + gap:] = ref_value
    full[det_size_full + gap + det_size:, starting:det_size_full] = ref_value
    full[det_size_full + gap + det_size:, det_size_full + gap:] = ref_value

    # For the specific field point and filter specified above:
    traces = []
    for stripe in stripes.values():
        if combined:
            full[stripe['y0']:stripe['y1'], stripe['x0']:stripe['x1']] += substripe_value
        else:
            trace = copy(full)
            trace[stripe['y0']:stripe['y1'], stripe['x0']:stripe['x1']] = substripe_value
            traces.append(trace)

    if plot:
        plt = figure(width=1000, height=1000)
        plt.image([full], x=0, y=0, dw=full.shape[0], dh=full.shape[1])
        for stripe in stripes.values():
            plt.line([stripe['x0'], stripe['x1']], [(stripe['y0']+stripe['y1'])/2.]*2)
        show(plt)

    return full[1064:3192, :] if combined else [trace[1064:3192] for trace in traces]


def NIRISS_SOSS_trace_mask(aperture, radius=20):
    """
    Construct a trace mask for SOSS data

    Parameters
    ----------
    radius: int
        The radius in pixels of the trace

    Returns
    -------
    np.ndarray
        The SOSS trace mask
    """
    traces = np.array([np.polyval(coeff, np.arange(2048)) for coeff in APERTURES[aperture]['coeffs']])
    ydim = APERTURES[aperture]['subarr_y'][2] - APERTURES[aperture]['subarr_y'][1]
    mask1 = np.zeros((ydim, 2048))
    mask2 = np.zeros((ydim, 2048))
    mask3 = np.zeros((ydim, 2048))
    for col in np.arange(2048):
        mask1[int(traces[0][col]) - radius: int(traces[0][col]) + radius, col] = 1
        mask2[int(traces[1][col]) - radius: int(traces[1][col]) + radius, col] = 1
        mask3[int(traces[2][col]) - radius: int(traces[2][col]) + radius, col] = 1

    # Right referecnce pixels
    mask1[:, :4] = 0

    # Left reference pixels
    mask1[:, -4:] = 0
    mask2[:, :4] = 0
    mask3[:, :4] = 0

    # Top reference pixels
    mask1[-5:, :] = 0
    mask2[-5:, :] = 0
    mask3[-5:, :] = 0

    # Order 3 cutoff
    mask3[:, 823:] = 0

    return mask1, mask2, mask3


def find_sources(ra, dec, width=7.5*u.arcmin, catalog='Gaia', target_date=Time.now(), verbose=False, pm_corr=True):
    """
    Find all the stars in the vicinity and estimate temperatures

    Parameters
    ----------
    ra : float
        The RA of the target in decimal degrees
    dec : float
        The Dec of the target in decimal degrees
    width: astropy.units.quantity
        The width of the square search box
    catalog: str
        The name of the catalog to use, ['2MASS', 'Gaia']
    target_date: Time, int, str
        The target epoch year of the observation, e.g. '2025'
    verbose: bool
        Print details
    pm_corr: bool
        Correct source coordinates based on their proper motion

    Returns
    -------
    astropy.table.Table
        The table of stars
    """
    # Converting to degrees and query for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=u.deg if isinstance(ra, float) and isinstance(dec, float) else (u.hour, u.deg))

    # Search Gaia for stars
    if catalog == 'Gaia':

        if verbose:
            print('Searching {} Catalog to find all stars within {} of RA={}, Dec={}...'.format(catalog, width, ra, dec))

        stars = Gaia.query_object_async(coordinate=targetcrd, width=width, height=width)

        # Perform XMatch between Gaia and SDSS DR16
        xmatch_result = XMatch.query(cat1=stars, cat2='vizier:V/154/sdss16', max_distance=2 * u.arcsec, colRA1='ra', colDec1='dec', colRA2='RA_ICRS', colDec2='DE_ICRS')

        try:
            # Join Gaia results with XMatch results based on source_id
            merged_results = join(stars, xmatch_result, keys='source_id', join_type='left')

            # Extract SDSS DR16 source types from XMatch results
            stars['type'] = ['STAR' if source_id not in xmatch_result['source_id'] else ('STAR' if sdss_type == '' else sdss_type) for source_id, sdss_type in zip(stars['source_id'], merged_results['spCl'])]

        except ValueError:
            pass

        # Or infer galaxy from parallax
        stars['type'] = [classify_source(row) for row in stars]

        # Derived from K. Volk
        stars['Teff'] = [GAIA_TEFFS[0][(np.abs(GAIA_TEFFS[1] - row['bp_rp'])).argmin()] for row in stars]

        # Calculate relative flux
        stars['fluxscale'] = stars['phot_g_mean_flux'] / stars['phot_g_mean_flux'][0]

        # Star names
        stars['name'] = [str(i) for i in stars['source_id']]

        # Catalog name
        cat = 'I/350/gaiaedr3'

    # Search 2MASS
    elif catalog == '2MASS':
        stars = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone', radius=width)

        jhMod = np.array([0.545, 0.561, 0.565, 0.583, 0.596, 0.611, 0.629, 0.642, 0.66, 0.679, 0.696, 0.71, 0.717, 0.715, 0.706, 0.688, 0.663, 0.631, 0.601, 0.568, 0.537, 0.51, 0.482, 0.457, 0.433, 0.411, 0.39, 0.37, 0.314, 0.279])
        hkMod = np.array([0.313, 0.299, 0.284, 0.268, 0.257, 0.247, 0.24, 0.236, 0.229, 0.217,0.203, 0.188, 0.173, 0.159, 0.148, 0.138, 0.13, 0.123, 0.116, 0.112, 0.107, 0.102, 0.098, 0.094, 0.09, 0.086, 0.083, 0.079, 0.07, 0.067])
        teffMod = np.array([2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5800, 6000])

        # Make sure colors are calculated
        stars['j_h'] = stars['j_m'] - stars['h_m']
        stars['h_k'] = stars['h_m'] - stars['k_m']
        stars['j_k'] = stars['j_m'] - stars['k_m']

        # Find Teff of each star from the color
        stars['Teff'] = [teffMod[np.argmin((row['j_h'] - jhMod) ** 2 + (row['h_k'] - hkMod) ** 2)] for row in stars]

        # Calculate relative flux
        stars['fluxscale'] = 10.0 ** (-0.4 * (stars['j_m'] - stars['j_m'][0]))

        # Determine source type (not sure how to do this with 2MASS)
        stars['source_type'] = ['STAR'] * len(stars)

        # Star names
        stars['name'] = [str(i) for i in stars['designation']]

        # Catalog name
        cat = 'II/246/out'

    # Add URL (before PM correcting coordinates)
    search_radius = 1
    urls = ['https://vizier.u-strasbg.fr/viz-bin/VizieR-5?-ref=VIZ62fa613b20f3fc&-out.add=.&-source={}&-c={}&eq=ICRS&rs={}&-out.orig=o'.format(cat, quote_plus(f"{ra_deg} {dec_deg}"), search_radius) for ra_deg, dec_deg in zip(stars['ra'], stars['dec'])]
    stars.add_column(urls, name='url')

    # Cope coordinates to new column
    stars.add_column(stars['ra'], name='obs_ra')
    stars.add_column(stars['dec'], name='obs_dec')

    # Update RA and Dec using proper motion data
    if pm_corr:
        for row in stars:
            new_ra, new_dec = calculate_current_coordinates(row['ra'], row['dec'], row['pmra'], row['pmdec'], row['ref_epoch'], target_date=target_date, verbose=verbose)

            if not hasattr(new_ra, 'mask'):
                row['ra'] = new_ra
            if not hasattr(new_dec, 'mask'):
                row['dec'] = new_dec

    # Find distance from target to each star
    sindRA = (stars['ra'][0] - stars['ra']) * np.cos(stars['dec'][0])
    cosdRA = stars['dec'][0] - stars['dec']
    stars.add_column(np.sqrt(sindRA ** 2 + cosdRA ** 2) * u.deg.to(u.arcsec), name='distance')
    stars.sort('distance')

    # Add detector location to the table
    stars.add_columns(np.zeros((10, len(stars))), names=['xtel', 'ytel', 'xdet', 'ydet', 'xsci', 'ysci', 'xord0', 'yord0', 'xord1', 'yord1'])

    return stars


def calculate_current_coordinates(ra, dec, pm_ra, pm_dec, epoch, target_date=Time.now(), verbose=False):
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
    verbose: bool
        Print extra stuff

    Returns
    -------
    new_ra, new_dec
        The corrected RA and Dec for the source
    """
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

    if verbose:
        print(f"delta_T = {dt_years} years")
        print(f'RA: {ra} + {pm_ra} mas/yr => {new_ra}')
        print(f'Dec: {dec} + {pm_dec} mas/yr => {new_dec}')

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


def classify_source(row, verbose=False):
    """
    Classify a Gaia EDR3 source as STAR or GALAXY based on proxy criteria.

    Parameters
    ----------
    row : astropy.table.Row
        A single row from an Astropy Table with Gaia EDR3 columns.
    verbose: bool
        Print helpful stuff

    Returns
    -------
    str
        'STAR' if the source appears to be stellar, 'GALAXY' otherwise.
    """
    # Default is STAR
    stype = 'STAR'

    # Check to see if GALAXY
    if hasattr(row['parallax'], 'mask'):

        if hasattr(row['astrometric_excess_noise'], 'mask'):
            if row['astrometric_excess_noise'] > 1.0:
                stype = 'GALAXY'

        if hasattr(row['phot_bp_rp_excess_factor'], 'mask') and hasattr(row['bp_rp'], 'mask'):
            if row['phot_bp_rp_excess_factor'] > (1.0 + 0.015 * row['bp_rp']**2):
                stype = 'GALAXY'

    else:
        if row['parallax'] < 0.5:
            stype = 'GALAXY'

    if verbose:
        print(stype, row['parallax'], row['astrometric_excess_noise'], row['phot_bp_rp_excess_factor'])

    return stype


def calc_v3pa(V3PA, stars, aperture, data=None, tilt=0, plot=False, verbose=False, logging=True, POM=False):
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
    verbose: bool
        Print statements
    POM: bool
        Show the pick off mirror footprint in the plot

    Returns
    -------
    targframe, starframe
        The frame containing the target trace and a frame containing all contaminating star traces
    """
    if logging:
        open('contam_tool.log', 'w').close()
        start_time = datetime.datetime.now()
        log_checkpoint(f'Logging in {log_file}...')

    if verbose:
        print("Checking PA={} with {} stars in the vicinity".format(V3PA, len(stars['ra'])))

    if isinstance(aperture, str):

        if verbose:
            print("Getting aperture info from pysiaf...")

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

    if logging:
        log_checkpoint('Loaded aperture info from pySIAF.')

    # Get APA from V3PA
    APA = V3PA + aperture.V3IdlYAngle
    if APA > 360:
        APA = APA - 360
    elif APA < 0:
        APA = APA + 360

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
    if verbose:
        print("Getting star locations for {} stars at PA={} from pysiaf...".format(len(stars), APA))

    for idx, star in enumerate(stars[1:]):

        # Get the TEL coordinates (V2, V3) of the star
        V2, V3 = pysiaf.utils.rotations.sky_to_tel(attitude, star['ra'], star['dec'])
        star['xtel'], star['ytel'] = V2.to(u.arcsec).value, V3.to(u.arcsec).value

        # Get the DET coordinates of the star
        star['xdet'], star['ydet'] = aperture.tel_to_det(star['xtel'], star['ytel'])

        # Get the DET coordinates of the star
        star['xsci'], star['ysci'] = aperture.det_to_sci(star['xdet'], star['ydet'])

        # Order 0 location relative to pysiaf SCI coordinates (with distortion corrections)
        star['xord0'] = int(star['xsci'] + aper['c0x0'] + aper['c1x0'] * (stars['xsci'][0] - star['xsci']))
        star['yord0'] = int(star['ysci'] + aper['c0y0'] + aper['c1y0'] * (stars['ysci'][0] - star['ysci']))

        # Order 1/2/3 location relative to order 0 location (with distortion corrections)
        x_shift = int(aper['c1x1'] * (stars[0]['xord0'] - star['xord0']))
        y_shift = int(aper['c1y1'] * (stars[0]['yord0'] - star['yord0'])) + int(aper['c2y1'] * (stars[0]['xord0'] - star['xord0']))
        star['xord1'] = star['xord0'] + aper['xord0to1'] + x_shift
        star['yord1'] = star['yord0'] + aper['yord0to1'] + y_shift

    if logging:
        log_checkpoint(f'Calculated target and {len(stars)-1} source sci coordinates.')

    # Just sources in FOV (Should always have at least 1, the target)
    lft, rgt, top, bot = aper['lft'], aper['rgt'], aper['top'], aper['bot']
    FOVstars = stars[(lft < stars['xord0']) & (stars['xord0'] < rgt) & (bot < stars['yord0']) & (stars['yord0'] < top)]

    # Get the traces for sources in the FOV and add the column to the source table
    star_traces = [get_trace(aperture.AperName, temp, typ, verbose=False) for temp, typ in zip(FOVstars['Teff'], FOVstars['type'])]
    FOVstars['traces'] = star_traces

    # Remove Teff value for GALAXY type
    FOVstars['Teff'] = [np.nan if t == 'GALAXY' else i for i, t in zip(FOVstars['Teff'], FOVstars['type'])]
    FOVstars['Teff_str'] = ['---' if t == 'GALAXY' else str(int(i)) for i, t in zip(FOVstars['Teff'], FOVstars['type'])]

    if logging:
        log_checkpoint(f'Found {len(FOVstars)} sources in the FOV.')

    if verbose:
        print("Calculating contamination from {} other sources in the FOV".format(len(FOVstars) - 1))

    # Make frame for the target and a frame for all the other stars
    n_traces = len(star_traces[0])
    targframes = [np.zeros((subY, subX))] * n_traces
    starframe = np.zeros((subY, subX))

    # Get order 0
    order0 = get_order0(aperture.AperName) * 1.5e8 # Scaling factor based on observations

    # Get trace masks
    trace_masks = NIRCam_DHS_trace_mask(aperture.AperName) if 'NRCA5' in aperture.AperName else NIRISS_SOSS_trace_mask(aperture.AperName)

    if logging:
        log_checkpoint(f'Fetched arrays, masks, and traces...')

    # Iterate over all stars in the FOV and add their scaled traces to the correct frame
    for idx, star in enumerate(FOVstars):

        # Scale the traces for this source
        traces = [trace * star['fluxscale'] * aper['empirical_scale'][n + 1] for n, trace in enumerate(star['traces'])]

        # Add each target trace to it's own frame
        if idx == 0:

            for n, trace in enumerate(traces):

                # Assumes the lower lft corner of the trace is in the lower left corner of the 'targframe' array
                targframes[n] = add_array_at_position(targframes[n], trace, 0, 0)

        # Add all orders to the same frame (if it is a STAR)
        else:

            # Scale the order 0 image and add it to the starframe
            scale0 = copy(order0) * star['fluxscale'] * aper['empirical_scale'][0]
            starframe = add_array_at_position(starframe, scale0, int(star['xord0'] - aper['subarr_x'][0]), int(star['yord0'] - aper['subarr_y'][1]), centered=True)

            # NOTE: Take this conditional out if you want to see galaxy traces!
            if star['type'] == 'STAR':
                for trace in traces:
                    starframe = add_array_at_position(starframe, trace, int(star['xord1'] - stars['xord1'][0]), int(star['yord1'] - stars['yord1'][0]))

    if logging:
        log_checkpoint(f'Added {len(FOVstars)} sources to the simulated frames.')

    # Adding frames together
    simframes = [tframe + starframe for tframe in targframes]
    simframe = np.sum(targframes, axis=0) + starframe
    pctframes = [np.divide(starframe, sframe, out=np.full_like(starframe, np.nan), where=(sframe != 0) & ~np.isnan(sframe)) for sframe in simframes]
    pctlines = []
    for i, (pframe, mask) in enumerate(zip(pctframes, trace_masks)):
        masked = pframe * mask
        with np.errstate(invalid='ignore', divide='ignore'):
            mean_line = np.nanmean(masked, axis=0)
        pctlines.append(mean_line)

    # Make results dict
    result = {'pa': V3PA, 'target': np.sum(targframes, axis=0), 'target_traces': targframes,
              'contaminants': starframe, 'sources': FOVstars, 'contam_levels': pctlines}

    if logging:
        log_checkpoint('Compiled final results.')

    if plot:

        # Make the plot
        tools = ['pan', 'reset', 'box_zoom', 'wheel_zoom', 'save']
        tips = [('Name', '@name'), ('Type', '@type'), ('RA', '@ra'), ('DEC', '@dec'), ('trace', '@trace'),
                ('scale', '@fluxscale'), ('Teff [K]', '@Teff_str'), ('distance [mas]', '@distance')]
        fig = figure(title='Generated FOV from Gaia EDR3', width=900, height=450 if inst['inst'] == 'NIRCam' else max(subY, 120), match_aspect=True, tools=tools)
        fig.title = '({}, {}) at PA={} in {}'.format(stars[0]['ra'], stars[0]['dec'], V3PA, aperture.AperName)

        # Plot config
        scale = 'log'
        color_map = 'Viridis256'

        # Plot the obs data if possible...
        if data is not None:
            simframe = data

        # Replace negatives
        simframe[simframe < 0] = 0

        # Rotate for PWCPOS
        simframe = rotate(simframe, tilt)

        # Plot the image data or simulation
        vmax = np.nanmax(simframe)
        mapper = LogColorMapper(palette=color_map, low=1, high=vmax) if scale == 'log' else LinearColorMapper(palette=color_map, low=0, high=vmax)
        imgsource = ColumnDataSource(data={'sim': [simframe]})
        fig.image(image='sim', x=aper['subarr_x'][0], dw=subX, y=aper['subarr_y'][1], dh=subY, source=imgsource, name="image", color_mapper=mapper)

        # Plot the detector gaps and reference pixels for visual inspection
        refframe = NIRCam_DHS_trace_mask(aperture.AperName, substripe_value=0, ref_value=1, gap_value=1, combined=True) if 'NRCA5' in aperture.AperName else np.zeros((subY, subX))
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

        if logging:
            log_checkpoint(f'Finished making plot of size {len(json.dumps(json_item(gp))) / 1e6:.2f}MB')
            log_checkpoint(f'Total calculation time: {(datetime.datetime.now() - start_time).total_seconds():.2f}s')

        return result, gp

    return result


def field_simulation(ra, dec, aperture, binComp=None, target_date=Time.now(), n_jobs=-1, interpolate=False, plot=False,
                     title='My Target', multi=True, verbose=True):
    """Produce a contamination field simulation at the given sky coordinates

    Parameters
    ----------
    ra : float
        The RA of the target
    dec : float
        The Dec of the target
    aperture: str
        The aperture to use, ['NIS_SUBSTRIP96', 'NIS_SUBSTRIP256', 'NRCA5_GRISM256_F444W', 'NRCA5_GRISM256_F322W2', 'MIRI_SLITLESSPRISM']
    binComp : dict
        A dictionary of parameters for a binary companion with keys {'name', 'ra', 'dec', 'fluxscale', 'teff'}
    target_date: Time, int, str
        The target epoch year of the observation, e.g. '2025'
    n_jobs: int
        Number of cores to use (-1 = All)
    interpolate: bool
        Skip every other PA and interpolate to speed up calculation
    plot: bool
        Return a plot

    Returns
    -------
    simuCube : np.ndarray
        The simulated data cube. Index 0 and 1 (axis=0) show the trace of
        the target for orders 1 and 2 (respectively). Index 2-362 show the trace
        of the target at every position angle (PA) of the instrument.
    plt: NoneType, bokeh.plotting.figure
        The plot of the contaminationas a function of PA

    Example
    -------
    from exoctk.contam_visibility import field_simulator as fs
    ra, dec = 91.872242, -25.594934
    targframe, starcube, results = fs.field_simulation(ra, dec, 'NIS_SUBSTRIP256')
    """
    # Initialize logging
    open('contam_tool.log', 'w').close()
    start_time = datetime.datetime.now()
    log_checkpoint(f'Logging in {log_file}...')

    # Check for contam tool data
    check_for_data('exoctk_contam')

    # Aperture names
    if aperture not in APERTURES:
        raise ValueError("Aperture '{}' not supported. Try {}".format(aperture, list(APERTURES.keys())))

    # Instantiate a pySIAF object
    if verbose:
        print('Getting info from pysiaf for {} aperture...'.format(aperture))

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
    stars = find_sources(ra, dec, target_date=target_date, verbose=False)

    # Add stars manually
    if isinstance(binComp, dict):
        stars = add_source(stars, **binComp)

    # Set the number of cores for multiprocessing
    max_cores = 8
    if n_jobs == -1 or n_jobs > max_cores:
        n_jobs = max_cores

    log_checkpoint(f'Found {len(stars)} sources in the neighborhood')

    # Get full list from ephemeris
    ra_hms, dec_dms = re.sub('[a-z]', ':', targetcrd.to_string('hmsdms')).split(' ')
    goodPAs = get_exoplanet_positions(ra_hms, dec_dms, in_FOR=True)

    # Get all observable PAs and convert to ints
    goodPA_vals = list(goodPAs[~goodPAs['{}_min_pa_angle'.format(inst['inst'].upper())].isna()]['{}_min_pa_angle'.format(inst['inst'].upper())]) + list(goodPAs[~goodPAs['{}_nominal_angle'.format(inst['inst'].upper())].isna()]['{}_nominal_angle'.format(inst['inst'].upper())]) + list(goodPAs[~goodPAs['{}_max_pa_angle'.format(inst['inst'].upper())].isna()]['{}_max_pa_angle'.format(inst['inst'].upper())])
    goodPA_ints = np.sort(np.unique(np.array(goodPA_vals).astype(int)))

    # Group good PAs to find gaps in visibility
    good_groups = []
    current_group = [goodPA_ints[0]]
    max_gap = 7 # Biggest PA gap considered to still be observable
    for i in range(1, len(goodPA_ints)):
        if goodPA_ints[i] - current_group[-1] <= max_gap:
            current_group.append(goodPA_ints[i])
        else:
            good_groups.append(current_group)
            current_group = [goodPA_ints[i]]

    good_groups.append(current_group)
    good_group_bounds = [(min(grp), max(grp)) for grp in good_groups]
    goodPA_list = np.concatenate([np.arange(grp[0], grp[1]+1) for grp in good_group_bounds]).ravel()

    # Try to speed it up by doing every other visible PA and then interpolating
    # if interpolate:


    log_checkpoint(f'Found {len(goodPA_ints)}/360 visible position angles to check')

    # Flatten list and check against 360 angles to get all bad PAs
    # badPA_list = [pa for pa in pa_list if pa not in goodPA_list]

    # Time it
    if verbose:
        print('Calculating target contamination from {} neighboring sources in position angle ranges {}...'.format(len(stars), [(int(rng[0]), int(rng[1])) for rng in good_group_bounds]))
        start = time.time()

    # Calculate contamination of all stars at each PA
    # -----------------------------------------------
    # To multiprocess, or not to multiprocess. That is the question.
    # Whether 'tis nobler in the code to suffer
    # The slings and arrows of outrageous list comprehensions,
    # Or to take arms against a sea of troubles,
    # And by multiprocessing end them?
    if multi:
        pl = pool.ThreadPool(n_jobs)
        func = partial(calc_v3pa, stars=stars, aperture=aper, plot=False, verbose=False, logging=False)
        results = pl.map(func, goodPA_list)
        pl.close()
        pl.join()

    else:
        results = []
        for pa in goodPA_list:
            result = calc_v3pa(pa, stars=stars, aperture=aper, plot=False, verbose=False, logging=False)
            results.append(result)

    log_checkpoint('Finished all calc_v3pa calculations')

    # We only need one target frame frames
    targframes = [np.asarray(trace) for trace in results[0]['target_traces']]

    # Make sure starcube is of shape (PA, rows, cols)
    starcube = np.zeros((360, targframes[0].shape[0], targframes[0].shape[1]))

    # Make the contamination plot
    for result in results:
        starcube[result['pa'], :, :] = result['contaminants']

    if verbose:
        print('Contamination calculation complete: {} {}'.format(round(time.time() - start, 3), 's'))

    log_checkpoint('Bundled up all the results.')

    # Make contam plot
    if plot:

        # Get bad PA list from missing angles between 0 and 360
        badPAs = [j for j in np.arange(0, 360) if j not in goodPA_list]

        # Make old contam plot
        starcube_targ = np.zeros((362, 2048, 96 if aperture == 'NIS_SUBSTRIP96' else 256))
        starcube_targ[0, :, :] = (targframes[0]).T[::-1, ::-1]
        starcube_targ[1, :, :] = (targframes[1]).T[::-1, ::-1]
        starcube_targ[2:, :, :] = starcube.swapaxes(1, 2)[:, ::-1, ::-1]
        contam_plot = cf.contam(starcube_targ, aperture, targetName=title, badPAs=badPAs)

        # # Slider plot
        # contam_plot = contam_slider_plot(results, plot=False)

        return targframes, starcube, results, contam_plot

    else:

        return targframes, starcube, results


def contam_slider_plot(contam_results, threshold=0.05, plot=False):
    """
    Make the contamination plot with a slider

    Parameters
    ----------
    contam_results: dict
        The dictionary of results from the field_simulation function
    plot: bool
        Show the plot if True

    Returns
    -------
    bokeh.layouts.column
        The column of plots
    """
    # Full PA list
    pa_list = np.arange(360)
    goodPA_list = [result['pa'] for result in contam_results]
    badPA_list = [pa for pa in pa_list if pa not in goodPA_list]

    # Grab one target frame
    targframe = np.asarray(contam_results[0]['target'])

    # Make the contamination plot
    order1_contam = np.zeros((360, targframe.shape[1]))
    order2_contam = np.zeros((360, targframe.shape[1]))
    order3_contam = np.zeros((360, targframe.shape[1]))
    for result in contam_results:
        order1_contam[result['pa'], :] = result['order1_contam']
        order2_contam[result['pa'], :] = result['order2_contam']
        order3_contam[result['pa'], :] = result['order3_contam']

    # Define data
    contam_dict = {'contam1_{}'.format(result['pa']): result['order1_contam'] for result in contam_results}
    contam_dict.update({'contam2_{}'.format(result['pa']): result['order2_contam'] for result in contam_results})
    contam_dict.update({'contam3_{}'.format(result['pa']): result['order3_contam'] for result in contam_results})

    # Wrap the data in two ColumnDataSources
    source_visible = ColumnDataSource(
        data=dict(col=np.arange(2048), zeros=np.zeros(2048), contam1=order1_contam[0], contam2=order2_contam[0],
                  contam3=order3_contam[0]))
    source_available = ColumnDataSource(data=contam_dict)

    # Define plot elements
    plt = figure(width=900, height=300, tools=['reset', 'save'])
    plt.line('col', 'contam1', source=source_visible, color='blue', line_width=2, line_alpha=0.6,
             legend_label='Order 1')
    plt.line('col', 'contam2', source=source_visible, color='red', line_width=2, line_alpha=0.6, legend_label='Order 2')
    plt.line('col', 'contam3', source=source_visible, color='green', line_width=2, line_alpha=0.6,
             legend_label='Order 3')
    glyph1 = VArea(x="col", y1="zeros", y2="contam1", fill_color="blue", fill_alpha=0.3)
    plt.add_glyph(source_visible, glyph1)
    glyph2 = VArea(x="col", y1="zeros", y2="contam2", fill_color="red", fill_alpha=0.3)
    plt.add_glyph(source_visible, glyph2)
    glyph3 = VArea(x="col", y1="zeros", y2="contam3", fill_color="green", fill_alpha=0.3)
    plt.add_glyph(source_visible, glyph3)
    plt.y_range = Range1d(0, min(1, max(np.nanmax(order1_contam), np.nanmax(order2_contam), np.nanmax(order3_contam))))
    plt.x_range = Range1d(0, 2048)
    plt.xaxis.axis_label = ''
    plt.yaxis.axis_label = 'Contamination / Target Flux'
    slider = Slider(title='Position Angle',
                    value=pa_list[0],
                    start=min(pa_list),
                    end=max(pa_list),
                    step=int((max(pa_list) - min(pa_list)) / (len(pa_list) - 1)),
                    sizing_mode='stretch_width')

    span = Span(line_width=2, location=slider.value, dimension='height')

    # Define CustomJS callback, which updates the plot based on selected function by updating the source_visible
    callback = CustomJS(
        args=dict(source_visible=source_visible, source_available=source_available, span=span), code="""
            var selected_pa = (cb_obj.value).toString();
            var data_visible = source_visible.data;
            var data_available = source_available.data;
            data_visible['contam1'] = data_available['contam1_' + selected_pa];
            data_visible['contam2'] = data_available['contam2_' + selected_pa];
            data_visible['contam3'] = data_available['contam3_' + selected_pa];
            span.location = cb_obj.value;
            source_visible.change.emit();
        """)

    # Make a guide that shows which PAs are unobservable
    viz_none = np.array([1 if i in badPA_list else 0 for i in pa_list])
    viz_ord1 = np.array([1 if i > threshold else 0 for i in np.nanmax(order1_contam, axis=1)])
    viz_ord2 = np.array([1 if i > threshold else 0 for i in np.nanmax(order2_contam, axis=1)])
    viz_ord3 = np.array([1 if i > threshold else 0 for i in np.nanmax(order3_contam, axis=1)])

    # Make the plot
    viz_plt = figure(width=900, height=300, x_range=Range1d(0, 359))
    viz_plt.step(np.arange(360), np.mean(order1_contam, axis=1), color='blue', mode="center")
    viz_plt.step(np.arange(360), np.mean(order2_contam, axis=1), color='red', mode="center")
    viz_plt.step(np.arange(360), np.mean(order3_contam, axis=1), color='green', mode="center")
    c1 = viz_plt.vbar(x=np.arange(360), top=viz_ord1, line_color=None, fill_color='blue', alpha=0.2)
    c2 = viz_plt.vbar(x=np.arange(360), top=viz_ord2, line_color=None, fill_color='red', alpha=0.2)
    c3 = viz_plt.vbar(x=np.arange(360), top=viz_ord3, line_color=None, fill_color='green', alpha=0.2)
    c0 = viz_plt.vbar(x=np.arange(360), top=viz_none, width=1.1, line_color=None, fill_color='#555555')
    viz_plt.x_range = Range1d(0, 359)
    viz_plt.y_range = Range1d(0, 1)
    viz_plt.add_layout(span)

    legend = Legend(items=[
        ('Ord 1 > {}% Contaminated'.format(threshold * 100), [c1]),
        ('Ord 2 > {}% Contaminated'.format(threshold * 100), [c2]),
        ('Ord 3 > {}% Contaminated'.format(threshold * 100), [c3]),
        ("Target not visible", [c0]),
    ], location=(50, 0), orientation='horizontal', border_line_alpha=0)
    viz_plt.add_layout(legend, 'below')

    # Put plot together
    slider.js_on_change('value', callback)
    layout = column(plt, slider, viz_plt)

    if plot:
        show(layout)

    return layout


def get_order0(aperture):
    """Get the order 0 image for the given aperture

    Parameters
    ----------
    aperture: str
        The aperture to use

    Returns
    -------
    np.ndarray
        The 2D order 0 image
    """
    # Get file
    # TODO: Add order 0 files for other modes
    filename = 'NIS_order0.npy' if 'NIS' in aperture else 'NIS_order0.npy'

    # Get the path to the trace files
    trace_path = os.path.join(os.environ['EXOCTK_DATA'], f'exoctk_contam/order0/{filename}')

    # Make frame
    trace = np.load(trace_path)

    return trace


def get_trace(aperture, teff, stype, verbose=False, plot=False):
    """Get the trace for the given aperture at the given temperature

    Parameters
    ----------
    aperture: str
        The aperture to use
    teff: int
        The temperature [K]
    stype: str
        The source type, ['STAR', 'GALAXY']
    verbose: bool
        Print statements
    plot: bool
        Plot the trace

    Returns
    -------
    np.ndarray
        The 2D trace
    """
    if 'DHS_F322W2' in aperture:
        aperpath = 'NRCA5_GRISM256_F322W2'
    elif 'DHS_F444W' in aperture:
        aperpath = 'NRCA5_GRISM256_F444W'
    elif 'NIS' in aperture:
        aperpath = 'NIS_SUBSTRIP256'
    else:
        aperpath = aperture

    # Get the path to the trace files
    traces_path = os.path.join(os.environ['EXOCTK_DATA'], f'exoctk_contam/traces/{aperpath}/*.fits')

    # Glob the file names
    trace_files = glob.glob(traces_path)

    # Get closest Teff
    teffs = np.array([int(os.path.basename(file).split('_')[-1][:-5]) for file in trace_files])
    file = trace_files[np.argmin((teffs - teff)**2)]
    if verbose:
        print('Fetching {} {}K trace from {}'.format(aperture, teff, file))

    # Get data
    if 'NIS' in aperture:
        # Orders stored separately just in case ;)
        traceo1 = fits.getdata(file, ext=0)
        traceo2 = fits.getdata(file, ext=1)
        traceo3 = fits.getdata(file, ext=2)

        # Zero out background
        traceo1[traceo1 < 1] = 0
        traceo2[traceo2 < 1] = 0
        traceo3[traceo3 < 1] = 0

        # 1.5 scaling based on comparison with observational data
        obs_scale = 1.5
        traces = [traceo1 * obs_scale, traceo2 * obs_scale, traceo3 * obs_scale]

        if stype == 'GALAXY':

            # Just mask trace area
            traces = NIRISS_SOSS_trace_mask(aperture)

    elif 'NRCA5' in aperture:

        # Get the trace and replace the NaN values
        trace = fits.getdata(file)[0]
        trace = replace_NaNs(trace)

        # Put the trace in each of the DHS trace positions
        traces = []
        for stripe in DHS_STRIPES[aperture].values():
            y, x = int((stripe['y1'] + stripe['y0']) / 2.), int((stripe['x1'] + stripe['x0']) / 2.)
            dhs_trace = add_array_at_position(np.zeros((4257, 4257)), trace, x, y, centered=True)
            traces.append(dhs_trace)

        if stype == 'GALAXY':
            traces = NIRCam_DHS_trace_mask(aperture, plot=False)

        traces = [trace[APERTURES[aperture]['subarr_y'][1]:APERTURES[aperture]['subarr_y'][2], :] for trace in traces]

    else:
        traces = [fits.getdata(file)]

    if plot:
        f = figure(width=900, height=450)
        final = np.sum(traces, axis=0)
        f.image([final], x=APERTURES[aperture]['subarr_x'][0], y=APERTURES[aperture]['subarr_y'][1], dw=final.shape[1], dh=final.shape[0])
        show(f)

    return traces


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

from copy import copy
from functools import partial
import glob
from multiprocessing import pool, cpu_count
import os
import re
import time
from pkg_resources import resource_filename

import astropy.coordinates as crd
from astropy.io import fits
from astroquery.irsa import Irsa
import astropy.units as u
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.irsa import Irsa
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
from astroquery.gaia import Gaia
from astropy.io import fits
from bokeh.plotting import figure, show
from bokeh.layouts import gridplot, column
from bokeh.models import Range1d, LinearColorMapper, LogColorMapper, Label, ColorBar, ColumnDataSource, HoverTool, Slider, CustomJS, VArea, CrosshairTool, TapTool, OpenURL, Span, Legend
from bokeh.palettes import PuBu, Spectral6
from bokeh.transform import linear_cmap
from hotsoss.plotting import plot_frame
from hotsoss.locate_trace import trace_polynomial
from scipy.ndimage.interpolation import rotate
import numpy as np
import pysiaf
import regions

from ..utils import get_env_variables, check_for_data
# from .visibilityPA import using_gtvt
from .new_vis_plot import build_visibility_plot, get_exoplanet_positions
from .contamination_figure import contam

Vizier.columns = ["**", "+_r"]
Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # DR2 is default catalog
Gaia.ROW_LIMIT = 100

from exoctk import utils

APERTURES = {'NIS_SOSSFULL': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.8, 2.8], 'subarr_x': [0, 2048, 2048, 0], 'subarr_y':[0, 0, 2048, 2048], 'trim': [127, 126, 252, 1]},
             'NIS_SUBSTRIP96': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.8, 2.8], 'subarr_x': [0, 2048, 2048, 0], 'subarr_y':[1792, 1792, 1888, 1888], 'trim': [47, 46, 0, 1]},
             'NIS_SUBSTRIP256': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.8, 2.8], 'subarr_x': [0, 2048, 2048, 0], 'subarr_y':[1792, 1792, 2048, 2048], 'trim': [127, 126, 0, 1]},
             'NRCA5_GRISM256_F277W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [2.395, 3.179], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F322W2': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [2.413, 4.083], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F356W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [3.100, 4.041], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F444W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [3.835, 5.084], 'trim': [0, 1, 1250, 1]},
             'MIRIM_SLITLESSPRISM': {'inst': 'MIRI', 'full': 'MIRIM_FULL', 'scale': 0.11, 'rad': 2.0, 'lam': [5, 12], 'trim': [6, 5, 0, 1]}}

# Gaia color-Teff relation
GAIA_TEFFS = np.asarray(np.genfromtxt(resource_filename('exoctk', 'data/contam_visibility/predicted_gaia_colour.txt'), unpack=True))

def SOSS_trace_mask(aperture, radius=20):
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
    # TODO: Implement order 3 trace in hotsoss.locate_trace.trace_polynomial then enable it here.
    traces = trace_polynomial(evaluate=True)
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


def trace_dict(teffs=None):
    """
    Load the trace data for all the given Teff values into a dictionary

    Parameters
    ----------
    teffs: sequence
        The teff values to fetch

    Returns
    -------
    dict
        The trace data for the given Teff values
    """
    teff_dict = {}

    if teffs is None:
        teffs = np.arange(2000, 12100, 100)

    # Make sure they're ints
    teffs = [int(teff) for teff in teffs]

    for teff in teffs:
        teff_dict[teff] = get_trace('NIS_SUBSTRIP256', teff, verbose=False)

    return teff_dict


def find_stars(ra, dec, width=5*u.arcmin, catalog='Gaia', verbose=False):
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

        # Derived from K. Volk
        # TODO: What to do for sources with no bp-rp color? Uses 2300K if missing.
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

        # Star names
        stars['name'] = [str(i) for i in stars['designation']]

        # Catalog name
        cat = 'II/246/out'

    # Find distance from target to each star
    sindRA = (stars['ra'][0] - stars['ra']) * np.cos(stars['dec'][0])
    cosdRA = stars['dec'][0] - stars['dec']
    stars.add_column(np.sqrt(sindRA ** 2 + cosdRA ** 2) * u.deg.to(u.arcsec), name='distance')
    stars.sort('distance')

    # Add detector location to the table
    stars.add_columns(np.zeros((10, len(stars))), names=['xtel', 'ytel', 'xdet', 'ydet', 'xsci', 'ysci', 'xord0', 'yord0', 'xord1', 'yord1'])

    # Add URL
    urls = ['https://vizier.u-strasbg.fr/viz-bin/VizieR-5?-ref=VIZ62fa613b20f3fc&-out.add=.&-source={}&-c={}%20%2b{},eq=ICRS,rs=2&-out.orig=o'.format(cat, ra_deg, dec_deg) for ra_deg, dec_deg in zip(stars['ra'], stars['dec'])]
    stars.add_column(urls, name='url')

    # Get traces
    traces = trace_dict(teffs=np.unique(stars['Teff']))
    stars.add_column([np.zeros((256, 2048))] * len(stars), name='trace_o1')
    stars.add_column([np.zeros((256, 2048))] * len(stars), name='trace_o2')
    stars.add_column([np.zeros((256, 2048))] * len(stars), name='trace_o3')
    for star in stars:
        star['trace_o1'], star['trace_o2'], star['trace_o3'] = traces[int(star['Teff'])]

    return stars


def add_star(startable, name, ra, dec, teff, fluxscale=None, delta_mag=None, dist=None, pa=None):
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

    # Get the trace
    trace = get_trace('NIS_SUBSTRIP256', teff, verbose=False)

    # Add the row to the table
    startable.add_row({'name': name, 'ra': ra, 'dec': dec, 'Teff': teff, 'fluxscale': fluxscale, 'distance': dist, 'trace_o1': trace[0], 'trace_o2': trace[1], 'trace_o3': trace[2]})
    startable.sort('distance')

    return startable


def calc_v3pa(V3PA, stars, aperture, data=None, c0x0=885, c0y0=1462, c1x0=-0.11, c1y0=0, c1y1=0, c1x1=0.02, tilt=0, ord0scale=1, ord1scale=1, plot=False, verbose=False):
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
    ref: str
        The reference frame to plot in, ['tel', 'det', 'sci']
    floor: float
        The noise floor to zero out
    plot: bool
        Plot the full frame and subarray bounds with all traces
    verbose: bool
        Print statements

    Returns
    -------
    targframe, starframe
        The frame containing the target trace and a frame containing all contaminating star traces
    """
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

    # Order 0 location relative to pysiaf SCI coordinates
    x_sweet = 2865
    y_sweet = 1720
    stars['xord0'][0] = int(stars['xsci'][0] + c0x0 + c1x0 * (stars['ysci'][0] + c0y0 - y_sweet))
    stars['yord0'][0] = int(stars['ysci'][0] + c0y0)
    stars['xord1'][0] = stars['xord0'][0] - x_sweet + aper['subarr_x'][0]
    stars['yord1'][0] = stars['yord0'][0] - y_sweet + aper['subarr_y'][1]

    # Get target's attitude matrix for each Position Angle
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

        # Order 0 location relative to pysiaf SCI coordinates
        star['xord0'] = int(star['xsci'] + c0x0 + c1x0 * (star['ysci'] + c0y0 - y_sweet))
        star['yord0'] = int(star['ysci'] + c0y0)

        # Order 1/2/3 location relative to order 0 location
        x_shift = int(c1x0 + c1x1 * (stars[0]['xord0'] - star['xord0']))
        y_shift = int(c1y0 + c1y1 * (stars[0]['yord0'] - star['yord0']) - c1x1 * (stars[0]['xord0'] - star['xord0']))
        star['xord1'] = star['xord0'] - x_sweet + aper['subarr_x'][0] + x_shift
        star['yord1'] = star['yord0'] - y_sweet + aper['subarr_y'][1] + y_shift

    # Just stars in FOV (Should always have at least 1, the target)
    lft, rgt, top, bot = 700, 5100, 1940, 1400
    FOVstars = stars[(lft < stars['xord0']) & (stars['xord0'] < rgt) & (bot < stars['yord0']) & (stars['yord0'] < top)]

    if verbose:
        print("Calculating contamination from {} other stars in the FOV".format(len(FOVstars) - 1))

    # Make frame for the target and a frame for all the other stars
    targframe_o1 = np.zeros((subY, subX))
    targframe_o2 = np.zeros((subY, subX))
    targframe_o3 = np.zeros((subY, subX))
    starframe = np.zeros((subY, subX))

    if plot:
        # Set up hover tool
        tips = [('Name', '@name'), ('RA', '@ra'), ('DEC', '@dec'), ('scale', '@fluxscale'), ('Teff', '@Teff'), ('ord0', '@xord0{int}, @yord0{int}')]
        hover = HoverTool(tooltips=tips, name='stars')
        crosshair = CrosshairTool(dimensions="height")
        taptool = TapTool(behavior='select', callback=OpenURL(url="@url"))

        # Make the plot
        tools = ['pan', crosshair, 'reset', 'box_zoom', 'wheel_zoom', 'save', taptool, hover]
        fig = figure(title='Generated FOV from Gaia EDR3', width=900, height=subY, match_aspect=True, tools=tools)
        fig.title = '({}, {}) at PA={} in {}'.format(stars[0]['ra'], stars[0]['dec'], V3PA, aperture.AperName)

        # Plot config
        scale = 'log'
        color_map = 'Viridis256'

        # Plot the obs data if possible
        if data is not None:
            vmax = np.nanmax(data)
            if scale == 'log':
                mapper = LogColorMapper(palette=color_map, low=1, high=vmax)
            else:
                mapper = LinearColorMapper(palette=color_map, low=0, high=vmax)
            data[data < 0] = 0
            data = rotate(data, tilt)
            fig.image([data], x=0, y=2048 - data.shape[0], dh=data.shape[0], dw=2048, color_mapper=mapper)

    # Get order 0
    order0 = get_order0(aperture.AperName) * 1.5e8 # Scaling factor based on observations

    # SOSS trace masks
    mask1, mask2, mask3 = SOSS_trace_mask(aperture.AperName)

    # Iterate over all stars in the FOV and add their scaled traces to the correct frame
    for idx, star in enumerate(FOVstars):

        # Scale the order 0 image and get dims
        scale0 = copy(order0) * star['fluxscale'] * ord0scale
        dim0y, dim0x = scale0.shape
        dim0y0 = int(dim0y / 2)
        dim0y1 = dim0y - dim0y0
        dim0x0 = int(dim0x / 2)
        dim0x1 = dim0x - dim0x0

        # Locations of the order 0 pixels on the subarray
        f0x0, f1x0 = int(max(aper['subarr_x'][0], star['xord0'] - dim0x0)), int(min(aper['subarr_x'][1], star['xord0'] + dim0x1))
        f0y0, f1y0 = int(max(aper['subarr_y'][1], star['yord0'] - dim0y0)), int(min(aper['subarr_y'][2], star['yord0'] + dim0y1))

        if 0 < f1x0 - f0x0 <= dim0x and 0 < f1y0 - f0y0 <= dim0y:

            # How many pixels of the order 0 image fall on the subarray
            t0x0 = dim0x - (f1x0 - f0x0) if f0x0 == aper['subarr_x'][0] else 0
            t1x0 = f1x0 - f0x0 if f1x0 == aper['subarr_x'][1] else dim0x
            t0y0 = dim0y - (f1y0 - f0y0) if f0y0 == aper['subarr_y'][0] else 0
            t1y0 = f1y0 - f0y0 if f1y0 == aper['subarr_y'][2] else dim0y

            # if verbose:
            #     print("{} x {} pixels of star {} order 0 fall on {}".format(t1y0 - t0y0, t1x0 - t0x0, idx, aperture.AperName))

            # Target order 0 is never on the subarray so add all order 0s to the starframe
            starframe[f0y0 - aper['subarr_y'][1]:f1y0 - aper['subarr_y'][1], f0x0 - aper['subarr_x'][1]:f1x0 - aper['subarr_x'][0]] += scale0[t0y0:t1y0, t0x0:t1x0]

        # Higher Orders ============================================================================

        # Get the appropriate trace
        trace_o1 = star['trace_o1']
        trace_o2 = star['trace_o2']
        trace_o3 = star['trace_o3']

        # Orient trace if need be
        if 'NIS' in aperture.AperName:
            trace_o1 = np.rot90(trace_o1.T[:, ::-1] * 1.5, k=1) # Scaling factor based on observations
            trace_o2 = np.rot90(trace_o2.T[:, ::-1] * 1.5, k=1) # Scaling factor based on observations
            trace_o3 = np.rot90(trace_o3.T[:, ::-1] * 1.5, k=1) # Scaling factor based on observations

            # Pad or trim SUBSTRIP256 simulation for SUBSTRIP96 or FULL frame
            if aperture.AperName == 'NIS_SOSSFULL':
                trace_o1 = np.pad(trace_o1, ((1792, 0), (0, 0)), 'constant')
                trace_o2 = np.pad(trace_o2, ((1792, 0), (0, 0)), 'constant')
                trace_o3 = np.pad(trace_o3, ((1792, 0), (0, 0)), 'constant')
            elif aperture.AperName == 'NIS_SUBSTRIP96':
                trace_o1 = trace_o1[:96, :]
                trace_o2 = trace_o2[:96, :]
                trace_o3 = trace_o3[:96, :]

        # Get the trace and shift into the correct subarray position
        # trace *= mask1 + mask2 + mask3

        # Scale the order 1, 2, 3 image and get dims
        trace_o1 *= star['fluxscale'] * ord1scale
        if aperture.AperName.startswith('NIS'):
            trace_o2 *= star['fluxscale'] * ord1scale
            trace_o3 *= star['fluxscale'] * ord1scale
        dimy, dimx = trace_o1.shape

        # Func to remove background
        # def bkg(trc, cutoff=0.00001):
        #     trc[trc < 0] = cutoff
        #     trc[np.isnan(trc)] = cutoff
        #     trc[np.isinf(trc)] = cutoff
        #     return trc

        # Trim background
        # trace_o1 = bkg(trace_o1)
        # if aperture.AperName.startswith('NIS'):
        #     trace_o2 = bkg(trace_o2)
        #     trace_o3 = bkg(trace_o3)

        # Location of full trace footprint
        fpx0 = int(star['xord1'])
        fpx1 = int(fpx0 + dimx)
        fpy0 = int(star['yord1'])
        fpy1 = int(fpy0 + dimy)

        # Locations of the trace pixels on the subarray
        f0x, f1x = max(aper['subarr_x'][0], fpx0), min(aper['subarr_x'][1], fpx1)
        f0y, f1y = max(aper['subarr_y'][1], fpy0), min(aper['subarr_y'][2], fpy1)

        # print(idx, f0x, f1x, f0y, f1y)
        if 0 < f1x - f0x <= dimx and 0 < f1y - f0y <= dimy:

            # How many pixels of the trace image fall on the subarray
            t0x = dimx - (f1x - f0x) if f0x == aper['subarr_x'][0] else 0
            t1x = f1x - f0x if f1x == aper['subarr_x'][1] else dimx
            t0y = dimy - (f1y - f0y) if f0y == aper['subarr_y'][0] else 0
            t1y = f1y - f0y if f1y == aper['subarr_y'][2] else dimy

            # if verbose:
            #     print("{} x {} pixels of star {} trace fall on {}".format(t1y - t0y, t1x - t0x, idx, aperture.AperName))

            # Box to show footprint of full trace
            # fig.patch([fpx0, fpx1, fpx1, fpx0], [fpy0, fpy0, fpy1, fpy1], line_color="black", fill_color='black', fill_alpha=0.1)

            # Add each order to it's own frame
            if idx == 0:

                targframe_o1[f0y - aper['subarr_y'][1]:f1y - aper['subarr_y'][1], f0x - aper['subarr_x'][1]:f1x - aper['subarr_x'][0]] += trace_o1[t0y:t1y, t0x:t1x]
                if aperture.AperName.startswith('NIS'):
                    targframe_o2[f0y - aper['subarr_y'][1]:f1y - aper['subarr_y'][1], f0x - aper['subarr_x'][1]:f1x - aper['subarr_x'][0]] += trace_o2[t0y:t1y, t0x:t1x]
                    targframe_o3[f0y - aper['subarr_y'][1]:f1y - aper['subarr_y'][1], f0x - aper['subarr_x'][1]:f1x - aper['subarr_x'][0]] += trace_o3[t0y:t1y, t0x:t1x]

            # Add all orders to the same frame
            else:

                starframe[f0y - aper['subarr_y'][1]:f1y - aper['subarr_y'][1], f0x - aper['subarr_x'][1]:f1x - aper['subarr_x'][0]] += trace_o1[t0y:t1y, t0x:t1x]
                if aperture.AperName.startswith('NIS'):
                    starframe[f0y - aper['subarr_y'][1]:f1y - aper['subarr_y'][1], f0x - aper['subarr_x'][1]:f1x - aper['subarr_x'][0]] += trace_o2[t0y:t1y, t0x:t1x]
                    starframe[f0y - aper['subarr_y'][1]:f1y - aper['subarr_y'][1], f0x - aper['subarr_x'][1]:f1x - aper['subarr_x'][0]] += trace_o3[t0y:t1y, t0x:t1x]

    # Contam per order
    simframe_o1 = targframe_o1 + starframe
    simframe_o2 = targframe_o2 + starframe
    simframe_o3 = targframe_o3 + starframe
    simframe = targframe_o1 + targframe_o2 + targframe_o3 + starframe

    # Calculate contam/total counts in each detector column
    pctframe_o1 = starframe / simframe_o1
    pctframe_o2 = starframe / simframe_o2
    pctframe_o3 = starframe / simframe_o3
    pctline_o1 = np.nanmean(pctframe_o1 * mask1, axis=0)
    pctline_o2 = np.nanmean(pctframe_o2 * mask2, axis=0)
    pctline_o3 = np.nanmean(pctframe_o3 * mask3, axis=0)

    result = {'pa': V3PA, 'target': targframe_o1 + targframe_o2 + targframe_o3, 'target_o1': targframe_o1, 'target_o2': targframe_o2, 'target_o3': targframe_o3,  'contaminants': starframe, 'sources': FOVstars, 'order1_contam': pctline_o1, 'order2_contam': pctline_o2, 'order3_contam': pctline_o3}

    if plot:

        # Plot the simulated frame
        vmax = np.nanmax(simframe)
        if scale == 'log':
            mapper = LogColorMapper(palette=color_map, low=1, high=vmax)
        else:
            mapper = LinearColorMapper(palette=color_map, low=0, high=vmax)

        # Only plot the simulation if no data is available to plot
        if data is None:
            fig.image(image=[simframe], x=aper['subarr_x'][0], dw=subX, y=aper['subarr_y'][1], dh=subY, color_mapper=mapper)

        mapper = linear_cmap(field_name='Teff', palette=Spectral6, low=np.nanmin(FOVstars['Teff']), high=np.nanmax(FOVstars['Teff']))

        # Plot order 0 locations
        fig.circle('xord0', 'yord0', color=mapper, size=15, line_width=3, fill_color=None, name='stars', source=dict(FOVstars[['Teff', 'xord0', 'yord0', 'ra', 'dec', 'name', 'url', 'fluxscale', 'xdet', 'ydet', 'xtel', 'ytel']]))
        fig.circle([x_sweet], [y_sweet], size=10, line_width=3, fill_color=None, line_color='black')

        fig = plot_traces(FOVstars, fig)

        # Show the figure
        fig.x_range = Range1d(aper['subarr_x'][0], aper['subarr_x'][1])
        fig.y_range = Range1d(aper['subarr_y'][1], aper['subarr_y'][2])

        # Source for ratio plot
        rsource = ColumnDataSource(data=dict(x=np.arange(subX), zeros=np.zeros(subX), o1=pctline_o1, o2=pctline_o2, o3=pctline_o3))

        # Make plot
        rfig = figure(title='Target Contamination', width=900, height=200, match_aspect=True, tools=tools, x_range=fig.x_range)
        rfig.line(np.arange(subX), pctline_o1, color='blue', legend_label='Order 1')
        glyph1 = VArea(x='x', y1='zeros', y2='o1', fill_color="blue", fill_alpha=0.3)
        rfig.add_glyph(rsource, glyph1)
        if aperture.AperName not in ['NIS_SUBSTRIP96']:
            rfig.line(np.arange(subX), pctline_o2, color='red', legend_label='Order 2')
            glyph2 = VArea(x='x', y1='zeros', y2='o2', fill_color="red", fill_alpha=0.3)
            rfig.add_glyph(rsource, glyph2)
            rfig.line(np.arange(subX), pctline_o3, color='green', legend_label='Order 3')
            glyph3 = VArea(x='x', y1='zeros', y2='o3', fill_color="green", fill_alpha=0.3)
            rfig.add_glyph(rsource, glyph3)
        rfig.y_range = Range1d(0, min(1, max(pctline_o1.max(), pctline_o2.max(), pctline_o3.max())))
        rfig.yaxis.axis_label = 'Contam / Total Counts'
        rfig.xaxis.axis_label = 'Detector Column'

        # Color bar
        # color_bar = ColorBar(color_mapper=mapper['transform'], width=10, location=(0, 0), title="Teff")
        # fig.add_layout(color_bar, 'right')

        # Plot grid
        gp = gridplot([[fig], [rfig]])

        return result, gp

    return result


def plot_traces(star_table, fig, color='red'):
    """
    PLot the trace locations of all the stars in the table

    Parameters
    ----------
    star_table: astropy.table.Table
        The table of stars
    fig: bokeh.plotting.figure.Figure
        The figure to plot on

    Returns
    -------
    fig
        The figure
    """

    # Trace extends in dispersion direction further than 2048 subarray edges
    blue_ext = 150
    red_ext = 200

    # Get the new x-ranges
    xr0 = np.linspace(-blue_ext, 2048 + red_ext, 1000)
    xr1 = np.linspace(-blue_ext, 1820 + red_ext, 1000)
    xr2 = np.linspace(-blue_ext, 1130 + red_ext, 1000)

    # Add the y-intercept to the c0 coefficient
    polys = trace_polynomial()
    yr0 = np.polyval(polys[0], xr0)
    yr1 = np.polyval(polys[1], xr1)
    yr2 = np.polyval(polys[2], xr2)

    for idx, star in enumerate(star_table):
        # Order 1/2/3 location relative to order 0
        fig.line(xr0 + star['xord1'], yr0 + star['yord1'], color=color, line_dash='solid' if idx == 0 else 'dashed')
        fig.line(xr1 + star['xord1'], yr1 + star['yord1'], color=color, line_dash='solid' if idx == 0 else 'dashed')
        fig.line(xr2 + star['xord1'], yr2 + star['yord1'], color=color, line_dash='solid' if idx == 0 else 'dashed')

    return fig


def field_simulation(ra, dec, aperture, binComp=None, n_jobs=-1, plot=False, multi=True, verbose=True):
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
    n_jobs: int
        Number of cores to use (-1 = All)

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
    stars = find_stars(ra, dec, verbose=verbose)

    # Add stars manually
    if isinstance(binComp, dict):
        stars = add_star(stars, **binComp)

    # Set the number of cores for multiprocessing
    max_cores = 8
    if n_jobs == -1 or n_jobs > max_cores:
        n_jobs = max_cores

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

    # Flatten list and check against 360 angles to get all bad PAs
    # badPA_list = [pa for pa in pa_list if pa not in goodPA_list]

    # Time it
    if verbose:
        print('Calculating target contamination from {} neighboring sources in position angle ranges {}...'.format(len(stars), good_group_bounds))
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
        func = partial(calc_v3pa, stars=stars, aperture=aper, plot=False, verbose=False)
        results = pl.map(func, goodPA_list)
        pl.close()
        pl.join()

    else:
        results = []
        for pa in goodPA_list:
            result = calc_v3pa(pa, stars=stars, aperture=aper, plot=False, verbose=False)
            results.append(result)

    # We only need one target frame frames
    targframe_o1 = np.asarray(results[0]['target_o1'])
    targframe_o2 = np.asarray(results[0]['target_o2'])
    targframe_o3 = np.asarray(results[0]['target_o3'])
    targframe = [targframe_o1, targframe_o2, targframe_o3]

    # Make sure starcube is of shape (PA, rows, cols)
    starcube = np.zeros((360, targframe_o1.shape[0], targframe_o1.shape[1]))

    # Make the contamination plot
    for result in results:
        starcube[result['pa'], :, :] = result['contaminants']

    if verbose:
        print('Contamination calculation complete: {} {}'.format(round(time.time() - start, 3), 's'))

    # Make contam plot
    if plot:
        contam_slider_plot(results, plot=plot)

    return targframe, starcube, results


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
    plt = figure(width=900, height=300, tools=['reset', 'box_zoom', 'wheel_zoom', 'save'])
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
                    step=int((max(pa_list) - min(pa_list)) / (len(pa_list) - 1)))

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
    viz_plt = figure(width=900, height=200, x_range=Range1d(0, 359))
    viz_plt.step(np.arange(360), np.mean(order1_contam, axis=1), color='blue', mode="center")
    viz_plt.step(np.arange(360), np.mean(order2_contam, axis=1), color='red', mode="center")
    viz_plt.step(np.arange(360), np.mean(order3_contam, axis=1), color='green', mode="center")
    c1 = viz_plt.vbar(x=np.arange(360), top=viz_ord1, line_color=None, fill_color='blue', alpha=0.2)
    c2 = viz_plt.vbar(x=np.arange(360), top=viz_ord2, line_color=None, fill_color='red', alpha=0.2)
    c3 = viz_plt.vbar(x=np.arange(360), top=viz_ord3, line_color=None, fill_color='green', alpha=0.2)
    c0 = viz_plt.vbar(x=np.arange(360), top=viz_none, color='black', alpha=0.6)
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
    if 'NIS' in aperture:
        filename = 'NIS_order0.npy'

    # Get the path to the trace files
    trace_path = os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/order0/{}'.format(filename))

    # Make frame
    trace = np.load(trace_path)

    return trace


def get_trace(aperture, teff, verbose=False):
    """Get the trace for the given aperture at the given temperature

    Parameters
    ----------
    aperture: str
        The aperture to use
    teff: int
        The temperature [K]

    Returns
    -------
    np.ndarray
        The 2D trace
    """
    # Get the path to the trace files
    traces_path = os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/traces/{}/*.fits'.format('NIS_SUBSTRIP256' if 'NIS' in aperture else aperture))

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

        trace = [traceo1, traceo2, traceo3]

    else:
        trace = [fits.getdata(file)]

    # # Expand to SUBSTRIP256 to FULL frame for NIS_SOSSFULL
    # if aperture == 'NIS_SOSSFULL':
    #     full_trace = np.zeros((2301, 2301))
    #     full_trace[:, -257:] = trace
    #     trace = full_trace

    return trace


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

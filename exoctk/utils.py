# !/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for utility funtions
"""
import itertools
import os
import re
import requests
import urllib

from astropy.io import fits
import bokeh.palettes as bpal
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import numpy as np
from svo_filters import svo

EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
MODELGRID_DIR = os.path.join(EXOCTK_DATA, 'modelgrid/')
FORTGRID_DIR = os.path.join(EXOCTK_DATA, 'fortney/')
EXOCTKLOG_DIR = os.path.join(EXOCTK_DATA, 'exoctk_log/')
GENERICGRID_DIR = os.path.join(EXOCTK_DATA, 'generic/')

# Supported profiles
PROFILES = ['uniform', 'linear', 'quadratic',
            'square-root', 'logarithmic', 'exponential',
            '3-parameter', '4-parameter']

# Supported filters
FILTERS = svo.filters()

# Set the version
VERSION = '0.2'


def color_gen(colormap='viridis', key=None, n=10):
    """Color generator for Bokeh plots

    Parameters
    ----------
    colormap: str, sequence
        The name of the color map

    Returns
    -------
    generator
        A generator for the color palette
    """
    if colormap in dir(bpal):
        palette = getattr(bpal, colormap)

        if isinstance(palette, dict):
            if key is None:
                key = list(palette.keys())[0]
            palette = palette[key]

        elif callable(palette):
            palette = palette(n)

        else:
            raise TypeError("pallette must be a bokeh palette name or a sequence of color hex values.")

    elif isinstance(colormap, (list, tuple)):
        palette = colormap

    else:
        raise TypeError("pallette must be a bokeh palette name or a sequence of color hex values.")

    yield from itertools.cycle(palette)


COLORS = color_gen('Category10')


def interp_flux(mu, flux, params, values):
    """
    Interpolate a cube of synthetic spectra for a
    given index of mu

    Parameters
    ----------
    mu: int
        The index of the (Teff, logg, FeH, *mu*, wavelength)
        data cube to interpolate
    flux: np.ndarray
        The 5D data array
    params: list
        A list of each free parameter range
    values: list
        A list of each free parameter values

    Returns
    -------
    tu
        The array of new flux values
    """
    # Iterate over each wavelength (-1 index of flux array)
    shp = flux.shape[-1]
    flx = np.zeros(shp)
    generators = []
    for lam in range(shp):
        interp_f = RegularGridInterpolator(params, flux[:, :, :, mu, lam])
        f, = interp_f(values)

        flx[lam] = f
        generators.append(interp_f)

    return flx, generators


def calc_zoom(R_f, arr):
    """
    Calculate the zoom factor required to make the given
    array into the given resolution

    Parameters
    ----------
    R_f: int
        The desired final resolution of the wavelength array
    arr: array-like
        The array to zoom
    """
    # Get initial resolution
    lam = arr[-1]-arr[0]
    d_lam_i = np.nanmean(np.diff(arr))
    # R_i = lam/d_lam_i

    # Calculate zoom
    d_lam_f = lam/R_f
    z = d_lam_i/d_lam_f

    return z


def rebin_spec(spec, wavnew, oversamp=100, plot=False):
    """
    Rebin a spectrum to a new wavelength array while preserving
    the total flux

    Parameters
    ----------
    spec: array-like
        The wavelength and flux to be binned
    wavenew: array-like
        The new wavelength array

    Returns
    -------
    np.ndarray
        The rebinned flux

    """
    wave, flux = spec
    nlam = len(wave)
    x0 = np.arange(nlam, dtype=float)
    x0int = np.arange((nlam-1.) * oversamp + 1., dtype=float)/oversamp
    w0int = np.interp(x0int, x0, wave)
    spec0int = np.interp(w0int, wave, flux)/oversamp

    # Set up the bin edges for down-binning
    maxdiffw1 = np.diff(wavnew).max()
    w1bins = np.concatenate(([wavnew[0]-maxdiffw1],
                             .5*(wavnew[1::]+wavnew[0: -1]),
                             [wavnew[-1]+maxdiffw1]))

    # Bin down the interpolated spectrum:
    w1bins = np.sort(w1bins)
    nbins = len(w1bins)-1
    specnew = np.zeros(nbins)
    inds2 = [[w0int.searchsorted(w1bins[ii], side='left'),
              w0int.searchsorted(w1bins[ii+1], side='left')]
             for ii in range(nbins)]

    for ii in range(nbins):
        specnew[ii] = np.sum(spec0int[inds2[ii][0]: inds2[ii][1]])

    return specnew


def writeFITS(filename, extensions, headers=()):
    """
    Write some data to a new FITS file

    Parameters
    ----------
    filename: str
        The filename of the output FITS file
    extensions: dict
        The extension name and associated data to include
        in the file
    headers: array-like
        The (keyword, value, comment) groups for the PRIMARY
        header extension

    """
    # Write the arrays to a FITS file
    prihdu = fits.PrimaryHDU()
    prihdu.name = 'PRIMARY'
    hdulist = fits.HDUList([prihdu])

    # Write the header to the PRIMARY HDU
    hdulist['PRIMARY'].header.extend(headers, end=True)

    # Write the data to the HDU
    for k, v in extensions.items():
        hdulist.append(fits.ImageHDU(data=v, name=k))

    # Write the file
    hdulist.writeto(filename, clobber=True)
    hdulist.close()

    # Insert END card to prevent header error
    # hdulist[0].header.tofile(filename, endcard=True, clobber=True)


def smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    Parameters
    ----------
    x: sequence
        The input signal
    window_len: int
        The dimension of the smoothing window
    window: str
        The type of window from 'flat', 'hanning', 'hamming', 'bartlett',
        'blackman'. 'flat' window will produce a moving average smoothing.

    Retruns
    -------
    np.ndarray
        The smoothed signal

    Example
    -------
    t = linspace(-2, 2, 0.1)
    x = sin(t)+randn(len(t))*0.1
    y = smooth(x)
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming',\
                          'bartlett', 'blackman'")

    s = np.r_[2*np.median(x[0: window_len/5])-x[window_len: 1: -1], x,
              2*np.median(x[-window_len/5:])-x[-1: -window_len: -1]]

    if window == 'flat':
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='same')

    return y[window_len-1: -window_len+1]


def medfilt(x, window_len):
    """
    Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.

    Parameters
    ----------
    x: np.array
        The 1D array to smooth
    window_len: int
        The size of the smoothing window

    Returns
    -------
    np.ndarray
        The smoothed 1D array
    """
    # assert x.ndim == 1, "Input must be one-dimensional."
    if window_len % 2 == 0:
        s1 = "Median filter length ("
        s2 = ") must be odd. Adding 1."
        print(s1+str(window_len)+s2)
        window_len += 1
    window_len = int(window_len)
    k2 = int((window_len - 1)//2)
    s = np.r_[2*np.median(x[0: int(window_len/5)])-x[window_len: 1: -1],
              x, 2*np.median(x[int(-window_len/5):])-x[-1: -window_len: -1]]
    y = np.zeros((len(s), window_len), dtype=s.dtype)

    y[:, k2] = s
    for i in range(k2):
        j = k2 - i
        y[j:, i] = s[:-j]
        y[: j, i] = s[0]
        y[: -j, -(i+1)] = s[j:]
        y[-j:, -(i+1)] = s[-1]
    return np.median(y[window_len-1: -window_len+1], axis=1)


def filter_table(table, **kwargs):
    """Retrieve the filtered rows

    Parameters
    ----------
    table: astropy.table.Table, pandas.DataFrame
        The table to filter
    param: str
        The parameter to filter by, e.g. 'Teff'
    value: str, float, int, sequence
        The criteria to filter by,
        which can be single valued like 1400
        or a range with operators [<,<=,>,>=],
        e.g. ('>1200','<=1400')

    Returns
    -------
    astropy.table.Table, pandas.DataFrame
        The filtered table
    """
    for param, value in kwargs.items():

        # Check it is a valid column
        if param not in table.colnames:
            raise KeyError("No column named {}".format(param))

        # Wildcard case
        if isinstance(value, (str, bytes)) and '*' in value:

            # Get column data
            data = list(map(str, table[param]))

            if not value.startswith('*'):
                value = '^'+value
            if not value.endswith('*'):
                value = value+'$'

            # Strip double quotes and decod
            value = value.replace("'", '')\
                         .replace('"', '')\
                         .replace('*', '(.*)')

            # Regex
            reg = re.compile(value, re.IGNORECASE)
            keep = list(filter(reg.findall, data))

            # Get indexes
            idx = np.where([i in keep for i in data])

            # Filter table
            table = table[idx]

        else:

            # Make single value string into conditions
            if isinstance(value, str):

                # Check for operator
                if any([value.startswith(o) for o in ['<', '>', '=']]):
                    value = [value]

                # Assume eqality if no operator
                else:
                    value = ['=='+value]

            # Turn numbers into strings
            if isinstance(value, (int, float, np.float16)):
                value = ["=={}".format(value)]

            # Iterate through multiple conditions
            for cond in value:

                # Equality
                if cond.startswith('='):
                    v = cond.replace('=', '')
                    if v.replace('.','',1).isdigit():
                        table = table[table[param] == eval(v)]
                    else:
                        table = table[table[param] == v]

                # Less than or equal
                elif cond.startswith('<='):
                    v = cond.replace('<=', '')
                    table = table[table[param] <= eval(v)]

                # Less than
                elif cond.startswith('<'):
                    v = cond.replace('<', '')
                    table = table[table[param] < eval(v)]

                # Greater than or equal
                elif cond.startswith('>='):
                    v = cond.replace('>=', '')
                    table = table[table[param] >= eval(v)]

                # Greater than
                elif cond.startswith('>'):
                    v = cond.replace('>', '')
                    table = table[table[param] > eval(v)]

                else:
                    raise ValueError("'{}' operator not valid.".format(cond))

    return table


def find_closest(axes, points, n=1, values=False):
    """
    Find the n-neighboring elements of a given value in an array

    Parameters
    ----------
    axes: list, np.array
        The array(s) to search
    points: array-like, float
        The point(s) to search for
    n: int
        The number of values to the left and right of the points
    Returns
    -------
    np.ndarray
        The n-values to the left and right of 'points' in 'axes'
    """
    results = []
    if not isinstance(axes, list):
        axes = [axes]
        points = [points]

    for i, (axis, point) in enumerate(zip(axes, points)):
        if point >= min(axis) and point <= max(axis):
            axis = np.asarray(axis)
            idx = np.clip(axis.searchsorted(point), 1, len(axis)-1)
            slc = slice(max(0, idx-n), min(idx+n, len(axis)))

            if values:
                result = axis[slc]
            else:

                result = np.arange(0, len(axis))[slc].astype(int)

            results.append(result)
        else:
            print('Point {} outside grid.'.format(point))
            return

    return results

def build_target_url(target_name):
    '''Build restful api url based on target name.

    Parameters
        ----------
        target_name : string
            The name of the target transit.

        Returns
        -------
        target_url : string
    '''
    # Encode the target name string.
    encode_target_name = urllib.parse.quote(target_name, encoding='utf-8')
    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/{}/properties/".format(encode_target_name)

    return target_url

def get_canonical_name(target_name):
    '''Get ExoMAST prefered name for exoplanet.

        Parameters
        ----------
        target_name : string
            The name of the target transit.

        Returns
        -------
        canonical_name : string
    '''

    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/identifiers/"

    # Create params dict for url parsing. Easier than trying to format yourself.
    params = {"name":target_name}

    r = requests.get(target_url, params=params)
    planetnames = r.json()
    canonical_name = planetnames['canonicalName']

    return canonical_name

def get_target_data(target_name):
    """
    Send request to exomast restful api for target information.

    Parameters
    ----------
    target_name : string
        The name of the target transit

    Returns
    -------
    target_data: json:
        json object with target data.
    """

    canonical_name = get_canonical_name(target_name)

    target_url = build_target_url(canonical_name)

    r = requests.get(target_url)

    if r.status_code == 200:
        target_data = r.json()
    else:
        print('Whoops, no data for this target!')

    # Some targets have multiple catalogs
    # nexsci is the first choice.
    if len(target_data) > 1:
        # Get catalog names from exomast and make then the keys of a dictionary
        # and the values are its position in the json object.
        catalog_dict = {data['catalog_name']: index for index, data in enumerate(target_data)}

        # Parse based on catalog accuracy.
        if 'nexsci' in list(catalog_dict.keys()):
            target_data = target_data[catalog_dict['nexsci']]
        elif 'exoplanets.org' in list(catalog_dict.keys()):
            target_data = target_data[catalog_dict['exoplanets.org']]
        else:
            target_data = target_data[0]
    else:
        target_data = target_data[0]

    # Strip spaces and non numeric or alphabetic characters and combine.
    url = 'https://exo.mast.stsci.edu/exomast_planet.html?planet={}'.format(re.sub(r'\W+', '', canonical_name))

    return target_data, url

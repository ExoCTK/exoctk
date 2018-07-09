# !/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for utility funtions
"""
from __future__ import print_function

from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import numpy as np


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

    if plot:
        plt.figure()
        plt.loglog(wave, flux, c='b')
        plt.loglog(wavnew, specnew, c='r')

    return specnew


def multiplot(rows, columns, ylabel='', xlabel='', sharey=True, sharex=True,
              fontsize=20, figsize=(15, 7), title='', **kwargs):
    """
    Creates subplots with given number or *rows* and *columns*.

    Parameters
    ----------
    rows: int
        The number of rows in the figure
    columns: int
        The number of columns in the figure
    ylabel: str, list
        The shared y-label or list of y-labels for each column
    xlabel: str, list
        The shared y-label or list of y-labels for each column
    sharey: bool
        Same y-axis limits
    sharex: bool
        Same x-axis limits
    fontsize: int
        The fontsize to use throughout the figure
    figsize: tuple, list
        The (x, y) dimenstions of the figure
    title: str
        The title of the figure
    Returns
    -------
    list
        A list of the figure and axes objects for the current figure

    Example
    -------
    >>> fig, (ax11, ax12, ax13), (ax21, ax22, ax23) = multiplot(2, 3)
    >>> ax11.plot(x, y, label='Row 1, Col 1 Plot')
    """
    # Initialize the plot
    fig, axes = plt.subplots(rows, columns, sharey=sharey, sharex=sharex,
                             figsize=figsize)
    plt.rc('text', usetex=True)
    plt.rc('font', size=fontsize)

    # Set the y-label(s)
    if ylabel:
        if isinstance(ylabel, str):
            fig.text(0.04, 0.54, ylabel, ha='center', va='center',
                     rotation='vertical', **kwargs)
        else:
            if columns > 1:
                for a, l in zip(axes, ylabel):
                    a[0].set_ylabel(l, fontsize=fontsize, labelpad=fontsize)
            else:
                for a, l in zip(axes, ylabel):
                    a.set_ylabel(l, fontsize=fontsize, labelpad=fontsize)

    # Set the x-label(s)
    if xlabel:
        if isinstance(xlabel, str):
            fig.text(0.54, 0.04, xlabel, ha='center', va='center', **kwargs)
        else:
            if rows > 1:
                for a, l in zip(axes, xlabel):
                    a[0].set_xlabel(l, fontsize=fontsize, labelpad=fontsize)
            else:
                for a, l in zip(axes, xlabel):
                    a.set_xlabel(l, fontsize=fontsize, labelpad=fontsize)

    # Plot formatting
    plt.suptitle(title)
    plt.subplots_adjust(right=0.96, top=0.93 if title else 0.96, bottom=0.15,
                        left=0.12, hspace=0, wspace=0)
    fig.canvas.draw()

    return [fig] + list(axes)


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

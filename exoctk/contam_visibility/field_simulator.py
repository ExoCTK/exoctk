import glob
import os
from multiprocessing import pool, cpu_count
import time
from functools import partial
import re

import astropy.coordinates as crd
import astropy.units as u
from astropy.table import Table
from bokeh.plotting import show
from bokeh.layouts import gridplot
from bokeh.plotting import figure
from bokeh.models import Range1d, LinearColorMapper, Label
from bokeh.models.widgets import Panel, Tabs
from bokeh.palettes import PuBu
import numpy as np
from astroquery.irsa import Irsa
from matplotlib import cm
from scipy.io import readsav
from astropy.io import fits
import pysiaf
from pysiaf.utils import rotations

from ..utils import get_env_variables, check_for_data
from .visibilityPA import using_gtvt


APERTURES = {'NIS_SUBSTRIP96': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.6, 1.4]},
             'NIS_SUBSTRIP256': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.8, 2.8]},
             'NRCA5_GRISM256_F444W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [3.063, 5.111]},
             'NRCA5_GRISM256_F322W2': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [2.369, 4.417]},
             'MIRIM_SLITLESSPRISM': {'inst': 'MIRI', 'full': 'MIRIM_FULL', 'scale': 0.11, 'rad': 2.0, 'lam': [5, 12]}}


def calc_v3pa(V3PA, stars, aper, targetIndex, pixel_pos):
    """
    Calculate the V3 position angle for each target at the given PA

    Parameters
    ----------
    V3PA: float
        The PA in V3
    stars: dict
        The stars in the FOV
    aper: pysiaf.aperture.JwstAperture
        The aperture object for the given mode
    targetIndex: int
        The index of the target in the stars dictionary
    pixel_pos: sequence
        The positionx of the detectro corners

    Returns
    -------

    """
    # Calculate corresponding V2/V3 (TEL) coordinates for Sweetspot
    xSweet, ySweet = aper.reference_point('det')
    v2targ, v3targ = aper.det_to_tel(xSweet, ySweet)
    v3angle = aper.V3IdlYAngle

    # Get APA from V3PA
    APA = V3PA + v3angle
    if APA > 360:
        APA = APA - 360
    elif APA < 0:
        APA = APA + 360

    # Get subarray dims
    subX, subY = aper.XSciSize, aper.YSciSize

    # Get target's attitude matrix for each Position Angle
    attitude = rotations.attitude_matrix(v2targ, v3targ, stars['RA'][targetIndex], stars['DEC'][targetIndex], APA)

    xdet, ydet = [], []
    xsci, ysci = [], []
    for starRA, starDEC in zip(stars['RA'], stars['DEC']):

        # Get the TEL coordinates of each star w attitude matrix
        V2, V3 = rotations.sky_to_tel(attitude, starRA, starDEC)

        # Convert to arcsec and turn into a float
        V2, V3 = V2.to(u.arcsec).value, V3.to(u.arcsec).value

        # Save detector coords
        XDET, YDET = aper.tel_to_det(V2, V3)
        XSCI, YSCI = aper.det_to_sci(XDET, YDET)
        xdet.append(XDET)
        ydet.append(YDET)
        xsci.append(XSCI)
        ysci.append(YSCI)

    # Record keeping
    stars['xdet'], stars['ydet'] = np.array(xdet), np.array(ydet)
    stars['xsci'], stars['ysci'] = np.array(xsci), np.array(ysci)

    # Target position
    sci_targx = stars['xsci'][targetIndex]
    sci_targy = stars['ysci'][targetIndex]

    # Determine stars in FOV
    inFOV = []
    for star in range(stars['RA'].size):

        x, y = stars['xdet'][star], stars['ydet'][star]
        if (pixel_pos[0] < y) & (y < pixel_pos[1]) & (pixel_pos[2] < x) & (x < pixel_pos[3]):
            inFOV.append(star)
    inFOV = np.array(inFOV)

    # Make simulated image
    simuCube = np.zeros((361, subY, subX))
    for idx in inFOV:

        # Get the Teff
        sci_dx = round(sci_targx - stars['xsci'][idx])
        sci_dy = round(sci_targy - stars['ysci'][idx])
        temp = stars['Temp'][idx]

        # Get the flux scaling
        fluxscale = 10.0**(-0.4 * (stars['Jmag'][idx] - stars['Jmag'][targetIndex]))

        # Get the trace and pad it
        trace = get_trace(aper.AperName, int(temp))
        # trace = np.rot90(trace, k=3)
        pad_trace = np.pad(trace, pad_width=5000, mode='constant', constant_values=0)

        # Determine the highest pixel value of trace
        maxY, maxX = np.where(pad_trace == pad_trace.max())
        peakY, peakX = maxY[0], maxX[0]

        # Use relative distances (sci_dx, sci_dy) to find target
        xTarg = peakX + sci_dx
        yTarg = peakY + sci_dy

        # Use the (xTarg, yTarg) coordinates to slice out subarray
        # remember X is columns, Y is rows
        dimX0, dimX1 = xTarg - sci_targx, xTarg + subX - sci_targx
        dimY0, dimY1 = yTarg - sci_targy, yTarg + subY - sci_targy

        if dimX0 < 0:
            dimX0 = 0
            dimX1 = subX
        if dimY0 < 0:
            dimY0 = 0
            dimY1 = subY

        traceY, traceX = pad_trace.shape
        if dimX1 > traceX:
            dimX1 = traceX
            dimX0 = traceX - subX
        if dimY1 > traceY:
            dimY1 = traceY
            dimY0 = traceY - subY

        if (dimX1 < 0) or (dimY1 < 0):
            continue

        # -1 because pySIAF is 1-indexed
        mx0, mx1 = int(dimX0) - 1, int(dimX1) - 1
        my0, my1 = int(dimY0) - 1, int(dimY1) - 1

        # Add trace to the image
        tr = pad_trace[my0:my1, mx0:mx1] * fluxscale
        trY, trX = tr.shape

        # Fleshing out index 0 of the simulation cube (trace of target)
        if (sci_dx == 0) & (sci_dy == 0):  # this is the target

            simuCube[0, 0:trY, 0:trX] = tr

        # Fleshing out indexes 1-361 of the simulation cube
        # (trace of neighboring stars at every position angle)
        else:

            simuCube[V3PA + 1, 0:trY, 0:trX] += tr

    return simuCube


def field_simulation(ra, dec, aperture, binComp='', n_jobs=-1, nPA=360, plot=True):
    """Produce a contamination field simulation at the given sky coordinates

    Parameters
    ----------
    ra : float
        The RA of the target
    dec : float
        The Dec of the target
    aperture: str
        The aperture to use, ['NIS_SUBSTRIP96', 'NIS_SUBSTRIP256', 'NRCA5_GRISM256_F444W', 'NRCA5_GRISM256_F322W2', 'MIRI_SLITLESSPRISM']
    binComp : sequence
        The parameters of a binary companion
    n_jobs: int
        Number of cores to use (-1 = All)
    nPA: int
        The number of position angles
    plot: bool
        Return contamination plot

    Returns
    -------
    simuCube : np.ndarray
        The simulated data cube. Index 0 and 1 (axis=0) show the trace of
        the target for orders 1 and 2 (respectively). Index 2-362 show the trace
        of the target at every position angle (PA) of the instrument.
    plt: NoneType, bokeh.plotting.figure
        The plot of the contaminationas a function of PA
    """
    # Check for contam tool data
    check_for_data('exoctk_contam')

    # Time it
    print('Calculating target contamination from neighboring sources...')
    start = time.time()

    # Aperture names
    if aperture not in APERTURES:
        raise ValueError("{}: Not a supported aperture. Try {}".format(aperture, APERTURES.keys()))

    # Instantiate a pySIAF object
    inst = APERTURES[aperture]
    siaf = pysiaf.Siaf(inst['inst'])

    # Get the full and subarray apertures
    full = siaf.apertures[inst['full']]
    aper = siaf.apertures[aperture]

    # Full frame pixel positions
    rows, cols = full.corners('det')
    pixel_pos = rows.min(), rows.max(), cols.min(), cols.max()

    # Converting to degrees
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=u.deg)
    targetRA = targetcrd.ra.value
    targetDEC = targetcrd.dec.value

    # Querying for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    info = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone', radius=inst['rad'] * u.arcmin)

    # Coordinates of all the stars in FOV, including target
    stars = {}
    stars['RA'] = info['ra'].data.data
    stars['DEC'] = info['dec'].data.data
    nStars = stars['DEC'].size

    # Distance to each FOV star
    sindRA = (targetRA - stars['RA']) * np.cos(targetDEC)
    cosdRA = targetDEC - stars['DEC']
    stars['distance'] = np.sqrt(sindRA**2 + cosdRA**2)

    # The target is the "closest" star from the Irsa query
    targetIndex = np.argmin(stars['distance'])

    # 2MASS - Teff relations
    jhMod = np.array([0.545, 0.561, 0.565, 0.583, 0.596, 0.611, 0.629, 0.642, 0.66, 0.679, 0.696, 0.71, 0.717, 0.715, 0.706, 0.688, 0.663, 0.631, 0.601, 0.568, 0.537, 0.51, 0.482, 0.457, 0.433, 0.411, 0.39, 0.37, 0.314, 0.279])
    hkMod = np.array([0.313, 0.299, 0.284, 0.268, 0.257, 0.247, 0.24, 0.236, 0.229, 0.217,0.203, 0.188, 0.173, 0.159, 0.148, 0.138, 0.13, 0.123, 0.116, 0.112, 0.107, 0.102, 0.098, 0.094, 0.09, 0.086, 0.083, 0.079, 0.07, 0.067])
    teffMod = np.array([2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5800, 6000])

    # JHK bands of all stars in FOV, including target
    stars['Jmag'] = info['j_m'].data.data
    stars['Hmag'] = info['h_m'].data.data
    stars['Kmag'] = info['k_m'].data.data

    # Add any missing companion
    if binComp != '':
        deg2rad = np.pi / 180
        bb = binComp[0] / 3600 / np.cos(allDEC[targetIndex] * deg2rad)
        stars['RA'] = np.append(stars['RA'], (stars['RA'][targetIndex] + bb))
        stars['DEC'] = np.append(stars['DEC'], (stars['DEC'][targetIndex] + binComp[1] / 3600))
        stars['Jmag'] = np.append(stars['Jmag'], binComp[2])
        stars['Hmag'] = np.append(stars['Hmag'], binComp[3])
        stars['Kmag'] = np.append(stars['Kmag'], binComp[4])

    # J-H band, H-K band. This will be used to derive the Teff
    J_Hobs = stars['Jmag'] - stars['Hmag']
    H_Kobs = stars['Hmag'] - stars['Kmag']

    # Find/assign Teff of each star
    stars['Temp'] = np.empty(nStars)
    for j in range(nStars):
        color_separation = (J_Hobs[j] - jhMod)**2 + (H_Kobs[j] - hkMod)**2
        min_separation_ind = np.argmin(color_separation)
        stars['Temp'][j] = teffMod[min_separation_ind]

    # Print star data
    table = Table(stars)
    print("{} stars in target FOV:".format(len(table)))
    table.pprint()

    # Set the number of cores for multiprocessing
    max_cores = cpu_count()
    if n_jobs == -1 or n_jobs > max_cores:
        n_jobs = max_cores

    # Calculate contamination at each PA
    pl = pool.ThreadPool(n_jobs)
    func = partial(calc_v3pa, stars=stars, aper=aper, targetIndex=targetIndex, pixel_pos=pixel_pos)
    cubes = zip(*pl.map(func, np.arange(0, nPA, 1)))
    pl.close()
    pl.join()

    # Sum cubes into a single cube
    simuCube = np.nansum(np.asarray(list(cubes)).swapaxes(0, 1), axis=0)
    print('Contamination calculation complete: {} {}'.format(round(time.time() - start, 3), 's'))

    plt = None
    if plot:

        # Convert coords to format gtvt can read
        ra_hms, dec_dms = re.sub('[a-z]', ':', targetcrd.to_string('hmsdms')).split(' ')

        # Determine PA visibility
        minPA, maxPA, _, _, _, badPAs = using_gtvt(ra_hms[:-1], dec_dms[:-1], inst['inst'])

        # Make the plot
        plt = plot_contamination(simuCube, inst['lam'], badPAs)

    return simuCube, plt


def get_trace(aperture, teff):
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
    traces_path = os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/traces/{}/*.fits'.format(aperture))

    # Glob the file names
    trace_files = glob.glob(traces_path)

    # Get closest Teff
    teffs = np.array([int(os.path.basename(file).split('_')[-1][:-5]) for file in trace_files])
    file = trace_files[np.argmin((teffs - teff)**2)]

    # Make into dict
    trace = fits.getdata(file)

    return trace.squeeze()


def plot_contamination(data, wlims, badPAs=[], minPA=0, maxPA=360, title=''):
    """
    Plot the contamination

    Parameters
    ----------
    data: np.ndarray
        The cube of data
    wlims: tuple
        The wavelength min and max
    badPAs: list
        The list of position angles with no visibility
    minPA: int
        The minimum position angle to plot
    maxPA: int
        The maximum position angle to plot

    Returns
    -------
    bokeh.layouts.gridplot
        The contamination figure
    """
    # Target star and interlopers
    targ = data[0, :, :]
    cube = data[1:, :, :]

    fig = figure()
    fig.image([targ], x=0, y=0, dw=targ.shape[1], dh=targ.shape[0])
    show(fig)

    # Data dimensions
    nPA, rows, cols = cube.shape
    dPA = 1
    PA = np.arange(minPA, maxPA, dPA)

    # Plot setup
    tools = 'pan, box_zoom, crosshair, reset, hover'
    xlim0, xlim1 = wlims
    ylim0 = minPA - 0.5
    ylim1 = maxPA + 0.5
    color_mapper = LinearColorMapper(palette=PuBu[8][::-1][2:], low=-4, high=1)
    color_mapper.low_color = 'white'
    color_mapper.high_color = 'black'

    # The width of the target trace
    peak = targ.max()
    low_lim_col = np.where(targ > 0.0001 * peak)[1].min()
    high_lim_col = np.where(targ > 0.0001 * peak)[1].max()

    # The length of the target trace
    targ_trace_start = np.where(targ > 0.0001 * peak)[0].min()
    targ_trace_stop = np.where(targ > 0.0001 * peak)[0].max()

    # Begin contam calculation at each channel (row) Y
    contam = np.zeros([rows, nPA])
    for Y in np.arange(rows):

        # Calculate limits
        peakX = np.argmax(targ[Y, :])
        LEFT, RIGHT = peakX - low_lim_col, peakX + high_lim_col

        # Calculate weights
        tr = targ[Y, LEFT:RIGHT]
        wt = tr / np.sum(tr**2)

        # Add to the contam figure
        contam[Y, :] = np.sum(cube[:, Y, LEFT:RIGHT] * wt, where=~np.isnan(cube[:, Y, LEFT:RIGHT] * wt), axis=1)

    # Trim the figure
    contam = contam[targ_trace_start:targ_trace_stop, :]

    # Prep figure data
    fig_data = np.log10(np.clip(contam, 1.e-10, 1.))

    # Trace plot
    trplot = figure(tools=tools, width=500, height=500, title=title, x_range=Range1d(xlim0, xlim1), y_range=Range1d(ylim0, ylim1))
    trplot.image([fig_data], x=xlim0, y=ylim0, dw=xlim1 - xlim0, dh=ylim1 - ylim0, color_mapper=color_mapper)
    trplot.xaxis.axis_label = 'Wavelength (um)'
    trplot.yaxis.axis_label = 'Aperture Position Angle (degrees)'

    # Shade bad position angles on the trace plot
    nbadPA = len(badPAs)
    if nbadPA > 0:
        tops = [np.max(badPA) for badPA in badPAs]
        bottoms = [np.min(badPA) for badPA in badPAs]
        left = [xlim0] * nbadPA
        right = [xlim1] * nbadPA
        trplot.quad(top=tops, bottom=bottoms, left=left, right=right, color='#555555', alpha=0.6)

    # Make a figure summing the contamination at a given PA
    sumplot = figure(tools=tools, width=150, height=500, x_range=Range1d(0, 100), y_range=trplot.y_range, title=None)
    sumplot.line(100 * np.sum(contam >= 0.001, axis=1) / rows, PA - dPA / 2, line_color='blue', legend_label='> 0.001')
    sumplot.line(100 * np.sum(contam >= 0.01, axis=1) / rows, PA - dPA / 2, line_color='green', legend_label='> 0.01')
    sumplot.xaxis.axis_label = '% channels contam.'
    sumplot.yaxis.major_label_text_font_size = '0pt'

    final = gridplot(children=[[trplot, sumplot]])

    return final


if __name__ == '__main__':
    ra, dec = "04 25 29.0162", "-30 36 01.603"  # Wasp 79
    field_simulation(ra, dec, 'NIS_SUBSTRIP256')

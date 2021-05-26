import glob
import os
from multiprocessing import pool, cpu_count
import time
from functools import partial

import astropy.coordinates as crd
import astropy.units as u
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


APERTURES = {'NIS_SUBSTRIP96': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'traces': 'NIRISS/modelOrder12_teff*', 'scale': 0.065, 'rad': 2.5, 'lam': [(0.8, 2.8), (0.6, 1.4)]},
             'NIS_SUBSTRIP256': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'traces': 'NIRISS/modelOrder12_teff*', 'scale': 0.065, 'rad': 2.5, 'lam': [(0.8, 2.8), (0.6, 1.4)]},
             'NRCA5_GRISM256_F444W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'traces': 'NIRCam_F444W/*', 'scale': 0.063, 'rad': 2.5, 'lam': [(3.063, 5.111)]},
             'NRCA5_GRISM256_F322W2': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'traces': 'NIRCam_F322W2/*', 'scale': 0.063, 'rad': 2.5, 'lam': [(2.369, 4.417)]},
             'MIRIM_SLITLESSPRISM': {'inst': 'MIRI', 'full': 'MIRIM_FULL', 'traces': 'MIRI/_*', 'scale': 0.11, 'rad': 2.0, 'lam': [(5, 12)]}}


def get_traces(aperture):
    """Get the traces for the given aperture

    Parameters
    ----------
    aperture: str
        The aperture to use

    Returns
    -------
    np.ndarray
        A data cube of the traces
    """
    # Get aperture info
    aper = APERTURES[aperture]

    # Get the path to the trace files
    traces_path = os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/traces/', aper['traces'])

    # Glob the file names
    trace_files = glob.glob(traces_path)

    # Make into dict
    trace_dict = {os.path.basename(file).split('_')[-1].replace('teff', '').split('.')[0]: fits.getdata(file)[0] if aper['inst'] == 'MIRI' else fits.getdata(file, 1)[0] if aper['inst'] == 'NIRCam' else readsav(file, verbose=False)['modelo12'] for file in trace_files}

    return trace_dict


def field_simulation(ra, dec, aperture, binComp='', n_jobs=-1, nPA=360):
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

    Returns
    -------
    simuCube : np.ndarray
        The simulated data cube. Index 0 and 1 (axis=0) show the trace of
        the target for orders 1 and 2 (respectively). Index 2-362 show the trace
        of the target at every position angle (PA) of the instrument.
    """
    # Check for contam tool data
    check_for_data('exoctk_contam')

    # Time it
    print('Starting contam cube...')
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
    subX, subY = aper.XSciSize, aper.YSciSize

    # Generate cube of field simulation at every degree of APA rotation
    V3PAs = np.arange(0, nPA, 1)
    xSweet, ySweet = aper.reference_point('det')
    add_to_v3pa = aper.V3IdlYAngle

    # Full Frame dimensions
    rows, cols = full.corners('det')
    minrow, maxrow = rows.min(), rows.max()
    mincol, maxcol = cols.min(), cols.max()

    # Determine the traces for the given instrument
    traces = get_traces(aperture)

    # Converting to degrees
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=u.deg)
    targetRA = targetcrd.ra.value
    targetDEC = targetcrd.dec.value

    # Querying for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    info = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone', radius=inst['rad'] * u.arcmin)

    # Coordinates of all the stars in FOV, including target
    allRA = info['ra'].data.data
    allDEC = info['dec'].data.data

    # Initiating a dictionary to hold all relevant star information
    stars = {}
    stars['RA'], stars['DEC'] = allRA, allDEC

    sindRA = (targetRA - stars['RA']) * np.cos(targetDEC)
    cosdRA = targetDEC - stars['DEC']
    distance = np.sqrt(sindRA**2 + cosdRA**2)
    # if np.min(distance) > 1.0*(10**-4):
    #     coords = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg)).to_string('decimal')
    #     ra, dec = coords.split(' ')[0], coords.split(' ')[1]
    #     raise Exception('Unable to detect a source with coordinates [RA: {}, DEC: {}] within IRSA`s 2MASS Point-Source Catalog. Please enter different coordinates or contact the JWST help desk.'.format(str(ra), str(dec)))

    # The target is the "closest" star from the Irsa query
    targetIndex = np.argmin(distance)

    # 2MASS - Teff relations
    jhMod = np.array([0.545, 0.561, 0.565, 0.583, 0.596, 0.611, 0.629, 0.642, 0.66, 0.679, 0.696, 0.71, 0.717, 0.715, 0.706, 0.688, 0.663, 0.631, 0.601, 0.568, 0.537, 0.51, 0.482, 0.457, 0.433, 0.411, 0.39, 0.37, 0.314, 0.279])
    hkMod = np.array([0.313, 0.299, 0.284, 0.268, 0.257, 0.247, 0.24, 0.236, 0.229, 0.217,0.203, 0.188, 0.173, 0.159, 0.148, 0.138, 0.13, 0.123, 0.116, 0.112, 0.107, 0.102, 0.098, 0.094, 0.09, 0.086, 0.083, 0.079, 0.07, 0.067])
    teffMod = np.array([2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5800, 6000])

    # JHK bands of all stars in FOV, including target
    Jmag = info['j_m'].data.data
    Hmag = info['h_m'].data.data
    Kmag = info['k_m'].data.data

    # J-H band, H-K band. This will be used to derive the Teff
    J_Hobs = Jmag - Hmag
    H_Kobs = Hmag - Kmag

    # Add any missing companion
    if binComp != '':
        deg2rad = np.pi / 180
        bb = binComp[0] / 3600 / np.cos(allDEC[targetIndex] * deg2rad)
        allRA = np.append(allRA, (allRA[targetIndex] + bb))
        allDEC = np.append(allDEC, (allDEC[targetIndex] + binComp[1] / 3600))
        Jmag = np.append(Jmag, binComp[2])
        Hmag = np.append(Kmag, binComp[3])
        Kmag = np.append(Kmag, binComp[4])
        J_Hobs = Jmag - Hmag
        H_Kobs = Hmag - Kmag

    # Number of stars
    nStars = stars['RA'].size

    # Find/assign Teff of each star
    starsT = np.empty(nStars)
    for j in range(nStars):
        color_separation = (J_Hobs[j] - jhMod)**2 + (H_Kobs[j] - hkMod)**2
        min_separation_ind = np.argmin(color_separation)
        starsT[j] = teffMod[min_separation_ind]

    # Record keeping
    stars['Temp'] = starsT
    stars['Jmag'] = Jmag

    # Calculate corresponding V2/V3 (TEL) coordinates for Sweetspot
    v2targ, v3targ = aper.det_to_tel(xSweet, ySweet)

    # Set the number of cores for multiprocessing
    max_cores = cpu_count()
    if n_jobs == -1 or n_jobs > max_cores:
        n_jobs = max_cores

    # Multiprocess PAs
    pl = pool.ThreadPool(n_jobs)
    func = partial(calc_v3pa, add_to_v3pa=add_to_v3pa, v2targ=v2targ, v3targ=v3targ, targetRA=targetRA, targetDEC=targetDEC, stars=stars, traces=traces, subY=subY, subX=subX, aper=aper, targetIndex=targetIndex, nStars=nStars, mincol=mincol, maxcol=maxcol, minrow=minrow, maxrow=maxrow)
    images = zip(*pl.map(func, V3PAs))
    pl.close()
    pl.join()

    # Make frames into a cube
    simuCube = np.asarray(list(images))
    print('Contam cube finished: {} {}'.format(round(time.time() - start, 3), 's'))

    return simuCube


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

    # The width of the trace
    peak = targ.max()
    low_lim_col = np.where(targ > 0.0001 * peak)[1].min()
    high_lim_col = np.where(targ > 0.0001 * peak)[1].max()

    # The length of the trace
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

    show(final)


def calc_v3pa(V3PA, add_to_v3pa, v2targ, v3targ, targetRA, targetDEC, stars, traces, subY, subX, aper, targetIndex, nStars, mincol, maxcol, minrow, maxrow):
    """
    Function to make a frame of the simulated cube for the given PA

    Returns
    -------
    np.ndarray
        The 2D array
    """
    # Get APA from V3PA
    APA = V3PA + add_to_v3pa
    if APA > 360:
        APA = APA-360
    elif APA < 0:
        APA = APA+360

    # Get target's attitude matrix for each Position Angle
    attitude = rotations.attitude_matrix(v2targ, v3targ, targetRA, targetDEC, APA)

    xdet, ydet = [], []
    xsci, ysci = [], []
    for starRA, starDEC in zip(stars['RA'], stars['DEC']):

        # Get the TEL coordinates of each star w attitude matrix
        V2, V3 = rotations.sky_to_tel(attitude, starRA, starDEC)

        # Convert to arcsec and turn to a float
        V2, V3 = V2.to(u.arcsec).value, V3.to(u.arcsec).value

        XDET, YDET = aper.tel_to_det(V2, V3)
        XSCI, YSCI = aper.det_to_sci(XDET, YDET)

        xdet.append(XDET)
        ydet.append(YDET)
        xsci.append(XSCI)
        ysci.append(YSCI)

    # Record keeping
    stars['xdet'], stars['ydet'] = np.array(xdet), np.array(ydet)
    stars['xsci'], stars['ysci'] = np.array(xsci), np.array(ysci)

    sci_targx = stars['xsci'][targetIndex]
    sci_targy = stars['ysci'][targetIndex]

    #########################STEP 5#####################################
    ####################################################################

    inFOV = []
    for star in range(0, nStars):

        x, y = stars['xdet'][star], stars['ydet'][star]
        if (mincol < x) & (x < maxcol) & (minrow < y) & (y < maxrow):
            inFOV.append(star)

    inFOV = np.array(inFOV)

    #########################STEP 6#####################################
    ####################################################################

    simuImage = np.zeros((subY, subX))
    for idx in inFOV:

        # Get the Teff
        sci_dx = round(sci_targx - stars['xsci'][idx])
        sci_dy = round(sci_targy - stars['ysci'][idx])
        temp = stars['Temp'][idx]

        # Get the flux scaling
        fluxscale = 10.0**(-0.4 * (stars['Jmag'][idx] - stars['Jmag'][targetIndex]))

        # Get the trace and pad it
        trace = traces[str(int(temp))]
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
        simuImage[0:trY, 0:trX] += tr

    return simuImage



if __name__ == '__main__':
    ra, dec = "04 25 29.0162", "-30 36 01.603"  # Wasp 79
    field_simulation(ra, dec, 'NIS_SUBSTRIP256')

from functools import partial
import glob
from multiprocessing import pool, cpu_count
import os
import re
import time

import astropy.coordinates as crd
import astropy.units as u
from astropy.table import Table
from astroquery.irsa import Irsa
from astropy.io import fits
from bokeh.plotting import figure, show
from bokeh.layouts import gridplot
from bokeh.models import Range1d, LinearColorMapper, Label, ColorBar, ColumnDataSource, HoverTool
from bokeh.models.widgets import Panel, Tabs
from bokeh.palettes import PuBu, Spectral6
from bokeh.transform import linear_cmap
from hotsoss.plotting import plot_frame
import numpy as np
import pysiaf

from ..utils import get_env_variables, check_for_data
from .visibilityPA import using_gtvt
from .contamination_figure import contam


APERTURES = {'NIS_SUBSTRIP96': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.8, 2.8], 'trim': [47, 46, 0, 1]},
             'NIS_SUBSTRIP256': {'inst': 'NIRISS', 'full': 'NIS_SOSSFULL', 'scale': 0.065, 'rad': 2.5, 'lam': [0.8, 2.8], 'trim': [127, 126, 0, 1]},
             'NRCA5_GRISM256_F277W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [2.395, 3.179], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F322W2': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [2.413, 4.083], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F356W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [3.100, 4.041], 'trim': [0, 1, 0, 1]},
             'NRCA5_GRISM256_F444W': {'inst': 'NIRCam', 'full': 'NRCA5_FULL', 'scale': 0.063, 'rad': 2.5, 'lam': [3.835, 5.084], 'trim': [0, 1, 1250, 1]},
             'MIRIM_SLITLESSPRISM': {'inst': 'MIRI', 'full': 'MIRIM_FULL', 'scale': 0.11, 'rad': 2.0, 'lam': [5, 12], 'trim': [6, 5, 0, 1]}}


def calc_v3pa(V3PA, stars, aperture, ref='sci', plot=True, verbose=True):
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

    # Get subarray dims
    subX, subY = aperture.XSciSize, aperture.YSciSize

    # Calculate corresponding V2/V3 (TEL) coordinates for Sweetspot
    stars['xdet'][0], stars['ydet'][0] = aperture.reference_point('det')
    stars['xtel'][0], stars['ytel'][0] = aperture.det_to_tel(stars['xdet'][0], stars['ydet'][0])
    stars['xsci'][0], stars['ysci'][0] = aperture.det_to_sci(stars['xdet'][0], stars['ydet'][0])

    # Get target's attitude matrix for each Position Angle
    attitude = pysiaf.utils.rotations.attitude_matrix(stars['xtel'][0], stars['ytel'][0], stars['ra'][0], stars['dec'][0], APA)

    # Get relative coordinates of the stars based on target attitude
    for idx, star in enumerate(stars[1:]):

        # Get the TEL coordinates (V2, V3) of the star
        V2, V3 = pysiaf.utils.rotations.sky_to_tel(attitude, star['ra'], star['dec'])
        star['xtel'], star['ytel'] = V2.to(u.arcsec).value, V3.to(u.arcsec).value

        # Get the DET coordinates of the star
        star['xdet'], star['ydet'] = aperture.tel_to_det(star['xtel'], star['ytel'])

        # Get the DET coordinates of the star
        star['xsci'], star['ysci'] = aperture.det_to_sci(star['xdet'], star['ydet'])

    # Aperture info
    aper = APERTURES[aperture.AperName]
    xleft, xright, ybot, ytop = aper['trim']

    # Just stars in FOV (Should always have at least 1, the target)
    FOVstars = stars[(aperture.mincol < stars['xdet']) & (stars['xdet'] < aperture.maxcol) & (aperture.minrow < stars['ydet']) & (stars['ydet'] < aperture.maxrow)]
    if verbose:
        print("Calculating contamination from {} other stars in the FOV".format(len(FOVstars) - 1))

    # Make frame for the target and a frame for all the other stars
    targframe = np.zeros((subY, subX))
    starframe = np.zeros((subY, subX))

    if plot:
        # Set up hover tool
        tips = [('Name', '@name'), ('RA', '@ra'), ('DEC', '@dec'), ('Jmag', '@j_m'), ('Hmag', '@h_m'), ('Kmag', '@k_m'), ('Teff', '@Teff'), ('Sci', '(@xsci{int}, @ysci{int})'), ('Det', '(@xdet{int}, @ydet{int})'), ('Tel', '(@xtel{int}, @ytel{int})')]
        hover = HoverTool(tooltips=tips, names=['stars'])

        # Make the plot
        tools = ['pan', 'reset', 'box_zoom', 'wheel_zoom', 'save', hover]
        fig = figure(title='Generated FOV from 2MASS IRSA fp_psc', match_aspect=True, tools=tools)
        fig.title = 'The FOV in {} coordinates at APA {} for {}'.format(ref.upper(), V3PA, aperture.AperName)
        fig.patch(full.corners(ref)[0], full.corners(ref)[1], color="black", alpha=0.1)
        fig.patch(aperture.corners(ref)[0], aperture.corners(ref)[1], line_color="blue", fill_color='blue', fill_alpha=0.1)

    # Iterate over all stars in the FOV and add their scaled traces to the correct frame
    for idx, star in enumerate(FOVstars):

        # Get the trace and shift into the correct subarray position
        trace = get_trace(aperture.AperName, star['Teff'])
        trace = trace[xleft:-xright, ybot:-ytop]
        x, y = int(star['x{}'.format(ref)]), int(star['y{}'.format(ref)])

        if 'NIS' in aperture.AperName:
            trace = trace.T[::-1]
            height, width = trace.shape
            x0 = x - width + 68
            y0 = y - height

        elif 'F322W2' in aperture.AperName:
            height, width = trace.shape
            x0 = x - width + 467  # 2048 - 1581
            y0 = y - round(height / 2)

        elif 'F356W' in aperture.AperName:
            height, width = trace.shape
            x0 = x - width + 467  # 2048 - 1581
            y0 = y - round(height / 2)

        elif 'F277W' in aperture.AperName:
            height, width = trace.shape
            x0 = x - width - 600
            y0 = y - round(height / 2)

        elif 'F444W' in aperture.AperName:
            trace = trace.T[:, ::-1]
            height, width = trace.shape
            x0 = x - width + 1096  # 2048 - 952
            y0 = y - round(height / 2)

        elif 'MIRI' in aperture.AperName:
            trace = trace[::-1]
            height, width = trace.shape
            x0 = x - round(width / 2)
            y0 = y - height + 113  # 529 - 416

        else:
            pass

        # print(star[['xdet', 'ydet', 'xsci', 'ysci']], width, height, x, y, subX, subY)

        # Scale flux
        trace[trace < np.nanmax(trace)/10000] = 0
        trace *= star['fluxscale']

        # Add target trace to target frame...
        if idx == 0:

            targframe[y0:y0 + height, x0:x0 + width] = trace

        # ..or star trace to star frame
        else:

            # TODO: Simulate 256 subarray trace for 96 and trim appropriately
            # Add just the portion of this star that falls on the target frame
            f0x, f1x = max(0, x0), min(subX, x0 + width)
            f0y, f1y = max(0, y0), min(subY, y0 + height)
            t0x, t1x = max(0, -x0), min(width, subX - x0)
            t0y, t1y = max(0, -y0), min(height, subY - y0)
            if t1y - t0y > 0 and t1x - t0x > 0:
                if verbose:
                    print("{} x {} pixels of star {} fall on the target frame".format(t1y-t0y, t1x-t0x, idx))
                starframe[f0y:f1y, f0x:f1x] += trace[t0y:t1y, t0x:t1x]

        if plot:
            fig.image(image=[trace], x=x0, dw=width, y=y0, dh=height, alpha=0.5)

    if plot:
        mapper = linear_cmap(field_name='Teff', palette=Spectral6, low=min(stars['Teff']), high=max(stars['Teff']))
        fig.star('x{}'.format(ref), 'y{}'.format(ref), color=mapper, size=12, name='stars', source=dict(stars[['Teff', 'x{}'.format(ref), 'y{}'.format(ref), 'ra', 'dec', 'name', 'j_m', 'h_m', 'k_m', 'xdet', 'ydet', 'xtel', 'ytel']]))
        fig.circle(FOVstars['x{}'.format(ref)][0], FOVstars['y{}'.format(ref)][0], size=12, fill_color=None, line_color='red')
        color_bar = ColorBar(color_mapper=mapper['transform'], width=10,  location=(0,0), title="Teff")
        fig.add_layout(color_bar, 'right')
        show(fig)

        show(plot_frame(targframe + starframe, title='Target'))
        # show(plot_frame(starframe, title='Contamination'))

    return targframe, starframe


def field_simulation(ra, dec, aperture, binComp='', n_jobs=-1, pa_list=np.arange(360), plot=True, multi=False):
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
    pa_list: sequence
        The position angles to calculate

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

    # Aperture names
    if aperture not in APERTURES:
        raise ValueError("Aperture '{}' not supported. Try {}".format(aperture, list(APERTURES.keys())))

    # Instantiate a pySIAF object
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
    aper.mincol, aper.maxcol = cols.min(), cols.max()

    # Find stars in the vicinity
    stars = find_stars(ra, dec, radius=inst['rad'] * u.arcmin, binComp=binComp)

    # Time it
    print('Calculating target contamination from {} neighboring sources...'.format(len(stars)))
    start = time.time()

    # Set the number of cores for multiprocessing
    max_cores = cpu_count()
    if n_jobs == -1 or n_jobs > max_cores:
        n_jobs = max_cores

    # Exclude PAs where target is not visible to speed up calculation
    ra_hms, dec_dms = re.sub('[a-z]', ':', targetcrd.to_string('hmsdms')).split(' ')
    minPA, maxPA, _, _, _, badPAs = using_gtvt(ra_hms[:-1], dec_dms[:-1], inst['inst'])
    pa_list = [pa for pa in pa_list if pa not in badPAs]

    # Calculate contamination of all stars at each PA
    # To multiprocess, or not to multiprocess. That is the question.
    if multi:
        pl = pool.ThreadPool(n_jobs)
        func = partial(calc_v3pa, stars=stars, aperture=aper, plot=False, verbose=False)
        targframes, starframes = zip(*pl.map(func, pa_list))
        pl.close()
        pl.join()

    else:
        targframes = []
        starframes = []
        for pa in pa_list:
            print(pa)
            tarf, starf = calc_v3pa(pa, stars=stars, aperture=aper, plot=False, verbose=False)
            targframes.append(tarf)
            starframes.append(starf)

    # We only need one target frame frames
    targframe = np.asarray(targframes[0])

    # Make sure starcube is of shape (PA, rows, cols)
    starcube = np.zeros((360, subY, subX))
    for pa, starframe in zip(pa_list, starframes):
        starcube[pa, :, :] = starframe

    print('Contamination calculation complete: {} {}'.format(round(time.time() - start, 3), 's'))

    # Make the contamination plot
    plt = None
    if plot:

        # Old plot
        # starcube = np.concatenate([targframe[None, :, :], targframe[None, :, :], starcube], axis=0)
        starcube = starcube.swapaxes(1, 2)
        targframe = targframe.swapaxes(0, 1)
        starcube = starcube[:, :2048, :256]
        targframe = targframe[:2048, :256]
        # plt = contam(starcube, 'NIS_SUBSTRIP256', targetName='Foo', badPAs=badPAs, fig='bokeh')

        # New plot
        plt = plot_contamination(targframe, starcube, [0, 2048], badPAs)

    return targframe, starcube, plt


def find_stars(ra, dec, radius, binComp=''):
    """
    Find all the stars in the vicinity and estimate temperatures

    Parameters
    ----------
    ra : float
        The RA of the target in decimal degrees
    dec : float
        The Dec of the target in decimal degrees
    radius: astropy.units.quantity
        The search radius

    Returns
    -------
    astropy.table.Table
        The table of stars
    """
    # Converting to degrees and query for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=u.deg if isinstance(ra, float) and isinstance(dec, float) else (u.hour, u.deg))
    stars = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone', radius=radius)

    jhMod = np.array([0.545, 0.561, 0.565, 0.583, 0.596, 0.611, 0.629, 0.642, 0.66, 0.679, 0.696, 0.71, 0.717, 0.715, 0.706, 0.688, 0.663, 0.631, 0.601, 0.568, 0.537, 0.51, 0.482, 0.457, 0.433, 0.411, 0.39, 0.37, 0.314, 0.279])
    hkMod = np.array([0.313, 0.299, 0.284, 0.268, 0.257, 0.247, 0.24, 0.236, 0.229, 0.217,0.203, 0.188, 0.173, 0.159, 0.148, 0.138, 0.13, 0.123, 0.116, 0.112, 0.107, 0.102, 0.098, 0.094, 0.09, 0.086, 0.083, 0.079, 0.07, 0.067])
    teffMod = np.array([2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5800, 6000])

    # Make sure colors are calculated
    stars['j_h'] = stars['j_m'] - stars['h_m']
    stars['h_k'] = stars['h_m'] - stars['k_m']
    stars['j_k'] = stars['j_m'] - stars['k_m']

    # Add any missing companion (ra, dec, J, H, K)
    if binComp != '':
        deg2rad = np.pi / 180
        bb = binComp[0] / 3600 / np.cos(stars[0]['ra'] * deg2rad)
        star = {'ra': stars['ra'][0] + bb, 'dec': stars['dec'][0] + binComp[1] / 3600, 'j_m': binComp[2], 'h_m': binComp[3], 'k_m': binComp[4], 'j_h': binComp[2] - binComp[3], 'h_k': binComp[3] - binComp[4]}
        star['Teff'] = teffMod[np.argmin((star['j_h'] - jhMod) ** 2 + (star['h_k'] - hkMod) ** 2)]
        stars.add_row(star)

    # Find distance from target to each star
    sindRA = (stars['ra'][0] - stars['ra']) * np.cos(stars['dec'][0])
    cosdRA = stars['dec'][0] - stars['dec']
    stars.add_column(np.sqrt(sindRA ** 2 + cosdRA ** 2), name='distance')
    stars.sort('distance')

    # Find Teff of each star from the color
    teffs = [teffMod[np.argmin((row['j_h'] - jhMod) ** 2 + (row['h_k'] - hkMod) ** 2)] for row in stars]

    # Add Teffs and detector location to the table
    stars.add_column(teffs, name='Teff')
    stars.add_columns(np.zeros((6, len(stars))), names=['xtel', 'ytel', 'xdet', 'ydet', 'xsci', 'ysci'])

    # Add names for tooltips
    stars.add_column([i.decode('UTF-8') for i in stars['designation']], name='name')

    # Calculate relative flux (and make plot marker size proportional)
    fluxscale = 10.0 ** (-0.4 * (stars['j_m'] - stars['j_m'][0]))
    stars.add_column(fluxscale, name='fluxscale')

    return stars


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

    # Get data
    trace = fits.getdata(file)

    return trace.squeeze()[::-1]


def plot_contamination(targframe, starcube, wlims, badPAs=[], title=''):
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
    minPA: int
        The minimum position angle to plot
    maxPA: int
        The maximum position angle to plot

    Returns
    -------
    bokeh.layouts.gridplot
        The contamination figure
    """
    # Data dimensions
    PAs, rows, cols = starcube.shape

    # Remove background values < 1 as it can blow up contamination
    targframe = np.where(targframe < 1, 0, targframe)

    # The width of the target trace
    peak = targframe.max()
    low_lim_col = np.where(targframe > 0.0001 * peak)[1].min()
    high_lim_col = np.where(targframe > 0.0001 * peak)[1].max()

    # The length of the target trace
    targ_trace_start = np.where(targframe > 0.0001 * peak)[0].min()
    targ_trace_stop = np.where(targframe > 0.0001 * peak)[0].max()

    # # Calculate limits of the target trace
    # cutoff = targframe.max() / 0.0001
    # row_starts = np.argmax(targframe > cutoff, axis=1)
    # row_stops = -np.argmax(targframe[:, ::-1] > cutoff, axis=1)
    #
    # # Iterate over rows
    # contam = np.zeros([rows, PAs])
    # for row, (start, stop) in enumerate(zip(row_starts, row_stops)):
    #
    #     # Calculate weights
    #     tr = targframe[row, start:stop]
    #     wt = tr / np.sum(tr**2)
    #     ww = np.tile(wt, PAs).reshape([PAs, tr.size])
    #
    #     # Add to contam figure
    #     contam[row, :] = np.sum(starcube[:, row, start:stop] * ww, axis=1)

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
    # contam = np.clip(contam, 1.e-10, 1.)

    # Hover tool
    hover = HoverTool(tooltips=[("Wavelength", "$x"), ("PA", "$y"), ('Value', '@data')], names=['contam'])
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

    # Make a figure summing the contamination at a given PA
    sumplot = figure(tools=tools, width=150, height=500, x_range=Range1d(0, 100), y_range=trplot.y_range, title=None)
    sumplot.line(100 * np.sum(contam >= 0.001, axis=1) / rows, np.arange(PAs) - 0.5, line_color='blue', legend_label='> 0.001')
    sumplot.line(100 * np.sum(contam >= 0.01, axis=1) / rows, np.arange(PAs) - 0.5, line_color='green', legend_label='> 0.01')
    sumplot.xaxis.axis_label = '% channels contam.'
    sumplot.yaxis.major_label_text_font_size = '0pt'

    return gridplot(children=[[trplot, sumplot]])


import glob
import os
import pysiaf

import astropy.coordinates as crd
from astropy.io import fits
from astroquery.irsa import Irsa
import astropy.units as u
import numpy as np
from pysiaf.utils import rotations
from scipy.io import readsav

from exoctk import utils

EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
TRACES_PATH = os.path.join(os.environ.get('EXOCTK_DATA'), 'exoctk_contam', 'traces')


def sossFieldSim(ra, dec, binComp='', dimX=256, frame=0):
    """ Produce a SOSS field simulation for a target.
    Parameters
    ----------
    ra: float
        The RA of the target.
    dec: float
        The Dec of the target.
    binComp: sequence
        The parameters of a binary companion.
    dimX: int
        The subarray size.
    Returns
    -------
    simuCub : np.ndarray
        The simulated data cube.
    """

    # STEP 1
    # Pulling stars from IRSA point-source catalog
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=u.deg if isinstance(ra, float) and isinstance(dec, float) else (u.hour, u.deg))
    targetRA = targetcrd.ra.deg
    targetDEC = targetcrd.dec.deg
    info = Irsa.query_region(targetcrd,
                             catalog='fp_psc',
                             spatial='Cone',
                             radius=2.5 * u.arcmin)

    # Coordinates of all stars in FOV, including target
    allRA = info['ra'].data.data
    allDEC = info['dec'].data.data
    Jmag = info['j_m'].data.data
    Hmag = info['h_m'].data.data
    Kmag = info['k_m'].data.data

    # J-H band, H-K band. This will be used to derive the stellar Temps later
    J_Hobs = Jmag - Hmag
    H_Kobs = Hmag - Kmag

    # Determining target index by calculating the relative distance between
    # each source and the target. The target will have the smallest distance
    # from itself (oof) so whatever that index is will be the targetIndex
    aa = ((targetRA - allRA) * np.cos(targetDEC))
    distance = np.sqrt(aa ** 2 + (targetDEC - allDEC) ** 2)
    targetIndex = np.argmin(distance)

    # Add any missing companion
    if binComp != '':
        binComp = [float(i) for i in binComp.split(',')]

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
    nStars = allRA.size

    # Restoring model parameters
    modelParam = readsav(os.path.join(TRACES_PATH, 'NIRISS', 'modelsInfo.sav'),
                         verbose=False)
    models = modelParam['models']
    modelPadX = modelParam['modelpadx']
    modelPadY = modelParam['modelpady']
    dimXmod = modelParam['dimxmod']
    dimYmod = modelParam['dimymod']
    jhMod = modelParam['jhmod']
    hkMod = modelParam['hkmod']
    teffMod = modelParam['teffmod']

    # Find/assign Teff of each star
    starsT = np.empty(nStars)
    for j in range(nStars):
        color_separation = (J_Hobs[j] - jhMod) ** 2 + (H_Kobs[j] - hkMod) ** 2
        min_separation_ind = np.argmin(color_separation)
        starsT[j] = teffMod[min_separation_ind]

    sweetSpot = dict(x=856, y=107, RA=allRA[targetIndex],
                     DEC=allDEC[targetIndex], jmag=Jmag[targetIndex])

    radeg = 180 / np.pi
    niriss_pixel_scale = 0.065  # arcsec
    # offset between all stars and target
    dRA = (allRA - sweetSpot['RA']) * np.cos(sweetSpot['DEC'] / radeg) * 3600
    dDEC = (allDEC - sweetSpot['DEC']) * 3600

    # Put field stars positions and magnitudes in structured array
    _ = dict(RA=allRA, DEC=allDEC, dRA=dRA, dDEC=dDEC, jmag=Jmag, T=starsT,
             x=np.empty(nStars), y=np.empty(nStars), dx=np.empty(nStars),
             dy=np.empty(nStars), distance=distance)
    stars = np.empty(nStars, dtype=[(key, val.dtype) for key, val in _.items()])
    for key, val in _.items():
        stars[key] = val

    # Initialize final fits cube that contains the modelled traces
    # with contamination
    PAmin = 0  # instrument PA, degrees
    PAmax = 360
    dPA = 1  # degrees

    # Set of IPA values to cover
    PAtab = np.arange(PAmin, PAmax, dPA)  # degrees
    nPA = len(PAtab)

    dimY = 2048
    # cube of trace simulation at every degree of field rotation,
    # +target at O1 and O2
    simuCube = np.zeros([nPA + 2, dimY, dimX])

    saveFiles = glob.glob(
        os.path.join(
            TRACES_PATH,
            'NIRISS',
            '*modelOrder12*.sav'))

    # Big loop to generate a simulation at each instrument PA

    # for kPA in [frame]:  # range(PAtab.size):
    for kPA in range(PAtab.size):
        APA = PAtab[kPA]
        print('Generating field at APA : {}'.format(str(APA)))

        sindx = np.sin((np.pi / 2) + APA / radeg) * stars['dDEC']
        cosdx = np.cos((np.pi / 2) + APA / radeg) * stars['dDEC']
        nps = niriss_pixel_scale
        stars['dx'] = (np.cos((np.pi / 2) + APA / radeg) * stars['dRA'] - sindx) / nps
        stars['dy'] = (np.sin((np.pi / 2) + APA / radeg) * stars['dRA'] + cosdx) / nps
        stars['x'] = stars['dx'] + sweetSpot['x']
        stars['y'] = stars['dy'] + sweetSpot['y']

        # Retain stars that are within the Direct Image NIRISS POM FOV
        ind, = np.where(
            (stars['x'] >= -162) & (stars['x'] <= 2047 + 185) & (stars['y'] >= -154) & (stars['y'] <= 2047 + 174))
        starsInFOV = stars[ind]

        for i in range(len(ind)):
            intx = round(starsInFOV['dx'][i])
            inty = round(starsInFOV['dy'][i])

            k = np.where(teffMod == starsInFOV['T'][i])[0][0]

            fluxscale = 10.0 ** (-0.4 * (starsInFOV['jmag'][i] - sweetSpot['jmag']))

            # deal with subection sizes.
            # these variables will determine where the
            # trace will land on the array based on the
            # neighbor's position relative to the target's position
            mx0 = int(modelPadX - intx)
            mx1 = int(modelPadX - intx + dimX)
            my0 = int(modelPadY - inty)
            my1 = int(modelPadY - inty + dimY)

            if (mx0 > dimXmod) or (my0 > dimYmod):
                continue
            if (mx1 < 0) or (my1 < 0):
                continue

            x0 = (mx0 < 0) * (-mx0)
            y0 = (my0 < 0) * (-my0)
            mx0 *= (mx0 >= 0)
            mx1 = dimXmod if mx1 > dimXmod else mx1
            my0 *= (my0 >= 0)
            my1 = dimYmod if my1 > dimYmod else my1

            # if target and first kPA, add target traces of order 1 and 2
            # in output cube
            if (intx == 0) & (inty == 0) & (kPA == 0):
                fNameModO12 = saveFiles[k]

                modelO12 = readsav(fNameModO12, verbose=False)['modelo12']
                ord1 = modelO12[0, my0:my1, mx0:mx1] * fluxscale
                ord2 = modelO12[1, my0:my1, mx0:mx1] * fluxscale
                simuCube[0, y0:y0 + my1 - my0, x0:x0 + mx1 - mx0] = ord1
                simuCube[1, y0:y0 + my1 - my0, x0:x0 + mx1 - mx0] = ord2

            if (intx != 0) or (inty != 0):
                mod = models[k, my0:my1, mx0:mx1]
                simuCube[kPA + 2, y0:y0 + my1 - my0, x0:x0 + mx1 - mx0] += mod * fluxscale

    # fra = simuCube[frame + 2, :, :]
    # tar = simuCube[0, :, :]
    # ff = plot_frame(fra.T + tar.T)
    # show(ff)

    return simuCube, Table(starsInFOV)


if __name__ == '__main__':
    ra, dec = "04 25 29.0162", "-30 36 01.603"  # Wasp 79
    field_simulation(ra, dec, 'NIS_SUBSTRIP256')

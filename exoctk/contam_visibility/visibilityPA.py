"""Produces a graph of the visibility & accessible position angles
for a given RA & DEC, and prints out corresponding information,
including the ranges of accessible and inaccessible PAs.

Usage: python visibilityPA.py RA DEC [targetName]
 if targetName is specified, then the figure is saved

-Created by David Lafreniere, March 2016
-makes use of (and hacks) several scripts created by Pierre Ferruit
 that are part of the JWST Python tools JWSTpylib and JWSTpytools
"""
import datetime
import math
import pkg_resources

from astropy.table import Table
from astropy.time import Time
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import HoverTool, ranges
from bokeh.models.widgets import Panel, Tabs
import matplotlib.dates as mdates
import numpy as np

from . import ephemeris_old2x as EPH
from jwst_gtvt.find_tgt_info import get_table


D2R = math.pi / 180.  # degrees to radians
R2D = 180. / math.pi  # radians to degrees


def checkVisPA(ra, dec, targetName=None, ephFileName=None, fig=None):
    """Check the visibility at a range of position angles.

    Parameters
    ----------
    ra: str
        The RA of the target in hh:mm:ss.s or dd:mm:ss.s or representing a float
    dec: str
        The Dec of the target in hh:mm:ss.s or dd:mm:ss.s or representing a float
    targetName: str
        The target name
    ephFileName: str
        The filename of the ephemeris file
    fig: bokeh.plotting.figure
        The figure to plot to

    Returns
    -------
    paGood : float
        The good position angle.
    paBad : float
        The bad position angle.
    gd : matplotlib.dates object
       The greogrian date.
    fig : bokeh.plotting.figure object
        The plotted figure.

    """
    if ephFileName is None:
        eph_file = 'data/contam_visibility/JWST_ephem_short.txt'
        ephFileName = pkg_resources.resource_filename('exoctk', eph_file)
    if ra.find(':') > -1:  # format is hh:mm:ss.s or  dd:mm:ss.s
        ra = convert_ddmmss_to_float(ra) * 15. * D2R
        dec = convert_ddmmss_to_float(dec) * D2R
    else:  # format is decimal
        ra = float(ra) * D2R
        dec = float(dec) * D2R

    # load ephemeris
    eclFlag = False
    eph = EPH.Ephemeris(ephFileName, eclFlag)

    # convert dates from MJD to Gregorian calendar dates
    mjd = np.array(eph.datelist)
    d = mdates.julian2num(mjd + 2400000.5)
    gd = mdates.num2date(d)

    # loop through dates and determine VIS and PAs (nominal, min, max)
    vis = np.empty(mjd.size, dtype=bool)
    paNom = np.empty(mjd.size)
    paMin = np.empty(mjd.size)
    paMax = np.empty(mjd.size)
    for i in range(mjd.size):

        # is it visible?
        vis[i] = eph.in_FOR(mjd[i], ra, dec)

        # nominal PA at this date
        pa = eph.normal_pa(mjd[i], ra, dec)

        # search for minimum PA allowed by roll
        pa0 = pa
        while eph.is_valid(mjd[i], ra, dec, pa0 - 0.002):
            pa0 -= 0.002

        # search for maximum PA allowed by roll
        pa1 = pa
        while eph.is_valid(mjd[i], ra, dec, pa1 + 0.002):
            pa1 += 0.002

        paNom[i] = (pa * R2D) % 360
        paMin[i] = (pa0 * R2D) % 360
        paMax[i] = (pa1 * R2D) % 360

    # does PA go through 360 deg?
    wrap = np.any(np.abs(np.diff(paNom[np.where(vis)[0]])) > 350)

    # Determine good and bad PA ranges
    # Good PAs
    i, = np.where(vis)
    pa = np.concatenate((paNom[i], paMin[i], paMax[i]))

    if wrap:
        pa = np.append(pa, (0., 360.))
    pa.sort()

    i1, = np.where(np.diff(pa) > 10)
    i0 = np.insert(i1 + 1, 0, 0)
    i1 = np.append(i1, -1)
    paGood = np.dstack((pa[i0], pa[i1])).round(1).reshape(-1, 2).tolist()

    # bad PAs (complement of the good PAs)
    paBad = []
    if paGood[0][0] > 0:
        paBad.append([0., paGood[0][0]])
    for i in range(1, len(paGood)):
        paBad.append([paGood[i - 1][1], paGood[i][0]])
    if paGood[-1][1] < 360.:
        paBad.append([paGood[-1][1], 360.])

    # Make a figure
    if fig is None or fig:
        tools = 'crosshair, reset, hover, save'
        radec = ', '.join([str(ra), str(dec)])
        fig = figure(tools=tools, plot_width=800, plot_height=400,
                     x_axis_type='datetime',
                     title=targetName or radec)

    # Do all figure calculations
    iBad, = np.where(vis == False)
    paMasked = np.copy(paNom)
    paMasked[iBad] = np.nan
    gdMasked = np.copy(gd)

    i = np.argmax(paNom)
    if paNom[i + 1] < 10:
        i += 1
    paMasked = np.insert(paMasked, i, np.nan)
    gdMasked = np.insert(gdMasked, i, gdMasked[i])

    i = np.argmax(paMin)
    goUp = paMin[i - 2] < paMin[i - 1]  # PA going up at wrap point?

    # Top part
    i0_top = 0 if goUp else i
    i1_top = i if goUp else paMin.size - 1
    paMaxTmp = np.copy(paMax)
    paMaxTmp[np.where(paMin > paMax)[0]] = 360

    # Bottom part
    i = np.argmin(paMax)
    i0_bot = i if goUp else 0
    i1_bot = paMin.size - 1 if goUp else i
    paMinTmp = np.copy(paMin)
    paMinTmp[np.where(paMin > paMax)[0]] = 0

    # Convert datetime to a number for Bokeh
    gdMaskednum = [datetime.date(2019, 6, 1) + datetime.timedelta(days=n)
                   for n, d in enumerate(gdMasked)]
    color = 'green'

    # Draw the curve and error
    try:
        fig.line(
            gdMaskednum,
            paMasked,
            legend_label='cutoff',
            line_color=color)
    except AttributeError:
        fig.line(gdMaskednum, paMasked, legend='cutoff', line_color=color)

    # Top
    terr_y = np.concatenate([paMin[i0_top:i1_top + 1],
                             paMaxTmp[i0_top:i1_top + 1][::-1]])
    terr_x = np.concatenate([gdMaskednum[i0_top:i1_top + 1],
                             gdMaskednum[i0_top:i1_top + 1][::-1]])
    fig.patch(terr_x, terr_y, color=color, fill_alpha=0.2, line_alpha=0)

    # Bottom
    berr_y = np.concatenate([paMinTmp[i0_bot:i1_bot + 1],
                             paMax[i0_bot:i1_bot + 1][::-1]])
    berr_x = np.concatenate([gdMaskednum[i0_bot:i1_bot + 1],
                             gdMaskednum[i0_bot:i1_bot + 1][::-1]])
    fig.patch(berr_x, berr_y, color='red', fill_alpha=0.2, line_alpha=0)
    from bokeh.plotting import show
    show(fig)

    # Plot formatting
    fig.xaxis.axis_label = 'Date'
    fig.yaxis.axis_label = 'Aperture Position Angle (degrees)'

    return paGood, paBad, gd, fig


def using_gtvt(ra, dec, instrument, targetName='noName', ephFileName=None, output='bokeh'):
    """Plot the visibility (at a range of position angles) against time.

    Parameters
    ----------
    ra : str
        The RA of the target (in degrees) hh:mm:ss.s or dd:mm:ss.s or representing a float
    dec : str
        The Dec of the target (in degrees) hh:mm:ss.s or dd:mm:ss.s or representing a float
    instrument : str
        Name of the instrument. Can either be (case-sensitive):
        'NIRISS', 'NIRCam', 'MIRI', 'FGS', or 'NIRSpec'
    ephFileName : str
        The filename of the ephemeris file.
    output : str
        Switches on plotting with Bokeh. Parameter value must be 'bokeh'.

    Returns
    -------
    paGood : float
        The good position angle.
    paBad : float
        The bad position angle.
    gd : matplotlib.dates object
       The gregorian date.
    fig : bokeh.plotting.figure object
        The plotted figure.

    """
    # Getting calculations from GTVT (General Target Visibility Tool)
    tab = get_table(ra, dec)

    gd = tab['Date']
    paMin = tab[str(instrument) + ' min']
    paMax = tab[str(instrument) + ' max']
    paNom = tab[str(instrument) + ' nom']
    v3min = tab['V3PA min']
    v3max = tab['V3PA max']

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOTE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Addressing NIRSpec issue*
    # *the issue that NIRSpec's angle goes beyond 360 degrees with some targs,
    # thus resetting back to 0 degrees, which can make the plot look weird
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    index = np.arange(0, len(paNom), 1)

    for idx in index:

        try:
            a1 = paNom[idx]
            b1 = paNom[idx + 1]

            if (np.isfinite(a1)) & (np.isfinite(b1)):
                delta = np.abs(a1 - b1)

                if delta > 250:

                    gd = np.insert(gd, idx + 1, np.nan)
                    paMin = np.insert(paMin, idx + 1, np.nan)
                    paMax = np.insert(paMax, idx + 1, np.nan)
                    paNom = np.insert(paNom, idx + 1, np.nan)
                    v3min = np.insert(v3min, idx + 1, np.nan)
                    v3max = np.insert(v3min, idx + 1, np.nan)
        except BaseException:
            pass

    # Setting up HoverTool parameters & other variables
    COLOR = 'green'
    TOOLS = 'pan, wheel_zoom, box_zoom, reset, save'
    SOURCE = ColumnDataSource(data=dict(pamin=paMin,
                                        panom=paNom,
                                        pamax=paMax,
                                        date=gd))
    TOOLTIPS = [('Date', '@date{%F}'),
                ('Maximum Aperture PA', '@pamax'),
                ('Nominal Aperture PA', '@panom'),
                ('Minimum Aperture PA', '@pamin')]

    # Time to plot
    if output == 'bokeh':
        fig = figure(tools=TOOLS,
                     plot_width=800,
                     plot_height=400,
                     x_axis_type='datetime',
                     title='{} Visibility with {}'.format(targetName,
                                                          instrument))

    # Draw the curve and PA min/max circles
    try:
        nom = fig.line('date', 'panom',
                       line_color=COLOR,
                       legend_label='Nominal Aperture PA',
                       alpha=.5,
                       source=SOURCE)
    except AttributeError:
        nom = fig.line('date', 'panom',
                       line_color=COLOR,
                       legend='Nominal Aperture PA',
                       alpha=.5,
                       source=SOURCE)

    fig.circle('date', 'pamin', color=COLOR, size=1, source=SOURCE)
    fig.circle('date', 'pamax', color=COLOR, size=1, source=SOURCE)

    # Adding HoverTool
    fig.add_tools(HoverTool(renderers=[nom],
                            tooltips=TOOLTIPS,
                            formatters={'date': 'datetime'},
                            mode='vline'))

    # Plot formatting
    fig.xaxis.axis_label = 'Date'
    fig.yaxis.axis_label = 'Aperture Position Angle (degrees)'
    fig.y_range = ranges.Range1d(0, 360)

    # Making the output table
    # Creating new lists w/o the NaN values
    v3minnan, v3maxnan, paNomnan, paMinnan, paMaxnan, gdnan, mjds = \
        [], [], [], [], [], [], []

    for vmin, vmax, pnom, pmin, pmax, date in zip(
            v3min, v3max, paNom, paMin, paMax, gd):
        if np.isfinite(pmin):
            v3minnan.append(vmin)
            v3maxnan.append(vmax)
            paNomnan.append(pnom)
            paMinnan.append(pmin)
            paMaxnan.append(pmax)
            gdnan.append(date)

    # Adding MJD column
    mjdnan = []
    for date in gdnan:
        t = Time(str(date), format='iso')
        mjd = t.mjd
        mjdnan.append(mjd)

    # Adding lists to a table object
    table = Table([v3minnan,
                   v3maxnan,
                   paMinnan,
                   paMaxnan,
                   paNomnan,
                   gdnan,
                   mjdnan],
                  names=('min_V3_PA',
                         'max_V3_PA',
                         'min_Aperture_PA',
                         'max_Aperture_PA',
                         'nom_Aperture_PA',
                         'Gregorian',
                         'MJD'))

    # Getting bad PAs
    allPAs = np.arange(0, 360, 1)
    badPAs = []

    for pa in allPAs:
        if (pa not in np.round(paMinnan)) & \
           (pa not in np.round(paMaxnan)) & \
           (pa not in np.round(paNomnan)):

            badPAs.append(pa)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOTE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This addresses a bokeh shading issue that accidentally shades
    # accessible PAs (e.g: trappist-1b)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    badPAs = select_badPAs_ge_paNomnan(badPAs, paNomnan)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOTE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Grouping the bad PAs into lists within the badPAs list.
    # This will make bad PA shading easier in the contamination Bokeh plot
    # (sossContamFig.py)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    badPAs = np.sort(badPAs)

    if len(badPAs > 0):
        grouped_badPAs = [[badPAs[0]]]

        for idx in range(1, len(badPAs)):

            if ((badPAs[idx - 1] + 1) == badPAs[idx]):

                grouped_badPAs[len(grouped_badPAs) - 1].append(badPAs[idx])

            elif ((badPAs[idx - 1] + 1) < badPAs[idx]):
                grouped_badPAs.append([badPAs[idx]])

        grouped_badPAs = np.asarray(grouped_badPAs)

    else:  # Accounting for targets with 100% visibility
        grouped_badPAs = np.asarray([])

    return paMin, paMax, gd, fig, table, grouped_badPAs


def select_badPAs_ge_paNomnan(badPAs, paNomnan, threshold=7):
    """Returns the absolute difference between each badPAs and paNomnan
    Should be greater than threshold (default=7)

    Parameters
    ----------
    badPAs: list
        The list of bad position angles
    paNomnan: list
        The list of nominal PAs

    Returns
    -------
    np.ndarray
        The array of PAs
    """
    # Reshaping
    badPAs_array = np.array(badPAs)[np.newaxis]  # (1, len(badPAs))
    paNomnan_array = np.array(paNomnan)[np.newaxis].T  # (len(paNomnan), 1)

    # elementwise absolute difference
    diff = np.abs(np.subtract(badPAs_array, paNomnan_array))  # (len(paNomnan), len(badPAs))

    # boolean array above threshold
    above_thresh = np.all(diff >= threshold, axis=0)  # (len(badPAs),)

    # index and return those that are above threshold
    return badPAs_array[0, above_thresh]
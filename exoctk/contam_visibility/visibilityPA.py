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

from bokeh.plotting import figure
import matplotlib.dates as mdates
import numpy as np

from . import ephemeris_old2x as EPH
from jwst_gtvt.find_tgt_info import get_table


D2R = math.pi/180.  # degrees to radians
R2D = 180./math.pi  # radians to degrees


def convert_ddmmss_to_float(astring):
    """Convert sexigesimal to decimal degrees

    Parameters
    ----------
    astring: str
        The sexigesimal coordinate.

    Returns
    -------
    hour_or_deg : float
        The converted coordinate.
    """
    aline = astring.split(':')
    d = float(aline[0])
    m = float(aline[1])
    s = float(aline[2])
    hour_or_deg = (s/60.+m)/60.+d
    return hour_or_deg


def checkVisPA(ra, dec, targetName=None, ephFileName=None, fig=None):
    """Check the visibility at a range of position angles

    Parameters
    ----------
    ra: float
        The RA of the target
    dec: float
        The Dec of the target
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
    d = mdates.julian2num(mjd+2400000.5)
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
        while eph.is_valid(mjd[i], ra, dec, pa0-0.002):
            pa0 -= 0.002

        # search for maximum PA allowed by roll
        pa1 = pa
        while eph.is_valid(mjd[i], ra, dec, pa1+0.002):
            pa1 += 0.002

        paNom[i] = (pa*R2D) % 360
        paMin[i] = (pa0*R2D) % 360
        paMax[i] = (pa1*R2D) % 360

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
    i0 = np.insert(i1+1, 0, 0)
    i1 = np.append(i1, -1)
    paGood = np.dstack((pa[i0], pa[i1])).round(1).reshape(-1, 2).tolist()

    # bad PAs (complement of the good PAs)
    paBad = []
    if paGood[0][0] > 0:
        paBad.append([0., paGood[0][0]])
    for i in range(1, len(paGood)):
        paBad.append([paGood[i-1][1], paGood[i][0]])
    if paGood[-1][1] < 360.:
        paBad.append([paGood[-1][1], 360.])

    # Make a figure
    if fig is None or fig == True:
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
    if paNom[i+1] < 10:
        i += 1
    paMasked = np.insert(paMasked, i, np.nan)
    gdMasked = np.insert(gdMasked, i, gdMasked[i])

    i = np.argmax(paMin)
    goUp = paMin[i-2] < paMin[i-1]  # PA going up at wrap point?

    # Top part
    i0_top = 0 if goUp else i
    i1_top = i if goUp else paMin.size-1
    paMaxTmp = np.copy(paMax)
    paMaxTmp[np.where(paMin > paMax)[0]] = 360

    # Bottom part
    i = np.argmin(paMax)
    i0_bot = i if goUp else 0
    i1_bot = paMin.size-1 if goUp else i
    paMinTmp = np.copy(paMin)
    paMinTmp[np.where(paMin > paMax)[0]] = 0

    # Convert datetime to a number for Bokeh
    gdMaskednum = [datetime.date(2019, 6, 1)+datetime.timedelta(days=n)
                   for n, d in enumerate(gdMasked)]
    color = 'green'

    # Draw the curve and error
    fig.line(gdMaskednum, paMasked, legend='cutoff', line_color=color)

    # Top
    terr_y = np.concatenate([paMin[i0_top:i1_top+1],
                            paMaxTmp[i0_top:i1_top+1][::-1]])
    terr_x = np.concatenate([gdMaskednum[i0_top:i1_top+1],
                            gdMaskednum[i0_top:i1_top+1][::-1]])
    fig.patch(terr_x, terr_y, color=color, fill_alpha=0.2, line_alpha=0)


    # Bottom
    berr_y = np.concatenate([paMinTmp[i0_bot:i1_bot+1],
                            paMax[i0_bot:i1_bot+1][::-1]])
    berr_x = np.concatenate([gdMaskednum[i0_bot:i1_bot+1],
                            gdMaskednum[i0_bot:i1_bot+1][::-1]])
    fig.patch(berr_x, berr_y, color='red', fill_alpha=0.2, line_alpha=0)
    from bokeh.plotting import show
    show(fig)


    # Plot formatting
    fig.xaxis.axis_label = 'Date'
    fig.yaxis.axis_label = 'Position Angle (degrees)'

    #return paGood, paBad, gd, fig
    return terr_x, terr_y, berr_x, berr_y

def using_gtvt(ra, dec, instrument, targetName=None, ephFileName=None, FIG='bokeh'):
    """Check the visibility at a range of position angles

    Parameters
    ----------
    ra: float
        The RA of the target
    dec: float
        The Dec of the target
    instrument: str
        Name of the instrument. Can either be (case-sensitive):
        NIRISS, NIRCam, MIRI, FGS, NIRSpec
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
       The gregorian date.
    fig : bokeh.plotting.figure object
        The plotted figure.

    """

    from bokeh.plotting import save, output_file, show
    from bokeh.io import reset_output

    output_file('/Users/jmedina/Desktop/visib_test4.html')
    # getting calculations from GTVT (General Target Visibility Tool)
    tab = get_table(ra, dec)

    gd = tab['Date']
    paMin = tab[str(instrument)+' min']
    paMax = tab[str(instrument)+' max']
    paNom = tab['V3PA']

    paMin, paMax, gd = np.asarray(paMin), np.asarray(paMax), np.asarray(gd)
    paMinnan, paMaxnan, gdnan = [], [], []


    for i, j, k in zip(paMin, paMax, gd):
        maskMin, maskMax = np.isnan(i), np.isnan(j)
        if (maskMin==False)&(maskMax==False):
            paMinnan.append(i)
            paMaxnan.append(j)
            gdnan.append(k)

    paMinnan = np.asarray(paMinnan)
    paMaxnan = np.asarray(paMaxnan)
    gdnan = np.asarray(gdnan)




    # Convert datetime to a number for Bokeh
    #gdMaskednum = [datetime.date(2019, 6, 1)+datetime.timedelta(days=n) for n,d in enumerate(gd)]
    color = 'green'
    from bokeh.models import HoverTool
    from bokeh.plotting import ColumnDataSource


    SOURCE = ColumnDataSource(data=dict(pamin=paMin, panom=paNom, pamax=paMax, date=gd))
    TOOLS = 'pan, wheel_zoom, reset, save'

    if FIG=='bokeh':
        fig = figure(tools=TOOLS, plot_width=800, plot_height=400, x_axis_type='datetime')

    # Draw the curve and error
    #fig.line(gd, paMax, legend='min PA', line_color='blue', source=SOURCE)
    min = fig.line('date', 'pamin', color='red', legend='min PA', source=SOURCE)
    nom = fig.line('date', 'panom', color='black', legend='nom PA', source=SOURCE)
    max = fig.line('date', 'pamax', color='blue', legend='max PA', source=SOURCE)

    fig.add_tools(HoverTool(renderers=[min], tooltips=[('Date','@date{%F}'), ('Minimum Angle', '@pamin')], formatters={'date':'datetime'}))
    fig.add_tools(HoverTool(renderers=[nom], tooltips=[('Date','@date{%F}'), ('Nominal Angle', '@panom')], formatters={'date':'datetime'}))
    fig.add_tools(HoverTool(renderers=[max], tooltips=[('Date','@date{%F}'), ('Maximum Angle', '@pamax')], formatters={'date':'datetime'}))

    #x = np.append(gdnan, gdnan[::-1])
    #y = np.append(paMinnan, paMaxnan[::-1])
    #fig.patch(x,y, color='purple', fill_alpha=0.2, line_alpha=0)
    #fig.patch(gd, paMin)
    #fig.patches([gd,gd],[paMax,paMin])



    """
    i = np.argmax(paMinnan)
    goUp = paMinnan[i-2] < paMinnan[i-1]  # PA going up at wrap point?

    # Top part
    i0_top = 0 if goUp else i
    i1_top = i if goUp else paMinnan.size-1
    paMaxTmp = np.copy(paMaxnan)
    paMaxTmp[np.where(paMinnan > paMaxnan)[0]] = 360

    # Bottom part
    i = np.argmin(paMaxnan)
    i0_bot = i if goUp else 0
    i1_bot = paMinnan.size-1 if goUp else i
    paMinTmp = np.copy(paMinnan)
    paMinTmp[np.where(paMinnan > paMaxnan)[0]] = 0

    # Top

    err_y = np.concatenate([paMinnan[i0_top:i1_top+1],
                            paMaxTmp[i0_top:i1_top+1][::-1]])
    err_x = np.concatenate([gdnan[i0_top:i1_top+1],
                            gdnan[i0_top:i1_top+1][::-1]])
    for i in zip(err_y, err_x):
        print(i)
    fig.patch(err_x, err_y, color=color, fill_alpha=0.2, line_alpha=0)

    # Bottom
    err_y = np.concatenate([paMinTmp[i0_bot:i1_bot+1],
                            paMaxnan[i0_bot:i1_bot+1][::-1]])
    err_x = np.concatenate([gdnan[i0_bot:i1_bot+1],
                            gdnan[i0_bot:i1_bot+1][::-1]])
    print('bottom')
    for i in zip(err_y, err_x):
        print(i)
    fig.patch(err_x, err_y, color=color, fill_alpha=0.2, line_alpha=0)
    """
    """
    # Do all figure calculations
    #iBad, = np.where(vis == False)
    paMasked = np.copy(paNom)
    #paMasked[iBad] = np.nan


    #i = np.argmax(paNom)
    #if paNom[i+1] < 10:
    #    i += 1
    #paMasked = np.insert(paMasked, i, np.nan)
    #gdMasked = np.insert(gdMasked, i, gdMasked[i])

    # future jenny: i think these 2 lines are the issue (w/ the ribbon plotting)
    i = np.argmax(paMinnan)
    goUp = paMinnan[i-2] < paMinnan[i-1]  # PA going up at wrap point?

    # Top part
    i0_top = 0 if goUp else i               # indexing
    i1_top = i if goUp else paMinnan.size-1 # indexing
    paMaxTmp = np.copy(paMaxnan)
    paMaxTmp[np.where(paMinnan > paMaxnan)[0]] = 360 # <- what is the purpose of this

    # Bottom part
    i = np.argmin(paMaxnan)
    i0_bot = i if goUp else 0
    i1_bot = paMinnan.size-1 if goUp else i
    paMinTmp = np.copy(paMinnan)
    paMinTmp[np.where(paMinnan > paMaxnan)[0]] = 0

    # Convert datetime to a number for Bokeh
    #gdMaskednum = [datetime.date(2019, 6, 1)+datetime.timedelta(days=n)
    #               for n, d in enumerate(gdMasked)]
    color = 'green'

    # Draw the curve and error
    #fig.line(gd, paMasked, legend='cutoff', line_color=color)

    # Top
    err_y = np.concatenate([paMinnan[i0_top:i1_top+1],
                            paMaxTmp[i0_top:i1_top+1][::-1]])
    err_x = np.concatenate([gdnan[i0_top:i1_top+1],
                            gdnan[i0_top:i1_top+1][::-1]])

    fig.patch(err_x, err_y, color=color, fill_alpha=0.2, line_alpha=0)

    print('top')
    for i in zip(err_x, err_y):
        print(i)
    # Bottom
    err_y = np.concatenate([paMinTmp[i0_bot:i1_bot+1],
                            paMaxnan[i0_bot:i1_bot+1][::-1]])
    err_x = np.concatenate([gdnan[i0_bot:i1_bot+1],
                            gdnan[i0_bot:i1_bot+1][::-1]])


    fig.patch(err_x, err_y, color=color, fill_alpha=0.2, line_alpha=0)
    print('bottom')
    for j in zip(err_x, err_y):
        print(j)
    """
    # Plot formatting
    fig.xaxis.axis_label = 'Date'
    fig.yaxis.axis_label = 'Position Angle (degrees)'


    save(fig)

    print('this is a test2s')
    return paMin, paMax, gd, paNom

import os
import pkg_resources
import sys

from astropy.io import fits
from bokeh.layouts import gridplot
from bokeh.plotting import figure
from bokeh.models import Range1d, LinearColorMapper, Label
from bokeh.models.widgets import Panel, Tabs
from bokeh.palettes import PuBu
import numpy as np

from . import visibilityPA as vpa
from ..utils import fill_between


EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
if not EXOCTK_DATA:
    print(
        'WARNING: The $EXOCTK_DATA environment variable is not set. Contamination overlap will not work. Please set the '
        'value of this variable to point to the location of the exoctk_data '
        'download folder.  Users may retreive this folder by clicking the '
        '"ExoCTK Data Download" button on the ExoCTK website, or by using '
        'the exoctk.utils.download_exoctk_data() function.')
    TRACES_PATH = None
else:
    TRACES_PATH = os.path.join(EXOCTK_DATA, 'exoctk_contam', 'traces')

disp_nircam = 0.001  # microns
lam0_nircam322w2 = 2.369
lam1_nircam322w2 = 4.417
lam0_nircam444w = 3.063
lam1_nircam444w = 5.111


def nirissContam(cube, paRange=[0, 360]):
    """ Generates the contamination figure that will be plotted on the website
    for NIRISS SOSS.
    """
    # Get data from FITS file
    if isinstance(cube, str):
        hdu = fits.open(cubeName)
        cube = hdu[0].data
        hdu.close()

    # Pull out the target trace and cube of neighbor traces
    trace1 = cube[0, :, :]
    trace2 = cube[1, :, :]
    cube = cube[2:, :, :]

    plotPAmin, plotPAmax = paRange

    # Start calculations
    if not TRACES_PATH:
        return None
    lam_file = os.path.join(TRACES_PATH, 'NIRISS', 'lambda_order1-2.txt')
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    nPA = cube.shape[0]
    rows = cube.shape[1]
    cols = cube.shape[2]
    print('cols ', cols)
    dPA = 360 // nPA
    PA = np.arange(nPA) * dPA

    contamO1 = np.zeros([rows, nPA])
    contamO2 = np.zeros([rows, nPA])

    low_lim_col = 20
    high_lim_col = 41

    for row in np.arange(rows):
        # Contamination for order 1 of target trace
        i = np.argmax(trace1[row, :])
        tr = trace1[row, i - low_lim_col:i + high_lim_col]
        w = tr / np.sum(tr**2)
        ww = np.tile(w, nPA).reshape([nPA, tr.size])
        contamO1[row, :] = np.sum(
            cube[:, row, i - low_lim_col:i + high_lim_col] * ww, axis=1)

        # Contamination for order 2 of target trace
        if lamO2[row] < 0.6:
            continue
        i = np.argmax(trace2[row, :])
        tr = trace2[row, i - 20:i + 41]
        w = tr / np.sum(tr**2)
        ww = np.tile(w, nPA).reshape([nPA, tr.size])
        contamO2[row, :] = np.sum(cube[:, row, i - 20:i + 41] * ww, axis=1)

    return contamO1, contamO2


def nircamContam(cube, instrument, paRange=[0, 360]):
    """ Generates the contamination figure that will be plotted on the website
    for NIRCam Grism Time Series mode.

    PARAMETERS
    ----------
    cube : arr or str
        A 3D array of the simulated field at every Aperture Position Angle (APA).
        The shape of the cube is (361, subY, subX).
        or
        The name of an HDU .fits file sthat has the cube.

    instrument : str
        The name of the instrument + what filter is being used. For NIRCam the
        options are: 'NIRCam F322W2', 'NIRCam F444W'

    RETURNS
    -------
    bokeh plot
    """
    # Get data from FITS file
    if isinstance(cube, str):
        hdu = fits.open(cubeName)
        cube = hdu[0].data
        hdu.close()

    # Pull out the target trace and cube of neighbor traces
    targ = cube[0, :, :]  # target star order 1 trace
    # neighbor star order 1 and 2 traces in all the angles
    cube = cube[1:, :, :]

    # Remove background values < 1 as it can blow up contamination
    targ = np.where(targ < 1, 0, targ)

    PAmin, PAmax = paRange[0], paRange[1]
    PArange = np.arange(PAmin, PAmax, 1)

    nPA, rows, cols = cube.shape[0], cube.shape[1], cube.shape[2]

    contamO1 = np.zeros([nPA, cols])

    # the width of the trace (in Y-direction for NIRCam GTS)
    peak = targ.max()
    low_lim_row = np.where(targ > 0.0001 * peak)[0].min()
    high_lim_row = np.where(targ > 0.0001 * peak)[0].max()

    # the length of the trace (in X-direction for NIRCam GTS)
    targ_trace_start = np.where(targ > 0.0001 * peak)[1].min()
    targ_trace_stop = np.where(targ > 0.0001 * peak)[1].max()

    # Begin contam calculation at each channel (column) X
    for X in np.arange(cols):
        if (X < targ_trace_start) or (X > targ_trace_stop):
            continue

        peakY = np.argmax(targ[:, X])
        TOP, BOT = peakY + high_lim_row, peakY - low_lim_row

        tr = targ[BOT:TOP, X]

        # calculate weights
        wt = tr / np.sum(tr**2)
        ww = np.tile(wt, nPA).reshape([nPA, tr.size])

        contamO1[:, X] = np.sum(cube[:, BOT:TOP, X] * ww, axis=1)

    contamO1 = contamO1[:, targ_trace_start:targ_trace_stop]
    return contamO1


def miriContam(cube, paRange=[0, 360]):
    """ Generates the contamination figure that will be plotted on the website
    for MIRI LRS.
    """
    # Get data from FITS file
    if isinstance(cube, str):
        hdu = fits.open(cubeName)
        cube = hdu[0].data
        hdu.close()

    # Pull out the target trace and cube of neighbor traces
    targ = cube[0, :, :]  # target star order 1 trace
    # neighbor star order 1 and 2 traces in all the angles
    cube = cube[1:, :, :]

    # Remove background values < 1 as it can blow up contamination
    targ = np.where(targ < 1, 0, targ)

    PAmin, PAmax = paRange[0], paRange[1]
    PArange = np.arange(PAmin, PAmax, 1)

    nPA, rows, cols = cube.shape[0], cube.shape[1], cube.shape[2]

    contamO1 = np.zeros([rows, nPA])

    # the width of the trace (in Y-direction for NIRCam GTS)
    peak = targ.max()

    low_lim_col = np.where(targ > 0.0001 * peak)[1].min()
    high_lim_col = np.where(targ > 0.0001 * peak)[1].max()

    # the length of the trace (in X-direction for NIRCam GTS)
    targ_trace_start = np.where(targ > 0.0001 * peak)[0].min()
    targ_trace_stop = np.where(targ > 0.0001 * peak)[0].max()
    # Begin contam calculation at each channel (row) Y
    for Y in np.arange(rows):
        if (Y < targ_trace_start) or (Y > targ_trace_stop):
            continue

        peakX = np.argmax(targ[Y, :])
        LEFT, RIGHT = peakX - low_lim_col, peakX + high_lim_col

        tr = targ[Y, LEFT:RIGHT]

        # calculate weights
        wt = tr / np.sum(tr**2)
        ww = np.tile(wt, nPA).reshape([nPA, tr.size])

        contamO1[Y, :] = np.sum(cube[:, Y, LEFT:RIGHT] * wt,
                                where=~np.isnan(cube[:, Y, LEFT:RIGHT] * wt),
                                axis=1)

        #target = np.sum(cube[0, Y, LEFT:RIGHT], axis=0)

        # contamO1[Y, :] = np.sum(cube[:, Y, LEFT:RIGHT]*ww,
        #                        where=~np.isnan(cube[:, Y, LEFT:RIGHT]),
        #                        axis=1)#/target
    contamO1 = contamO1[targ_trace_start:targ_trace_stop, :]
    return contamO1


def contam(cube, instrument, targetName='noName', paRange=[0, 360],
           badPAs=np.asarray([]), tmpDir="", fig='', to_html=True):

    lam_file = os.path.join(TRACES_PATH, 'NIRISS', 'lambda_order1-2.txt')
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    nPA, rows, cols = cube.shape[0], cube.shape[1], cube.shape[2]

    PAmin, PAmax = paRange[0], paRange[1]
    PA = np.arange(PAmin, PAmax, 1)

    # Generate the contam figure
    if instrument == 'NIRISS':
        contamO1, contamO2 = nirissContam(cube)
    elif (instrument == 'NIRCam F322W2') or (instrument == 'NIRCam F444W'):
        contamO1 = nircamContam(cube, instrument)
    elif instrument == 'MIRI':
        contamO1 = miriContam(cube)

    TOOLS = 'pan, box_zoom, crosshair, reset, hover'

    y = np.array([0., 0.])
    y1 = 0.07
    y2 = 0.12
    y3 = 0.17
    y4 = 0.23
    bad_PA_color = '#dddddd'
    bad_PA_alpha = 0.7
    dPA = 1

    # Order 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Contam plot
    if instrument == 'NIRISS':
        xlim0 = lamO1.min()
        xlim1 = lamO1.max()
    elif instrument == 'NIRCam F322W2':
        xlim0 = lam0_nircam322w2
        xlim1 = lam1_nircam322w2
    elif instrument == 'NIRCam F444W':
        xlim0 = lam0_nircam444w
        xlim1 = lam1_nircam444w
    elif instrument == 'MIRI':
        xlim0 = 5
        xlim1 = 12

    ylim0 = PAmin - 0.5
    ylim1 = PAmax + 0.5
    color_mapper = LinearColorMapper(palette=PuBu[8][::-1][2:],
                                     low=-4, high=1)
    color_mapper.low_color = 'white'
    color_mapper.high_color = 'black'

    orders = 'Orders 1 & 2' if instrument == 'NIRCam' else 'Order 1'
    s2 = figure(
        tools=TOOLS, width=500, height=500, title='{} {} Contamination with {}'.format(
            orders, targetName, instrument), x_range=Range1d(
            xlim0, xlim1), y_range=Range1d(
                ylim0, ylim1))

    contamO1 = contamO1 if 'NIRCam' in instrument else contamO1.T
    contamO1 = np.fliplr(contamO1) if (
        instrument == 'MIRI') or (
        instrument == 'NIRCam F322W2') else contamO1
    fig_data = np.log10(np.clip(contamO1, 1.e-10, 1.)
                        )  # [:, :361] # might this
    # index have somethig to
    # do w the choppiness
    # of o1 in all instruments

    X = xlim1 if (instrument == 'MIRI') or (
        instrument == 'NIRCam F322W2') else xlim0
    DW = xlim0 - xlim1 if (instrument == 'MIRI') or (instrument ==
                                                     'NIRCam F322W2') else xlim1 - xlim0

    # Begin plotting ~~~~~~~~~~~~~~~~~~~~~~~~

    s2.image([fig_data], x=xlim0, y=ylim0, dw=xlim1 - xlim0, dh=ylim1 - ylim0,
             color_mapper=color_mapper)
    s2.xaxis.axis_label = 'Wavelength (um)'
    if instrument != 'NIRISS':
        s2.yaxis.axis_label = 'Aperture Position Angle (degrees)'

    # Add bad PAs
    bad_PA_color = '#555555'
    bad_PA_alpha = 0.6
    # for ybad0, ybad1 in badPA:
    if len(badPAs) > 0:

        tops, bottoms, lefts, rights = [], [], [], []
        for idx in range(0, len(badPAs)):
            PAgroup = badPAs[idx]
            top_idx = np.max(PAgroup)
            bot_idx = np.min(PAgroup)

            tops.append(top_idx)
            bottoms.append(bot_idx)
            lefts.append(xlim0)
            rights.append(xlim1)

        s2.quad(top=tops, bottom=bottoms,
                left=lefts, right=rights,
                color=bad_PA_color, alpha=bad_PA_alpha)

    # Line plot
    #ax = 1 if 'NIRCam' in instrument else 0
    channels = cols if 'NIRCam' in instrument else rows
    s3 = figure(tools=TOOLS, width=150, height=500,
                x_range=Range1d(0, 100), y_range=s2.y_range, title=None)

    try:
        s3.line(100 * np.sum(contamO1 >= 0.001, axis=1) / channels, PA - dPA / 2,
                line_color='blue', legend_label='> 0.001')
        s3.line(100 * np.sum(contamO1 >= 0.01, axis=1) / channels, PA - dPA / 2,
                line_color='green', legend_label='> 0.01')
    except AttributeError:
        s3.line(100 * np.sum(contamO1 >= 0.001, axis=1) / channels, PA - dPA / 2,
                line_color='blue', legend='> 0.001')
        s3.line(100 * np.sum(contamO1 >= 0.01, axis=1) / channels, PA - dPA / 2,
                line_color='green', legend='> 0.01')

    s3.xaxis.axis_label = '% channels contam.'
    s3.yaxis.major_label_text_font_size = '0pt'

    # ~~~~~~ Order 2 ~~~~~~
    # Contam plot
    if instrument == 'NIRISS':
        xlim0 = lamO2.min()
        xlim1 = lamO2.max()
        ylim0 = PA.min() - 0.5 * dPA
        ylim1 = PA.max() + 0.5 * dPA
        xlim0 = 0.614
        s5 = figure(
            tools=TOOLS,
            width=500,
            height=500,
            title='Order 2 {} Contamination with {}'.format(
                targetName,
                instrument),
            x_range=Range1d(
                xlim0,
                xlim1),
            y_range=s2.y_range)
        fig_data = np.log10(np.clip(contamO2.T, 1.e-10, 1.))[:, 300:]
        s5.image(
            [fig_data],
            x=xlim0,
            y=ylim0,
            dw=xlim1 - xlim0,
            dh=ylim1 - ylim0,
            color_mapper=color_mapper)
        #s5.yaxis.major_label_text_font_size = '0pt'
        s5.xaxis.axis_label = 'Wavelength (um)'
        s5.yaxis.axis_label = 'Aperture Position Angle (degrees)'

        if len(badPAs) > 0:

            tops, bottoms, lefts, rights = [], [], [], []
            for idx in range(0, len(badPAs)):
                PAgroup = badPAs[idx]
                top_idx = np.max(PAgroup)
                bot_idx = np.min(PAgroup)

                tops.append(top_idx)
                bottoms.append(bot_idx)
                lefts.append(xlim0)
                rights.append(xlim1)

            s5.quad(top=tops, bottom=bottoms,
                    left=lefts, right=rights,
                    color=bad_PA_color, alpha=bad_PA_alpha)

        # Line plot
        s6 = figure(tools=TOOLS, width=150, height=500, y_range=s2.y_range,
                    x_range=Range1d(100, 0), title=None)

        try:
            s6.line(100 * np.sum(contamO2 >= 0.001, axis=0) / rows, PA - dPA / 2,
                    line_color='blue', legend_label='> 0.001')
            s6.line(100 * np.sum(contamO2 >= 0.01, axis=0) / rows, PA - dPA / 2,
                    line_color='green', legend_label='> 0.01')
        except AttributeError:
            s6.line(100 * np.sum(contamO2 >= 0.001, axis=0) / rows, PA - dPA / 2,
                    line_color='blue', legend='> 0.001')
            s6.line(100 * np.sum(contamO2 >= 0.01, axis=0) / rows, PA - dPA / 2,
                    line_color='green', legend='> 0.01')

        s6.xaxis.axis_label = '% channels contam.'
        s6.yaxis.major_label_text_font_size = '0pt'

    if len(badPAs) > 0:

        tops, bottoms, lefts, rights = [], [], [], []
        for idx in range(0, len(badPAs)):
            PAgroup = badPAs[idx]
            top_idx = np.max(PAgroup)
            bot_idx = np.min(PAgroup)

            tops.append(top_idx)
            bottoms.append(bot_idx)
            lefts.append(0)
            rights.append(100)

        s3.quad(top=tops, bottom=bottoms,
                left=lefts, right=rights,
                color=bad_PA_color, alpha=bad_PA_alpha)
        if instrument == 'NIRISS':
            s6.quad(top=tops, bottom=bottoms,
                    left=rights, right=lefts,
                    color=bad_PA_color, alpha=bad_PA_alpha)

    # ~~~~~~ Plotting ~~~~~~

    if instrument != 'NIRISS':
        fig = gridplot(children=[[s2, s3]])
    else:
        fig = gridplot(children=[[s6, s5, s2, s3]])

    return fig  # , contamO1


if __name__ == "__main__":
    # arguments RA & DEC, conversion to radians
    argv = sys.argv

    ra = argv[1]
    dec = argv[2]
    cubeNameSuf = argv[3]

    pamin = 0 if len(argv) < 5 else int(argv[4])
    pamax = 360 if len(argv) < 6 else int(argv[5])

    cubeName = argv[6]
    targetName = None if len(argv) < 8 else argv[7]
    save = False if len(argv) < 8 else True  # if name provided -> save
    tmpDir = "." if len(argv) < 9 else argv[8]
    os.makedirs(tmpDir, exist_ok=True)

    goodPA, badPA, _ = vpa.checkVisPA(ra, dec, targetName)

    contam(cubeName, targetName=targetName, paRange=[pamin, pamax],
           badPA=badPA, tmpDir=tmpDir)

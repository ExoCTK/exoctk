import os
import pkg_resources
import sys

from astropy.io import fits
from bokeh.layouts import gridplot
from bokeh.plotting import figure
from bokeh.models import Range1d, LinearColorMapper, Label
from bokeh.models.widgets import Panel, Tabs
from bokeh.palettes import inferno
import numpy as np

from . import visibilityPA as vpa

TRACES_PATH = os.path.join(os.environ.get('EXOCTK_DATA'),  'exoctk_contam', 'traces')

disp_nircam = 0.001 # microns
lam0_nircam322w2 = 2.369
lam1_nircam322w2 = 4.417
lam0_nircam444w = 3.063
lam1_nircam444w = 5.111

def contam(cube, instrument, targetName='noName', paRange=[0, 360],
           badPAs=np.asarray([]), tmpDir="", fig='', to_html=True):

    # Get data from FITS file
    if isinstance(cube, str):
        hdu = fits.open(cubeName)
        cube = hdu[0].data
        hdu.close()

    if instrument != 'NIRISS':
        trace1 = cube[0, :, :] # target star order 1 trace
        cube = cube[1:, :, :] # neighbor star order 1 and 2 traces in all the angles

    elif instrument=='NIRISS': # only NIRISS has both orders accounted for
        trace1 = cube[0, :, :]
        trace2 = cube[1, :, :]
        cube = cube[2:, :, :]

    plotPAmin, plotPAmax = paRange

    # start calculations
    lam_file = os.path.join(TRACES_PATH, 'NIRISS', 'lambda_order1-2.txt')
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    nPA = cube.shape[0]
    rows = cube.shape[1]
    cols = cube.shape[2]
    print('cols ', cols)
    dPA = 360//nPA
    PA = np.arange(nPA)*dPA

    contamO1 = np.zeros([rows, nPA])
    if instrument=='NIRISS':
        contamO2 = np.zeros([rows, nPA])

    for row in np.arange(rows):
        i = np.argmax(trace1[row, :])
        tr = trace1[row, i-20:i+41]
        #tr = trace1[row, i-100:i+100]
        w = tr/np.sum(tr**2)
        ww = np.tile(w, nPA).reshape([nPA, tr.size])

        contamO1[row, :] = np.sum(cube[:, row, i-20:i+41]*ww, axis=1)

        if instrument=='NIRISS':
            if lamO2[row] < 0.6:
                continue
            i = np.argmax(trace2[row, :])
            tr = trace2[row, i-20:i+41]
            w = tr/np.sum(tr**2)
            ww = np.tile(w, nPA).reshape([nPA, tr.size])
            contamO2[row, :] = np.sum(cube[:, row, i-20:i+41]*ww, axis=1)


    TOOLS = 'pan, box_zoom, crosshair, reset, hover'

    y = np.array([0., 0.])
    y1 = 0.07
    y2 = 0.12
    y3 = 0.17
    y4 = 0.23
    bad_PA_color = '#dddddd'
    bad_PA_alpha = 0.7

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

    ylim0 = PA.min()-0.5*dPA
    ylim1 = PA.max()+0.5*dPA
    color_mapper = LinearColorMapper(palette=inferno(8)[::-1],
                                     low=-4, high=1)
    color_mapper.low_color = 'white'
    color_mapper.high_color = 'black'
    s2 = figure(tools=TOOLS, width=500, height=500,
                title='Order 1 {} Contamination with {}'.format(targetName, instrument),
                x_range=Range1d(xlim0, xlim1),
                y_range=Range1d(ylim0, ylim1))
    fig_data = np.log10(np.clip(contamO1.T, 1.e-10, 1.))
    s2.image([fig_data], x=xlim0, y=ylim0, dw=xlim1-xlim0, dh=ylim1-ylim0,
             color_mapper=color_mapper)
    s2.xaxis.axis_label = 'Wavelength (um)'
    if instrument != 'NIRISS':
        s2.yaxis.axis_label = 'Aperture Position Angle (degrees)'

    # Add bad PAs
    bad_PA_color = '#555555'
    bad_PA_alpha = 0.6
    #for ybad0, ybad1 in badPA:
    if len(badPAs)>0:

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
    s3 = figure(tools=TOOLS, width=150, height=500,
                x_range=Range1d(0, 100), y_range=s2.y_range, title=None)
    s3.line(100*np.sum(contamO1 >= 0.001, axis=0)/rows, PA-dPA/2,
            line_color='blue', legend='> 0.001')
    s3.line(100*np.sum(contamO1 >= 0.01, axis=0)/rows, PA-dPA/2,
            line_color='green', legend='> 0.01')
    s3.xaxis.axis_label = '% channels contam.'
    s3.yaxis.major_label_text_font_size = '0pt'

    # ~~~~~~ Order 2 ~~~~~~
    # Contam plot
    if instrument=='NIRISS':
        xlim0 = lamO2.min()
        xlim1 = lamO2.max()
        ylim0 = PA.min()-0.5*dPA
        ylim1 = PA.max()+0.5*dPA
        xlim0 = 0.614
        s5 = figure(tools=TOOLS, width=500, height=500,
                    title='Order 2 {} Contamination with {}'.format(targetName, instrument),
                    x_range=Range1d(xlim0, xlim1), y_range=s2.y_range)
        fig_data = np.log10(np.clip(contamO2.T, 1.e-10, 1.))[:, 300:]
        s5.image([fig_data], x=xlim0, y=ylim0, dw=xlim1-xlim0, dh=ylim1-ylim0,
                 color_mapper=color_mapper)
        #s5.yaxis.major_label_text_font_size = '0pt'
        s5.xaxis.axis_label = 'Wavelength (um)'
        s5.yaxis.axis_label = 'Aperture Position Angle (degrees)'

        if len(badPAs)>0:

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
        s6.line(100*np.sum(contamO2 >= 0.001, axis=0)/rows, PA-dPA/2,
                line_color='blue', legend='> 0.001')
        s6.line(100*np.sum(contamO2 >= 0.01, axis=0)/rows, PA-dPA/2,
                line_color='green', legend='> 0.01')
        s6.xaxis.axis_label = '% channels contam.'
        s6.yaxis.major_label_text_font_size = '0pt'

    if len(badPAs)>0:

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
        if instrument=='NIRISS':
            s6.quad(top=tops, bottom=bottoms,
                    left=rights, right=lefts,
                     color=bad_PA_color, alpha=bad_PA_alpha)

    # ~~~~~~ Plotting ~~~~~~

    if instrument!='NIRISS':
        fig = gridplot(children=[[s2, s3]])
    else:
        fig = gridplot(children=[[s6, s5, s2, s3]])

    return fig#, contamO1


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

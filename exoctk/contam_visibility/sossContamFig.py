import os
import pkg_resources
import sys

from astropy.io import fits
from bokeh.layouts import gridplot
from bokeh.plotting import figure
from bokeh.models import Range1d, LinearColorMapper, Label
from bokeh.palettes import inferno
import numpy as np

from . import visibilityPA as vpa


def contam(cube, targetName='noName', paRange=[0, 360], badPA=[], tmpDir="",
           fig='', to_html=True):
    """
    Generate the contamination plot.

    Parameters
    ----------
    cube: array-like, str
        The data cube or FITS filename containing the data.
    targetName: str
        The name of the target.
    paRange: sequence
        The position angle range to consider.
    badPA: sequence
        Position angles to exclude.
    tmpDir: str
        A directory to write the files to.
    fig: matplotlib.figure, bokeh.figure
        A figure to add the plots to.
    to_html: bool
        Return the image as bytes for HTML.

    Returns
    -------
    fig : matplotlib.pyplot or bokeh object
        The populated matplotlib or bokeh plot.
    """
    # Get data from FITS file
    if isinstance(cube, str):
        # hdu = fits.open('cube_'+target+'.fits')
        hdu = fits.open(cubeName)
        cube = hdu[0].data
        hdu.close()

    trace2dO1 = cube[0, :, :]  # order 1
    trace2dO2 = cube[1, :, :]  # order 2
    cube = cube[2:, :, :]  # all the angles

    plotPAmin, plotPAmax = paRange

    # start calculations
    loc = 'data/contam_visibility/lambda_order1-2.txt'
    lam_file = pkg_resources.resource_filename('exoctk', loc)
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    ny = trace2dO1.shape[0]
    nPA = cube.shape[0]
    dPA = 360//nPA
    PA = np.arange(nPA)*dPA

    contamO1 = np.zeros([ny, nPA])
    contamO2 = np.zeros([ny, nPA])
    for y in np.arange(ny):
        i = np.argmax(trace2dO1[y, :])
        tr = trace2dO1[y, i-20:i+41]
        w = tr/np.sum(tr**2)
        ww = np.tile(w, nPA).reshape([nPA, tr.size])
        contamO1[y, :] = np.sum(cube[:, y, i-20:i+41]*ww, axis=1)

        if lamO2[y] < 0.6:
            continue
        i = np.argmax(trace2dO2[y, :])
        tr = trace2dO2[y, i-20:i+41]
        w = tr/np.sum(tr**2)
        ww = np.tile(w, nPA).reshape([nPA, tr.size])
        contamO2[y, :] = np.sum(cube[:, y, i-20:i+41]*ww, axis=1)

    # Otherwise, it's a Bokeh plot
    if fig:

        TOOLS = 'pan, wheel_zoom, box_zoom, crosshair, reset, hover'

        y = np.array([0., 0.])
        y1 = 0.07
        y2 = 0.12
        y3 = 0.17
        y4 = 0.23
        bad_PA_color = '#dddddd'
        bad_PA_alpha = 0.7

        # Order 1

        # Contam plot
        xlim0 = lamO1.min()
        xlim1 = lamO1.max()
        ylim0 = PA.min()-0.5*dPA
        ylim1 = PA.max()+0.5*dPA
        color_mapper = LinearColorMapper(palette=inferno(8)[::-1],
                                         low=-4, high=1)
        color_mapper.low_color = 'white'
        color_mapper.high_color = 'black'
        s2 = figure(tools=TOOLS, width=500, height=500, title=None,
                    x_range=Range1d(xlim0, xlim1),
                    y_range=Range1d(ylim0, ylim1))
        fig_data = np.log10(np.clip(contamO1.T, 1.e-10, 1.))
        s2.image([fig_data], x=xlim0, y=ylim0, dw=xlim1-xlim0, dh=ylim1-ylim0,
                 color_mapper=color_mapper)
        s2.xaxis.axis_label = 'Wavelength (um)'
        s2.yaxis.axis_label = 'Position Angle (degrees)'

        # Line plot
        s3 = figure(tools=TOOLS, width=150, height=500,
                    x_range=Range1d(0, 100), y_range=s2.y_range, title=None)
        s3.line(100*np.sum(contamO1 >= 0.001, axis=0)/ny, PA-dPA/2,
                line_color='blue', legend='> 0.001')
        s3.line(100*np.sum(contamO1 >= 0.01, axis=0)/ny, PA-dPA/2,
                line_color='green', legend='> 0.01')
        s3.xaxis.axis_label = '% channels contam.'
        s3.yaxis.major_label_text_font_size = '0pt'

        # Add bad PAs
        for ybad0, ybad1 in badPA:
            s2.patch([xlim0, xlim1, xlim1, xlim0],
                     [ybad1, ybad1, ybad0, ybad0],
                     color=bad_PA_color, alpha=bad_PA_alpha)
            s3.patch([0, 100, 100, 0], [ybad1, ybad1, ybad0, ybad0],
                     color=bad_PA_color, alpha=bad_PA_alpha, legend='Bad PA')

        # Line list
        s1 = figure(tools=TOOLS, width=500, plot_height=100,
                    x_range=s2.x_range, title=None)

        l = np.array([0.89, 0.99])
        s1.line(l, y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                      text='H2O', render_mode='css', text_font_size='8pt'))

        l = np.array([1.09, 1.2])
        s1.line(l, y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                      text='H2O', render_mode='css', text_font_size='8pt'))

        l = np.array([1.1, 1.24])
        s1.line(l, y+y2, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data',
                      text='CH4', render_mode='css', text_font_size='8pt'))

        l = np.array([1.3, 1.51])
        s1.line(l, y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                      text='H2O', render_mode='css', text_font_size='8pt'))

        l = np.array([1.6, 1.8])
        s1.line(l, y+y2, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data',
                      text='CH4', render_mode='css', text_font_size='8pt'))

        l = np.array([1.75, 2.05])
        s1.line(l, y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                      text='H2O', render_mode='css', text_font_size='8pt'))

        l = np.array([2.3, lamO1.max()])
        s1.line(l, y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                      text='H2O', render_mode='css', text_font_size='8pt'))

        l = np.array([2.15, 2.5])
        s1.line(l, y+y2, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data',
                      text='CH4', render_mode='css', text_font_size='8pt'))

        l = np.array([1.1692, 1.1778])
        s1.line(l[0], [y3, y3+0.02], line_color='black')
        s1.line(l[1], [y3, y3+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                      y_units='data', text='K', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([1.2437, 1.2529])
        s1.line(l[0], [y3, y3+0.02], line_color='black')
        s1.line(l[1], [y3, y3+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                      y_units='data', text='K', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([1.5168])
        s1.line(l[0], [y3, y3+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                      y_units='data', text='K', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([1.1384, 1.1409])
        s1.line(l[0], [y4, y4+0.02], line_color='black')
        s1.line(l[1], [y4, y4+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                      y_units='data', text='Na', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([1.2682])
        s1.line(l[0], [y4, y4+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                      y_units='data', text='Na', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([2.2063, 2.2090])
        s1.line(l[0], [y4, y4+0.02], line_color='black')
        s1.line(l[1], [y4, y4+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                      y_units='data', text='Na', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([2.2935, 2.3227, 2.3525, 2.3830, 2.4141])
        s1.line(l[0], [y3, y3+0.02], line_color='black')
        s1.line(l[1], [y3, y3+0.02], line_color='black')
        s1.line(l[2], [y3, y3+0.02], line_color='black')
        s1.line(l[3], [y3, y3+0.02], line_color='black')
        s1.line(l[4], [y3, y3+0.02], line_color='black')
        s1.line(l[[0, -1]], y+y3+0.02, line_color='black', line_width=1)
        s1.add_layout(Label(x=l[[0, -1]].mean(), y=y3+0.02, x_units='data',
                      y_units='data', text='CO', render_mode='css',
                      text_font_size='8pt'))

        s1.xaxis.major_label_text_font_size = '0pt'
        s1.yaxis.major_label_text_font_size = '0pt'

        # Order 2

        # Contam plot
        xlim0 = lamO2.min()
        xlim1 = lamO2.max()
        ylim0 = PA.min()-0.5*dPA
        ylim1 = PA.max()+0.5*dPA
        xlim0 = 0.614
        s5 = figure(tools=TOOLS, width=250, height=500, title=None,
                    x_range=Range1d(xlim0, xlim1), y_range=s2.y_range)
        fig_data = np.log10(np.clip(contamO2.T, 1.e-10, 1.))[:, 300:]
        s5.image([fig_data], x=xlim0, y=ylim0, dw=xlim1-xlim0, dh=ylim1-ylim0,
                 color_mapper=color_mapper)
        s5.yaxis.major_label_text_font_size = '0pt'
        s5.xaxis.axis_label = 'Wavelength (um)'

        # Line plot
        s6 = figure(tools=TOOLS, width=150, height=500, y_range=s2.y_range,
                    x_range=Range1d(100, 0), title=None)
        s6.line(100*np.sum(contamO2 >= 0.001, axis=0)/ny, PA-dPA/2,
                line_color='blue', legend='> 0.001')
        s6.line(100*np.sum(contamO2 >= 0.01, axis=0)/ny, PA-dPA/2,
                line_color='green', legend='> 0.01')
        s6.xaxis.axis_label = '% channels contam.'
        s6.yaxis.major_label_text_font_size = '0pt'

        # Dummy plots for nice spacing
        s0 = figure(tools=TOOLS, width=150, plot_height=100, title=None)
        s0.outline_line_color = "white"
        s7 = figure(tools=TOOLS, width=150, plot_height=100, title=targetName)
        s7.outline_line_color = "white"

        # Add bad PAs
        for ybad0, ybad1 in badPA:
            s5.patch([xlim0, xlim1, xlim1, xlim0],
                     [ybad1, ybad1, ybad0, ybad0],
                     color=bad_PA_color, alpha=bad_PA_alpha)
            s6.patch([0, 100, 100, 0], [ybad1, ybad1, ybad0, ybad0],
                     color=bad_PA_color, alpha=bad_PA_alpha, legend='Bad PA')

        # Line list
        s4 = figure(tools=TOOLS, width=250, plot_height=100,
                    x_range=s5.x_range, title=None)
        l = np.array([0.89, 0.99])
        s4.line(l, y+y1, line_color='black', line_width=1.5)
        s4.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                      text='H2O', render_mode='css', text_font_size='8pt'))

        l = np.array([1.09, 1.2])
        s4.line(l, y+y1, line_color='black', line_width=1.5)
        s4.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                      text='H2O', render_mode='css', text_font_size='8pt'))

        l = np.array([1.1, 1.24])
        s4.line(l, y+y2, line_color='black', line_width=1.5)
        s4.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data',
                      text='CH4', render_mode='css', text_font_size='8pt'))

        l = np.array([1.3, lamO2.max()])
        s4.line(l, y+y1, line_color='black', line_width=1.5)
        s4.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                      text='H2O', render_mode='css', text_font_size='8pt'))

        l = np.array([0.7665, 0.7699])
        s4.line(l[0], [y3, y3+0.02], line_color='black')
        s4.line(l[1], [y3, y3+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                      y_units='data', text='K', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([1.1692, 1.1778])
        s4.line(l[0], [y3, y3+0.02], line_color='black')
        s4.line(l[1], [y3, y3+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                      y_units='data', text='K', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([1.2437, 1.2529])
        s4.line(l[0], [y3, y3+0.02], line_color='black')
        s4.line(l[1], [y3, y3+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                      y_units='data', text='K', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([1.1384, 1.1409])
        s4.line(l[0], [y4, y4+0.02], line_color='black')
        s4.line(l[1], [y4, y4+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                      y_units='data', text='Na', render_mode='css',
                      text_font_size='8pt'))

        l = np.array([1.2682])
        s4.line(l[0], [y4, y4+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                      y_units='data', text='Na', render_mode='css',
                      text_font_size='8pt'))

        s4.xaxis.major_label_text_font_size = '0pt'
        s4.yaxis.major_label_text_font_size = '0pt'

        # put all the plots in a grid layout
        fig = gridplot(children=[[s7, s4, s1, s0], [s6, s5, s2, s3]])

    return fig


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

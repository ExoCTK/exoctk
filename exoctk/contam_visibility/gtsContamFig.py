import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import math as math

from astroquery.irsa import Irsa
from astroquery.simbad import Simbad
import astropy.units as u

from tqdm import tqdm

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


D2R = math.pi / 180.0
def gaussian1D(center, width, height = None, offset = None):
    """
        Written by Nate Lust
        Edited  by Jonathan Fraine

        Returns a 1D gaussian function with the given parameters

        center  = center of gaussian profile

        width   = width  of gaussian profile

        height  = height of gaussian profile
                    -- defaults to `1 / np.sqrt(2.*pi*sigma**2.)`

        offset  = background, lower limit value for gaussian
                    -- defaults to 0.0
    """

    if height == None:
        height  = np.sqrt(2.*np.pi*width**2.)
        height  = 1.0/height

    if offset == None:
        offset = 0.0

    width   = float(width)

    return lambda x: height*np.exp(-(((center - x)/width)**2)/2) + offset

def gaussian2D(center_y, center_x, width_y, width_x = None, height = None, offset = None):
    """
        Written by Nate Lust
        Edited  by Jonathan Fraine

        Returns a 2D gaussian function with the given parameters

        center_y, center_x  = center position of 2D gaussian profile

        width_y , width_x   = widths of 2D gaussian profile (if width_y != width_x, then gaussian crossection = ellipse)

        height  = height of gaussian profile
                    -- defaults to `1 / np.sqrt(2.*pi*sigma**2.)`

        offset  = background, lower limit value for gaussian
                    -- defaults to 0.0
    """

    if width_x == None:
        width_x = width_y

    if height == None:
        height = np.sqrt(2*np.pi*(width_x**2 + width_y**2))
        height = 1./height

    if offset == None:
        offset = 0.0

    width_x = float(width_x)
    width_y = float(width_y)

    return lambda y,x: height*np.exp(-(((center_x-x)/width_x)**2 + ( (center_y-y)/width_y)**2)/2)+offset

def conv1D(arr1, arr2):
    '''
        Convolve 2 arrays together
            -- used by `smooth_gaussconv`
    '''
    fft1    = np.fft.fft(arr1)
    fft2    = np.fft.fft(arr2)

    conv    = np.fft.ifft(fft1*fft2)

    return np.real(np.fft.fftshift(conv))

def conv2D(arr1, arr2):
    '''
        Convolve 2 arrays together
            -- used by `smooth_gaussconv`
    '''
    fft1    = np.fft.fft2(arr1)
    fft2    = np.fft.fft2(arr2)

    conv    = np.fft.ifft2(fft1*fft2)

    return np.real(np.fft.fftshift(conv))
def smooth_gaussconv(arr, sigma):
    '''
        Gaussian smooths `arr` over a width of `sigma`
            -- uses `conv1D` and `conv2D`, which both use np.fft
    '''
    if len(arr.shape) == 1:
        gs1D    = gaussian1D(arr.size/2, sigma)(np.arange(arr.size))
        return conv1D(gs1D, arr) / gs1D.sum()

    if len(arr.shape) == 2:
        gs2D    = gaussian2D(arr.shape[0]/2., arr.shape[1]/2., sigma)(*np.indices(arr.shape))
        return conv2D(gs2D, arr) / gs2D.sum()

    dra   = sources['ra'][(sources['k_m'] - K0) < dmag_limit] - RAd0
    ddec  = sources['dec'][(sources['k_m'] - K0) < dmag_limit] - Decd0
    dmag  = sources['k_m'][(sources['k_m'] - K0) < dmag_limit] - K0
    msize = 10**(-0.1*dmag)*15
    nlimit= sum((sources['k_m'] - K0) < dmag_limit)
    col   = np.zeros(nlimit) # No spectral collision by default

    #LWA R Field Point positions
    X_F322W2 = 479 # X Field point (undeviated); old val 514
    X_F444W = 1161 # X Field point undeviated); old val 1081
    X_buff = 10 # Keep other spectra this far away from xmin1 and xxmax1 (dist from edge)
    Y_field = 32  # Y Field point for locating object
    Y_buff = 1.0 # Keep other spectra this far away from Y_field (dist from center)

    # sweet spot ^

    xmaxWL_F322W2 = 4.013 # xmax1 wavelength of F322W2 in microns per https://jwst-docs.stsci.edu/display/JTI/NIRCam+Filters
    MINWL_F444W = 3.880 #minimum wavelength of F444W in microns per https://jwst-docs.stsci.edu/display/JTI/NIRCam+Filters

    wavelen_undev = 3.97
    disp = -0.001 #microns / pxl
    dx_F444W = 1106 # length of F444W spectrum in pixels
    dx_F322W2m1 = 1583 # length of F322W2 m=1 spectrum in pixels
    dx_F322W2m2 = math.fabs((xmaxWL_F322W2 - 2.4) * 2.0 / disp) # length of F322W2 m=2 spectrum in pixels

    xmin1_F322W2m1 = (xmaxWL_F322W2 - wavelen_undev)/disp + X_F322W2
    xmax_F322W2 = xmin1_F322W2m1 + dx_F322W2m1
    xmin1_F322W2m2 = (xmaxWL_F322W2*2.0 - wavelen_undev)/disp + X_F322W2
    xmax_F322W2m2 = xmin1_F322W2m2 + dx_F322W2m2
    # print(xmin1_F322W2m2, xmax_F322W2m2, dx_F322W2m2

    xxmax1_F444W = (MINWL_F444W -  wavelen_undev)/disp + X_F444W
    xmin1_F444W = xxmax1_F444W - dx_F444W
    # print(xmin1_F322W2m1, xmax_F322W2, xmin1_F444W, xxmax1_F444W

    # LWA field extremes for pickoff mirror POM
    xmin1_xmax1y = -142
    xmin1_miny = -132
    xxmax1_miny = 2254
    xxmax1_xmax1y = 2198

    yxmax1_minx = 2822
    ymin_minx = -174
    ymin_xmax1x = -180
    yxmax1_xmax1x = 2800

    nyquistWidth = 2. / 2.3548
    # contaminationF322W2 = np.zeros((nPA, nWaves))

    normalHeight   = 1/np.sqrt(2*np.pi*nyquistWidth*nyquistWidth)
    targetGaussian = gaussian1D(center=Y_field, width=nyquistWidth, height=normalHeight, offset=0)(y_gs+Y_field)

    xbox = [-0, -0, 4040, 4040, -0]
    ybox = [-0, 4040, 4040, -0, -0]

def plot_f322w2(contaminationF322W2_order1, contaminationF322W2_order2, badPA = []):
    """This function will use the contamination output for NIRCam's F322W2 GTS
    data to create a bokeh plot that will be outputted in ExoCTK.stsci.edu's
    Contamination Overlap page (when the user opts to plot the contamination
    plot, in addition to the visibility plot).

    Credits
    -------
    Written by ?, 20??
    Modified by Joseph Filipazo, 201?
    Modified by Jennifer V. Medina, 2018

    Parameters
    ----------

    Returns
    -------
    A bokeh plot
    """
    TOOLS = 'pan, box_zoom, crosshair, reset, hover, save'

    y = np.array([0., 0.])
    y1 = 0.07
    y2 = 0.12
    y3 = 0.17
    y4 = 0.23
    bad_PA_color = '#dddddd'
    bad_PA_alpha = 0.7

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Initializations
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    loc = 'data/contam_visibility/lambda_order1-2.txt'
    lam_file = pkg_resources.resource_filename('exoctk', loc)
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    contamO1 = contaminationF322W2_order1
    contamO2 = contaminationF322W2_order2

    ny = contamO1.shape[1]
    nPA = contamO1.shape[0]
    dPA = 362//nPA
    PA = np.arange(nPA)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Order 1 (plots on the far right)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ylim0 = 0
    ylim1 = 360
    xlim0 = 0
    xlim1 = contamO1.shape[1]

    o1data = smooth_gaussconv(contamO1, 0.5)+1e-10

    color_mapper = LinearColorMapper(palette=inferno(8)[::-1],
                                     low=0.1, high=0.2)
    color_mapper.low_color = 'white'
    color_mapper.high_color = 'black'
    s2 = figure(tools=TOOLS, width=500, height=500, title=None,
                x_range=Range1d(xlim0, xlim1),
                y_range=Range1d(ylim0, ylim1))

    s2.image([o1data], x=0, y=0, dw=xlim1-xlim0, dh=ylim1-ylim0, \
            color_mapper=color_mapper)

    disp = 0.001 # microns
    lam0 = 2.369
    lam1 = 4.417
    def tolam(x):
        dlam = disp*x
        lam = 2.369+dlam

        return lam

    # Change x-ticks from column # --to--> wavelength (microns)
    s2.xaxis.ticker = [0, 500, 1000, 1500, 2048]
    s2.xaxis.major_label_overrides = {0: str(lam0) , 500: str(tolam(500)), \
            1000: str(tolam(1000)), 1500:str(tolam(1500)) , 2048: str(lam1)}

    # Labeling
    s2.xaxis.axis_label = 'Wavelength (um)'
    s2.yaxis.axis_label = 'Position Angle (degrees)'



    # Line plot
    s3 = figure(tools=TOOLS, width=150, height=500,
                x_range=Range1d(0, 100), y_range=s2.y_range, title="order 1")
    s3.line(100*np.sum(contamO1 > 0.001, axis=1)/ny, PA-dPA/2,
            line_color='blue', legend='> 0.001')
    s3.line(100*np.sum(contamO1 > 0.01, axis=1)/ny, PA-dPA/2,
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






    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Order 2 (plots on the far left)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Contam plot (PA vs. wavelength)
    ylim0 = 0
    ylim1 = 360
    xlim0 = 0
    xlim1 = contamO2.shape[1]
    o2data = smooth_gaussconv(contamO2, 0.5)+1e-10

    s5 = figure(tools=TOOLS, width=500, height=500, title="NIRCam GTS F322W2 order 2",
                x_range=Range1d(xlim0, xlim1), y_range=s2.y_range)

    s5.image([o2data], x=0, y=0, dw=xlim1-xlim0, dh=ylim1-ylim0, \
            color_mapper=color_mapper)

    s5.xaxis.ticker = [0, 500, 1000, 1500, 2048]
    s5.xaxis.major_label_overrides = {0: str(lam0) , 500: str(tolam(500)), \
            1000: str(tolam(1000)), 1500:str(tolam(1500)) , 2048: str(lam1)}

    s5.yaxis.major_label_text_font_size = '0pt'
    s5.xaxis.axis_label = 'Wavelength (um)'




    # Line plot
    s6 = figure(tools=TOOLS, width=150, height=500, y_range=s2.y_range,
                x_range=Range1d(100, 0), title=None)
    s6.line(100*np.sum(contamO2 > 0.001, axis=1)/ny, PA-dPA/2,
            line_color='blue', legend='> 0.001')
    s6.line(100*np.sum(contamO2 > 0.01, axis=1)/ny, PA-dPA/2,
            line_color='green', legend='> 0.01')
    s6.xaxis.axis_label = '% channels contam.'
    s6.yaxis.major_label_text_font_size = '0pt'

    # Add bad PAs
    for ybad0, ybad1 in badPA:
        s5.patch([xlim0, xlim1, xlim1, xlim0],
                 [ybad1, ybad1, ybad0, ybad0],
                 color=bad_PA_color, alpha=bad_PA_alpha)
        s6.patch([0, 100, 100, 0], [ybad1, ybad1, ybad0, ybad0],
                 color=bad_PA_color, alpha=bad_PA_alpha, legend='Bad PA')

    # put all the plots in a grid layout
    fig = gridplot(children=[[s6, s5, s2, s3]])

    from bokeh.plotting import save, output_file


    output_file('/Users/jmedina/Desktop/contam_test.html')
    save(fig)
    #return fig




"""
return fig
"""










def plot_f444W(contaminationF444W_order1, badPA = []):
    """This function will use the contamination output for NIRCam's F444W GTS
    data to create a bokeh plot that will be outputted in ExoCTK.stsci.edu's
    Contamination Overlap page (when the user opts to plot the contamination
    plot, in addition to the visibility plot).

    Credits
    -------
    Written by ?, 20??
    Modified by Joseph Filipazo, 201?
    Modified by Jennifer V. Medina, 2018

    Parameters
    ----------

    Returns
    -------
    A bokeh plot
    """
    TOOLS = 'pan, box_zoom, crosshair, reset, hover, save'

    y = np.array([0., 0.])
    y1 = 0.07
    y2 = 0.12
    y3 = 0.17
    y4 = 0.23
    bad_PA_color = '#dddddd'
    bad_PA_alpha = 0.7

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Initializations
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    loc = 'data/contam_visibility/lambda_order1-2.txt'
    lam_file = pkg_resources.resource_filename('exoctk', loc)
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    contamO1 = contaminationF444W_order1
    print(np.shape(contamO1))
    print(contamO1)
    ny = contamO1.shape[1]
    nPA = contamO1.shape[0]
    dPA = 362//nPA
    PA = np.arange(nPA)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Order 1 (plots on the far right)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ylim0 = 0
    ylim1 = 360
    xlim0 = 0
    xlim1 = contamO1.shape[1]

    o1data = smooth_gaussconv(contamO1, 0.5)+1e-10

    color_mapper = LinearColorMapper(palette=inferno(8)[::-1],
                                     low=0.1, high=0.2)
    color_mapper.low_color = 'white'
    color_mapper.high_color = 'black'
    s2 = figure(tools=TOOLS, width=500, height=500, title="NIRCam GTS F444W order 1",
                x_range=Range1d(xlim0, xlim1),
                y_range=Range1d(ylim0, ylim1))

    s2.image([o1data], x=0, y=0, dw=xlim1-xlim0, dh=ylim1-ylim0, \
            color_mapper=color_mapper)

    disp = 0.001 # microns
    lam0 = 3.063
    lam1 = 5.11
    def tolam(x, instrument, filter=''):
        """x is column # of subarray"""
        dlam = disp*x
        if 'nircam' or 'NIRCam' in instrument:
            if '322' in filter:
                lam = 2.369+dlam
            elif '444' in filter:
                lam = 3.063+dlam

        if 'MIRI' or 'miri' in instrument:
            pass


        return lam

    # Change x-ticks from column # --to--> wavelength (microns)
    s2.xaxis.ticker = [0, 500, 1000, 1500, 2048]
    s2.xaxis.major_label_overrides = {0: str(lam0) , 500: str(tolam(500, 'nircam', '444')), \
            1000: str(tolam(1000, 'nircam', '444')), 1500:str(tolam(1500, 'nircam', '444')) , 2048: str(lam1)}

    # Labeling
    s2.xaxis.axis_label = 'Wavelength (um)'
    s2.yaxis.axis_label = 'Position Angle (degrees)'



    # Line plot
    s3 = figure(tools=TOOLS, width=150, height=500,
                x_range=Range1d(0, 100), y_range=s2.y_range, title=None)
    s3.line(100*np.sum(contamO1 > 0.001, axis=1)/ny, PA-dPA/2,
            line_color='blue', legend='> 0.001')
    s3.line(100*np.sum(contamO1 > 0.01, axis=1)/ny, PA-dPA/2,
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








    # put all the plots in a grid layout
    fig = gridplot(children=[[s2, s3]])

    from bokeh.plotting import save, output_file


    output_file('/Users/jmedina/Desktop/contam_testf444w.html')
    save(fig)






def plot_lrs(contamination_lrs, badPA = []):
    """This function will use the contamination output for MIRI's LRS mode
    data to create a bokeh plot that will be outputted in ExoCTK.stsci.edu's
    Contamination Overlap page (when the user opts to plot the contamination
    plot, in addition to the visibility plot).

    Credits
    -------
    Written by ?, 20??
    Modified by Joseph Filipazo, 201?
    Modified by Jennifer V. Medina, 2018

    Parameters
    ----------

    Returns
    -------
    A bokeh plot
    """
    TOOLS = 'pan, box_zoom, crosshair, reset, hover, save'

    y = np.array([0., 0.])
    y1 = 0.07
    y2 = 0.12
    y3 = 0.17
    y4 = 0.23
    bad_PA_color = '#dddddd'
    bad_PA_alpha = 0.7

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Initializations
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    loc = 'data/contam_visibility/lambda_order1-2.txt'
    lam_file = pkg_resources.resource_filename('exoctk', loc)
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    contamO1 = contamination_lrs
    print(np.shape(contamO1))
    print(contamO1)
    ny = contamO1.shape[1]
    nPA = contamO1.shape[0]
    dPA = 362//nPA
    PA = np.arange(nPA)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Order 1 (plots on the far right)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ylim0 = 0
    ylim1 = 360
    xlim0 = 0
    xlim1 = contamO1.shape[1]

    o1data = smooth_gaussconv(contamO1, 0.5)+1e-10

    color_mapper = LinearColorMapper(palette=inferno(8)[::-1],
                                     low=0.1, high=0.2)
    color_mapper.low_color = 'white'
    color_mapper.high_color = 'black'
    s2 = figure(tools=TOOLS, width=500, height=500, title="NIRCam GTS F444W order 1",
                x_range=Range1d(xlim0, xlim1),
                y_range=Range1d(ylim0, ylim1))

    s2.image([o1data], x=0, y=0, dw=xlim1-xlim0, dh=ylim1-ylim0, \
            color_mapper=color_mapper)

    disp = 0.001 # microns
    lam0 = 3.063
    lam1 = 5.11
    def tolam(x, instrument, filter=''):
        """x is column # of subarray"""
        dlam = disp*x
        if 'nircam' or 'NIRCam' in instrument:
            if '322' in filter:
                lam = 2.369+dlam
            elif '444' in filter:
                lam = 3.063+dlam

        if 'MIRI' or 'miri' in instrument:
            pass


        return lam

    # Change x-ticks from column # --to--> wavelength (microns)
    s2.xaxis.ticker = [0, 500, 1000, 1500, 2048]
    s2.xaxis.major_label_overrides = {0: str(lam0) , 500: str(tolam(500, 'nircam', '444')), \
            1000: str(tolam(1000, 'nircam', '444')), 1500:str(tolam(1500, 'nircam', '444')) , 2048: str(lam1)}

    # Labeling
    s2.xaxis.axis_label = 'Wavelength (um)'
    s2.yaxis.axis_label = 'Position Angle (degrees)'



    # Line plot
    s3 = figure(tools=TOOLS, width=150, height=500,
                x_range=Range1d(0, 100), y_range=s2.y_range, title=None)
    s3.line(100*np.sum(contamO1 > 0.001, axis=1)/ny, PA-dPA/2,
            line_color='blue', legend='> 0.001')
    s3.line(100*np.sum(contamO1 > 0.01, axis=1)/ny, PA-dPA/2,
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








    # put all the plots in a grid layout
    fig = gridplot(children=[[s2, s3]])

    from bokeh.plotting import save, output_file


    output_file('/Users/jmedina/Desktop/contam_testf444w.html')
    save(fig)



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

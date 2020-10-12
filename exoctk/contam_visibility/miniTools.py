""" The contamVerify mini tool will be a companion to ExoCTK's Contamination
Overlap tool, as it will visualize the Contaminaton Bokeh plots on the website.

Functions are:
plotTemps    - Plots the temperatures of stars according to color.
traceLengths - Fine-tunes the trace lengths in the plot.
contamVerify - The main mini tool. Outputs a .pdf file with one or more figures
               showing the FOV in the science frame according to the input
               Aperture Position Angle(s) it is fed.

Author(s)
---------
Jennifer V. Medina, 2020
"""
import astropy.coordinates as crd
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os
import pysiaf

from astropy.io import fits
from astroquery.irsa import Irsa
from matplotlib import cm
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
from scipy.io import readsav

EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
if not EXOCTK_DATA:
    print(
        'WARNING: The $EXOCTK_DATA environment variable is not set. Contamination overlap will not work. Please set the '
        'value of this variable to point to the location of the exoctk_data '
        'download folder.  Users may retreive this folder by clicking the '
        '"ExoCTK Data Download" button on the ExoCTK website, or by using '
        'the exoctk.utils.download_exoctk_data() function.')
    TRACES_PATH = None

TRACES_PATH = os.path.join(EXOCTK_DATA, 'exoctk_contam', 'traces')


def plotTemps(TEMPS, allRA, allDEC):
    """ The stars' colors in the plot will be a function of effective stellar
    temperatures when plotting with this function. """

    # Getting the color palette
    colors = cm.get_cmap('viridis', len(TEMPS))
    colors_0 = np.asarray(colors.colors)

    # Assigning index arrays to TEMPS array
    i = TEMPS.argsort()
    ii = TEMPS.argsort().argsort()

    # Matching the colors to the corresponding magnitude
    colors = colors_0
    starsx, starsy = allRA[i], allDEC[i]

    plt.style.use('dark_background')
    for x, y, c in zip(starsx, starsy, colors):
        plt.scatter(
            x,
            y,
            marker='*',
            s=150,
            color=c,
            picker=True,
            lw=0.5,
            edgecolor='white')

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis,
                               norm=plt.Normalize(vmin=TEMPS.min(),
                                                  vmax=TEMPS.max()))
    sm._A = []

    cbar = plt.colorbar(sm, fraction=0.046, pad=0.04)
    cbar.set_label('effective Temperature (K)', fontsize=20)


def traceLength(inst):
    """ For fine-tuning the trace lengths in the contamVerify output figures """

    # Getting example trace to calculate rough estimate of trace lengths
    if 'NIRCam' in inst:
        FILE = 'rot_o1_6000.0.fits'
    elif 'MIRI' in inst:
        FILE = 'LOWbg_6000.0.fits'
    elif 'NIRISS' in inst:
        FILE = 'modelOrder12_teff6000.sav'

    trFile = os.path.join(TRACES_PATH, inst.replace(' ', '_'), FILE)
    trData = readsav(trFile)['modelo12'] if 'NIRISS' in inst \
        else fits.getdata(trFile, 1)
    trData = trData[0]
    print(np.shape(trData))
    ax = 1 if 'NIRCam' in inst else 0
    peak = trData.max()

    # the length of the trace
    targ_trace_start = np.where(trData > 0.0001 * peak)[ax].min()
    targ_trace_stop = np.where(trData > 0.0001 * peak)[ax].max()

    return targ_trace_start, targ_trace_stop


def contamVerify(RA, DEC, INSTRUMENT, APAlist, binComp=[], PDF='', web=False):
    """ Generates a PDF file of figures displaying a simulation
    of the science image for any given observation using the parameters provided.

    Parameter(s)
    ------------
    RA  : str
        The Right Ascension of your target in HH:MM:SS
    DEC : str
        The Declination of your target in DD:MM:SS
    INSTRUMENT : str
        The instrument you are observing with (case-sensitive).
        The software currently supports:
        'MIRI', 'NIRISS', 'NIRCam F322W2', 'NIRCam F444W'
    APAlist : list
        A list of Aperture Position Angle(s). Element(s) must be in integers.
        Example 1:
        [1, 25, 181, 205]
        Example 2:
        [25]
    binComp : list
        A list containing parameters of a missing companion that is not
        in the 2MASS IRSA point-source catalog. The format is:
        [RA (arcseconds), DEC (arcseconds), J mag, H mag, K mag]
        [string, string, integer, integer, integer]
    PDF : string
        The path to where the PDF file will be saved. If left blank, the PDF
        file will be saved in your current working directory.
        Example:
        'path/to/my/file.pdf'
    web : boolean
        Makes it easier to integrate it onto the website. Leave this as false,
        unless you're running this in app_exoctk.py

    Returns
    -------
    A .PDF file containing a simulation of the FOV of your target in the
    science coordinate system. Some things to consider when reading the figures:

    1. The target is circled in red
    2. Stellar temperatures of all sources are plotted by color
    3. The gray region oulined in blue represents the aperture for the given
       instrument.
    4. The blue square represents the readout region, or the "origin"

    """
    print('Generating FOV...')
    # Decimal degrees --> HMSDMS for Irsa.query_region()
    targetcrd = crd.SkyCoord(ra=RA, dec=DEC, unit='deg').to_string('hmsdms')
    targetRA, targetDEC = RA, DEC

    # Querying for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    rad = 2.5
    print('Querying for point-sources within {} arcminutes...'.format(str(rad)))
    info = Irsa.query_region(
        targetcrd,
        catalog='fp_psc',
        spatial='Cone',
        radius=rad * u.arcmin)

    # Coordinates of all stars in FOV, including target
    allRA = info['ra'].data.data
    allDEC = info['dec'].data.data

    # Initiating a dictionary to hold all relevant star information
    stars = {}
    stars['RA'], stars['DEC'] = allRA, allDEC

    print('Total point-sources found in region: {}'.format(len(stars['RA'])))
    # Finding the target using relative distances
    sindRA = (targetRA - stars['RA']) * np.cos(targetDEC)
    cosdRA = targetDEC - stars['DEC']
    distance = np.sqrt(sindRA**2 + cosdRA**2)
    targetIndex = np.argmin(distance)

    # Appending missing companion to the above lists (if any)
    if binComp != []:
        print('Adding missing companion...')
        bb = binComp[0] / 3600 / np.cos(allDEC[targetIndex] * deg2rad)
        allRA = np.append(allRA, (allRA[targetIndex] + bb))
        allDEC = np.append(allDEC, (allDEC[targetIndex] + binComp[1] / 3600))
        Jmag = np.append(Jmag, binComp[2])
        Hmag = np.append(Kmag, binComp[3])
        Kmag = np.append(Kmag, binComp[4])
        J_Hobs = Jmag - Hmag
        H_Kobs = Hmag - Kmag

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

    # JHK bands of all stars in FOV, including target
    Jmag = info['j_m'].data.data
    Hmag = info['h_m'].data.data
    Kmag = info['k_m'].data.data
    # J-H band, H-K band. This will be used to derive the stellar Temps later
    J_Hobs = Jmag - Hmag
    H_Kobs = Hmag - Kmag

    # Number of stars
    nStars = stars['RA'].size

    # Find/assign Teff of each star
    print('Calculating effective temperatures...')
    starsT = np.empty(nStars)
    for j in range(nStars):
        color_separation = (J_Hobs[j] - jhMod)**2 + (H_Kobs[j] - hkMod)**2
        min_separation_ind = np.argmin(color_separation)
        starsT[j] = teffMod[min_separation_ind]

    # Record keeping
    stars['Temp'] = starsT

    # Initiating a dictionary for customizability
    apertures = {}
    apertures['NIRISS'] = ['NIS_SOSSFULL', 'NIS_SOSSFULL']
    apertures['NIRCam F444W'] = ['NRCA5_GRISM256_F444W', 'NRCA5_FULL']
    apertures['NIRCam F322W2'] = ['NRCA5_GRISM256_F322W2', 'NRCA5_FULL']
    apertures['MIRI'] = ['MIRIM_SLITLESSPRISM', 'MIRIM_FULL']

    # Instantiate SIAF object
    siaf = pysiaf.Siaf(INSTRUMENT.split(' ')[0])

    aper = siaf.apertures[apertures[INSTRUMENT][0]]
    full = siaf.apertures[apertures[INSTRUMENT][1]]

    # DET_targ -> TEL_targ -> get attitude matrix for target
    # -> TEL_neighbor -> DET_neighbor -> SCI_neighbor
    print('Converting Sky --> Science coordinates...')
    xSweet, ySweet = aper.reference_point('det')

    v2targ, v3targ = aper.det_to_tel(xSweet, ySweet)

    contam = {}

    if not web:
        filename = 'contam_{}_{}_{}.pdf'.format(RA, DEC, INSTRUMENT)
        defaultPDF = os.path.join(os.getcwd(), filename).replace(' ', '_')
        PDF = defaultPDF if PDF == '' else PDF
    elif web:
        filename = 'contam_{}_{}_{}.pdf'.format(RA, DEC, INSTRUMENT)
        PDF = os.path.join(TRACES_PATH, filename)

    print('Saving figures to: {}'.format(PDF))
    print('This will take a second...')
    pdfobj = PdfPages(PDF)
    for APA in APAlist:

        attitude = pysiaf.utils.rotations.attitude_matrix(
            v2targ, v3targ, targetRA, targetDEC, APA)

        xdet, ydet = [], []
        xsci, ysci = [], []

        for starRA, starDEC in zip(stars['RA'], stars['DEC']):
            # Get the TEL coordinates of each star using the attitude
            # matrix of the target
            V2, V3 = pysiaf.utils.rotations.sky_to_tel(
                attitude, starRA, starDEC)
            # Convert to arcsec and turn to a float
            V2, V3 = V2.to(u.arcsec).value, V3.to(u.arcsec).value

            XDET, YDET = aper.tel_to_det(V2, V3)
            XSCI, YSCI = aper.det_to_sci(XDET, YDET)

            xdet.append(XDET)
            ydet.append(YDET)
            xsci.append(XSCI)
            ysci.append(YSCI)

        XDET, YDET = np.array(xdet), np.array(ydet)
        XSCI, YSCI = np.array(xsci), np.array(ysci)

        starsAPA = {'xdet': XDET, 'ydet': YDET, 'xsci': XSCI, 'ysci': YSCI}

        # Finding indexes of neighbor sources that land on detector
        rows, cols = full.corners('det')

        minrow, maxrow = rows.min(), rows.max()
        mincol, maxcol = cols.min(), cols.max()

        inFOV = []
        for star in range(0, nStars):

            x, y = starsAPA['xdet'][star], starsAPA['ydet'][star]
            if (mincol < x) & (x < maxcol) & (minrow < y) & (y < maxrow):
                inFOV.append(star)

        inFOV = np.array(inFOV)

        # Making final plot
        fig = plt.figure(figsize=(15, 15))
        aper.plot(frame='sci', fill_color='gray', color='blue')
        plt.scatter(
            XSCI[targetIndex],
            YSCI[targetIndex],
            s=400,
            lw=1.5,
            facecolor='gray',
            edgecolor='red')
        plotTemps(starsT[inFOV], XSCI[inFOV], YSCI[inFOV])
        aper.plot_frame_origin(frame='sci', which='sci')

        # Fine-tune trace lengths
        start, stop = traceLength(INSTRUMENT)

        # Plotting the trace footprints
        for x, y in zip(XSCI[inFOV], YSCI[inFOV]):

            if 'F322W2' in INSTRUMENT:
                plt.plot([x - stop, x + start], [y, y],
                         lw=40, color='white', alpha=0.2)
                plt.plot([x - stop, x + start],
                         [y, y], lw=2., color='white')
            elif 'F444W' in INSTRUMENT:
                plt.plot([x - start, x + stop], [y, y],
                         lw=40, color='white', alpha=0.2)
                plt.plot([x - start, x + stop],
                         [y, y], lw=2., color='white')
            else:
                plt.plot([x, x], [y - stop, y + start],
                         lw=40, color='white', alpha=0.2)
                plt.plot([x, x], [y - stop, y + start],
                         lw=2., color='white')

        # Labeling
        aperstr = str(aper.AperName.replace('_', ' '))
        tx, ty = str(
            round(
                XSCI[targetIndex])), str(
            round(
                YSCI[targetIndex]))
        plt.title(
            'The FOV in SCIENCE coordinates at APA {}$^o$'.format(
                str(APA)) +
            '\n' +
            '{}'.format(aperstr) +
            '\n' +
            'Target (X,Y): {}, {}'.format(
                tx,
                ty),
            fontsize=20)

        # Adding to PDF
        pdfobj.savefig(fig, bbox_inches='tight')

    pdfobj.close()

    if web:
        return PDF

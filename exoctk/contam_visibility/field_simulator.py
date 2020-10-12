import glob
import os
import pysiaf

import astropy.coordinates as crd
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astroquery.irsa import Irsa
from matplotlib import cm
from scipy.io import readsav
from astropy.io import fits
from exoctk.utils import get_env_variables
from pysiaf.utils import rotations

EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
if not EXOCTK_DATA:
    print('WARNING: The $EXOCTK_DATA environment variable is not set. '
          'Contamination overlap will not work. Please set the '
          'value of this variable to point to the location of the exoctk_data '
          'download folder.  Users may retreive this folder by clicking the '
          '"ExoCTK Data Download" button on the ExoCTK website, or by using '
          'the exoctk.utils.download_exoctk_data() function.'
          )
    TRACES_PATH = None
else:
    TRACES_PATH = os.path.join(EXOCTK_DATA, 'exoctk_contam', 'traces')


def sossFieldSim(ra, dec, binComp='', dimX=256):
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
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg))
    targetRA = targetcrd.ra.value
    targetDEC = targetcrd.dec.value
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
    distance = np.sqrt(aa**2 + (targetDEC - allDEC)**2)
    targetIndex = np.argmin(distance)

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
        color_separation = (J_Hobs[j] - jhMod)**2 + (H_Kobs[j] - hkMod)**2
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
             dy=np.empty(nStars))
    stars = np.empty(nStars,
                     dtype=[(key, val.dtype) for key, val in _.items()])
    for key, val in _.items():
        stars[key] = val

    # Initialize final fits cube that contains the modelled traces
    # with contamination
    PAmin = 0  # instrument PA, degrees
    PAmax = 360
    dPA = 1  # degrees

    # Set of IPA values to cover
    PAtab = np.arange(PAmin, PAmax, dPA)    # degrees
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

    for kPA in range(PAtab.size):
        APA = PAtab[kPA]
        print('Generating field at APA : {}'.format(str(APA)))

        V3PA = APA + 0.57  # from APT

        sindx = np.sin((np.pi / 2) + APA / radeg) * stars['dDEC']
        cosdx = np.cos((np.pi / 2) + APA / radeg) * stars['dDEC']
        nps = niriss_pixel_scale
        stars['dx'] = (np.cos((np.pi / 2) + APA / radeg)
                       * stars['dRA'] - sindx) / nps
        stars['dy'] = (np.sin((np.pi / 2) + APA / radeg)
                       * stars['dRA'] + cosdx) / nps
        stars['x'] = stars['dx'] + sweetSpot['x']
        stars['y'] = stars['dy'] + sweetSpot['y']

        # Retain stars that are within the Direct Image NIRISS POM FOV
        ind, = np.where((stars['x'] >= -162) & (stars['x'] <= 2047 + 185) &
                        (stars['y'] >= -154) & (stars['y'] <= 2047 + 174))
        starsInFOV = stars[ind]

        for i in range(len(ind)):
            intx = round(starsInFOV['dx'][i])
            inty = round(starsInFOV['dy'][i])

            k = np.where(teffMod == starsInFOV['T'][i])[0][0]

            fluxscale = 10.0**(-0.4 *
                               (starsInFOV['jmag'][i] - sweetSpot['jmag']))

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
                simuCube[kPA + 2, y0:y0 + my1 - my0,
                         x0:x0 + mx1 - mx0] += mod * fluxscale
    return simuCube


def gtsFieldSim(ra, dec, filter, binComp=''):
    """ Produce a Grism Time Series field simulation for a target.
    Parameters
    ----------
    ra : float
        The RA of the target.
    dec : float
        The Dec of the target.
    filter : str
        The NIRCam filter being used. Can either be:
        'F444W' or 'F322W2' (case-sensitive)
    binComp : sequence
        The parameters of a binary companion.

    Returns
    -------
    simuCube : np.ndarray
        The simulated data cube. Index 0 and 1 (axis=0) show the trace of
        the target for orders 1 and 2 (respectively). Index 2-362 show the trace
        of the target at every position angle (PA) of the instrument.
    """
    # Instantiate a pySIAF object
    siaf = pysiaf.Siaf('NIRCam')

    full = siaf.apertures['NRCA5_FULL']
    if filter == 'F444W':
        aper = siaf.apertures['NRCA5_GRISM256_F444W']
    elif filter == 'F322W2':
        aper = siaf.apertures['NRCA5_GRISM256_F322W2']

    # Calling the variables
    deg2rad = np.pi / 180
    subX, subY = aper.XSciSize, aper.YSciSize
    rad = 2.5  # arcmins
    pixel_scale = 0.063  # arsec/pixel
    V3PAs = np.arange(0, 360, 1)
    nPA = len(V3PAs)
    # Generate cube of field simulation at every degree of APA rotation
    simuCube = np.zeros([nPA + 1, subY, subX])
    xSweet, ySweet = aper.reference_point('det')
    add_to_v3pa = aper.V3IdlYAngle
    # NIRCam Full Frame dimensions
    rows, cols = full.corners('det')
    minrow, maxrow = rows.min(), rows.max()
    mincol, maxcol = cols.min(), cols.max()

    #############################STEP 1#####################################
    ########################################################################
    # Converting to degrees
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg))
    targetRA = targetcrd.ra.value
    targetDEC = targetcrd.dec.value

    # Querying for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    info = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone',
                             radius=rad * u.arcmin)

    # Coordinates of all the stars in FOV, including target
    allRA = info['ra'].data.data
    allDEC = info['dec'].data.data

    # Initiating a dictionary to hold all relevant star information
    stars = {}
    stars['RA'], stars['DEC'] = allRA, allDEC

    #############################STEP 2#####################################
    ########################################################################
    sindRA = (targetRA - stars['RA']) * np.cos(targetDEC)
    cosdRA = targetDEC - stars['DEC']
    distance = np.sqrt(sindRA**2 + cosdRA**2)
    if np.min(distance) > 1.0*(10**-4):
        coords = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg)).to_string('decimal')
        ra, dec = coords.split(' ')[0], coords.split(' ')[1]
        raise Exception('Unable to detect a source with coordinates [RA: {}, DEC: {}] within IRSA`s 2MASS Point-Source Catalog. Please enter different coordinates or contact the JWST help desk.'.format(str(ra), str(dec)))

    targetIndex = np.argmin(distance)

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

    #############################STEP 3#####################################
    ########################################################################
    # JHK bands of all stars in FOV, including target
    Jmag = info['j_m'].data.data
    Hmag = info['h_m'].data.data
    Kmag = info['k_m'].data.data
    # J-H band, H-K band. This will be used to derive the Teff
    J_Hobs = Jmag - Hmag
    H_Kobs = Hmag - Kmag

    # Add any missing companion
    if binComp != '':
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

    #############################STEP 4#####################################
    ########################################################################
    # Calculate corresponding V2/V3 (TEL) coordinates for Sweetspot
    v2targ, v3targ = aper.det_to_tel(xSweet, ySweet)

    for V3PA in range(0, nPA, 1):
        # Get APA from V3PA
        APA = V3PA + add_to_v3pa
        if APA > 360:
            APA = APA-360
        elif APA < 0:
            APA = APA+360

        print('Generating field at APA : {}'.format(str(APA)))

        # Get target's attitude matrix for each Position Angle
        attitude = rotations.attitude_matrix(v2targ, v3targ,
                                             targetRA, targetDEC,
                                             APA)

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

        sci_targx, sci_targy = stars['xsci'][targetIndex],\
            stars['ysci'][targetIndex]

    #############################STEP 5#####################################
    ########################################################################
        inFOV = []
        for star in range(0, nStars):

            x, y = stars['xdet'][star], stars['ydet'][star]
            if (mincol < x) & (x < maxcol) & (minrow < y) & (y < maxrow):
                inFOV.append(star)

        inFOV = np.array(inFOV)

    #############################STEP 6#####################################
    ########################################################################
        nircam_path = 'NIRCam_F444W' if filter == 'F444W' else 'NIRCam_F322W2'
        fitsFiles = glob.glob(
            os.path.join(
                TRACES_PATH,
                nircam_path,
                'rot*.fits'))
        fitsFiles = np.sort(fitsFiles)

        for idx in inFOV:

            sci_dx = round(sci_targx - stars['xsci'][idx])
            sci_dy = round(sci_targy - stars['ysci'][idx])
            temp = stars['Temp'][idx]

            for file in fitsFiles:
                if str(temp) in file:
                    trace = fits.getdata(file, 1)[0]

            fluxscale = 10.0**(-0.4 * \
                               (stars['Jmag'][idx] - stars['Jmag'][targetIndex]))

            # Padding array
            pad_trace = np.pad(trace, pad_width=5000, mode='constant',
                               constant_values=0)

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

            traceX, traceY = np.shape(pad_trace)[1], np.shape(pad_trace)[0]
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

            # Fleshing out index 0 of the simulation cube (trace of target)
            if (sci_dx == 0) & (sci_dy == 0):  # this is the target

                tr = pad_trace[my0:my1, mx0:mx1] * fluxscale
                trX, trY = np.shape(tr)[1], np.shape(tr)[0]

                simuCube[0, 0:trY, 0:trX] = tr

            # Fleshing out indexes 1-361 of the simulation cube
            # (trace of neighboring stars at every position angle)
            else:

                tr = pad_trace[my0:my1, mx0:mx1] * fluxscale
                trX, trY = np.shape(tr)[1], np.shape(tr)[0]
                simuCube[V3PA + 1, 0:trY, 0:trX] += tr

    return simuCube


def lrsFieldSim(ra, dec, binComp=''):
    """ Produce a Grism Time Series field simulation for a target.
    Parameters
    ----------
    ra : float
        The RA of the target.
    dec : float
        The Dec of the target.
    binComp : sequence
        The parameters of a binary companion.

    Returns
    -------
    simuCube : np.ndarray
        The simulated data cube. Index 0 and 1 (axis=0) show the trace of
        the target for orders 1 and 2 (respectively). Index 2-362 show the trace
        of the target at every position angle (PA) of the instrument.
    """
    #############################INSTRUMENT PARAMETERS######################
    ########################################################################
    # Instantiate a pySIAF object
    siaf = pysiaf.Siaf('MIRI')
    aper = siaf.apertures['MIRIM_SLITLESSPRISM']
    full = siaf.apertures['MIRIM_FULL']

    # Calling the variables
    deg2rad = np.pi / 180
    subX, subY = aper.XSciSize, aper.YSciSize
    rad = 2.0  # arcmins
    pixel_scale = 0.11  # arsec/pixel
    V3PAs = np.arange(0, 360, 1)
    nPA = len(V3PAs)
    # Generate cube of field simulation at every degree of APA rotation
    simuCube = np.zeros([nPA + 1, subY, subX])
    xSweet, ySweet = aper.reference_point('det')
    add_to_v3pa = aper.V3IdlYAngle
    # MIRI Full Frame dimensions
    rows, cols = full.corners('det')
    minrow, maxrow = rows.min(), rows.max()
    mincol, maxcol = cols.min(), cols.max()

    #############################STEP 1#####################################
    ########################################################################
    # Converting to degrees
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg))
    targetRA = targetcrd.ra.value
    targetDEC = targetcrd.dec.value

    # Querying for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    info = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone',
                             radius=rad * u.arcmin)

    # Coordinates of all the stars in FOV, including target
    allRA = info['ra'].data.data
    allDEC = info['dec'].data.data

    # Initiating a dictionary to hold all relevant star information
    stars = {}
    stars['RA'], stars['DEC'] = allRA, allDEC

    #############################STEP 2#####################################
    ########################################################################
    sindRA = (targetRA - stars['RA']) * np.cos(targetDEC)
    cosdRA = targetDEC - stars['DEC']
    distance = np.sqrt(sindRA**2 + cosdRA**2)
    if np.min(distance) > 1.0*(10**-4):
        coords = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg)).to_string('decimal')
        ra, dec = coords.split(' ')[0], coords.split(' ')[1]
        raise Exception('Unable to detect a source with coordinates [RA: {}, DEC: {}] within IRSA`s 2MASS Point-Source Catalog. Please enter different coordinates or contact the JWST help desk.'.format(str(ra), str(dec)))

    targetIndex = np.argmin(distance)

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

    #############################STEP 3#####################################
    ########################################################################
    # JHK bands of all stars in FOV, including target
    Jmag = info['j_m'].data.data
    Hmag = info['h_m'].data.data
    Kmag = info['k_m'].data.data
    # J-H band, H-K band. This will be used to derive the Teff
    J_Hobs = Jmag - Hmag
    H_Kobs = Hmag - Kmag

    # Add any missing companion
    if binComp != '':
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

    #############################STEP 4#####################################
    ########################################################################
    # Calculate corresponding V2/V3 (TEL) coordinates for Sweetspot
    v2targ, v3targ = aper.det_to_tel(xSweet, ySweet)

    for V3PA in range(0, nPA, 1):
        # Get APA from V3PA
        APA = V3PA + add_to_v3pa
        if APA > 360:
            APA = APA-360
        elif APA < 0:
            APA = APA+360

        print('Generating field at APA : {}'.format(str(APA)))

        # Get target's attitude matrix for each Position Angle
        attitude = rotations.attitude_matrix(v2targ, v3targ,
                                             targetRA, targetDEC,
                                             APA)

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

        sci_targx, sci_targy = stars['xsci'][targetIndex],\
            stars['ysci'][targetIndex]
    #############################STEP 5#####################################
    ########################################################################
        inFOV = []
        for star in range(0, nStars):

            x, y = stars['xdet'][star], stars['ydet'][star]
            if (mincol < x) & (x < maxcol) & (minrow < y) & (y < maxrow):
                inFOV.append(star)

        inFOV = np.array(inFOV)

    #############################STEP 6#####################################
    ########################################################################
        fitsFiles = glob.glob(os.path.join(TRACES_PATH, 'MIRI', 'LOW*.fits'))
        fitsFiles = np.sort(fitsFiles)

        for idx in inFOV:

            sci_dx = round(sci_targx - stars['xsci'][idx])
            sci_dy = round(sci_targy - stars['ysci'][idx])

            temp = stars['Temp'][idx]

            for file in fitsFiles:
                if str(temp) in file:
                    trace = fits.getdata(file)[0]

            fluxscale = 10.0**(-0.4 * \
                               (stars['Jmag'][idx] - stars['Jmag'][targetIndex]))

            # Padding array
            pad_trace = np.pad(trace, pad_width=5000, mode='constant',
                               constant_values=0)

            # Determine the highest pixel value of trace
            maxY, maxX = np.where(pad_trace == pad_trace.max())
            peakY, peakX = maxY[0], maxX[0]

            # Use relative distances (sci_dx, sci_dy) to find target
            # xTarg,yTarg are essentially the "sweetspot" in the PADDED arr
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

            traceX, traceY = np.shape(pad_trace)[1], np.shape(pad_trace)[0]
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

            # Fleshing out index 0 of the simulation cube (trace of target)
            if (sci_dx == 0) & (sci_dy == 0):  # this is the target

                tr = pad_trace[my0:my1, mx0:mx1] * fluxscale
                trX, trY = np.shape(tr)[1], np.shape(tr)[0]

                simuCube[0, 0:trY, 0:trX] = tr

            # Fleshing out indexes 1-361 of the simulation cube
            # (trace of neighboring stars at every position angle)
            else:

                tr = pad_trace[my0:my1, mx0:mx1] * fluxscale
                trX, trY = np.shape(tr)[1], np.shape(tr)[0]
                simuCube[V3PA + 1, 0:trY, 0:trX] += tr

    return simuCube


def fieldSim(ra, dec, instrument, binComp='', testing=False):
    """ Wraps ``sossFieldSim``, ``gtsFieldSim``, and ``lrsFieldSim`` together.
    Produces a field simulation for a target using any instrument (NIRISS,
    NIRCam, or MIRI).

    Parameters
    ----------
    ra : float
        The RA of the target.
    dec : float
        The Dec of the target.
    instrument : str
        The instrument the contamination is being calculated for.
        Can either be (case-sensitive):
        'NIRISS', 'NIRCam F322W2', 'NIRCam F444W', 'MIRI'
    binComp : sequence
        The parameters of a binary companion.
    testing : bool
        Shoud be ``True`` if running fieldSim for testing / troubleshooting
        purposes. This will generate a matplotlib figure showing the target
        FOV. The neighboring stars in this FOV will be included in the
        contamination calculation (contamFig.py).

    Returns
    -------
    simuCube : np.ndarray
        The simulated data cube. Index 0 and 1 (axis=0) show the trace of
        the target for orders 1 and 2 (respectively). Index 2-362 show the trace
        of the target at every position angle (PA) of the instrument.
    plt.plot() : matplotlib object
        A plot. Only if `testing` parameter is set to True.
    """

    # Calling the variables which depend on what instrument you use
    if instrument == 'NIRISS':
        simuCube = sossFieldSim(ra, dec, binComp)

    elif instrument == 'NIRCam F444W':
        simuCube = gtsFieldSim(ra, dec, 'F444W', binComp)

    elif instrument == 'NIRCam F322W2':
        simuCube = gtsFieldSim(ra, dec, 'F322W2', binComp)

    elif instrument == 'MIRI':
        simuCube = lrsFieldSim(ra, dec, binComp)

    return simuCube


if __name__ == '__main__':
    ra, dec = "04 25 29.0162", "-30 36 01.603"  # Wasp 79
    #sossFieldSim(ra, dec)
    if EXOCTK_DATA:
        fieldSim(ra, dec, instrument='NIRISS')

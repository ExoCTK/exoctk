import glob
import os
import pysiaf
from multiprocessing import pool, cpu_count
import time
from functools import partial

import astropy.coordinates as crd
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astroquery.irsa import Irsa
from matplotlib import cm
from scipy.io import readsav
from astropy.io import fits
from exoctk.utils import get_env_variables, check_for_data
from pysiaf.utils import rotations


APERTURES = {'NIS_SUBSTRIP96': ['NIRISS', 'NIS_SOSSFULL', 'NIRISS/modelOrder12_teff*'],
             'NIS_SUBSTRIP256': ['NIRISS', 'NIS_SOSSFULL', 'NIRISS/modelOrder12_teff*'],
             'NRCA5_GRISM256_F444W': ['NIRCam', 'NRCA5_FULL', 'NIRCam_F444W/*'],
             'NRCA5_GRISM256_F322W2': ['NIRCam', 'NRCA5_FULL', 'NIRCam_F322W2/*'],
             'MIRI_SLITLESSPRISM': ['MIRI', 'MIRIM_FULL', 'MIRI/_*']}

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
    instrument, full, tracedir = APERTURES[aperture]

    # Get the path to the trace files
    traces_path = os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/traces/', tracedir)

    # Glob the file names
    trace_files = glob.glob(traces_path)

    # Make into dict
    trace_dict = {os.path.basename(file).split('_')[-1].replace('teff', '').split('.')[0]: fits.getdata(file)[0] if instrument == 'MIRI' else fits.getdata(file, 1)[0] if instrument == 'NIRCam' else readsav(file, verbose=False)['modelo12'] for file in trace_files}

    return trace_dict


def field_simulation(ra, dec, aperture, binComp='', n_jobs=-1):
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
    instrument, full_aper, trace_dir = APERTURES[aperture]
    siaf = pysiaf.Siaf(instrument)

    # Get the full and subarray apertures
    full = siaf.apertures[full_aper]
    aper = siaf.apertures[aperture]

    # Calling the variables
    deg2rad = np.pi / 180
    subX, subY = aper.XSciSize, aper.YSciSize
    rad = 2.5  # arcmins
    pixel_scale = 0.063  # arsec/pixel
    V3PAs = np.arange(0, 360, 1)
    nPA = len(V3PAs)

    # Generate cube of field simulation at every degree of APA rotation
    xSweet, ySweet = aper.reference_point('det')
    add_to_v3pa = aper.V3IdlYAngle

    # Full Frame dimensions
    rows, cols = full.corners('det')
    minrow, maxrow = rows.min(), rows.max()
    mincol, maxcol = cols.min(), cols.max()

    # Determine the traces for the given instrument
    traces = get_traces(aperture)

    #############################STEP 1#####################################
    ########################################################################

    # Converting to degrees
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=u.deg)
    targetRA = targetcrd.ra.value
    targetDEC = targetcrd.dec.value

    # Querying for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    info = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone', radius=rad * u.arcmin)

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
    # if np.min(distance) > 1.0*(10**-4):
    #     coords = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg)).to_string('decimal')
    #     ra, dec = coords.split(' ')[0], coords.split(' ')[1]
    #     raise Exception('Unable to detect a source with coordinates [RA: {}, DEC: {}] within IRSA`s 2MASS Point-Source Catalog. Please enter different coordinates or contact the JWST help desk.'.format(str(ra), str(dec)))

    targetIndex = np.argmin(distance)

    # Restoring model parameters
    modpath = os.path.join(os.environ['EXOCTK_DATA'], 'exoctk_contam/traces/NIRISS/modelsInfo.sav')
    modelParam = readsav(modpath, verbose=False)
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
    V3PAs = list(range(0, nPA, 1))

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

        # Get the trace
        trace = traces[str(int(temp))]
        if trace.ndim == 3:
            trace = np.sum()

        # Get the flux scaling
        fluxscale = 10.0**(-0.4 * (stars['Jmag'][idx] - stars['Jmag'][targetIndex]))

        # Padding array
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

        # Add trace to the image
        tr = pad_trace[my0:my1, mx0:mx1] * fluxscale
        trX, trY = np.shape(tr)[1], np.shape(tr)[0]
        simuImage[0:trY, 0:trX] += tr

        # # Fleshing out index 0 of the simulation cube (trace of target)
        # if (sci_dx == 0) & (sci_dy == 0):  # this is the target
        #
        #     tr = pad_trace[my0:my1, mx0:mx1] * fluxscale
        #     trX, trY = np.shape(tr)[1], np.shape(tr)[0]
        #
        #     simuImage[0:trY, 0:trX] = tr
        #
        # # Fleshing out indexes 1-361 of the simulation cube
        # # (trace of neighboring stars at every position angle)
        # else:
        #
        #     tr = pad_trace[my0:my1, mx0:mx1] * fluxscale
        #     trX, trY = np.shape(tr)[1], np.shape(tr)[0]
        #     simuImage[0:trY, 0:trX] += tr

    return simuImage



if __name__ == '__main__':
    ra, dec = "04 25 29.0162", "-30 36 01.603"  # Wasp 79
    field_simulation(ra, dec, 'NIS_SUBSTRIP256')

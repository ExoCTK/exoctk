import glob
import os

from astroquery.irsa import Irsa
import astropy.coordinates as crd
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav
from astropy.io import fits
from exoctk.utils import get_env_variables

TRACES_PATH = os.path.join(os.environ.get('EXOCTK_DATA'),  'exoctk_contam', 'traces')
if TRACES_PATH == '':
    raise NameError("You need to have an exported 'EXOCTK_DATA' environment \
                     variable and data set up before we can continue.")

def sossFieldSim(ra, dec, binComp='', dimX=256):
    """ Produce a SOSS field simulation for a target

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
    # stars in large field around target
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg))
    targetRA = targetcrd.ra.value
    targetDEC = targetcrd.dec.value
    info = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone',
                             radius=2.5*u.arcmin)

    # coordinates of all stars in FOV, including target
    allRA = info['ra'].data.data
    allDEC = info['dec'].data.data
    Jmag = info['j_m'].data.data
    Hmag = info['h_m'].data.data
    Kmag = info['k_m'].data.data
    J_Hobs = Jmag-Hmag
    H_Kobs = Hmag-Kmag

    # target coords
    aa = ((targetRA-allRA)*np.cos(targetDEC))
    distance = np.sqrt(aa**2 + (targetDEC-allDEC)**2)
    targetIndex = np.argmin(distance)  # the target

    # add any missing companion
    if binComp != '':
        deg2rad = np.pi/180
        bb = binComp[0]/3600/np.cos(allDEC[targetIndex]*deg2rad)
        allRA = np.append(allRA, (allRA[targetIndex] + bb))
        allDEC = np.append(allDEC, (allDEC[targetIndex] + binComp[1]/3600))
        Jmag = np.append(Jmag, binComp[2])
        Hmag = np.append(Kmag, binComp[3])
        Kmag = np.append(Kmag, binComp[4])
        J_Hobs = Jmag-Hmag
        H_Kobs = Hmag-Kmag

    # number of stars
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

    # find/assign Teff of each star
    starsT = np.empty(nStars)
    for j in range(nStars):
        color_separation = (J_Hobs[j]-jhMod)**2+(H_Kobs[j]-hkMod)**2
        min_separation_ind = np.argmin(color_separation)
        starsT[j] = teffMod[min_separation_ind]

    radeg = 180/np.pi
    niriss_pixel_scale = 0.065  # arcsec
    sweetSpot = dict(x=856, y=107, RA=allRA[targetIndex],
                     DEC=allDEC[targetIndex], jmag=Jmag[targetIndex])

    # offset between all stars and target
    dRA = (allRA - sweetSpot['RA'])*np.cos(sweetSpot['DEC']/radeg)*3600
    dDEC = (allDEC - sweetSpot['DEC'])*3600

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

    # dimX=256 #2048 #########now as argument, with default to 256
    dimY = 2048
    # cube of trace simulation at every degree of field rotation,
    # +target at O1 and O2
    simuCube = np.zeros([nPA+2, dimY, dimX])

    # saveFiles = glob.glob('idlSaveFiles/*.sav')[:-1]
    saveFiles = glob.glob(os.path.join(TRACES_PATH, 'NIRISS', '*.sav'))[:-1]
    # pdb.set_trace()

    # Big loop to generate a simulation at each instrument PA
    for kPA in range(PAtab.size):
        APA = PAtab[kPA]
        V3PA = APA+0.57  # from APT
        sindx = np.sin(np.pi/2+APA/radeg)*stars['dDEC']
        cosdx = np.cos(np.pi/2+APA/radeg)*stars['dDEC']
        nps = niriss_pixel_scale
        stars['dx'] = (np.cos(np.pi/2+APA/radeg)*stars['dRA']-sindx)/nps
        stars['dy'] = (np.sin(np.pi/2+APA/radeg)*stars['dRA']+cosdx)/nps
        stars['x'] = stars['dx']+sweetSpot['x']
        stars['y'] = stars['dy']+sweetSpot['y']

        # Display the star field (blue), target (red), subarray (green),
        # full array (blue), and axes
        if (kPA == 0 and nStars > 1) and False:
            plt.plot([0, 2047, 2047, 0, 0], [0, 0, 2047, 2047, 0], 'b')
            plt.plot([0, 255, 255, 0, 0], [0, 0, 2047, 2047, 0], 'g')
            # the order 1 & 2 traces
            path = '/Users/david/Documents/work/jwst/niriss/soss/data/'
            t1 = np.loadtxt(path+'trace_order1.txt', unpack=True)
            plt.plot(t1[0], t1[1], 'r')
            t2 = np.loadtxt(path+'trace_order2.txt', unpack=True)
            plt.plot(t2[0], t2[1], 'r')
            plt.plot(stars['x'], stars['y'], 'b*')
            plt.plot(sweetSpot['x'], sweetSpot['y'], 'r*')
            plt.title("APA= {} (V3PA={})".format(APA, V3PA))
            ax = plt.gca()

            # add V2 & V3 axes
            l, hw, hl = 250, 50, 50
            adx, ady = -l*np.cos(-0.57 / radeg), -l*np.sin(-0.57 / radeg)
            ax.arrow(2500, 1800, adx, ady, head_width=hw, head_length=hl,
                     length_includes_head=True, fc='k')  # V3
            plt.text(2500+1.4*adx, 1800+1.4*ady, "V3", va='center',
                     ha='center')
            adx, ady = -l*np.cos((-0.57-90)/radeg), -l*np.sin((-0.57-90)/radeg)
            ax.arrow(2500, 1800, adx, ady, head_width=hw, head_length=hl,
                     length_includes_head=True, fc='k')  # V2
            plt.text(2500+1.4*adx, 1800+1.4*ady, "V2", va='center',
                     ha='center')
            # add North and East
            adx, ady = -l*np.cos(APA/radeg), -l*np.sin(APA/radeg)
            ax.arrow(2500, 1300, adx, ady, head_width=hw, head_length=hl,
                     length_includes_head=True, fc='k')  # N
            plt.text(2500+1.4*adx, 1300+1.4*ady, "N", va='center', ha='center')
            adx, ady = -l*np.cos((APA-90)/radeg), -l*np.sin((APA-90)/radeg)
            ax.arrow(2500, 1300, adx, ady, head_width=hw, head_length=hl,
                     length_includes_head=True, fc='k')  # E
            plt.text(2500+1.4*adx, 1300+1.4*ady, "E", va='center', ha='center')

            ax.set_xlim(-400, 2047+800)
            ax.set_ylim(-400, 2047+400)
            ax.set_aspect('equal')
            plt.show()

        # Retain stars that are within the Direct Image NIRISS POM FOV
        ind, = np.where((stars['x'] >= -162) & (stars['x'] <= 2047+185) &
                        (stars['y'] >= -154) & (stars['y'] <= 2047+174))
        starsInFOV = stars[ind]

        for i in range(len(ind)):
            intx = round(starsInFOV['dx'][i])
            inty = round(starsInFOV['dy'][i])

            k = np.where(teffMod == starsInFOV['T'][i])[0][0]

            fluxscale = 10.0**(-0.4*(starsInFOV['jmag'][i]-sweetSpot['jmag']))

            # deal with subection sizes
            mx0 = int(modelPadX-intx)
            mx1 = int(modelPadX-intx+dimX)
            my0 = int(modelPadY-inty)
            my1 = int(modelPadY-inty+dimY)

            if (mx0 > dimXmod) or (my0 > dimYmod):
                continue
            if (mx1 < 0) or (my1 < 0):
                continue

            x0 = (mx0 < 0)*(-mx0)
            y0 = (my0 < 0)*(-my0)
            mx0 *= (mx0 >= 0)
            mx1 = dimXmod if mx1 > dimXmod else mx1
            my0 *= (my0 >= 0)
            my1 = dimYmod if my1 > dimYmod else my1

            # if target and first kPA, add target traces of order 1 and 2
            # in output cube
            if (intx == 0) & (inty == 0) & (kPA == 0):
                fNameModO12 = saveFiles[k]
                modelO12 = readsav(fNameModO12, verbose=False)['modelo12']
                ord1 = modelO12[0, my0:my1, mx0:mx1]*fluxscale
                ord2 = modelO12[1, my0:my1, mx0:mx1]*fluxscale
                simuCube[0, y0:y0+my1-my0, x0:x0+mx1-mx0] = ord1
                simuCube[1, y0:y0+my1-my0, x0:x0+mx1-mx0] = ord2

            if (intx != 0) or (inty != 0):
                mod = models[k, my0:my1, mx0:mx1]
                simuCube[kPA+2, y0:y0+my1-my0, x0:x0+mx1-mx0] += mod*fluxscale
                #mod = models[k, 1500:1500+2048, 100:100+256]
                #simuCube[kPA+2, 0:dimY, 0:dimX] += mod*fluxscale
    return simuCube

def fieldSim(ra, dec, instrument, binComp=''):
    """Produce a field simulation for a target
    Parameters
    ----------
    ra : float
        The RA of the target.
    dec : float
        The Dec of the target.
    instrument : str
        The instrument the contamination is being calculated for.
        Can either be (case-sensitive):
        'NIRISS', 'NIRCam', 'MIRI', 'NIRSpec'
    binComp : sequence
        The parameters of a binary companion.

    Returns
    -------
    simuCube : np.ndarray
        The simulated data cube. Index 0 and 1 (axis=0) show the trace of
        the target for orders 1 and 2 (respectively). Index 2-362 show the trace
        of the target at every position angle (PA) of the instrument.
    """

    # Calling the variables which depend on what instrument you use
    if instrument=='NIRISS':
        simuCube = sossFieldSim(ra, dec)

    else:
        if instrument=='NIRCam F444W':
            dimX = 51
            dimY = 1343 # <- Q: should be conservative w/ sub-array size?
            rad = 2.5
            pixel_scale = 0.063 # arsec
            xval, yval = 1096.9968649303112, 34.99693173255946 # got from PYSIAF
            add_to_apa = 0.0265 # got from jwst_gtvt/find_tgt_info.py
        elif instrument=='NIRCam F322W2':
            dimX = 51
            dimY = 1823 # <- Q: should be conservative w/ sub-array size?
            rad = 2.5
            pixel_scale = 0.063 # arsec
            xval, yval = 468.0140991987737, 35.007956285677665
            add_to_apa = 0.0265 # got from jwst_gtvt/find_tgt_info.py
        elif instrument=='MIRI':
            dimX = 55
            dimY = 427
            rad = 2.5
            pixel_scale = 0.11 # arsec
            xval, yval = 38.5, 829.0
            add_to_apa = 4.83425324

        # stars in large field around target
        targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg))
        targetRA = targetcrd.ra.value
        targetDEC = targetcrd.dec.value
        info = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone',
                                 radius=rad*u.arcmin)

        # coordinates of all stars in FOV, including target
        allRA = info['ra'].data.data
        allDEC = info['dec'].data.data
        Jmag = info['j_m'].data.data
        Hmag = info['h_m'].data.data
        Kmag = info['k_m'].data.data
        J_Hobs = Jmag-Hmag
        H_Kobs = Hmag-Kmag

        # target coords
        aa = ((targetRA-allRA)*np.cos(targetDEC))
        distance = np.sqrt(aa**2 + (targetDEC-allDEC)**2)
        targetIndex = np.argmin(distance)  # the target

        # add any missing companion
        if binComp != '':
            deg2rad = np.pi/180
            bb = binComp[0]/3600/np.cos(allDEC[targetIndex]*deg2rad)
            allRA = np.append(allRA, (allRA[targetIndex] + bb))
            allDEC = np.append(allDEC, (allDEC[targetIndex] + binComp[1]/3600))
            Jmag = np.append(Jmag, binComp[2])
            Hmag = np.append(Kmag, binComp[3])
            Kmag = np.append(Kmag, binComp[4])
            J_Hobs = Jmag-Hmag
            H_Kobs = Hmag-Kmag

        # number of stars
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
        #teffMod = modelParam['teffmod']
        teffMod = np.linspace(2000, 6000, 41)

        # find/assign Teff of each star
        starsT = np.empty(nStars)
        for j in range(nStars):
            color_separation = (J_Hobs[j]-jhMod)**2+(H_Kobs[j]-hkMod)**2
            min_separation_ind = np.argmin(color_separation)
            starsT[j] = teffMod[min_separation_ind]

        radeg = 180/np.pi
        sweetSpot = dict(x=xval, y=yval, RA=allRA[targetIndex],
                         DEC=allDEC[targetIndex], jmag=Jmag[targetIndex])

        # offset between all stars and target
        dRA = (allRA - sweetSpot['RA'])*np.cos(sweetSpot['DEC']/radeg)*3600
        dDEC = (allDEC - sweetSpot['DEC'])*3600

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

        # cube of trace simulation at every degree of field rotation,
        # +target at O1 and O2

        if instrument == 'NIRISS':
            simuCube = np.zeros([nPA+2, dimY+1, dimX+1])
            savFiles = glob.glob(os.path.join(TRACES_PATH, instrument.replace(' ', '_'), 'modelOrder*.sav'))
            savFiles = np.sort(savFiles) # **organizes the elements from lowest to
                                           # highest temp (2000->6000K)
        elif instrument == 'NIRCam F444W':
            simuCube = np.zeros([nPA+1, dimY+1, dimX+1])
            fitsFiles = glob.glob(os.path.join(TRACES_PATH, instrument.replace(' ', '_'), 'o1*.0.fits'))
            fitsFiles = np.sort(fitsFiles)
            print(TRACES_PATH)

        elif instrument == 'NIRCam F322W2':
            simuCube = np.zeros([nPA+1, dimY+1, dimX+1])
            fitsFiles = glob.glob(os.path.join(TRACES_PATH, instrument.replace(' ', '_'), 'o1*.0.fits'))
            fitsFiles = np.sort(fitsFiles)

        elif instrument == 'MIRI':
            simuCube = np.zeros([nPA+1, dimY+1, dimX+1])
            print(TRACES_PATH)
            fitsFiles = glob.glob(os.path.join(TRACES_PATH, instrument, '_*.fits'))
            fitsFiles = np.sort(fitsFiles)


        # Big loop to generate a simulation at each instrument PA
        for kPA in range(PAtab.size):
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print('ANGLE: ', kPA)
            APA = PAtab[kPA] # Aperture Position Angle (PA of instrument)
            V3PA = APA+add_to_apa  # from APT (PA of V3 axis)
            sindx = np.sin(np.pi/2+APA/radeg)*stars['dDEC']
            cosdx = np.cos(np.pi/2+APA/radeg)*stars['dDEC']
            ps = pixel_scale
            stars['dx'] = (np.cos(np.pi/2+APA/radeg)*stars['dRA']-sindx)/ps
            stars['dy'] = (np.sin(np.pi/2+APA/radeg)*stars['dRA']+cosdx)/ps
            stars['x'] = stars['dx']+sweetSpot['x']
            stars['y'] = stars['dy']+sweetSpot['y']

            # Retain stars that are within the Direct Image NIRISS POM FOV
            # I think this extends the subarray edges to the detector edges.
            # like it keeps the stars that fall out of the subarray but still
            # fall into the detector.
            #ind, = np.where((stars['x'] >= -162) & (stars['x'] <= dimY+185) &
            #                (stars['y'] >= -154) & (stars['y'] <= dimY+174))
            ind, = np.where((stars['x'] >= -8000) & (stars['x'] <= dimY+8000) &
                            (stars['y'] >= -8000) & (stars['y'] <= dimY+8000))
            starsInFOV = stars[ind]

            #print('INDDDDD')
            #print(ind)
            #print(stars['x'])

            for i in range(len(ind)):
                # are these the coordinates of the stars
                # on the detector (w/ undeviated wavelength)?
                intx = round(starsInFOV['dx'][i])
                inty = round(starsInFOV['dy'][i])

                # **this indexing assumes that teffMod is
                # sorted the same way fitsFiles was sorted
                k = np.where(teffMod == starsInFOV['T'][i])[0][0]


                fluxscale = 10.0**(-0.4*(starsInFOV['jmag'][i]-sweetSpot['jmag']))

                # deal with subection sizes
                if instrument!='NIRISS':
                    modelPadX = 0
                    modelPadY = 0
                else:
                    print(modelPadX, modelPadY)

                mx0 = int(modelPadX-intx)
                mx1 = int(modelPadX-intx+dimX)
                my0 = int(modelPadY-inty)
                my1 = int(modelPadY-inty+dimY)
                print('intx,y ',intx, inty)

                # this is where MIRI breaks
                if (mx0 > dimX) or (my0 > dimY):
                    print('broken here')
                    print(mx0, dimX, my0, dimY)
                    continue
                if (mx1 < 0) or (my1 < 0):
                    print('or broken here')
                    print(mx1, my1)
                    continue

                print('made it here')
                print(mx0, dimX, my0, dimY)
                print(mx1, my1)

                x0 = (mx0 < 0)*(-mx0)
                y0 = (my0 < 0)*(-my0)
                mx0 *= (mx0 >= 0)
                mx1 = dimX if mx1 > dimX else mx1
                my0 *= (my0 >= 0)
                my1 = dimY if my1 > dimY else my1

                #print('y0 ', y0)
                #print('my1 ', my1)
                #print('my0 ', my0)
                #print('x0 ', x0)
                #print('mx1 ', mx1)
                #print('mx0 ', mx0)
                # if target and first kPA, add target traces of order 1 and 2
                # in output cube
                # the target will have intx = 0, inty = 0

                if (intx == 0) & (inty == 0) & (kPA == 0):
                    print('fuuuuuuuuuuu')
                    if instrument=='NIRISS':
                        fNameModO12 = savFiles[k]
                        modelO12 = readsav(fNameModO12, verbose=False)['modelo12']
                        ord1 = modelO12[0, my0:my1, mx0:mx1]*fluxscale
                        ord2 = modelO12[1, my0:my1, mx0:mx1]*fluxscale
                        simuCube[0, y0:y0+my1-my0, x0:x0+mx1-mx0] = ord1
                        simuCube[1, y0:y0+my1-my0, x0:x0+mx1-mx0] = ord2

                    else:
                        fNameModO12 = fitsFiles[k]
                        modelO1 = fits.getdata(fNameModO12, 1)
                        ord1 = modelO1[0, my0:my1, mx0:mx1]*fluxscale
                        simuCube[0, y0:y0+my1-my0, x0:x0+mx1-mx0] = ord1



                if (intx != 0) or (inty != 0):
                    print('this is a milestone')
                    if instrument=='NIRISS':
                        mod = models[k, my0:my1, mx0:mx1]
                        simuCube[kPA+2, y0:y0+my1-my0, x0:x0+mx1-mx0] += mod*fluxscale
                    else:
                        print("MIRI IS WORKING???????? ")
                        fNameModO12 = fitsFiles[k]
                        modelO12 = fits.getdata(fNameModO12)
                        simuCube[kPA+1, y0:y0+my1-my0, x0:x0+mx1-mx0] += modelO12[0, my0:my1, mx0:mx1]*fluxscale

    return simuCube

def _fsnircam(ra, dec, instrument, binComp=''):
    """Produce a field simulation for a target
    Parameters
    ----------
    ra : float
        The RA of the target.
    dec : float
        The Dec of the target.
    instrument : str
        The instrument the contamination is being calculated for.
        Can either be (case-sensitive):
        'NIRISS', 'NIRCam', 'MIRI', 'NIRSpec'
    binComp : sequence
        The parameters of a binary companion.

    Returns
    -------
    simuCube : np.ndarray
        The simulated data cube. Index 0 and 1 (axis=0) show the trace of
        the target for orders 1 and 2 (respectively). Index 2-362 show the trace
        of the target at every position angle (PA) of the instrument.
    """

    # Calling the variables which depend on what instrument you use
    if instrument=='NIRISS':
        #dimX = 256
        #dimY = 2048
        dimX = 257 # columns in your simulation array (PandExo)
        dimY = 2301 # rows in your simulation array (PandExo)
        rad = 2.5 # radius (in arcmins) to query for neighboring stars
        pixel_scale = 0.065 # arcsec/pixel
        xval, yval = 856, 107 # <- Q: sweetSpot can be off the sub-array?
        add_to_apa = 0.57
    elif instrument=='NIRCam F444W':
        dimX = 51
        dimY = 1343 # <- Q: should be conservative w/ sub-array size?
        rad = 2.5
        pixel_scale = 0.065
        xval, yval = 16, 16
        add_to_apa = 0.0265 # got from jwst_gtvt/find_tgt_info.py
    elif instrument=='NIRCam F322W2':
        dimX = 51
        dimY = 1823 # <- Q: should be conservative w/ sub-array size?
        rad = 2.5
        pixel_scale = 0.065
        xval, yval = 16, 16
        add_to_apa = 0.0265 # got from jwst_gtvt/find_tgt_info.py
    #elif instrument=='MIRI':
    #    dimX = 256
    #    dimY =
    #    rad =
    #    pixel_scale =
    #    xval, yval =
    #    add_to_apa =
    #elif instrument=='NIRSpec':
    #    dimX = 256
    #    dimY =
    #    rad =
    #    pixel_scale = 0.065
    #    xval, yval =
    #    add_to_apa =

    # stars in large field around target
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg))
    targetRA = targetcrd.ra.value
    targetDEC = targetcrd.dec.value
    info = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone',
                             radius=rad*u.arcmin)

    # coordinates of all stars in FOV, including target
    allRA = info['ra'].data.data
    allDEC = info['dec'].data.data
    Jmag = info['j_m'].data.data
    Hmag = info['h_m'].data.data
    Kmag = info['k_m'].data.data
    J_Hobs = Jmag-Hmag
    H_Kobs = Hmag-Kmag

    # target coords
    aa = ((targetRA-allRA)*np.cos(targetDEC))
    distance = np.sqrt(aa**2 + (targetDEC-allDEC)**2)
    targetIndex = np.argmin(distance)  # the target

    # add any missing companion
    if binComp != '':
        deg2rad = np.pi/180
        bb = binComp[0]/3600/np.cos(allDEC[targetIndex]*deg2rad)
        allRA = np.append(allRA, (allRA[targetIndex] + bb))
        allDEC = np.append(allDEC, (allDEC[targetIndex] + binComp[1]/3600))
        Jmag = np.append(Jmag, binComp[2])
        Hmag = np.append(Kmag, binComp[3])
        Kmag = np.append(Kmag, binComp[4])
        J_Hobs = Jmag-Hmag
        H_Kobs = Hmag-Kmag

    # number of stars
    nStars = allRA.size

    # Restoring model parameters
    modelParam = readsav(os.path.join(IDLSAVE_PATH, 'modelsInfo.sav'),
                         verbose=False)
    models = modelParam['models']
    modelPadX = modelParam['modelpadx']
    modelPadY = modelParam['modelpady']
    dimXmod = modelParam['dimxmod']
    dimYmod = modelParam['dimymod']
    jhMod = modelParam['jhmod']
    hkMod = modelParam['hkmod']
    #teffMod = modelParam['teffmod']
    teffMod = np.linspace(2000, 6000, 41)

    # find/assign Teff of each star
    starsT = np.empty(nStars)
    for j in range(nStars):
        color_separation = (J_Hobs[j]-jhMod)**2+(H_Kobs[j]-hkMod)**2
        min_separation_ind = np.argmin(color_separation)
        starsT[j] = teffMod[min_separation_ind]

    radeg = 180/np.pi
    sweetSpot = dict(x=xval, y=yval, RA=allRA[targetIndex],
                     DEC=allDEC[targetIndex], jmag=Jmag[targetIndex])

    # offset between all stars and target
    dRA = (allRA - sweetSpot['RA'])*np.cos(sweetSpot['DEC']/radeg)*3600
    dDEC = (allDEC - sweetSpot['DEC'])*3600

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

    # cube of trace simulation at every degree of field rotation,
    # +target at O1 and O2
    simuCube = np.zeros([nPA+2, dimY+1, dimX+1])

    mypath = '/user/jmedina/Forks/exoctk/exoctk/data/contam_visibility/traces/nircam_f322w2/'
    #mypath = '/user/jmedina/Forks/exoctk/exoctk/data/contam_visibility/traces/nircam_f444w/'
    #mypath = '/Users/jmedina/Desktop/pandexo_niriss_traces'
    fitsFiles = glob.glob(os.path.join(mypath, 'o1*.0.fits'))
    fitsFiles = np.sort(fitsFiles) # **organizes the elements from lowest to
                                   # highest temp (2000->6000K)
    #mypath = '/user/jfilippazzo/Modules/ExoCTK/TOR/idlSaveFiles'
    #fitsFiles = glob.glob(os.path.join(mypath, 'modelO*.sav'))
    #fitsFiles = np.sort(fitsFiles)

    # Big loop to generate a simulation at each instrument PA
    for kPA in range(PAtab.size):
        print('ANGLE: ', kPA)
        APA = PAtab[kPA]
        V3PA = APA+add_to_apa  # from APT
        sindx = np.sin(np.pi/2+APA/radeg)*stars['dDEC']
        cosdx = np.cos(np.pi/2+APA/radeg)*stars['dDEC']
        ps = pixel_scale
        stars['dx'] = (np.cos(np.pi/2+APA/radeg)*stars['dRA']-sindx)/ps
        stars['dy'] = (np.sin(np.pi/2+APA/radeg)*stars['dRA']+cosdx)/ps
        stars['x'] = stars['dx']+sweetSpot['x']
        stars['y'] = stars['dy']+sweetSpot['y']

        # Retain stars that are within the Direct Image NIRISS POM FOV
        ind, = np.where((stars['x'] >= -162) & (stars['x'] <= dimY+185) &
                        (stars['y'] >= -154) & (stars['y'] <= dimY+174))
        starsInFOV = stars[ind]

        for i in range(len(ind)):
            # are these the coordinates of the stars
            # on the detector (w/ undeviated wavelength)?
            intx = round(starsInFOV['dx'][i])
            inty = round(starsInFOV['dy'][i])

            # **this indexing assumes that teffMod is
            # sorted the same way fitsFiles was sorted
            k = np.where(teffMod == starsInFOV['T'][i])[0][0]
            print(starsInFOV['T'][i], teffMod)
            print(k)

            fluxscale = 10.0**(-0.4*(starsInFOV['jmag'][i]-sweetSpot['jmag']))

            # deal with subection sizes
            modelPadX = 0 # or maybe dimX-256 ?
            modelPadY = 0 # or maybe dimY-2048 ?
            mx0 = int(modelPadX-intx)
            mx1 = int(modelPadX-intx+dimX)
            my0 = int(modelPadY-inty)
            my1 = int(modelPadY-inty+dimY)
            print('intx,y ',intx, inty)

            if (mx0 > dimX) or (my0 > dimY):
                continue
            if (mx1 < 0) or (my1 < 0):
                continue

            x0 = (mx0 < 0)*(-mx0)
            y0 = (my0 < 0)*(-my0)
            mx0 *= (mx0 >= 0)
            mx1 = dimX if mx1 > dimX else mx1
            my0 *= (my0 >= 0)
            my1 = dimY if my1 > dimY else my1

            print('y0 ', y0)
            print('my1 ', my1)
            print('my0 ', my0)
            print('x0 ', x0)
            print('mx1 ', mx1)
            print('mx0 ', mx0)
            # if target and first kPA, add target traces of order 1 and 2
            # in output cube
            # the target will have intx = 0, inty = 0

            if (intx == 0) & (inty == 0) & (kPA == 0):
                fNameModO12 = fitsFiles[k]
                #print('GENERATE THIS TRACE O12 SEPARATELY IN PANDEXO: ', fNameModO12)
                print(fNameModO12)
                modelO1 = fits.getdata(fNameModO12, 1)
                #modelO2 = fits.getdata(fNameModO12, 2)
                ord1 = modelO1[0, my0:my1, mx0:mx1]*fluxscale
                #ord2 = modelO1[0, my0:my1, mx0:mx1]*fluxscale
                simuCube[0, y0:y0+my1-my0, x0:x0+mx1-mx0] = ord1
                #simuCube[1, y0:y0+my1-my0, x0:x0+mx1-mx0] = ord2

                print('fluxscale ', fluxscale)


            if (intx != 0) or (inty != 0):
                #mod = models[k, my0:my1, mx0:mx1]
                #simuCube[kPA+2, y0:y0+my1-my0, x0:x0+mx1-mx0] += mod*fluxscale
                fNameModO12 = fitsFiles[k]
                modelO12 = fits.getdata(fNameModO12)
                #modelO12 = readsav(fNameModO12)
                print('BLAH2')
                print(simuCube.shape, modelO12.shape)
                print(y0+my1-my0-y0, x0+mx1-mx0-x0)
                print(my1-my0, mx1-mx0)


                simuCube[kPA+1, y0:y0+my1-my0, x0:x0+mx1-mx0] += modelO12[0, my0:my1, mx0:mx1]*fluxscale

    return simuCube


if __name__ == '__main__':
    ra, dec = "04 25 29.0162", "-30 36 01.603"  # Wasp 79
    #sossFieldSim(ra, dec)
    fieldSim(ra, dec, instrument='NIRISS')

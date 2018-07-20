import pysynphot as S
import numpy as np


def colorMod():
    from scipy.interpolate import CubicSpline
    # Load the Filters
    # load the J filter
    JBandFilterWave, JBandFilterThru = np.loadtxt('2MASS-2MASS.J.dat').T
    # load the H filter
    HBandFilterWave, HBandFilterThru = np.loadtxt('2MASS-2MASS.H.dat').T
    # load the K filter
    KBandFilterWave, KBandFilterThru = np.loadtxt('2MASS-2MASS.Ks.dat').T

    # Interpolate the Filters onto the Stellar Model Wavelength Grid
    # we are going to use this to setup the filters in the SYNPHOT framework
    tempModel = S.Icat('phoenix', 5500, 0.0, 4.5)

    JBandFilterSpline = CubicSpline(JBandFilterWave, JBandFilterThru)
    HBandFilterSpline = CubicSpline(HBandFilterWave, HBandFilterThru)
    KBandFilterSpline = CubicSpline(KBandFilterWave, KBandFilterThru)

    # For some reason pysnphot stores the Vega spectrum with different
    # number of flux points
    w = S.Vega.wave
    # isolate model region in JBand
    JWaveVRange = (w >= JBandFilterWave.min()) * (w <= JBandFilterWave.max())
    # isolate model region in HBand
    HWaveVRange = (w >= HBandFilterWave.min()) * (w <= HBandFilterWave.max())
    # isolate model region in KBand
    KWaveVRange = (w >= KBandFilterWave.min()) * (w <= KBandFilterWave.max())

    f = S.Vega.flux
    VegaJFlux = np.trapz(f[JWaveVRange] * JBandFilterSpline(w[JWaveVRange]))
    VegaHFlux = np.trapz(f[HWaveVRange] * HBandFilterSpline(w[HWaveVRange]))
    VegaKFlux = np.trapz(f[KWaveVRange] * KBandFilterSpline(w[KWaveVRange]))

    # Now to compare the JHK filters to the PHOENIX models
    jmin = tempModel.wave >= JBandFilterWave.min()
    jmax = tempModel.wave <= JBandFilterWave.max()
    hmin = tempModel.wave >= HBandFilterWave.min()
    hmax = tempModel.wave <= HBandFilterWave.max()
    kmin = tempModel.wave >= KBandFilterWave.min()
    kmax = tempModel.wave <= KBandFilterWave.max()
    JWaveRange = jmin * jmax  # isolate model region in JBand
    HWaveRange = hmin * hmax  # isolate model region in HBand
    KWaveRange = kmin * kmax  # isolate model region in KBand

    # Loop over the 30 models
    nModels = 30
    Jmags = np.zeros(nModels)
    Hmags = np.zeros(nModels)
    Kmags = np.zeros(nModels)

    teff_list = [x for x in range(2800, 5500+100, 100)] + [5800] + [6000]

    # check if I did that write
    assert(len(teff_list) == nModels)

    for kt, teff in enumerate(teff_list):
        modelTeff = S.Icat('phoenix', teff, 0.0, 4.5)  # load the model
        # integrate the flux in the J bandpass
        jmod = modelTeff.flux[JWaveRange]
        jsp = JBandFilterSpline(tempModel.wave[JWaveRange])
        Jmags[kt] = -2.5*np.log10(np.trapz(jmod * jsp / VegaJFlux))
        # integrate the flux in the H bandpass
        hmod = modelTeff.flux[HWaveRange]
        hsp = HBandFilterSpline(tempModel.wave[HWaveRange])
        Hmags[kt] = -2.5*np.log10(np.trapz(hmod * hsp / VegaHFlux))
        # integrate the flux in the K bandpass
        kmod = modelTeff.flux[KWaveRange]
        ksp = KBandFilterSpline(tempModel.wave[KWaveRange])
        Kmags[kt] = -2.5*np.log10(np.trapz(kmod * ksp / VegaKFlux))

    jhmod = Jmags - Hmags
    hkmod = Hmags - Kmags

    return jhmod, hkmod, teff_list

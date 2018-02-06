from pylab import *;ion()
import pysynphot as S

from scipy.interpolate import CubicSpline

def colorMod():
	# Load the Filters
	JBandFilterWave, JBandFilterThru = np.loadtxt('2MASS-2MASS.J.dat').T # load the J filter
	HBandFilterWave, HBandFilterThru = np.loadtxt('2MASS-2MASS.H.dat').T # load the H filter
	KBandFilterWave, KBandFilterThru = np.loadtxt('2MASS-2MASS.Ks.dat').T # load the K filter

	# Interpolate the Filters onto the Stellar Model Wavelength Grid
	tempModel = S.Icat('phoenix', 5500, 0.0, 4.5) # we are going to use this to setup the filters in the SYNPHOT framework

	JBandFilterSpline = CubicSpline(JBandFilterWave, JBandFilterThru)
	HBandFilterSpline = CubicSpline(HBandFilterWave, HBandFilterThru)
	KBandFilterSpline = CubicSpline(KBandFilterWave, KBandFilterThru)

	# For some reason pysnphot stores the Vega spectrum with different number of flux points
	JWaveVRange= (S.Vega.wave >= JBandFilterWave.min()) * (S.Vega.wave <= JBandFilterWave.max()) # isolate model region in JBand
	HWaveVRange= (S.Vega.wave >= HBandFilterWave.min()) * (S.Vega.wave <= HBandFilterWave.max()) # isolate model region in HBand
	KWaveVRange= (S.Vega.wave >= KBandFilterWave.min()) * (S.Vega.wave <= KBandFilterWave.max()) # isolate model region in KBand

	VegaJFlux = np.trapz(S.Vega.flux[JWaveVRange] * JBandFilterSpline(S.Vega.wave[JWaveVRange]))
	VegaHFlux = np.trapz(S.Vega.flux[HWaveVRange] * HBandFilterSpline(S.Vega.wave[HWaveVRange]))
	VegaKFlux = np.trapz(S.Vega.flux[KWaveVRange] * KBandFilterSpline(S.Vega.wave[KWaveVRange]))

	# Now to compare the JHK filters to the PHOENIX models
	JWaveRange= (tempModel.wave >= JBandFilterWave.min()) * (tempModel.wave <= JBandFilterWave.max()) # isolate model region in JBand
	HWaveRange= (tempModel.wave >= HBandFilterWave.min()) * (tempModel.wave <= HBandFilterWave.max()) # isolate model region in HBand
	KWaveRange= (tempModel.wave >= KBandFilterWave.min()) * (tempModel.wave <= KBandFilterWave.max()) # isolate model region in KBand


	# Loop over the 30 models
	nModels   = 30
	Jmags     = np.zeros(nModels)
	Hmags     = np.zeros(nModels)
	Kmags     = np.zeros(nModels)

	teff_list = [x for x in range(2800,5500+100,100)] + [5800] + [6000]

	# check if I did that write
	assert(len(teff_list) == nModels)

	for kt, teff in enumerate(teff_list):
		modelTeff = S.Icat('phoenix', teff, 0.0, 4.5) # load the model
		Jmags[kt] = -2.5*np.log10(np.trapz(modelTeff.flux[JWaveRange] * JBandFilterSpline(tempModel.wave[JWaveRange]) / VegaJFlux)) # integrate the flux in the J bandpass
		Hmags[kt] = -2.5*np.log10(np.trapz(modelTeff.flux[HWaveRange] * HBandFilterSpline(tempModel.wave[HWaveRange]) / VegaHFlux)) # integrate the flux in the H bandpass
		Kmags[kt] = -2.5*np.log10(np.trapz(modelTeff.flux[KWaveRange] * KBandFilterSpline(tempModel.wave[KWaveRange]) / VegaKFlux)) # integrate the flux in the K bandpass

	jhmod = Jmags - Hmags
	hkmod = Hmags - Kmags
	
	return jhmod, hkmod, teff_list
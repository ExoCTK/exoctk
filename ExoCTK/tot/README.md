# PandExo_HST
PandExo: A community tool for transiting exoplanet science with JWST & HST

Similar to an exposure time calculator (ETC), PandExo is a transiting exoplanet noise simulator. It is based on Pandeia, the ETC for JWST, and has been expanded to include HST's WFC3 instrument.

The included Jupyter Notebook tutorial demonstrates how to use PandExo to:

    1. Optimize WFC3's NSAMP and SAMP-SEQ parameters,
    2. Predict transmission/emission spectrum uncertainties, and
    3. Determine the observation start window
    
for any system observed with HST/WFC3 using the G102 & G141 grisms. The current implementation scales the measured flux, variance, and exposure time values from previously-observed systems, computes the expected rms per spectrophotometric channel, and estimates the transit/eclipse depth error based on the anticipated number of valid in- and out-of-transit data points.  The uncertainty estimates depend on the orbital properties of the system, instrument configuration, and observation duration.  The code assumes Gaussian-distributed white noise and uniform uncertainties over the G102 and G141 grisms; both assumptions are consistent with published results.


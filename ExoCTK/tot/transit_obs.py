from ExoCTK import core
import numpy as np
import os, json
import pickle as pkl
import pkg_resources

class SetDefaultModes():
    """
    This class contains functionality for setting observing modes for exoplanet observations 
    This is NOT a complete set of observing modes. Instead, if offers a starting point for 
    choosing one instrument specification. There is one function for each instrument. 
    
    Possible inputs are: 
       "WFC3 G102"
       "WFC3 G141"
    """

    def __init__(self, inst):
        inst, mode = inst.lower().split(' ')
        self.instrument = inst
        self.config = mode
    
    def pick(self): 
        try: 
            return getattr(self, self.instrument)()
        except: 
            print("INVALID INSTRUMENT NAME")
            return 
                                   
    def wfc3(self):
        #wfc3_input
        json_file = pkg_resources.resource_filename('ExoCTK', 'data/tot/wfc3_input.json')
        with open(json_file) as data_file:
            pandeia_data = json.load(data_file)
            pandeia_data["configuration"]["instrument"]["disperser"] = self.config
        return pandeia_data

def load_exo_dict():
    """
    Function loads in empty exoplanet dictionary for pandexo input
    """
    pandexo_input = {
        "star":{
            "type" : "user or phoenix", 
            "starpath" : "file path",
            "w_unit": "Angs,cm,um,cm or Hz",
            "f_unit": "W/m2/um, FLAM, Jy, or erg/s/cm2/Hz",
            "mag": "magnitude",
            "ref_wave": "corresponding ref wave",
            "temp": "Only if phoenix, (in K)",
            "metal": "in log Fe/H",
            "logg": "cgs"
            },

        "planet" :{
	        "type": "user",
            "exopath" : "file path",
            "w_unit": "Angs,cm,um,cm or Hz",
            "f_unit": "rp/r* or fp/f*"
            },

        "observation": {
            "sat_level": "in % sat level",
            "transit_duration": "in seconds",
            "noccultations": "num transits",
            "wave_bin": "in micron",
            "fraction": "time in/out", 
            "noise_floor":"constant number or file name"
            }
    }
    print("Replace all inputs before calling run_pandexo.")
    return pandexo_input   
    
def load_mode_dict(inst):
    '''
    Small function to pull in correct instrument dictionary 
    '''
    return SetDefaultModes(inst).pick()

def run_pandexo(pandexo_input, pandeia_input, save_file=True, output_path=os.getcwd(), output_file = ''):
    """
    This program contains the functionality for running PandExo 
    
    Parameters
    ---------- 
    exo : exoplanet input dictionary 
    inst : instrument input dictionary
    save_file : default saves file output but you can spefiy no file output
    output_path : path for output files. default will print out files in current working directory 
    output_file : file name, default = single.p for single run 
    
    Returns
    -------
    result : For single run output will just be a single PandExo output dictionary
            
    Attributes
    ----------
    load_mode_dict : to call instrument dictionaries based off keys 
    run_inst_space : to run instrument parameters space in parallel 
    run_param_space : to run exoplanet parameter space in parallel
    """
    
    print("Running Single Case w/ User Instrument Dict.")

    #define the calculation we'll be doing 
    calculation = pandexo_input['calculation'].lower()
    
    if calculation == 'scale':
        hmag            = pandexo_input['star']['hmag']
        trdur           = pandexo_input['observation']['transit_duration']
        numTr           = pandexo_input['observation']['noccultations']
        nchan           = pandexo_input['observation']['nchan']
        scanDirection   = pandexo_input['observation']['scanDirection']
        norbits         = pandexo_input['observation']['norbits']
        schedulability  = pandexo_input['observation']['schedulability']
        disperser       = pandeia_input['configuration']['instrument']['disperser'].lower()
        subarray        = pandeia_input['configuration']['detector']['subarray'].lower()
        nsamp           = pandeia_input['configuration']['detector']['nsamp']
        samp_seq        = pandeia_input['configuration']['detector']['samp_seq']
        results         = wfc3_TExoNS(hmag, trdur, numTr, nchan, disperser, scanDirection, subarray, nsamp, samp_seq, norbits, schedulability)
    else:
        print("****HALTED: Unknown calculation: %s" % calculation)
        return
        
    if output_file == '':
        output_file = 'singlerun.p'
    if save_file: 
        pkl.dump(results, open(os.path.join(output_path,output_file), "wb"))
    
    return results
    
def wfc3_GuessNOrbits(trdur):
    '''
    Determine number of HST orbits needed for transit observation when not provided by the user.
    
    PARAMETERS
    ----------
    trdur           : float, transit duration in seconds
    
    RETURNS
    -------
    norbits         : float, number of HST orbits per transit (including discarded thermal-settling orbit)
    
    HISTORY
    -------
    Written by Kevin Stevenson      October 2016
    '''
    # Compute number of HST orbits during transit
    # ~96 minutes per HST orbit
    orbitsTr    = trdur/96./60.
    if orbitsTr <= 1.5:
        norbits = 4.
    elif orbitsTr <= 2.0:
        norbits = 5.
    else:
        norbits = np.ceil(orbitsTr*2+1)
    
    return norbits

def wfc3_GuessParams(hmag, disperser, scanDirection, subarray, obsTime, maxScanHeight=180., maxExptime=150.):
    '''
    Determine optimal nsamp and samp_seq values when not provided by the user.
    
    PARAMETERS
    ----------
    hmag            : float, H-band magnitude
    disperser       : string, grism (g102 or g141)
    scanDirection   : string, spatial scan direction (Forward or Round Trip)
    subarray        : string, subarray aperture (grism256 or grism512)
    obsTime         : float, available observing time per HST orbit in seconds
    maxScanHeight   : (optional) float, maximum scan height in pixels
    maxExptime      : (optional) float, maximum exposure time in seconds
    
    RETURNS
    -------
    bestnsamp       : float, number of up-the-ramp samples (1..15)
    bestsamp_seq    : string, time between non-destructive reads (SPARS5, SPARS10, or SPARS25)
    
    HISTORY
    -------
    Written by Kevin Stevenson      October 2016
    '''
    
    allnsamp    = np.arange(1,16)
    allsampseq  = ['spars5', 'spars10', 'spars25']
    maxDutyCycle= 0
    for samp_seq in allsampseq:
        for nsamp in allnsamp:
            exptime, tottime, scanRate, scanHeight, fluence = wfc3_obs(hmag, disperser, scanDirection, 
                                                                       subarray, nsamp, samp_seq)
            # Compute duty cycle and compare
            # Exposure time should be less than 2.5 minutes to achieve good time resolution
            ptsOrbit    = np.floor(obsTime/tottime)
            dutyCycle   = (exptime*(ptsOrbit-1))/50./60*100
            if (dutyCycle > maxDutyCycle) and (exptime < maxExptime) and (scanHeight < maxScanHeight):
                maxDutyCycle    = dutyCycle
                bestsampseq     = samp_seq
                bestnsamp       = nsamp
    
    return bestnsamp, bestsampseq

def wfc3_obs(hmag, disperser, scanDirection, subarray, nsamp, samp_seq):
    '''
    Determine the recommended exposure time, scan rate, scan height, and overheads.
    
    PARAMETERS
    ----------
    hmag            : float, H-band magnitude
    disperser       : string, grism (g102 or g141)
    scanDirection   : string, spatial scan direction (Forward or Round Trip)
    subarray        : string, subarray aperture (grism256 or grism512)
    nsamp           : float, number of up-the-ramp samples (1..15)
    samp_seq        : string, time between non-destructive reads (SPARS5, SPARS10, or SPARS25)
    
    RETURNS
    -------
    exptime         : float, exposure time in seconds
    tottime         : float, total frame time including overheads in seconds
    scanRate        : float, recommended scan rate in arcsec/s
    scanHeight      : float, scan height in pixels
    fluence         : float, maximum pixel fluence in electrons
    
    HISTORY
    -------
    Written by Kevin Stevenson      October 2016
    '''
    
    # Estimate exposure time
    if subarray == 'grism512':
        # GRISM512
        if samp_seq == 'spars5':
            exptime     = 0.853 + (nsamp-1)*2.9215  #SPARS5
        elif samp_seq == 'spars10':
            exptime     = 0.853 + (nsamp-1)*7.9217  #SPARS10
        elif samp_seq == 'spars25':
            exptime     = 0.853 + (nsamp-1)*22.9213 #SPARS25
        else:
            print("****HALTED: Unknown SAMP_SEQ: %s" % samp_seq)
            return
    else:
        # GRISM256
        if samp_seq == 'spars5':
            exptime     = 0.280 + (nsamp-1)*2.349   #SPARS5
        elif samp_seq == 'spars10':
            exptime     = 0.278 + (nsamp-1)*7.3465  #SPARS10
        elif samp_seq == 'spars25':
            exptime     = 0.278 + (nsamp-1)*22.346  #SPARS25
        else:
            print("****HALTED: Unknown SAMP_SEQ: %s" % samp_seq)
            return
    
    # Recommended scan rate
    scanRate    = np.round(1.9*10**(-0.4*(hmag-5.9)),3)     #arcsec/s
    if disperser == 'g102':
        # G102/G141 flux ratio is ~0.8
        scanRate *= 0.8
    # Max fluence in electrons/pixel
    fluence     = (5.5/scanRate)*10**(-0.4*(hmag-15))*2.4   #electrons
    if disperser == 'g102':
        # WFC3_ISR_2012-08 states that the G102/G141 scale factor is 0.96 DN/electron
        fluence *= 0.96
        # G102/G141 flux ratio is ~0.8
        fluence *= 0.8
    # Scan height in pixels
    scanRatep   = scanRate/0.121                            #pixels/s
    scanHeight  = scanRatep*exptime                         #pixels
    '''
    #Quadratic correlation between scanRate and read overhead
    foo     = np.array([0.0,0.1,0.3,0.5,1.0,2.0,3.0,4.0])/0.121
    bar     = np.array([ 40, 40, 41, 42, 45, 50, 56, 61])
    c       = np.polyfit(foo,bar,2)
    model   = c[2] + c[1]*foo + c[0]*foo**2
    #c = [  6.12243227e-04,   6.31621064e-01,   3.96040946e+01]
    '''
    # Define instrument overheads (in seconds)
    c       = [6.12243227e-04, 6.31621064e-01, 3.96040946e+01]
    read    = c[2] + c[1]*scanRatep + c[0]*scanRatep**2
    # Correlation between scanHeight/scanRate and pointing overhead was determined elsewhere
    if scanDirection == 'Round Trip':
        # Round Trip scan direction doesn't have to return to starting point, therefore no overhead
        pointing    = 0.
    elif scanDirection == 'Forward':
        c           = [  3.18485340e+01,   3.32968829e-02,   1.65687590e-02,   
                         7.65510038e-01,  -6.24504499e+01,   5.51452028e-03]
        pointing    = c[0]*(1 - np.exp(-c[2]*(scanHeight-c[4]))) + c[1]*scanHeight + c[3]*scanRatep +c[5]*scanRatep**2
    else:
        print("****HALTED: Unknown scan direction: %s" % scanDirection)
        return
    # Estimate total frame time including overheads
    tottime     = exptime+read+pointing         #seconds
    
    return exptime, tottime, scanRate, scanHeight, fluence

def wfc3_TExoNS(hmag, trdur, numTr, nchan, disperser, scanDirection='Forward', subarray='grism256', nsamp=0, samp_seq=None, norbits=None, schedulability='100'):
    '''
    Compute transit/eclipse depth uncertainty for a defined system and number of spectrophotometric channels.
    
    PARAMETERS
    ----------
    hmag            : float, H-band magnitude
    trdur           : array, transit duration in seconds
    numTr           : array, number of observed transits/eclipses
    nchan           : float, number of spectrophotometric channels
    disperser       : string, grism (g102 or g141)
    scanDirection   : string, spatial scan direction (Forward or Round Trip)
    subarray        : string, subarray aperture (grism256 or grism512)
    nsamp           : float, number of up-the-ramp samples (1..15)
    samp_seq        : string, time between non-destructive reads (SPARS5, SPARS10, or SPARS25)
    norbits         : float, number of requested orbits per transit (including discarded thermal-settling orbit)
    schedulability  : string, % time HST can observe target (30 or 100)
    
    RETURNS
    -------
    deptherr        : float, transit/eclipse depth uncertainty per spectrophotometric channel
    
    HISTORY
    -------
    Written by Kevin Stevenson      October 2016
    '''
    
    disperser   = disperser.lower()
    subarray    = subarray.lower()
    if isinstance(samp_seq, str):
        samp_seq    = samp_seq.lower()
    
    if disperser == 'g141':
        # Define reference Hmag, flux, variance, and exposure time for GJ1214
        refmag      = 9.094
        refflux     = 2.32e8
        refvar      = 2.99e8
        refexptime  = 88.436
    elif disperser == 'g102':
        # Define reference Hmag, flux, variance, and exposure time for WASP12
        refmag      = 10.228
        refflux     = 8.26e7
        refvar      = 9.75e7
        refexptime  = 103.129
    else:
        print("****HALTED: Unknown disperser: %s" % disperser)
        return
    
    # Define maximum recommended scan height
    if subarray == 'grism512':
        maxScanHeight = 430
    elif subarray == 'grism256':
        maxScanHeight = 180
    else:
        print("****HALTED: Unknown subarray aperture: %s" % subarray)
        return
    
    # Define maximum frame time
    maxExptime  = 150.
    
    # Define available observing time per HST orbit in seconds
    if schedulability == '30':
        obsTime   = 51.3*60
    elif schedulability == '100':
        obsTime   = 46.3*60
    else:
        print("****HALTED: Unknown schedulability: %s" % schedulability)
        return
    
    # Compute recommended number of HST orbits and compare to user specified value
    guessorbits = wfc3_GuessNOrbits(trdur)
    if norbits == None:
        norbits = guessorbits
    elif norbits != guessorbits:
        print("****WARNING: Number of specified HST orbits does not match number of recommended orbits: %0.0f" % guessorbits)
    
    if nsamp == 0 or samp_seq == None or samp_seq == "None":
        # Determine optimal nsamp and samp_seq values
        nsamp, samp_seq = wfc3_GuessParams(hmag, disperser, scanDirection, subarray, 
                                           obsTime, maxScanHeight, maxExptime)
    # Calculate observation parameters
    exptime, tottime, scanRate, scanHeight, fluence = wfc3_obs(hmag, disperser, scanDirection, 
                                                               subarray, nsamp, samp_seq)
    if scanHeight > maxScanHeight:
        print("****WARNING: Computed scan height exceeds maximum recommended height of %0.0f pixels." % maxScanHeight)
    if exptime > maxExptime:
        print("****WARNING: Computed frame time (%0.0f seconds) exceeds maximum recommended duration of %0.0f seconds." % (exptime,maxExptime))
    
    # Compute number of data points (frames) per orbit
    ptsOrbit    = np.floor(obsTime/tottime)
    # First point (frame) is always low, ignore when computing duty cycle
    dutyCycle   = (exptime*(ptsOrbit-1))/50./60*100
    
    # Compute number of non-destructive reads per orbit
    readsOrbit  = ptsOrbit*(nsamp+1)
    
    # Look for mid-orbit buffer dumps
    if (subarray == 'grism256') and (readsOrbit >= 300) and (exptime <= 43):
        print("****WARNING: Observing plan may incur mid-orbit buffer dumps.  Check with APT.")
    if (subarray == 'grism512') and (readsOrbit >= 120) and (exptime <= 100):
        print("****WARNING: Observing plan may incur mid-orbit buffer dumps.  Check with APT.")
    
    # Compute number of HST orbits per transit
    # ~96 minutes per HST orbit
    orbitsTr    = trdur/96./60.
    
    # Estimate number of good points during planet transit
    # First point in each HST orbit is flagged as bad; therefore, subtract from total
    if orbitsTr < 0.5:
        # Entire transit fits within one HST orbit
        ptsInTr = ptsOrbit * orbitsTr/0.5 - 1
    elif orbitsTr <= 1.5:
        # Assume one orbit centered on mid-transit time
        ptsInTr = ptsOrbit - 1
    elif orbitsTr < 2.:
        # Assume one orbit during full transit and one orbit during ingress/egress
        ptsInTr = ptsOrbit * (np.floor(orbitsTr) + np.min((1,np.remainder(orbitsTr-np.floor(orbitsTr)-0.5,1)/0.5))) - 2
    else:
        # Assume transit contains 2+ orbits timed to maximize # of data points.
        ptsInTr = ptsOrbit * (np.floor(orbitsTr) + np.min((1,np.remainder(orbitsTr-np.floor(orbitsTr),1)/0.5))) - np.ceil(orbitsTr)
    
    # Estimate number of good points outside of transit
    # Discard first HST orbit
    ptsOutTr    = (ptsOrbit-1) * (norbits-1) - ptsInTr
    
    # Compute transit depth uncertainty per spectrophotometric channel
    ratio       = 10**((refmag - hmag)/2.5)
    flux        = ratio*refflux*exptime/refexptime
    fluxvar     = ratio*refvar*exptime/refexptime
    chanflux    = flux/nchan
    chanvar     = fluxvar/nchan
    chanrms     = np.sqrt(chanvar)/chanflux*1e6     #ppm
    inTrrms     = chanrms/np.sqrt(ptsInTr*numTr)    #ppm
    outTrrms    = chanrms/np.sqrt(ptsOutTr*numTr)   #ppm
    deptherr    = np.sqrt(inTrrms**2 + outTrrms**2) #ppm
    
    print("Number of HST orbits: %0.0f" % norbits)
    print("WFC3 parameters: NSAMP = %0.0f, SAMP_SEQ = %s" %(nsamp,samp_seq.upper()))
    print("Recommended scan rate: %0.3f arcsec/s" % scanRate)
    print("Scan height: %0.1f pixels" % scanHeight)
    print("Maximum pixel fluence: %0.0f electrons" % fluence)
    print("Estimated duty cycle (outside of Earth occultation): %0.1f%%" % dutyCycle)
    print("Transit depth uncertainty: %0.1f ppm for each of %0.0f channel(s)" % (deptherr, nchan))
    
    return deptherr/1e6, chanrms/1e6, ptsOrbit

def calc_StartWindow(eventType, rms, ptsOrbit, numOrbits, depth, inc, aRs, period, windowSize, ecc=0, w=90., duration=None, offset=0.):
    '''
    Plot earliest and latest possible spectroscopic light curves for given start window size
    
    PARAMETERS
    ----------
    eventType       : str, 'transit' or 'eclipse'
    rms             : float, light curve root-mean-square
    ptsOrbit        : float, number of frames per HST orbit
    numOrbits       : float, number of HST orbits per visit
    depth           : float, transit/eclipse depth
    inc             : float, orbital inclination in degrees
    aRs             : float, Semi-major axis in units of stellar radii (a/R*)
    period          : float, orbital period in days
    windowSize      : float, observation start window size in minutes
    ecc             : (Optional) float, eccentricity (default is 0)
    w               : (Optional) float, longitude of periastron (default is 90 degrees)
    duration        : (Optional) float, full transit/eclipse duration in days
    offset          : (Optional) float, manual offset in observation start time, in minutes
    
    RETURNS
    -------
    minphase        : float, earliest observation start phase
    maxphase        : float, latest observation start phase
    
    HISTORY
    -------
    Written by Kevin Stevenson          October 2016
    Added eccentric orbit handling      December 2016
    '''
    import matplotlib.pyplot as plt
    import batman
    
    hstperiod   = 96./60/24                 # HST orbital period, in days
    punc        = windowSize/120./24/period # Half start window size, in phase
    cosi        = np.cos(inc*np.pi/180)     # Cosine of the inclination
    rprs        = np.sqrt(depth)            # Planet-star radius ratio
        
    params          = batman.TransitParams()
    if eventType == 'transit':
        midpt       = period
        b           = aRs*cosi*(1-ecc**2)/(1+ecc*np.sin(w*np.pi/180)) # Impact parameter
        sfactor     = np.sqrt(1-ecc**2)/(1+ecc*np.sin(w*np.pi/180))   # Account for planet speed on eccentric orbits
        params.u    = [0.1, 0.1]                                      # limb darkening coefficients
    elif eventType == 'eclipse':
        midpt       = period/2*(1+4*ecc*np.cos(w*np.pi/180)/np.pi)
        b           = aRs*cosi*(1-ecc**2)/(1-ecc*np.sin(w*np.pi/180)) # Impact parameter
        sfactor     = np.sqrt(1-ecc**2)/(1-ecc*np.sin(w*np.pi/180))   # Account for planet speed on eccentric orbits
        params.u    = [0.0, 0.0]                                      # limb darkening coefficients
    else:
        print("****HALTED: Unknown event type: %s" % eventType)
        return
    params.t0       = midpt/period          # phase of transit/eclipse
    params.per      = 1.                    # orbital period, units are orbital phase
    params.rp       = rprs                  # planet radius (in units of stellar radii)
    params.a        = aRs                   # semi-major axis (in units of stellar radii)
    params.inc      = inc                   # orbital inclination (in degrees)
    params.ecc      = ecc                   # eccentricity
    params.w        = w                     # longitude of periastron (in degrees)
    params.limb_dark= "quadratic"           # limb darkening model
    
    if duration == None:                    # Transit/eclipse duration
        duration = period/np.pi*np.arcsin(1./aRs*np.sqrt(((1+rprs)**2-(aRs*cosi)**2)/(1-cosi**2)))*sfactor
    phase1      = (midpt + duration/2. - hstperiod*(numOrbits-2) - hstperiod/2 + offset/24./60)/period
    phase2      = (midpt - duration/2. - hstperiod*2 + offset/24./60)/period
    minphase    = (phase1+phase2)/2-punc
    maxphase    = (phase1+phase2)/2+punc

    #Plot extremes of possible HST observations
    npts        = 4 * ptsOrbit * numOrbits
    phdur       = duration/period
    phase1      = np.linspace(minphase+hstperiod/period,minphase+hstperiod/period*(numOrbits-1)+hstperiod/period/2,npts)
    phase2      = np.linspace(maxphase+hstperiod/period,maxphase+hstperiod/period*(numOrbits-1)+hstperiod/period/2,npts)
    m           = batman.TransitModel(params, phase1)
    trmodel1    = m.light_curve(params)
    m           = batman.TransitModel(params, phase2)
    trmodel2    = m.light_curve(params)
    obsphase1   = []
    obsphase2   = []
    for i in range(numOrbits):
        obsphase1   = np.r_[obsphase1, np.linspace(minphase+hstperiod/period*i,minphase+hstperiod/period*i+hstperiod/period/2,ptsOrbit)]
        obsphase2   = np.r_[obsphase2, np.linspace(maxphase+hstperiod/period*i,maxphase+hstperiod/period*i+hstperiod/period/2,ptsOrbit)]
    m           = batman.TransitModel(params, obsphase1)
    obstr1      = m.light_curve(params) + np.random.normal(0, rms, obsphase1.shape)
    m           = batman.TransitModel(params, obsphase2)
    obstr2      = m.light_curve(params) + np.random.normal(0, rms, obsphase2.shape)
    
    plt.figure(None, figsize=(12,4))
    plt.clf()
    a=plt.subplot(211)
    a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%.4f'))
    plt.title('Earliest Start Phase', size=12)
    plt.errorbar(obsphase1, obstr1, rms, fmt='go')
    plt.plot(phase1, trmodel1, 'b-', lw=2)
    ylim1   = plt.ylim()
    xlim1   = plt.xlim()
    plt.ylabel("Normalized Flux", size=12)
    plt.xlabel("Orbital Phase", size=12)
    b=plt.subplot(212)
    b.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%.4f'))
    plt.title('Latest Start Phase', size=12)
    plt.errorbar(obsphase2, obstr2, rms, fmt='ro')
    plt.plot(phase2, trmodel2, 'b-', lw=2)
    ylim2   = plt.ylim()
    xlim2   = plt.xlim()
    plt.ylabel("Normalized Flux", size=12)
    plt.xlabel("Orbital Phase", size=12)
    #Put both subplots onto same x,y scale
    ylim    = [np.min((ylim1[0],ylim2[0])), np.max((ylim1[1],ylim2[1]))]
    xlim    = [np.min((xlim1[0],xlim2[0])), np.max((xlim1[1],xlim2[1]))]
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.subplot(211)
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.tight_layout()
    
    return minphase, maxphase

def plot_PlanSpec(specfile, w_unit, disperser, deptherr, nchan, smooth=None, labels=None):
    '''
    Plot exoplanet transmission/emission spectrum
    
    PARAMETERS
    ----------
    specfile        : string, filename for model spectrum [wavelength, flux]
    w_unit          : string, wavelength unit (um or nm)
    disperser       : string, grism (g102 or g141)
    deptherr        : float or array of length nchan, simulated transit/eclipse depth uncertainty
    nchan           : float, number of spectrophotometric channels
    smooth          : (Optional) float, length of smoothing kernel
    labels          : (Optional) Legend labels: ['Model name', 'Simulated obs. name']
        
    HISTORY
    -------
    Written by Kevin Stevenson      October 2016
    '''
    import matplotlib.pyplot as plt
    # Load model wavelengths and spectrum
    mwave, mspec = np.loadtxt(specfile, unpack=True)
    # Convert wavelength to microns
    if w_unit == 'um':
        pass
    elif w_unit == 'nm':
        mwave /= 1000.
    else:
        print("****HALTED: Unrecognized wavelength unit '%s'" % w_unit)
        return
    
    # Smooth model spectrum (optional)
    if smooth != None:
        mspec = core.smooth(mspec, smooth)
    
    # Determine disperser wavelength boundaries
    if disperser == 'g141':
        wmin = 1.125
        wmax = 1.650
    elif disperser == 'g102':
        wmin = 0.84
        wmax = 1.13
    else:
        print("WARNING: Unrecognized disperser name '%s'" % disperser)
    
    # Determine wavelength bins
    binsize     = (wmax - wmin)/nchan
    wave_low    = np.round([i for i in np.linspace(wmin, wmax-binsize, nchan)],3)
    wave_hi     = np.round([i for i in np.linspace(wmin+binsize, wmax, nchan)],3)
    binwave     = (wave_low + wave_hi)/2.
    
    # Create simulated spectrum by binning model spectrum and addding uncertainty
    binspec     = np.zeros(nchan)
    for i in range(nchan):
        ispec       = np.where((mwave >= wave_low[i])*(mwave <= wave_hi[i]))
        binspec[i]  = np.mean(mspec[ispec])
    binspec    += np.random.normal(0,deptherr,nchan)
    
    plt.figure(None, figsize=(12,4))
    plt.clf()
    plt.plot(mwave, mspec, '-k')
    plt.errorbar(binwave, binspec, deptherr, fmt='bo', ms=8)
    plt.xlim(wmin, wmax)
    plt.ylim(np.min(binspec)-2*deptherr, np.max(binspec)+2*deptherr)
    if labels is not None: plt.legend(labels, loc='upper left')
    plt.xlabel("Wavelength ($\mu m$)",size=12)
    plt.ylabel("Depth (ppm)",size=12)
    plt.tight_layout()
    
    return binspec

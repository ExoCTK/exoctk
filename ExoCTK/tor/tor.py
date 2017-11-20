## EXOCTK LEGGGGGGGO

"""This is essentially demystify_apt.py -- it aims to take in parameters about the transit
at hand and spit out parameters that can go straight (lol) into APT.

It has two main pieces:
    
    1. Convert the literal amount of time you need for a transit into groups and integrations.
    2. Make sure that the exposure isn't over-saturated. 

Author : 
    Jules Fowler, April 2017

Use : 
    Part of ExoCTK, but tbh tbd?

Outputs :
    Also tbd.
"""

## -- IMPORTS
import json
import math
import os

from astropy.io import ascii
import numpy as np
from scipy import interpolate
from scipy.integrate import quad

## -- FUNCTIONS
# - I think this is going to have a few different sections up in hereeee...

## -- PHOTOMETRY/THROUHGPUT ESTIMATES

def convert_mag_to_flux(mag):
    """Converts the magnitude to flux. Assumes ST zeropoints. 

    Parameters
    ----------
    mag : float
        The magnitude at a given lambda.
    
    Returns
    -------
    flux : float
        The flux at a given lambda.
    """
    
    flux = 10**((mag + 21.1)/(-2.5))
    return flux


def calc_cr(mag, band, temp, px_size, throughput, filt):
    """Calculates the maximum countrate for one pixel.

    Parameters
    ----------
    mag : float
        The magnitude throughout whichever band.
    band : str
        The magnitude band at hand. 
    temp : float
        The temperature of the object (in K.)
    px_size : float
        The pixels/m^2 for the given detector.
    throughput : str
        Whether to extimate with maximum (max), median (med),
        or mean (mean) throughput. 

    Returns
    -------
    countrate : float
        The countrate per pixel per second for the filter/detector.
    """
    
    # Assign band/filter wavelengths.
    bands, filters = create_band_filter_dicts()
    band_min, band_max, band_zpt = bands[band]
    
    filter_min, filter_max = filters[filt]['min_max']
    throughput_factor = filters[filt][throughput]

    # Calculate flux through given band -- based on input magnitude
    band_flux = ((10**((mag - band_zpt)/(-2.5)))*px_size)*(1e19)*(filter_max-filter_min)
    print('Band Flux : ' + str(band_flux))
    # Assign it a blackbody and transform the flux for the filter at hand.
    I_band = quad(estimate_blackbody, band_min, band_max, args=(temp))[0]
    
    I_filter = quad(estimate_blackbody, filter_min, filter_max, args=(temp))[0]
    print(I_filter)
    print(I_band)
    # Calculate the ration of filter to band flux
    I_ratio = I_filter/I_band
    
    print('Flux ratio : ' + str(I_ratio))
    filter_flux = I_ratio*band_flux
#    print('Flux through filter : ' + str(filter_flux))
    
    # Find out how much goes through the filter.
    countrate = throughput_factor*filter_flux
    print(countrate) 
    return countrate

def create_band_filter_dicts():
    """Create the dictionary structures that will come in handy
    for the other funcitons here.

    Returns
    -------
    bands : dict
        Dictionary full of min/max wavelengths for each band.
    filters : dict
        Nested dictionary full of min/max wavelengths and throughput vals.
    """
    
    # Initialize dicts
    bands, filters = {}, {}

    # Magnitude bands -- in angstrom
    bands['U'] = (3650-330, 3650+330, -20.94)
    bands['B'] = (4450-470, 4450+470, -20.45)
    bands['V'] = (5510-440, 5510+440, -21.12)
    bands['R'] = (6580-690, 6580+690, -21.61)
    bands['I'] = (8060-755, 8060+755, -22.27)
    bands['J'] = (12200-565, 12200+565, -23.80)
    bands['H'] = (16300-1535, 16300+1535, -24.80)
    bands['K'] = (21900-1950, 21900+1950, -26.00)
    bands['L'] = (34500-2360, 34500+2360, -27.87)
#    bands['M'] = (47500-2300, 47600+2300)
#    bands['N'] = (105000-12500, 105000+12500)
#    bands['Q'] = (210000-24900, 210000+24900)
    
    ## MIRI filters -- in angstrom
    filters['LRS'] = {'min_max': (50710.4818218, 118739.11626), 'max': 0.602303853487, 'med': 0.552979004444, 'mean': 0.504712002351}

    # NIRCam filters -- in angstrom
    filters['F322W2'] = {'min_max': (26245.490982, 40147.8957916), 'max': 0.572233552821, 'med': 0.479853694536, 'mean': 0.459312105383}
    filters['F444W'] = {'min_max': (38794.8296593, 49745.8917836), 'max': 0.545676241333, 'med': 0.486952905631, 'mean': 0.465809117579}
    filters['F277W'] = {'min_max': (25402.4048096, 31318.2364729), 'max': 0.41290462156, 'med': 0.321649982683, 'mean': 0.318531860247}
    
    # NIRISS
    filters['GR700XD'] = {'min_max': (9000, 23500), 'max': 0.35, 'med': 0.31, 'mean': 0.31}

    # NIRSpec
    filters['G140H'] = {'min_max': (10000.0, 12700.0), 'max': 0.744763372884, 'med': 0.627978608462, 'mean': 0.626210220963}
    filters['G140M'] = {'min_max': (10000.0, 12700.0), 'max': 0.789124407879, 'med': 0.685992244672, 'mean': 0.667456122386}
    filters['G235H'] = {'min_max': (17000.0, 30703.3), 'max': 0.750483467307, 'med': 0.657618898472, 'mean': 0.646282919026}
    filters['G235M'] = {'min_max': (17000.0, 31000.0), 'max': 0.80118579198, 'med': 0.692135889694, 'mean': 0.691814768153}
    filters['G395H'] = {'min_max': (29000.0, 51776.6), 'max': 0.863484315687, 'med': 0.76015947863, 'mean': 0.727039127476}
    filters['G395M'] = {'min_max': (29000.0, 51787.1), 'max': 0.921712341837, 'med': 0.82187746764, 'mean': 0.778701455106}

    return bands, filters
    
def estimate_blackbody(wavelength, temp):
    """Creates the integrand for a Planck blackbody.
    This is techically a self contained function that will calculate 
    the blackbody spectral radiance, but its purpose is be integrated.

    Parameters
    ----------
    wavelength : float
        The of the source (in microns.)
    temp : float
        The temperature of teh source (in K.)

    Returns
    -------
    integrand : float
        Final expression for spectral radiance.
    """
    # Constants
    h = 6.626e-34 # m^2 kg/s 
    k = 1.3806e-23 # m^2kg/s^2K
    c = 3e8 # m/s\
    # Convert wavelength to m and then frequency
    wavelength = wavelength*(1e-10)
    nu = c/(wavelength*1e6) # in Hz
    
    # Expression for blackbody
    integrand = (2*h*(nu**3)/(c**2))/(np.e**((h*nu)/(k*temp)) - 1)
    return integrand

## -- I GAVE UP AND KIND OF USE PANDEIA FUNCTIONS

def calc_groups_from_exp_time(max_exptime_per_int, t_frame):
    """
    Given the maximum saturation time, calculates the number
    of frames per group.

    Parameters
    ----------
    max_exptime_per_int : float
        The maximum number of seconds an integration can last
        before it's oversaturated. 
    t_frame : float
        The time per frame.
    """

    groups = max_exptime_per_int/t_frame
    return np.floor(groups)

def create_pandeia_dicts(infile):
    """
    WIP -- How to initialize these dicts quick w/o reading a 
    text file? RN the answer is hard coding...
    Time will tell.

    Returns
    -------
    pandeia_dicts : dict of dict of array
    """
    
    dat = ascii.read(infile)
    pandeia_dict = dict(dat)
    return pandeia_dict


def interpolate_from_dat(mag, ins, filt, sub, mod, band, t_frame, sat_lvl, infile, ta=False):
    """
    Interpolates the precalculated pandeia data to estimate the saturation limit.

    Parameters
    ----------
    mag : float
        The magnitude of the source. (Takes between 4.5-12.5)
    ins : str
       The instrument, allowable "miri", "niriss", "nirspec", "nircam". 
    filt : str
        Filter.
    sub : str   
        Subarray.
    mod : str
        Phoenix model key.
    band : str
        Magnitude band -- unsed rn.
    t_frame : float
        Frame time.
    sat_lvl : float 
        The maximum fullwell saturation we'll allow.
    infile : str
        The data file to use.
    ta : bool, optional
        Whether or not we're running this for TA.

    Returns
    -------
    n_group : int
        The number of groups that won't oversaturate the detector.
    """
    
    # Create the dictionaries for each filter and select out the prerun data
    with open(infile) as f:
        dat = json.load(f)
    
    ta_or_tor = 'tor'
    if ta:
        ta_or_tor = 'ta_tor'

    mags = dat['mags']
    sat = dat[ta_or_tor][ins][filt][sub][mod]
    
    # Inter/olate the given magnitude
    log_sat = np.log10(sat)
    func_log = interpolate.interp1d(mags, log_sat)
    func = interpolate.interp1d(mags, sat)
    max_log_sat = func_log(mag)
    max_sat_test = func(mag)
    max_sat = 10**(max_log_sat)

    # Figure out what it means in wake of the given sat lvl
    max_exptime = sat_lvl/max_sat

    # Calculate the nearest number of groups
    n_group = calc_groups_from_exp_time(max_exptime, t_frame)
    
    return n_group, max_sat


def map_to_ta_modes(ins, max_group, min_group):
    """Turns the min/max groups into the closest allowable
    TA group mode.

    Parameters
    ----------
    ins : str
        Instrument.
    max_group : int
        The maximum number of groups without oversaturating.
    
    min_group : int
        The groups needed to hit the target SNR.
    
    Returns
    -------
    ta_groups : int
        The suggestions for ta groups.
    """

    # Allowable group modes for each ins
    groups = {'miri': [3, 5, 9, 15, 23, 33, 45, 59, 75, 93, 113, 135, 159, 185, 243, 275, 513],
              'niriss': [3, 5, 7, 9, 1, 13, 15, 17, 19],
              'nirspec': [3], 
              'nircam': [3, 5, 9, 17, 33, 65]
              }
    
    # Split the diff and match 'em up. 
    avg_groups = np.mean([max_group, min_group])
    allowable_groups = groups[ins]
    ta_groups = min(allowable_groups, key=lambda x:abs(x-avg_groups))
    
    # Unless it was oversaturated from the get-go OR there aren't enough groups
    # for SNR
    if min_group == 0:
        ta_groups = 0
    if min_group > max(allowable_groups):
        ta_groups = -1
    return ta_groups

def min_groups(mag, ins, filt, sub, mod, band, infile):
    """Estimates the minimum number of groups to reach 
    target acq sat requirements.

    Parameters
    ----------
    mag : float 
        Magnitude of star.
    ins : str
        Instrument.
    filt : str
        Filter.
    sub : str
        Subarray.
    mod : str
        Phoenix model key.
    band : str, currently unused?
        The band -- right now only k sooo?
    infile : str
        The file with the pandeia data.
    
    Returns
    -------
    min_groups : int
        The minimum number of groups to reach target snr.
    """

    with open(infile) as f:
        dat = json.load(f)
    
    # Match to closest magnitude
    mags = dat['mags']
    closest_mag = min(mags, key=lambda x:abs(x-mag))
    print(closest_mag)
    index = mags.index(closest_mag)
    
    # Match to data 
    min_groups = dat['ta_snr'][ins][filt][sub][mod][index]
    print(min_groups)
    return min_groups


## -- SIMPLE LOGIC TIME CALCULATIONS

def calc_n_group(ins, countrate, sat_max, t_frame, n_frame, n_skip):
    """Calculates number of groups per integration.

    Parameters
    ----------
    countrate : float
        The counts per pixel per second for your source/detector/filter.
    sat_max : float
        The maximum saturation limit we'll reach. 
    t_frame : float
        The frame time (in seconds.)
    n_frame : int
        Frames per group -- always 1 except maybe brown dwarves.
    n_skip : int
        The skips per integration -- always 1 except maybe brown dwarves.
    
    Returns
    -------
    n_group : int
        The minimum integer of groups per integration that will NOT excede the 
        saturatin limit.
    """
    if ins == 'NIRCam':
        n_group = 4
        while (n_group*n_frame - n_skip)*(t_frame*countrate) > sat_max:
            n_group += -1
    
    if ins == 'NIRSpec':
        n_group = 5
        if compare_sat(n_group, countrate, sat_max, t_frame, n_frame, n_skip):
            n_group = 4
        if compare_sat(n_group, countrate, sat_max, t_frame, n_frame, n_skip):
            n_group = 1
        
    if ins == 'NIRISS':
        n_group = 30 
        while compare_sat(n_group, countrate, sat_max, t_frame, n_frame, n_skip):
            n_group += -1
        if n_group < 1:
            n_group = 1
            message = 'Your source may be too bright for this filter!'

    if ins == 'MIRI':
        n_group = ((sat_max)/(t_frame*countrate) + n_skip)/n_frame
        n_group = math.floor(n_group)
    
    return n_group
    
def compare_sat(n_group, countrate, sat_max, t_frame, n_frame, n_skip):
    
    print('Sat max ' + str(sat_max))
    print('Cr ' + str(countrate))
    print('frame # ' + str(n_frame))
    print('frame time ' + str(t_frame))

    comp = (n_group*n_frame - n_skip)*(t_frame*countrate) > sat_max
    print(comp)
    return comp

def calc_n_int(transit_time, n_group, n_reset, t_frame, n_frame):
    """Calculates number of integrations required.

    Parameters
    ----------
    transit_time : float
        The time of your planeting transit (in hours.)
    n_group : int
        Groups per integration.
    n_reset : int
        Reset frames per integration.
    t_frame : float
        The frame time (in seconds).
    n_frame : int
        Frames per group -- always 1 except maybe brown dwarves.

    Returns
    -------
    n_ints : float
        The required number of integraitions.
    """

    hour = 3600
    n_ints = (transit_time*hour)/(t_frame*(n_group*n_frame+n_reset))
    print('N frame : ' + str(n_frame))
    print('N reset : ' + str(n_reset))
    return math.ceil(n_ints)


def calc_obs_efficiency(t_exp, t_duration):
    """Calculates the observation efficiency.

    Parameters
    ----------
    t_exp : float
        Exposure time (in seconds).
    t_duration : float
        Duration time (in seconds).
    
    Returns
    -------
    obs_eff : float
        Observation efficiency.
    """

    obs_eff = t_exp/t_duration
    return obs_eff


def calc_t_duration(n_group, n_int, n_reset, t_frame, n_frame=1):
    """Calculates duration time (or exposure duration as told by APT.)

    Parameters
    ----------
    n_group : int
        Groups per integration.
    n_int : int
        Integrations per exposure.
    n_reset : int
        Reset frames per integration.
    t_frame : float
        Frame time (in seconds).
    n_frame : Frames per group -- always one except brown dwarves.

    Returns
    -------
    t_duration : float
        Duration time (in seconds).
    """
    
    t_duration = t_frame*(n_group*n_frame + n_reset)*n_int
    return t_duration


def calc_t_exp(n_int, t_ramp):
    """Calculates exposure time (or photon collection duration as told by APT.)

    Parameters
    ----------
    n_int : int
        Integrations per exposure.
    t_ramp : float
        Ramp time (in seconds).

    Returns
    -------
    t_exp : float
        Exposure time (in seconds).
    """

    t_exp = n_int*t_ramp
    return t_exp


def calc_t_frame(n_col, n_row, n_amp, ins):
    """Calculates the frame time for a given ins/readmode/subarray.

    Parameters
    ----------
    n_col : int
        Number of columns.
    n_row : int
        Number of rows.
    n_amp : int
        Amplifiers reading data.

    Returns:
    t_frame : float
        The frame time (in seconds).
    """
    n_col, n_amp, n_row = int(n_col), int(n_amp), int(n_row)
    
    if ins == 'nirspec':
        n = 2
    if ins in ['nircam', 'niriss']:
        n = 1

    t_frame = (n_col/n_amp + 12)*(n_row + n)*(1e-5)
    print(n_row)
    print(n_col)
    print(t_frame)
    return t_frame


def calc_t_int(n_group, t_frame, n_frame, n_skip):
    """Calculates the integration time.
    
    Parameters
    ----------
    n_group : int
        Groups per integration.
    t_frame : float
        Frame time (in seconds).
    n_frame : int
        Frames per group -- always 1 except maybe brown dwarves.
    n_skip : int 
        Skips per integration -- always 0 except maybe brown dwarves.

    Returns
    -------
    t_int : float
        Integration time (in seconds). 
    """
    
    t_int = (n_group*(n_frame + n_skip) - n_skip)*t_frame
    print(t_int)
    print(n_group)
    print(n_frame)
    print(n_skip)
    print(t_frame)
    return t_int


def calc_t_ramp(t_int, n_reset, t_frame):
    """Calculates the ramp time -- or the integration time plus overhead for resets.

    Parameters
    ----------
    t_int : float
        Integration time (in seconds).
    n_reset : int
        Rest frames per integration
    t_frame : float
        Frame time (in seconds).

    Returns
    -------
    t_ramp : float
        Ramp time (in seconds).
    """

    t_ramp = t_int + (n_reset - 1)*t_frame
    return t_ramp


def create_tor_dict(params, n_frame=1, n_skip=0):
    """Calculates all of the tor things and puts them in a dictionary for easy access. 

    Parameters
    ----------
    params : dict
        Dictionary of all the needed parameters. Must include:
        obs_time, n_group, mag, mod, band, filt, filt_ta,
        ins, subarray, subarray_ta, sat_mode, sat_max, infile
    n_frame :int, optional
        The number of frames -- almost always 1.
    n_skip: int, optional
        Number of skips -- almost always 0
    Returns
    -------
    tor_dict : dict, str
        Dictionary of values related to tor calculations. If 
        the calculation throws an error it will return a string 
        error message instead.

    Returns
    -------
    tor_dict : dict
        Dictionary of values related to tor calculations.
    """
    
    ## -- TARGET ACQ
    ta_frame_time = set_t_frame(params['ins'], params['subarray_ta'], ta=True)

    if ta_frame_time == 'Error!':
        return 'Looks like you mismatched your instrument and subarray. Please try again.'

    max_group, sat_rate_ta = interpolate_from_dat(params['mag'], params['ins'], params['filt_ta'], 
                                     params['subarray_ta'], params['mod'], params['band'],
                                     ta_frame_time, params['sat_max'], params['infile'], ta=True)
    min_group = min_groups(params['mag'], params['ins'], params['filt_ta'], params['subarray_ta'], 
                           params['mod'], params['band'], params['infile'])
    ta_groups = map_to_ta_modes(params['ins'], max_group, min_group)
    t_duration_ta = calc_t_duration(ta_groups, 1, 1, ta_frame_time)

    ## -- THE OTHER STUFF I GUESS??
    
    # Figure out the rows/cols/amps/px_size and sat
    ins_params = set_params_from_ins(params['ins'], params['subarray'])
    frame_time = set_t_frame(params['ins'], params['subarray'])

    # Test for an empty set of parameters
    if ins_params == 'Error!' or frame_time == 'Error':
        return 'Looks like you mismatched your instrument, filter, and subarray. Please try again.'
    
    # Or continue as normal
    else:
        n_row, n_col, n_amp, px_size, t_frame, n_reset = ins_params
        band_ins = '{0}_{1}'.format(params['band'], params['filt'])
        
        params['sat_max'] = convert_sat(params['sat_max'], params['sat_mode'], params['ins'])

        # Calculate countrate and n_groups if it isn't supplied
        n_group, sat_rate = interpolate_from_dat(params['mag'],
            params['ins'], params['filt'], params['subarray'], params['mod'],
            params['band'], frame_time, params['sat_max'], params['infile'])
        
        if str(params['n_group']) == 'optimize':
            params['n_group'] = int(n_group)
        else:
            params['n_group'] = int(float(params['n_group']))
        
        
        # Calculate times/ramps/etc
        t_int = calc_t_int(params['n_group'], frame_time, n_frame, n_skip)
        t_ramp = calc_t_ramp(t_int, n_reset, frame_time)
    
        # Calculate nubmer of integrations (THE MEAT)
        n_int = calc_n_int(params['obs_time'], params['n_group'], n_reset, frame_time, n_frame)
    
        # Other things that may come in handy who knows?
        t_exp = calc_t_exp(n_int, t_ramp)
        t_duration = calc_t_duration(params['n_group'], n_int, n_reset, frame_time, n_frame)
        obs_eff = calc_obs_efficiency(t_exp, t_duration)

        # Update params with new friends
        params['n_col'] = n_col
        params['n_row'] = n_row
        params['n_amp'] = n_amp
        params['n_reset'] = n_reset
        params['n_frame'] = n_frame
        params['n_skip'] = n_skip
        params['t_frame'] = round(frame_time, 3)
        params['t_int'] = round(t_int, 3)
        params['t_ramp'] = round(t_ramp, 3)
        params['n_int'] = n_int
        params['t_exp'] = round(t_exp/3600, 3)
        params['t_duration'] = round(t_duration/3600, 3)
        params['obs_eff'] = round(obs_eff, 3)
        params['ta_t_frame'] = ta_frame_time
        params['ta_groups'] = int(ta_groups)
        params['t_duration_ta'] = t_duration_ta
        params['max_sat_prediction'] = round(sat_rate*frame_time*params['n_group'], 3)
        params['max_sat_ta'] = round(sat_rate_ta*ta_frame_time*ta_groups, 3)
        return params 


## -- INS CONVERSION THINGS

def convert_sat(sat_max, sat_mode, ins):
     
    if sat_mode == 'well':
        ins_dict = {'nirspec': 65500, 'miri': 250000, 'nircam': 60000, 'niriss': 75000}
        sat_max = sat_max*ins_dict[ins]

    print(sat_max)
    return sat_max
    

def set_t_frame(ins, sub, ta=False):
    """ Assign the appropriate frame time based on the ins
    and subarray. For now, modes are implied.

    Parameters
    ----------
    ins : str
        The instrument : 'miri', 'niriss', 'nirspec', or 
        'nircam'.
    sub : str
        The subarray -- too lazy to write out the options 
        here.
    ta : bool   
        Whether this is for TA or not.

    Returns
    -------
    t_frame : float
        The frame time for this ins/sub combo.
    """
    
    # This dictionary holds all the frame times.
    frame_time = {'miri': {'ta': {'slitlessprism': 0.16}, 'slitlessprism': 0.16},
                  'niriss': {'ta': {'subtaami': 0.045}, 'substrip96': 2.219, 'substrip256': 5.49},
                  'nircam': {'ta': {'sub32tats': 0.015}, 
                             'full': 10.73666, 'subgrism256': 1.34666, 'subgrism128': 0.676, 'subgrism64': 0.3406},
                  'nirspec': {'ta': {'full': 10.736, 'sub32': 0.015},
                              'sub2048': 0.0916, 'sub1024a': 0.451, 'sub1024b': 0.451, 'sub512': 0.2255}
                  }
    try:

        if ta:
            print(ins, sub)
            t_frame = frame_time[ins]['ta'][sub]
        else:
            print(ins, sub)
            t_frame = frame_time[ins][sub]

    except (KeyError, IndexError) as e:
        print('RIP')
        t_frame = 'Error!'

    return t_frame


def set_params_from_ins(ins, subarray):
    """Sets/collects the running parameters from the instrument.
    
    Parameters
    ----------
    ins : str
        Instrument, options are nircam, niriss, nirpec, and miri.
    subarray : str
        Subarray mode.
    Returns
    -------
    rows : int  
        The number of pixels per row.
    cols : int
        The number of columns per row.
    """
    
    try:
        n_reset = 1
 
        if ins == 'nirspec':
            px_size = (40e-4)**2

            if subarray == 'sub2048':
                rows, cols = 2048, 32
            elif subarray in ['sub1024a', 'sub1024b']:
                rows, cols = 1024, 32
            elif subarray == 'sub512':
                rows, cols = 512, 32
            
            amps = 1 # 4 if not NRSRAPID????
            ft = calc_t_frame(cols, rows, amps, ins)
  
        elif ins == 'nircam':
            px_size = (18e-4)**2
            
            if subarray == 'full':
                rows, cols, amps = 2048, 2048, 4
            elif subarray == 'subgrism256':
                rows, cols, amps = 256, 256, 1
            elif subarray == 'subgrism128':
                rows, cols, amps = 128, 2048, 1
            elif subarray == 'subgrism64':
                rows, cols, amps = 64, 2048, 1
            
            ft = calc_t_frame(cols, rows, amps, ins)
  
        elif ins == 'miri':
            px_size = (25e-4)**2
            
            if subarray == 'slitlessprism':
                rows, cols, ft = 416, 72, .159
            
            amps = 4
            n_reset = 0 
            
        elif ins == 'niriss':
            px_size = (40e-4)**2
            
            if subarray == 'substrip96':
                rows, cols = 2048, 96
            elif subarray == 'substrip256':
                rows, cols = 2048, 256
            
            amps = 1
            ft = calc_t_frame(cols, rows, amps, ins)
        
        return rows, cols, amps, px_size, ft, n_reset
    
    except NameError:
        return 'Error!'

## -- RUN

if __name__ == "__main__":
    

    transit_time, n_col, n_row, n_reset, n_amp = 5, 96, 2048, 1, 1
    mag, band, temp, px_size, throughput, filt = 10.3, 'K',  2550, 3e9, 'med', 'G235M'
    sat_max = 20000

    tor_dict = create_tor_dict(transit_time, mag, band, temp, sat_max, px_size, throughput, filt, n_col, n_row, n_amp, n_reset)

    

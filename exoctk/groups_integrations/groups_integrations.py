"""
This is a module for calcuating groups and integrations with JWST on the fly. 
the main function (perform_calculation) takes a dictionary of inputs 
(modeled around how the web tool takes inputs) which must include : 
obs_time, n_group, mag, mod, band, filt, filt_ta, ins, subarray, subarray_ta, 
sat_mode, sat_max, and infile. It produces and dictionary of outputs that includes 
all of the original information as well as groups, integrations, saturation levels, 
and observation time estimates for target acquisition and science observations
with JWST. 


Authors
-------
    Jules Fowler, April 2017

Use
---
    This is mostly a module to be used by ExoCTKWeb -- but the main function
    can be run standalone.

"""

## -- IMPORTS
import json
import math
import os
from decimal import Decimal

from astropy.io import ascii
import numpy as np
from scipy import interpolate
from scipy.integrate import quad

## -- FUNCTIONS

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

    Returns
    -------
    groups : int
        The required number of groups.
    """

    groups = max_exptime_per_int/t_frame
    return np.floor(groups)


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
    n_ints = (float(transit_time)*hour)/(t_frame*(n_group*n_frame+n_reset))
    
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
    n_frame : int, optional
        Frames per group -- always one except brown dwarves.

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
    ins : str
        The instrument key.

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


def convert_sat(sat_max, sat_mode, ins, infile, ta=False):
    """ Converts full well fraction to a saturation in counts OR provides the 
    max fullwell for TA mode.

    Parameters
    ----------
    sat_max : float
        Either a full well fraction or counts.
    sat_mode : str
        'well' or 'counts'.
    ins : str
        The instrument.
    infile : str
        The path to the data file.
    ta : bool, optional
        Whether or not it's TA mode.
    
    Returns
    -------
    sat_max : float
        The fullwell to use in counts.
    """

    with open(infile) as f:
        dat = json.load(f)

    ins_dict = dat['fullwell'] 
    
    if sat_mode == 'well':
        sat_max = float(sat_max)*float(ins_dict[ins])
    
    if ta:
        sat_max = ins_dict[ins]
    
    return sat_max
    

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
    max_sat : int
        The maximum saturation level reached by that number of groups.
    """
    # Create the dictionaries for each filter and select out the prerun data
    with open(infile) as f:
        dat = json.load(f)

    ta_or_sci = 'sci_sat'

    if ta:
        ta_or_sci = 'ta_sat'

    # The data
    mags = np.array(dat['mags'])
    sat = dat[ta_or_sci][ins][filt][sub][mod]
    log_sat = np.log10(sat)

    # Interpolate the given magnitude
    func_log = interpolate.interp1d(mags, log_sat)
    max_log_sat = func_log(float(mag))
    max_sat = 10**(max_log_sat)

    # Figure out what it means in wake of the given sat lvl
    max_exptime = float(sat_lvl)/max_sat

    # Calculate the nearest number of groups
    n_group = calc_groups_from_exp_time(max_exptime, t_frame)

    # Can't have zero groups
    n_group = n_group or 1

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
    min_ta_groups : int
        The min possible groups to hit target SNR.
    max_ta_groups : int
        The max possible groups before saturation.
    """

    # Allowable group modes for each ins
    groups = {'miri': [3, 5, 9, 15, 23, 33, 45, 59, 75, 93, 113, 135, 159, 185, 243, 275, 513],
              'niriss': [3, 5, 7, 9, 1, 13, 15, 17, 19],
              'nirspec': [3], 
              'nircam': [3, 5, 9, 17, 33, 65]
              }
    
    # Match the literal min and max groups to the nearest mode. 
    allowable_groups = groups[ins]
    min_ta_groups = min(allowable_groups, key=lambda x:abs(x-min_group))
    max_ta_groups = min(allowable_groups, key=lambda x:abs(x-max_group))
    
    # Unless it was oversaturated from the get-go OR there aren't enough groups
    # for SNR
    if min_group == 0:
        min_ta_groups = 0
        max_ta_groups = 0
    if min_group > max(allowable_groups):
        min_ta_groups = -1
        max_ta_groups = 0
    
    # BOTH ARE FLIPPED RN -- I WILL FLIP BOTH BACK SOON...
    return max_ta_groups, min_ta_groups


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
    mags = [float(i) for i in dat['mags']]
    closest_mag = min(mags, key=lambda x:abs(x-float(mag)))
    index = mags.index(closest_mag)

    # Match to data 
    min_groups = dat['ta_snr'][ins][filt][sub][mod][index]

    return min_groups


def perform_calculation(params, n_frame=1, n_skip=0):
    """Calculates all of the outputs and puts them in a dictionary for easy access. 

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
    params : dict, str
        Dictionary of outputs and inputs. If 
        the calculation throws an error it will return a string 
        error message instead.
    """
    
    ## -- TARGET ACQ
    ta_frame_time = set_t_frame(params['infile'], params['ins'], params['subarray_ta'], ta=True)
    
    max_sat_ta_lvl = convert_sat(params['sat_max'], params['sat_mode'],
            params['ins'], params['infile'], ta=True)

    max_group, sat_rate_ta = interpolate_from_dat(params['mag'], params['ins'], params['filt_ta'], 
                                     params['subarray_ta'], params['mod'], params['band'],
                                     ta_frame_time, max_sat_ta_lvl, params['infile'], ta=True)
    min_group = min_groups(params['mag'], params['ins'], params['filt_ta'], params['subarray_ta'], 
                           params['mod'], params['band'], params['infile'])
    min_ta_groups, max_ta_groups = map_to_ta_modes(params['ins'], max_group, min_group)
    t_duration_ta_min = calc_t_duration(min_ta_groups, 1, 1, ta_frame_time)
    t_duration_ta_max = calc_t_duration(max_ta_groups, 1, 1, ta_frame_time)

    ## -- Science obs
    
    # Figure out the rows/cols/amps/px_size and sat
    ins_params = set_params_from_ins(params['ins'], params['subarray'])
    frame_time = set_t_frame(params['infile'], params['ins'], params['subarray'])

    # Run all the calculations
    n_row, n_col, n_amp, px_size, t_frame, n_reset = ins_params
    band_ins = '{0}_{1}'.format(params['band'], params['filt'])
    
    params['sat_max'] = convert_sat(params['sat_max'], params['sat_mode'], params['ins'], params['infile'])

    # Calculate countrate and n_groups if it isn't supplied
    n_group, sat_rate = interpolate_from_dat(params['mag'],
        params['ins'], params['filt'], params['subarray'], params['mod'],
        params['band'], frame_time, params['sat_max'], params['infile'])

    if str(params['n_group']) == 'optimize':
        params['n_group'] = int(n_group)
    else:
        params['n_group'] = int(float(params['n_group']))
    
    ## -- Aditional helpful params

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
    params['min_ta_groups'] = int(min_ta_groups)
    params['max_ta_groups'] = int(max_ta_groups)
    params['t_duration_ta_min'] = t_duration_ta_min
    params['t_duration_ta_max'] = t_duration_ta_max
    params['max_sat_prediction'] = round(sat_rate*frame_time*params['n_group'], 3)
    params['max_sat_ta'] = round(sat_rate_ta*ta_frame_time*max_ta_groups, 3)
    params['min_sat_ta'] = round(sat_rate_ta*ta_frame_time*min_ta_groups, 3)
    
    return params 


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
    

def set_t_frame(infile, ins, sub, ta=False):
    """ Assign the appropriate frame time based on the ins
    and subarray. For now, modes are implied.

    Parameters
    ----------
    infile: str
        The path to the data file.
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
    
    # Read in dict with frame times
    with open(infile) as f:
        frame_time = json.load(f)['frame_time']

    if ta:
        t_frame = frame_time[ins]['ta'][sub]
    else:
        t_frame = frame_time[ins][sub]

    return t_frame


## -- RUN

# There's no function call right now because these inputs are ornery. 
# I can add one if someone feels strongly.

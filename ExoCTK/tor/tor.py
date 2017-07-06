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
import math

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
    return np.round(groups)

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


def interpolate_from_dat(mag, ins, filt, sub, band, t_frame, sat_lvl, infile):
    """
    Interpolates the precalculated pandeia data to estimate the saturation limit.

    Parameters
    ----------
    mag : float
        The magnitude of the source.
    band_ins : str
        The band in which the magnitude is -- only does : 
    t_frame : float
        The seconds per frame for the instrument/subarray.

    Returns
    -------
    n_group : int
        The number of groups that won't oversaturate the detector.
    """
    
    # Create the dictionaries for each filter and select out the prerun data
    dict = create_pandeia_dicts(infile)
    print(ins, filt, sub, band)
    print('DO THESE PARAMS LOOKS RIGHT???')
    
    mag_dat, exptime = dict['input_mag'], dict[ins + '_'+ filt + '_' + sub + '_' + band]

    # Interpolate the given magnitude
    func = interpolate.interp1d(mag_dat, exptime)
    max_sat = func(mag)
    
    # Figure out what it means in wake of the given sat lvl
    max_exptime = sat_lvl/max_sat

    # Calculate the nearest number of groups
    n_group = calc_groups_from_exp_time(max_exptime, t_frame)
    
    return n_group



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


def calc_t_duration(n_group, n_int, n_reset, t_frame, n_frame):
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
    print('Dem ints doe.') 
    n_col, n_amp, n_row = int(n_col), int(n_amp), int(n_row)
    
    if ins == 'NIRSpec':
        n = 2
    if ins in ['NIRCam', 'NIRISS']:
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


def create_tor_dict(transit_time, n_group, mag, band, filt, ins, subarray, sat_mode, sat_max, n_reset, infile, n_frame=1, n_skip=0):
    """Calculates all of the tor things and puts them in a dictionary for easy access. 

    Parameters
    ---------
    transit_time : float
        How long your transit will take (in hours.)
    ins : str
        The instrument at hand.
    subarray : str
        The subarray at hand.
    n_reset : int
        Reset frames per integration.
    n_frame : int
        Frames per group -- always 1 expect maybe brown dwarves.
    n_skip : int
        Skips per integrations -- always 0 except maybe brown dwarves.

    Returns
    -------
    tor_dict : dict
        Dictionary of values related to tor calculations.
    """
    
    # Figure out the rows/cols/amps/px_size and sat
    params = set_params_from_ins(ins, subarray)
    
    # Test for an empty set of parameters
    if params == 'Error!':
        return 'Looks like you mismatched your instrument, filter, and subarray. Please try again.'
    
    # Or continue as normal
    else:
        n_row, n_col, n_amp, px_size, t_frame, MIRI_n_reset = params
        band_ins = str(band) + '_' + str(filt)
        
#        print(px_size)
        sat_max = convert_sat(sat_max, sat_mode, ins)

        # Calculate countrate and n_groups if it isn't supplied
        if n_group == 'optimize':
            n_group = interpolate_from_dat(mag, ins, filt, subarray, band, t_frame, sat_max, infile)
#            countrate = calc_cr(mag, band, temp, px_size, throughput, filt)
#            print('Countrate : ' + str(countrate))
#            n_group = calc_n_group(ins, countrate, sat_max, t_frame, n_frame, n_skip)
        
        n_group = int(float(n_group))
        # Calculate times/ramps/etc
        t_int = calc_t_int(n_group, t_frame, n_frame, n_skip)
        t_ramp = calc_t_ramp(t_int, n_reset, t_frame)
    
        # Calculate nubmer of integrations (THE MEAT)
        if ins == 'MIRI':
            n_reset = MIRI_n_reset
        n_int = calc_n_int(transit_time, n_group, n_reset, t_frame, n_frame)
    
        # Other things that may come in handy who knows?
        t_exp = calc_t_exp(n_int, t_ramp)
        t_duration = calc_t_duration(n_group, n_int, n_reset, t_frame, n_frame)
        obs_eff = calc_obs_efficiency(t_exp, t_duration)

        # Write out dict
        tor_dict = {'n_col': n_col, 'n_row': n_row, 'n_amp': n_amp, 'n_group': n_group, 'n_reset': n_reset, 'sat_max': sat_max,
            'n_frame': n_frame, 'n_skip': n_skip, 'obs_time': transit_time, 't_frame': round(t_frame, 3), 't_int': round(t_int, 3), 't_ramp': t_ramp, 
            'n_int': n_int, 't_exp': round(t_exp/3600, 3), 't_duration': round(t_duration/3600, 3), 'obs_eff': obs_eff}
        return tor_dict


## -- INS CONVERSION THINGS

def convert_sat(sat_max, sat_mode, ins):
     
    if sat_mode == 'well':
        ins_dict = {'NIRSpec': 60000, 'MIRI': 250000, 'NIRCam': 90000, 'NIRISS': 75000}
        sat_max = sat_max*ins_dict[ins]

    print(sat_max)
    return sat_max
    

def set_params_from_ins(ins, subarray):
    """Sets/collects the running parameters from the instrument.
    
    Parameters
    ----------
    ins : str
        Instrument, options are NIRCam, NIRISS, NIRSpec, and MIRI.
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
 
        if ins == 'NIRSpec':
            px_size = (40e-4)**2

            if subarray == 'SUB2048':
                rows, cols = 2048, 32
            if subarray in ['SUB1024A', 'SUB1024B']:
                rows, cols = 1024, 32
            if subarray == 'SUB512':
                rows, cols = 512, 32
            if subarray == 'SUB512S':
                rows, cols = 512, 16
            
            amps = 1 # 4 if not NRSRAPID????
            ft = calc_t_frame(cols, rows, amps, ins)
  
        if ins == 'NIRCam':
            px_size = (18e-4)**2
            
            if subarray == 'FULL':
                rows, cols, amps = 2048, 2048, 4
            if subarray == 'SUB640':
                 rows, cols, amps = 640, 640, 1
            if subarray == 'SUB320':
                rows, cols, amps = 320, 320, 1
            if subarray in ['SUB160', 'SUB160P']:
                rows, cols, amps = 160, 160, 1
            if subarray == 'SUB400':
                rows, cols, amps = 400, 400, 1
            if subarray == 'SUB64P':
                rows, cols, amps = 64, 64, 1
  
            if subarray == 'SUBGRISM256':
                rows, cols, amps = 256, 256, 1
            if subarray == 'SUBGRISM128':
                rows, cols, amps = 128, 2048, 1
            if subarray == 'SUBGRISM64':
                rows, cols, amps = 64, 2048, 1
            
            ft = calc_t_frame(cols, rows, amps, ins)
  
        if ins == 'MIRI':
            px_size = (25e-4)**2
            
            if subarray == 'SLITLESSPRISM':
                rows, cols, ft = 416, 72, .159
            if subarray == 'FULL':
                rows, cols, ft = 1024, 1032, 2.775
            
            amps = 4
            n_reset = 0 
            
        if ins == 'NIRISS':
            px_size = (40e-4)**2
            
            if subarray == 'SUBSTRIP96':
                rows, cols = 2048, 96
            if subarray == 'SUBSTRIP256':
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
    for key in tor_dict:
        print(key, tor_dict[key])

    

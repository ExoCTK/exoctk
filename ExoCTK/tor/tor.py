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
import numpy as np


## -- FUNCTIONS

# - Here I want something that selects all of the right parameters based on ins???

def calc_n_groups(temp, mag, band, ins):
    return temp

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

    return n_ints


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


def calc_t_frame(n_col, n_row, n_amp):
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
    
    t_frame = (n_col/n_amp + 12)*(n_row + 1)*(1e-5)
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

def create_tor_dict(transit_time, n_col, n_row, n_amp, n_group, n_reset, n_frame=1, n_skip=0):
    """Calculates all of the tor things and puts them in a dictionary for easy access. 

    Parameters
    ---------
    transit_time : float
        How long your transit will take (in hours.)
    n_col : int
        Number of columns.
    n_row : int
        Number of rows.
    n_amp : int
        Number of amps.
    n_group : int
        Number of groups per integration.
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
    
    # Calculate frame time and other initial properties
    t_frame = calc_t_frame(n_col, n_row, n_amp)
    t_int = calc_t_int(n_group, t_frame, n_frame, n_skip)
    t_ramp = calc_t_ramp(t_int, n_reset, t_frame)
    
    # Calculate nubmer of integrations (THE MEAT)
    n_int = calc_n_int(transit_time, n_group, n_reset, t_frame, n_frame)
    
    # Other things that may come in handy who knows?
    t_exp = calc_t_exp(n_int, t_ramp)
    t_duration = calc_t_duration(n_group, n_int, n_reset, t_frame, n_frame)
    obs_eff = calc_obs_efficiency(t_exp, t_duration)

    # Write out dict
    tor_dict = {'n_col': n_col, 'n_row': n_row, 'n_amp': n_amp, 'n_group': n_group, 'n_reset': n_reset,
                'n_frame': n_frame, 'n_skip': n_skip, 't_frame': t_frame, 't_int': t_int, 't_ramp': t_ramp, 
                'n_int': n_int, 't_exp': t_exp, 't_duration': t_duration, 'obs_eff': obs_eff}
    return tor_dict

def set_params_from_ins(ins, n_group, n_col, n_row):
    """Sets/collects the running parameters from the instrument.
    
    Parameters
    ----------
    ins : str
        Instrument, options are NIRCam, NIRISS, NIRSpec, MIRI, and SOSS.
    n_group : int
        Groups per integration.
    n_col : int
        Number of detector pixel columns.
    n_row : int
        Number of detector row columns.

    Returns
    -------
    params : dict
        Dictionary of all of the parameters needed to time the TOR.
    """
    
    if ins == 'NIRCam':
        n_reset = 2
    return ins

## -- RUN

if __name__ == "__main__":
    
    transit_time, n_col, n_row, n_amp, n_group, n_reset = 5, 96, 2048, 1, 8, 1
    tor_dict = create_tor_dict(transit_time, n_col, n_row, n_amp, n_group, n_reset)
    for key in tor_dict:
        print(key, tor_dict[key])

    

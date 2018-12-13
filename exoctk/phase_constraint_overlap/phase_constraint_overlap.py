"""Phase contraint overlap tool. This tool calculates the minimum and maximum phase of
the transit based on parameters provided by the user.

Usage:
  calculate_constraint <target_name> [--t_start=<t0>] [--period=<p>] [--obs_duration=<obs_dur>] [--transit_duration=<trans_dur>] [--window_size=<win_size>]
  
Arguments:
  <target_name>                     Name of target
Options:
  -h --help                         Show this screen.
  --version                         Show version.
  --t_start=<t0>                    The starting time of the transit in BJD or HJD.
  --period=<p>                      The period of the transit in days.
  --obs_duration=<obs_dur>          The duration of the observation in hours.
  --transit_duration=<trans_dur>    The duration of the transit in hours.
  --window_size=<win_size>          The window size of the transit in hours.
"""

import math
import os

import argparse
from docopt import docopt
import numpy as np
import requests
import urllib


def build_target_url(target_name):
    '''build restful api url based on target name
    '''
    # Encode the target name string.
    encode_target_name = urllib.parse.quote(target_name, encoding='utf-8')
    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/{}/properties/".format(encode_target_name)
    return target_url

def calculate_phase(period, t0, obs_duration, window_size=1.0):
    print("CALCULATING PHASE CONSTRATIONS WITH VALUES: Period: {}, t0: {}, obs_duration: {}, window_size: {}".format(period, t0, obs_duration, window_size))
    minphase = 1.0 -((obs_duration + window_size)/2.0/24/period)
    maxphase = 1.0 -((obs_duration - window_size)/2.0/24/period)
    
    return minphase, maxphase

def get_canonical_name(target_name):
    """Get ExoMAST prefered name for exoplanet.
    """

    # Generalized target url
    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/identifiers/"
    
    # Create params dict for url parsing. Easier than trying to format yourself.
    params = {"name":target_name}
    
    r = requests.get(target_url, params=params)
    planetnames = r.json()
    canonical_name = planetnames['canonicalName']
    
    return canonical_name

def get_transit_details(target_name):
    '''send request to exomast restful api for target information.
    '''
    canonical_name = get_canonical_name(target_name)

    target_url = build_target_url(canonical_name)
    
    r = requests.get(target_url)
    
    if r.status_code == 200:
        target_data = r.json()
    else:
        print('Whoops, no data for this target!')

    # Some exoplanets have multiple catalog entries
    # Temporary... Need to write a parser to select catalog
    target_data = target_data[0]
    
    period = target_data['orbital_period']
    t0 = target_data['transit_time']
    transit_duration = target_data['transit_duration']
    obs_duration = np.min((6, 3*transit_duration + 1))

    return period, t0, obs_duration


def phase_overlap_constraint(target_name, period=None, t0=None, obs_duration=None, window_size=None):
    ''' The main function to calculate the phase overlap constraints.
        We will update to allow a user to just plug in the target_name 
        and get the other variables.
        
        Parameters
        ----------
        period : float
            The period of the transit in days. 
        t0 : float
            The start time in BJD or HJD.
        obs_duration : float
            The duration of the observation in hours.
        winSize : float
            The window size of transit in hours. Default is 1 hour.
        target_name : string
            The name of the target transit. 
        
        Returns
        -------
        minphase : float
            The minimum phase constraint.
        maxphase : float
            The maximum phase constraint.'''

    # If user only gave target_name.
    if all(arg is None for arg in [period, t0, obs_duration, window_size]):
        period, t0, obs_duration = get_transit_details(target_name)
        minphase, maxphase = calculate_phase(period, t0, obs_duration)
    # If user doesn't provide observation details but wants to changes window size
    elif all(arg is None for arg in [period, t0, obs_duration]):
        period, t0, obs_duration = get_transit_details(target_name)
        minphase, maxphase = calculate_phase(period, t0, obs_duration, window_size)
    # Make sure users are passing all transit information if providing their own.
    # Can't supply period without t0 or obs_duration or any combination of these variables!
    elif any(arg is None for arg in [period, t0, obs_duration]):
        raise ValueError('If passing transit information, must include all values! period: {}, t0: {}, obs_duration: {}'.format(period, t0, obs_duration))
    # Else user needs to provide all of the data, calculate using their values.
    else:
        minphase, maxphase = calculate_phase(period, t0, obs_duration, window_size)
            
    print("MINIMUM PHASE: {}, MAXIMUM PHASE: {}".format(minphase, maxphase))


def parse_args():
    """Parses command line arguments.
    Returns
    -------
    args : obj
        An ``argparse`` object containing all of the added arguments.
    """

    # Create help string
    period_help = 'The period of the transit in days.'
    t0_help = 'The starting time of the transit in BJD or HJD.'
    obsDur_help = 'The duration of the observation in hours.'
    transitDur_help = 'The duration of the transit in hours.'
    winSize_help = 'The window size of the transit in hours.'
    target_name_help = 'The name of the transiting target.'

    # Add time arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-n --target_name',
        dest='target_name',
        action='store',
        type=str,
        required=True,
        help=target_name_help)
    
    parser.add_argument('-t --t_start',
        dest='t0',
        action='store',
        type=float,
        help=t0_help)

    parser.add_argument('-p --period',
        dest='period',
        action='store',
        type=float,
        help=period_help)
    
    parser.add_argument('-odur --obsDur',
        dest='obsDur',
        action='store',
        type=float,
        help=obsDur_help)
    
    parser.add_argument('-tdur --transitDur',
        dest='transitDur',
        action='store',
        type=float,
        help=transitDur_help)
        
    parser.add_argument('-w --winSize',
        dest='winSize',
        action='store',
        type=float,
        required=False,
        help=winSize_help)

    
    # Parse args
    args = parser.parse_args()

    return args
    

if __name__ == '__main__':
    # args = parse_args()
    args = docopt(__doc__, version='0.1')

    phase_overlap_constraint(args['<target_name>'], args['--period'], 
                             args['--t_start'], args['--transit_duration'], 
                             args['--window_size'])
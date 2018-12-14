#! /usr/bin/env python
"""Phase contraint overlap tool. This tool calculates the minimum and maximum phase of
the transit based on parameters provided by the user.

Authors:
    Catherine Martlin, 2018
    Mees Fix, 2018

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
  --window_size=<win_size>          The window size of the transit in hours [default: 1.0]
"""

import math
import os

import argparse
from docopt import docopt
import numpy as np
import requests
import urllib


def build_target_url(target_name):
    '''Build restful api url based on target name.

    Parameters
        ----------
        target_name : string
            The name of the target transit. 

        Returns
        -------
        target_url : string
    '''
    # Encode the target name string.
    encode_target_name = urllib.parse.quote(target_name, encoding='utf-8')
    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/{}/properties/".format(encode_target_name)

    return target_url

def calculate_phase(period, obsDur, winSize):
    ''' Function to calculate the min and max phase. 

        Parameters
        ----------
        period : float
            The period of the transit in days. 
        obsdur : float
            The duration of the observation in hours.
        winSize : float
            The window size of transit in hours. Default is 1 hour.

        Returns
        -------
        minphase : float
            The minimum phase constraint.
        maxphase : float
            The maximum phase constraint. '''

    minphase = 1.0 - ((obsDur + winSize)/2.0/24/period)
    maxphase = 1.0 - ((obsDur - winSize)/2.0/24/period)
    
    return minphase, maxphase

def calculate_obsDur(transitDur):
    ''' Function to calculate the min and max phase. 

        Parameters
        ----------
        transitDur : float
            The duration of the transit in hours.

        Returns
        -------
        obsdur : float
            The duration of the observation in hours. '''

    obsDur = np.min((6, 3*transitDur+1))

    return obsDur

def get_canonical_name(target_name):
    '''Get ExoMAST prefered name for exoplanet.

        Parameters
        ----------
        target_name : string
            The name of the target transit. 

        Returns
        -------
        canonical_name : string
    '''

    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/identifiers/"
    
    # Create params dict for url parsing. Easier than trying to format yourself.
    params = {"name":target_name}
    
    r = requests.get(target_url, params=params)
    planetnames = r.json()
    canonical_name = planetnames['canonicalName']
    
    return canonical_name

def get_transit_details(target_name):
    '''Send request to exomast restful api for target information.
        
        Parameters
        ----------
        target_name : string
            The name of the target transit. 

        Returns
        -------
        period : float
            The period of the transit in days. 
        transitDur : float
            The duration of the transit in hours. 

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
            The maximum phase constraint. '''

    if obs_duration == None:
        if period == None:
            period, transit_dur, t0 = get_transit_details(target_name)
        obs_duration = calculate_obsDur(transit_dur)

    minphase, maxphase = calculate_phase(period, obs_duration, window_size)
    
    # Is this the return that we want? Do we need to use t0 for something? 
    print('MINIMUM PHASE: {}, MAXIMUM PHASE: {}'.format(minphase, maxphase))

if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')

    # Ugh, docopt datatypes are funky.
    # This converts entries from strs to floats
    for k,v in args.items():
        try:
            args[k] = float(v)
        except (ValueError, TypeError):
            # Handles None and char strings.
            continue
    
    phase_overlap_constraint(args['<target_name>'], args['--period'], 
                             args['--t_start'], args['--transit_duration'], 
                             args['--window_size'])

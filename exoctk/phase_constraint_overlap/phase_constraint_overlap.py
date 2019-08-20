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

from exoctk.utils import get_target_data

def build_target_url(target_name):
    '''build restful api url based on target name
    '''
    # Building the URL will be tricky, we should use
    # as much source code as we can from the archives group..
    # Name matching, filtering, etc.
    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/Wasp-6%20b/properties/"
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

    minphase = 1.0 - ((0.5 + obsDur + winSize)/2.0/24/period)
    maxphase = 1.0 - ((0.5 + obsDur - winSize)/2.0/24/period)
    
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
            The duration of the observation in hours. Maximum of 6 hours. '''

    obsDur = np.min((6, 3*transitDur+1))

    return obsDur


def get_transit_details(target_name):
    '''send request to exomast restful api for target information.
    '''

    target_url = build_target_url(target_name)
    r = requests.get(target_url)

    if r.status_code == 200:
        target_data = r.json()
    else:
        print('Whoops, no data for this target!')

    # Requests from restapi can have multiple catalogs.
    # Discussion about which one we should use should 
    # be had in the future.
    for item in target_data:
        print(item['catalog_name'], item['orbital_period'])

    # are t0, obsDur, and winSize available via ExoMAST api?
    return period, t0, obsDur, winSize 


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
            data = get_target_data(target_name)
            
            period = data['orbital_period']
            transit_dur = data['transit_duration']
            t0 = data['transit_time']

        obs_duration = calculate_obsDur(transit_dur)

    minphase, maxphase = calculate_phase(period, obs_duration, window_size)
    
    # Is this the return that we want? Do we need to use t0 for something? 
    print('MINIMUM PHASE: {}, MAXIMUM PHASE: {}'.format(minphase, maxphase))

# Need to make entry point for this!
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

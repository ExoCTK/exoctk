#! /usr/bin/env python
"""Phase contraint overlap tool. This tool calculates the minimum and maximum phase of
the transit based on parameters provided by the user.

Authors:
    Catherine Martlin, 2018
    Mees Fix, 2018

Usage:
  calculate_constraint <target_name> [--period=<p>] [--obs_duration=<obs_dur>] [--transit_duration=<trans_dur>] [--window_size=<win_size>]
  
Arguments:
  <target_name>                     Name of target
Options:
  -h --help                         Show this screen.
  --version                         Show version.
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

    # are obsDur, and winSize available via ExoMAST api?
    return period, obsDur, winSize 


def phase_overlap_constraint_main(target_name, period=None, obs_duration=None, window_size=None):
    ''' The main function to calculate the phase overlap constraints.
        We will update to allow a user to just plug in the target_name 
        and get the other variables.
        
        Parameters
        ----------
        period : float
            The period of the transit in days. 
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
    
    data = get_target_data(target_name)
    transitDur = data['transit_duration'] 

    if obs_duration == None:
        obs_duration = calculate_obsDur(transitDur)
    else:
        obs_dur_min = transitDur * 2
        if obs_duration < obs_dur_min:
            print('WARNING: Your observation duration is less than twice the transit duration. You need to increase it.')
    
    if period == None:  
        period = data['orbital_period']

    if window_size == None: 
        window_size = 1.0

    minphase, maxphase = calculate_phase(period, obs_duration, window_size)
    
    # Is this the return that we want?
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
                             args['--transit_duration'], args['--window_size'])

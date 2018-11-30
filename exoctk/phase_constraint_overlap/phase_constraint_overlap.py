import math
import os

import argparse
import numpy as np
import requests


def build_target_url(target_name):
    '''build restful api url based on target name
    '''
    # Building the URL will be tricky, we should use
    # as much source code as we can from the archives group..
    # Name matching, filtering, etc.
    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/Wasp-6%20b/properties/"
    return target_url

def calculate_phase(period, t0, obsDur, winSize):
    minphase = 1.0 -((obsDur + winSize)/2.0/24/period)
    maxphase = 1.0 -((obsDur - winSize)/2.0/24/period)
    
    return minphase, maxphase

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


def phase_overlap_constraint(period, t0, obsDur, winSize, target_name):
    ''' The main function to calculate the phase overlap constraints.
        We will update to allow a user to just plug in the target_name 
        and get the other variables.
        
        Parameters
        ----------
        period : float
            The period of the transit in days. 
        t0 : float
            The start time in BJD or HJD.
        obsdur : float
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

    if target_name == parser.get_default('target_name'):
        print("Using the default values:")
        print(target_name)
        print(parser.get_default('target_name'))
        minphase, maxphase = calculate_phase(period, t0, obsDur, winSize)

    elif target_name != parser.get_default('target_name') and period != parser.get_default('period') or t0 != parser.get_default('t0') or obsDur != parser.get_default('obsDur') or winSize != parser.get_default('winSize'):
        print(target_name)
        print(parser.get_default('target_name'))
        print("Using your input values for calculations:")
        minphase, maxphase = calculate_phase(period, t0, obsDur, winSize)

    else:
        print("Have a different Target Name, but same defaults. Retriving variables from exoplanets.org for that target name:")
        #period,t0,obsDur,winSize = get_transit_details(target_name)
        minphase, maxphase = calculate_phase(period, t0, obsDur, winSize)
    
    print(minphase,maxphase)

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
    parser.add_argument('-p --period',
        dest='period',
        action='store',
        type=float,
        help=period_help,
        default=3.7354845)
    
    parser.add_argument('-t --t0',
        dest='t0',
        action='store',
        type=float,
        help=t0_help,
        default= 2454592.80154)

    parser.add_argument('-odur --obsDur',
        dest='obsDur',
        action='store',
        type=float,
        help=obsDur_help,
        default= 34447.38/3600.)
    
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
        help=winSize_help,
        default=1.0)

    parser.add_argument('-n --target_name',
        dest='target_name',
        action='store',
        type=str,
        help=target_name_help,
        default='WASP-17b')

    # Parse args
    args = parser.parse_args()

    return args
    

if __name__ == '__main__':
    
    args = parse_args()

    phase_overlap_constraint(args.period, args.t0, args.obsDur, args.winSize, args.target_name)
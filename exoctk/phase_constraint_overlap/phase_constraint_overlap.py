#! /usr/bin/env python

""" Provide phase overlap constraints for transits. 

This script will provide the min and max phases for a transit based
on the target name provided. 

Authors
-------
    Catherine Martlin, 2018
    Mees Fix, 2018

Use
---
    This script is can be executed as such:
    ::
        python phase_constraint_overlap.py 

    However, it will mostly be used through the webapp interface. 
"""

import math
import os

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

    for item in target_data:
         print(item['catalog_name'], item['orbital_period'])

    # NExSci has more up-to-date options for exoplanet data so we will use that
    # catalog by default.
    # cats = [d['catalog_name'] for d in target_data if 'nexsci' in d]
    # print(cats)

    # are t0, obsDur, and winSize available via ExoMAST api?
    # return period, transitDur, t0


def phase_constraint_overlap(period, t0, obsDur, winSize, target_name):
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

    if winSize == None:
        winSize = 1.0  # Default value.

    if obsDur == None:
        if period == None:
            period, transitDur, t0 = get_transit_details(target_name)
        obsDur = calculate_obsDur(transitDur)

    minphase, maxphase = calculate_phase(period, t0, obsDur, winSize)
    
    print(minphase, maxphase) # Is this the return that we want? Do we need to use t0 for something? 

if __name__ == '__main__':

    target_name = None
    period = None
    t0 = None 
    obsDur= None 
    winSize= None 

    phase_constraint_overlap(target_name, period, t0, obsDur, winSize)
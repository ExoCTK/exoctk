"""
This is a module for the NCS -- NIRCam Coronagraphy
Simulation package of ExoCTK. It's based on a script
written by Laurent Pueyo and Jason Wang -- and highly
dependent upon WebbPSF and the pandeia_coronagraphy packages.

This is a wrapper of the pandeia_coronagraphy package meant
to simulate NIRCam iterating through multiple roll angles. 
Use it to see if your transit will be visible! (Or don't. 
I can't tell you what to do.)

Authors
-------
Jules Fowler, Februrary 2018
Laurent Pueyo
Jason Wang

Use 
---
IDK what parameters will make it to the final, sooo...


Notes
-----
This is sloooow AF. Someday, we will do the interpolate?
"""

## -- IMPORTS
from copy import deepcopy
import json
import os
import sys

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

from pandeia_coronagraphy import analysis
from pandeia_coronagraphy import engine
from pandeia_coronagraphy import scene
from pandeia_coronagraphy import transformations

## -- FUNCTIONS

def main(filters, masks, mode='run'):
    """ Overall function to run through one roll angle?
    set of a given scene or test that all modes work."""
    
    if mode == 'run':
        scene = intialize_scene()
        for planet in planets:
            add_planet(planet)
        run_calc(sene)

    elif mode == 'test':
        
    else:
        print('I do not understand what you want from me.')


        
        
        


def build_scene

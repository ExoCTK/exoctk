"""UNDER CONSTRUCTION"""

import glob
import os
import pysiaf

import astropy.coordinates as crd
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astroquery.irsa import Irsa
from matplotlib import cm
from scipy.io import readsav
from astropy.io import fits
from exoctk.utils import get_env_variables

EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
if not EXOCTK_DATA:
    print('WARNING: The $EXOCTK_DATA environment variable is not set. '
          'Contamination overlap will not work. Please set the '
          'value of this variable to point to the location of the exoctk_data '
          'download folder.  Users may retreive this folder by clicking the '
          '"ExoCTK Data Download" button on the ExoCTK website, or by using '
          'the exoctk.utils.download_exoctk_data() function.'
          )
    TRACES_PATH = None
else:
    TRACES_PATH = os.path.join(EXOCTK_DATA,  'exoctk_contam', 'traces')


class fieldSim:
    def __init__(self, RA, DEC, INSTRUMENT):
        self.ra = RA
        self.dec = DEC
        self.instrument = INSTRUMENT

        if 'NIRCam' in self.instrument:
            siaf = pysiaf.Siaf('NIRCam')
            dimX, dimY = 51, 1343
            rad = 2.5
            pixel_scale = siaf.
            xSweet, ySweet = siaf. , siaf.
        elif 'MIRI' in self.instrument:
            siaf = pysiaf.Siaf('MIRI')
            dimX, dimY = 51, 1343
            rad = 2.5
            pixel_scale = siaf.
            xSweet, ySweet = siaf. , siaf.

    def getStars(self, ):

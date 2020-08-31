""" A couple functions we will use to generate example plots
in the Jupyter notebooks. Makes them look a little cleaner.
"""
import os
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from scipy.io import readsav
from matplotlib import cm

EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
if not EXOCTK_DATA:
    print('WARNING: The $EXOCTK_DATA environment variable is not set. Contamination overlap will not work. Please set the '
          'value of this variable to point to the location of the exoctk_data '
            'download folder.  Users may retreive this folder by clicking the '
            '"ExoCTK Data Download" button on the ExoCTK website, or by using '
            'the exoctk.utils.download_exoctk_data() function.'
          )
    TRACES_PATH = None

TRACES_PATH = os.path.join(EXOCTK_DATA,  'exoctk_contam', 'traces')

#def distorCorr():

def plotTemps(TEMPS, allRA, allDEC):
    # Getting the color palette
    colors = cm.get_cmap('viridis', len(TEMPS))
    colors_0 = np.asarray(colors.colors)

    # Assigning index arrays to TEMPS array
    i = TEMPS.argsort()
    ii = TEMPS.argsort().argsort()

    # Matching the colors to the corresponding magnitude
    colors = colors_0
    starsx, starsy = allRA[i], allDEC[i]

    for x, y, c in zip(starsx, starsy, colors):
        plt.scatter(x, y, marker='*', s=100, color=c, picker=True, lw=0.5, edgecolor='white')

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis,
                               norm=plt.Normalize(vmin=TEMPS.min(),
                                                  vmax=TEMPS.max()))
    sm._A = []
    plt.colorbar(sm)

def traceLength(inst):
    # Getting example trace to calculate rough estimate of trace lengths
    if 'NIRCam' in inst:
        FILE = 'rot_o1_6000.0.fits'
    elif 'MIRI' in inst:
        FILE = 'LOWbg_6000.0.fits'
    elif 'NIRISS' in inst:
        FILE = 'modelOrder12_teff6000.sav'

    trFile = os.path.join(TRACES_PATH, inst.replace(' ', '_'), FILE)
    trData = readsav(trFile)['modelo12'] if 'NIRISS' in inst \
                                         else fits.getdata(trFile, 1)
    trData = trData[0]
    print(np.shape(trData))
    ax = 1 if 'NIRCam' in inst else 0
    peak = trData.max()

    # the length of the trace
    targ_trace_start = np.where(trData > 0.0001*peak)[ax].min()
    targ_trace_stop = np.where(trData > 0.0001*peak)[ax].max()

    return targ_trace_start, targ_trace_stop

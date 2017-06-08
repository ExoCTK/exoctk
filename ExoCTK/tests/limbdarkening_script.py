#! /usr/bin/env python3

import os, sys
sys.path.append('/Users/gbruno/python/source/ExoCTK/')
from ExoCTK import core
from ExoCTK.ldc import ldcfit as lf
from ExoCTK.ldc import ldcplot as lp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import pdb

# Path for the models
#fits_files = '/Users/gbruno/idl/projects/limb_test/phxinten/'
#fits_files = '/Users/gbruno/python/projects/ExoCTK/kurucz/'
fits_files = '/user/jfilippazzo/Models/ATLAS9/'

def ld_bins():

    # Define wavelenght bins
    bins = np.arange(1.125, 1.650, 0.00453128*9)
    # Stellar parameters
    teff, logg, FeH = 5250, 4.0, 0.0
    #teff, logg, FeH = 2500, 4.5, 0.0
    # Initialize grid
    model_grid = core.ModelGrid(fits_files)

    # Filter?
    K_band = core.Filter('Kepler.K')
    #plt.plot(*K_band.rsr)

    # Open text file to save results
    fname = 'test1_atlas_ldc_ongrid_Kband_quadratic.dat'
    ldfile = open(fname, 'w')
    ldfile.write('bin_i\tbin_f\tu1\tu2\n')
    ldfile.write('-----\t-----\t--\t--\n')
    for i in range(1, len(bins)):
        # Compute LD coefficients
        model_grid.customize(Teff_rng=(5000,5400), logg_rng=(3.5,4.5), FeH_rng=(-.5,.5), wave_rng=(bins[i - 1], bins[i]))
        interpolation = lf.ldc(teff, logg, FeH, model_grid, 'quadratic', mu_min = 0.05, plot=False, bandpass = True)
        plt.close('all')
        coeffs = interpolation['quadratic']['coeffs']
        #coeffs, mu, radius = lf.ldc(teff, logg, FeH, model_grid, '4-parameter', plot=False)    # Old version

        # Write to file
        #ldfile.write('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' %(bins[i - 1], bins[i], coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
        ldfile.write('%.3f\t%.3f\t%.3f\t%.3f\n' %(bins[i - 1], bins[i], coeffs[0], coeffs[1]))
    # Whole range
    model_grid.customize(Teff_rng=(5000,5400), logg_rng=(3.5, 4.5), FeH_rng=(-0.5,0.5), wave_rng=(0.3, 1.0))
    #coeffs, mu, radius = lf.ldc(teff, logg, FeH, model_grid, '4-parameter', plot=False, bandpass = K_band)
    style = ['--', '-']
    for ii, laws in enumerate(['quadratic']):
        #coeffs, mu, radius = lf.ldc(teff, logg, FeH, model_grid, '4-parameter', plot=True)
        interpolation = lf.ldc(teff, logg, FeH, model_grid, laws, ls = style[ii], plot=True, mu_min = 0.05, bandpass = True)
        coeffs = interpolation[laws]['coeffs']
        #print('0.3', '1.0', coeffs, interpolation[laws]['err'])#, coeffs[1], coeffs[2], coeffs[3])
        #ldfile.write('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' %(float('0.3'), float('1.0'), coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
        ldfile.write('%.3f\t%.3f\t%.3f\t%.3f\n' %(float('0.3'), float('1.0'), coeffs[0], coeffs[1]))
    ldfile.close()
    mu = interpolation['mu']
    #mu = np.linspace(0.0, 1.0, 17)
    # From Sing's fits...
    #c1, c2, c3, c4 = 0.71364385, -0.68152905, 1.3952835, -0.62918711
    #c1, c2, c3, c4 = 0.71637808, -0.70077065, 1.4193015, -0.63346020
    a, b = 0.50412642, 0.19013756
    #a = 0.65774273

    # Comparison plot
    #sing4 = 1. - c1*(1. - mu**0.5) - c2*(1. - mu) - c3*(1. - mu**1.5) - c4*(1. - mu**2)
    #ldc4 = 1. - coeffs[0]*(1. - mu**0.5) - coeffs[1]*(1. - mu) - coeffs[2]*(1. - mu**1.5) - coeffs[3]*(1. - mu**2)
    #sing1 = 1 - a*(1. - mu)
    sing2 = 1. - a*(1. - mu) - b*(1. - mu)**2
    #ldc2 = 1. - coeffs[0]*(1. - mu) - coeffs[1]*(1. - mu)**2
    #plt.plot(mu, sing1, 'g-', label = 'Sing linear')
    plt.plot(mu, sing2, 'g-', label = 'Sing quadratic')
    #plt.plot(mu, sing4, 'r-', label = 'Sing 4 coeffs')
    #plt.plot(mu, ldc4, 'g-', label = 'LDC from coeffs')
    plt.xlabel('$\mu$', fontsize = 18)
    plt.ylabel('$I(\mu)/I(\mu = 1)$', fontsize = 18)
    plt.title('Teff = ' + str(teff) + ', log g = ' + str(logg) + ', [Fe/H] = ' + str(FeH))
    plt.legend(loc = 'best')
    plt.savefig('LDCvsSing_phoenix_ongrid_quadratic.png')

    return

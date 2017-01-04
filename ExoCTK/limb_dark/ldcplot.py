#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module of plotting tools for the limb darkening subpackage.
"""
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import os
from astropy.io import fits
try:
    from ExoCTK import core
    from ExoCTK.limb_dark import ldcfit
except ImportError:
    from ExoCTK.ExoCTK import core
    from ExoCTK.ExoCTK.limb_dark import ldcfit

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def test_plot():
    """
    Reproduce Jeff's plot from the TRAPPIST-1 paper
    """
    spec_files = '/Users/jfilippazzo/Documents/Modules/limb_dark_jeff/limb/specint/'
    G = core.ModelGrid(spec_files, 'foobar')
    G.customize(Teff_rng=(2500,2600), logg_rng=(5,5.5), FeH_rng=(-0.5,0.5),  wave_rng=(1.1,1.7))
    
    _ = ldcfit.ldc(2600,5.5,0,G,'nonlinear', plot=plt.gcf(), **{'ls':'-', 'c':'b'})
    _ = ldcfit.ldc(2600,5.5,0,G,'nonlinear', plot=plt.gcf(), **{'ls':'-', 'c':'b'})
    _ = ldcfit.ldc(2600,5.,0,G,'nonlinear', plot=plt.gcf(), **{'ls':'-', 'c':'g'})
    _ = ldcfit.ldc(2600,5.,0,G,'quadratic', plot=plt.gcf(), **{'ls':'--', 'c':'g'})
    _ = ldcfit.ldc(2500,5.,0,G,'quadratic', plot=plt.gcf(), **{'ls':'--', 'c':'m'})
    _ = ldcfit.ldc(2500,5.,0,G,'nonlinear', plot=plt.gcf(), **{'ls':'-', 'c':'m'})
    _ = ldcfit.ldc(2500,5.5,0,G,'nonlinear', plot=plt.gcf(), **{'ls':'-', 'c':'r'})
    _ = ldcfit.ldc(2500,5.5,0,G,'quadratic', plot=plt.gcf(), **{'ls':'--', 'c':'r'})
    _ = ldcfit.ldc(2550,5.22,0,G,'quadratic', plot=plt.gcf(), **{'ls':'--', 'c':'k', 'lw':3})
    _ = ldcfit.ldc(2550,5.22,0,G,'nonlinear', plot=plt.gcf(), **{'ls':'-', 'c':'k', 'lw':3})

def ld_v_mu(data, params):
    """
    Plot the limb darkening versus mu
    
    Parameters
    ----------
    data: array-like, str
        The list of (coefficient, mu) values or the path
        to a FITS file with COEFFS and MU extensions
    
    """    
    try:
    
        # Get the data directly...
        if isinstance(data, (list,tuple,np.ndarray)):
            coeffs, mu = data[:2]

        # ...or from a FITS file
        else:
            coeffs = fits.getdata(path, extname='COEFFS')
            mu = fits.getdata(path, extname='MU')
        
        # Create the plot
        plt.figure()
        plt.xlabel(r'$\mu$m')
        plt.ylabel=('Limb Darkening')
        
        # Get the param ranges
        
        for c,m in zip(coeffs, mu):
            plt.scatter(coeffs, mu)
    
    except:
        print('Sorry, that input is not plottable!')
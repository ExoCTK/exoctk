from __future__ import absolute_import, unicode_literals, print_function
import pymultinest
import math, os
from fm import *
import pdb
import numpy as np
import pickle
from matplotlib.pyplot import *
#output path--note this will be name of folder
#for multinest output and also the generated
#pickle 
outpath="./OUTPUT/hp41_wfc3_cc"
if not os.path.exists(outpath): os.mkdir(outpath)  #creating folder to dump MultiNest output


#load crosssections between wnomin and wnomax
xsects=xsects_HST(1500,27900) # 1500 cm-1 (6.66 um) to 28000 cm-1 (0.358 um)

# log-likelihood
def loglike(cube, ndim, nparams):

    #setting default parameters---will be fixed to these values unless replaced with 'theta'
    #planet/star system params--typically not free parameters in retrieval
    Rp = 1.685 #0.930#*x[4]# Planet radius in Jupiter Radii--this will be forced to be 10 bar radius--arbitrary (scaling to this is free par)
    Rstar = 1.786 #0.598   #Stellar Radius in Solar Radii
    M = 0.795 #1.78    #Mass in Jupiter Masses

    #TP profile params (3--Guillot 2010, Parmentier & Guillot 2013--see Line et al. 2013a for implementation)
    #Tirr= 2091 #1500#x[0]#544.54 #terminator **isothermal** temperature--if full redistribution this is equilibrium temp
    #logKir= -1.5  #TP profile IR opacity controlls the "vertical" location of the gradient
    #logg1= -0.7     #single channel Vis/IR opacity. Controls the delta T between deep T and TOA T
    Tint= 200.
    
    #Composition parameters---assumes "chemically consistnat model" described in Kreidberg et al. 2015
    #logMet= 0.0 #x[1]#1.5742E-2 #.   #Metallicity relative to solar log--solar is 0, 10x=1, 0.1x = -1 used -1.01*log10(M)+0.6
    #logCtoO= -0.26#x[2]#-1.97  #log C-to-O ratio: log solar is -0.26
    logPQCarbon = -5.5  #CH4, CO, H2O Qunech pressure--forces CH4, CO, and H2O to constant value at quench pressure value
    logPQNitrogen = -5.5  #N2, NH3 Quench pressure--forces N2 and NH3 to ""  --ad hoc for chemical kinetics--reasonable assumption
    
    #A&M Cloud parameters--includes full multiple scattering (for realzz) in both reflected and emitted light
    #logKzz=7 #log Rayleigh Haze Amplitude (relative to H2)
    #fsed=3.0 #haze slope--4 is Rayeigh, 0 is "gray" or flat.  
    #logPbase=-1.0  #gray "large particle" cloud opacity (-35 - -25)
    #logCldVMR=-25.0 #cloud fraction
    
    #simple 'grey+rayleigh' parameters--non scattering--just pure extinction
    logKcld = -40
    logRayAmp = -30
    RaySlope = 0

    #unpacking parameters to retrieve
    Tirr, logMet, logCtoO, logKzz, fsed ,logPbase,logCldVMR,xRp, logKir, logg1 = cube[0],cube[1],cube[2],cube[3],cube[4],cube[5],\
                                                                          cube[6],cube[7], cube[8], cube[9]

    ##all values required by forward model go here--even if they are fixed
    x=np.array([Tirr, logKir,logg1, Tint,logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp*xRp, Rstar, M, logKzz, fsed,logPbase,logCldVMR, logKcld, logRayAmp, RaySlope])
    gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1., 1., 1.]) #can be made free params if desired (won't affect mmw)#can be made free
    foo=fx_trans(x,wlgrid,gas_scale,xsects)
    y_binned=foo[0]

    loglikelihood=-0.5*np.sum((y_meas-y_binned)**2/err**2)  #nothing fancy here
    return loglikelihood


# prior transform
Tmin, Tmax = 1000., 2500.      # Tirr: 1000 - 25000
MHmin, MHmax = -2., 3.        # [M/H]: -2.0 - 3.0 (0.01x - 1000x) 
logCOmin, logCOmax = -2., 0.3  # log(C/O): -2 - 0.3 (0.01 to 2.0 )
logKZZmin, logKZZmax = 5., 11. # log(Kzz): 5 - 11 (1E5 - 1E11 cm2/s)
FSEDmin, FSEDmax = 0.5, 6.0
logPBASEmin, logPBASEmax = -6., 1.5 # logPbase: -6.0 - 1.5 (1 ubar - 30 bar)
logCldVMRmin, logCldVMRmax = -15., -2. # logCldVMR: -15 - -2
xRPmin, xRPmax = 0.5, 1.5 # xRp: 0.5 - 1.5 (multiplicative factor to "fiducial" 10 bar radius)
logKIRmin, logKIRmax = -3, 0.
logg1min, logg1max = -3., 0.


def prior(cube,ndim,nparams):
    #prior ranges...
    cube[0] = (Tmax - Tmin) * cube[0] + Tmin                             # Tirr: 400 - 1800
    cube[1] = (MHmax - MHmin)*cube[1] + MHmin                            # [M/H]: -2.0 - 3.0 (0.01x - 1000x)
    cube[2] = (logCOmax - logCOmin)*cube[2] + logCOmin                   # log(C/O): -2 - 0.3 (0.01 to 2.0 )
    cube[3] = (logKZZmax - logKZZmin)*cube[3] + logKZZmin                 # log(Kzz): 5 - 11 (1E5 - 1E11 cm2/s)
    cube[4] = (FSEDmax - FSEDmin)*cube[4] + FSEDmin                 # fsed: 0.5 - 6.0
    cube[5] = (logPBASEmax - logPBASEmin)*cube[5] + logPBASEmin                 # logPbase: -6.0 - 1.5 (1 ubar - 30 bar)
    cube[6] = (logCldVMRmax - logCldVMRmin)*cube[6] + logCldVMRmin                   # logCldVMR: -15 - -2
    cube[7] = (xRPmax - xRPmin)*cube[7] + xRPmin                   # xRp: 0.5 - 1.5 (multiplicative factor to "fiducial" 10 bar radius)
    cube[8] = (logKIRmax - logKIRmin)*cube[8] + logKIRmin
    cube[9] = (logg1max - logg1min)*cube[9] + logg1min

#####loading in data##########
wmin, wmax, d, derr = np.loadtxt('hat-p-41-data.dat',unpack=True)
wlgrid = (wmin + wmax)*0.5
y_meas = d*1e-6
err = derr*1e-6
idx = np.argsort(wlgrid)
wlgrid, y_meas, err = wlgrid[idx], y_meas[idx], err[idx]
#wlgrid, y_meas, err=np.loadtxt('w43b_trans.txt').T
outname=outpath+'.pic'  #dynesty output file name (saved as a pickle)
Nparam=10  #number of parameters--make sure it is the same as what is in prior and loglike
Nlive=500 #number of nested sampling live points

#calling pymultinest
pymultinest.run(loglike, prior, Nparam, outputfiles_basename=outpath+'/template_',resume=False, verbose=True,n_live_points=Nlive, importance_nested_sampling=False)

#converting pymultinest output into pickle format (weighting sampling points by volumes to get "true" posterior)
a = pymultinest.Analyzer(n_params = Nparam, outputfiles_basename=outpath+'/template_')
s = a.get_stats()
output=a.get_equal_weighted_posterior()
pickle.dump(output,open(outname,"wb"))

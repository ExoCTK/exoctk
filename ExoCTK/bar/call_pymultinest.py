from __future__ import absolute_import, unicode_literals, print_function
import pymultinest
import math, os
from fm import *
import pdb
import numpy as np
import pickle
if not os.path.exists("chains"): os.mkdir("chains")

#prior--only use this to transform "unit (0-1)" coordinates to
#physical values...
#for instance to transform something to have a value between -6 and 6
#would require -6+12*cube[0]. When cube[0]=0, then the transformed 
#cube[0] will equal -6+12.*0, or -6. When cube[0]=1, then the transformed cube[0]
#will be -6+12*1, or 6.
def prior(cube, ndim, nparams):
	#TP profile params
        cube[0]=300+2700*cube[0] #Irradiation Temp (linear)   --  300-3000 K
        cube[1]=-3+3*cube[1] #TP profile grey IR opacity (log) -- 1E-3 - 0 cm2/g (Freedman+2014)
        cube[2]=-3+4*cube[2] #vis/IR opacity (log) -- 1E-3 to 10
	#Composition parameters
        cube[3]=-2+5*cube[3] #metallicity x solar (log), 0=solar, 1=10x, -1=0.1x etc.--0.01 - 1000 x solar
        cube[4]=-2+4*cube[4] #C/O (log)-- -0.26 = solar (C/O=0.55)  0.01 - 100  
	cube[5]=1.5-8*cube[5] #carbon qunech pressure (log) 30bars - 0.3 ubar 
	cube[6]=1.5-8*cube[6] #nitrogen qunech pressure (log) 30bars - 0.3 ubar
	#Generic Mie Cloud params
	cube[7]=-2+6*cube[7] #Ray Haze Amplitude
	cube[8]=4*cube[8] #ray haze amp
	cube[9]=-7+9.5*cube[9]  #logPc
	#10 bar radius scaling
	cube[10]=0.5+1.*cube[10] #xRp

#loglikelihood
def loglike(cube, ndim, nparams):
    Tirr, logKir,logg1, logMet, logCtoO, logPQCarbon,logPQNitrogen=cube[0],cube[1],cube[2],cube[3],cube[4],cube[5],cube[6]
    RayAmp,RaySlp,logPc,xRp=cube[7],cube[8],cube[9],cube[10]
   
    #for fixing unnessary parameters
    Rp= 1.359# Planet radius in Jupiter Radii--this will be forced to be 10 bar radius--arbitrary (scaling to this is free par)
    Rstar=1.155   #Stellar Radius in Solar Radii
    M = 0.690    #Mass in Jupiter Masses

    ##all values required by forward model go here--even if they are fixed
    x=np.array([Tirr, logKir,logg1, logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp*xRp, Rstar, M,RayAmp,RaySlp,logPc,xRp])
    gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.])
    y=fx(x,gas_scale) 

    y_mod=y[0]
  

    #y_mod=y_meas+x[1]
    loglikelihood=-0.5*np.sum((y_meas-y_mod)**2/err**2)
    return loglikelihood

########################
#reading in the data from text file
########################
#MCMC output file name--saved to pickle
outfile='MCMC.pic'

########################
#reading in the data from text file
########################
data=pickle.load(open("Data.pic",'rb'))

wlgrid=data[0]
y_meas=data[1]
err=data[2]

n_params=11


pymultinest.run(loglike, prior, n_params, outputfiles_basename='./chains/template_',resume=False, verbose=True,n_live_points=1000,importance_nested_sampling=False)
a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename='./chains/template_')
s = a.get_stats()


output=a.get_equal_weighted_posterior()
pickle.dump(output,open(outfile,"wb"))



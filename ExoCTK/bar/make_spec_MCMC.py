from .fm import *
import pickle
from joblib import Parallel,delayed


#load data
data=pickle.load(open("Data.pic",'rb'))
wlgrid=data[0]
y_meas=data[1]
err=data[2]


#load chains*************************************
fname1='MCMC.pic'
pic=pickle.load(open(fname1,'rb'))
chain=pic[:,:-1]
lnprob=pic[:,-1]


chi2=-2.*lnprob/wlgrid.shape[0]

Npars=chain.shape[1]

samples=chain
chi2flat=chi2
chi2best=np.min(chi2flat)
locbest=np.array((chi2flat==chi2best).nonzero())[0,0]
xbest=samples[locbest,:]

y_mod_arr=[]
y_hires_arr=[]
Parr=[]
Tarr=[]
H2O=[]
CH4=[]
CO=[]
CO2=[]
NH3=[]
N2=[]
HCN=[]
H2S=[]
PH3=[]
C2H2=[]
C2H6=[]
Na=[]
K=[]
TiO=[]
VO=[]
FeH=[]
H=[]
H2=[]
He=[]
Cloud=[]


#
NN=500#len(samples)
draws=np.random.randint(len(samples),size=NN)
#this is dumb...have to input all of the gas values here 
#b/c pymultinest won't allow us to do the fun fix parameters
#thing like in emcee
Tirr=samples[draws,0]
logKir=samples[draws,1]
logg1=samples[draws,2]
logMet=samples[draws,3]
logCtoO=samples[draws,4]
logPQCarbon=samples[draws,5]
logPQNitrogen=samples[draws,6]
Rp=samples[draws,0]*0.+1.359
Rstar=samples[draws,0]*0.+1.155
M=samples[draws,0]*0.+0.690
RayAmp=samples[draws,7]
RaySlp=samples[draws,8]
logPc=samples[draws,9]
xRp=samples[draws,10]


##all values required by forward model go here--even if they are fixed
x=np.array([Tirr, logKir,logg1, logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp*xRp, Rstar, M, RayAmp,RaySlp,logPc])


gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]) #can be made free params if desired (won't affect mmw)
fxarr=Parallel(n_jobs=5)(delayed(fx)(x[:,i],gas_scale) for i in range(NN))
#
wnocrop=fxarr[0][1]
for i in range(NN):
    print(i)
    y_mod=fxarr[i][0]
    P,T, H2Oarr, CH4arr, COarr,CO2arr, NH3arr, N2arr, H2Sarr, HCNarr, C2H2arr, C2H6arr, H2arr, Naarr, Karr=fxarr[i][2]
    y_mod_arr=np.concatenate([y_mod_arr,y_mod])
    Parr=np.concatenate([Parr,P])
    Tarr=np.concatenate([Tarr,T])
    H2O=np.concatenate([H2O,H2Oarr])
    CH4=np.concatenate([CH4,CH4arr])
    CO=np.concatenate([CO,COarr])
    CO2=np.concatenate([CO2,CO2arr])
    NH3=np.concatenate([NH3,NH3arr])
    HCN=np.concatenate([HCN,HCNarr])
    C2H2=np.concatenate([C2H2,C2H2arr])


y_mod_arr=y_mod_arr.reshape(NN,wlgrid.shape[0])

Parr=Parr.reshape(NN,P.shape[0])
Tarr=Tarr.reshape(NN,P.shape[0])
H2O=H2O.reshape(NN,P.shape[0])
CH4=CH4.reshape(NN,P.shape[0])
CO=CO.reshape(NN,P.shape[0])
CO2=CO2.reshape(NN,P.shape[0])
NH3=NH3.reshape(NN,P.shape[0])
HCN=HCN.reshape(NN,P.shape[0])
C2H2=C2H2.reshape(NN,P.shape[0])



xarr=np.zeros((NN,samples.shape[1]))
for i in range(NN): xarr[i,:]=samples[draws[i],:]

array=[Parr, Tarr, H2O,CH4,CO,CO2,NH3,HCN,C2H2]
pickle.dump([y_mod_arr,wnocrop,xarr,array], open("y_arr_full.pic","wb") )


pdb.set_trace()





pdb.set_trace()







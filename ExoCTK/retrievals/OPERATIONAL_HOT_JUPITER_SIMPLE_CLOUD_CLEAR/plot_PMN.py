from numpy.random import * 
import numpy as np
import pickle
import pylab
from matplotlib.pyplot import *
from scipy.io.idl import readsav #USAGE: >>> readsave('file.save')
import pdb
import corner
from scipy import interp
rc('font',family='serif')


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


fname='HOT_JUPITER_CLEAR_RayHz_Cld_Top'

chi2=-2.*lnprob/wlgrid.shape[0]

Npars=chain.shape[1]

samples=chain
chi2flat=chi2
chi2best=np.min(chi2flat)
locbest=np.array((chi2flat==chi2best).nonzero())[0,0]
xbest=samples[locbest,:]


#plotting stair-pairs plot via traingle.py******************            
priorlow=np.array( [300. ,-3.,  -3., -2.,-2., -6.5, -6.5,-2,0,-6,0 ])
priorhigh=np.array([3000., 0.,   1.,  3., 2.,  1.5,  1.5, 4,4, 1.5,1])

Npars=len(priorlow)
ext=np.zeros([2,Npars])
for i in range(Npars):
	med=np.percentile(samples[:,i],50)
	stdev=np.std(samples[:,i])
	ext[0,i]=med-4.*stdev
	if ext[0,i] < priorlow[i]:
		ext[0,i]=priorlow[i]
	ext[1,i]=med+4.*stdev
	if ext[1,i] > priorhigh[i]:
		ext[1,i]=priorhigh[i]
ext=ext.T
#ext[:,0]=priorlow
#ext[:,1]=priorhigh


corner.corner(samples,labels=['Tirr', 'logKir','logg1', 'logMet', 'logCtoO', 'logPQC','logPQN','RayAmp','RaySlp','logPc','xRp' ], bins=25,plot_datapoints='False',quantiles=[.16,0.5,.84],show_titles='True',plot_contours='True',truths=[1200,-1.5,-1,0,-0.26,-5,-5,0,4,1.5,1],levels=(1.-np.exp(-(1)**2/2.),1.-np.exp(-(2)**2/2.),1.-np.exp(-(3)**2/2.)))

savefig(fname+"_stair_pairs.pdf",format='pdf')
close()


#Plotting spectra spread**********************
y_arr=pickle.load(open('y_arr_full.pic','rb'))
y_mod_arr=y_arr[0]
wnocrop=y_arr[1]
Parr, Tarr, H2O,CH4,CO,CO2,NH3,HCN,C2H2=y_arr[3]
NN=y_mod_arr.shape[0]

y_median=np.zeros(wnocrop.shape[0])
y_high_1sig=np.zeros(wnocrop.shape[0])
y_high_2sig=np.zeros(wnocrop.shape[0])
y_low_1sig=np.zeros(wnocrop.shape[0])
y_low_2sig=np.zeros(wnocrop.shape[0])

for i in range(wnocrop.shape[0]):
    percentiles=np.percentile(y_mod_arr[:,i],[4.55, 15.9, 50, 84.1, 95.45])
    y_low_2sig[i]=percentiles[0]*100
    y_low_1sig[i]=percentiles[1]*100
    y_median[i]=percentiles[2]*100
    y_high_1sig[i]=percentiles[3]*100
    y_high_2sig[i]=percentiles[4]*100


fig1, ax=subplots()
ymax=1.01*np.max(y_meas)*100
ymin=0.99*np.min(y_meas)*100
xlabel('$\lambda$ ($\mu$m)',size='xx-large')
ylabel('(R$_{p}$/R$_{\star}$)$^2$ [$\%$]',size='xx-large')
ax.fill_between(1E4/wnocrop,y_low_2sig,y_high_2sig,facecolor='r',alpha=0.1,edgecolor='None')  
ax.fill_between(1E4/wnocrop,y_low_1sig,y_high_1sig,facecolor='r',alpha=1.,edgecolor='None')  
ax.set_xscale('log')
ax.set_xticks([0.6,0.8,1,3,5,7,9,11])
ax.axis([0.6,11,ymin*0.98,ymax*1.02])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(length=10,width=1,labelsize='large',which='major')
ax.plot(1E4/wnocrop,y_median,color='b')
minorticks_on()
ax.errorbar(1E4/wnocrop, y_meas*100, yerr=err*100,xerr=None,fmt='D',color='black',alpha=0.5)
savefig(fname+'_spectra.pdf',format='pdf')
close()


#plotting TP and cloud
#TP
P=Parr[0,:]
T_median=np.zeros(P.shape[0])
T_high_1sig=np.zeros(P.shape[0])
T_high_2sig=np.zeros(P.shape[0])
T_low_1sig=np.zeros(P.shape[0])
T_low_2sig=np.zeros(P.shape[0])

for i in range(P.shape[0]):
    percentiles=np.percentile(Tarr[:,i],[4.55, 15.9, 50, 84.1, 95.45])
    T_low_2sig[i]=percentiles[0]
    T_low_1sig[i]=percentiles[1]
    T_median[i]=percentiles[2]
    T_high_1sig[i]=percentiles[3]
    T_high_2sig[i]=percentiles[4]

#TP
fill_betweenx(P,T_low_2sig,T_high_2sig,facecolor='r',edgecolor='None',alpha=0.1)
fill_betweenx(P,T_low_1sig,T_high_1sig,facecolor='r',edgecolor='None',alpha=1.)
semilogy(T_median,P,'b-',lw='2',label='TP')
xlabel('Temperature [K]',fontsize=20)
axis([0,3000,30,1E-7])
ylabel('Pressure [bar]',fontsize=20)
savefig(fname+'_TP.pdf',format='pdf')
close()

#plotting chemistry
axis([1E-8,1,30,1E-7])
semilogy()
for i in range(H2O.shape[0]):semilogx(H2O[i,:],P,alpha=0.05,color='blue')
for i in range(H2O.shape[0]):semilogx(CH4[i,:],P,alpha=0.05,color='black')
for i in range(H2O.shape[0]):semilogx(CO[i,:],P,alpha=0.05,color='green')
for i in range(H2O.shape[0]):semilogx(CO2[i,:],P,alpha=0.05,color='orange')
for i in range(H2O.shape[0]):semilogx(NH3[i,:],P,alpha=0.05,color='red')
semilogx( np.median(H2O,axis=0),P,alpha=0.5,color='blue',label='H2O')
semilogx( np.median(CH4,axis=0),P,alpha=0.5,color='black',label='CH4')
semilogx( np.median(CO,axis=0),P,alpha=0.5,color='green',label='CO')
semilogx( np.median(CO2,axis=0),P,alpha=0.5,color='orange',label='CO2')
semilogx( np.median(NH3,axis=0),P,alpha=0.5,color='red',label='NH3')
legend()
xlabel('Gas Volume Mixing Ratio',size='xx-large')
ylabel('Pressure [bar]',size='xx-large')
savefig(fname+'_Gas_VMR.pdf',format='pdf')
close()

pdb.set_trace()

from fm import *
import pickle
from matplotlib.ticker import FormatStrFormatter
rc('font',family='serif')



def spectrum(Rp= 1.93, Rstar=1.2, M = 0.5, Tirr=1700, logKir=-1.5, logg1=-1, logMet=0., logCtoO=-0.26, logPQCarbon=-5, logPQNitrogen=-5, RayAmp=-5., RaySlp=4., logPc=1.5):
	# #the parameters
	# #planet/star system params--typically not free parameters in retrieval
	# Rp= 1.93# Planet radius in Jupiter Radii--this will be forced to be 10 bar radius--arbitrary (scaling to this is free par)
	# Rstar=1.2   #Stellar Radius in Solar Radii
	# M = 0.5    #Mass in Jupiter Masses
	# #TP profile params (3--Guillot 2010, Parmentier & Guillot 2013--see Line et al. 2013a for implementation)
	# Tirr=1700 #terminator **isothermal** temperature--if full redistribution this is equilibrium temp
	# logKir=-1.5  #TP profile IR opacity controlls the "vertical" location of the gradient
	# logg1=-1     #single channel Vis/IR opacity. Controls the delta T between deep T and TOA T
	# #Composition parameters---assumes "chemically consistnat model" described in Kreidberg et al. 2015
	# logMet=0. #.   #Metallicity relative to solar log--solar is 0, 10x=1, 0.1x = -1 used -1.01*log10(M)+0.6
	# logCtoO=-0.26  #log C-to-O ratio: log solar is -0.26
	# logPQCarbon=-5  #CH4, CO, H2O Qunech pressure--forces CH4, CO, and H2O to constant value at quench pressure value
	# logPQNitrogen=-5  #N2, NH3 Quench pressure--forces N2 and NH3 to ""  --ad hoc for chemical kinetics--reasonable assumption
	# #generic cloud parameters--Rayleigh Haze(amplitude, slope) and Hard Cloud Top Pressure
	# RayAmp=-5. #log Rayleigh Haze Amplitude (relative to H2)
	# RaySlp=4. #haze slope--4 is Rayeigh, 0 is "gray" or flat.
	# logPc=1.5 #log of the hard gray cloud top pressure
	
	
	
	#seting up input state vector. Must be in this order as indicies are hard wired in fx inside fm
	           #  0    1        2       3     4      5              6           7     8    9     10     11     12        
	x=np.array([Tirr, logKir,logg1, logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp, Rstar, M, RayAmp, RaySlp,logPc])
	
	
	#calling forward model
	#thermochemical gas profile scaling factors
	# 0   1    2    3   4    5    6     7    8    9   10    11   12   13    14   15   16   17   18
	#H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He
	gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]) #can be made free params if desired (won't affect mmw)
	y_mod,wnocrop,atm=fx(x,gas_scale)  #returns model spectrum, wavenumber grid, and vertical abundance profiles from chemistry
	
	
	#read in from noise model here
	err=np.zeros(len(wnocrop))
	err[:]=50*1E-6
	y_meas=np.zeros(len(err))
	
	#adding gaussian noise 
	for i in range(y_meas.shape[0]): y_meas[i]=y_mod[i]+np.random.randn(1)*err[i]
	
	#computing chi-square of random noise instance
	print np.sum((y_meas-y_mod)**2/err**2)/len(y_meas)
	
	
	#dumping pickles--model then noised up data data
	output=[1E4/wnocrop, y_mod]
	pickle.dump(output,open("Model.pic","wb"))  #spectral model to be noised up by instrument noise model
	
	output=[1E4/wnocrop, y_meas,err]  #noised up "synthetic" spectrum
	pickle.dump(output,open("Data.pic","wb"))
	
	
	#plotting-------------------------------------------------------------------
	ymin=0.98*np.min(y_mod)*100
	ymax=1.02*np.max(y_mod)*100
	#ymin=2.9
	#ymax=3.05
	
	
	#plotting structure (TP, mixing ratios, cloud)-----------------------
	# atm=np.array([P,T, H2Oarr, CH4arr, COarr,CO2arr, NH3arr, N2arr, H2Sarr, HCNarr, C2H2arr, C2H6arr, H2arr, Naarr, Karr])
	P=atm[0]
	T=atm[1]
	H2O=atm[2]
	CH4=atm[3]
	CO=atm[4]
	CO2=atm[5]
	NH3=atm[6]
	HCN=atm[9]
	C2H2=atm[10]
	Cld=atm[-1]
	
	#mixing ratios
	fig2, ax1=subplots()
	ax1.semilogx(H2O,P,'b',ls='--',lw=2,label='H2O')
	ax1.semilogx(CH4,P,'black',ls='--',lw=2,label='CH4')
	ax1.semilogx(CO,P,'g',ls='--',lw=2,label='CO')
	ax1.semilogx(CO2,P,'orange',ls='--',lw=2,label='CO2')
	ax1.semilogx(NH3,P,'darkblue',ls='--',lw=2,label='NH3')
	ax1.semilogx(HCN,P,'r',ls='--',lw=2,label='HCN')
	ax1.semilogx(C2H2,P,'c',ls='--',lw=2,label='C2H2')
	ax1.set_xlabel('Mixing Ratio',fontsize=20)
	ax1.axis([1E-9,1,100,1E-6])
	
	#TP
	ax2=ax1.twiny()
	ax2.semilogy(T,P,'r-',lw='2',label='TP')
	ax2.set_xlabel('Temperature [K]',color='r',fontsize=20)
	ax2.axis([0,3000,100,1E-6])
	for tl in ax2.get_xticklabels(): tl.set_color('r')
	ax1.set_ylabel('Pressure [bar]',fontsize=20)
	ax1.legend(loc=4,frameon=False)
	
	
	#plotting spectrum-------------------------------------
	fig1, ax=subplots()
	ax.plot(1E4/wnocrop, y_mod*100.)
	xlabel('$\lambda$ ($\mu$m)',fontsize=18)
	ylabel('(R$_{p}$/R$_{*}$)$^{2} \%$',fontsize=18)
	minorticks_on()
	ax.errorbar(1E4/wnocrop, y_meas*100, yerr=err*100,xerr=None,fmt='D',color='black',alpha=0.5)
	ax.set_xscale('log')
	ax.set_xticks([0.6,0.8,1,3,5,7,9,11])
	ax.axis([0.6,11,ymin*0.98,ymax*1.02])
	ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.tick_params(length=10,width=1,labelsize='large',which='major')
	show()
	
	
	
	
	##########################TESTING, PLAYING TO GAIN INTUITION######################################
	#------twiddeling--testing cloud param values---can test other params to by forcing x index or gas_scale values to new value
	# 0   1    2    3   4    5    6     7    8    9   10    11   12   13    14   15   16   17   18  #19 gases (PH3, N2, and H are not opacity sources)
	#H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He
	# gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.])*0. #all gases zeroed out--just want to test effect of cloud
	gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.])*0. #all gases zeroed out--just want to test effect of cloud
	
	
	##plotting spectrum-------------------------------------
	#fig1, ax=subplots()
	#ax.plot(1E4/wnocrop, y_mod*100.)
	#
	## pdb.set_trace()
	#
	## Clean this plot up! Choose opacity sources or just plot top 5
	#for n,i in enumerate(gas_scale):  #looping trhough 10 param values
	#    gas_scale[n]=1
	#    y_mod2,wnocrop2,junk=fx(x,gas_scale)
	#    ax.plot(1E4/wnocrop2, y_mod2*100.,color='blue',alpha=0.2+0.2*i)
	#    gas_scale[n] = 0.
	#    # y_mod2,wnocrop2,junk=fx(x,gas_scale)
	#    # ax.plot(1E4/wnocrop2, y_mod2*100.,color='blue',alpha=0.2+0.2*i)
	#
	#xlabel('$\lambda$ ($\mu$m)',fontsize=18)
	#ylabel('(R$_{p}$/R$_{*}$)$^{2} \%$',fontsize=18)
	#minorticks_on()
	#ax.errorbar(1E4/wnocrop, y_meas*100, yerr=err*100,xerr=None,fmt='D',color='black',alpha=0.5)
	#ax.set_xscale('log')
	#ax.set_xticks([0.6,0.8,1,3,5,7,9,11])
	#ax.axis([0.6,11,ymin*0.98,ymax*1.02])
	#ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	#ax.tick_params(length=10,width=1,labelsize='large',which='major')
	#show()
	#plt.savefig('./test.png')
	#
	#pdb.set_trace()
	
	#truth vector
	#pickle.dump(x,open('truth.pic','wb'))
	
	
	
	
	
	
	
	

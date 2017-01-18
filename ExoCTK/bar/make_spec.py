from fm import *
import pickle
from matplotlib.ticker import FormatStrFormatter

def spectrum(Rp= 1.93, Rstar=1.2, Mp=0.5, Tirr=1700, logKir=-1.5, logg1=-1, logMet=0., logCtoO=-0.26, 
             logPQC=-5, logPQN=-5, RayAmp=-5., RaySlp=4., logPc=1.5, plot=True):
    """
    Description...
    
    Parameters
    ----------
    Rp: float
        The planet radius [RJup]
    Rstar: float
        The stellar radius [RSun]
    Mp: float
        The planet mass [MJup]
    Tirr: int
        The isothermal temperature at the terminator [K].
        If full redistribution, this is equilibrium temperature.
    logKir: float
        The TP profile IR opacity, the "vertical" location of the gradient
    logg1: float
        The single channel Vis/IR opacity. Controls the delta T between deep T and TOA T
    logMet: float
        The metallicity relative to solar log--solar is 0, 10x=1, 0.1x = -1 used -1.01*log10(M)+0.6
    logCtoO: float
        The log(C/O) ratio. Log solar is -0.26
    logPQC: float
        CH4, CO, H2O Qunech pressure--forces CH4, CO, and H2O to 
        constant value at quench pressure value
    logPQN: float
        N2, NH3 Quench pressure--forces N2 and NH3 to ""  
        --ad hoc for chemical kinetics-- reasonable assumption
    RayAmp: float
        log Rayleigh Haze Amplitude (relative to H2)
    RaySlp: float
        haze slope--4 is Rayeigh, 0 is "gray" or flat
    logPc: float
        log of the hard gray cloud top pressure
    plot: bool
        Plot the PT profiles
    
	"""
	#seting up input state vector. Must be in this order as indicies are hard wired in fx inside fm
	           #  0    1        2       3     4      5              6           7     8    9     10     11     12        
	x = np.array([Tirr, logKir,logg1, logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp, Rstar, M, RayAmp, RaySlp,logPc])
	
	
	#calling forward model
	#thermochemical gas profile scaling factors
	# 0   1    2    3   4    5    6     7    8    9   10    11   12   13    14   15   16   17   18
	#H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He
	gas_scale = np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]) #can be made free params if desired (won't affect mmw)
	
    # Return model spectrum, wavenumber grid, and vertical abundance profiles from chemistry
    y_mod, wnocrop, atm = fx(x,gas_scale)
	
	# Read in from noise model here
	err = np.zeros(len(wnocrop))
	err[:] = 50*1E-6
	y_meas = np.zeros(len(err))
	
	# Adding gaussian noise 
	for i in range(y_meas.shape[0]): 
        y_meas[i] = y_mod[i]+np.random.randn(1)*err[i]
	
	# Compute chi-square of random noise instance
	print(np.sum((y_meas-y_mod)**2/err**2)/len(y_meas))
	
	# Dump noised up data
	output = [1E4/wnocrop, y_mod]
	pickle.dump(output,open("Model.pic","wb"))
	
    # Dump noised up "synthetic" spectrum
	output=[1E4/wnocrop, y_meas, err]
	pickle.dump(output,open("Data.pic","wb"))
	
	# Plots
    if plot:
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
	#gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.])*0. #all gases zeroed out--just want to test effect of cloud
	
	
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
	
	
	
	
	
	
	
	

import h5py
import json
import numpy as np
import os
import pysynphot as psyn
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import time

class PlanetarySystem():
    """Class object defining the planetary system you wish to model with Chimera"""

    def __init__(self, json_file=None):
        """Initialize the class object."""

        if json_file:
            with open(json_file, "r") as read_file:
                self.json_data = json.load(read_file)
            for key in self.json_data:
                setattr(self, key, self.json_data[key])
            self.gas_scale = np.array(list(self.gas_abundances.values()))
        else:
            self.observatory = None
            self.wavenumber_min = None
            self.wavenumber_max = None
            self.chemistry_parameters = {}
            self.stellar_parameters = {}
            self.planetary_parameters = {}
            self.gas_abundances = {}
            self.cloud_parameters = {}

class GenerateModel():

    def __init__(self, cross_sections):
        """Initialize the class object."""

        self.cross_sections = cross_sections
        self.make_stellar_model()
        self.load_transmission_spectrum('/Users/mfix/Desktop/exoctk/exoctk/atmospheric_retrievals/w43b_trans.txt')
        self.fx_trans_free()


    def load_transmission_spectrum(self, transmission_file):
        self.wlgrid, self.y_meas, self.err=np.loadtxt(transmission_file).T


    def make_plot(self):
        """Plot Model"""
        P,T, H2O, CH4,CO,CO2,NH3,Na,K,TiO,VO,C2H2,HCN,H2S,FeH,H2,He,H,e, Hm,qc,r_eff,f_r=atm

        fig2, ax1=plt.subplots()
        #feel free to plot whatever you want here....
        ax1.semilogx(H2O,P,'b',ls='--',lw=2,label='H2O')
        ax1.semilogx(CH4,P,'black',ls='--',lw=2,label='CH4')
        ax1.semilogx(CO,P,'g',ls='--',lw=2,label='CO')
        ax1.semilogx(CO2,P,'orange',ls='--',lw=2,label='CO2')
        ax1.semilogx(NH3,P,'darkblue',ls='--',lw=2,label='NH3')
        ax1.semilogx(Na,P,'b',lw=2,label='Na')
        ax1.semilogx(K,P,'g',lw=2,label='K')
        ax1.semilogx(TiO,P,'k',lw=2,label='TiO')
        ax1.semilogx(VO,P,'orange',lw=2,label='VO')
        ax1.semilogx(qc,P,'gray',lw=1,ls='--',label='Cond. VMR.')  #<---- A&M Cloud Condensate VMR profile (not droplets)

        ax1.set_xlabel('Mixing Ratio',fontsize=20)
        ax1.set_ylabel('Pressure [bar]',fontsize=20)
        ax1.semilogy()
        ax1.legend(loc=4,frameon=False)
        ax1.axis([1E-9,1,100,1E-7])

        #plotting TP profile on other x-axis
        ax2=ax1.twiny()
        ax2.semilogy(T,P,'r-',lw='4',label='TP')
        ax2.set_xlabel('Temperature [K]',color='r',fontsize=20)
        ax2.axis([0.8*T.min(),1.2*T.max(),100,1E-6])
        for tl in ax2.get_xticklabels(): tl.set_color('r')
        ax2.legend(loc=1,frameon=False)

        plt.savefig('./plots/atmosphere_transmission_WFC3_FREE.pdf',fmt='pdf')
        plt.show()
        plt.close()


    def retrieve(self, method):
        """Perform the atmopsheric retrieval"""
        return None


    def make_stellar_model(self):
        """
        Make stellar using pysynphot

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        database = self.cross_sections.stellar_parameters['stellar_db']
        temp = self.cross_sections.stellar_parameters['temp']
        logMH = self.cross_sections.stellar_parameters['logmh']
        logg = self.cross_sections.stellar_parameters['logg']
        outfile = self.cross_sections.stellar_parameters['stellar_file']

        sp = psyn.Icat(database, temp, logMH, logg)
        sp.convert('flam') # initial units erg/cm^2/s/ Angst
        wave=sp.wave*1e-10  # in meters
        fluxmks = sp.flux*1e-7*1e4*1e10 # in W/m2/m
        f = h5py.File(outfile,'w')
        f.create_dataset('lambdastar',data=wave)
        f.create_dataset('Fstar0',data=fluxmks)
        f.close()


    def fx_trans_free(self):
        """Transmission spectrscopy with free chemistry
        Parameters
        ----------
        x : list 
            See tutorials for description of list and order of list 
        wlgrid : ndarray
            Array to regrid final specturm on (micron)
        gas_scale : ndarray
            array to scale mixing ratio of gases 
        xsects : list 
            cross section array from `xsecs` function 
        Returns
        -------
        y_binned,F,wno,chemarr
        binned array of spectrum, high res spectrum, og wavenumber grid, chemistry array
        which includes: 
        chemarr = [P,T, H2Oarr, CH4arr,COarr,CO2arr,NH3arr,Naarr,Karr,TiOarr,VOarr,C2H2arr,HCNarr,H2Sarr,FeHarr,H2arr,Hearr,Harr, earr, Hmarr,qc,r_eff,f_r])
        """

        # Unpacking Guillot 2010 TP profile params (3 params)
        Tirr = self.cross_sections.planetary_parameters['Tirr']
        logKir = self.cross_sections.planetary_parameters['logKir']
        logg1 = self.cross_sections.planetary_parameters['logg1']
        Tint = self.cross_sections.planetary_parameters['Tint']

        # Unpacking Chemistry Parms
        Met=10.** self.cross_sections.chemistry_parameters['Metallicity']  #metallicity
        CtoO=10.** self.cross_sections.chemistry_parameters['CtoO'] #C/O
        logPQC = self.cross_sections.chemistry_parameters['logPQC']  #carbon quench pressure
        logPQN = self.cross_sections.chemistry_parameters['logPQN']  #nitrogen quench pressure

        # unpacking planet params
        Rp = self.cross_sections.planetary_parameters['Rp']  #planet radius (in jupiter)
        Rstar = self.cross_sections.stellar_parameters['Rstar']   #stellar radius (in solar)
        M = self.cross_sections.planetary_parameters['M']   #planet mass (in jupiter)

        # unpacking and converting A&M cloud params
        Kzz=10** self.cross_sections.cloud_parameters['logKzz'] *1E-4  #Kzz for A&M cloud
        fsed = self.cross_sections.cloud_parameters['fsed']  #sedimentation factor for A&M cloud
        Pbase=10.** self.cross_sections.cloud_parameters['logPbase']  #cloud top pressure
        Cld_VMR=10** self.cross_sections.cloud_parameters['logCldVMR']  #Cloud Base Condensate Mixing ratio

        # unpacking and converting simple cloud params
        CldOpac=10** self.cross_sections.cloud_parameters['logKcld']
        RayAmp=10** self.cross_sections.cloud_parameters['logRayAmp']
        RaySlp = self.cross_sections.cloud_parameters['RaySlope']
        
        # Setting up atmosphere grid
        self.atmosphere_grid = {}
        self.atmosphere_grid['logP'] = np.arange(-6.8,1.5,0.1)+0.1
        self.atmosphere_grid['P'] = 10.0** self.atmosphere_grid['logP']
        self.atmosphere_grid['g0'] = 6.67384E-11*M*1.898E27/(Rp*71492.*1.E3)**2
        self.atmosphere_grid['kv'] = 10.**(logg1+logKir)
        self.atmosphere_grid['kth'] = 10.**logKir
        self.atmosphere_grid['alpha'] = 0.5

        tempurature, pressure = self.profile_tempeture_and_pressure(self.atmosphere_grid['kv'], self.atmosphere_grid['kv'], self.atmosphere_grid['kth'])
        self.T = sp.interp(self.atmosphere_grid['logP'],np.log10(pressure),tempurature)
        t1 = time.time()

        Tavg = 0.5*(self.T[1:]+self.T[:-1])
        Pavg = 0.5*(self.atmosphere_grid['P'][1:] + self.atmosphere_grid['P'][:-1])

        #interpolation chem
        # No need to unpack, trying to avoid passing big arrays around and just use class members in assignment.
        # logCtoO, logMet, Tarr, logParr, loggas = self.cross_sections.xsecarr[9:]  #see xsects_HST/JWST routine...    
        Ngas = self.cross_sections.loggas.shape[-2]
        gas = np.zeros((Ngas,len(self.atmosphere_grid['P'])))+1E-20
        #             H2O   CH4   CO    CO2    NH3   N2    HCN   H2S   PH3 C2H2   C2H6   Na     K    TiO    VO     FeH     H    H2   He   e  H-
        mu = np.array([18.02,16.04,28.01,44.01,17.03,28.01,27.02,34.08,  34.,26.04, 30.07,22.99, 39.1, 63.87, 66.94, 56.85, 1.01, 2.02, 4.0,0.,1.01,0 ])

        self.cross_sections.gas_scale*=1E0
        for i in range(Ngas): gas[i,:] = 10**self.cross_sections.gas_scale[i]

        H2Oarr, CH4arr, COarr, CO2arr, NH3arr, N2arr, HCNarr, H2Sarr,PH3arr,C2H2arr, C2H6arr, Naarr, Karr, TiOarr, VOarr,FeHarr, Harr,H2arr, Hearr,earr, Hmarr,mmw = gas

        H2He = 1.-np.sum(gas,axis=0)
        frac = 0.176471
        H2arr = H2He/(1.+frac)
        Hearr = frac*H2arr
        gas[-5] = H2arr
        gas[-4] = Hearr

        mmw[:] = gas.T.dot(mu)

        #ackerman & Marley cloud model here
        mmw_cond = 100.39 # molecular weight of condensate (in AMU)  MgSiO3=100.39
        rho_cond = 3250 # density of condensate (in kg/m3)           MgSiO3=3250.
        rr = 10** (np.arange(-2,2.6,0.1)) # Droplet radii to compute on: MUST BE SAME AS MIE COEFF ARRAYS!!!!!!!!! IF YOU CHANGE THIS IT WILL BREAK
        qc = self.cloud_profile(fsed, Cld_VMR, self.atmosphere_grid['P'], Pbase)
        r_sed, r_eff, r_g, f_r = self.particle_radius(fsed, Kzz, mmw, self.T, self.atmosphere_grid['P'], 
                                            self.atmosphere_grid['g0'], rho_cond,mmw_cond, qc, rr*1E-6)

        #10.1  #reference pressure bar-keep fixed
        #computing transmission spectrum-----------
        Pref=1.1
        
        self.spec = self.tran(self.cross_sections.xsecarr, self.T, self.atmosphere_grid['P'], mmw, Pref,
                CldOpac, H2Oarr, CH4arr,COarr,CO2arr,NH3arr,Naarr,Karr,TiOarr,VOarr,C2H2arr,HCNarr,H2Sarr,
                FeHarr,Harr,earr,Hmarr,H2arr,Hearr,RayAmp,RaySlp,f_r, M, Rstar, Rp)

        self.wno = self.spec[0]
        self.F = self.spec[1]

        self.y_binned, _ = self.instrument_tran_non_uniform(self.wlgrid, self.wno, self.F)

        self.chemarr = np.array([self.atmosphere_grid['P'], self.T, H2Oarr, CH4arr, 
                        COarr, CO2arr, NH3arr, Naarr, Karr, TiOarr, VOarr, C2H2arr, 
                        HCNarr, H2Sarr, FeHarr, H2arr, Hearr, Harr, earr, Hmarr, qc, r_eff, f_r])


    def instrument_tran_non_uniform(self, wlgrid,wno, Fp):
        """
        Rebins transmission spectra on new wavelength grid 
        Parameters
        ----------
        wlgrid : ndarray
            New wavelength grid 
        wno : ndarray
            Old wavenumber grid 
        Fp : ndarray
            (rp/rs)^2
        Returns
        -------
        Fratio_int
            new regridded spectrum
        Fp
            old high res spectrum
        """
        if isinstance(wlgrid,int):
            return Fp, Fp

        else:
            szmod=wlgrid.shape[0]
            delta=np.zeros(szmod)
            Fratio=np.zeros(szmod)
            for i in range(szmod-1):
                delta[i]=wlgrid[i+1]-wlgrid[i]  
            delta[szmod-1]=delta[szmod-2] 

            for i in range(szmod-1):
                i=i+1
                loc=np.where((1E4/wno > wlgrid[i]-0.5*delta[i-1]) & (1E4/wno < wlgrid[i]+0.5*delta[i]))
                Fratio[i]=np.mean(Fp[loc])

            loc=np.where((1E4/wno > wlgrid[0]-0.5*delta[0]) & (1E4/wno < wlgrid[0]+0.5*delta[0]))
            Fratio[0]=np.mean(Fp[loc])

            Fratio_int=Fratio
            return Fratio_int, Fp


    def kcoeff_interp(self, logPgrid, logTgrid, logPatm, logTatm, wnogrid, kcoeff):
        """
        This routine interpolates the correlated-K tables
        to the appropriate atmospheric P & T for each wavenumber and 
        g-ordinate for each gas. It uses a standard bi-linear 
        interpolation scheme.
        
        Parameters
        ---------- 
        logPgrid : ndarray
            pressure grid (log10) on which the CK coeff's are pre-computed (Npressure points)
        logTgrid : ndarray
            temperature grid (log10) on which the CK coeffs are pre-computed (Ntemperature points)
        logPatm : ndarray 
            atmospheric pressure grid (log10) 
        logTatm : ndarray 
            atmospheric temperature grid (log10)
        wnogrid : ndarray 
            CK wavenumber grid (Nwavenumber points) (actually, this doesn't need to be passed...it does nothing here...)
        kcoeff : ndarray 
            massive CK coefficient array (in log10)--Ngas x Npressure x Ntemperature x Nwavenumbers x Ngordinates
        
        Returns
        ------- 
        kcoeff_int 
            the interpolated-to-atmosphere CK coefficients (in log). 
            This will be Nlayers x Nwavenumber x Ngas x Ngordiantes
        """

        Ng, NP, NT, Nwno, Nord = self.cross_sections.xsecarr.shape

        Natm=len(logTatm)

        kcoeff_int=np.zeros((Natm,Nwno,Ng,Nord))

        for i in range(Natm):  #looping through atmospheric layers

            y=logPatm[i]
            x=logTatm[i]

            p_ind_hi=np.where(logPgrid>=y)[0][0]
            p_ind_low=np.where(logPgrid<y)[0][-1]
            T_ind_hi=np.where(logTgrid>=x)[0][0]
            T_ind_low=np.where(logTgrid<x)[0][-1]

            y2=logPgrid[p_ind_hi]
            y1=logPgrid[p_ind_low]
            x2=logTgrid[T_ind_hi]
            x1=logTgrid[T_ind_low]
        
            for j in range(Ng): #looping through gases
                for k in range(Nwno): #looping through wavenumber
                    for l in range(Nord): #looping through g-ord
                        arr=kcoeff[j,:,:,k,l]
                        Q11=arr[p_ind_low,T_ind_low]
                        Q12=arr[p_ind_hi,T_ind_low]
                        Q22=arr[p_ind_hi,T_ind_hi]
                        Q21=arr[p_ind_low,T_ind_hi]
                        fxy1=(x2-x)/(x2-x1)*Q11+(x-x1)/(x2-x1)*Q21
                        fxy2=(x2-x)/(x2-x1)*Q12+(x-x1)/(x2-x1)*Q22
                        fxy=(y2-y)/(y2-y1)*fxy1 + (y-y1)/(y2-y1)*fxy2
                        kcoeff_int[i,k,j,l]=fxy

        return kcoeff_int

    def tran(self, xsects, T, P, mmw,Ps,CldOpac,alphaH2O,alphaCH4,alphaCO,alphaCO2,alphaNH3,alphaNa,alphaK,alphaTiO,alphaVO, alphaC2H2, alphaHCN, alphaH2S,alphaFeH,fH,fe,fHm,fH2,fHe,amp,power,f_r,M,Rstar,Rp):
        """Runs all requisite routins for transmisison spectroscopy
        Returns
        -------
        wno
            Wavenumber grid 
        F
            (rp/rs)^2
        Z
            height above ref pressure
        """
        t1=time.time()

        #renaming variables, bbecause why not
        fH2=fH2
        fHe=fHe
        fH2O=alphaH2O
        fCH4=alphaCH4
        fCO=alphaCO
        fCO2=alphaCO2
        fNH3=alphaNH3
        fNa=alphaNa
        fK=alphaK
        fTiO=alphaTiO
        fVO=alphaVO
        fC2H2=alphaC2H2
        fHCN=alphaHCN
        fH2S=alphaH2S
        fFeH=alphaFeH
        mmw=mmw
        #pdb.set_trace()
        #Na and K are fixed in this model but can be made free parameter if desired
        #If T < 800 K set these equal to 0!!!--they condense out below this temperature (roughly)
    
        
        Fractions = np.array([fH2*fH2,fHe*fH2,fH2O, fCH4, fCO, fCO2, fNH3,fNa,fK,fTiO,fVO,fC2H2,fHCN,fH2S,fFeH,fH*fe, fHm])  #gas mole fraction profiles
                            #H2Ray, HeRay  Ray General,
        Frac_Cont = np.array([fH2,fHe,fH2*0.+1.,fH2*0.+1])  #continuum mole fraction profiles
        Frac_Cont=np.concatenate((Frac_Cont, f_r),axis=0)

        #Load measured cross-sectional values and their corresponding
        #T,P,and wno grids on which they were measured
        
        # Stop unpacking array, use class members for variable assignment
        #Pgrid, Tgrid, wno, gord, wts, xsecarr, radius, Mies, Fstar=xsects[0:9]
        Pgrid = self.cross_sections.P
        Tgrid = self.cross_sections.T
        wno = self.cross_sections.wno
        gord = self.cross_sections.g
        wts = self.cross_sections.wts
        xsecarr = self.cross_sections.xsecarr
        radius = self.cross_sections.radius
        Mies = self.cross_sections.mies_arr
        Fstar = self.cross_sections.Fstar
        
        #Calculate Temperature, Pressure and Height grids on which
        #transmissivity will be computed
        n = len(P)
        nv = len(wno)
        
        
        Z=np.zeros(n)  #level altitudes
        dZ=np.zeros(n)  #layer thickness array
        r0=Rp*71492.*1.E3  #converting planet radius to meters
        mmw=mmw*1.660539E-27  #converting mmw to Kg
        kb=1.38E-23
        G=6.67428E-11
        M=M*1.89852E27

        
        #Compute avg Temperature at each grid
        Tavg = np.array([0.0]*(n-1))
        Pavg = np.array([0.0]*(n-1))
        for z in range(n-1):
            Pavg[z] = np.sqrt(P[z]*P[z+1])
            Tavg[z] = sp.interp(np.log10(Pavg[z]),sp.log10(P),T)
        #create hydrostatic altitutde grid from P and T
        Phigh=P.compress((P>Ps).flat)  #deeper than reference pressure
        Plow=P.compress((P<=Ps).flat)   #shallower than reference pressure
        for i in range(Phigh.shape[0]):  #looping over levels above ref pressure
            i=i+Plow.shape[0]-1
            g=G*M/(r0+Z[i])**2#g0*(Rp/(Rp+Z[i]/(69911.*1E3)))**2
            H=kb*Tavg[i]/(mmw[i]*g)  #scale height
            dZ[i]=H*np.log(P[i+1]/P[i]) #layer thickness, dZ is negative
            Z[i+1]=Z[i]-dZ[i]   #level altitude
            #print(P[i], H/1000, Z[i]/1000, g)
        for i in range(Plow.shape[0]-1):  #looping over levels below ref pressure
            i=Plow.shape[0]-i-1
            g=G*M/(r0+Z[i])**2#g0*(Rp/(Rp+Z[i]/(69911.*1E3)))**2
            H=kb*Tavg[i]/(mmw[i]*g)
            dZ[i]=H*np.log(P[i+1]/P[i])
            Z[i-1]=Z[i]+dZ[i]
            #print(P[i], H/1000., Z[i]/1000, g)

        #pdb.set_trace()
        #Interpolate values of measured cross-sections at their respective
        #temperatures pressures to the temperature and pressure of the
        #levels on which the optical depth will be computed
        t2=time.time()
        #print('Setup', t2-t1)
        #make sure   200 <T <4000 otherwise off cross section grid
        TT=np.zeros(len(Tavg))
        TT[:]=Tavg
        TT[Tavg < 500] = 500.
        TT[Tavg > 3000] = 3000.
        PP=np.zeros(len(Pavg))
        PP[:]=Pavg
        PP[Pavg < 3E-6]=3E-6
        PP[Pavg >=300 ]=300


        kcoeffs_interp=10**self.kcoeff_interp(np.log10(Pgrid), np.log10(Tgrid), np.log10(PP), np.log10(TT), wno, xsecarr)
        t3=time.time()
        #print('Kcoeff Interp', t3-t2)
        #continuum opacities (nlayers x nwnobins x ncont)***********
        xsec_cont=kcoeffs_interp[:,:,0,0]
        wave = (1/wno)*1E8
        sigmaH2 = xsec_cont*0.+1*((8.14E-13)*(wave**(-4.))*(1+(1.572E6)*(wave**(-2.))+(1.981E12)*(wave**(-4.))))*1E-4  #H2 gas Ray
        sigmaHe = xsec_cont*0.+1*((5.484E-14)*(wave**(-4.))*(1+(2.44E5)*(wave**(-2.))))*1E-4   #He gas Ray
        #Rayleigh Haze from des Etangs 2008
        wno0=1E4/0.43
        sigmaRay=xsec_cont*0.+2.E-27*amp*(wno/wno0)**power*1E-4
        #grey cloud opacity
        sigmaCld=xsec_cont*0.+CldOpac

        #mie scattering 
        xsecMie=Mies[0,0,:,:].T
        sigmaMie=np.repeat(xsecMie[np.newaxis,:,:],len(Pavg),axis=0)

        xsecContinuum=np.array([sigmaH2.T,sigmaHe.T,sigmaRay.T,sigmaCld.T]).T #building continuum xsec array (same order as cont_fracs)
        xsecContinuum=np.concatenate((xsecContinuum, sigmaMie),axis=2)
        #(add more continuum opacities here and in fractions)
        t4=time.time()
        #print("Continuum Xsec Setup ", t4-t3)
        #********************************************
        #Calculate transmissivity as a function of
        #wavenumber and height in the atmosphere
        t=self.CalcTauXsecCK(kcoeffs_interp,Z,Pavg,Tavg, Fractions, r0,gord,wts,Frac_Cont,xsecContinuum)
        t5=time.time()
        #print('Transmittance', t5-t4)    

        #Compute Integral to get (Rp/Rstar)^2 (equation in brown 2001, or tinetti 2012)
        F=((r0+np.min(Z))/(Rstar*6.95508E8))**2+2./(Rstar*6.95508E8)**2.*np.dot((1.-t),(r0+Z)*dZ)
        t6=time.time()
        #print('Total in Trans', t6-t1)

        return wno, F, Z#, TauOne


    def CalcTauXsecCK(self, kcoeffs,Z,Pavg,Tavg, Fractions, r0,gord, wts, Fractions_Continuum, xsecContinuum):
        """
        Calculate opacity using correlated-k tables 
        
        Parameters
        ----------
        kcoeffs : ndarray 
            correlated-k tables 
        Z : ndarray
            height above ref pressure
        Pavg : ndarray
            pressure grid 
        Tavg : ndarray
            temperature grid 
        Fractions : ndarray
            volume mixing ratios of species 
        Fractions_Continuum : ndarray
            volume mixing ratios of continuum species 
        xsecContinuum : ndarray
            cross sections for continuum
        Returns
        -------
        transmission
        """
        ngas=Fractions.shape[0]
        nlevels=len(Z)
        nwno=kcoeffs.shape[1]
        trans=np.zeros((nwno, nlevels))+1.
        dlarr=np.zeros((nlevels,nlevels))
        ncont=xsecContinuum.shape[-1]
        uarr=np.zeros((nlevels,nlevels))
        kb=1.38E-23
        kbTavg=kb*Tavg
        Pavg_pascal=1E5*Pavg
        for i in range(nlevels):
            for j in range(i):
                index=i-j-1
                r1=r0+Z[i]
                r2=r0+Z[i-j]
                r3=r0+Z[index]
                dlarr[i,j]=(r3**2-r1**2)**0.5-(r2**2-r1**2)**0.5
                uarr[i,j]=dlarr[i,j]*Pavg_pascal[index]/kbTavg[index]
        
        for v in range(nwno):
            for i in range(nlevels):
                transfull=1.
                #for CK gases--try to do ALL gases as CK b/c of common interpolation
                for k in range(ngas):
                    transtmp=0.
                    for l in range(len(wts)):
                        tautmp=0.
                        for j in range(i):
                            index=i-j-1
                            tautmp+=2.*Fractions[k,index]*kcoeffs[index,v,k,l]*uarr[i,j]
                        transtmp+=np.exp(-tautmp)*wts[l]/2.
                    transfull*=transtmp
                #for continuum aborbers (gas rayligh, condensate scattering etc.--nlayers x nwno x ncont
                #'''
                for k in range(ncont):
                    tautmp=0.
                    for j in range(i):

                        index=i-j-1
                        tautmp+=2.*Fractions_Continuum[k,index]*xsecContinuum[index,v,k]*uarr[i,j]

                    transfull*=np.exp(-tautmp)
                #'''
                trans[v,i]=transfull
        return trans

    def cloud_profile(self, fsed,cloud_VMR, Pavg, Pbase):
        """
        A&M2001 eq 7., but in P-coordinates (using hydrostatic) and definition of f_sed
        Parameters
        ----------
        fsed : float 
            sedimentation efficiency 
        cloud_VMR : float 
            cloud base volume mixing ratio 
        Pavg : ndarray
            pressure grid 
        Pbase : float 
            location of cloud base
        Returns
        -------
        Condensate mixing ratio as a function of pressure
        """
        cond=cloud_VMR
        loc0=np.where(Pbase >= Pavg)[0][-1]
        cond_mix=np.zeros(len(Pavg))+1E-50
        cond_mix[0:loc0+1]=cond*(Pavg[0:loc0+1]/Pavg[loc0])**fsed
        return cond_mix


    def particle_radius(self, fsed, Kzz, mmw, Tavg, Pavg, g, rho_c, mmw_c, qc, rr):
        """
        Computes particle radius based on A&M2011
        Parameters
        ----------
        fsed : float 
            sedimentation efficiency 
        Kzz : float 
            vertical mixing 
        mmw : float 
            mean molecular weight of the atmosphere
        Tavg : float
            temperature 
        Pavg : float 
            pressure 
        g : float 
            gravity 
        rho_c : float 
            density of condensate 
        mmw_c : float 
            molecular weight of condensate 
        qc : float 
            mass mixing ratio of condensate 
        rr : float 
        """
        dlnr=np.abs(np.log(rr[1])-np.log(rr[0]))
        kb=1.38E-23  #boltzman constant
        mu0=1.66E-27  #a.m.u.
        d=2.827E-10  #bath gas molecule diameter (m)
        alpha=1.4  #alpha factor from A&M 2001 (don't need to change this)
        sig_eff=2  #log-normal particle size distribution width
        
        #atmosphere properties
        H=kb*Tavg/(mmw*mu0*g)  #scale height
        rho_a=Pavg*mmw*mu0*1E5/(kb*Tavg)  #atmospheric mass density
        
        wmix=Kzz/H  #vertical mixing velocity
        mfp=kb*Tavg/(2**0.5*np.pi*d**2*Pavg*1E5)   #mean free path
        eta=5./16.*np.sqrt(np.pi*2.3*mu0*kb*Tavg)*(Tavg/59.7)**.16/(1.22*np.pi*d**2) #dynamic viscosity of bath gas
        
        #computing varius radius profiles
        r_sed=2./3.*mfp*((1.+10.125*eta*wmix*fsed/(g*(rho_c-rho_a)*mfp**2))**.5-1.)  #sedimentation radius
        r_eff=r_sed*np.exp(-0.5*(alpha+1)*np.log(sig_eff)**2)  #A&M2011 equation 17 effective radius
        r_g=r_sed*np.exp(-0.5*(alpha+6.)*np.log(sig_eff)**2) #A&M formula (13)--lognormal mean (USE THIS FOR RAD)
        
        #droplet VMR
        f_drop=3.*mmw_c*mu0*qc/(4.*np.pi*rho_c*r_g**3)*np.exp(-4.5*np.log(sig_eff)**2)  #
        prob_lnr=np.zeros((len(rr),len(r_g)))
        for i in range(len(prob_lnr)): prob_lnr[i,:]=1./((2.*np.pi)**0.5*np.log(sig_eff))*np.exp(-0.5*np.log(rr[i]/r_g)**2/np.log(sig_eff)**2)*dlnr
        f_r=prob_lnr*f_drop
        
        return r_sed, r_eff, r_g, f_r

    def profile_tempeture_and_pressure(self, kv1, kv2, kth, alpha=0.5):
        """Guillot 2010 PT profile parameterization
        Returns
        -------
        temperature , pressure
        """
        Teff = self.cross_sections.planetary_parameters['Tint']
        f = 1.0  # solar re-radiation factor
        A = 0.0  # planetary albedo
        g0 = self.atmosphere_grid['g0']
        
        # Compute equilibrium temperature and set up gamma's
        T0 = self.cross_sections.planetary_parameters['Tirr']
        gamma1 = kv1/kth
        gamma2 = kv2/kth
        
        # Initialize arrays
        logtau =np.arange(-10,20,.1)
        tau =10**logtau
        
        #computing temperature
        T4ir = 0.75*(Teff**(4.))*(tau+(2.0/3.0))
        f1 = 2.0/3.0 + 2.0/(3.0*gamma1)*(1.+(gamma1*tau/2.0-1.0)*sp.exp(-gamma1*tau))+2.0*gamma1/3.0*(1.0-tau**2.0/2.0)*sp.special.expn(2.0,gamma1*tau)
        f2 = 2.0/3.0 + 2.0/(3.0*gamma2)*(1.+(gamma2*tau/2.0-1.0)*sp.exp(-gamma2*tau))+2.0*gamma2/3.0*(1.0-tau**2.0/2.0)*sp.special.expn(2.0,gamma2*tau)
        T4v1=f*0.75*T0**4.0*(1.0-alpha)*f1
        T4v2=f*0.75*T0**4.0*alpha*f2
        T=(T4ir+T4v1+T4v2)**(0.25)
        P=tau*g0/(kth*0.1)/1.E5
        
        
        # Return TP profile
        return T, P

class LoadCrossSections(PlanetarySystem):
    """Get correlated-K opacities, stellar spectrum, and chemistry 
    Routine that loads in the correlated-K opacities
    applicable to JWST. Here we assume the data will be binned 
    to an R=100 and span 50 - 28000 cm-1 (200 - 0.3 um)
    Also loads in Mie scattering properties, stellar spectrum
    for emission, and grid chemistry.  Note, to change Mie condensate, 
    change the name in this routine.  Have a look in ../ABSCOEFF_CK/MIES/
    for a current list of condensates.  Feel free to also
    load in a different stellar spectrum (for emission). 
    Does not matter where you get it from just so long as you can
    save it as a .h5 file in wavelength (from low to high)
    and flux, both in MKS units (m, W/m2/m)
    Parameters
    ----------
    wnomin : float 
        minimum wavenumber for computation (min 50 for JWST, 2000 for HST)
    wnomax : float 
        maximum wavenumber for computation (max 28000 for JWST, 30000 for HST)
    observatory : str 
        Options are 'JWST' or 'HST'
    directory : str
        Directory path to zip file from 
        https://www.dropbox.com/sh/o4p3f8ukpfl0wg6/AADBeGuOfFLo38MGWZ8oFDX2a?dl=0&lst=
    stellar_file : str, optional 
        Stellar file for emission spectra
    cond_name : str , optional
        Name of condensate to compute cloud, if using Ackerman Marley cloud scheme
    gauss_pts : str , optional
        Number of gauss points for ck opacity tables. Default is to use 
        10 for HST and 20 for JWST. 
        Options are : ['default', '10', '20']
    Returns
    -------
    P
        pressure grid for C-K-coefficient relaevant properties
    T
        temperature grid for C-K-coefficient relaevant properties (not the same as chemistry)
    wno
        wavenumber grid of xsecs (here, R=100)
    g
        CK gauss quad g-ordinate
    wts
        CK gauss quad weights
    xsecarr 
        Pre-computed grid of CK coefficients as a function of
        Gases x Pressure x Temperature x Wavenumber x GaussQuad Points. 
    radius
        Mie scattering/condensate relevant properties. 
        Condensate particle sizes in micron (0.01 - 316 um)
    mies_arr
        mie properties array.  It contains the total extinction
        cross-section (first element--Qext*pi*r^2), single scatter albedo (second,
        Qs/Qext),and assymetry parameter (third) as a function of condensate (in
        this case, just one), particle radius (0.01 - 316 um), and wavenumber. These
        were generated offline with the "pymiecoated" routine using
        the indicies of refraction given in Wakeford & Sing 2017.
        (Condensates x MieProperties x size bins x wavenumber)
    Fstar 
        Stellar Flux spectrum binned to cross-section wavenumber grid
    logCtoO
        chemistry C/O grid in log10 (-2.0 - +0.3, -0.26 = solar)
    logMet
        chemistry metallicity grid in log10 (-2 - +3, 0=solar)
    Tarr
        Chemistry Temperature grid (400 - 3400K)
    logParr
        Chemistry Pressure grid (-7.0 - 2.4, log10 in bar)
    gases
        log of Chemistry grid of molecular gas abundances (and 
        mean molecular weight) as a function of Metallicity, C/O, T, and P.
        (CtoO x logMet x Tarr x Gases x logParr). Gases: H2O  CH4  CO  CO2 NH3  N2  
        HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He   e- h-  mmw. 
        This grid was generated with NASA CEA2 routine and assumes pure
        equilibrium+condensate chemistry (no rainout here).
    """
    ### Read in CK arrays (can switch between 10 & 20 CK gp's )

    def __init__(self, data_directory, planetary_system_file=None):
        """Initialize the class object."""

        PlanetarySystem.__init__(self, json_file=planetary_system_file)

        self.directory = data_directory
        self.cond_name = 'MgSiO3'
        self.gauss_pts = 'default'

        available_opacities = ['H2H2', 'H2He', 'H2O','CH4','CO','CO2','NH3',
                            'Na_allard', 'K_allard','TiO','VO','C2H2','HCN',
                            'H2S','FeH','HMFF','HMBF']

        #set up filename tag 
        if self.observatory == 'HST':
            if self.gauss_pts == 'default': 
                gp = '10'
            elif self.gauss_pts in ['10','20']:
                gp = self.gauss_pts
            else: 
                raise Exception('Invalid gauss point choice.')

            filename = '_CK_STIS_WFC3_'+gp+'gp_2000_30000wno.h5'

        elif self.observatory =='JWST':
            if self.gauss_pts == 'default': 
                gp = '20'
            elif self.gauss_pts in ['10','20']:
                gp = self.gauss_pts
            else: 
                raise Exception('Invalid gauss point choice.')
            filename = '_CK_R100_'+gp+'gp_50_30000wno.h5'
        else: 
            raise Exception ('Pick a valid observatory: HST or JWST')

        #set up cross section array list 
        xsecarr = []
        for mol in available_opacities: 
            file = os.path.join(self.directory,mol+filename)

            #open up h5 file
            hf=h5py.File(file, 'r')

            #get parameters that are uniform for everything
            if mol == 'H2H2':
                wno=np.array(hf['wno'])
                self.T=np.array(hf['T'])
                self.P=np.array(hf['P'])
                self.g=np.array(hf['g'])
                self.wts=np.array(hf['wts'])

            #get kcoefficients for each molecule
            kcoeff=np.array(hf['kcoeff'])
            xsecarr_mol=10**(kcoeff-4.)
            hf.close()

            #add it to master list 
            xsecarr += [xsecarr_mol]

        #super big array stuffing all the gases in
        xsecarr = np.log10(np.array(xsecarr))



        #stellar flux file------------------------------------
        if self.stellar_parameters['stellar_file'] != None:
            hf=h5py.File(self.stellar_parameters['stellar_file'], 'r')  #user should generate their own stellar spectrum from whatever database they choose, save as .h5 
            lambdastar=np.array(hf['lambdastar'])  #in MKS (meters)
            Fstar0=np.array(hf['Fstar0'])  #in MKS (W/m2/m)
            hf.close()

            lambdastar=lambdastar*1E6
            loc=np.where((lambdastar >= 1E4/wno[-1]) & (lambdastar <=1E4/wno[0]))
            lambdastar=lambdastar[loc]
            lambdastar_hi=np.arange(lambdastar.min(),lambdastar.max(),0.0001)
            Fstar0=Fstar0[loc]
            Fstar0=sp.interp(np.log10(lambdastar_hi), np.log10(lambdastar), Fstar0) 

            #smooth stellar spectrum to CK bins
            szmod=len(wno)
            Fstar_smooth=np.zeros(szmod)
            dwno=wno[1:]-wno[:-1]
            for i in range(szmod-1):
                i=i+1
                loc=np.where((1E4/lambdastar_hi >= wno[i]-0.5*dwno[i-1]) & (1E4/lambdastar_hi < wno[i]+0.5*dwno[i-1]))
                Fstar_smooth[i]=np.mean(Fstar0[loc])

            Fstar_smooth[0]=Fstar_smooth[1]
            Fstar_smooth[-1]=Fstar_smooth[-2]
            Fstar=Fstar_smooth
        else: 
            Fstar = [np.nan]

        #loading mie coefficients-----------------------------
        file=os.path.join(self.directory, 'MIE_COEFFS',self.cond_name+'_r_0.01_300um_wl_0.3_200um_interp_STIS_WFC3_2000_30000wno.h5')
        hf=h5py.File(file, 'r')
        # wno_M=np.array(hf['wno_M'])
        radius=np.array(hf['radius'])
        Mies=np.array(hf['Mies'])
        hf.close()

        SSA=Mies[1,:,:]/Mies[0,:,:]#single scatter albedo
        Mies[1,:,:]=SSA  #single scatter albedo
        Mg2SiO4=Mies  #Mies = Qext, Qs, asym
        xxsec=Mg2SiO4[0,:,:].T*np.pi*radius**2*1E-12 #scattering cross-section
        Mg2SiO4[0,:,:]=xxsec.T
        mies_arr = np.array([Mg2SiO4])

        #cropping in wavenumber 
        loc=np.where((wno <= self.wavenumber_max) & (wno >= self.wavenumber_min))[0]
        self.wno=wno[loc]
        self.xsecarr=xsecarr[:,:,:,loc,:]
        self.mies_arr=mies_arr[:,:,:,loc]
        if self.stellar_parameters['stellar_file'] != None: 
            self.Fstar=Fstar[loc]
        else: 
            self.Fstar = np.array(Fstar*len(loc))

        #loading interpolatedable chemistry grid as a function of C/O, Metallicity, T, and P (for many gases)
        hf=h5py.File(os.path.join(self.directory,'CHEM','chem_grid.h5'), 'r')
        self.logCtoO = np.array(hf['logCtoO'])  #-2.0 - 0.3
        self.logMet = np.array(hf['logMet']) #-2.0 - 3.0
        self.Tarr = np.array(hf['Tarr'])  #400 - 3400
        self.logParr = np.array(hf['logParr'])   #-7.0 - 2.4

        self.gases = np.array(hf['gases'])  ##H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He   e- h-  mmw

        self.loggas = np.log10(self.gases)
        hf.close()

        self.radius = radius*1E-6
        print('Cross-sections Loaded')

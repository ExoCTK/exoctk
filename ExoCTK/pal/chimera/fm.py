from . import ctran
from .thermo import *

# import os
import sys
import math
import numpy as np
import scipy as sp
# from array import *
# from scipy import interpolate
from scipy import signal
from scipy import special
from scipy import interp
from scipy import ndimage
# import pdb
# from matplotlib.pyplot import *
# import datetime
from pickle import *
from numba import jit

@jit
def CalcTauXsecCK(kcoeffs,Z,Pavg,Tavg, Fractions, r0,gord, wts, Fractions_Continuum, xsecContinuum):
    """
    Calculate transmission with correlated-k method.

    Parameters
    ----------
    kcoeffs
    Z
    Pavg
    Tavg
    Fractions
    r0
    gord
    wts
    Fractions_Continuum
    xsecContinuum

    Returns
    -------

    """
    ngas=Fractions.shape[0]
    nlevels=len(Z)
    nwno=kcoeffs.shape[1]
    trans=np.zeros((nwno, nlevels))+1.
    dlarr=np.zeros((nlevels,nlevels))
    ncont=xsecContinuum.shape[-1]
    for i in range(nlevels-2):
        for j in range(i):
            r1=r0+Z[i]
            r2=r0+Z[i-j]
            r3=r0+Z[i-j-1]
            dlarr[i,j]=np.sqrt(r3**2-r1**2)-np.sqrt(r2**2-r1**2)

    kb=1.38E-23
    for v in range(nwno):
        for i in range(nlevels-2):
            transfull=1.
            #for CK gases--try to do ALL gases as CK b/c of common interpolation
            for k in range(ngas):
                transtmp=0.
                for l in range(len(wts)):
                    tautmp=0.
                    for j in range(i):
                        dl=dlarr[i,j]
                        curlevel=i-j-1
                        u=dl*1E5*Pavg[curlevel]/(kb*Tavg[curlevel])
                        tautmp+=2.*Fractions[k,curlevel]*kcoeffs[curlevel,v,k,l]*u
                    transtmp+=np.exp(-tautmp)*wts[l]/2.
                transfull*=transtmp
            #for continuum aborbers (gas rayligh, condensate scattering etc.--nlayers x nwno x ncont
            #'''
            for k in range(ncont):
                tautmp=0.
                for j in range(i):
                    dl=dlarr[i,j]
                    curlevel=i-j-1
                    u=dl*1E5*Pavg[curlevel]/(kb*Tavg[curlevel])*Fractions_Continuum[k,curlevel]
                    #print(Fractions_Continuum[k,curlevel])
                    tautmp+=2.*xsecContinuum[curlevel,v,k]*u
                transfull*=np.exp(-tautmp)
            #'''
            trans[v,i]=transfull
    return trans

def tran(T, P, mmw, Ps, Pc, alphaH2O, alphaCH4, alphaCO, alphaCO2, alphaNH3, alphaNaK, alphaTiO,
         alphaVO, alphaC2H2, alphaHCN, alphaH2S, alphaFeH, fH2, fHe, amp, power, M, Rstar, Rp, 
         wnomin, wnomax, path):
    """
    Calculate the transmission.

    Parameters
    ----------
    T
    P
    mmw
    Ps
    Pc
    alphaH2O
    alphaCH4
    alphaCO
    alphaCO2
    alphaNH3
    alphaNaK
    alphaTiO
    alphaVO
    alphaC2H2
    alphaHCN
    alphaH2S
    alphaFeH
    fH2
    fHe
    amp
    power
    M
    Rstar
    Rp
    wnomin
    wnomax
    path

    Returns
    -------

    """
    #print "Starting tran at ", datetime.datetime.now().time()

    #Convert parameters to proper units
    #converting molar mixing ratios in ppm to just mole fraction
    fH2=fH2
    fHe=fHe
    fH2O=alphaH2O
    fCH4=alphaCH4
    fCO=alphaCO
    fCO2=alphaCO2
    fNH3=alphaNH3
    fNaK=alphaNaK
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
    fNa=fNaK
    fK=fNa*0.05625
    
    Fractions = np.array([fH2*fH2,fHe*fH2,fH2O, fCH4, fCO, fCO2, fNH3,fNa,fK,fTiO,fVO,fC2H2,fHCN,fH2S,fFeH])  #gas mole fraction profiles
                        #H2Ray, HeRay  Ray General,
    Frac_Cont = np.array([fH2,fHe,fH2*0.+1.])  #continuum mole fraction profiles
    #Load measured cross-sectional values and their corresponding
    #T,P,and wno grids on which they were measured
    Pgrid, Tgrid, wno, gord, wts, xsecarr = xsects(path)
    # Pgrid = restore.xsects[0]
    # Tgrid = restore.xsects[1]
    # wno = restore.xsects[2]
    # gord=restore.xsects[3]
    # wts=restore.xsects[4]
    # xsecarr = restore.xsects[5]
    
    #Calculate Temperature, Pressure and Height grids on which
    #transmissivity will be computed
    n = len(P)
    nv = len(wno)
    
    
    Z=np.zeros(n)  #level altitudes
    dZ=np.zeros(n)  #layer thickness array
    r0=Rp*69911.*1.E3  #converting planet radius to meters
    mmw=mmw*1.660539E-27  #converting mmw to Kg
    kb=1.38E-23
    G=6.67384E-11
    M=M*1.898E27

    #cropping the wavenumber grid over selected range wnomin to wnomax
    loc=np.where((wno <= wnomax) & (wno >= wnomin))
    loc=loc[0]
    wno_offset = loc[0]
    wnocrop=wno[loc]
    count=len(wnocrop)
    
    #Compute avg Temperature at each grid
    Tavg = np.array([0.0]*(n-1))
    Pavg = np.array([0.0]*(n-1))
    for z in range(n-1):
        Pavg[z] = np.sqrt(P[z]*P[z+1])
        Tavg[z] = interp(np.log10(Pavg[z]),sp.log10(P),T)
    #create hydrostatic altitutde grid from P and T
    Phigh=P.compress((P>Ps).flat)  #deeper than reference pressure
    Plow=P.compress((P<=Ps).flat)   #shallower than reference pressure
    for i in range(Phigh.shape[0]):  #looping over levels above ref pressure
        i=i+Plow.shape[0]-1
        g=G*M/(r0+Z[i])**2#g0*(Rp/(Rp+Z[i]/(69911.*1E3)))**2
        H=kb*Tavg[i]/(mmw[i]*g)  #scale height
        dZ[i]=H*np.log(P[i+1]/P[i]) #layer thickness, dZ is negative
        Z[i+1]=Z[i]-dZ[i]   #level altitude
    for i in range(Plow.shape[0]-1):  #looping over levels below ref pressure
        i=Plow.shape[0]-i-1
        g=G*M/(r0+Z[i])**2#g0*(Rp/(Rp+Z[i]/(69911.*1E3)))**2
        H=kb*Tavg[i]/(mmw[i]*g)
        dZ[i]=H*np.log(P[i+1]/P[i])
        Z[i-1]=Z[i]+dZ[i]
    #Interpolate values of measured cross-sections at their respective
    #temperatures pressures to the temperature and pressure of the
    #levels on which the optical depth will be computed
    #print "Interpolating cross-sections at ", datetime.datetime.now().time()
    #make sure   200 <T <4000 otherwise off cross section grid
    TT=np.zeros(len(Tavg))
    TT[:]=Tavg
    TT[Tavg < 250] = 250.
    TT[Tavg > 3000] = 3000.
    PP=np.zeros(len(Pavg))
    PP[:]=Pavg
    PP[Pavg < 1E-6]=1E-6
    kcoeffs_interp=np.zeros((len(Pavg),len(wnocrop),xsecarr.shape[0],len(gord)))
    for i in range(len(gord)):
        xsecarr_inter, xsecarr_interRayleigh = ctran.InitXsects(xsecarr[:,:,:,:,i], Tgrid, Pgrid, TT, PP, wnocrop, wno_offset, amp,power )
        kcoeffs_interp[:,:,:,i]=xsecarr_inter  #nlayers x nwnobins x ngas x ngauss
    #continuum opacities (nlayers x nwnobins x ncont)***********
    xsec_cont=xsecarr_inter[:,:,0]
    wave = (1/wnocrop)*1E8
    sigmaH2 = xsec_cont*0.+1*((8.14E-13)*(wave**(-4.))*(1+(1.572E6)*(wave**(-2.))+(1.981E12)*(wave**(-4.))))*1E-4  #H2 gas Ray
    sigmaHe = xsec_cont*0.+1*((5.484E-14)*(wave**(-4.))*(1+(2.44E5)*(wave**(-2.))))*1E-4   #He gas Ray
    #Rayleigh Haze from des Etangs 2008
    wno0=1E4/0.43
    sigmaRay=xsec_cont*0.+2.E-27*amp*(wnocrop/wno0)**power*1E-4
    xsecContinuum=np.array([sigmaH2.T,sigmaHe.T,sigmaRay.T]).T #building continuum xsec array (same order as cont_fracs)
    #(add more continuum opacities here and in fractions)
    #********************************************
    #Calculate transmissivity as a function of
    #wavenumber and height in the atmosphere
    #print "Computing Transmittance ", datetime.datetime.now().time()
    t=CalcTauXsecCK(kcoeffs_interp,Z,Pavg,Tavg, Fractions, r0,gord,wts,Frac_Cont,xsecContinuum)
    locPc=np.where(P >= Pc)
    t[:,locPc]=0.
    #pdb.set_trace()
    #Compute Integral to get (Rp/Rstar)^2 (equation in brown 2001, or tinetti 2012)
    F=((r0+np.min(Z[:-1]))/(Rstar*6.955E8))**2+2./(Rstar*6.955E8)**2.*np.dot((1.-t),(r0+Z)*dZ)
    #print "Ending Tran at ", datetime.datetime.now().time()
    return wnocrop, F, Z#, TauOne


def xsects(path):
    """
    Read in cross-section arrays

    Parameters
    ----------
    path: str
        Path to directory with cross-section files

    Returns
    -------
    P: np.array
        Pressure points
    T: np.ndarray
        Temperature points
    wno: np.ndarray
        Wavenumber points
    g: np.ndarray
        g
    wts: np.ndarray
        wts
    chemarray: np.ndarray
        cross-sections
    """

    if sys.version_info.major >= 3:
        def load_pickle(file):
            # Python 3 tries to encode with ascii by default
            return load(file, encoding='bytes')
    else:
        def load_pickle(file):
            return load(file)

    CK = ['H2H2', 'H2He', 'H2O', 'CH4', 'CO', 'CO2', 'NH3', 'Na', 
          'K', 'TiO', 'VO', 'C2H2', 'HCN', 'H2S', 'FeH']

    CKarrs = []
    for species in CK:
        file = '{}CKarr{}_R100_900_16500.pic'.format(path,species)
        P, T, wno, kcoeff, g, wts = load_pickle(open(file,'rb'))
        CKarrs.append(10**(kcoeff-4.))
        
    
    return P, T, wno, g, wts, np.log10(np.array(CKarrs))


def TP(Teq, Teeff, g00, kv1, kv2, kth, alpha):
    """
    This function takes stellar, planetary, and atmospheric parameters and
    returns the temperature-pressure profile.

    Parameters
    ----------
    Teq
    Teeff
    g00
    kv1
    kv2
    kth
    alpha

    Returns
    -------
    T: np.ndarray
        The Temperature in Kelvin
    P: np.ndarray
        The Pressure in bar
    """
    
    
    Teff = Teeff
    f = 1.0  # solar re-radiation factor
    A = 0.0  # planetary albedo
    g0 = g00
    
    # Compute equilibrium temperature and set up gamma's
    T0 = Teq
    gamma1 = kv1/kth
    gamma2 = kv2/kth
    
    # Initialize arrays
    logtau =np.arange(-10,20,.1)
    tau =10**logtau
    
    
    #computing temperature
    T4ir = 0.75*(Teff**(4.))*(tau+(2.0/3.0))
    f1 = 2.0/3.0 + 2.0/(3.0*gamma1)*(1.+(gamma1*tau/2.0-1.0)*sp.exp(-gamma1*tau))+2.0*gamma1/3.0*(1.0-tau**2.0/2.0)*special.expn(2.0,gamma1*tau)
    f2 = 2.0/3.0 + 2.0/(3.0*gamma2)*(1.+(gamma2*tau/2.0-1.0)*sp.exp(-gamma2*tau))+2.0*gamma2/3.0*(1.0-tau**2.0/2.0)*special.expn(2.0,gamma2*tau)
    T4v1=f*0.75*T0**4.0*(1.0-alpha)*f1
    T4v2=f*0.75*T0**4.0*alpha*f2
    T=(T4ir+T4v1+T4v2)**(0.25)
    P=tau*g0/(kth*0.1)/1.E5
    
    
    # Return TP profile
    return T, P


def fx(x, gas_scale, path, cea_path=None):
    """
    Forward model--takes in state vector and returns the binned model points to
    compare directly to data.

    Be sure to change planet parameters when going to a new planet!

    Parameters
    ----------
    x: np.ndarray
        State vector [Tiso, logKir,logg1, logMet, logCtoO, logPQCarbon,
        logPQNitrogen, Rp, Rstar, M, RayAmp, RaySlp, logPc]
    gas_scale: np.ndarray
        Scaling factors for each gas
    path: str
        Path to cross-section files
    cea_path: str
        Path to a CEA executable to compute abundances on the fly

    Returns
    -------
    transmission: np.ndarray
        The spectrum
    wavenumber: np.ndarray
        The wavenumber bins
    chemarray: np.ndarray
        The TP-abundance profile for each chemical species.


    """
    #print "Entering Fx ", datetime.datetime.now().time()
    #  0    1        2       3     4          5           6          7     8    9       10        11           12        
    #Tiso, logKir,logg1, logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp, Rstar, M,   RayAmp    RaySlp        logPc
    
    #Guillot 2010 TP profile (3 params)
    Tiso=x[0]
    logKir=x[1]
    logg1=x[2]
    #Chemistry
    Met=10.**x[3]  #metallicity
    CtoO=10.**x[4] #C/O
    logPQC=x[5]  #carbon quench pressure
    logPQN=x[6]  #nitrogen quench pressure
    #planet params
    Rp=x[7]  #planet radius (in jupiter)
    Rstar=x[8]   #stellar radius (in solar)
    M=x[9]   #planet mass (in jupiter)
    #cloud params
    RayAmp=10**x[10]
    RaySlp=x[11]
    Pc=10.**x[12]

  
 
    
    #Setting up atmosphere grid****************************************
    logP = np.arange(-7,1.5,0.1)+0.1
    P = 10.0**logP
    g0=6.67384E-11*M*1.898E27/(Rp*69911.*1.E3)**2
    kv=10.**(logg1+logKir)
    kth=10.**logKir
    tp=TP(Tiso, 100,g0 , kv, kv, kth, 0.5)

    if cea_path is not None:
        T = interp(logP,np.log10(tp[1]),tp[0])
        #making courser chemistry grid
        Pchem=np.append(P[::4],P[-1])
        Tchem=np.append(T[::4],T[-1])
        H2Oarr, CH4arr, COarr, CO2arr, NH3arr, N2arr, HCNarr, H2Sarr, PH3arr, C2H2arr, C2H6arr, Naarr, Karr, TiOarr, VOarr, FeHarr, Harr, H2arr, Hearr, MMWarr=thermo(Met, CtoO, Tchem, Pchem,'foo_'+str(random()), cea_path)

        # Interoplating back to full grid
        H2Oarr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(H2Oarr))*gas_scale[0]
        CH4arr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(CH4arr))*gas_scale[1]
        COarr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(COarr))*gas_scale[2]
        CO2arr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(CO2arr))*gas_scale[3]
        NH3arr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(NH3arr))*gas_scale[4]
        N2arr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(N2arr))*gas_scale[5]
        HCNarr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(HCNarr))*gas_scale[6]
        H2Sarr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(H2Sarr))*gas_scale[7]
        PH3arr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(PH3arr))*gas_scale[8]
        C2H2arr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(C2H2arr))*gas_scale[9]
        C2H6arr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(C2H6arr))*gas_scale[10]
        Naarr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(Naarr))*gas_scale[11]
        Karr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(Karr))*gas_scale[12]
        TiOarr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(TiOarr))*gas_scale[13]
        VOarr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(VOarr))*gas_scale[14]
        FeHarr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(FeHarr))*gas_scale[15]
        Harr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(Harr))*gas_scale[16]
        H2arr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(H2arr))*gas_scale[17]
        Hearr = 10.**interp(np.log10(P),np.log10(Pchem),np.log10(Hearr))*gas_scale[18]
        MMWarr = interp(np.log10(P),np.log10(Pchem), MMWarr)

    
    #poor mans rain-out
    #TiO/VO rainout hack
    if np.min(TiOarr) < 1E-12:
        loc=np.where(TiOarr < 1E-12)
        TiOarr[P < P[loc[0][-1]]]=1E-14
    if np.min(VOarr) < 1E-12:
        loc=np.where(VOarr < 1E-12)
        VOarr[P < P[loc[0][-1]]]=1E-14

    #condensing out Na/K
    #Na/K rainout hack
    if np.min(Naarr) < 1E-12:
        loc=np.where(Naarr < 1E-12)
        Naarr[P < P[loc[0][-1]]]=1E-14
    if np.min(Karr) < 1E-12:
        loc=np.where(Karr < 1E-12)
        Karr[P < P[loc[0][-1]]]=1E-14

    #condensing out FeH
    if np.min(FeHarr) < 1E-12:
        loc=np.where(FeHarr < 1E-12)
        FeHarr[P < P[loc[0][-1]]]=1E-14

    PQC=10.**logPQC
    loc=np.where(P <= PQC)
    CH4arr[loc]=CH4arr[loc][-1]
    COarr[loc]=COarr[loc][-1]
    H2Oarr[loc]=H2Oarr[loc][-1]
    CO2arr[loc]=CO2arr[loc][-1]


    PQN=10.**logPQN
    loc=np.where(P <= PQN)
    NH3arr[loc]=NH3arr[loc][-1]
    N2arr[loc]=N2arr[loc][-1]

    mmw = MMWarr

    #wavenumber range over which to compute spectrum (min 500, max 30000)
    wnomin =910#6000
    wnomax =15500. #9800.


    Pref=10.1  #reference pressure bar-keep fixed
    #cloud profile
    fcarr=np.zeros(len(H2arr))
    Qc=0.
    Rc=0.
    #computing transmission spectrum
    spec = tran(T, P, mmw, Pref, Pc, H2Oarr, CH4arr, COarr, CO2arr, NH3arr, Naarr+Karr, TiOarr, 
                VOarr, C2H2arr, HCNarr, H2Sarr, FeHarr, H2arr, Hearr, RayAmp, RaySlp, M, Rstar, 
                Rp, wnomin, wnomax, path)
    wnocrop = spec[0]
    F = spec[1]
    #print "Exiting Fx ", datetime.datetime.now().time()
    chemarr=np.array([P,T, H2Oarr, CH4arr, COarr,CO2arr, NH3arr, N2arr, H2Sarr, HCNarr, C2H2arr, C2H6arr, H2arr, Naarr, Karr])

    return F, wnocrop, chemarr



"""Utility and numba functions for chimera.py
"""

from audioop import cross
from copy import deepcopy
import numpy as np
import scipy as sp
from numba import jit


@jit(nopython=True)
def cloud_profile(fsed,cloud_VMR, Pavg, Pbase):
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
    cond_mix[0:loc0+1]=cond*(Pavg[0:loc0+1]/Pavg[loc0])**fsed  #
    return cond_mix


@jit(nopython=True)
def particle_radius(fsed,Kzz,mmw,Tavg, Pavg,g, rho_c,mmw_c, qc,rr):
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


@jit(nopython=True)
def CalcTauXsecCK(kcoeffs,Z,Pavg,Tavg, Fractions, r0, gord, wts, Fractions_Continuum, xsecContinuum):
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


@jit(nopython=True)
def kcoeff_interp(logPgrid, logTgrid, logPatm, logTatm, wnogrid, kcoeff):
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

    Ng, NP, NT, Nwno, Nord = kcoeff.shape

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

@jit(nopython=True)
def mix_two_gas_CK(k1,k2,VMR1,VMR2,gord, wts):
    """
    Key function that properly mixes the CK coefficients
    for two individual gases via the "resort-rebin" procedure
    as descrbied in Lacis & Oinas 1991, Molliere et al. 2015 and 
    Amundsen et al. 2017. Each pair of gases can be treated as a
    new "hybrid" gas that can then be mixed again with another
    gas. That is the resort-rebin magic.  This is all for a *single*
    wavenumber bin for a single pair of gases.
    Parameters
    ---------- 
    k1 : ndarray 
        k-coeffs for gas 1 (on Nk ordinates)
    k2 : ndarray 
        k-coeffs for gas 2
    VMR1 : ndarray 
        volume mixing ratio of gas 1
    VMR2 : ndarray 
        volume mixing ratio for gas 2
    gord : ndarray
        g-ordinate array for gauss. quad.
    wts : ndarray 
        gauss quadrature wts--same for both gases
    Returns
    -------
    kmix_bin
        mixed CK coefficients for the given pair of gases
    VMR
        volume mixing ratio of "mixed gas".
    """
    VMR=VMR1+VMR2  #"new" VMR is sum of individual VMR's
    Nk=len(wts)
    kmix=np.zeros(Nk**2)   #Nk^2 mixed k-coeff array
    wtsmix=np.zeros(Nk**2) #Nk^2 mixed weights array
    #mixing two gases weighting by their relative VMR
    for i in range(Nk):
        for j in range(Nk):
            kmix[i*Nk+j]=(VMR1*k1[i]+VMR2*k2[j])/VMR #equation 9 Amundsen 2017 (equation 20 Mollier 2015)
            wtsmix[i*Nk+j]=wts[i]*wts[j]    #equation 10 Amundsen 2017

    #resort-rebin procedure--see Amundsen et al. 2016 or section B.2.1 in Molliere et al. 2015
    sort_indicies=np.argsort(kmix)  #sort new "mixed" k-coeff's from low to high--these are indicies
    kmix_sort=kmix[sort_indicies]  #sort k-coeffs from low to high
    wtsmix_sort=wtsmix[sort_indicies]  #sort mixed weights using same indicie mapping from sorted mixed k-coeffs
    #combining w/weights--see description on Molliere et al. 2015--not sure why this works..similar to Amundson 2016 weighted avg?
    int=np.cumsum(wtsmix_sort)
    x=int/np.max(int)*2.-1
    logkmix=np.log10(kmix_sort)
    #kmix_bin=10**np.interp(gord,x,logkmix)  #interpolating via cumulative sum of sorted weights...
    kmix_bin=np.zeros(len(gord))
    for i in range(len(gord)):
        loc=np.where(x >= gord[i])[0][0]
        kmix_bin[i]=10**logkmix[loc]

    return kmix_bin, VMR


@jit(nopython=True)
def mix_multi_gas_CK(CK,VMR,gord, wts):
    """
    Key function that properly mixes the CK coefficients
    for multiple gases by treating a pair of gases at a time.
    Each pair becomes a "hybrid" gas that can be mixed in a pair
    with another gas, succesively. This is performed at a given
    wavenumber and atmospheric layer.
    
    Parameters
    ----------
    CK : ndarray 
        array of CK-coeffs for each gas: Ngas x nordinates at a given wavenumber and pressure level
    VMR : ndarray 
        array of mixing ratios for Ngas.
    gord : ndarray 
        g-ordinates
    wts : ndarray 
        gauss quadrature wts--same for both gases
    
    Returns
    -------
    kmix_bin 
        mixed CK coefficients for the given pair of gases
    VMR 
        Volume mixing ratio of "mixed gas".
    """
    ngas=CK.shape[0]
    #begin by mixing first two gases
    kmix,VMRmix=mix_two_gas_CK(CK[0,:],CK[1,:],VMR[0],VMR[1],gord,wts)
    loc1=np.where((VMR > 1E-12) & (CK[:,-1] > 1E-50))[0]
    ngas=len(loc1)
    #mixing in rest of gases inside a loop
    for j in range(2,ngas):
        kmix,VMRmix=mix_two_gas_CK(kmix,CK[loc1[j],:],VMRmix,VMR[loc1[j]],gord,wts)
    #kmix,VMRmix=mix_two_gas_CK(kmix,CK[j,:],VMRmix,VMR[j],gord,wts)
    
    
    return kmix, VMRmix


@jit(nopython=True)
def compute_tau(CK,xsecContinuum,Mies, mass_path, Fractions,Fractions_Continuum, Fractions_Cond, gord, wts):  #this is the bottleneck right here...
    """
    Key function that computes the layer optical depths
    at each wavenumber,and g-ordiante. It also does the confusing mixing
    for the single-scatter abledo and asymetry parameter by
    appropriately weighting each by the scattering/extincition
    optical depths.  Each g-bin is treated like a psuedo"wavenumber" bin.
    Check out: https://spacescience.arc.nasa.gov/mars-climate-modeling-group/brief.html
    Parameters
    ----------
    CK : ndarray 
        array of interpolated-to-atmosphere grid CK-coeffs for each gas (Nlayers x Nwavenumbers x Ngas x Ngordinates)
    xsecContinuum : ndarray 
        CK coefficients for continuum gases, here just the rayleigh scattering opacities. Each g-bin is 'flat'
    Mies : ndarray 
        Condensate Mie scattering extenction cross sections (e.g., Qe*pi*r^2) for each condensate, particle size, wavenumber
    wts : ndarray 
        gauss quadrature wts--same for both gases
    Returns
    -------
    kmix_bin 
        mixed CK coefficients for the given pair of gases
    VMR
        Volume mixing ratio of "mixed gas".
    """
    Nlay=CK.shape[0]
    Nord=CK.shape[2]
    dtau_gas=np.zeros((Nlay,Nord))
    dtau_cont=np.zeros((Nlay,Nord))
    dtau_cond=np.zeros((Nlay,Nord))
    ssa=np.zeros((Nlay,Nord))
    asym=np.zeros((Nlay,Nord))
    
    for j in range(Nlay):
        k, VMR=mix_multi_gas_CK(CK[j,:,:],Fractions[:,j],gord,wts)
        dtau_gas[j,:]=VMR*k*mass_path[j]
        
        #add continuum opacities here--just add to dtau linearly
        weighted_xsec=np.sum(xsecContinuum[j,:]*Fractions_Continuum[:,j])  #summing xsecs x abundances
        xsec_cont=np.zeros(Nord)+weighted_xsec  #crating array that is size of n-ordiantes to fill with continuum xsecs (flat k-dist)
        dtau_cont[j,:]=xsec_cont*mass_path[j]   #computing continuum optical depth
        
        #everything condensate scattering here!!
        #total extinction cross-section of condensates (in g-space)
        weighted_cond_xsec=np.sum(Fractions_Cond[:,j]*Mies[0,0,:])
        xsec_cond=np.zeros(Nord)+weighted_cond_xsec  #crating array that is size of n-ordiantes to fill with continuum xsecs (flat k-dist)
        dtau_cond[j,:]=xsec_cond*mass_path[j]
        
        #weighted ssa and asym--this is hard--espceially after 2 beers....it's a weird "weighted in g-space ssa"
        #ssa
        weighted_ssa=np.sum(Fractions_Cond[:,j]*Mies[0,1,:]*Mies[0,0,:])  #what did I do here???
        ssa_cond=np.zeros(Nord)+weighted_ssa
        ssa[j,:]=(dtau_cont[j,:]*0.999999+ssa_cond*mass_path[j])/(dtau_cont[j,:]+dtau_cond[j,:]+dtau_gas[j,:])
        #(gotta keep that ray scattering under control..thus the 0.99999)
        #asym, aka <cos(theta)>
        weighted_asym=np.sum(Fractions_Cond[:,j]*Mies[0,0,:]*Mies[0,2,:]*Mies[0,1,:])
        asym_cond=np.zeros(Nord)+weighted_asym
        asym[j,:]=(asym_cond*mass_path[j])/(ssa_cond*mass_path[j]+dtau_cont[j,:])
    
    #print '--------------'
    #print ssa[:,0]
    
    dtau=dtau_cont+dtau_gas+dtau_cond*1.
    return dtau,dtau_cond,dtau_gas, ssa*1., asym*1.

@jit(nopython=True)
def blackbody(T,wl):
    """
    This function takes in a temperature (T) and a
    wavelength grid (wl) and returns a blackbody flux grid.
    This is used to compute the layer/slab "emission". All in MKS units.
    Parameters
    ---------- 
    T : float 
        temperature of blackbody in kelvin
    wl : ndarray 
        wavelength grid in meters
    Returns
    -------
    B : ndarray 
        an array of blackbody fluxes (W/m2/m/ster) at each wavelength (size Nwavelengths)
    """
    # Define constants used in calculation
    h = 6.626E-34
    c = 3.0E8
    k = 1.38E-23
    
    # Calculate Blackbody Flux (B) at each wavelength point (wl)
    B = ((2.0*h*c**2.0)/(wl**5.0))*(1.0/(sp.exp((h*c)/(wl*k*T)) - 1.0))
    
    # Return blackbody flux
    return B


@jit(nopython = True)
def tri_diag_solve(l, a, b, c, d):
    '''
    Tri-diagnoal matrix inversion solver. This is used
    for the two-stream radiative transver matrix inversions
    to solve for the boundary-condition coefficents on the layer
    interfaces.  
    A, B, C and D refer to: A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = D(I)
    Parameters
    ---------- 
    L : int 
        array size
    A : array or list
    B : array or list
    C : array or list
    C : array or list
    Returns
    -------
    solution coefficients, XK
    '''   
    AS, DS, CS, DS,XK = np.zeros(l), np.zeros(l), np.zeros(l), np.zeros(l), np.zeros(l) # copy arrays
    
    AS[-1] = a[-1]/b[-1]
    DS[-1] = d[-1]/b[-1]
    
    for i in range(l-2, -1, -1):
        x = 1.0 / (b[i] - c[i] * AS[i+1])
        AS[i] = a[i] * x
        DS[i] = (d[i]-c[i] * DS[i+1]) * x
    XK[0] = DS[0]
    for i in range(1,l):
        XK[i] = DS[i] - AS[i] * XK[i-1]
    return XK


@jit(nopython=True)
def src_func_loop(B,tautop,Bsurf,ubari,nlay,lam,dtau,taump,taut,ssa,hg,k1,k2,B0,B1,uarr,w):
    '''
    Two stream source function technique described in 
    Toon, McKay, & Ackerman, 1989, JGR, 94.
    Parameters
    ----------
    B : ndarray 
        planck function array (with optical depth)
    tautop : ndarray 
        optical depth of top most layer
    Bsurf : ndarray 
        "surface" blackbody
    ubari : float 
        0.5 for hemispheric mean approximation
    nlay : int 
        number of model layers
    lam : narray
        Wavelength 
    dtau : ndarray
        optical depth per layer
    taump : ndarray 
    taut : ndarray 
    ssa : ndarray 
        single scattering albedo 
    hg : ndarray
        asymmetry 
    k1 : ndarray 
    k2 : ndarray 
    B0 : ndarray 
    B1 : ndarray 
    uarr : ndarray 
    w : ndarray
    Returns
    -------
    Fu 
        layer interface monochromatic upward flux
    Fd
        layer interface monochromatic downward flux
    '''
    twopi=2.*np.pi

    Fd=np.zeros(nlay+1)
    Fu=np.zeros(nlay+1)
 
    alphax=((1.-ssa)/(1.-ssa*hg))**0.5
    g=twopi*ssa*k1*(1+hg*alphax)/(1.+alphax)
    h=twopi*ssa*k2*(1-hg*alphax)/(1.+alphax)
    xj=twopi*ssa*k1*(1-hg*alphax)/(1.+alphax)
    xk=twopi*ssa*k2*(1+hg*alphax)/(1.+alphax)
    alpha1=twopi*(B0+B1*(ubari*ssa*hg/(1.-ssa*hg)))
    alpha2=twopi*B1
    sigma1=twopi*(B0-B1*(ubari*ssa*hg/(1.-ssa*hg)))
    sigma2=alpha2
    
    #so overflow/underflows don't happen
    g[ssa < 0.01]=0.
    h[ssa < 0.01]=0.
    xj[ssa < 0.01]=0.
    xk[ssa < 0.01]=0.
    alpha1[ssa < 0.01]=twopi*B0[ssa < 0.01]
    alpha2[ssa < 0.01]=twopi*B1[ssa < 0.01]
    sigma1[ssa < 0.01]=alpha1[ssa < 0.01]
    sigma2[ssa < 0.01]=alpha2[ssa < 0.01]
    
    #more array definitions
    fpt=np.zeros(nlay+1)
    fmt=np.zeros(nlay+1)
   
    em=np.exp(-lam*dtau)
    obj=lam*dtau
    obj[obj > 35.]=35.
    epp=np.exp(obj)
    em4_mp=np.exp(-lam*dtau)
    obj2=0.5*lam*dtau
    obj2[obj2 > 35.]=35.
    epp2=np.exp(obj2)
    ngauss=len(uarr)
    
    for i in range(ngauss):
        ugauss=uarr[i]
        fpt[:]=0.
        fmt[:]=0.
       
        fpt[-1]=twopi*(Bsurf+B1[-1]*ugauss)  #bottom BD
        fmt[0]=twopi*(1.-np.exp(-tautop/ugauss))*B[0]  #top BC
        em2=np.exp(-dtau/ugauss)
        em3=em*em2
        #ray tracing intensities from bottom to top (upwards intensity, fpt) and from top to bottom (downwards intensity, fmt)     
        for j in range(nlay):  #j is from TOA "down"
        	#downards emission intensity from TOA
            fmt[j+1]=fmt[j]*em2[j]+(xj[j]/(lam[j]*ugauss+1.))*(epp[j]-em2[j])+(xk[j]/(lam[j]*ugauss-1.))*(em2[j]-em[j])+sigma1[j]*(1.-em2[j])+sigma2[j]*(ugauss*em2[j]+dtau[j]-ugauss)
            
            #upwards emission intensity from bottom of atmosphere (nlay-1)
            z=nlay-1-j  #flipping indicies to integrate from bottom up
            fpt[z]=fpt[z+1]*em2[z]+(g[z]/(lam[z]*ugauss-1))*(epp[z]*em2[z]-1)+(h[z]/(lam[z]*ugauss+1))*(1.-em3[z])+alpha1[z]*(1.-em2[z])+alpha2[z]*(ugauss-(dtau[z]+ugauss)*em2[z])
        
        #gauss quadrature integration for fluxes
        Fu=Fu+fpt*uarr[i]*w[i]
        Fd=Fd+fmt*uarr[i]*w[i]
    
    return Fu, Fd



@jit(nopython=True)
def toon(dtau1, ssa1, hg1, B):
    """
    Monochromatic two-stream radiative transfer solver described in
    Toon, McKay, & Ackerman, 1989, JGR, 94.  Modified from
    https://github.com/adamkovics/atmosphere/blob/master/atmosphere/rt/twostream.py
    (described in Adamkovics et al. 2016) and the Ames Mars GCM radiative transfer 
    from https://spacescience.arc.nasa.gov/mars-climate-modeling-group/models.html,
    Hollingsworth et al. This then calls the src_func_loop which recomputes the l
    ayer intensities at select cos(ZA) using the two stream solution as the intial
    source function.
    Parameters
    ----------
    dtau : ndarray
        layer/slab optical depths
    ssa : ndarray 
        layer/slab single scatter albedo
    hg  : ndarray 
        layer/slab asymmetry parameter
    B : ndarray 
        planck function at each *level* (N levels, N-1 layers/slabs)
    Returns
    -------
    Fup : ndarray
        monochromatic upwards flux at each layer interface
    Fdown : ndarray
        monochromatic downwards flux at each layer interface
    """
    dtau1[dtau1 < 1E-5]=1E-5

    #delta eddington correction for peaky scatterers
    ssa=(1.-hg1**2)*ssa1/(1.-ssa1*hg1**2)
    dtau=(1.-ssa1*hg1**2)*dtau1
    hg=hg1/(1.+hg1)

    ubari=0.5#1./np.sqrt(3)#0.5
    nlay = len(dtau)
    
    taub = np.cumsum(dtau)    # Cumulative optical depth at layer bottoms
    taut=np.zeros(len(dtau))      
    taut[1:]=taub[:-1] # Cumulative optical depth at layer tops (right? the bottom of one is the top of the one below it =)
    
    taump=taut+0.5*dtau #midpoint optical depths
    
    twopi = np.pi+np.pi  #2pi
    
    
    #AMES MARS CODE equations--Hemispheric Mean Approximation for plankian source (ubari=0.5 in IR)
    #see also Table 1 in Toon et al. 1989, plus some algebra
    alpha = ((1.-ssa)/(1.-ssa*hg))**0.5
    lam = alpha*(1.-ssa*hg)/ubari
    gamma = (1.-alpha)/(1.+alpha)
    term = ubari/(1.-ssa*hg)
    
    #computing linearized planck function (the glorious linear-in-tau)
    B0=B[0:-1]
    B1=(B[1:]-B[:-1])/dtau
    loc=np.where(dtau <= 3E-6)[0]
    B1[loc]=0.
    B0[loc][:-1]=0.5*(B0[loc][1:]+B0[loc][:-1])
    
    # Cpm1 and Cmm1 are the C+ and C- terms evaluated at the top of the layer (at dtau=0).
    Cpm1 =B0+B1*term  #ames code
    Cmm1 =B0-B1*term
    
    # Cp and Cm are the C+ and C- terms evaluated at the bottom of the layer.
    Cp =B0+B1*dtau+B1*term #ames code
    Cm =B0+B1*dtau-B1*term
    
    #
    tautop=dtau[0]*np.exp(-1)
    Btop=(1.-np.exp(-tautop/ubari))*B[0]
    Bsurf=B[-1]
    bottom=Bsurf+B1[-1]*ubari
    
    # Solve for the coefficients of system of equations using boundary conditions
    exptrm = lam*dtau
    exptrm[exptrm>35] = 35 # clipped so that exponential doesn't explode
    Ep = np.exp(exptrm)
    Em = 1./Ep
    
    E1 = Ep + gamma*Em
    E2 = Ep - gamma*Em
    E3 = gamma*Ep + Em
    E4 = gamma*Ep - Em
    
    L = nlay+nlay
    Af = np.zeros(L)
    Bf = np.zeros(L)
    Cf = np.zeros(L)
    Df = np.zeros(L)
    
    # First Term
    Af[0] = 0.0
    Bf[0] = gamma[0] + 1.
    Cf[0] = gamma[0] - 1.
    Df[0] = Btop - Cmm1[0]
    
    AA = (E1[:-1]+E3[:-1])*(gamma[1:]-1)
    BB = (E2[:-1]+E4[:-1])*(gamma[1:]-1)
    CC = 2.*(1.-gamma[1:]*gamma[1:])
    DD = (gamma[1:]-1) * (Cpm1[1:] - Cp[:-1]) + (1-gamma[1:]) * (Cm[:-1]-Cmm1[1:])
    Af[1:-1:2]=AA
    Bf[1:-1:2]=BB
    Cf[1:-1:2]=CC
    Df[1:-1:2]=DD
    
    AA = 2.*(1.-gamma[:-1]*gamma[:-1])
    BB = (E1[:-1]-E3[:-1])*(gamma[1:]+1.)
    CC = (E1[:-1]+E3[:-1])*(gamma[1:]-1.)
    DD = E3[:-1]*(Cpm1[1:] - Cp[:-1]) + E1[:-1]*(Cm[:-1] - Cmm1[1:])
    Af[2::2]=AA
    Bf[2::2]=BB
    Cf[2::2]=CC
    Df[2::2]=DD
    
    # Last term:
    rsf=0
    Af[-1] = E1[-1]-rsf*E3[-1]
    Bf[-1] = E2[-1]-rsf*E4[-1]
    Cf[-1] = 0.0
    Df[-1] = Bsurf - Cp[-1]+rsf*Cm[-1]
    
    k=tri_diag_solve(L, Af, Bf, Cf, Df)
    
    # Unmix coefficients
    even = np.arange(0,2*nlay,2)
    odd  = even+1
    k1 = k[even] + k[odd]
    k2 = k[even] - k[odd]
    
    #this would be the raw two stream solution Fluxes, but we
    #won't use these. Will use source function technique
    Fupraw = np.pi*(k1*Ep + gamma*k2*Em + Cpm1) #
    Fdownraw=np.pi*(k1*Ep*gamma + k2*Em + Cmm1)
    
    #source function part
    uarr=np.array([0.1834346,0.5255324,0.7966665,0.9602899])  #angluar gauss quad cos zenith angles
    w=np.array([0.3626838,0.3137066, 0.2223810, 0.1012885 ])  #gauss quad weights
   
    Fup, Fdown=src_func_loop(B,tautop,Bsurf,ubari,nlay,lam,dtau,taump,taut, ssa,hg,k1,k2,B0,B1,uarr,w)
    
    return Fup,Fdown#, Fupraw, Fdownraw


@jit(nopython=True)
def toon_solar(dtau1, ssa1, hg1, rsurf, MU0, F0PI, BTOP):
    """
    Monochromatic two-stream radiative transfer solver described in
    Toon, McKay, & Ackerman, 1989, JGR, 94.  Modified from
    https://github.com/adamkovics/atmosphere/blob/master/atmosphere/rt/twostream.py
    (described in Adamkovics et al. 2016) and the Ames Mars GCM radiative transfer 
    from https://spacescience.arc.nasa.gov/mars-climate-modeling-group/models.html,
    Hollingsworth et al. This then calls the src_func_loop which recomputes the l
    ayer intensities at select cos(ZA) using the two stream solution as the intial
    source function.
    Parameters
    ----------
    dtau : ndarray 
        layer/slab optical depths
    ssa : ndarray 
        layer/slab single scatter albedo
    hg  : ndarray 
        layer/slab asymmetry parameter
    rsurf : ndarray 
        surface reflectivity
    MU0 : ndarray 
        cos(theta) where theta is both the incidence an emission angle.
    F0PI : ndarray 
        top of atmosphere incident flux * pi (the direct beam)
    BTOP : ndarray 
        top of the atmosphere diffuse flux.
    Returns
    -------
    Fup 
        monochromatic upwards reflected stellar flux at each layer interface
    Fdn
        monochromatic downwards direct and diffuse stellar flux at each layer interface
    """
    #delta eddington scaling
    dtau1[dtau1 < 1E-5]=1E-5
    ssa=(1.-hg1**2)*ssa1/(1.-ssa1*hg1**2)
    dtau=(1.-ssa1*hg1**2)*dtau1
    hg=hg1/(1.+hg1)
    
    
    nlay = len(dtau)
    
    # Cumulative optical depth
    taub = np.cumsum(dtau)
    taut=np.zeros(len(dtau))
    taut[1:]=taub[:-1]
    
    taump=taut+0.5*dtau
    
    # Surface reflectance and lower boundary condition
    bsurf = rsurf * MU0 * F0PI * np.exp(-taub[-1]/MU0)
    
    twopi = np.pi+np.pi
    #
    g1 = 0.86602540378 * (2.-ssa*(1+hg))
    g2 = (1.7320508075688772*ssa/2.) * (1-hg)
    g2[g2 == 0.] = 1E-10
    g3 = (1.-1.7320508075688772*hg*MU0)/2
    g4 = 1. - g3
    
    lam = np.sqrt(g1*g1 - g2*g2)
    gamma = (g1-lam)/g2
    alpha = np.sqrt( (1.-ssa) / (1.-ssa*hg) )
    
    Am = F0PI * ssa *(g4 * (g1 + 1./MU0) + g2*g3 )/ (lam*lam - 1./(MU0*MU0))
    Ap = F0PI * ssa *(g3 * (g1 - 1./MU0) + g2*g4 )/ (lam*lam - 1./(MU0*MU0))
    
    # Cpm1 and Cmm1 are the C+ and C- terms evaluated at the top of the layer.
    Cpm1 = Ap * np.exp(-taut/MU0)
    Cmm1 = Am * np.exp(-taut/MU0)
    # Cp and Cm are the C+ and C- terms evaluated at the bottom of the layer.
    Cp = Ap * np.exp(-taub/MU0)
    Cm = Am * np.exp(-taub/MU0)
    
    #  Solve for the coefficients of system of equations using boundary conditions
    # Exponential terms:
    exptrm = lam*dtau
    exptrm[exptrm>35] = 35 # clipped so that exponential doesn't explode
    Ep = np.exp(exptrm)
    Em = 1./Ep
    
    E1 = Ep + gamma*Em
    E2 = Ep - gamma*Em
    E3 = gamma*Ep + Em
    E4 = gamma*Ep - Em
    
    L = nlay+nlay
    Af = np.empty(L)
    Bf = np.empty(L)
    Cf = np.empty(L)
    Df = np.empty(L)
    
    # First Term
    Af[0] = 0.0
    Bf[0] = gamma[0] + 1.
    Cf[0] = gamma[0] - 1.
    Df[0] = BTOP - Cmm1[0]
    
    AA = (E1[:-1]+E3[:-1])*(gamma[1:]-1)
    BB = (E2[:-1]+E4[:-1])*(gamma[1:]-1)
    CC = 2.*(1.-gamma[1:]*gamma[1:])
    DD = (gamma[1:]-1) * (Cpm1[1:] - Cp[:-1]) + (1-gamma[1:]) * (Cm[:-1]-Cmm1[1:])
    Af[1:-1:2]=AA
    Bf[1:-1:2]=BB
    Cf[1:-1:2]=CC
    Df[1:-1:2]=DD
    
    AA = 2.*(1.-gamma[:-1]*gamma[:-1])
    BB = (E1[:-1]-E3[:-1])*(gamma[1:]+1.)
    CC = (E1[:-1]+E3[:-1])*(gamma[1:]-1.)
    DD = E3[:-1]*(Cpm1[1:] - Cp[:-1]) + E1[:-1]*(Cm[:-1] - Cmm1[1:])
    Af[2::2]=AA
    Bf[2::2]=BB
    Cf[2::2]=CC
    Df[2::2]=DD
    # Last term:
    Af[-1] = E1[-1] - rsurf*E3[-1]
    Bf[-1] = E2[-1] - rsurf*E4[-1]
    Cf[-1] = 0.0
    Df[-1] = bsurf - Cp[-1] + rsurf*Cm[-1]
    
    k=tri_diag_solve(L, Af, Bf, Cf, Df)
    
    # Unmix coefficients
    even = np.arange(0,2*nlay,2)
    odd  = even+1
    k1 = k[even] + k[odd]
    k2 = k[even] - k[odd]
    
    Fup0=k1*Ep+gamma*k2*Em+Cpm1
    Fdn0=k1*Ep*gamma+k2*Em+Cmm1
    
    Cpmid = Ap * np.exp(-taump/MU0)
    Cmmid = Am * np.exp(-taump/MU0)
    Fup_diffuse=k1*Ep+gamma*k2*Em+Cpmid
    Fdn_diffuse=k1*Ep*gamma+k2*Em+Cmmid
    
    Fdn=np.zeros(nlay+1)
    Fup=np.zeros(nlay+1)
    Fdn[0:-1]=Fdn_diffuse+MU0*F0PI*np.exp(-taump/MU0)
    Fup[0:-1]=Fup_diffuse
    
    return Fup, Fdn




@jit(nopython=True)
def compute_RT(wnocrop,T,kcoeffs_interp,xsecContinuum,Mies,mass_path,
    Fractions,Fractions_Continuum,Fractions_Cond, gord,wts, Fstar):
    """
    Calls all requisite radiative transfer routines for emission
    Parameters
    ----------
    wnocrop : ndarray 
        wavenumber grid 
    T : ndarray 
        temperatures 
    kcoeffs_interp : ndarray 
        Interpolated k-coefficients 
    xsecContinuum : ndarray 
        cross sections for continuum opacity 
    Mies : ndarray 
        Mie scattering parameters 
    mass_path : ndarray 
    Fractions : ndarray
        volume mixing ratios of all absorbers (excluding continuum)
    Fractions_Continuum : ndarray
        volume mixing ratios for continuum 
    Fractions_Cond : ndarray
        volume mixing ratios for condensates 
    gord : ndarray
        g-oordinates for k tables
    wts : ndarray 
        weights for k tables 
    Fstar : ndarray
        Flux of parent star
    Returns
    -------
    """
    mu0=0.5#1./np.sqrt(3.) #what should I do with this???
    Nwno=len(wnocrop)
    Nlay=kcoeffs_interp.shape[0]
    Fuparr=np.zeros((Nwno, Nlay+1))
    Fdnarr=np.zeros((Nwno, Nlay+1))
    Fdnarr_star=np.zeros((Nwno, Nlay+1))
    Fuparr_star=np.zeros((Nwno, Nlay+1))
    dtauarr=np.zeros((Nlay,len(gord),Nwno))
    dtau_condarr=np.zeros((Nlay,len(gord),Nwno))
    dtau_gasarr=np.zeros((Nlay,len(gord),Nwno))
    ssaarr=np.zeros((Nlay,len(gord),Nwno))
    asymarr=np.zeros((Nlay,len(gord),Nwno))
    for i in range(Nwno): #looping over wavenumber
        #print 'BEGIN compute_RT WNO Loop: ', datetime.datetime.now().time()
        B=blackbody(T,1E4/wnocrop[i]*1E-6)
        #print 'Compute tau: ', datetime.datetime.now().time()
        dtau,dtau_cond,dtau_gas,ssa, asym=compute_tau(kcoeffs_interp[:,i,:,:],xsecContinuum[:,i,:],Mies[:,:,:,i], mass_path, Fractions,Fractions_Continuum, Fractions_Cond,gord, wts)
        dtauarr[:,:,i]=dtau
        ssaarr[:,:,i]=ssa
        asymarr[:,:,i]=asym
        dtau_condarr[:,:,i]=dtau_cond
        dtau_gasarr[:,:,i]=dtau_gas
    
        #print 'Looping over Gords: ', datetime.datetime.now().time()
        for j in range(len(gord)): #looping over g-space
            Fup_k,Fdn_k=toon(dtau[:,j], ssa[:,j]*1. ,asym[:,j]*1. , B)  #toon
            Fup_star_k, Fdn_star_k=stellar_flux(dtau[:,j], ssa[:,j]*1. ,asym[:,j]*1.,Fstar[i],mu0)
            Fuparr[i,:]+=Fup_k*wts[j]
            Fdnarr[i,:]+=Fdn_k*wts[j]
            Fdnarr_star[i,:]+=Fdn_star_k*wts[j]*1.
            Fuparr_star[i,:]+=Fup_star_k*wts[j]*1.

    return 0.5*Fuparr+0.5*Fuparr_star, 0.5*Fdnarr+0.5*Fdnarr_star, 0.5*Fuparr,0.5*Fuparr_star , dtauarr, ssaarr, asymarr #0.5 is b/c g-ordinates go from -1 - 1 instead of 0 -1


@jit(nopython=True)
def stellar_flux(dtau,ssa,hg, Fstar0, mu0):  #need to give this dtau and mu0
    """ 
    """
    rsfc=0.0
    Fup_s,Fdn_s=toon_solar(dtau, ssa, hg, rsfc,mu0,Fstar0,0.)
    return 0.5*Fup_s, 0.5*Fdn_s


def evaluate_model(wlgrid, model_parameters, model, cross_sections):
    # from copy import deepcopy

    def _assign_priors(model_params, cross_sections):
        """Programatically assign the values to cross sections given a dictionary
        of key value pairs.

        cross_sections : exoctk.chimera.cross_sections
            Cross sections object generated by LoadCrossSections
        model_params : dict
            Dictionary of key, value pairs to update cross sections with.
        """

        cross_sections_copy = deepcopy(cross_sections)
        
        for pname, pvalue in model_params.items():
            # Check the attributes of cross sections object 
            for attribute in vars(cross_sections_copy):
                # The priors defined are nested into dictionaires broken up by the type of parameter.
                # planetary_parameters, cloud_parameters, chemistry_parameters etc so here we are just
                # checking to see if the attributes are dictionaries because thats where there keys for
                # the priors are located.
                if isinstance(vars(cross_sections_copy)[attribute], dict):
                    # if attribute is a dictionary, see if the prior key is an option in the dictionary.
                    if pname in vars(cross_sections_copy)[attribute].keys():
                        # if prior is in dictionary keys, assign the prior to transformed value.
                        vars(cross_sections_copy)[attribute][pname] = pvalue
                    else:
                        continue
        
        return cross_sections_copy


    # Make a copy of the cross_sections obj to manipulate to evaluate model.
    temp_cross_sections = _assign_priors(model_parameters, cross_sections)
    
    temp_model = deepcopy(model)
    temp_model.wlgrid = wlgrid
    temp_model.fx_trans_free()

    spec = temp_model.tran(temp_cross_sections, temp_model.T, temp_model.atmosphere_grid['P'], temp_model.mmw, temp_model.Pref, temp_model.CldOpac, temp_model.H2Oarr, temp_model.CH4arr, 
                            temp_model.COarr, temp_model.CO2arr, temp_model.NH3arr, temp_model.Naarr, temp_model.Karr, temp_model.TiOarr, temp_model.VOarr,
                            temp_model.C2H2arr, temp_model.HCNarr, temp_model.H2Sarr, temp_model.FeHarr, temp_model.Harr, temp_model.earr, temp_model.Hmarr, 
                            temp_model.H2arr, temp_model.Hearr, temp_model.RayAmp, temp_model.RaySlp, temp_model.f_r, temp_cross_sections.planetary_parameters['M'], temp_cross_sections.stellar_parameters['Rstar'], 
                            temp_cross_sections.planetary_parameters['Rp'])

    wno = spec[0]
    F = spec[1]

    y_binned, _ = model.instrument_tran_non_uniform(wlgrid, wno, F)

    return y_binned, temp_model
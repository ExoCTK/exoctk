import h5py
import numpy as np
import os
import pysynphot as psyn
from scipy import interp


class PlanetarySystem():
    """Class object defining the planetary system you wish to model with Chimera"""

    def __init__(self):
        """Initialize the class object."""

        self.directory = None
        self.planetary_parameters = None
        self.output_plot = 'model.png'
        self.observatory = None
        self.stellar_parameters = None
        self.wavenumber_min = None
        self.wavenumber_max = None


    def make_plot(self):
        """Plot Model"""
        return None


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

        database = self.stellar_parameters['stellar_db']
        temp = self.stellar_parameters['temp']
        logMH = self.stellar_parameters['logmh']
        logg = self.stellar_parameters['logg']
        outfile = self.stellar_parameters['stellar_file']

        sp = psyn.Icat(database, temp, logMH, logg)
        sp.convert('flam') # initial units erg/cm^2/s/ Angst
        wave=sp.wave*1e-10  # in meters
        fluxmks = sp.flux*1e-7*1e4*1e10 # in W/m2/m
        f = h5py.File(outfile,'w')
        f.create_dataset('lambdastar',data=wave)
        f.create_dataset('Fstar0',data=fluxmks)
        f.close()


    def calculate_cross_sections(self):
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
            Fstar0=interp(np.log10(lambdastar_hi), np.log10(lambdastar), Fstar0) 

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
        file=os.path.join(self.directory, 'MIE_COEFFS',self.cond_name+'_r_0.01_300um_wl_0.3_200um_interp_R100_20gp_50_30000wno.h5')
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
        wno=wno[loc]
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

        gases = np.array(hf['gases'])  ##H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He   e- h-  mmw

        self.gases = np.log10(gases)
        hf.close()

        self.radius = radius*1E-6
        print('Cross-sections Loaded')

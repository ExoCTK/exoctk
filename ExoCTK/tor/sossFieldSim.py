from astroquery.irsa import Irsa
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.coordinates as crd
import astropy.units as u
import numpy as np
import os
import glob
import idlsave
import pdb

idlsave_path = os.environ.get('EXOCTK_CONTAM_DIR')

def sossFieldSim(ra, dec, binComp='', dimX=256):
    # binComp: [deltaRA,deltaDEC,J,H,K]

    # stars in large field around target
    targetcrd = crd.SkyCoord(ra = ra, dec = dec, unit=(u.hour, u.deg))
    targetRA  = targetcrd.ra.value
    targetDEC = targetcrd.dec.value
    info      = Irsa.query_region(targetcrd, catalog = 'fp_psc', spatial = 'Cone', radius = 2.5*u.arcmin) 
	
    # coordinates of all stars in FOV, including target
    allRA   = info['ra'].data.data 
    allDEC  = info['dec'].data.data
    Jmag    = info['j_m'].data.data
    Hmag    = info['h_m'].data.data
    Kmag    = info['k_m'].data.data
    J_Hobs  = Jmag-Hmag
    H_Kobs  = Hmag-Kmag
    
    # target coords
    distance    = np.sqrt( ((targetRA-allRA)*np.cos(targetDEC))**2 + (targetDEC-allDEC)**2 )
    targetIndex = np.argmin(distance) # the target
    
    # add any missing companion
    cubeNameSuf=''
    if binComp!='':
        deg2rad     = np.pi/180
        allRA       = np.append(allRA, (allRA[targetIndex] + binComp[0]/3600/np.cos(allDEC[targetIndex]*deg2rad)))
        allDEC      = np.append(allDEC, (allDEC[targetIndex] + binComp[1]/3600))
        Jmag        = np.append(Jmag,binComp[2])
        Hmag        = np.append(Kmag,binComp[3])
        Kmag        = np.append(Kmag,binComp[4])
        J_Hobs      = Jmag-Hmag
        H_Kobs      = Hmag-Kmag
        cubeNameSuf ='_custom'          

    #number of stars
    nStars=allRA.size

    cooTar=crd.SkyCoord(ra=allRA[targetIndex],dec=allDEC[targetIndex], unit=(u.deg, u.deg))
    # cubeName = 'cubes/cube_RA_DEC_'+cooTar.to_string('hmsdms',precision=1).replace(' ','')+'.fits'
    #
    # print('Cube name:')
    # print(cubeName)
    # print('cubNameSuf:')
    # print(cubeNameSuf)
    # print('Target coordinates:',cooTar.to_string('hmsdms',precision=1))
    #
    # if os.path.exists(cubeName) and (binComp is None):
    #     return

    """
    # the trace models
    if 0:
        modelDir  = '/Users/david/Documents/work/jwst/niriss/soss/data'

        # synthetic mag & colors of models
        readcol,modelDir+'/'+'magnitude_modele_atm.txt',Ttheo,logg,MH,aH,Jtheo,Htheo,Ktheo,jhtheo,jktheo,hktheo,format='f,f,f,f,f,f,f,f,f,f'
        remove,where(logg ne 4.50 or Ttheo lt 2800),Ttheo,logg,jhtheo,jktheo,hktheo 

        # load model traces
        fNameMod=file_search(modelDir,'simu_lte*-4.5-0.0.BT-Settl_LLNL_UDEMPSF_cropped.fits.gz',count=nMod)
        teffMod=fix(strmid(fNameMod,strpos(fNameMod[0],'simu_lte')+8,3))*100
        isort=sort(teffMod)
        teffMod=teffMod[isort]
        fNameMod=fNameMod[isort]

        # cropped model traces are cubes whose third axis are: 0->sum of all orders, 1->order 1, 2->order 2
        h0=headfits(fNameMod[0])
        dimXmod=sxpar(h0,'naxis1')
        dimYmod=sxpar(h0,'naxis2')
        models=fltarr(dimXmod,dimYmod,nMod)
        jhMod=fltarr(nMod)
        hkMod=fltarr(nMod)
        for i=0,nMod-1 do begin
            tmp=readfits(fNameMod[i])
            models[*,*,i]=tmp[*,*,0] # keep only the sum of all orders

            modelO12=tmp[*,*,1:2]
            fNameModO12='modelOrder12_teff'+strtrim(teffMod[i],2)+'.sav'
            save,modelO12,file=fNameModO12
        
            j=where(Ttheo eq teffMod[i])
            jhMod[i]=jhTheo[j]
            hkMod[i]=hkTheo[j]  
        endfor
        modelPadX = 2000-1900
        modelPadY = 2000

        # save the model traces
        save,models,fNameModO12,teffMod,jhMod,hkMod,dimXmod,dimYmod,modelPadX,modelPadY,Ttheo,logg,jhtheo,jktheo,hktheo,file='modelsInfo.sav'
    endif else begin
        restore,'modelsInfo.sav'
    endelse
    """
    
    #Restoring model parameters 
    modelParam = idlsave.read(os.path.join(idlsave_path,'modelsInfo.sav'),verbose=False) 
    models     = modelParam['models']
    modelPadX  = modelParam['modelpadx']
    modelPadY  = modelParam['modelpady'] 
    dimXmod    = modelParam['dimxmod']
    dimYmod    = modelParam['dimymod']
    jhMod      = modelParam['jhmod']
    hkMod      = modelParam['hkmod']
    teffMod    = modelParam['teffmod'] 

    # find/assign Teff of each star
    starsT=np.empty(nStars)
    for j in range(nStars):
        color_separation = (J_Hobs[j]-jhMod)**2+(H_Kobs[j]-hkMod)**2
        min_separation_ind = np.argmin(color_separation)
        starsT[j]=teffMod[min_separation_ind]

    radeg = 180/np.pi
    niriss_pixel_scale = 0.065  # arcsec
    sweetSpot = dict(x=856,y=107,RA=allRA[targetIndex],DEC=allDEC[targetIndex],jmag=Jmag[targetIndex])

    #offset between all stars and target
    dRA=(allRA - sweetSpot['RA'])*np.cos(sweetSpot['DEC']/radeg)*3600
    dDEC=(allDEC - sweetSpot['DEC'])*3600
    
    # Put field stars positions and magnitudes in structured array
    _ = dict(RA=allRA, DEC=allDEC, dRA=dRA, dDEC=dDEC, jmag=Jmag, T=starsT,
             x=np.empty(nStars), y=np.empty(nStars), dx=np.empty(nStars), dy=np.empty(nStars))
    stars=np.empty(nStars,dtype=[(key,val.dtype) for key,val in _.items()])
    for key,val in _.items(): stars[key]=val
   
    # Initialize final fits cube that contains the modelled traces with contamination
    PAmin = 0 #instrument PA, degrees
    PAmax = 360
    dPA = 1 # degrees
    # Set of IPA values to cover
    PAtab = np.arange(PAmin, PAmax, dPA)    # degrees
    nPA    = len(PAtab)

    # dimX=256 #2048 #########now as argument, with default to 256
    dimY=2048
    simuCube=np.zeros([nPA+2,dimY, dimX])  # cube of trace simulation at every degree of field rotation, +target at O1 and O2

    # saveFiles = glob.glob('idlSaveFiles/*.sav')[:-1]
    saveFiles = glob.glob(os.path.join(idlsave_path,'*.sav'))[:-1]
    #pdb.set_trace()
    
    # Big loop to generate a simulation at each instrument PA
    for kPA in range(PAtab.size):
        APA= PAtab[kPA]
        V3PA=APA+0.57 #from APT
    
        stars['dx']= (np.cos(np.pi/2+APA/radeg)*stars['dRA']-np.sin(np.pi/2+APA/radeg)*stars['dDEC'])/niriss_pixel_scale
        stars['dy']= (np.sin(np.pi/2+APA/radeg)*stars['dRA']+np.cos(np.pi/2+APA/radeg)*stars['dDEC'])/niriss_pixel_scale
        stars['x'] = stars['dx']+sweetSpot['x']
        stars['y'] = stars['dy']+sweetSpot['y']

        # Display the star field (blue), target (red), subarray (green), full array (blue), and axes
        if (kPA==0 and nStars > 1) and False:
            print(kPA)
            plt.plot([0,2047,2047,0,0],[0,0,2047,2047,0], 'b')
            plt.plot([0,255,255,0,0],[0,0,2047,2047,0], 'g')
            #the order 1 & 2 traces
            t1=np.loadtxt('/Users/david/Documents/work/jwst/niriss/soss/data/trace_order1.txt',unpack=True)
            plt.plot(t1[0],t1[1],'r')
            t2=np.loadtxt('/Users/david/Documents/work/jwst/niriss/soss/data/trace_order2.txt',unpack=True)
            plt.plot(t2[0],t2[1],'r')
            plt.plot(stars['x'], stars['y'], 'b*')
            plt.plot(sweetSpot['x'], sweetSpot['y'], 'r*')
            plt.title("APA= {} (V3PA={})".format(APA,V3PA))
            ax=plt.gca()

            #add V2 & V3 axes
            l,hw,hl=250,50,50
            adx,ady=-l*np.cos(-0.57/radeg),-l*np.sin(-0.57/radeg)
            ax.arrow(2500, 1800, adx,ady, head_width=hw, head_length=hl, length_includes_head=True, fc='k') #V3
            plt.text(2500+1.4*adx,1800+1.4*ady,"V3",va='center',ha='center')
            adx,ady=-l*np.cos((-0.57-90)/radeg),-l*np.sin((-0.57-90)/radeg)
            ax.arrow(2500, 1800, adx, ady, head_width=hw, head_length=hl, length_includes_head=True, fc='k') #V2
            plt.text(2500+1.4*adx,1800+1.4*ady,"V2",va='center',ha='center')
            #add North and East
            adx,ady=-l*np.cos(APA/radeg),-l*np.sin(APA/radeg)
            ax.arrow(2500, 1300, adx, ady, head_width=hw, head_length=hl, length_includes_head=True, fc='k') #N
            plt.text(2500+1.4*adx,1300+1.4*ady,"N",va='center',ha='center')
            adx,ady=-l*np.cos((APA-90)/radeg),-l*np.sin((APA-90)/radeg)
            ax.arrow(2500, 1300, adx, ady, head_width=hw, head_length=hl, length_includes_head=True, fc='k') #E
            plt.text(2500+1.4*adx,1300+1.4*ady,"E",va='center',ha='center')

            ax.set_xlim(-400,2047+800)
            ax.set_ylim(-400,2047+400)
            ax.set_aspect('equal')
            plt.show()
 
            # email from Michael Maszkiewicz on 2016 May 17:
            #  POM edge is:
            #     1.67-2.01 mm left of detector in x
            #     2.15-2.22 mm above detector in y
            #     2.22-2.43 mm right of detector in x
            #     1.49-1.87 mm below detector in y
            #  we take the larger value in each direction and add a 50 pixels margin
            #  given the uncertainty of the POM location wrt to detector FOV

        # Retain stars that are within the Direct Image NIRISS POM FOV
        ind, = np.where((stars['x'] >= -162) & (stars['x'] <= 2047+185) & (stars['y'] >= -154) & (stars['y'] <= 2047+174))
        starsInFOV=stars[ind]
        
    #   print(' Number of stars in FOV:',ind.size)
    #   print(' J magnitude of those stars:', starsInFOV['Jmag'])
    #   print(' Temperature of those stars:', starsInFOV['T'])

        for i in range(len(ind)):
            intx = round(starsInFOV['dx'][i])
            inty = round(starsInFOV['dy'][i])
            # print(intx,inty)

            k=np.where(teffMod == starsInFOV['T'][i])[0][0]
            
            fluxscale = 10.0**(-0.4*(starsInFOV['jmag'][i] - sweetSpot['jmag']))

            #deal with subection sizes
            mx0=int(modelPadX-intx)
            mx1=int(modelPadX-intx+dimX)
            my0=int(modelPadY-inty)
            my1=int(modelPadY-inty+dimY)
            
            if (mx0 > dimXmod) or (my0 > dimYmod):
                continue
            if (mx1 < 0) or (my1 < 0):
                continue
            
            x0  =(mx0<0)*(-mx0)
            y0  =(my0<0)*(-my0)
            mx0 *=(mx0 >= 0)
            mx1 = dimXmod if mx1>dimXmod else mx1
            my0 *=(my0 >= 0)
            my1 =dimYmod if my1>dimYmod else my1

            # if target and first kPA, add target traces of order 1 and 2 in output cube
            if (intx == 0) & (inty == 0) & (kPA == 0):              
                fNameModO12 = saveFiles[k]
                modelO12 = idlsave.read(fNameModO12,verbose=False)['modelo12']
                simuCube[0, y0:y0+my1-my0, x0:x0+mx1-mx0] = modelO12[0, my0:my1, mx0:mx1] * fluxscale # order 1
                simuCube[1, y0:y0+my1-my0, x0:x0+mx1-mx0] = modelO12[1, my0:my1, mx0:mx1] * fluxscale # order 2
                
            if (intx != 0) or (inty != 0): #field star
                simuCube[kPA+2, y0:y0+my1-my0, x0:x0+mx1-mx0] += models[k, my0:my1, mx0:mx1] * fluxscale
    
    # fits.writeto(cubeName, simuCube, clobber = True)
    # print(cubeName)
    
    return simuCube


if __name__=='__main__':
    ra,dec="04 25 29.0162","-30 36 01.603" #Wasp 79
    sossFieldSim(ra,dec)

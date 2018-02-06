#Produces a graph of the visibility & accessible position angles
#for a given RA & DEC, and prints out corresponding information,
#including the ranges of accessible and inaccessible PAs.
#
#Usage: python visibilityPA.py RA DEC [targetName]
# if targetName is specified, then the figure is saved
#
#-Created by David Lafreniere, March 2016
#-makes use of (and hacks) several scripts created by Pierre Ferruit
# that are part of the JWST Python tools JWSTpylib and JWSTpytools


import sys
#import ExoCTK
#import ExoCTK.tor.contam_tool.ephemeris_old2x as EPH
import ephemeris_old2x as EPH
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
from astropy.io import ascii
import os, base64, io, mpld3


def calc_vis(RA, DEC, targetName=None):
	D2R = math.pi / 180.  #degrees to radians
	R2D = 180. / math.pi #radians to degrees 

	def convert_ddmmss_to_float(astring):
		aline = astring.split(':')
		d= float(aline[0])
		m= float(aline[1])
		s= float(aline[2])
		hour_or_deg = (s/60.+m)/60.+d
		return hour_or_deg

	def checkVisPA(ra,dec,targetName=None,save=False):
		if ra.find(':')>-1:  #format is hh:mm:ss.s or  dd:mm:ss.s  
			ra  = convert_ddmmss_to_float(ra) * 15. * D2R
			dec   = convert_ddmmss_to_float(dec) * D2R
		else: #format is decimal
			ra  = float(ra) * D2R
			dec   = float(dec) * D2R

		#load ephemeris
		#ephFileName,eclFlag=os.path.join(os.path.dirname(ExoCTK.__file__),'data/tor/JWST_ephem_short.txt'),False
		ephFileName,eclFlag='JWST_ephem_short.txt',False
		eph = EPH.Ephemeris(ephFileName, eclFlag)
		#convert dates from MJD to Gregorian calendar dates
		mjd=np.array(eph.datelist)
		d=mdates.julian2num(mjd+2400000.5)
		gd=mdates.num2date(d)

		#loop through dates and determine VIS and PAs (nominal, min, max)
		vis=np.empty(mjd.size,dtype=bool)
		paNom,paMin,paMax=np.empty(mjd.size),np.empty(mjd.size),np.empty(mjd.size)
		for i in range(mjd.size):
			vis[i]=eph.in_FOR(mjd[i],ra,dec) #is it visible?
			pa=eph.normal_pa(mjd[i],ra,dec) #nominal PA at this date
	
			#search for minimum PA allowed by roll
			pa0=pa
			while eph.is_valid(mjd[i],ra,dec,pa0-0.002):
				pa0-=0.002
			#search for maximum PA allowed by roll
			pa1=pa
			while eph.is_valid(mjd[i],ra,dec,pa1+0.002):
				pa1+=0.002
	
			paNom[i]=(pa*R2D)%360
			paMin[i]=(pa0*R2D)%360
			paMax[i]=(pa1*R2D)%360

		#does PA go through 360 deg?
		wrap=np.any(np.abs(np.diff(paNom[np.where(vis)[0]])) > 350)

		#Determine good and bad PA ranges
		#Good PAs
		i,=np.where(vis)
		pa=np.concatenate((paNom[i],paMin[i],paMax[i]))
		if wrap: pa=np.append(pa,(0.,360.))
		pa.sort()
		i1,=np.where(np.diff(pa)>10)
		i0=np.insert(i1+1,0,0)
		i1=np.append(i1,-1)
		paGood=np.dstack((pa[i0],pa[i1])).round(1).reshape(-1,2).tolist()

		#bad PAs (complement of the good PAs)
		paBad=[]
		if paGood[0][0]>0: paBad.append([0.,paGood[0][0]])
		for i in range(1,len(paGood)):
			paBad.append([paGood[i-1][1],paGood[i][0]])
		if paGood[-1][1]<360.: paBad.append([paGood[-1][1],360.])

		#print results to file
		"""
		if save:
			fName='visibilityPA-'+targetName+'.txt'
			fic=open(fName,'w')

			fic.write('#Date    MJD          VIS?  PAnom   PArange\n')
			for i in range(vis.size):
				tmp1='{:7.3f}'.format(paNom[i]) if vis[i] else 7*'-'
				tmp2='{:7.3f}--{:7.3f}'.format(paMin[i],paMax[i]) if vis[i] else 16*'-'
				#fic.write(gd[i].strftime("%y-%m-%d")+' {:f} {:5s} {:7.3f} {:7.3f}--{:7.3f} \n'.format(mjd[i],str(vis[i]),paNom[i],paMin[i],paMax[i]))
				fic.write(gd[i].strftime("%y-%m-%d")+' {:f} {:5s} {} {} \n'.format(mjd[i],str(vis[i]),tmp1,tmp2))

			fic.write("\n")
			fic.write("Accessible PA ranges: ")
			fic.write(','.join([str(x) for x in paGood]))
			fic.write("\n")
			fic.write("Non-accessible PA ranges: ")
			fic.write(','.join([str(x) for x in paBad]))
			fic.write("\n")
			fic.close()
		"""

		#make the plot
		fig=plt.figure()
		ax=plt.axes()
		plt.title(targetName)
	

		#plot nominal PA
		iBad,=np.where(vis==False)
		paMasked=np.copy(paNom)
		paMasked[iBad]=np.nan
		gdMasked=np.copy(gd)
		if wrap:
			i=np.argmax(paNom)
			if paNom[i+1]<10: i+=1
			paMasked=np.insert(paMasked,i,np.nan)
			gdMasked=np.insert(gdMasked,i,gdMasked[i])
		plt.plot(gdMasked,paMasked,color='k')
	# 		plt.xticks(rotation = 'vertical')

		#plot ranges allowed through roll
		if wrap:
			i=np.argmax(paMin)
			goUp=paMin[i-2]<paMin[i-1] #PA going up at wrap point?
			#top part
			i0=0 if goUp else i
			i1=i if goUp else paMin.size-1
			paMaxTmp=np.copy(paMax)
			paMaxTmp[np.where(paMin>paMax)[0]]=360
			plt.fill_between(gd[i0:i1+1],paMin[i0:i1+1],paMaxTmp[i0:i1+1],where=vis[i0:i1+1],lw=0,facecolor='k',alpha=0.5)
			#bottom part
			i=np.argmin(paMax)
			i0=i if goUp else 0
			i1=paMin.size-1 if goUp else i 
			paMinTmp=np.copy(paMin)
			paMinTmp[np.where(paMin>paMax)[0]]=0
			plt.fill_between(gd[i0:i1+1],paMinTmp[i0:i1+1],paMax[i0:i1+1],where=vis[i0:i1+1],lw=0,facecolor='k',alpha=0.5)
		else:
			plt.fill_between(gd,paMin,paMax,where=vis,lw=0,facecolor='k',alpha=0.5)

		plt.ylabel('Position Angle (degrees)')
		plt.xlim(min(gd),max(gd))
		ax.xaxis.set_major_locator(mdates.MonthLocator())
		ax.xaxis.set_major_formatter(mdates.DateFormatter("%b '%y"))
		ax.xaxis.set_minor_locator(mdates.DayLocator(list(range(1,32,5))))
		plt.ylim(0,360)
		ax.yaxis.set_major_locator(MultipleLocator(25))
		ax.yaxis.set_minor_locator(MultipleLocator(5))
		plt.grid()
		for label in ax.get_xticklabels():
			label.set_rotation(45)
		if save:
	# 			fname=tmpDir+'/visibilityPA-'+targetName+'.png'
			png = mpld3.fig_to_html(fig)
		

		return png, paGood, paBad

	#arguments RA & DEC, conversion to radians
	save=False if targetName==None else True
# 	os.makedirs(tmpDir, exist_ok=True)

	png, pG, pB = checkVisPA(RA, DEC,targetName=targetName,save=save)

	return pG, pB, png

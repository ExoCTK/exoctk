"""
This file is meant to be used as a module for the ExoCTK website i.e. 
app_exoctk.py. It is part of the ExoCTK package. It produces a graph of 
the visibility & accessible position angles for a given RA & DEC, and 
prints out corresponding information, including the ranges of accessible 
and inaccessible PAs.

Authors:
	Rafia Bushra, University of Arizona
	
	David Lafreniere, University de Montreal

Usage:
	* From python script:
		from ExoCTK.tor.tor2.visibilityPA import calc_vis
		calc_vis(ra, dec, targetName)
		
	* From terminal:
		This function returns a figure in form of bytes. So I didn't keep a terminal 
		version, unless you want to use ipython from terminal and import calc_vis from 
		visibilityPA 
		
	
Requirements:
	* ExoCTK package must be installed
	* ephemeris_old2x.py must be in the same folder
	* JWST_ephem_short.txt must be in ExoCTK/data/tor
"""


import sys
import ExoCTK
from . import ephemeris_old2x as EPH
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
from bokeh import mpl
from bokeh.plotting import show
from bokeh.models import DatetimeTickFormatter 
import os
import base64
import io


def calc_vis(RA, DEC, targetName=None):

	"""
	This function prduces the visibility plot and saves it to memory
	Args:
		ra (str)        : Right Ascension of target in hh:mm:ss
		dec (str)       : Declination of target in ddd:mm:ss
		targetName(str) : Name of target
		
	Returns:
		paGood (list)   : list os position angle where target is accessible
		paBad (list)    : list os position angle where target is not accessible
		png (bytes)     : visibility plot in bytes form that should be passed on to an 
		   				  html file
	"""
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
		ephFileName,eclFlag=os.path.join(os.path.dirname(ExoCTK.__file__),'data/tor/JWST_ephem_short.txt'),False
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
# 			png = mpld3.fig_to_html(fig)
			# bk = mpl.to_bokeh(fig)
# 			bk.plot_width=400
# 			show(bk)
			buff = io.BytesIO()
			plt.savefig(buff, format = 'png')
			buff.seek(0)
			figdata_png = base64.b64encode(buff.getvalue()).decode('ascii')

		return paGood, paBad, figdata_png

	save=False if targetName==None else True

	paGood, paBad, png = checkVisPA(RA, DEC, targetName=targetName, save=save)

	return paGood, paBad, png

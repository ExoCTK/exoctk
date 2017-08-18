"""
This file is meant to be used as a module for the ExoCTK website i.e. 
app_exoctk.py. It is part of the ExoCTK package. It runs the simulation 
for two filters of NIRCam and  LRS mode of MIRI. In essence, it prodeuces 
figures that are similar to those produced by sossContamFig.py.

Authors:
	Rafia Bushra, University of Arizona
	
	Jonathan Fraine, Space Telescope Science Institute
	
	Tom Greene, NASA Ames Research Center
	
Usage:
	* From python script:
		from ExoCTK.tor.tor2.ncam_miri import *
		ncam_miri(target_name, instrument)
		
	* From terminal:
		cd into ExoCTK/ExoCTK/tor/tor2 directory
		$ python ncam_miri.py target_name instrument
		
	
Requirements:
	* ExoCTK package must be installed
"""




# T. Greene 04 April 2017 tom.greene@nasa.gov; copied Simbad.query & Irsa.query from J Fraine

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm
from sys import argv
import numpy as np
import math as math

from astroquery.irsa import Irsa
from astroquery.simbad import Simbad
import astropy.units as u
from tqdm import tqdm
import io
import base64

D2R = math.pi / 180.0

# Input targetNamename and V3PA here
def ncam_miri(targetName, ins = 'NIRCam_F322W2'):
	
	"""
	This function runs simulation for spectral overlap and prodeuces a contamination figure
	
	Args:
		targetName (str) : Name of Target
		ins (str)        : Name of instrument. Four choices:
					NIRCam_F322W2 (default)
					NIRCam_f444W
					MIRI
		
	
	Returns:
		figdata_png (bytes) : Contamination plot in bytes form that should be passed on to a html file
	"""
	
	# PA     = 7  #position angle of field rel. to V3 for given observation date/visibility


	fontsize = 15

	# PA     = 50+90  #position angle of field re. to V3 for given observation date/visibility
	minPA  = 0
	xmax1PA  = 360
	spacePA= 1
	nPA    = (xmax1PA - minPA) // spacePA + 2

	setPA  = np.arange(minPA-spacePA, xmax1PA + spacePA, spacePA)

	dmag_limit = 7.5 # delta K mag limit for source consideration

	nWaves        = 2048

	y_gs= np.arange(-2048,2048)
	x_gs= np.arange(-2048,2048)



	# Get coordinates from SIMBAD and then find nearby 2MASS point sources

	target_Info  = Simbad.query_object(targetName)
	target_Data  = Irsa.query_region(targetName, catalog="fp_psc", spatial='Box',width=1 * u.arcsec)

	target_RA  = target_Info['RA'][0].split()
	target_DEC = target_Info['DEC'][0].split()

	RAd0  = (float(target_RA[0]) + float(target_RA[1])/60. + float(target_RA[2])/3600.)*15
	Decd0 = math.fabs(float(target_DEC[0])) + math.fabs(float(target_DEC[1])/60.) + math.fabs(float(target_DEC[2])/3600.)

	if (float(target_DEC[0]) < 0 or float(target_DEC[1]) < 0 or float(target_DEC[2]) < 0):
		Decd0 = -1.0 * Decd0

	K0     = np.sort(target_Data['k_m'])[0]#5.6 #K mag of host star

	sources = Irsa.query_region(targetName, catalog="fp_psc", spatial='Box',width=4 * u.arcmin)
	nsources = len(sources)


	# Set up arrays for x offset, yoffset, mag difference, plot size, and spectral collision
	dx = np.zeros(nsources)
	dy = np.zeros(nsources)
	dmag = np.zeros(nsources)
	msize = np.zeros(nsources)
	col = np.zeros(nsources)

	print(target_RA , target_DEC, RAd0, Decd0)

	dmag_limit = 7.5 # delta K mag limit for source contamination consideration


	def gaussian1D(center, width, height = None, offset = None):
		"""
			Written by Nate Lust
			Edited  by Jonathan Fraine

			Returns a 1D gaussian function with the given parameters

			center  = center of gaussian profile

			width   = width  of gaussian profile

			height  = height of gaussian profile
						-- defaults to `1 / np.sqrt(2.*pi*sigma**2.)`

			offset  = background, lower limit value for gaussian
						-- defaults to 0.0
		"""

		if height == None:
			height  = np.sqrt(2.*np.pi*width**2.)
			height  = 1.0/height

		if offset == None:
			offset = 0.0

		width   = float(width)

		return lambda x: height*np.exp(-(((center - x)/width)**2)/2) + offset

	def gaussian2D(center_y, center_x, width_y, width_x = None, height = None, offset = None):
		"""
			Written by Nate Lust
			Edited  by Jonathan Fraine

			Returns a 2D gaussian function with the given parameters

			center_y, center_x  = center position of 2D gaussian profile

			width_y , width_x   = widths of 2D gaussian profile (if width_y != width_x, then gaussian crossection = ellipse)

			height  = height of gaussian profile
						-- defaults to `1 / np.sqrt(2.*pi*sigma**2.)`

			offset  = background, lower limit value for gaussian
						-- defaults to 0.0
		"""

		if width_x == None:
			width_x = width_y

		if height == None:
			height = np.sqrt(2*np.pi*(width_x**2 + width_y**2))
			height = 1./height

		if offset == None:
			offset = 0.0

		width_x = float(width_x)
		width_y = float(width_y)

		return lambda y,x: height*np.exp(-(((center_x-x)/width_x)**2 + ( (center_y-y)/width_y)**2)/2)+offset

	def conv1D(arr1, arr2):
		'''
			Convolve 2 arrays together
				-- used by `smooth_gaussconv`
		'''
		fft1    = np.fft.fft(arr1)
		fft2    = np.fft.fft(arr2)

		conv    = np.fft.ifft(fft1*fft2)

		return np.real(np.fft.fftshift(conv))

	def conv2D(arr1, arr2):
		'''
			Convolve 2 arrays together
				-- used by `smooth_gaussconv`
		'''
		fft1    = np.fft.fft2(arr1)
		fft2    = np.fft.fft2(arr2)

		conv    = np.fft.ifft2(fft1*fft2)

		return np.real(np.fft.fftshift(conv))

	def smooth_gaussconv(arr, sigma):
		'''
			Gaussian smooths `arr` over a width of `sigma`
				-- uses `conv1D` and `conv2D`, which both use np.fft
		'''
		if len(arr.shape) == 1:
			gs1D    = gaussian1D(arr.size/2, sigma)(np.arange(arr.size))
			return conv1D(gs1D, arr) / gs1D.sum()

		if len(arr.shape) == 2:
			gs2D    = gaussian2D(arr.shape[0]/2., arr.shape[1]/2., sigma)(*np.indices(arr.shape))
			return conv2D(gs2D, arr) / gs2D.sum()


	dra   = sources['ra'][(sources['k_m'] - K0) < dmag_limit] - RAd0
	ddec  = sources['dec'][(sources['k_m'] - K0) < dmag_limit] - Decd0
	dmag  = sources['k_m'][(sources['k_m'] - K0) < dmag_limit] - K0
	msize = 10**(-0.1*dmag)*15
	nlimit= sum((sources['k_m'] - K0) < dmag_limit)
	col   = np.zeros(nlimit) # No spectral collision by default

	#change to 1 if you want to see a plot of the sources
	if 0:
		plt.plot(dra*3600,ddec*3600,'o')
		plt.plot(0,0,'o')

	"""	
	Set up for NIRCam, visualize R grism + F322W2
	Set Constants for Geometric Location of Targets on the Detector for F322W2
	"""

	#LWA R Field Point positions
	X_F322W2 = 479 # X Field point (undeviated); old val 514
	X_F444W = 1161 # X Field point undeviated); old val 1081
	X_buff = 20 # Keep other spectra this far away from xmin1 and xxmax1 (dist from edge)
	Y_field = 32  # Y Field point for locating object
	Y_buff = 50 # Keep other spectra this far away from Y_field (dist from center)

	xmaxWL_F322W2 = 4.013 # xmax1 wavelength of F322W2 in microns per https://jwst-docs.stsci.edu/display/JTI/NIRCam+Filters
	MINWL_F444W = 3.880 #minimum wavelength of F444W in microns per https://jwst-docs.stsci.edu/display/JTI/NIRCam+Filters

	wavelen_undev = 3.97
	disp = -0.001 #microns / pxl
	dx_F444W = 1106 # length of F444W spectrum in pixels
	dx_F322W2 = 1583 # length of F322W2 m=1 spectrum in pixels
	dx_F322W2m2 = math.fabs((xmaxWL_F322W2 - 2.4) * 2.0 / disp) # length of F322W2 m=2 spectrum in pixels

	xmin1_F322W2 = (xmaxWL_F322W2 - wavelen_undev)/disp + X_F322W2
	xmax_F322W2 = xmin1_F322W2 + dx_F322W2
	xmin1_F322W2m2 = (xmaxWL_F322W2*2.0 - wavelen_undev)/disp + X_F322W2
	xmax_F322W2m2 = xmin1_F322W2m2 + dx_F322W2m2
	# print(xmin1_F322W2m2, xmax_F322W2m2, dx_F322W2m2

	xxmax1_F444W = (MINWL_F444W -  wavelen_undev)/disp + X_F444W
	xmin1_F444W = xxmax1_F444W - dx_F444W
	# print(xmin1_F322W2, xmax_F322W2, xmin1_F444W, xxmax1_F444W

	# LWA field extremes for pickoff mirror POM
	xmin1_xmax1y = -142
	xmin1_miny = -132
	xxmax1_miny = 2254
	xxmax1_xmax1y = 2198

	yxmax1_minx = 2822
	ymin_minx = -174
	ymin_xmax1x = -180
	yxmax1_xmax1x = 2800


	nyquistWidth = 2. / 2.3548
	# contaminationF322W2 = np.zeros((nPA, nWaves))

	normalHeight   = 1/np.sqrt(2*np.pi*nyquistWidth*nyquistWidth)
	targetGaussian = gaussian1D(center=Y_field, width=nyquistWidth, height=normalHeight, offset=0)(y_gs+Y_field)


	if ins == 'NIRCam_F322W2':
		#Plot for F322W2
		xbox = [-0, -0, 4040, 4040, -0]
		ybox = [-0, 4040, 4040, -0, -0]

		def plot_spectral_overlap_F322W2(col, dra, ddec, Decd0, PA, msize, dmag, targetName):
			title = "%s, NCam LWAR F322W2, V3PA =%3.0f$^\circ$" % (targetName, PA)
	
			dx = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.cos((PA-0) * D2R) - ddec * 3600.0 / 0.065 * math.sin((PA-0) *D2R)
			dy = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.sin((PA-0) * D2R) + ddec * 3600.0 / 0.065 * math.cos((PA-0) *D2R)

			fig, ax = plt.subplots(figsize=(6,6))
			ax.set_title(title)
			majorLocator = MultipleLocator(500)
			minorLocator = MultipleLocator(100)

			ax.set_xlim([-500, 2300])
			ax.set_ylim([-500, 2300])
			ax.xaxis.set_major_locator(majorLocator)
			ax.xaxis.set_minor_locator(minorLocator)
			ax.yaxis.set_major_locator(majorLocator)
			ax.yaxis.set_minor_locator(minorLocator)
			plt.gca().invert_xaxis()  # make X axis increase to left, along +V2 direction
			ax.set_xlabel('Detector pixel (x)')
			ax.set_ylabel('Detector pixel (y)')
			plt.text(2000, 2100, '$\Delta$K<%4.1f mag' % (dmag_limit), color='black', fontsize=10)

			#draw V2 & V3 axes
			ax.arrow(-250, -250, 0, 1500, head_width=75, head_length=125, fc='k', ec='k')
			plt.text(-300, 500, 'V3', color='black', fontsize=14)
			ax.arrow(-250, -250, 1500, 0, head_width=75, head_length=125, fc='k', ec='k')
			plt.text(500, -400, 'V2', color='black', fontsize=14)


			#draw N & E axes
			dx_N = 200 * math.sin(-PA * D2R)
			dy_N = 200 * math.cos(-PA * D2R) 
			dx_E = 150 * math.sin((-PA+90) * D2R)
			dy_E = 150 * math.cos((-PA+90) * D2R) 
			ax.arrow(1600, 1600, dx_N, dy_N, head_width=50, head_length=80, fc='r', ec='r')
			ax.arrow(1600, 1600, dx_E, dy_E, head_width=50, head_length=80, fc='r', ec='r')
			plt.text(1600+dx_N, 1600+dy_N, 'N', color='black', fontsize=14)
			plt.text(1600+dx_E, 1600+dy_E, 'E', color='black', fontsize=14)

			ax.plot(xbox,ybox) # detector bounds
			#ax.plot([-142, -132, 2254, 2198], [2822, -174, -180, 2800], color='r') # field area that can be imaged on detector
			ax.plot([xmin1_xmax1y, xmin1_miny, xxmax1_miny, xxmax1_xmax1y], [yxmax1_minx, ymin_minx, ymin_xmax1x, yxmax1_xmax1x], color='r') # field area that can be imaged on detector
			ax.plot(X_F322W2, Y_field, color='g', marker='o')
			ax.plot([xmin1_F322W2, xmax_F322W2], [Y_field, Y_field], linewidth=2.0, color='g')
			ax.text(xmax_F322W2-200, Y_field+50, 'F322W2 spectrum', color='g')
			#ax.plot(X_F322W2+dx, Y_field+dy, color='black', linestyle='none', marker='o')
			for i in range(nlimit):
				#ax.plot(X_F322W2+dx[i], Y_field+dy[i],  color='black', linestyle='none', marker='o', markersize=msize[i])
				ax.plot(X_F322W2+dx[i], Y_field+dy[i],  color='black', linestyle='none', marker='o', markersize=(K0-dmag[i])*3)
				if (col[i] == 1 or col[i] == 3):
					ax.plot(X_F322W2+dx[i], Y_field+dy[i],  color='red', linestyle='none', marker='o', markersize=msize[i])
					ax.plot([xmin1_F322W2+dx[i], xmax_F322W2+dx[i]], [Y_field+dy[i], Y_field+dy[i]], linewidth=2.0, color='r')
				if (col[i] == 2 or col[i] == 3):
					ax.plot(X_F322W2+dx[i], Y_field+dy[i],  color='orange', linestyle='none', marker='o', markersize=msize[i])
					ax.plot([xmin1_F322W2m2+dx[i], xmax_F322W2m2+dx[i]], [Y_field+dy[i], Y_field+dy[i]], linewidth=2.0, color='orange')
				col[i] = 0 # reset collision value for next filter
			title = "%s_F322W2_PA%03.0f" % (targetName, PA)
		#     plt.ylim(-30,30)
			#fig.savefig(title)
	

		#NIRCam source offsets
		contaminationF322W2_order1 = np.zeros((nPA, nWaves))
		contaminationF322W2_order2 = np.zeros((nPA, nWaves))

		for kPA, PA in tqdm(enumerate(setPA), total=nPA):
			#NIRCam source offsets
			dx = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.cos((PA-0) * D2R) - ddec * 3600.0 / 0.065 * math.sin((PA-0) *D2R)
			dy = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.sin((PA-0) * D2R) + ddec * 3600.0 / 0.065 * math.cos((PA-0) *D2R)
			col_F322W2= np.zeros(nlimit)
	
			# Check for F322W2 spectral collisitons
			for i in range(nlimit):
				xmin1 = X_F322W2+dx[i] + (xmaxWL_F322W2 - wavelen_undev)/disp # m = 1 spectrum limits
				xmax1 = xmin1 + dx_F322W2
				x     = X_F322W2 + dx[i]
				y     = Y_field + dy[i]
				xmin2 = X_F322W2+dx[i] + (xmaxWL_F322W2 - wavelen_undev)*2.0 /disp # m = 2 spectrum limits
				xmax2 = xmin1 + dx_F322W2m2
		
				# First check if 1st order spectrum companion steps on 1st order spectrum of target 
				if dmag[i] != 0.0:
					if (((xmin1 > (xmin1_F322W2-X_buff) and xmin1 < (xmax_F322W2+X_buff)) or (xmax1 > (xmin1_F322W2-X_buff) and (xmax1 < (xmax_F322W2+X_buff)))) and y > (Y_field-Y_buff) and y < (Y_field+Y_buff)):
						# Next check to see that offending source is withn POM FOV and not the target itself:
						# if ((x+X_buff) > xmin1_xmax1y and (x-X_buff) < xmax1_miny and math.fabs(dx[i]) > 3.0 and math.fabs(dy[i]) > 3.0): # don't flag the target itself
						col[i] = 1
						# if dmag[i] is not 0:
						bgGaussian     = gaussian1D(center=y, width=nyquistWidth, height=normalHeight, offset=0)(y_gs+Y_field)
						contNow= sum(targetGaussian*bgGaussian)
						for kx in range(int(round(xmin1)), int(round(xmax1+1))):
						# for kx in range(int(round(xmin1_F322W2+dx[i])), int(round(xmax_F322W2))):
							# for kx in range(int(round(xmin1_F322W2)), int(round(xmax_F322W2))):
							if (kx >= 0) and (kx < xmax_F322W2):
								contaminationF322W2_order1[kPA,kx] += contNow
						# print("m = 1: dx {0:.1f}, dy {1:.1f}, dK = {2:.1f}".format(dx[i], dy[i], dmag[i]))
					# Now check if 2nd order spectrum of companion steps on 1st order spectrum of target
					if (((xmin2 > (xmin1_F322W2-X_buff) and xmin2 < (xmax_F322W2+X_buff)) or (xmax2 > (xmin1_F322W2-X_buff) and (xmax2 < (xmax_F322W2+X_buff)))) and y > (Y_field-Y_buff) and y < (Y_field+Y_buff)):
						# if ((x+X_buff) > xmin1_xmax1y and (x-X_buff) < xmax1_miny and math.fabs(dx[i]) > 3.0 and math.fabs(dy[i]) > 3.0): # don't flag the target itself
						col[i] = col[i] + 2 
						# if dmag[i] is not 0:
						bgGaussian     = gaussian1D(center=y, width=nyquistWidth, height=normalHeight, offset=0)(y_gs+Y_field)
						contNow= sum(targetGaussian*bgGaussian)
						for kx in range(int(round(xmin1)), int(round(xmax1+1))):
						# for kx in range(int(round(xmin1_F322W2+dx[i])), int(round(xmax_F322W2))):
							# for kx in range(int(round(xmin1_F322W2)), int(round(xmax_F322W2))):
							if (kx >= xmin1_F322W2m2) and (kx < 2048):
								contaminationF322W2_order2[kPA,kx] += contNow
						#print("m = 2: dx {0:.1f}, dy {1:.1f}, dK = {2:.1f}".format(dx[i], dy[i], dmag[i]))

	
		plt.figure(figsize=(20,10))
		#F322 1st order
		plt.subplot(121)
		plt.imshow(smooth_gaussconv(contaminationF322W2_order1, 0.5) + 1e-10, extent=[2.4,5.0,0,360], cmap=plt.cm.Reds, norm=LogNorm(), aspect='auto')
		plt.title("%s - NIRCam LWAR %s %s Order" % (targetName, 'F322W2', '1st'),fontsize=fontsize)
		plt.ylabel('Position Angle (deg)',fontsize=fontsize)
		plt.xlabel('Wavelength $(\mu m)$',fontsize=fontsize)
		plt.colorbar()

		#F322 2nd order
		plt.subplot(122)
		plt.imshow(smooth_gaussconv(contaminationF322W2_order2, 0.5) + 1e-10, extent=[2.4,5.0,0,360], cmap=plt.cm.Reds, norm=LogNorm(), aspect='auto')
		plt.title("%s - NIRCam LWAR %s 2nd Order" % (targetName, 'F322W2'),fontsize=fontsize)
		plt.ylabel('Position Angle (deg)',fontsize=fontsize)
		plt.xlabel('Wavelength $(\mu m)$',fontsize=fontsize)
		plt.colorbar()

		buff = io.BytesIO()
		plt.savefig(buff, format = 'png')
		buff.seek(0)
		figdata_png = base64.b64encode(buff.getvalue()).decode('ascii')
		
		
	elif ins == 'NIRCam_F444W':
		xbox = [-0, -0, 2040, 2040, -0]
		ybox = [-0, 2040, 2040, -0, -0]

		def plot_spectral_overlap_F444W(col, dra, ddec, Decd0, PA, msize, dmag, targetName):
	
			dx = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.cos((PA-0) * D2R) - ddec * 3600.0 / 0.065 * math.sin((PA-0) *D2R)
			dy = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.sin((PA-0) * D2R) + ddec * 3600.0 / 0.065 * math.cos((PA-0) *D2R)
	
			title = "%s, NCam LWAR F444W, V3PA =%3.0f$^\circ$" % (targetName, PA)
			fig, ax = plt.subplots(figsize=(6,6))
			ax.set_title(title)
			majorLocator = MultipleLocator(500)
			minorLocator = MultipleLocator(100)

			ax.set_xlim([-500, 2300])
			ax.set_ylim([-500, 2300])
			ax.xaxis.set_major_locator(majorLocator)
			ax.xaxis.set_minor_locator(minorLocator)
			ax.yaxis.set_major_locator(majorLocator)
			ax.yaxis.set_minor_locator(minorLocator)
			plt.gca().invert_xaxis()  # make X axis increase to left, along +V2 direction
			ax.set_xlabel('Detector pixel (x)')
			ax.set_ylabel('Detector pixel (y)')
			plt.text(2000, 2100, '$\Delta$K<%4.1f mag' % (dmag_limit), color='black', fontsize=10)

			#draw V2 & V3 axes
			ax.arrow(-250, -250, 0, 1500, head_width=75, head_length=125, fc='k', ec='k')
			plt.text(-300, 500, 'V3', color='black', fontsize=14)
			ax.arrow(-250, -250, 1500, 0, head_width=75, head_length=125, fc='k', ec='k')
			plt.text(500, -400, 'V2', color='black', fontsize=14)


			#draw N & E axes
			dx_N = 200 * math.sin(-PA * D2R)
			dy_N = 200 * math.cos(-PA * D2R) 
			dx_E = 150 * math.sin((-PA+90) * D2R)
			dy_E = 150 * math.cos((-PA+90) * D2R) 
			ax.arrow(1600, 1600, dx_N, dy_N, head_width=50, head_length=80, fc='r', ec='r')
			ax.arrow(1600, 1600, dx_E, dy_E, head_width=50, head_length=80, fc='r', ec='r')
			plt.text(1600+dx_N, 1600+dy_N, 'N', color='black', fontsize=14)
			plt.text(1600+dx_E, 1600+dy_E, 'E', color='black', fontsize=14)

			ax.plot(xbox,ybox) # detector bounds
			#ax.plot([-142, -132, 2254, 2198], [2822, -174, -180, 2800], color='r') # field area that can be imaged on detector
			ax.plot([xmin1_xmax1y, xmin1_miny, xxmax1_miny, xxmax1_xmax1y], [yxmax1_minx, ymin_minx, ymin_xmax1x, yxmax1_xmax1x], color='r') # field area that can be imaged on detector
			ax.plot(X_F444W, Y_field, color='g', marker='o')
			ax.plot([xmin1_F444W, xxmax1_F444W], [Y_field, Y_field], linewidth=2.0, color='g')
			ax.text(xxmax1_F444W-200, Y_field+50, 'F444W spectrum', color='g')
			for i in range(nlimit):
				ax.plot(X_F444W+dx[i], Y_field+dy[i],  color='black', linestyle='none', marker='o', markersize=msize[i])
				#ax.plot(X_F444W+dx[i], Y_field+dy[i],  color='black', linestyle='none', marker='o', markersize=(K0-dmag[i])*3)
				if (col[i] == 1):
					ax.plot(X_F444W+dx[i], Y_field+dy[i],  color='red', linestyle='none', marker='o', markersize=10)
					ax.plot([xmin1_F444W+dx[i], xxmax1_F444W+dx[i]], [Y_field+dy[i], Y_field+dy[i]], linewidth=2.0, color='r')
				col[i] = 0 # reset collision value for next filter

			title = "%s_F444W_PA%03.0f" % (targetName, PA)
			#fig.savefig(title)
			
		# Check for spectral collisitons
		contaminationF444W_order1 = np.zeros((nPA, nWaves))
		for kPA, PA in tqdm(enumerate(setPA), total=nPA):
			#NIRCam source offsets
			dx = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.cos((PA-0) * D2R) - ddec * 3600.0 / 0.065 * math.sin((PA-0) *D2R)
			dy = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.sin((PA-0) * D2R) + ddec * 3600.0 / 0.065 * math.cos((PA-0) *D2R)
			col= np.zeros(nlimit)

			for i in range(nlimit):
				xxmax1 = (MINWL_F444W -  wavelen_undev)/disp + X_F444W + dx[i]
				xmin1 = xxmax1_F444W - dx_F444W 
				x    = X_F444W + dx[i]
				y    = Y_field + dy[i]
				if dmag[i] != 0:
					if (((xmin1 > (xmin1_F444W-X_buff) and xmin1 < (xxmax1_F444W+X_buff)) or (xxmax1 > (xxmax1_F444W-X_buff) and xxmax1 < (xxmax1_F444W+X_buff))) and y > (Y_field-Y_buff) and y < (Y_field+Y_buff)):
						# Next check to see that offending source is withn POM FOV and not the target itself:
						# if ((x+X_buff) > xmin1_xmax1y and (x-X_buff) < xxmax1_miny and math.fabs(dx[i]) > 3.0 and math.fabs(dy[i]) > 3.0): # don't flag the target itself
							col[i] = 1
							# if dmag[i] is not 0:
							bgGaussian = gaussian1D(center=y, width=nyquistWidth, height=normalHeight, offset=0)(y_gs+Y_field)
							contNow    = sum(targetGaussian*bgGaussian)
							for kx in range(int(round(xmin1)), int(round(xxmax1+1))):
							# for kx in range(int(round(xmin1_F322W2+dx[i])), int(round(xmax_F322W2))):
								# for kx in range(int(round(xmin1_F322W2)), int(round(xmax_F322W2))):
								if (kx >= 0) and (kx < xmax_F322W2):
									contaminationF444W_order1[kPA,kx] += contNow
		
		plt.figure(figsize=(10,10))
		plt.imshow(smooth_gaussconv(contaminationF444W_order1, 0.5) + 1e-10, extent=[2.4, 5.0, 0, 360], cmap=plt.cm.Reds, norm=LogNorm(), aspect='auto')
		plt.title("%s, NIRCam LWAR %s 1st Order" % (targetName, 'F444W'),fontsize=fontsize)
		plt.ylabel('Position Angle (deg)',fontsize=fontsize)
		plt.xlabel('Wavelength $(\mu m)$',fontsize=fontsize)
		plt.colorbar()		
		
		buff = io.BytesIO()
		plt.savefig(buff, format = 'png')
		buff.seek(0)
		figdata_png = base64.b64encode(buff.getvalue()).decode('ascii')
		
	elif ins == 'MIRI':
		#MIRI: SLITLESSPRISM subarray is 72 x 416, with bottom at 1, 529
		#MIRI -Y axis PA is -5deg rel to V3 
		# https://jwst-docs.stsci.edu/display/JTI/MIRI+Overview

		X_LRS = 36 #WAG star position
		Y_LRS = 820 #WAG star position
		Y_xmax1 = 912 #WAG for 5 micron location
		dy_LRS = 370 # length of 5 - 12 micron LRS spectrum in pixels
		Y_MIN = Y_xmax1 - dy_LRS #WAG for 12 micron location


		MIRIPA = 5 # deg
		MIRIscale = 0.11 #arcsec / pxlMIRIscale = 0.11 #arcsec / pxl

		X_buff = 10 # Keep other spectra this far away from X_LRS (dist from center)
		Y_buff = 10 # Keep other spectra this far away from xmin1 and xxmax1 (dist from edge)		
		
		x_det = [-0, -0, 1024, 1024, -0]
		y_det = [-0, 1020, 1020, -0, -0]
		x_sub = [0, 0, 72, 72, -0]
		y_sub = [529, 529+416, 529+416, 529, 529]

		def plot_spectral_overlap_LRS(col, dx, dy, msize, dmag, targetName, PA):
			title = "%s, MIRI Slitless LRS, V3 PA =%3.0f$^\circ$" % (targetName, PA)
			fig, ax = plt.subplots(figsize=(6,6))
			ax.set_title(title)
			majorLocator = MultipleLocator(500)
			minorLocator = MultipleLocator(100)

			ax.set_xlim([-200, 1200])
			ax.set_ylim([-200, 1200])
			ax.xaxis.set_major_locator(majorLocator)
			ax.xaxis.set_minor_locator(minorLocator)
			ax.yaxis.set_major_locator(majorLocator)
			ax.yaxis.set_minor_locator(minorLocator)


			ax.set_xlabel('Detector pixel (x)')
			ax.set_ylabel('Detector pixel (y)')

			ax.plot(x_det,y_det)
			ax.plot(x_sub,y_sub)

			#draw V2 & V3 axes
			axlen = 800
			ax.arrow(1060, -100, axlen*math.sin(MIRIPA * D2R), axlen*math.cos(MIRIPA * D2R), head_width=40, head_length=60, fc='k', ec='k')
			plt.text(1030, 500, 'V3', color='black', fontsize=14)
			ax.arrow(1060, -100, -axlen*math.cos(MIRIPA * D2R), axlen*math.sin(MIRIPA * D2R), head_width=40, head_length=60, fc='k', ec='k')
			plt.text(500, -140, 'V2', color='black', fontsize=14)

			#draw N & E axes
			dx_N = 125 * math.sin((PA+MIRIPA) * D2R)
			dy_N = 125 * math.cos((PA+MIRIPA) * D2R) 
			dx_E = 100 * math.sin((PA+MIRIPA-90) * D2R)
			dy_E = 100 * math.cos((PA+MIRIPA-90) * D2R) 
			ax.arrow(800, 800, dx_N, dy_N, head_width=25, head_length=40, fc='r', ec='r')
			ax.arrow(800, 800, dx_E, dy_E, head_width=25, head_length=40, fc='r', ec='r')
			plt.text(800+dx_N, 800+dy_N, 'N', color='black', fontsize=14)
			plt.text(800+dx_E, 800+dy_E, 'E', color='black', fontsize=14)

			ax.plot(X_LRS, Y_LRS, color='g', marker='o')
			ax.plot([X_LRS, X_LRS], [Y_MIN, Y_xmax1], linewidth=2.0, color='g')
			for i in range(nlimit):
				ax.plot(X_LRS+dx[i], Y_LRS+dy[i],  color='black', linestyle='none', marker='o', markersize=msize[i])
				if (col[i] == 1):
					ax.plot(X_LRS+dx[i], Y_LRS+dy[i],  color='red', linestyle='none', marker='o', markersize=10)
					ax.plot([X_LRS+dx[i], X_LRS+dx[i]], [Y_MIN+dy[i], Y_xmax1+dy[i]], linewidth=2.0, color='r')
				col[i] = 0 # reset collision value for next filter 	

		contaminationLRS_order1 = np.zeros((nPA, nWaves))

		for kPA, PA in tqdm(enumerate(setPA), total=nPA):
			#MIRI source offsets
			dx  = -(dra * 3600.0 / MIRIscale) * math.cos(Decd0 * D2R) * math.cos((PA+MIRIPA) * D2R) + ddec * 3600.0 / MIRIscale * math.sin((PA-MIRIPA) *D2R)
			dy  = (dra * 3600.0 / MIRIscale) * math.cos(Decd0 * D2R) * math.sin((PA+MIRIPA) * D2R) + ddec * 3600.0 / MIRIscale * math.cos((PA-MIRIPA) *D2R)
			col = np.zeros(nlimit) # No spectral collision by default

			# Check for spectral collisitons
			for i in range(nlimit):
				yxmax1 = dy[i] + Y_LRS
				ymin = yxmax1 - dy_LRS
				x = X_LRS + dx[i]
				y = Y_LRS + dy[i]
				if dmag[i] != 0:
					if (((ymin > (Y_MIN-Y_buff) and ymin < (Y_xmax1+Y_buff)) or (yxmax1 > (Y_MIN-Y_buff) and (yxmax1 < (Y_xmax1+Y_buff)))) and x > (X_LRS-X_buff) and x < (X_LRS+X_buff)):
						# Next check to see that offending source is withn POM FOV and not the target itself:
						# if ((x+X_buff) > xmin1_xmax1y and (x-X_buff) < xxmax1_miny and math.fabs(dx[i]) > 1 and math.fabs(dy[i]) > 1): # don't flag the target itself
							col[i] = 1         
							# if dmag[i] is not 0:
							bgGaussian  = gaussian1D(center=x, width=nyquistWidth, height=normalHeight, offset=0)(x_gs+X_LRS)
							contNow     = sum(targetGaussian*bgGaussian)
							for kx in range(int(round(ymin)), int(round(yxmax1+1))):
								if (kx >= 0) and (kx < 2048):
									contaminationLRS_order1[kPA,kx] += contNow
							# print("dx {0:.1f}, dy {1:.1f}, dK = {2:.1f}".format(dx[i], dy[i], dmag[i]))		
							

		plt.figure(figsize=(10,10))
		plt.imshow(smooth_gaussconv(contaminationLRS_order1, 0.5) + 1e-10, extent=[5,12,0,360], cmap=plt.cm.Reds, norm=LogNorm(), aspect='auto')
		plt.title("%s - MIRI %s" % (targetName, 'LRS'),fontsize=fontsize)
		plt.ylabel('Position Angle (deg)',fontsize=fontsize)
		plt.xlabel('Wavelength $(\mu m)$',fontsize=fontsize)
		plt.colorbar()
		
		buff = io.BytesIO()
		plt.savefig(buff, format = 'png')
		buff.seek(0)
		figdata_png = base64.b64encode(buff.getvalue()).decode('ascii')				

	else:
		print('Instrument not recognized')
		print('Your options are:')
		print('NIRCam_F322W2')
		print('NIRCam_F444W')
		print('MIRI')
		
	return figdata_png	
		
if __name__ == '__main__':
	tname = argv[1]
	inst  = 'ncam-f322' if len(argv)<3 else argv[2]
	ncam_miri(targetName = tname, ins = inst)		

		

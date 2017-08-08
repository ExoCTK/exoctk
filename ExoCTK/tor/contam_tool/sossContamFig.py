import ExoCTK
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator,MaxNLocator
from . import visibilityPA as vis
import os, io, base64

#Hack to override default hatch linewidth for PDF
import matplotlib
import six
from matplotlib.path import Path
from matplotlib.backends.backend_pdf import Name, Op
from matplotlib.transforms import Affine2D
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
# from flask import make_response


def setCustomHatchWidth(customWidth):
	def _writeHatches(self):
		hatchDict = dict()
		sidelen = 72.0
		for hatch_style, name in six.iteritems(self.hatchPatterns):
			ob = self.reserveObject('hatch pattern')
			hatchDict[name] = ob
			res = {'Procsets':
				   [Name(x) for x in "PDF Text ImageB ImageC ImageI".split()]}
			self.beginStream(
				ob.id, None,
				{'Type': Name('Pattern'),
				 'PatternType': 1, 'PaintType': 1, 'TilingType': 1,
				 'BBox': [0, 0, sidelen, sidelen],
				 'XStep': sidelen, 'YStep': sidelen,
				 'Resources': res})

			# lst is a tuple of stroke color, fill color,
			# number of - lines, number of / lines,
			# number of | lines, number of \ lines
			rgb = hatch_style[0]
			self.output(rgb[0], rgb[1], rgb[2], Op.setrgb_stroke)
			if hatch_style[1] is not None:
				rgb = hatch_style[1]
				self.output(rgb[0], rgb[1], rgb[2], Op.setrgb_nonstroke,
							0, 0, sidelen, sidelen, Op.rectangle,
							Op.fill)

			self.output(customWidth, Op.setlinewidth) ###the new width###

			# TODO: We could make this dpi-dependent, but that would be
			# an API change
			self.output(*self.pathOperations(
				Path.hatch(hatch_style[2]),
				Affine2D().scale(sidelen),
				simplify=False))
			self.output(Op.stroke)

			self.endStream()
		self.writeObject(self.hatchObject, hatchDict)

	matplotlib.backends.backend_pdf.PdfFile.writeHatches = _writeHatches

def cmap_discretize(cmap, N):
	"""Return a discrete colormap from the continuous colormap cmap.

		cmap: colormap instance, eg. cm.jet. 
		N: number of colors.

	Example
		x = resize(arange(100), (5,100))
		djet = cmap_discretize(cm.jet, 5)
		imshow(x, cmap=djet)
	"""

	if type(cmap) == str:
		cmap = matplotlib.cm.get_cmap(cmap)
	colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
	colors_rgba = cmap(colors_i)
	indices = np.linspace(0, 1., N+1)
	cdict = {}
	for ki,key in enumerate(('red','green','blue')):
		cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in range(N+1) ]
	# Return colormap object.
	return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

def contam(ra, dec, targetName, cubeName, pamin = 0, pamax = 360):

	if (type(pamin) == str) | (type(pamax) == str):
		try:
			pamin, pamax = int(pamin), int(pamax)
		except:
			pamin, pamax = 0, 360
	
	paRange = [pamin, pamax]
	goodPA, badPA, _ =vis.calc_vis(ra, dec, targetName=targetName)
	
	plotPAmin,plotPAmax=paRange
	suffix='_PA'+str(plotPAmin)+'-'+str(plotPAmax)

	#start calculations
	ypix,lamO1,lamO2=np.loadtxt(os.path.join(os.path.dirname(ExoCTK.__file__),'data/tor/lambda_order1-2.txt'),unpack=True)

#    hdu=fits.open('cube_'+target+'.fits')
	hdu=fits.open(cubeName)
	trace2dO1=hdu[0].data[0,:,:] #order 1
	trace2dO2=hdu[0].data[1,:,:] #order 2
	cube=hdu[0].data[2:,:,:] #all the angles
	hdu.close()

	ny=trace2dO1.shape[0]
	nPA=cube.shape[0]
	dPA=360//nPA
	PA=np.arange(nPA)*dPA

	contamO1=np.zeros([ny,nPA])
	contamO2=np.zeros([ny,nPA])
	for y in np.arange(ny):
		i=np.argmax(trace2dO1[y,:])
		tr=trace2dO1[y,i-20:i+41]
		w=tr/np.sum(tr**2)
		ww=np.tile(w,nPA).reshape([nPA,tr.size])
		contamO1[y,:]=np.sum(cube[:,y,i-20:i+41]*ww,axis=1)

		if lamO2[y]<0.6: continue
		i=np.argmax(trace2dO2[y,:])
		tr=trace2dO2[y,i-20:i+41]
		w=tr/np.sum(tr**2)
		ww=np.tile(w,nPA).reshape([nPA,tr.size])
		contamO2[y,:]=np.sum(cube[:,y,i-20:i+41]*ww,axis=1)

	####make the plot####
	fig=plt.figure(figsize=(10,5))

	#main panel order 1
	ax1=plt.axes([0.425,0.11,0.4,0.7])
	ax1.set_xticks(np.arange(1,3,0.25))
	ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
	plt.ylim(plotPAmin-0.5*dPA,plotPAmax+0.5*dPA)
	if plotPAmax-plotPAmin>=200:
		minTicksMult=5
	elif plotPAmax-plotPAmin>=90:
		minTicksMult=2
	else:
		minTicksMult=1
	ax1.yaxis.set_minor_locator(MultipleLocator(minTicksMult))
	#ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
	plt.grid()
	plt.xlabel('Wavelength (microns)',fontsize='small')
	plt.ylabel('Position Angle (degrees)',fontsize='small')
	#contamCmap=cmap_discretize('spectral_r',8)
	contamCmap=cmap_discretize('CMRmap_r',8)
	contamCmap.set_under('w')
	plt.imshow(np.log10(np.clip(contamO1.T,1.e-10,1.)),extent=(lamO1.min(),lamO1.max(),PA.min()-0.5*dPA,PA.max()+0.5*dPA),
		vmin=-4,vmax=0,aspect='auto',origin='lower',interpolation='nearest',cmap=contamCmap)

	for p in badPA:
		plt.fill_between([lamO1.min(),lamO1.max()],p[0],p[1],hatch='xx',facecolors='none',edgecolor='k',lw=0)
	setCustomHatchWidth(0.4)

	#main panel order 2
	ax2=plt.axes([0.18,0.11,0.18,0.7])
	ax2.set_xticks(np.arange(0.6,1.4,0.25))
	ax2.xaxis.set_minor_locator(MultipleLocator(0.05))
	plt.ylim(plotPAmin-0.5*dPA,plotPAmax+0.5*dPA)
	ax2.yaxis.set_minor_locator(MultipleLocator(minTicksMult))
	ax2.set_yticklabels([])
	plt.grid()
	plt.xlabel('Wavelength (microns)',fontsize='small')
	yminO2=np.argmin(np.abs(lamO2-0.6))
	plt.imshow(np.log10(np.clip(contamO2[yminO2:,:].T,1.e-10,1.)),extent=(lamO2[yminO2],lamO2.max(),PA.min()-0.5*dPA,PA.max()+0.5*dPA),
		vmin=-4,vmax=0,aspect='auto',origin='lower',interpolation='nearest',cmap=contamCmap)

	for p in badPA:
		plt.fill_between([lamO2[yminO2],lamO2.max()],p[0],p[1],hatch='xx',facecolors='none',edgecolor='k',lw=0)


	#scale bar associated with main panel orders 1 and 2
	ax4=plt.axes([0.85,0.9,0.11,0.04])
	cbar=np.tile(np.arange(-4,0.01,0.1),2).reshape([2,-1])
	plt.imshow(cbar,vmin=-4,vmax=0,aspect='auto',origin='lower',interpolation='nearest',
		extent=(-4,0,0,1),cmap=contamCmap)
	ax4.set_yticks([])
	plt.xticks(np.arange(-4,0.1),fontsize='small')
	plt.xlabel('log(contam. level)',fontsize='small')

	#right panel, fraction of contaminated wavelength channels vs PA, order 1
	ax1a=plt.axes([0.84,0.11,0.13,0.7])
	plt.xlim(0,95)
	ax1a.xaxis.set_minor_locator(MultipleLocator(5))
	plt.ylim(plotPAmin-0.5*dPA,plotPAmax+0.5*dPA)
	ax1a.yaxis.set_minor_locator(MultipleLocator(minTicksMult))
	ax1a.set_yticklabels([])
	plt.xlabel('\% channels contam. \n above threshold',fontsize='small')
	plt.grid()
	plt.plot(100*np.sum(contamO1 >= 0.001,axis=0)/ny,PA-dPA/2,ls='steps',linewidth=1.5,label='>0.001')
	plt.plot(100*np.sum(contamO1 >= 0.01,axis=0)/ny,PA-dPA/2,ls='steps',linewidth=1.5,label='>0.01')
	#plt.plot(np.sum(contam >= 0.05,axis=0)/ny,PA,ls='steps',linewidth=1.5)
	for p in badPA:
		plt.fill_between([0,100],p[0],p[1],hatch='xx',facecolors='none',edgecolor='k',lw=0)
	leg=plt.legend(fontsize='xx-small',loc='best',framealpha=0.75,
		handlelength=1.5,handletextpad=0.2,labelspacing=0.3)
	leg.get_frame().set_linewidth(0.0)

	#left panel, fraction of contaminated wavelength channels vs PA, order 2
	ax2a=plt.axes([0.03,0.11,0.13,0.7])
	plt.xlim(0,95)
	ax2a.xaxis.set_minor_locator(MultipleLocator(5))
	plt.ylim(plotPAmin-0.5*dPA,plotPAmax+0.5*dPA)
	ax2a.yaxis.set_minor_locator(MultipleLocator(minTicksMult))
	ax2a.set_yticklabels([])
	plt.xlabel('\% channels contam. \n above threshold',fontsize='small') 
	plt.grid()
	plt.plot(100*np.sum(contamO2 >= 0.001,axis=0)/ny,PA-dPA/2,ls='steps',linewidth=1.5,label='>0.001')
	plt.plot(100*np.sum(contamO2 >= 0.01,axis=0)/ny,PA-dPA/2,ls='steps',linewidth=1.5,label='>0.01')
	#plt.plot(np.sum(contam >= 0.05,axis=0)/ny,PA,ls='steps',linewidth=1.5)
	for p in badPA:
		plt.fill_between([0,100],p[0],p[1],hatch='xx',facecolors='none',edgecolor='k',lw=0)


	#top panel, Annotation of spectral features, order 1
	ax1b=plt.axes([0.425,0.83,0.4,0.14])
	plt.xlim(lamO1.min(),lamO1.max())
	plt.ylim(0.05,0.3)
	ax1b.set_xticks(np.arange(1,3,0.25))
	ax1b.xaxis.set_minor_locator(MultipleLocator(0.05))
	ax1b.set_yticks([])
	ax1b.set_xticklabels([])

	plt.grid(axis='x')
	y=np.array([0.,0.])
	y1=0.07
	y2=0.12
	y3=0.17
	y4=0.23

	l=np.array([0.89,0.99])
	plt.plot(l,y+y1,color='k',linewidth=1.5)
	plt.text(l.mean(),y1,r'H$_2$O',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.09,1.2])
	plt.plot(l,y+y1,color='k',linewidth=1.5)
	plt.text(l.mean(),y1,r'H$_2$O',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.1,1.24])
	plt.plot(l,y+y2,color='k',linewidth=1.5)
	plt.text(l.mean(),y2,r'CH$_4$',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.3,1.51])
	plt.plot(l,y+y1,color='k',linewidth=1.5)
	plt.text(l.mean(),y1,r'H$_2$O',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.6,1.8])
	plt.plot(l,y+y2,color='k',linewidth=1.5)
	plt.text(l.mean(),y2,r'CH$_4$',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.75,2.05])
	plt.plot(l,y+y1,color='k',linewidth=1.5)
	plt.text(l.mean(),y1,r'H$_2$O',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([2.3,lamO1.max()])
	plt.plot(l,y+y1,color='k',linewidth=1.5)
	plt.text(l.mean(),y1,r'H$_2$O',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([2.15,2.5])
	plt.plot(l,y+y2,color='k',linewidth=1.5)
	plt.text(l.mean(),y2,r'CH$_4$',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.1692,1.1778])
	plt.vlines(l,y3,y3+0.02,color='k')
	plt.text(l.mean(),y3+0.02,'K',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.2437,1.2529])
	plt.vlines(l,y3,y3+0.02,color='k')
	plt.text(l.mean(),y3+0.02,'K',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.5168])
	plt.vlines(l,y3,y3+0.02,color='k')
	plt.text(l.mean(),y3+0.02,'K',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.1384,1.1409])
	plt.vlines(l,y4,y4+0.02,color='k')
	plt.text(l.mean(),y4+0.02,'Na',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.2682])
	plt.vlines(l,y4,y4+0.02,color='k')
	plt.text(l.mean(),y4+0.02,'Na',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([2.2063,2.2090])
	plt.vlines(l,y4,y4+0.02,color='k')
	plt.text(l.mean(),y4+0.02,'Na',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([2.2935,2.3227,2.3525,2.3830,2.4141])
	plt.vlines(l,y3,y3+0.02,color='k')
	plt.plot(l[[0,-1]],y+y3+0.02,color='k',linewidth=1)
	plt.text(l[[0,-1]].mean(),y3+0.02,'CO',ha='center',va='bottom',fontsize='xx-small')

	#top panel, Annotation of spectral features, order 2
	ax2b=plt.axes([0.18,0.83,0.18,0.14])
	plt.xlim(lamO2[yminO2],lamO2.max())
	plt.ylim(0.05,0.3)
	ax2b.set_xticks(np.arange(0.6,1.4,0.25))
	ax2b.xaxis.set_minor_locator(MultipleLocator(0.05))
	ax2b.set_yticks([])
	ax2b.set_xticklabels([])

	plt.grid(axis='x')

	l=np.array([0.89,0.99])
	plt.plot(l,y+y1,color='k',linewidth=1.5)
	plt.text(l.mean(),y1,r'H$_2$O',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.09,1.2])
	plt.plot(l,y+y1,color='k',linewidth=1.5)
	plt.text(l.mean(),y1,r'H$_2$O',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.1,1.24])
	plt.plot(l,y+y2,color='k',linewidth=1.5)
	plt.text(l.mean(),y2,r'CH$_4$',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.3,lamO2.max()])
	plt.plot(l,y+y1,color='k',linewidth=1.5)
	plt.text(l.mean(),y1,r'H$_2$O',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([0.7665,0.7699])
	plt.vlines(l,y3,y3+0.02,color='k')
	plt.text(l.mean(),y3+0.02,'K',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.1692,1.1778])
	plt.vlines(l,y3,y3+0.02,color='k')
	plt.text(l.mean(),y3+0.02,'K',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.2437,1.2529])
	plt.vlines(l,y3,y3+0.02,color='k')
	plt.text(l.mean(),y3+0.02,'K',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.1384,1.1409])
	plt.vlines(l,y4,y4+0.02,color='k')
	plt.text(l.mean(),y4+0.02,'Na',ha='center',va='bottom',fontsize='xx-small')

	l=np.array([1.2682])
	plt.vlines(l,y4,y4+0.02,color='k')
	plt.text(l.mean(),y4+0.02,'Na',ha='center',va='bottom',fontsize='xx-small')

	#add target name
	plt.figtext(0.04,0.9,targetName,fontsize='large')
	suffix="_PA{}-{}".format(*paRange)
	buff = io.BytesIO()
	plt.savefig(buff, format = 'png')
	buff.seek(0)
	figdata_png = base64.b64encode(buff.getvalue()).decode('ascii')
	return figdata_png
		

if __name__ == "__main__":
    #arguments RA & DEC, conversion to radians
    argv = sys.argv

    ra=argv[1]
    dec=argv[2]
    cubeName=argv[3]

    pamin=0 if len(argv)<5 else int(argv[4])
    pamax=360 if len(argv)<6 else int(argv[5])

    targetName=None if len(argv)<7 else argv[6]
    save=False if len(argv)<7 else True #if name provided -> save
    tmpDir="." if len(argv)<8 else argv[7]
    os.makedirs(tmpDir, exist_ok=True)
#     pdb.set_trace()
    goodPA, badPA, _ =calc_vis(ra, dec, targetName=targetName)

    contam(cubeName,targetName=targetName,paRange=[pamin,pamax],badPA=badPA)
    


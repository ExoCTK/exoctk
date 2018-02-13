import numpy as np
from astropy.io import fits
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator,MaxNLocator
from . import visibilityPA as vpa
import os
import matplotlib
import six
from matplotlib.path import Path
from matplotlib.backends.backend_pdf import Name, Op
from matplotlib.transforms import Affine2D
import pkg_resources
import base64
import io
from bokeh.io import gridplot, show
from bokeh.plotting import figure
from bokeh.models import Range1d, LinearColorMapper, LogColorMapper, LogTicker, ColorBar, Label
from bokeh.palettes import inferno

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

def contam(cube, targetName='noName', paRange=[0,360], badPA=[], tmpDir="", fig='', to_html=True):
    """
    Generate the contamination plot
    
    Parameters
    ----------
    cube: array-like, str
        The data cube or FITS filename containing the data
    targetName: str
        The name of the target
    paRange: sequence
        The position angle range to consider
    badPA: sequence
        Position angles to exclude
    tmpDir: str
        A directory to write the files to
    fig: matplotlib.figure, bokeh.figure
        A figure to add the plots to
    to_html: bool
        Return the image as bytes for HTML
    
    Returns
    -------
    fig
        The populated matplotlib or bokeh plot
    """
    # Get data from FITS file
    if isinstance(cube, str):
        # hdu = fits.open('cube_'+target+'.fits')
        hdu = fits.open(cubeName)
        cube = hdu[0].data
        hdu.close()
    
    trace2dO1 = cube[0,:,:] #order 1
    trace2dO2 = cube[1,:,:] #order 2
    cube = cube[2:,:,:] #all the angles
    
    plotPAmin,plotPAmax=paRange
    suffix='_PA'+str(plotPAmin)+'-'+str(plotPAmax)

    #start calculations
    lam_file = pkg_resources.resource_filename('ExoCTK', 'data/tor/lambda_order1-2.txt')
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    ny = trace2dO1.shape[0]
    nPA = cube.shape[0]
    dPA = 360//nPA
    PA = np.arange(nPA)*dPA

    contamO1 = np.zeros([ny,nPA])
    contamO2 = np.zeros([ny,nPA])
    for y in np.arange(ny):
        i = np.argmax(trace2dO1[y,:])
        tr = trace2dO1[y,i-20:i+41]
        w = tr/np.sum(tr**2)
        ww = np.tile(w,nPA).reshape([nPA,tr.size])
        contamO1[y,:] = np.sum(cube[:,y,i-20:i+41]*ww,axis = 1)

        if lamO2[y]<0.6: continue
        i = np.argmax(trace2dO2[y,:])
        tr = trace2dO2[y,i-20:i+41]
        w = tr/np.sum(tr**2)
        ww = np.tile(w,nPA).reshape([nPA,tr.size])
        contamO2[y,:] = np.sum(cube[:,y,i-20:i+41]*ww,axis = 1)
        
    if fig=='':
        fig = plt.figure(figsize=(10,5))
        
    # Make the default plot matplotlib
    if isinstance(fig, matplotlib.figure.Figure):
        
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
        # contamCmap=cmap_discretize('spectral_',8)
        contamCmap=cmap_discretize('CMRmap_r',8)
        contamCmap.set_under('w')
        plt.imshow(np.log10(np.clip(contamO1.T,1.e-10,1.)),extent=(lamO1.min(),lamO1.max(),PA.min()-0.5*dPA,PA.max()+0.5*dPA),
            vmin=-4,vmax=0,aspect='auto',origin='lowe',interpolation='nearest',cmap=contamCmap)

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
            vmin=-4,vmax=0,aspect='auto',origin='lowe',interpolation='nearest',cmap=contamCmap)

        for p in badPA:
            plt.fill_between([lamO2[yminO2],lamO2.max()],p[0],p[1],hatch='xx',facecolors='none',edgecolor='k',lw=0)


        #scale bar associated with main panel orders 1 and 2
        ax4=plt.axes([0.85,0.9,0.11,0.04])
        cbar=np.tile(np.arange(-4,0.01,0.1),2).reshape([2,-1])
        plt.imshow(cbar,vmin=-4,vmax=0,aspect='auto',origin='lowe',interpolation='nearest',
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
        plt.xlabel('Pct. channels contam. \n above threshold',fontsize='small')
        plt.grid()
        plt.plot(100*np.sum(contamO1 >= 0.001,axis=0)/ny,PA-dPA/2,ls='steps',linewidth=1.5,label='> 0.001')
        plt.plot(100*np.sum(contamO1 >= 0.01,axis=0)/ny,PA-dPA/2,ls='steps',linewidth=1.5,label='> 0.01')
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
        plt.xlabel('Pct. channels contam. \n above threshold',fontsize='small')
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
        plt.text(l.mean(),y1,'H2O',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([1.09,1.2])
        plt.plot(l,y+y1,color='k',linewidth=1.5)
        plt.text(l.mean(),y1,'H2O',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([1.1,1.24])
        plt.plot(l,y+y2,color='k',linewidth=1.5)
        plt.text(l.mean(),y2,'CH4',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([1.3,1.51])
        plt.plot(l,y+y1,color='k',linewidth=1.5)
        plt.text(l.mean(),y1,'H2O',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([1.6,1.8])
        plt.plot(l,y+y2,color='k',linewidth=1.5)
        plt.text(l.mean(),y2,'CH4',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([1.75,2.05])
        plt.plot(l,y+y1,color='k',linewidth=1.5)
        plt.text(l.mean(),y1,'H2O',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([2.3,lamO1.max()])
        plt.plot(l,y+y1,color='k',linewidth=1.5)
        plt.text(l.mean(),y1,'H2O',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([2.15,2.5])
        plt.plot(l,y+y2,color='k',linewidth=1.5)
        plt.text(l.mean(),y2,'CH4',ha='center',va='bottom',fontsize='xx-small')

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
        plt.text(l.mean(),y1,'H2O',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([1.09,1.2])
        plt.plot(l,y+y1,color='k',linewidth=1.5)
        plt.text(l.mean(),y1,'H2O',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([1.1,1.24])
        plt.plot(l,y+y2,color='k',linewidth=1.5)
        plt.text(l.mean(),y2,'CH4',ha='center',va='bottom',fontsize='xx-small')

        l=np.array([1.3,lamO2.max()])
        plt.plot(l,y+y1,color='k',linewidth=1.5)
        plt.text(l.mean(),y1,'H2O',ha='center',va='bottom',fontsize='xx-small')

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
        # suffix="_PA{}-{}".format(*paRange)
        # plt.savefig(tmpDir+'/contamination-'+targetName+suffix+'.pdf',bbox_inches='tight')
        # plt.savefig(tmpDir+'/contamination-'+targetName+suffix+'.png',bbox_inches='tight')
        # plt.show()
        
        if to_html:
            
            #save plot to memory in form of bytes
            buff = io.BytesIO()
            plt.savefig(buff, format = 'png')
            buff.seek(0)
            figdata_png = base64.b64encode(buff.getvalue()).decode('ascii')
        
            return figdata_png
        
    # Otherwise, it's a Bokeh plot
    else:
        
        TOOLS = 'crosshair,reset,hover,save'
        
        # ==================================================================================================
        # Order 1 ==========================================================================================
        # ==================================================================================================
        
        # Line list
        s1 = figure(tools=TOOLS, width=500, plot_height=100, title=None)
        y=np.array([0.,0.])
        y1=0.07
        y2=0.12
        y3=0.17
        y4=0.23

        l=np.array([0.89,0.99])
        s1.line(l,y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data', text='H2O', render_mode='css', text_font_size='8pt'))

        l=np.array([1.09,1.2])
        s1.line(l,y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data', text='H2O', render_mode='css', text_font_size='8pt'))

        l=np.array([1.1,1.24])
        s1.line(l,y+y2, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data', text='CH4', render_mode='css', text_font_size='8pt'))

        l=np.array([1.3,1.51])
        s1.line(l,y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data', text='H2O', render_mode='css', text_font_size='8pt'))

        l=np.array([1.6,1.8])
        s1.line(l,y+y2, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data', text='CH4', render_mode='css', text_font_size='8pt'))

        l=np.array([1.75,2.05])
        s1.line(l,y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data', text='H2O', render_mode='css', text_font_size='8pt'))

        l=np.array([2.3,lamO1.max()])
        s1.line(l,y+y1, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data', text='H2O', render_mode='css', text_font_size='8pt'))

        l=np.array([2.15,2.5])
        s1.line(l,y+y2, line_color='black', line_width=1.5)
        s1.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data', text='CH4', render_mode='css', text_font_size='8pt'))

        l=np.array([1.1692,1.1778])
        s1.line(l[0], [y3,y3+0.02], line_color='black')
        s1.line(l[1], [y3,y3+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data', y_units='data', text='K', render_mode='css', text_font_size='8pt'))

        l=np.array([1.2437,1.2529])
        s1.line(l[0], [y3,y3+0.02], line_color='black')
        s1.line(l[1], [y3,y3+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data', y_units='data', text='K', render_mode='css', text_font_size='8pt'))
       
        l=np.array([1.5168])
        s1.line(l[0], [y3,y3+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data', y_units='data', text='K', render_mode='css', text_font_size='8pt'))
        #
        l=np.array([1.1384,1.1409])
        s1.line(l[0], [y4,y4+0.02], line_color='black')
        s1.line(l[1], [y4,y4+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data', y_units='data', text='Na', render_mode='css', text_font_size='8pt'))
       
        l=np.array([1.2682])
        s1.line(l[0], [y4,y4+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data', y_units='data', text='Na', render_mode='css', text_font_size='8pt'))
       
        l=np.array([2.2063,2.2090])
        s1.line(l[0], [y4,y4+0.02], line_color='black')
        s1.line(l[1], [y4,y4+0.02], line_color='black')
        s1.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data', y_units='data', text='Na', render_mode='css', text_font_size='8pt'))
       
        l=np.array([2.2935,2.3227,2.3525,2.3830,2.4141])
        s1.line(l[0], [y3,y3+0.02], line_color='black')
        s1.line(l[1], [y3,y3+0.02], line_color='black')
        s1.line(l[2], [y3,y3+0.02], line_color='black')
        s1.line(l[3], [y3,y3+0.02], line_color='black')
        s1.line(l[4], [y3,y3+0.02], line_color='black')
        s1.line(l[[0,-1]],y+y3+0.02, line_color='black', line_width=1)
        s1.add_layout(Label(x=l[[0,-1]].mean(), y=y3+0.02, x_units='data', y_units='data', text='CO', render_mode='css', text_font_size='8pt'))
        
        s1.xaxis.major_label_text_font_size = '0pt'
        s1.yaxis.major_label_text_font_size = '0pt'

        # Contam plot
        xlim0, xlim1, ylim0, ylim1 = lamO1.min(), lamO1.max(), PA.min()-0.5*dPA, PA.max()+0.5*dPA
        color_mapper = LinearColorMapper(palette=inferno(8)[::-1], low=-4, high=1)
        color_mapper.low_color = 'white'
        color_mapper.high_color = 'black'
        s2 = figure(tools=TOOLS, width=500, height=500, title=None, x_range=Range1d(xlim0, xlim1), y_range=Range1d(ylim0, ylim1))
        fig_data = np.log10(np.clip(contamO1.T,1.e-10,1.))
        s2.image([fig_data], x=xlim0, y=ylim0, dw=xlim1-xlim0, dh=ylim1-ylim0, color_mapper=color_mapper)
        # color_bar = ColorBar(color_mapper=color_mapper, location=(0,0))
        s2.xaxis.axis_label = 'Wavelength (um)'
        s2.yaxis.axis_label = 'Position Angle (degrees)'
        # s2.add_layout(color_bar, 'below')

        # Line plot
        s3 = figure(tools=TOOLS, width=150, height=500, title=None)
        s3.line(100*np.sum(contamO1 >= 0.001,axis=0)/ny, PA-dPA/2, line_color='blue', legend='> 0.001')
        s3.line(100*np.sum(contamO1 >= 0.01,axis=0)/ny, PA-dPA/2, line_color='green', legend='> 0.01')
        s3.xaxis.axis_label = '% channels contam.'
        s3.yaxis.major_label_text_font_size = '0pt'
        s3.x_range = Range1d(0, 100)
        s3.y_range = Range1d(0, 360)

        # Add bad PAs
        for ybad0,ybad1 in badPA:
            s2.patch([xlim0,xlim1,xlim1,xlim0], [ybad1,ybad1,ybad0,ybad0], color='white', alpha=0.7)
            s3.patch([0,100,100,0], [ybad1,ybad1,ybad0,ybad0], color='white', alpha=0.7)

        # ==================================================================================================
        # Order 2 ==========================================================================================
        # ==================================================================================================
        
        # Line list
        s4 = figure(tools=TOOLS, width=250, plot_height=100, title=None)
        l=np.array([0.89,0.99])
        s4.line(l,y+y1, line_color='black', line_width=1.5)
        s4.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data', text='H2O', render_mode='css', text_font_size='8pt'))

        l=np.array([1.09,1.2])
        s4.line(l,y+y1, line_color='black', line_width=1.5)
        s4.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data', text='H2O', render_mode='css', text_font_size='8pt'))

        l=np.array([1.1,1.24])
        s4.line(l,y+y2, line_color='black', line_width=1.5)
        s4.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data', text='CH4', render_mode='css', text_font_size='8pt'))

        l=np.array([1.3,lamO2.max()])
        s4.line(l,y+y1, line_color='black', line_width=1.5)
        s4.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data', text='H2O', render_mode='css', text_font_size='8pt'))

        l=np.array([0.7665,0.7699])
        s4.line(l[0], [y3,y3+0.02], line_color='black')
        s4.line(l[1], [y3,y3+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data', y_units='data', text='K', render_mode='css', text_font_size='8pt'))

        l=np.array([1.1692,1.1778])
        s4.line(l[0], [y3,y3+0.02], line_color='black')
        s4.line(l[1], [y3,y3+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data', y_units='data', text='K', render_mode='css', text_font_size='8pt'))

        l=np.array([1.2437,1.2529])
        s4.line(l[0], [y3,y3+0.02], line_color='black')
        s4.line(l[1], [y3,y3+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data', y_units='data', text='K', render_mode='css', text_font_size='8pt'))

        l=np.array([1.1384,1.1409])
        s4.line(l[0], [y4,y4+0.02], line_color='black')
        s4.line(l[1], [y4,y4+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data', y_units='data', text='Na', render_mode='css', text_font_size='8pt'))

        l=np.array([1.2682])
        s4.line(l[0], [y4,y4+0.02], line_color='black')
        s4.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data', y_units='data', text='Na', render_mode='css', text_font_size='8pt'))
        
        s4.xaxis.major_label_text_font_size = '0pt'
        s4.yaxis.major_label_text_font_size = '0pt'
        
        # Contam plot
        xlim0, xlim1, ylim0, ylim1 = lamO2.min(), lamO2.max(), PA.min()-0.5*dPA, PA.max()+0.5*dPA
        xlim0 = 0.614
        s5 = figure(tools=TOOLS, width=250, height=500, title=None, x_range=Range1d(xlim0, xlim1), y_range=Range1d(ylim0, ylim1))
        fig_data = np.log10(np.clip(contamO2.T,1.e-10,1.))[:,300:]
        s5.image([fig_data], x=xlim0, y=ylim0, dw=xlim1-xlim0, dh=ylim1-ylim0, color_mapper=color_mapper)
        s5.yaxis.major_label_text_font_size = '0pt'
        s5.xaxis.axis_label = 'Wavelength (um)'

        # Line plot
        s6 = figure(tools=TOOLS, width=150, height=500, title=None)
        s6.line(100*np.sum(contamO2 >= 0.001,axis=0)/ny, PA-dPA/2, line_color='blue', legend='> 0.001')
        s6.line(100*np.sum(contamO2 >= 0.01,axis=0)/ny, PA-dPA/2, line_color='green', legend='> 0.01')
        s6.xaxis.axis_label = '% channels contam.'
        s6.yaxis.major_label_text_font_size = '0pt'
        s6.x_range = Range1d(100, 0)
        s6.y_range = Range1d(0, 360)
        
        # Dummy plots for nice spacing
        s0 = figure(tools=TOOLS, width=150, plot_height=100, title=None)
        s0.outline_line_color = "white"
        s7 = figure(tools=TOOLS, width=150, plot_height=100, title=None)
        s7.outline_line_color = "white"
        
        # Add bad PAs
        for ybad0,ybad1 in badPA:
            s5.patch([xlim0,xlim1,xlim1,xlim0], [ybad1,ybad1,ybad0,ybad0], color='white', alpha=0.7)
            s6.patch([0,100,100,0], [ybad1,ybad1,ybad0,ybad0], color='white', alpha=0.7)
        
        # put all the plots in a grid layout
        fig = gridplot(children=[[s7, s4, s1, s0], [s6, s5, s2, s3]])
        
        
    return fig

if __name__ == "__main__":
    #arguments RA & DEC, conversion to radians
    argv = sys.argv

    ra=argv[1]
    dec=argv[2]
    cubeNameSuf=argv[3]

    pamin=0 if len(argv)<5 else int(argv[4])
    pamax=360 if len(argv)<6 else int(argv[5])

    cubeName = argv[6] 
    targetName=None if len(argv)<8 else argv[7]
    save=False if len(argv)<8 else True #if name provided -> save
    tmpDir="." if len(argv)<9 else argv[8]
    os.makedirs(tmpDir, exist_ok=True)

    goodPA, badPA, _ = vpa.checkVisPA(ra, dec, targetName)

    #cubeName='cubes/cube_RA'+ra+'DEC'+dec+cubeNameSuf+'.fits'
    contam(cubeName,targetName=targetName,paRange=[pamin,pamax],badPA=badPA,tmpDir=tmpDir)
    


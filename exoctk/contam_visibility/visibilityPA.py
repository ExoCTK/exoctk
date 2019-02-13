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

<<<<<<< HEAD
import os
import pkg_resources
import datetime
import sys
=======

import sys
from . import ephemeris_old2x as EPH
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
<<<<<<< HEAD

from . import ephemeris_old2x as EPH
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from astropy.io import ascii
=======
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from astropy.io import ascii
import os
import pkg_resources
import datetime
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017
from jwst_gtvt.find_tgt_info import get_table

D2R = math.pi/180.  #degrees to radians
R2D = 180./math.pi #radians to degrees

def convert_ddmmss_to_float(astring):
    aline = astring.split(':')
    d = float(aline[0])
    m = float(aline[1])
    s = float(aline[2])
    hour_or_deg = (s/60.+m)/60.+d
    return hour_or_deg

def checkVisPA(ra, dec, targetName=None, ephFileName=pkg_resources.resource_filename('exoctk', 'data/contam_visibility/JWST_ephem_short.txt'), save=False, fig=''):

    if ra.find(':')>-1:  #format is hh:mm:ss.s or  dd:mm:ss.s
        ra = convert_ddmmss_to_float(ra) * 15. * D2R
        dec = convert_ddmmss_to_float(dec) * D2R
    else: #format is decimal
        ra = float(ra) * D2R
        dec = float(dec) * D2R

    #load ephemeris
    eclFlag = False
    eph = EPH.Ephemeris(ephFileName, eclFlag)

    #convert dates from MJD to Gregorian calendar dates
    mjd = np.array(eph.datelist)
    d = mdates.julian2num(mjd+2400000.5)
    gd = mdates.num2date(d)

    #loop through dates and determine VIS and PAs (nominal, min, max)
    vis = np.empty(mjd.size,dtype=bool)
    paNom, paMin, paMax = np.empty(mjd.size), np.empty(mjd.size), np.empty(mjd.size)
    for i in range(mjd.size):

        #is it visible?
        vis[i] = eph.in_FOR(mjd[i],ra,dec)

        #nominal PA at this date
        pa = eph.normal_pa(mjd[i],ra,dec)

        #search for minimum PA allowed by roll
        pa0 = pa
        while eph.is_valid(mjd[i],ra,dec,pa0-0.002):
            pa0 -= 0.002

        #search for maximum PA allowed by roll
        pa1 = pa
        while eph.is_valid(mjd[i],ra,dec,pa1+0.002):
            pa1 += 0.002

        paNom[i] = (pa*R2D)%360
        paMin[i] = (pa0*R2D)%360
        paMax[i] = (pa1*R2D)%360

    #does PA go through 360 deg?
    wrap = np.any(np.abs(np.diff(paNom[np.where(vis)[0]])) > 350)

    #Determine good and bad PA ranges
    #Good PAs
    i, = np.where(vis)
    pa = np.concatenate((paNom[i],paMin[i],paMax[i]))

    if wrap:
        pa = np.append(pa,(0.,360.))
    pa.sort()

    i1, = np.where(np.diff(pa)>10)
    i0 = np.insert(i1+1,0,0)
    i1 = np.append(i1,-1)
    paGood = np.dstack((pa[i0],pa[i1])).round(1).reshape(-1,2).tolist()

    #bad PAs (complement of the good PAs)
    paBad = []
    if paGood[0][0]>0:
        paBad.append([0.,paGood[0][0]])
    for i in range(1,len(paGood)):
        paBad.append([paGood[i-1][1],paGood[i][0]])
    if paGood[-1][1]<360.:
        paBad.append([paGood[-1][1],360.])

    #print results to file

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

    # Make a figure
    if not fig or fig==True:
        fig = plt.gcf()

    # Do all figure calculations
    iBad, = np.where(vis==False)
    paMasked = np.copy(paNom)
    paMasked[iBad] = np.nan
    gdMasked = np.copy(gd)

    i = np.argmax(paNom)
    if paNom[i+1]<10:
        i+=1
    paMasked = np.insert(paMasked,i,np.nan)
    gdMasked = np.insert(gdMasked,i,gdMasked[i])

    i = np.argmax(paMin)
    goUp = paMin[i-2]<paMin[i-1] #PA going up at wrap point?

    # Top part
    i0_top = 0 if goUp else i
    i1_top = i if goUp else paMin.size-1
    paMaxTmp = np.copy(paMax)
    paMaxTmp[np.where(paMin>paMax)[0]] = 360

    # Bottom part
    i = np.argmin(paMax)
    i0_bot = i if goUp else 0
    i1_bot = paMin.size-1 if goUp else i
    paMinTmp = np.copy(paMin)
    paMinTmp[np.where(paMin>paMax)[0]] = 0

    # Add fits to matplotlib
    if isinstance(fig, matplotlib.figure.Figure):

        # Make axes
        ax = plt.axes()
        plt.title(targetName)

        #plot nominal PA
        plt.plot(gdMasked,paMasked,color='k')

        #plot ranges allowed through roll
        if wrap:
            i = np.argmax(paMin)
            goUp = paMin[i-2]<paMin[i-1] #PA going up at wrap point?

            #top part
            plt.fill_between(gd[i0_top:i1_top+1],paMin[i0_top:i1_top+1],paMaxTmp[i0_top:i1_top+1],where=vis[i0_top:i1_top+1],lw=0,facecolor='k',alpha=0.5)

            #bottom part
            plt.fill_between(gd[i0_bot:i1_bot+1],paMinTmp[i0_bot:i1_bot+1],paMax[i0_bot:i1_bot+1],where=vis[i0_bot:i1_bot+1],lw=0,facecolor='k',alpha=0.5)

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

    # Or to bokeh!
    else:

        # Convert datetime to a number for Bokeh
        gdMaskednum = [datetime.date(2019, 6, 1)+datetime.timedelta(days=n) for n,d in enumerate(gdMasked)]
        color = 'green'

        # Draw the curve and error
        fig.line(gdMaskednum, paMasked, legend='cutoff', line_color=color)

        # Top
        err_y = np.concatenate([paMin[i0_top:i1_top+1],paMaxTmp[i0_top:i1_top+1][::-1]])
        # err_x = np.concatenate([[d.timestamp() for d in gd[i0_top:i1_top+1]],[d.timestamp() for d in gd[i0_top:i1_top+1]][::-1]])
        err_x = np.concatenate([gdMaskednum[i0_top:i1_top+1],gdMaskednum[i0_top:i1_top+1][::-1]])
<<<<<<< HEAD

=======
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017
        fig.patch(err_x, err_y, color=color, fill_alpha=0.2, line_alpha=0)

        # Bottom
        err_y = np.concatenate([paMinTmp[i0_bot:i1_bot+1],paMax[i0_bot:i1_bot+1][::-1]])
        # err_x = np.concatenate([[d.timestamp() for d in gd[i0_bot:i1_bot+1]],[d.timestamp() for d in gd[i0_bot:i1_bot+1]][::-1]])
        err_x = np.concatenate([gdMaskednum[i0_bot:i1_bot+1],gdMaskednum[i0_bot:i1_bot+1][::-1]])
        fig.patch(err_x, err_y, color=color, fill_alpha=0.2, line_alpha=0)

        print('err_y and err_x')
        print(err_y)
        print(err_x)
        print('shapes')
        print(np.shape(err_y))
        print(np.shape(err_x))

        # Plot formatting
        fig.xaxis.axis_label = 'Date'
        fig.yaxis.axis_label = 'Position Angle (degrees)'

    return paGood, paBad, gd, fig

<<<<<<< HEAD
def using_gtvt(ra, dec, instrumentName, targetName=None, save=False, \
    ephFileName=pkg_resources.resource_filename('exoctk', 'data/contam_visibility/JWST_ephem_short.txt'), \
    fig=''):
=======
def using_gtvt(ra, dec, instrumentName, targetName=None, save=False, ephFileName=pkg_resources.resource_filename('exoctk', 'data/contam_visibility/JWST_ephem_short.txt'), fig=''):
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017
    """using gtvt to find PAmin and PAmax for NIRISS
    yay

    test

    parameters
    ----------
    ra : str
        right ascencion
    dec : str
        declination
    instrumentName : str
        name of the instrument. can either be (case-sensitive):
        NIRISS, NIRCam, MIRI, FGS, NIRSpec
    """
    # getting calculations from GTVT (General Target Visibility Tool)
    tab = get_table(ra, dec)
    gd = tab['Date']
<<<<<<< HEAD
    print(gd)
    print('stop')
=======
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017
    paMin = tab[str(instrumentName)+' min']
    paMax = tab[str(instrumentName)+' max']
    paNom = tab['V3PA']

    #loop through dates and determine VIS
<<<<<<< HEAD
    #load ephemeris
    eclFlag = False
    eph = EPH.Ephemeris(ephFileName, eclFlag)

    mjd = np.array(eph.datelist)
    #print(mjd.size)
    #for i in range(mjd.size):

        #is it visible?

    #    vis[i] = eph.in_FOR(mjd[i],float(ra),float(dec))

    #does PA go through 360 deg?
    #wrap = np.any(np.abs(np.diff(paNom[np.where(vis)[0]])) > 350)
=======
    mjd = np.array(eph.datelist)
    for i in range(mjd.size):

        #is it visible?
        vis[i] = eph.in_FOR(mjd[i],ra,dec)

    #does PA go through 360 deg?
    wrap = np.any(np.abs(np.diff(paNom[np.where(vis)[0]])) > 350)
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017

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
    # Make a figure
    if not fig or fig==True:
        fig = plt.gcf()



    # Top part
    #i0_top = 0 if goUp else i
    #i1_top = i if goUp else paMin.size-1
    #paMaxTmp = np.copy(paMax)
    #paMaxTmp[np.where(paMin>paMax)[0]] = 360

    # Bottom part
    #i = np.argmin(paMax)
    #i0_bot = i if goUp else 0
    #i1_bot = paMin.size-1 if goUp else i
    #paMinTmp = np.copy(paMin)
    #paMinTmp[np.where(paMin>paMax)[0]] = 0

    # Add fits to matplotlib
    tab = get_table(ra, dec)
    gd = tab['Date'] # gd: gregorian date
    paMin = tab[str(instrumentName)+' min']
    paMax = tab[str(instrumentName)+' max']
    paNom = tab['V3PA']

    pas = np.arange(1, 361) # all possible angles (stops at 360)
<<<<<<< HEAD


    paGood, paBad= [], []
    for min, max in zip(paMin, paMax):

        good = pas[(pas >= min) & (pas <= max)] # good angles: pamin < pa < pamax
        bad = pas[(pas < min) | (pas > max)] # bad angles: pamax < pa < pamin


        if (len(good) > 0) & (len(bad > 0)):
            paGood.extend(good)
            paBad.extend(bad)


    print('paGood')
    print(paGood)

    #for arr in paGood:

        #newgd = gd[np.where(np.in1d(paMin, arr)==True)]
        #gds.append(newgd)

    # Make a figure
    if not fig or fig==True:
        fig = plt.gcf()


=======
    paGood = pas[(pas >= paMin) & (pas <= paMax)] # good angles: pamin < pa < pamax
    paBad = pas[(pas < paMin) | (pas > paMax)] # bad angles: pamax < pa < pamin
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017

    if isinstance(fig, matplotlib.figure.Figure):

        # Make axes
        ax = plt.axes()

        # Add grid
        ax.grid(b=True, which='major', color='black', alpha=0.1, linestyle='--')

        # plot nominal PA
        plt.plot(gd, paNom, color='k', lw=0.7, label='nominal V3PA angle (Deg)')

        # plot ranges allowed through roll
        #if wrap:
            #i = np.argmax(paMin)
            #goUp = paMin[i-2]<paMin[i-1] #PA going up at wrap point?

            #top part
            #plt.fill_between(gd[i0_top:i1_top+1],paMin[i0_top:i1_top+1],paMaxTmp[i0_top:i1_top+1],where=vis[i0_top:i1_top+1],lw=0,facecolor='k',alpha=0.5)

        #    #bottom part
            #plt.fill_between(gd[i0_bot:i1_bot+1],paMinTmp[i0_bot:i1_bot+1],paMax[i0_bot:i1_bot+1],where=vis[i0_bot:i1_bot+1],lw=0,facecolor='k',alpha=0.5)

        #else:

<<<<<<< HEAD
        plt.fill_between(gd, paMin, paMax, where=paMax>paMin, lw=1.0, edgecolor=(0,0,0,.3), facecolor=(0,0,0,.5))
=======
        plt.fill_between(gd, paMin, paMax, where=paMax>paMin, lw=1.0, edgecolor=(0,0,0,.5), facecolor=(0,0,0,.5))
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017





        plt.ylabel('Position Angle (degrees)')
        plt.title('Target: Trappist-1'+'\n'+'Visibility with '+str(instrumentName)+' instrument - calculated using GTVT')
<<<<<<< HEAD
        #plt.xlim(min(gd),max(gd))
=======
        plt.xlim(min(gd),max(gd))
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b '%y"))
        ax.xaxis.set_minor_locator(mdates.DayLocator(list(range(1,32,5))))
        plt.ylim(0,360)
        ax.yaxis.set_major_locator(MultipleLocator(25))
        ax.yaxis.set_minor_locator(MultipleLocator(5))

        for label in ax.get_xticklabels():
            label.set_rotation(45)

        ax.set_xticks(ax.get_xticks()[::2])
        plt.legend()
<<<<<<< HEAD

    # Or to bokeh!
    #paMax = np.max(paGood)
    #paMin = np.min(paGood)

    else:
        #paMax = np.max(paGood)
        #paMin = np.min(paGood)
        #print('test')
        #print(paMax)
=======
    """
    # Or to bokeh!
    else:

>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017
        # Convert datetime to a number for Bokeh
        gd = [datetime.date(2019, 6, 1)+datetime.timedelta(days=n) for n,d in enumerate(gd)]
        color = 'green'

        # Draw the curve and error
        fig.line(gd, paNom, legend='cutoff', line_color=color)

        # Top
<<<<<<< HEAD
        #err_y = np.concatenate(np.asarray(paMin),np.asarray(paMax))
        # err_x = np.concatenate([[d.timestamp() for d in gd[i0_top:i1_top+1]],[d.timestamp() for d in gd[i0_top:i1_top+1]][::-1]])
        #err_x = np.concatenate([gdMaskednum[i0_top:i1_top+1],gdMaskednum[i0_top:i1_top+1][::-1]])
        #paMin, paMax = list(paMin), list(paMax)
        #PAs = paMin.extend(paMax)
        #gd = list(gd)
        #fig.patch(gd.extend(gd), PAs, color='green')

        # Bottom
        #err_y = np.concatenate([paMinTmp[i0_bot:i1_bot+1],paMax[i0_bot:i1_bot+1][::-1]])
        # err_x = np.concatenate([[d.timestamp() for d in gd[i0_bot:i1_bot+1]],[d.timestamp() for d in gd[i0_bot:i1_bot+1]][::-1]])
        #err_x = np.concatenate([gdMaskednum[i0_bot:i1_bot+1],gdMaskednum[i0_bot:i1_bot+1][::-1]])
        #fig.patch(g, err_y, color=color, fill_alpha=0.2, line_alpha=0)

        # Plot formatting
        fig.xaxis.axis_label = 'update test'
        fig.yaxis.axis_label = 'Position Angle (degrees)'


        from bokeh.plotting import save, output_file
        from bokeh.io import reset_output
        output_file('/Users/jmedina/Desktop/visib_test.html')
        save(fig)

    return paGood, paBad, gd
=======
        err_y = np.concatenate([paMin[i0_top:i1_top+1],paMaxTmp[i0_top:i1_top+1][::-1]])
        # err_x = np.concatenate([[d.timestamp() for d in gd[i0_top:i1_top+1]],[d.timestamp() for d in gd[i0_top:i1_top+1]][::-1]])
        err_x = np.concatenate([gdMaskednum[i0_top:i1_top+1],gdMaskednum[i0_top:i1_top+1][::-1]])
        fig.patch(err_x, err_y, color=color, fill_alpha=0.2, line_alpha=0)

        # Bottom
        err_y = np.concatenate([paMinTmp[i0_bot:i1_bot+1],paMax[i0_bot:i1_bot+1][::-1]])
        # err_x = np.concatenate([[d.timestamp() for d in gd[i0_bot:i1_bot+1]],[d.timestamp() for d in gd[i0_bot:i1_bot+1]][::-1]])
        err_x = np.concatenate([gdMaskednum[i0_bot:i1_bot+1],gdMaskednum[i0_bot:i1_bot+1][::-1]])
        fig.patch(err_x, err_y, color=color, fill_alpha=0.2, line_alpha=0)

        # Plot formatting
        fig.xaxis.axis_label = 'Date'
        fig.yaxis.axis_label = 'Position Angle (degrees)'
    """
    return paGood, paBad, gd, fig
>>>>>>> 56e76344f1b13e1346e36d5790c345f9ccbbb017

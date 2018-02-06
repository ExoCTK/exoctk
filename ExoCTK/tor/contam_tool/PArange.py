#! /usr/bin/env python
# #######################################################################
# IDENT       p_computeRangePA.py
# LANGUAGE    Python
# AUTHOR      P.FERRUIT from code created by WAYNE KINZEL      
# PURPOSE     Procedure to compute and display the range of PA
#             available for a given (RA,DEC) and over a given
#             search interval.
#
# Procedure derived from the code of Wayne Kinzel provided by Jeff Valenti
# Extract from the e-mail of Wayne Kinzel:
# As before, the code is not officially tested, nor is it an official STScI product.
# Users should be warned that the apparent position of the Sun changes ~+/-0.2 degrees
# depending upon where JWST is in its orbit.   So do not rely strongly on these results
# if the target is within ~0.2 degrees of |ecliptic latitude| 45 degrees or 85 degrees.
# For example if a target is at 84.9 degrees latitude and the tool says it is CVZ, it
# may not be with the operational orbit.
#
# SYNTAX
# usage: p_computeRangePA.py ra dec [ephemeris] [mjdmin] [mjdmax]
#
# positional arguments:
#    ra             Input right ascension coordinate (equatorial coordinates;
#                        dd:mm:ss.s or degrees).   
#    dec            Input declination coordinate (equatorial coordinates;
#                        dd:mm:ss.s or degrees).
#    out            Name of the output figure file.
# 
# optional arguments:
#    -e, --ephemeris   Input ephemeris file containing JWST position as a function
#                        of time. Defaulted to: JWST_elv4_20181031130000_ephem.txt
#    -mjdmin, --mjdmin Beginning of the search interval (in modified Julian date).
#                        Must be inside the period covered by the ephemeris.
#                        Defaulted to None in which case the ephemeris start date
#                        is used.
#    -mjdmax, --mjdmax Beginning of the search interval (in modified Julian date).
#                        Must be inside the period covered by the ephemeris.
#                        Defaulted to None in which case the ephemeris start date
#                        is used.
#
# VERSION
# 1.0.0 12.04.2015 PF Creation from the code of Wayne Kinzel
#
#########################################################################
version = '1.0.0'

import sys
import argparse
import datetime
import numpy
import string
import math
import matplotlib.pyplot as plt

import f_visibilityPeriods as f_visibilityPeriods
import ephemeris_old2x as EPH
import astropy.time

D2R = math.pi / 180.  #degrees to radians
R2D = 180. / math.pi #radians to degrees 
PI2 = 2. * math.pi   # 2 pi

def convert_ddmmss_to_float(astring):
    aline = astring.split(':')
    d= float(aline[0])
    m= float(aline[1])
    s= float(aline[2])
    hour_or_deg = (s/60.+m)/60.+d
    return hour_or_deg

# =======================================================================
# Declaring and loading the input arguments
# =======================================================================
parser = argparse.ArgumentParser()
parser.add_argument('ra', type=str, help='Input right ascension coordinate (equatorial coordinates; dd:mm:ss.s or degrees).')
parser.add_argument('dec', type=str, help='Input declination coordinate (equatorial coordinates; dd:mm:ss.s or degrees).')
parser.add_argument('-e', '--ephemeris', help='Input ephemeris file containing JWST position as a function of time.', default='JWST_elv4_20181031130000_ephem.txt')
parser.add_argument('-min', '--mjdmin', help='Start of the search interval (in modified Julian date).', default=None)
parser.add_argument('-max', '--mjdmax', help='End of the search interval (in modified Julian date).', default=None)

args = parser.parse_args()
argv = sys.argv
narg = len(argv)-1
print("=====================================")
print(argv[0])
print("Version: ", version)
print((datetime.datetime.now()).isoformat())
print("=====================================")

# =======================================================================
# Parsing the input arguments
# =======================================================================
# Input coordinates (two possible input formats)
if ':' in args.ra:  #format is hh:mm:ss.s or  dd:mm:ss.s  
  ra  = convert_ddmmss_to_float(args.ra) * 15. * D2R
  dec   = convert_ddmmss_to_float(args.dec) * D2R
else: #format is decimal
  ra  = float(args.ra) * D2R
  dec   = float(args.dec) * D2R
# Optional ephemeris filename
ephemerisFilename = args.ephemeris
mjdmin = args.mjdmin
mjdmax = args.mjdmax
    
print("###############################################")
print("{:s}".format(sys.argv[0]))
print("# RA  = {:12.8f} rad = {:12.8f} deg".format(ra, ra / D2R))
print("# DEC = {:12.8f} rad = {:12.8f} deg".format(dec, dec / D2R))
print("###############################################")


# =====================================================================================
# Loading the ephemeris table from file 
# =====================================================================================
# Ephemeris is in equatorial coordinates, setting the ecliptic coordinate flag to
# False
ECL = False
print("# Loading the ephemeris for JWST [{:s}].".format(ephemerisFilename))
A_eph = EPH.Ephemeris(ephemerisFilename, ECL)
rangeOfTimeEphemeris = astropy.time.Time([A_eph.amin, A_eph.amax], format='mjd')
print(rangeOfTimeEphemeris.iso[0])
print("# Start date of the ephemeris:   {:s} [{:f} MJD]".format(rangeOfTimeEphemeris.iso[0], A_eph.amin))
print("# End date of the ephemeris:     {:s} [{:f} MJD]".format(rangeOfTimeEphemeris.iso[1], A_eph.amax))
if ((mjdmin == None) or (float(mjdmin) < A_eph.amin)):
    mjdmin = A_eph.amin
    #print "# The start of the search interval has been set to the start of the period covered by the ephemeris."
else:
    mjdmin = float(mjdmin)
if ((mjdmax == None) or (float(mjdmax) > A_eph.amax)):
    mjdmax = A_eph.amax
    #print "# The end of the search interval has been set to the end of the period covered by the ephemeris."
else:
    mjdmax = float(mjdmax)
rangeOfTime = astropy.time.Time([mjdmin, mjdmax], format='mjd')
print("# Start date of the search  interval:   {:s} [{:f} MJD]".format(rangeOfTime.iso[0], mjdmin))
print("# End date of the search interval:     {:s} [{:f} MJD]".format(rangeOfTime.iso[1], mjdmax))
print("###############################################")

# =====================================================================================
# Searching the visibility periods
# =====================================================================================
startList, endList, statusList = f_visibilityPeriods.f_computeVisibilityPeriods(A_eph, mjdmin, mjdmax, ra, dec)

# =====================================================================================
# Printing the results
# =====================================================================================
print("|                               Observation window(s)                                         |         Normal V3 PA [deg] | Equatorial coordinates of the object")
print("Start date              End date                  Start[MJD]    End[MJD]          Dur.   Type   Start         End            coord1        coord2")
indexWindow = 0
for currentStatus in statusList:
    currentWStart = startList[indexWindow]
    currentWEnd = endList[indexWindow]
    currentDuration = currentWEnd - currentWStart
    rangeOfTime = astropy.time.Time([currentWStart, currentWEnd], format='mjd')
    pa_start = A_eph.normal_pa(currentWStart, ra, dec)
    pa_end   = A_eph.normal_pa(currentWEnd, ra, dec)
    print("%s %s %13.5f %13.5f %11.2f  %d %13.5f %13.5f %13.5f %13.5f " % (rangeOfTime.iso[0], rangeOfTime.iso[1], currentWStart, currentWEnd, currentDuration, currentStatus, pa_start*R2D,pa_end*R2D,ra*R2D,dec*R2D))
    indexWindow += 1



import numpy as np
from astroquery.irsa import Irsa
from astroquery.simbad import Simbad
import astropy.units as u	

def resolve_target(targetName):
	try:
		target_info = Simbad.query_object(targetName)
		targetRA = target_info['RA']
		targetDEC = target_info['DEC']
		ra = (targetRA[0].replace(' ', ':'))
		dec = (targetDEC[0].replace(' ', ':'))
		return ra, dec
		
	except:
		ra = 'unresolved'
		dec = 'unresolved'
		return ra, dec
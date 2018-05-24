#! /usr/bin/env python
# Version 1. August 2, 2010
# Version 2. August 3, 2010
#   Got rid of degrees trig functions

from math import *

D2R = pi/180.
R2D = 180./pi
PI2 = 2. * pi
epsilon = 23.43929 * D2R #obliquity of the ecliptic J2000
unit_limit = lambda x: min(max(-1.,x),1.)


def pa(tgt_c1,tgt_c2,obj_c1,obj_c2):
    """calculates position angle of object at tgt position."""
    y = cos(obj_c2)*sin(obj_c1-tgt_c1)
    x = (sin(obj_c2)*cos(tgt_c2)-cos(obj_c2)*sin(tgt_c2)*cos(obj_c1-tgt_c1))
    p = atan2(y,x)
    if p < 0.: p += PI2
    if p >= PI2: p -= PI2
    return p

def delta_pa_no_roll(pos1_c1,pos1_c2,pos2_c1,pos2_c2):
    """Calculates the change in position angle between two positions with no roll about V1"""
    u = (sin(pos1_c2) + sin(pos2_c2)) * sin(pos2_c1 - pos1_c1)
    v = cos(pos2_c1 - pos1_c1) + cos(pos1_c2)*cos(pos2_c2)+ sin(pos1_c2)*sin(pos2_c2)*cos(pos2_c1 - pos1_c1)
    return atan2(u,v)

def dist(obj1_c1,obj1_c2,obj2_c1,obj2_c2):
    """angular distance betrween two objects, positions specified in spherical coordinates."""
    x = cos(obj2_c2)*cos(obj1_c2)*cos(obj2_c1-obj1_c1) + sin(obj2_c2)*sin(obj1_c2)
    return acos(unit_limit(x))
    
def JWST_same_ori(tgt0_c1,tgt0_c2,p0,tgt_c1,tgt_c2):
    """Calculates normal orientation of second target, given first target's orientation is normal.
This is in Ecliptic coordinates!""" 
    long_sun = atan2(-sind(p0)*sin(tgt0_c2),cosd(p0)) + tgt0_c1
    pp = atan2(-sin(long_sun-tgt_c1),cos(long_sun-tgt_c1)*sin(tgt_c2))
    if pp<0.:
        pp += PI2
    return pp


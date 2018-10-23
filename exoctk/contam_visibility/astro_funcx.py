#! /usr/bin/env python
# Version 1. August 2, 2010
# Version 2. August 3, 2010
#   Got rid of degrees trig functions
# Joe Filippazzo, 2018/06/26, Made PEP8 compliant

from math import cos, sin, atan2, acos, pi

from .math_extensionsx import sind, cosd

D2R = pi/180.
R2D = 180./pi
PI2 = 2. * pi
epsilon = 23.43929 * D2R  # obliquity of the ecliptic J2000


def pa(tgt_c1, tgt_c2, obj_c1, obj_c2):
    """Calculates position angle of object at tgt position.

    Parameters
    ----------
    tgt_c1: float
        The RA of the target.
    tgt_c2: float
        The Dec of the target.
    obj_c1: float
        The RA of the reference.
    obj_c2: float
        The Dec of the reference.

    Returns
    -------
    float
        The position angle.
    """
    y = cos(obj_c2) * sin(obj_c1-tgt_c1)
    c = cos(obj_c2) * sin(tgt_c2) * cos(obj_c1-tgt_c1)
    x = sin(obj_c2) * cos(tgt_c2) - c
    p = atan2(y, x)
    if p < 0.:
        p += PI2
    if p >= PI2:
        p -= PI2
    return p


def delta_pa_no_roll(pos1_c1, pos1_c2, pos2_c1, pos2_c2):
    """Calculates the change in position angle between two positions with no
    roll about V1

    Parameters
    ----------
    pos1_c1: float
        The RA of the first position.
    pos1_c2: float
        The Dec of the first position.
    pos2_c1: float
        The RA of the second position.
    pos2_c2: float
        The Dec of the second position.

    Returns
    -------
    float
        The change in position angle.
    """
    u = (sin(pos1_c2) + sin(pos2_c2)) * sin(pos2_c1 - pos1_c1)
    a = sin(pos2_c2) * cos(pos2_c1 - pos1_c1)
    v = cos(pos2_c1 - pos1_c1) + cos(pos1_c2) * cos(pos2_c2) + sin(pos1_c2) * a
    return atan2(u, v)


def dist(obj1_c1, obj1_c2, obj2_c1, obj2_c2):
    """Calculates the angular distance betrween two objects, positions
    specified in spherical coordinates

    Parameters
    ----------
    obj1_c1: float
        The RA of the first object.
    obj1_c2: float
        The Dec of the first object.
    obj2_c1: float
        The RA of the second object.
    obj2_c2: float
        The Dec of the second object.

    Returns
    -------
    float
        The distance between the objects.
    """
    b = sin(obj2_c2) * sin(obj1_c2)
    x = cos(obj2_c2) * cos(obj1_c2) * cos(obj2_c1-obj1_c1) + b
    return acos(unit_limit(x))


def JWST_same_ori(tgt0_c1, tgt0_c2, p0, tgt_c1, tgt_c2):
    """Calculates normal orientation of second target, given first target's
    orientation is normal. This is in Ecliptic coordinates!

    Parameters
    ----------
    tgt0_c1: float
        The RA of the first target.
    tgt0_c2: float
        The Dec of the first target.
    p0: float
        The origin.
    tgt1_c1: float
        The RA of the second target.
    tgt1_c2: float
        The Dec of the second target.

    Returns
    -------
    float
        The normal orientation.
    """
    long_sun = atan2(-sind(p0) * sin(tgt0_c2), cosd(p0)) + tgt0_c1
    pp = atan2(-sin(long_sun-tgt_c1), cos(long_sun-tgt_c1) * sin(tgt_c2))
    if pp < 0.:
        pp += PI2
    return pp


def unit_limit(x):
    """ Forces value to be in [-1, 1]

    Parameters
    ----------
    x: float, int
        The value to adjust.

    Retruns
    -------
    float
        The adjusted value.
    """
    return min(max(-1., x), 1.)

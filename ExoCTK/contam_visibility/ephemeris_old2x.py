# ! /usr/bin/env python
# Module ephemeris.py
import math
import sys
import time

from . import astro_funcx as astro_func
from . import quaternionx as qx
from . import time_extensionsx as time2

D2R = math.pi/180.  # degrees to radians
R2D = 180. / math.pi  # radians to degrees
PI2 = 2. * math.pi  # 2 math.pi
MIN_SUN_ANGLE = 84.8 * D2R  # minimum Sun angle, in radians
MAX_SUN_ANGLE = 135.0 * D2R  # maximum Sun angle, in radians
SUN_ANGLE_PAD = 0.5 * D2R   # pad away from Sun angle limits

Qecl2eci = qx.QX(23.439291*D2R)  # At J2000 equinox


class Ephemeris:
    """A class for the ephemeris of an observation.

    History
    -------
    07/31/2008 Added functions for sun position
    12/05/2008 Switched to spline fit.
    07/29/2010 Added OP window calc
    08/03/2010 Got rid of degrees trig functions
              Removed in_FOR, is_valid_att etc. as those are S/C dependent
    """
    def __init__(self, ephem_file, cnvrt=False):
        """Ephemeris constructor

        Parameters
        ----------
        ephem_file: str
            The path to the ephemeris file.
        cnvrt: bool, optional
            Converts into Ecliptic frame.
        """
        if cnvrt:
            print("Using Ecliptic Coordinates")
        else:
            print("Using Equatorial Coordinates")
        self.datelist = []
        self.xlist = []
        self.ylist = []
        self.zlist = []
        self.amin = 0.
        self.amax = 0.
        aV = qx.Vector(0., 0., 0.)
        if ephem_file.find("l2_halo_FDF_060619.trh") > -1:
            ascale = 0.001
        else:
            ascale = 1.0
        fin = open(ephem_file, 'r')
        finLines = fin.readlines()
        for item in finLines[2:]:
            item = item.strip()
            item = item.split()
            adate = time2.mjd_from_string(item[0])  # represent dates as mjds
            x = float(item[1])*ascale
            y = float(item[2])*ascale
            z = float(item[3])*ascale
            if cnvrt:
                aV.set_eq(x, y, z)
                ll = aV.length()
                aV = aV/ll
                aV = Qecl2eci.inv_cnvrt(aV)
                aV = aV*ll
                x = aV.rx()
                y = aV.ry()
                z = aV.rz()
            self.datelist.append(adate)
            self.xlist.append(x)
            self.ylist.append(y)
            self.zlist.append(z)

            if self.amin == 0.:
                self.amin = adate
        self.amax = adate
        fin.close()

        
    def report_ephemeris(self, limit=100000, pathname=None):
        """Prints a formatted report of the ephemeris.

        Parameters
        ----------
        limit: int, optional
            The number of records to report.
        pathname: str, optional
            The path to a file to hold the report.
        """
        num_to_report = min(limit, len(self.datelist))

        if (pathname):
            dest = open(pathname, 'w')
            print(('# Generated %s\n' % (time.ctime())), file=dest)
        else:
            dest = sys.stdout  # defaults to standard output

        print(('%17s  %14s  %14s  %14s\n' % ('DATE      ', 'X (KM)   ',
               'Y (KM)   ', 'Z (KM)   ')), file=dest)

        for num in range(num_to_report):
            date = self.datelist[num]
            x = self.xlist[num]
            y = self.ylist[num]
            z = self.zlist[num]

            fmt = time2.display_date(date), x, y, z
            print(('%17s  %14.3f  %14.3f  %14.3f' % fmt), file=dest)

        if (pathname):
            dest.close()   # Clean up

            
    def pos(self, adate):
        """Computes the position of the telescope at a given date using the
        grid of positions of the ephemeris as a starting point and
        applying a linear interpolation between the ephemeris grid points

        Parameters
        ----------
        adate: datetime.datetime object
            The date of the observation.

        Returns
        -------
        qx.Vector object
            The position of the telescope as a Vector object
        """
        cal_days = adate - self.datelist[0]
        index = int(cal_days)
        if (index == len(self.datelist)-1):
            index = index - 1
        frac = cal_days - index
        x = (self.xlist[index+1] - self.xlist[index])*frac + self.xlist[index]
        y = (self.ylist[index+1] - self.ylist[index])*frac + self.ylist[index]
        z = (self.zlist[index+1] - self.zlist[index])*frac + self.zlist[index]
        return qx.Vector(x, y, z)

    
    def Vsun_pos(self, adate):
        """The vector of the sun at the given date

        Parameters
        ----------
        adate: datetime
            The date of the observation

        Returns
        -------
        Vector
            The position of the sun as a Vector object
        """
        Vsun = -1. * self.pos(adate)
        Vsun = Vsun / Vsun.length()
        return Vsun

    
    def sun_pos(self, adate):
        """The coordinates of the sun at the given date

        Parameters
        ----------
        adate: datetime.datetime object
            The date of the observation.

        Returns
        -------
        tuple
            The coordinates of the sun.
        """
        Vsun = -1. * self.pos(adate)
        Vsun = Vsun / Vsun.length()
        coord2 = math.asin(unit_limit(Vsun.z))
        coord1 = math.atan2(Vsun.y, Vsun.x)
        if coord1 < 0.:
            coord1 += PI2
        return (coord1, coord2)

    
    def normal_pa(self, adate, tgt_c1, tgt_c2):
        """Calculate the V3 position

        Parameters
        ----------
        adate: datetime.datetime object
            The date of the observation.
        tgt_c1: float
            The RA in radians.
        tgt_c2: float
            The Dec in radians.

        Returns
        -------
        float
            The V3 position.
        """
        (sun_c1, sun_c2) = self.sun_pos(adate)
        sun_pa = astro_func.pa(tgt_c1, tgt_c2, sun_c1, sun_c2)
        V3_pa = sun_pa + math.pi  # We want -V3 pointed towards sun.
        if V3_pa < 0.:
            V3_pa += PI2
        if V3_pa >= PI2:
            V3_pa -= PI2
        return V3_pa

    
    def long_term_attitude(self, date):
        """Defines a long-term safe attitude as of a given date.

        Parameters
        ----------
        date: float
             The date of computation, as an mjd.

        Returns
        -------
        Attitude
            The Attitude object at the given date.
        """
        # Retrieve Sun's position and transform to ecliptic coordinates.
        (sun_ra, sun_dec) = self.sun_pos(date)   # RA range 0-PI2
        vSun = qx.CelestialVector(sun_ra, sun_dec, degrees=False)
        vSun = vSun.transform_frame('ec')

        # Now subtract the minimum Sun angle plus a pad from the
        # ecliptic longitude.
        # Sun angle steadily decreases as the Earth (and JWST with it)
        # revolve counterclockwise, so this should maximize the duration
        # when the attitude is within the FOR.  Set the latitude to
        # 0 (ecliptic).
        # Normalize to 0-360 degrees, and convert back into equatorial
        # coordinates.
        # Then set the normal PA, which should be valid at the ecliptic for
        # the duration of the visibility window.
        longitude = vSun.ra - MIN_SUN_ANGLE - SUN_ANGLE_PAD

        if (longitude < 0):
            longitude = longitude + PI2

        vec1 = qx.CelestialVector(longitude, 0.0, frame='ec', degrees=False)
        vec1 = vec1.transform_frame('eq')
        pa = self.normal_pa(date, vec1.ra, vec1.dec)
        return(qx.Attitude(vec1.ra, vec1.dec, pa, degrees=False))

    
    def is_valid(self, date, ngc_1, ngc_2, V3pa):
        """Indicates whether an attitude is valid at a given date.

        Parameters
        ----------
        date: float
            The date of the observation.
        ngc_1: flaot
            The RA of the reference in radians.
        ngc_2: float
            The Dec of the reference in radians.
        V3pa: float
            The V3 position of the telescope.

        Returns
        -------
        bool
            Is it a valid PA.
        """
        # First check that the date is within the time interval of
        # the ephemeris.
        if ((date < self.amin) or (date > self.amax)):
            return False

        (sun_1, sun_2) = self.sun_pos(date)
        d = astro_func.dist(ngc_1, ngc_2, sun_1, sun_2)
        vehicle_pitch = math.pi/2 - d   # see JI memo from May 2006
        # sun pitch is always equal or greater than sun angle (V1 to sun)
        if (d < MIN_SUN_ANGLE or d > MAX_SUN_ANGLE):
            return False
        # now checking the roll and pitch angle combination
        pa = astro_func.pa(ngc_1, ngc_2, sun_1, sun_2) + math.pi
        roll = math.acos(math.cos(V3pa - pa))
        sun_roll = math.asin(math.sin(roll) * math.cos(vehicle_pitch))
        if (abs(sun_roll) <= 5.2*D2R):
            sun_pitch = math.atan2(math.tan(vehicle_pitch), math.cos(roll))
            if (sun_pitch <= 5.2*D2R and sun_pitch >= -45.*D2R):
                return True
        return False

    
    def in_FOR(self, date, ngc_1, ngc_2):
        """Test if in the FOR

        Parameters
        ----------
        date: float
            The date of the observation.
        ngc_1: flaot
            The RA of the reference in radians.
        ngc_2: float
            The Dec of the reference in radians.

        Returns
        -------
        bool
            Is it in the FOR.
        """
        (sun_1, sun_2) = self.sun_pos(date)
        d = astro_func.dist(ngc_1, ngc_2, sun_1, sun_2)
        # sun pitch is always equal or greater than sun angle (V1 to sun)
        if (d < MIN_SUN_ANGLE or d > MAX_SUN_ANGLE):
            return False
        return True

    
    def bisect_by_FOR(self, in_date, out_date, ngc_1, ngc_2):
        """Find the midpoint in time between in and out of FOR,
        assumes only one "root" in interval

        Parameters
        ----------
        in_date: float
            The in date of the observation.
        out_date: float
            The out date of the observation.
        ngc_1: flaot
            The RA of the reference in radians.
        ngc_2: float
            The Dec of the reference in radians.

        Returns
        -------
        float
            The midpoint in time.
        """
        delta_days = 200.
        mid_date = (in_date+out_date)/2.
        while delta_days > 0.000001:
            (sun_1, sun_2) = self.sun_pos(mid_date)
            d = astro_func.dist(ngc_1, ngc_2, sun_1, sun_2)
            if (d > MAX_SUN_ANGLE or d < MIN_SUN_ANGLE):
                out_date = mid_date
            else:
                in_date = mid_date
            mid_date = (in_date+out_date)/2.
            delta_days = abs(in_date-out_date)/2.
            # print "UU", mid_date
        # ensure returned date always in FOR
        if in_date > out_date:
            mid_date = mid_date + 0.000001
        else:
            mid_date = mid_date - 0.000001
        return mid_date

    
    def bisect_by_attitude(self, in_date, out_date, ngc_1, ngc_2, pa):
        """Find the midpoint in time between in and out of FOR,
        assumes only one "root" in interval

        Parameters
        ----------
        in_date: float
            The in date of the observation.
        out_date: float
            The out date of the observation.
        ngc_1: flaot
            The RA of the reference in radians.
        ngc_2: float
            The Dec of the reference in radians.
        pa: float
            The position angle.

        Returns
        -------
        float
            The midpoint in time.
        """
        icount = 0
        delta_days = 200.
        mid_date = (in_date+out_date)/2.
        # print "bisect >", in_date, out_date, abs(in_date-out_date )
        while delta_days > 0.000001:
            if self.is_valid(mid_date, ngc_1, ngc_2, pa):
                in_date = mid_date
            else:
                out_date = mid_date
            mid_date = (in_date+out_date)/2.
            delta_days = abs(in_date-out_date)/2.
            # print "UU", mid_date
            icount = icount + 1
        # print " bisected >", icount
        return mid_date

    
    def OP_window(self, adate, ngc_1, ngc_2, pa, mdelta, pdelta):
        """Attitude at adate must be valid, else returns (0, 0).
        If valid, returns (adate-mdelta, adate+pdelta) or the constraint
        window, which ever is smaller.

        Parameters
        ----------
        date: float
            The date of the observation.
        ngc_1: flaot
            The RA of the reference in radians.
        ngc_2: float
            The Dec of the reference in radians.
        pa: float
            The position angle.
        mdelta: float
            Delta time on low end.
        pdelta: float
            Delta time on high end.

        Returns
        -------
        tuple
            The OP window.
        """
        if self.is_valid(adate, ngc_1, ngc_2, pa):
            if self.is_valid(adate-mdelta, ngc_1, ngc_2, pa):
                OP_min = adate-mdelta
            else:
                OP_min = self.bisect_by_attitude(adate, adate-mdelta,
                                                 ngc_1, ngc_2, pa)
            if self.is_valid(adate+pdelta, ngc_1, ngc_2, pa):
                OP_max = adate+pdelta
            else:
                OP_max = self.bisect_by_attitude(adate, adate+pdelta,
                                                 ngc_1, ngc_2, pa)
        else:
            OP_min = 0.
            OP_max = 0.
        return (OP_min, OP_max)


def unit_limit(x):
    """ Forces value to be in [-1, 1].

    Parameters
    ----------
    x: float, int
        The value to adjust.

    Returns
    -------
    float
        The adjusted value.
    """
    return min(max(-1., x), 1.)

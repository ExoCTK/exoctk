"""
Series of functions to compute the visibility periods for a given (RA,DEC)
with in some cases the possibility to select a PA value.

Functions derived from the code of Wayne Kinzel provided by Jeff Valenti
Extract from the e-mail of Wayne Kinzel:
As before, the code is not officially tested, nor is it an official STScI
product. Users should be warned that the apparent position of the Sun changes
~+/-0.2 degrees epending upon where JWST is in its orbit.
So do not rely strongly on these results if the target is within ~0.2 degrees
of |ecliptic latitude| 45 degrees or 85 degrees.
For example if a target is at 84.9 degrees latitude and the tool says it is
CVZ, it may not be with the operational orbit.
"""
import math

D2R = math.pi / 180.  # degrees to radians
R2D = 180. / math.pi  # radians to degrees
PI2 = 2. * math.pi   # 2 pi


def f_computeVisibilityPeriods(ephemeris, mjdmin, mjdmax, ra, dec):
    """
    # -----------------------------------------------------------
    # METHOD         f_computeVisibilityPeriods()
    # TYPE           function
    #
    # DESCRIPTION    function that will compute the visibility
    #                periods for a given (RA,DEC) over a given
    #                time period.
    #
    # SYNTAX    f_computeVisibilityPeriods(ephemeris, mjdmin,
    #                    mjdmax, ra, dec)
    #
    # ephemeris: input ephemeris object
    # mjdmin: beginning of the search interval (modified
    #    Julian date). It must be covered by the ephemeris.
    # mjdmax: end of the search interval (modified
    #    Julian date). It must be covered by the ephemeris.
    # ra: input RA coordinate (equatorial coordinate, in rad)
    # dec: input DEC coordinate (equatorial coordinate, in rad)
    #
    # Returns two lists containing the start end end of each
    # visibility period and a list containing a status flag:
    # flag = 0 visibility period fully in the search interval
    # flag = -1 start of the visibility period truncated by
    # the start of the search interval
    # flag = -2 end of the visibility period truncated by
    # the end of the search interval
    # flag = +1 the search interval is fully included in
    # the visibility period
    #
    # -----------------------------------------------------------
    """
    # ===========================================================
    # Paranoid checks
    # ===========================================================
    # print "# RA  = {:12.8f} rad = {:12.8f} deg".format(ra, ra / D2R)
    # print "# DEC = {:12.8f} rad = {:12.8f} deg".format(dec, dec / D2R)
    # print"# No constraint on the PA."
    if (ephemeris.amin > mjdmin):
        print("""f_computeVisibilityPeriods(): the start of the search\
                 interval is not covered by the ephemeris.""")
        print("""Ephemeris start date (modified Julian date):\
                 {:8.5f}""".format(ephemeris.amin))
        print("""Search interval start date (modified Julian date):\
                 {:8.5f}""".format(mjdmin))
        raise ValueError
    if (ephemeris.amax < mjdmax):
        print("""f_computeVisibilityPeriods(): the end of the search interval\
                 is not covered by the ephemeris.""")
        print("""Ephemeris end date (modified Julian date):\
                 {:8.5f}""".format(ephemeris.amax))
        print("""Search interval end date (modified Julian date):\
                 {:8.5f}""".format(mjdmax))
        raise ValueError
    # ===========================================================
    # Scanning the search period
    # ===========================================================
    # Flag used to track the beginning and the end of a
    # visibility period
    iflip = False
    wstart = mjdmin
    startList = []
    endList = []
    statusList = []
    # Scannning step size (must be small enough to make sure that
    # it cannot contain a full vsibility period (we would miss
    # it)
    scanningStepSize = 0.1
    span = int((mjdmax - mjdmin) / scanningStepSize)
    # Initialisation (first step of the scan is outside from the
    # loop
    iflag_old = ephemeris.in_FOR(mjdmin, ra, dec)
    for i in range(span):
        # Current date (the last step may be partial to remain
        # within the search interval
        currentdate = mjdmin + (i + 1) * scanningStepSize
        if (currentdate >= mjdmax):
            currentdate = mjdmax
        iflag = ephemeris.in_FOR(currentdate, ra, dec)
        # Checking if we are reaching the beginning or the end of a
        # visibility period
        # (in which case the iflag value will change)
        if iflag != iflag_old:
            # Setting the iflip flag to True to keep track of the change
            # (in order to
            # detect CVZ object which are permanenetly visible)
            # If iflag = True we are starting a visibility period and use
            # a bisection method
            # to find the exact transition date. This assumes that there is
            # a single
            # transition in the interval => it looks like a step size of
            # 0.1 day is
            # sufficient to ensure that.
            if (iflag):
                step = currentdate-scanningStepSize
                wstart = ephemeris.bisect_by_FOR(currentdate, step, ra, dec)
            # IF iflag = False we are reaching the end of a visibility period.
            # Like for the previous case a bisection method is used to locate
            # accurately the end of the visibility period.
            else:
                step = currentdate-scanningStepSize
                wend = ephemeris.bisect_by_FOR(step, currentdate, ra, dec)
                startList.append(wstart)
                endList.append(wend)
                if (iflip):
                    statusList.append(0)
                else:
                    statusList.append(-1)
            iflip = True
            iflag_old = iflag

    # If there was a transition and we end up with a valid date, we close the
    # interval with the
    # end of the search interval
    if (iflag and iflip):
        startList.append(wstart)
        endList.append(currentdate)
        statusList.append(-2)
    # There is also the case were the visibility period covers the complete
    # search interval
    if (iflag and (not iflip)):
        startList.append(mjdmin)
        endList.append(mjdmax)
        statusList.append(1)
    # End of the function
    return startList, endList, statusList


def f_computeVisibilityPeriodsWithPA(ephemeris, mjdmin, mjdmax, ra, dec, pa):
    """
    # -----------------------------------------------------------
    # METHOD         f_computeVisibilityPeriodsWithPA()
    # TYPE           function
    #
    # DESCRIPTION    function that will compute the visibility
    #                periods for a given (RA,DEC), a given PA and
    #                over a given time period.
    #
    # SYNTAX    f_computeVisibilityPeriodsWithPA(ephemeris, mjdmin,
    #                    mjdmax, ra, dec, pa)
    #
    # ephemeris: input ephemeris object
    # mjdmin: beginning of the search interval (modified
    #    Julian date). It must be covered by the ephemeris.
    # mjdmax: end of the search interval (modified
    #    Julian date). It must be covered by the ephemeris.
    # ra: input RA coordinate (equatorial coordinate, in rad)
    # dec: input DEC coordinate (equatorial coordinate, in rad)
    # pa: input PA (in rad)
    #
    # Returns two lists containing the start end end of each
    # visibility period and a list containing a status flag:
    # flag = 0 visibility period fully in the search interval
    # flag = -1 start of the visibility period truncated by
    # the start of the search interval
    # flag = -2 end of the visibility period truncated by
    # the end of the search interval
    # flag = +1 the search interval is fully included in
    # the visibility period
    #
    # -----------------------------------------------------------
    """
    # ===========================================================
    # Paranoid checks
    # ===========================================================
    # print "# RA  = {:12.8f} rad = {:12.8f} deg".format(ra, ra / D2R)
    # print "# DEC = {:12.8f} rad = {:12.8f} deg".format(dec, dec / D2R)
    # print"# No constraint on the PA."
    if (ephemeris.amin > mjdmin):
        print("""f_computeVisibilityPeriodsWithPA(): the start of the search\
                 interval is not covered by the ephemeris.""")
        print("""Ephemeris start date (modified Julian date):\
                 {:8.5f}""".format(ephemeris.amin))
        print("""Search interval start date (modified Julian date):\
                 {:8.5f}""".format(mjdmin))
        raise ValueError
    if (ephemeris.amax < mjdmax):
        print("""f_computeVisibilityPeriodsWithPA(): the end of the search\
                 interval is not covered by the ephemeris.""")
        print("""Ephemeris end date (modified Julian date):\
                 {:8.5f}""".format(ephemeris.amax))
        print("""Search interval end date (modified Julian date):\
                 {:8.5f}""".format(mjdmax))
        raise ValueError
    # ===========================================================
    # Scanning the search period
    # ===========================================================
    # Flag used to track the beginning and the end of a
    # visibility period
    iflip = False
    wstart = mjdmin
    startList = []
    endList = []
    statusList = []
    # Scannning step size (must be small enough to make sure that
    # it cannot contain a full vsibility period (we would miss
    # it)
    scanningStepSize = 0.1
    span = int((mjdmax - mjdmin) / scanningStepSize)
    # Initialisation (first step of the scan is outside from the
    # loop
    iflag_old = ephemeris.is_valid(mjdmin, ra, dec, pa)
    for i in range(span):
        # Current date (the last step may be partial to remain
        # within the search interval
        currentdate = mjdmin + (i + 1) * scanningStepSize
        if (currentdate >= mjdmax):
            currentdate = mjdmax
        iflag = ephemeris.is_valid(currentdate, ra, dec, pa)
        # Checking if we are reaching the beginning or the end of a
        # visibility period
        # (in which case the iflag value will change)
        if iflag != iflag_old:
            # Setting the iflip flag to True to keep track of the change
            # (in order to
            # detect CVZ object which are permanenetly visible)
            # If iflag = True we are starting a visibility period and use
            # a bisection method
            # to find the exact transition date. This assumes that there
            # is a single
            # transition in the interval => it looks like a step size of
            # 0.1 day is
            # sufficient to ensure that.
            if (iflag):
                step = currentdate-scanningStepSize
                wstart = ephemeris.bisect_by_attitude(currentdate, step,
                                                      ra, dec, pa)
            # IF iflag = False we are reaching the end of a visibility period.
            # Like for the previous case a bisection method is used to locate
            # accurately the end of the visibility period.
            else:
                step = currentdate-scanningStepSize
                wend = ephemeris.bisect_by_attitude(step, currentdate,
                                                    ra, dec, pa)
                startList.append(wstart)
                endList.append(wend)
                if (iflip):
                    statusList.append(0)
                else:
                    statusList.append(-1)
            iflip = True
            iflag_old = iflag

    # If there was a transition and we end up with a valid date, we close
    # the interval with the
    # end of the search interval
    if (iflag and iflip):
        startList.append(wstart)
        endList.append(currentdate)
        statusList.append(-2)
    # There is also the case were the visibility period covers the
    # complete search interval
    if (iflag and (not iflip)):
        startList.append(mjdmin)
        endList.append(mjdmax)
        statusList.append(1)
    # End of the function
    return startList, endList, statusList


def f_computeDurationOfVisibilityPeriodWithPA(ephemeris,  mjdmin, mjdmax,
                                              ra, dec, pa, mjdc):
    """
    # -----------------------------------------------------------
    # METHOD         f_computeDurationOfVisibilityPeriodWithPA()
    # TYPE           function
    #
    # DESCRIPTION    function that will compute the duration of
    #                a specific visibility period associated to
    #                a given (RA,DEC), a given PA and given
    #                date.
    #
    # SYNTAX    f_computeDurationOfVisibilityPeriodWithPA(ephemeris,
    #                     mjdmin, mjdmax, ra, dec, pa, mjdc)
    #
    # ephemeris: input ephemeris object
    # mjdmin: beginning of the search interval (modified
    #    Julian date). It must be covered by the ephemeris.
    # mjdmax: end of the search interval (modified
    #    Julian date). It must be covered by the ephemeris.
    # ra: input RA coordinate (equatorial coordinate, in rad)
    # dec: input DEC coordinate (equatorial coordinate, in rad)
    # pa: input PA (in rad)
    # mjdc: date within the visibility period (i.e. compatible
    #    with (RA,DEC) and PA.
    #
    # Returns start,end,status
    # Status flag:
    # flag = 0 visibility period fully in the search interval
    # flag = -1 start of the visibility period truncated by
    # the start of the search interval
    # flag = -2 end of the visibility period truncated by
    # the end of the search interval
    # flag = +1 the search interval is fully included in
    # the visibility period
    #
    # -----------------------------------------------------------
    """
    # ===========================================================
    # Paranoid checks
    # ===========================================================
    # print "# RA  = {:12.8f} rad = {:12.8f} deg".format(ra, ra / D2R)
    # print "# DEC = {:12.8f} rad = {:12.8f} deg".format(dec, dec / D2R)
    # print"# No constraint on the PA."
    if (ephemeris.amin > mjdmin):
        print("""f_computeDurationOfVisibilityPeriodWithPA(): the start of\
                 thesearch interval is not covered by the ephemeris.""")
        print("""Ephemeris start date (modified Julian date):\
                 {:8.5f}""".format(ephemeris.amin))
        print("""Search interval start date (modified Julian date):\
                 {:8.5f}""".format(mjdmin))
        raise ValueError
    if (ephemeris.amax < mjdmax):
        print("""f_computeDurationOfVisibilityPeriodWithPA(): the end of the\
                 search interval is not covered by the ephemeris.""")
        print("""Ephemeris end date (modified Julian date):\
                 {:8.5f}""".format(ephemeris.amax))
        print("""Search interval end date (modified Julian date):\
                 {:8.5f}""".format(mjdmax))
        raise ValueError
    if (mjdmin > mjdc):
        print("""f_computeDurationOfVisibilityPeriodWithPA():\
                 initial date is not included in the search interval.""")
        print("""Search interval start date (modified Julian date):\
                 {:8.5f}""".format(mjdmin))
        print("Initial date (modified Julian date): {:8.5f}".format(mjdc))
        raise ValueError
    if (mjdmax < mjdc):
        print("""f_computeDurationOfVisibilityPeriodWithPA(): initial date is\
                 not included in the search interval.""")
        print("""Search interval end date (modified Julian date):\
                 {:8.5f}""".format(mjdmax))
        print("Initial date (modified Julian date): {:8.5f}".format(mjdc))
        raise ValueError

    iflag = ephemeris.is_valid(mjdc, ra, dec, pa)
    if (not iflag):
        print("""f_computeDurationOfVisibilityPeriodWithPA(): invalid date\
                 (not in a vsibility period).""")
        print("Date (modified Julian date): {:8.5f}".format(mjdc))
        raise ValueError

    # ===========================================================
    # Lookign for the start of the visibility period
    # ===========================================================
    scanningStepSize = 0.1
    iflipLeft = False
    currentmjd = mjdc
    continueFlag = True
    boundaryFlag = False
    while (continueFlag):
        currentmjd -= scanningStepSize
        if (currentmjd < mjdmin):
            currentmjd = mjdmin
            boundaryFlag = True
            continueFlag = False
        iflag = ephemeris.is_valid(currentmjd, ra, dec, pa)
        if (not iflag):
            wstart = ephemeris.bisect_by_attitude(currentmjd,
                                                  currentmjd+scanningStepSize,
                                                  ra, dec, pa)
            iflipLeft = True
            continueFlag = False
        elif (boundaryFlag):
            wstart = mjdmin

    iflipRight = False
    currentmjd = mjdc
    boundaryFlag = False
    continueFlag = True
    while (continueFlag):
        currentmjd += scanningStepSize
        if (currentmjd > mjdmax):
            currentmjd = mjdmax
            boundaryFlag = True
            continueFlag = False
        iflag = ephemeris.is_valid(currentmjd, ra, dec, pa)
        if (not iflag):
            wend = ephemeris.bisect_by_attitude(currentmjd-scanningStepSize,
                                                currentmjd, ra, dec, pa)
            iflipRight = True
            continueFlag = False
        elif (boundaryFlag):
            wend = mjdmax

    if ((not iflipLeft) and (not iflipRight)):
        status = 1
    elif (not iflipLeft):
        status = -1
    elif (not iflipRight):
        status = -2
    else:
        status = 0

    # End of the function
    return wstart, wend, status

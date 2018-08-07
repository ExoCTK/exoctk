"""Module containing library functions for time manipulation.
Standard for time representation in this project is fractional days.
Dates are represented as modified Julian dates (mjd).
An mjd gives the number of days since midnight on November 17, 1858.
"""
import string
from math import ceil, floor

# Constant for converting Julian dates to modified Julian dates
MJD_BASELINE = 2400000.5


def is_leap_year(year):
    """Returns True if the year is a leap year, False otherwise.

    Parameters
    ----------
    year: int
        The year to check.

    Returns
    -------
    bool
        Is the year a leap year?
    """
    return (((year % 4 == 0) and ((year % 100 > 0) or (year % 400 == 0))))


def days_in_year(year):
    """Returns the number of days in a year.

    Parameters
    ----------
    year: int
        The year to search.

    Returns
    -------
    days : int
        The number of days that year.
    """
    days = 365

    if (is_leap_year(year)):
        days += 1

    return(days)


def leap_years(year1, year2):
    """Returns the number of leap years between year1 and year2,
    non-inclusive.

    year1 and year2 must be integers, with year2 > year1

    Parameters
    ----------
    year1: int
        The start year.
    year2: int
        The end year.

    Returns
    -------
    int
        The number of leap years between year1 and year2.
    """

    # Find next years after year1 that are divisible by 4, 100, and 400
    next_div4 = int(4 * ceil(year1/4.0))
    next_div100 = int(100 * ceil(year1/100.0))
    next_div400 = int(400 * ceil(year1/400.0))

    # Now compute number of years between year1 and year2 that are
    # evenly divisible by 4, 100, 400
    div4_years = int(ceil((year2 - next_div4)/4.0))
    div100_years = int(ceil((year2 - next_div100)/100.0))
    div400_years = int(ceil((year2 - next_div400)/400.0))

    # Leap years are years divisible by 4, except for years
    # divisible by 100 that are not divisible by 400
    return(div4_years - (div100_years - div400_years))


def integer_days(time):
    """Takes a time in fractional days and returns integer component.

    Parameters
    ----------
    time: float
        The float time.

    Returns
    -------
    int
        The integer time.
    """
    # If time is negative, integer days is a larger negative number
    return(int(floor(time)))


def seconds_into_day(time):
    """Takes a time in fractional days and returns number of seconds since
    the start of the current day.

    Parameters
    ----------
    time: float
        The time as a float.

    Returns
    -------
    int
        The day's duration in seconds.
    """
    return(int(round(86400.0 * (time % 1))))


def days_to_seconds(days):
    """Takes a time in fractional days and converts it into integer
    seconds.

    Parameters
    ----------
    days: float
        The number of days as a float.

    Returns
    -------
    int
        The number of seconds in as many days.
    """
    return(int(round(86400 * days)))


def seconds_to_days(seconds):
    """Takes a time in integer seconds and converts it into fractional
    days.

    Parameters
    ----------
    seconds: int
        The number of seconds.

    Returns
    -------
    float
        The number of days as a float.
    """
    return(seconds / 86400.0)


def round_to_second(time):
    """Rounds a time in days to the nearest second.

    Parameters
    ----------
    time: int
        The number of days as a float.

    Returns
    -------
    float
        The number of seconds in as many days.
    """
    return(round(time * 86400)/86400.0)


def display_time(time, force_hours=False):
    """Returns a string representation of a time specified in fractional
    days.

    Parameters
    ----------
    time: float
        The time as a float.
    force_hours: bool
        Force the hour calculation.

    Returns
    -------
    str
        The time as a string.
    """
    # round to nearest second before extracting fields
    time = round_to_second(time)
    # if time is negative, print a minus sign and display absolute value
    if (time < 0):
        neg_string = '-'
        time = abs(time)
    else:
        neg_string = ''

    days = integer_days(time)
    day_string = hour_string = min_string = ''
    secs_within_day = seconds_into_day(time)
    hours_within_day = int(secs_within_day / 3600)
    secs_within_hour = secs_within_day % 3600
    mins_within_hour = int(secs_within_hour / 60)
    secs_within_min = secs_within_hour % 60

    # Unless force_hours is specified, only print a field if it or a
    # higher field is nonzero.  Fill with leading zeros.
    # The force_hours option is useful because Excel can get confused
    # when reading short time strings.
    if (days != 0):
        day_string = '%s:' % ((str(days)).zfill(3))

    if ((days != 0) or (hours_within_day > 0) or force_hours):
        hour_string = '%s:' % ((str(hours_within_day)).zfill(2))

    if ((days != 0) or (secs_within_day >= 60) or force_hours):
        min_string = '%s:' % ((str(mins_within_hour)).zfill(2))

    if ((days == 0) and (hours_within_day == 0) and (mins_within_hour == 0)):
        # avoid zero fill when there are only seconds
        sec_string = '%s' % ((str(secs_within_min)))
    else:
        sec_string = '%s' % ((str(secs_within_min)).zfill(2))

    return(neg_string + day_string + hour_string + min_string + sec_string)


def time_from_string(time_string):
    """Takes a string of the form ddd:hh:mm:ss and converts it to fractional
    days. All subfields above seconds are optional and may be omitted if the
    subfield and all higher-order ones are zero.

    Parameters
    ----------
    time_string: str
        The time as a string.

    Returns
    -------
    float
        The fractional days.
    """
    # extract fields
    fields = (string.split(time_string, ':'))
    seconds = int(fields[-1])
    num_fields = len(fields)
    # default to zero if not provided
    minutes = hours = days = 0

    if (num_fields > 1):
        minutes = int(fields[-2])
        if (num_fields > 2):
            hours = int(fields[-3])
            if (num_fields > 3):
                days = int(fields[-4])

    total_seconds = seconds + 60 * minutes + 3600 * hours + 86400 * days
    return(seconds_to_days(total_seconds))


def display_date(mjd):
    """Returns a string representation of the date represented by a
    modified Julian date.

    Parameters
    ----------
    mjd: float
        The modified julian day.

    Returns
    -------
    str
        The MJD as a string.
    """
    # adjust to number of days since Dec. 31, 1857
    int_days = int(floor(321.0 + mjd))
    # seconds_in_day = seconds_into_day(mjd)
    fractional_day = mjd % 1

    # First compute year and day without allowing for leap years, then adjust
    year = 1858 + int_days/365
    day_of_year = int_days % 365 - leap_years(1858, year)
    # handle case where leap year adjustment has made day negative
    while (day_of_year < 1):
        year -= 1
        day_of_year = day_of_year + days_in_year(year)

    year_string = '%s:' % (year)
    return(year_string + display_time(day_of_year + fractional_day))


def compute_mjd(year, day_of_year, hour, minute, second):
    """Computes a modified Julian date from a date specified as a year,
    day of year, hour, minute, and second.
    Arguments should be integers.

    Parameters
    ----------
    year: int
        The year.
    day_of_year: int
        The day.
    hour: int
        The hour.
    minute: int
        The minute.
    second: int
        The second.

    Returns
    -------
    float
        The modified julian day.
    """
    fractional_days = (hour * 3600 + minute * 60 + second)/86400.0
    mjd_years = year - 1859
    num_leaps = leap_years(1858, year)  # number of leap years since 1858

    # Add 45 days from Nov. 17 to end of 1858
    return((365*mjd_years)+num_leaps+45+(day_of_year-1)+fractional_days)


def mjd_from_string(time_string):
    """Takes a string of the form yyyy.ddd:hh:mm:ss and returns an mjd.

    Parameters
    ----------
    time_string: str
        The MJD as a string.

    Returns
    -------
    float
        The modified julian day.
    """
    years = int(time_string[0:4])
    days = int(time_string[5:8])
    hours = int(time_string[9:11])
    minutes = int(time_string[12:14])
    seconds = int(time_string[15:17])

    return(compute_mjd(years, days, hours, minutes, seconds))


def mjd_to_jd(mjd):
    """Converts a modified Julian date to a true Julian date.

    Parameters
    ----------
    mjd: float
        The modified julian day.

    Returns
    -------
    float
        The true Julian day.
    """
    return(MJD_BASELINE + mjd)


def jd_to_mjd(jd):
    """Converts a Julian date to a modified Julian date.

    Parameters
    ----------
    jd: float
        The true Julian day.

    Returns
    -------
    float
        The modified Julian day.
    """
    return (jd - MJD_BASELINE)


class Interval(object):
    """Class to represent a simple temporal interval.
    """
    def __init__(self, start, end):
        """Constructor for an interval.

        Parameters
        ----------
        start: float
            The start time.
        end: float
            The end time.
        """
        self.start = start
        self.end = end

    def __str__(self):
        """Returns a string representation of the interval."""
        return('Interval: start: %s, end: %s' % (display_date(self.start),
                                                 display_date(self.end)))

    def start_time(self):
        """Returns the start of the interval."""
        return(self.start)

    def end_time(self):
        """Returns the end of the interval."""
        return(self.end)

    def duration(self):
        """Returns the duration of an interval in fractional days."""
        return(self.end_time() - self.start_time())

    def temporal_relationship(self, time):
        """Returns the temporal relationship between an interval and an
        absolute time.

        Returns 'before' if the interval ends at or before the time,
        'after' if the interval begins at or after the time,
        'includes' if the time occurs during the interval.

        Parameters
        ----------
        time: float
            The time.

        Returns
        -------
        rel : str
            The temporal relationship.
        """
        if (self.end_time() <= time):
            rel = 'before'
        elif (self.start_time() >= time):
            rel = 'after'
        else:
            rel = 'includes'

        return(rel)


class FlexibleInterval(Interval):
    """Class to represent an interval with flexibility on when it can
    start and end.
    """
    def __init__(self, est, lst, let):
        """Constructor for a FlexibileInterval.

        Parameters
        ----------
        est: float
            Earliest start time (mjd).
        lst: float
            Latest start time (mjd).
        let: float
            Latest end time (mjd).
        """
        self.est = est
        self.lst = lst
        self.let = let

    def __str__(self):
        """Returns a string representation of the FlexibleInterval."""
        txt = (display_date(self.est), display_date(self.lst),
               display_date(self.let))
        return('FlexibleInterval: EST: %s, LST: %s, LET: %s' % txt)

    def start_time(self):
        """Returns the start of the FlexibleInterval."""

        return(self.est)

    def end_time(self):
        """Returns the end of the FlexibleInterval."""

        return(self.let)

    def flexibility(self):
        """Returns the flexibility of the FlexibleInterval, in
        fractional days."""

        return(self.lst - self.est)

    def maximum_duration(self):
        """Returns the maximum duration of the FlexibleInterval, in
        fractional days."""

        return(self.let - self.lst)

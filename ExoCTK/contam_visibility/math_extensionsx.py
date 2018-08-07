"""This module provides simple extensions to the Python mathematical
library.
2018/06/26 Made PEP8 compliant and added
sind() and cosd() - Joe Filippazzo
Version 1 August 23, 2010 RLH - Added OBLIQUITY.
Version 0 August 6, 2010 RLH  - Created
"""
from copy import deepcopy
from math import radians, acos, asin, cos, sin, sqrt, pi, ceil, exp, atan2

R2D = 180.0/pi
D2R = 1/R2D
EPSILON = 1.0e-10   # Use for safe comparisons of floating-point numbers
OBLIQUITY = 23.43929 * D2R  # Obliquity of Earth's orbit, in radians


def sind(x):
    """Return the sin in degrees.

    Parameters
    ----------
    x: float
        The evaluand.

    Returns
    -------
    float
        The sin of x in degrees.
    """
    return sin(radians(x))


def cosd(x):
    """Return the cos in degrees

    Parameters
    ----------
    x: float
        The evaluand.

    Returns
    -------
    float
        The cos of x in degrees.
    """
    return cos(radians(x))


def atan2d(x):
    """Return the arctan in degrees

    Parameters
    ----------
    x: float
        The evaluand.

    Returns
    -------
    float
        The arctan of x in degrees.
    """
    return atan2(radians(x))


def really_less_than(x, y):
    """Safe less-than function that returns true if and only if x is
    "significantly" less than y.

    Parameters
    ----------
    x: float
        The first number.
    y: float
        The second number.

    Returns
    -------
    bool
        True if x is less, else False.
    """
    return(x < y - EPSILON)


def really_greater_than(x, y):
    """Safe greater-than function that returns true if and only if x is
    "significantly" greater than y

    Parameters
    ----------
    x: float
        The first number.
    y: float
        The second number.

    Returns
    -------
    bool
        True if x is greater, else False.
    """
    return(x > y + EPSILON)


def asin2(val):
    """Safe version of asin that handles invalid arguments.

    Arguments greater than 1 are truncated to 1; arguments less than -1 are
    set to -1.

    Parameters
    ----------
    val: float
        The evaluand.

    Returns
    -------
    float
        The arcsin of the value.
    """
    return(asin(max(-1.0, min(1.0, val))))


def acos2(val):
    """Safe version of acos that handles invalid arguments in the same way as
    asin2

    Parameters
    ----------
    val: float
        The evaluand.

    Returns
    -------
    float
        The arccos of the value.
    """
    return(acos(max(-1.0, min(1.0, val))))


def avg(l):
    """Returns the average of a list of numbers

    Parameters
    ----------
    l: sequence
        The list of numbers.

    Returns
    -------
    float
        The average.
    """
    return(sum(l) / float(len(l)))


def avg2(num1, num2):
    """Returns the average of two numbers

    Parameters
    ----------
    num1: float
        The first number.
    num2: float
        The second number.

    Returns
    -------
    float
        The average.
    """
    return((num1 + num2)/2.0)


def output_as_percentage(num, fractional_digits=1):
    """Output a percentage neatly.

    Parameters
    ----------
    num: float
        The number to make into a percentage.
    fractional_digits: int
        Number of digits to output as fractions of a percent.
        If not supplied, there is no reduction in precision.

    Returns
    -------
    str
        The percentage.
    """
    if (fractional_digits is not None):
        format_str = '%%.%.df' % (fractional_digits)  # creates format string
    else:
        format_str = '%f'

    return('%s%%' % (format_str % (num)))


def percent_str(num, fractional_digits=1):
    """Output a number as a percentage.

    Parameters
    ----------
    num: float
        The number to make into a percentage.
    fractional_digits: int
        Number of digits to output as fractions of a percent.
        If not supplied, there is no reduction in precision.

    Returns
    -------
    str
        The percentage.
    """
    return(output_as_percentage(100 * num, fractional_digits))


def variance(l):
    """Variance of a list of numbers that represent sample values

    Parameters
    ----------
    l: sequence
        The list to take the variance of.

    Returns
    -------
    float
        The variance.
    """
    # Returns the sample variance (n-1 formula).
    mean = avg(l)
    sumsq = sum([(i - mean)**2 for i in l])
    return(sumsq / float(len(l) - 1))


def stdev(l):
    """Standard deviation of a list of numbers that represent sample values

    Parameters
    ----------
    l: sequence
        The list to take the standard deviation of.

    Returns
    -------
    float
        The standard deviation.
    """
    # Simply take the square root of the sample variance.
    return(sqrt(variance(l)))


def factorial(num):
    """Returns the factorial of a nonnegative integer.
    This function is provided in math module starting with Python 2.6,
    but implement anyway for compatibility with older systems.

    Parameters
    ----------
    num: int
        The number to factorialize.

    Returns
    -------
    int
        The factorial.
    """
    result = 1   # factorial of 0 is defined as 1

    for i in range(2, num + 1):
        result = result * i
    return(result)


def conditional_probability(p_joint, p_B):
    """Returns probability of event A given event B.

    Parameters
    ----------
    p_joint: float
        P(A,B).
    p_B: float
        Probability of event B.

    Returns
    -------
    float
        The probability.
    """
    return(p_joint / p_B)


class Polynomial(object):
    """Class to represent a polynomial."""
    def __init__(self, coefficients):
        """Constructor for a polynomial.
        Coefficients = a list of coefficients, starting with order 0 and
        increasing.

        Parameters
        ----------
        coefficients: sequence
            The list of coefficients, starting with order 0 and increasing.
        """
        self.coefficients = coefficients

    def __str__(self):
        """Returns a string representation of a polynomial."""
        order = len(self.coefficients) - 1
        return_string = 'Polynomial: order % d' % (order)

        for index in range(order + 1):  # print list of coefficients
            c = index, self.coefficients[index]
            return_string = return_string + '\nCoefficient % d = % .3f' % c

        return(return_string)

    def apply(self, value):
        """Returns the result of applying a polynomial to an input value

        Parameters
        ----------
        value: float
            The evaluand.

        Returns
        -------
        float
            The result of the evaluated equation.
        """
        result = 0

        # For each index, raise the input to the power given by that index
        # and multiply by the coefficient
        for index in range(len(self.coefficients)):
            result = result + self.coefficients[index] * value**index

        return (result)


class LinearEquation(Polynomial):
    """Subclass of Polynomial for linear equations.
    This implementation is three times faster, so Polynomial should be reserved
    for higher orders."""
    def __init__(self, coeff0, coeff1):
        """Constructor for a linear equation to provide a more 'natural'
        interface without using a list.

        Parameters
        ----------
        coeff0: float
            Additive constant.
        coeff1: float
            Multiplicative coefficient.
        """
        self.coefficients = [coeff0, coeff1]

    def apply(self, value):
        """Applies a linear equation to an input value.

        This is intended to be faster than the more general method with
        Polynomial.

        Parameters
        ----------
        value: float
            The evaluand.

        Returns
        -------
        float
            The result of the evaluated equation.
        """
        return(value * self.coefficients[1] + self.coefficients[0])


class HistogramBin(object):
    """Class to represent a bin within a histogram."""
    def store_items(self, num_items=1):
        """Stores a given number of items in the bin.

        Parameters
        ----------
        num_items: int
            Number of items to store (default 1).
        """
        self.count += num_items


class DiscreteBin(HistogramBin):
    """Class to represent a bin with a fixed value."""
    def __init__(self, bin_value):
        """Constructor for a fixed-value bin

        Parameters
        ----------
        bin_value: float
            The value of the bin.
        """
        self.bin_value = bin_value
        self.count = 0

    def __str__(self):
        """Returns a printed representation of the bin.
        """
        # Don't assume a type for the count, as it may not be an integer
        # if normalized.
        return('%s: % s items' % (self.bin_value, self.count))

    def ismatch(self, value):
        """Returns True if the value matches the bin, False otherwise

        Parameters
        ----------
        value: float
            The value to compare.

        Returns
        -------
        bool
            True if matches, else False.
        """
        return(value == self.bin_value)


class RangeBin(HistogramBin):
    """Class to represent a bin with a range."""
    def __init__(self, min_value=None, max_value=None, lower_inclusive=False,
                 upper_inclusive=True):
        """Constructor for a range bin.

        Parameters
        ----------
        min_value: float
            Minimum value for the bin.
        max_value: float
            Maximum value for the bin.
        lower_inclusive: bool
            True if min_value is inclusive, else False (default).
        upper_inclusive: bool
            True if max_value is inclusive, else False (default).
        """
        self.min_value = min_value
        self.max_value = max_value
        self.lower_inclusive = lower_inclusive
        self.upper_inclusive = upper_inclusive
        self.count = 0

    def describe_limits(self, precision=2):
        """Returns a printed representation of the limits of the bin.

        Parameters
        ----------
        precision: int
            Number of digits to print after the decimal point.

        Returns
        -------
        str
            The limits of the bin.
        """
        if(self.min_value is None):
            if(self.upper_inclusive):
                result = '<= % . * f:' % (precision, self.max_value)
            else:
                result = '< % . * f:' % (precision, self.max_value)
        elif(self.max_value is None):
            if(self.lower_inclusive):
                result = '>= % . * f:' % (precision, self.min_value)
            else:
                result = '> % . * f:' % (precision, self.min_value)
        else:
            v = (precision, self.min_value, precision, self.max_value)
            result = '%. * f to % . * f:' % v

        return(result)

    def __str__(self):
        """Returns a printed representation of the bin
        """
        # Don't assume a type for the count, as it may not be an integer
        # if normalized.
        return('%s % s items' % (self.describe_limits(), self.count))

    def istoo_high(self, value):
        """Returns True if the specified value is too high for the bin.
        Assumes the bin has an upper limit.

        Parameters
        ----------
        value: float
            The value to compare.

        Returns
        -------
        bool
            True if too high, else False.
        """
        # value equal to limit not consideredtoo high
        if(self.upper_inclusive):
            result = (value > self.max_value)
        else:
            result = (value >= self.max_value)

        return(result)

    def ismatch(self, value):
        """Indicates whether the bin matches the value.

        Parameters
        ----------
        value: float
            The value to compare.

        Returns
        -------
        bool
            True if matches, else False.
        """
        return(False)   # not really applicable, so always fail


class Histogram(object):
    """Class to represent a histogram."""

    def retrieve_count(self, bin_index):
        """Returns the number of items stored in a given bin of the histogram.

        Parameters
        ----------
        bin_index: int
            The index to use (starts with 1).

        Returns
        -------
        int
            The number of items in the bin.
        """
        return(self.bins[bin_index - 1].count)

    def num_items(self):
        """Returns the total number of items stored in the histogram
        """
        return(sum([bin.count for bin in self.bins]))

    def __str__(self):
        """Returns a printed representation of the histogram.
        """
        # Don't assume a type for the total, as it may not be an integer
        # if normalized.
        v = (len(self.bins), self.num_items())
        result = 'Histogram: % d bins, % s items\n' % v

        for bin in self.bins:
            result = result + '%s\n' % (bin.__str__())

        return(result)

    def normalize(self, total=None):
        """Takes a histogram and returns a new histogram that normalizes all
        its values.

        Parameters
        ----------
        total: int
            Number of items to divide each bin by for the normalization.
            If not supplied, it defaults to the total in the histogram.

        Returns
        -------
        new_histogram : Histogram
            The new histogram.
        """
        if (total is None):
            total = self.num_items()

        new_histogram = deepcopy(self)   # make a deep copy

        for bin in new_histogram.bins:
            bin.count = (float(bin.count)) / total

        return(new_histogram)


class DiscreteHistogram(Histogram):
    """Class to represent a histogram with discrete values."""

    def __init__(self, values):
        """Initializes a histogram with discrete values.

        Parameters
        ----------
        values: sequence
            List of the discrete values.
        """
        self.bins = []

        # Create a DiscreteBin object for each value and add it to the list.
        for value in values:
            self.bins.append(DiscreteBin(value))

    def retrieve_values(self):
        """Returns the list of bin values of a discrete histogram
        """
        return([bin.bin_value for bin in self.bins])

    def retrieve_count_by_value(self, value):
        """Returns the count matching a certain value.  If not found,
        return None

        Parameters
        ----------
        value: float
            The value to retrieve.
        """
        bin_index = 0
        result = None

        while ((not result) and (bin_index < len(self.bins))):
            bin = self.bins[bin_index]

            if (bin.ismatch(value)):
                result = bin.count

            bin_index += 1

        return(result)

    def store_items(self, value, count=1):
        """Stores a value in the discrete histogram if it matches one of
        the bin values.

        Count = number of items with that value to store (default 1).

        Returns True if a match was found and the value could be stored,
        False otherwise

        Parameters
        ----------
        value: float
            The value to store.
        count: int
            Number of items with that value to store (default 1).
        
        Returns
        -------
        found : bool
            Whether or not a match was found.
        """
        bin_index = 0
        found = False

        while ((not found) and (bin_index < len(self.bins))):
            bin = self.bins[bin_index]

            if (bin.ismatch(value)):
                bin.store_items(count)
                found = True

            bin_index += 1

        return(found)


class ContinuousHistogram(Histogram):
    """Class to represent a histogram with continuous values."""

    def __init__(self, boundaries, highest_inclusive=False):
        """Initializes a continuous histogram.

        Default behavior with highest_inclusive = False:
           Bin 0 is defined by x <= boundaries[0].
           For i > 0, bin i is defined by boundaries[i-1] < x <= boundaries[i].
           Bin n+1 is defined by x > boundaries[n-1].

        Behavior with highest_exclusive = True:
           Bins below n are defined in the same way as above.
           Bin n is defined by boundaries[n-2] < x < boundaries[n-1].
           Bin n+1 is defined by x >= boundaries[n-1].

        Parameters
        ----------
        boundaries: sequence
            List of numbers that separate the bins, in increasing order.
        highest_inclusive: bool
            True if highest bin includes the last boundary,
            False (default) otherwise.
        """
        self.bins = []
        self.highest_inclusive = highest_inclusive
        lower_lim = None

        # A histogram is a list of bins.
        # For the first bin, the lower limit is None (unbounded).
        # For the last bin, the upper limit is None.
        # Rely on default limits behavior in HistogramBin constructor.

        for index in range(len(boundaries)):
            upper_lim = boundaries[index]
            self.bins.append(RangeBin(min_value=lower_lim,
                                      max_value=upper_lim))
            lower_lim = upper_lim

        self.bins.append(RangeBin(min_value=lower_lim))   # Add highest bin

        # If highest_inclusive is set, change limit behavior of two
        # highest bins.
        if (highest_inclusive):
            self.bins[-2].upper_inclusive = False
            self.bins[-1].lower_inclusive = True

    def retrieve_boundaries(self):
        """Returns the list of boundaries of a continuous histogram.
        """
        # Return upper limits of all bins except the last.
        return([bin.max_value for bin in self.bins[:-1]])

    def store_items(self, value, count=1):
        """Stores a value in the continuous histogram.

        Parameters
        ----------
        value: float
            The value to store.
        count: int
            Number of items with that value to store (default 1).
        """
        bin_index = 0
        found = False

        # Search all bins up to the next-highest.  If value is not too high
        # for the bin, store the items there.
        while ((not found) and (bin_index < len(self.bins) - 1)):
            bin = self.bins[bin_index]

            if (not (bin.istoo_high(value))):
                found = True
                bin.store_items(count)

            bin_index += 1

        # The value must belong in the highest bin if not found earlier.
        if (not found):
            self.bins[-1].store_items(count)


def combine_histograms(histograms):
    """Takes a list of histograms and returns a new Histogram object that sums
    the values in each bin.

    All histograms in the list must be identical except for the count.

    Parameters
    ----------
    histograms: sequence
        A lst of Histogram objects to combine.

    Returns
    -------
    Histogram
        The combined histogram.
    """
    # Initialize the new histogram with properties of the first histogram
    # in the list.
    if (isinstance(histograms[0], ContinuousHistogram)):
        hist_bounds = histograms[0].retrieve_boundaries()
        hist_high = histograms[0].highest_inclusive
        new_histogram = ContinuousHistogram(hist_bounds, hist_high)
    else:  # assume discrete histogram
        new_histogram = DiscreteHistogram(histograms[0].retrieve_values())

    for bin_index in range(len(new_histogram.bins)):
        total_items = sum([hist.bins[bin_index].count for hist in histograms])
        new_histogram.bins[bin_index].store_items(total_items)

    return(new_histogram)


def average_histograms(histograms):
    """Takes a list of histogram objects and simply averages all the bin values.

    All histograms in the list must be identical except for the count.

    Parameters
    ----------
    histograms: sequence
        A lst of Histogram objects to combine.

    Returns
    -------
    new_histogram : Histogram
        The averaged histogram.
    """
    # Make a copy of the first histogram in the list.
    new_histogram = deepcopy(histograms[0])

    for bin_index in range(len(new_histogram.bins)):
        new_histogram.bins[bin_index].count = avg([hist.bins[bin_index].count
                                                   for hist in histograms])

    return(new_histogram)


class PoissonDistribution(DiscreteHistogram):
    """Class to represent a Poisson distribution."""

    def probability(self, k):
        """Computes the probability that the Poisson distribution takes on
        the value k.

        Value must be a nonnegative integer.

        Parameters
        ----------
        k: float
            The value to compute.

        Returns
        -------
        float
            The probability.
        """
        u = self.mean

        return((u**k * exp(-u)) / factorial(k))

    def generate_distribution(self):
        """Populates a Poisson distribution up to the maximum bin.
        """
        cum = 0

        # For each value from 0 up to the next-to-last, compute the
        # Poisson distribution for that value, populate the bin, and keep
        # track of the cumulative total.
        for value in range(len(self.bins) - 1):
            self.bins[value].count = p = self.probability(value)
            cum += p

        # remainder of distribution goes in the last bin
        self.bins[-1].count = 1.0 - cum

    def __init__(self, mean, max_boundary):
        """Constructor function for the Poisson distribution.

        Parameters
        ----------
        mean: float
            Mean parameter for the Poisson distribution.
        max_boundary: float
            The largest parameter for which the probability is to
            be computed. All values larger than max_boundary will be lumped
            into the highest bin.
        """
        self.mean = mean
        self.bins = []

        # Create discrete bins for values from 0 to max_boundary, inclusive.
        for value in range(max_boundary + 1):
            self.bins.append(DiscreteBin(value))

        # Use a RangeBin object for the highest bin.
        self.bins.append(RangeBin(min_value=max_boundary))

        # Now generate the distribution.
        self.generate_distribution()

    def __str__(self):
        """Inspector function for Poisson distribution."""
        poisson_info = 'PoissonDistribution: Mean: % .2f\n' % (self.mean)
        generic_info = super(self.__class__, self).__str__()

        return(poisson_info + generic_info)

    def retrieve_values(self):
        """Returns the list of bin values for the Poisson distribution."""
        return(list(range(len(self.bins) - 1)))   # leave out last bin

    def retrieve_count_by_value(self, value):
        """Returns the number of items in the histogram that have the
        designated value.

        Value must be an integer between 0 and max_boundary.

        Parameters
        ----------
        value: float
            The value to retrieve.

        Returns
        -------
        result : int
            The number of items with the given value.
        """
        # If the value exceeds max_boundary, return None.  Otherwise just
        # use the value as an index into the list of bins.
        if (value > len(self.bins) - 2):
            result = None
        else:
            result = self.bins[value].count

        return(result)

    def cumulative_probability(self, value):
        """Returns the probability that a random variable will have a value no
        greater than the one specified.

        Parameters
        ----------
        value: float
            The value between 0 and the max_boundary of the distribution.

        Returns
        -------
        float
            The probability.
        """
        return(sum([self.bins[i].count for i in range(value + 1)]))


class StatisticalList(list):
    """Numeric list class with statistical attributes."""

    def __init__(self, data=None):
        """Initializes a statistical list.

        Parameters
        ----------
        data: sequence
            List of inputs to list.
        """
        # If data were provided, copy into the list.
        if (data is not None):
            for i in range(len(data)):
                self.append(data[i])

    def compute_variance(self):
        """Computes the variance of a statistical list.
        """
        mean = self.mean

        # Variance is defined as the sum of the squares of the differences
        # between each data point and the mean, divided by the number of
        # degrees of freedom.  For now assume the list is a sample, so
        # DOF = n-1.
        return(sum([(n - mean)**2 for n in self]) / (len(self) - 1))

    def compute_rms(self):
        """Computes the rms value of a statistical list.
        """
        # Simply sum the squares, divide by n, and take the square root.
        return(sqrt(avg([n**2 for n in self])))

    def compute_statistics(self, min_value=None, max_value=None,
                           max_bins=None):
        """Computes statistics for a StatisticalList object; must contain at
        least one element.

        Parameters
        ----------
        min_value: float
            Minimum value for cutoff of histogram (defaults to
            minimum in list).
        max_value: float
            Maximum value for cutoff of histogram (defaults to
            maximum in list).
        max_bins: int
            Maximum number of bins in histogram.
        """
        # first sort the list in increasing order -- note this is destructive
        self.sort()

        # Compute min, max, mean, median, variance, standard deviation,
        # and rms value.
        num_elements = len(self)
        self.min = self[0]
        self.max = self[-1]
        self.mean = avg(self)
        self.median = self[(num_elements/2)]
        self.variance = self.compute_variance()
        self.stdev = sqrt(self.variance)
        self.rms_value = self.compute_rms()

        # Create a list of percentiles at 2% steps from 0 to 100%.
        # At each percentile, find the nearest index and add an entry.
        # Each entry is a tuple with (<percentile>, <value>).
        self.percentiles = []

        for i in range(51):
            index = int(round((num_elements - 1) * (i/50.0)))
            self.percentiles.append((i * 2, self[index]))

        # Define histogram parameters.
        # Number of bins defaults to one-fourth the number of elements.
        # If max_bins is specified, limit number of bins accordingly.
        num_bins = max(3, int(ceil(num_elements/4.0)))   # minimum of 3 bins

        if (max_bins is not None):
            num_bins = min(num_bins, max_bins)

        bin_step = (self.max - self.min) / (float(num_bins))
        min_bin_value = self.min + bin_step
        max_bin_value = self.max - bin_step

        # If min_value and/or max_value is specified, limit min/max bin
        # values accordingly.
        if (min_value is not None):
            min_bin_value = min_value

        if (max_value is not None):
            max_bin_value = max_value

        # recalculate bin step
        bin_step = (max_bin_value - min_bin_value) / (float(num_bins - 2))

        # Create a histogram with the specified number of evenly spaced
        # bounds. Number of bounds is equal to number of bins - 1.
        boundary_value = min_bin_value
        bounds = []

        for i in range(num_bins - 1):
            bounds.append(boundary_value)
            boundary_value = boundary_value + bin_step

        self.histogram = ContinuousHistogram(bounds, highest_inclusive=False)

        # Now store the data in the histogram
        for value in self:
            self.histogram.store_items(value)

    def __str__(self):
        """Prints data on a statistical list after statistics are generated."""
        vals = [len(self), self.min, self.max, self.mean, self.median]
        vals += [self.variance, self.stdev, self.rms_value]
        string = """StatisticalList: % d elements, min = % .4f, max = % .4f,
                    mean = % .4f, median = % .4f, variance = % .4f,
                    stdev = % .4f, rms = % .4f""" % vals

        return string


class Circle(object):
    """Class to represent a circle."""

    def __init__(self, radius):
        """Initialize a circle with a specified radius.

        Parameters
        ----------
        radius: float
            The radius of the circle.
        """
        self.radius = radius

    def __str__(self):
        """Inspector method for the circle.
        """
        return('Circle: radius = % .2f' % (self.radius))

    def area(self):
        """Returns the area of the circle.
        """
        return(pi * self.radius**2)


class Rectangle(object):
    """Class to represent a rectangle."""

    def __init__(self, length, width):
        """Initialize a rectangle with a specified length and width.

        Parameters
        ----------
        length: float
            The length of the rectangle.
        width: float
            The width of the rectangle.
        """
        self.length = length
        self.width = width

    def __str__(self):
        """Inspector method for the rectangle.
        """
        dims = (self.length, self.width)
        return('Rectangle: length = % .2f, width = % .2f' % dims)

    def area(self):
        """Returns the area of the rectangle.
        """
        return(self.length * self.width)

    def motion_tolerant_area(self, motion_length, motion_angle):
        """Returns the area within a rectangle that can tolerate a motion in
        a known direction while remaining within the rectangle.

        Parameters
        ----------
        motion_length: float
            Distance of motion (same units as rectangle length and width).
        motion_angle: float
            Angle in radians between the direction of motion and long
            direction of rectangle.

        Returns
        -------
        float
            The area.
        """
        # Compute the x and y distances to the edge.  A position
        # within the rectangle is only motion-tolerant if it exceeds
        # both of these distances.  The effective length is thus reduced
        # by delta_x and the effective width by delta_y.
        delta_x = motion_length * cos(motion_angle)  # lengthwise direction
        delta_y = motion_length * sin(motion_angle)
        return ((self.length - delta_x) * (self.width - delta_y))


class Square(Rectangle):
    """Class to represent a square."""
    def __init__(self, side):
        """Initialize a square with a specified side length.

        Parameters
        ----------
        side: float
            The length of the square side.
        """
        self.side = side

        # call superclass method with length and width
        super(self.__class__, self).__init__(side, side)

    def __str__(self):
        """Inspector method for the square.
        """
        return('Square: side = % .2f' % (self.side))

    def inner_area(self, excluded_width):
        """Returns the area of the square after removing a strip of specified
        width along each edge.

        Parameters
        ----------
        excluded_width: float
            The width of the strip to remove.

        Returns
        -------
        float
            The area.
        """
        # Return the area of a square that is reduced in side length by twice
        # the specified width, because both sides are affected.
        return(Square(self.side - 2 * excluded_width).area())

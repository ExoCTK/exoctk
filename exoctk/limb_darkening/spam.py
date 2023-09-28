"""
Module to calculate SPAM coefficients in LDC tool

Author: Nestor Espinoza
"""
import numpy as np
from scipy.optimize import minimize
from batman import TransitModel, TransitParams

from .. import utils


def init_batman(t, ld_law, nresampling=None, etresampling=None):
    """  
    This function initializes the batman lightcurve generator object.

    Parameters
    ----------

    t: array
      Array containing the times at which the lightcurve will be evaluated. Assumes units of days.
    ld_law: string
      Limb-darkening law used to compute the model. Available ld laws: uniform, linear, quadratic, 
      logarithmic, exponential, squareroot, nonlinear, power2.
    nresampling: int
      Number of resampled points in case resampling is wanted.
    etresampling: float
      Exposure time of the resampling, same units as input time.

    Returns
    -------
    
    params: batman object
      Object containing the parameters of the lightcurve model.
    m: batman object
      Object that enables the lightcurve calculation.

    """

    params = TransitParams()
    params.t0 = 0. 
    params.per = 1. 
    params.rp = 0.1
    params.a = 15.
    params.inc = 87.
    params.ecc = 0. 
    params.w = 90.

    if ld_law == 'linear':
        params.u = [0.5]

    elif ld_law == 'nonlinear':
        params.u = [0.1, 0.1, 0.1, 0.1]

    else:       
        params.u = [0.1,0.3]

    params.limb_dark = ld_law

    if nresampling is None or etresampling is None:
        m = TransitModel(params, t)

    else:
        m = TransitModel(params, t, supersample_factor=nresampling, exp_time=etresampling)

    return params, m


def spam_objective_function(theta, params, m, params_twop, m_twop):
    """
    Objective function that the SPAM algorithm is looking to minimize. 
    """

    u1, u2 = theta
    params_twop.u = [u1, u2]

    return np.sum((m.light_curve(params) - m_twop.light_curve(params_twop))**2)


def transform_coefficients(c1, c2, c3, c4, planet_name=None, planet_data=None, ld_law='quadratic', ndatapoints=1000, method='BFGS', u1_guess=0.5, u2_guess=0.5):
    """
    Given a set of non-linear limb-darkening coefficients (c1, c2, c3 and c4) and either a planet name ('planet_name') or a dictionary with 
    the planet's data ('planet_data'), this function returns the Synthetic-Photometry/Atmosphere-Model (SPAM; https://arxiv.org/abs/1106.4659) 
    limb-darkening coefficients of the star. 

    Reference:
        Howarth, I., 2011, "On stellar limb darkening and exoplanetary transits", MNRAS, 418, 1165
        https://ui.adsabs.harvard.edu/abs/2011MNRAS.418.1165H/abstract

    Parameters
    ----------

    c1: float
      First limb-darkening coefficient of the non-linear law.
    c2: float
      Same as c1, but second.
    c3: float
      Same as c1, but third.
    c4: float
      Same as c1, but fourth.
    planet_name: string
      String with the name of the input planet (e.g., 'WASP-19b'); this will be used to query the planet properties from MAST.
    planet_data: dict
      Dictionary containing the planet properties. In particular, this dictionary should contain the keys 
      'transit_duration', 'orbital_period' (days), 'Rp/Rs', 'a/Rs', 'inclination' (degrees), 'eccentricity' and 'omega' (degrees) 
      for the algorithm to work. Properties in this dictionary take prescedence over the ones retrieved from MAST.
    ld_law: string
      Limb-darkening law for which SPAM coefficients are wanted. Default is 'quadratic', but can also be 'squareroot' or 'logarithmic'. 
    ndatapoints: int
      Number of datapoints that will be used for the lightcurve simulations to extract back the SPAM coefficients.
    method: string
      Minimization method to match lightcurves. Default is 'BFGS', but can be any of the ones available for scipy.optimize.minimize.
      Details: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
    u1_guess: float
      Guess starting value for u1. Default is 0.5
    u2_guess: flat
      Guess starting value for u2. Default is 0.5.

    Returns
    -------
    
    u1: float
      The first limb-darkening coefficient of the selected two-parameter law.
    u2: flat
      Same as u1, but for the second coefficient.

    """

    planet_properties = ['transit_duration', 'orbital_period', 'Rp/Rs', 'a/Rs', 'inclination', 'eccentricity', 'omega']

    # Check if planet name is given
    if planet_name is not None:

        # If planet_name is given, retrieve MAST properties
        mast_planet_data, url = utils.get_target_data(planet_name)

        # Merge the dictionaries, prioritizing the manual input
        mast_planet_data.update(planet_data or {})
        planet_data = mast_planet_data

    if planet_data is None:

        raise Exception("User must input either 'planet_name' and/or 'planet_data' for SPAM to work. See details by doing exoctk.limb_darkening.spam.transform_coefficients?.")

    # Check that all properties exist in the input dictionary
    missing = [planet_property for planet_property in planet_properties if planet_data.get(planet_property) is None]
    if len(missing) > 0:

        # Print current data for user
        print('{} properties'.format(planet_name or 'Planet'))
        for planet_property in planet_properties:
            print('{}: {}'.format(planet_property, planet_data.get(planet_property)))

        raise ValueError("{} missing planet propert{} needed for SPAM to work. Please include this in the 'planet_data' input dictionary".format(len(missing), 'y is' if len(missing) == 1 else 'ies are'))

    # User inputs check done. Now jump into the algorithm. First, define times around transit
    times = np.linspace(-planet_data['transit_duration']/2., planet_data['transit_duration']/2., ndatapoints)

    # Now initialize models
    params, m = init_batman(times, 'nonlinear')
    params_twop, m_twop = init_batman(times, ld_law)

    # Define params according to the planet_data dictionary
    for par in [params, params_twop]:

        par.per = planet_data['orbital_period']
        par.rp = planet_data['Rp/Rs']
        par.a = planet_data['a/Rs']
        par.inc = planet_data['inclination']
        par.ecc = planet_data['eccentricity']
        par.w = planet_data['omega']

    # Set non-linear law coefficients
    params.u = [c1, c2, c3, c4]

    # Allright, now given these models, optimize u1 and u2 such that they match as good as possible the lightcurves of 
    # the non-linear law. We use scipy.optimize for this
    results = minimize(spam_objective_function, [u1_guess, u2_guess], args=(params, m, params_twop, m_twop), method=method.lower())

    return results.x, {prop: val for prop, val in planet_data.items() if prop in planet_properties}

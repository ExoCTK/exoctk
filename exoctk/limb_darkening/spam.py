import numpy as np
from scipy.optimize import minimize

try:
    import batman

except:
    print('Batman library not installed. Install it by doing "pip install batman-package" to use exoctk.limb_darkening.spam.')

from .. import utils

def init_batman(t, ld_law, nresampling = None, etresampling = None):
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

    params = batman.TransitParams()
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
        m = batman.TransitModel(params, t)

    else:
        m = batman.TransitModel(params, t, supersample_factor=nresampling, exp_time=etresampling)

    return params,m

def spam_objective_function(theta, params, m, params_twop, m_twop):
    """
    Objective function that the SPAM algorithm is looking to minimize. 
    """

    u1, u2 = theta
    params_twop.u = [u1, u2]

    return np.sum((m.light_curve(params) - m_twop.light_curve(params_twop))**2)

def transform_coefficients(c1, c2, c3, c4, planet_name = '', planet_data = None, ld_law = 'quadratic', ndatapoints = 1000, method = 'BFGS', u1_guess = 0.5, u2_guess = 0.5):
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

    if planet_name != '':

        # If planet_name is given, retrieve MAST properties:
        mast_planet_data, url = utils.get_target_data(planet_name)

        # If planet properties also given, add the ones not in that dictinoary:
        if planet_data is not None:

            for planet_property in planet_properties:

                if planet_property not in list(planet_data.keys()):

                    planet_data[planet_property] = mast_planet_data[planet_property]

        else:

            planet_data = mast_planet_data

    else:
    
        if planet_data is None:
            
            raise Exception("User must input either 'planet_name' and/or 'planet_data' for SPAM to work. See details by doing exoctk.limb_darkening.spam.transform_coefficients?.")

        # Check that all properties exist in the input dictionary:
        for planet_property in planet_properties:

            if planet_property not in list(planet_data.keys()):

                raise Exception("Input 'planet_data' does not have a '"+planet_property+"' key. This is a needed key for SPAM to work.")

    # User inputs check done. Now jump into the algorithm. First, define times around transit:
    times = np.linspace(-planet_data['transit_duration']/2., planet_data['transit_duration']/2., ndatapoints)

    # Now initialize models:
    params, m = init_batman(times, 'nonlinear')
    params_twop, m_twop = init_batman(times, ld_law)

    # Define params according to the planet_data dictionary:
    for par in [params, params_twop]:

        par.per = planet_data['orbital_period']
        par.rp = planet_data['Rp/Rs']
        par.a = planet_data['a/Rs']
        par.inc = planet_data['inclination']
        par.ecc = planet_data['eccentricity']
        par.w = planet_data['omega']

    # Set non-linear law coefficients:
    params.u = [c1, c2, c3, c4]

    # Allright, now given these models, optimize u1 and u2 such that they match as good as possible the lightcurves of 
    # the non-linear law. We use scipy.optimize for this:
    results = minimize(spam_objective_function, [u1_guess, u2_guess], args=(params, m, params_twop, m_twop)) 

    # Return SPAM coefficients:
    return results.x

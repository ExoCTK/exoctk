"""Light Curve Fitting

Generally:
1. Start with raw data of light curve
2. Load the known planet/stellar orbital and atmospheric parameters to create a theoretical light curve
3. 


Kevin
1. Create parameters file replete with systematic models
2. Load the obs data
3. Generate a model of initial conditions
4. LMFIT to determine best fit
5. MCMC routine to determine uncertainties with best fit priors


Jonathan
1. Load the obs data
2. Load planet params f(P, Tc, a/Rs, i(b), d, ldcs, e, w)
3. Plot initial lightcurve
4. Fidget with planet params if Tc can't be found
5. Use LMFIT with models to determine params


ExoCTK Version
1. Create base class instance with default parameters, also allow a parameter file argument
"""
import os
import numpy as np
import batman
import lmfit
import json
import pandas as pd


class LightCurveFitter:
    def __init__(self, model_class):
        self.flux = np.ones(100)
        self.time = np.arange(100)
        self.results = pd.DataFrame(names=('fit_number','wavelength', 'P', 'Tc', 'a/Rs', 'b', 'd', 'ldcs', 'e', 'w','model_name','chi2'))


    def run(self):
        """Run the model fits"""
        pass


    # Method to return sliced results table
    def master_slicer(self, value, param_name='wavelength'):
        return self.results.iloc[self.results[param_name]==value]


class LightCurve:
    def __init__(self, time, flux, parameters=None, fitter=None):
        """
        A class to store the actual light curve 
        
        Parameters
        ----------
        time: sequence
            The time axis in days
        flux: sequence
            The flux at each time
        parameters: str, object (optional)
            The orbital parameters of the star/planet system,
            may be a path to a JSON file or a parameter object
        fitter: 
            The instance used to perform the fit
        """
        # Check data
        if len(time)!=len(flux):
            assert ValueError('Time and flux axes must be the same length.')

        self.time = time
        self.flux = flux
        self._parameters = parameters


    @property
    def parameters(self):
        """A getter for the parameters"""
        return self._parameters


    @parameters.setter
    def parameters(self, params):
        """A setter for the parameters"""
        # Process if it is a parameters file
        if isinstance(params, str) and os.file.exists(params):
            params = json.load(params)


    def fit(self, fitter):
        self.fit_results = self.sampler.do_it()


    def plot(self):
        pass


class Parameters:
    """A class to hold the planetary parameters
    """
    def __init__(self, param_file=None, **kwargs):
        """Initialize the parameter object
        
        Parameters
        ----------
        param_file: str
            A text file of the parameters to parse
        
        """
        # Make an empty params dict
        params = {}

        # semi-major axis (in units of stellar radii)
        self._a = None

        # eccentricity
        self._ecc = None

        # inclunaition in degrees
        self._inc = None

        # limb darkening model
        self._limb_dark = None

        # orbital period
        self._per = None

        # planet radius (in units of stellar radii)
        self._rp = None

        # time of inferior conjunction
        self._t0 = None

        # limb darkening coefficients
        self._u = None

        # longitude of periastron (in degrees)
        self._w = None

        # If a param_file is given, make sure it exists
        if param_file is not None and os.path.exists(param_file):

            # Parse the ASCII file
            if param_file.endswith('.txt'):

                # Add the params to a dict
                data = np.genfromtxt(param_file)
                params = {i:j for i,j in data}

            # Parse the JSON file
            elif param_file.endswith('.json'):

                with open(param_file) as json_data:
                    params = json.load(json_data)

        # Add any kwargs to the parameter dict
        params.update(kwargs)

        # Try to store each as an attribute
        for param, value in params.items():
            setattr(self, param, value)


    @property
    def a(self):
        "A getter for semi-major axis"
        return self._a


    @a.setter
    def a(self, value):
        "A setter for semi-major axis (in units of stellar radii)"
        if not isinstance(value, float) or a<=0:
            assert ValueError("a must be a positive float")

        self._a = value


    @property
    def ecc(self):
        "A getter for eccentricity"
        return self._ecc


    @ecc.setter
    def ecc(self, value):
        "A setter for eccentricity"
        if not isinstance(value, float) or 0>value>1:
            assert ValueError("ecc must be a float between 0 and 1")

        self._ecc = value


    @property
    def inc(self):
        "A getter for inclination"
        return self._inc


    @inc.setter
    def inc(self, value):
        "A setter for inc"
        if not isinstance(value, float) or 0>value>90:
            assert ValueError("inc must be a float between 0 and 90")

        self._inc = value


    @property
    def limb_dark(self):
        "A getter for limb darkening profile"
        return self._limb_dark


    @limb_dark.setter
    def limb_dark(self, value):
        "A setter for limb darkening profile"
        profiles = ['uniform','linear','quadratic','square-root',
                    'logarithmic','exponential','3-parameter','4-parameter']

        if value not in profiles:
            assert NameError("limb_dark must be one of the following:",','.join(profiles))

        self.limb_dark = value


    @property
    def per(self):
        "A getter for period"
        return self._per


    @per.setter
    def per(self, value):
        "A setter for period"
        if not isinstance(value, float) or value<=0:
            assert ValueError("per must be a positive float")

        self._per = value


    @property
    def rp(self):
        "A getter for planet radius over stellar radius"
        return self._rp


    @rp.setter
    def rp(self, value):
        "A setter for planet radius over stellar radius"
        if not isinstance(value, float) or value<=0:
            assert ValueError("rp must be a positive float")

        self._rp = value


    @property
    def t0(self):
        "A getter for t0"
        return self._t0


    @t0.setter
    def t0(self, value):
        "A setter for t0"
        if not isinstance(value, float):
            assert ValueError("t0 must be a float")

        self.t0 = value


    @property
    def u(self):
        "A getter for limb darkening coefficients"
        return self._u


    @u.setter
    def u(self, value):
        "A setter for limb darkening coefficients"
        if not isinstance(value, (list, np.ndarray, tuple)):
            assert ValueError("u must be a sequence of limb darkening coefficients")

        self._u = value


    @property
    def w(self):
        "A getter for w"
        return self._w


    @w.setter
    def w(self, value):
        "A setter for w"
        if not isinstance(value, float) or 0>value>360:
            assert ValueError("w must be a float between 0 and 360")

        self._w = value


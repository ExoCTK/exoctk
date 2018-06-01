"""Base and child classes to handle orbital parameters

Author: Joe Filippazzo
Email: jfilippazzo@stsci.edu
"""
import os

class Parameters:
    """A class to hold the planetary parameters
    """
    def __init__(self, param_file=None, **kwargs):
        """Initialize the parameter object
        
        Parameters
        ----------
        param_file: str
            A text file of the parameters to parse
        
        Example
        -------
        params = lightcurve.Parameters(a=20, ecc=0.1, inc=89, limb_dark='quadratic')
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

        self._limb_dark = value


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
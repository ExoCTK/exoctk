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
import matplotlib.pyplot as plt

from .parameters import Parameters

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
    def __init__(self, time, flux, unc=None, parameters=None, fitter=None):
        """
        A class to store the actual light curve 
        
        Parameters
        ----------
        time: sequence
            The time axis in days
        flux: sequence
            The flux at each time
        unc: sequence
            The uncertainty on the flux
        parameters: str, object (optional)
            The orbital parameters of the star/planet system,
            may be a path to a JSON file or a parameter object
        fitter: 
            The instance used to perform the fit
        
        Example
        -------
        from ExoCTK.lightcurve_fitting import lightcurve
        time = np.arange(100)
        flux = np.random.normal([0.9 if 25<i<75 else 1 for i in range(100)], scale=0.01)
        unc = np.random.normal(size=100, scale=0.02)
        params = lightcurve.Parameters(a=20, ecc=0.1, inc=89, limb_dark='quadratic')
        lc = lightcurve.LightCurve(time, flux, unc, params)
        """
        # Check data
        if len(time)!=len(flux):
            assert ValueError('Time and flux axes must be the same length.')

        # Set the data arrays
        self.time = time
        self.flux = flux
        if unc is not None:
            if len(unc)!=len(time):
                raise ValueError('Time and unc axes must be the same length.')
                
            self.unc = unc
            
        else:
            self.unc = np.array([np.nan]*len(self.time))
        
        # Store the orbital parameters
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
            params = Parameters(params)
            
        # Or a Parameters instance
        if not isinstance(params, Parameters):
            raise TypeError("'params' argument must be a JSON file, ascii file, or ExoCTK.lightcurve_fitting.lightcurve.Parameters instance.")
            
        # Set the parameters attribute
        self._parameters = params


    def fit(self, model, fitter):
        """Fit the model to the lightcurve
        
        Parameters
        ----------
        model: ExoCTK.lightcurve_fitter.models.Model
            The model to fit to the data
        fitter: 
        """
        # Interpolate model to data
        mod = model.model
        
        
        self.fit_results = self.sampler.do_it()


    def plot(self):
        """Plot the light curve with all available fits"""
        plt.figure()
        
        plt.errorbar(self.time, self.flux, yerr=self.unc, marker='o', ls='none')
        
        plt.xlabel('Days')
        plt.ylabel('Rp/R*')


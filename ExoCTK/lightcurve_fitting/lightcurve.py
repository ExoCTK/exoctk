"""Base and child classes to handle light curve fitting

Author: Joe Filippazzo
Email: jfilippazzo@stsci.edu
"""
import os
import numpy as np
import batman
import lmfit
import json
import pandas as pd
import matplotlib.pyplot as plt

from .parameters import Parameters
from .models import Model

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


class LightCurve(Model):
    def __init__(self, time, flux, unc=None, parameters=None, fitter=None, units='day'):
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
        # Initialize the model
        super().__init__()
        
        # Check data
        if len(time)!=len(flux):
            raise ValueError('Time and flux axes must be the same length.')

        # Set the data arrays
        if unc is not None:
            if len(unc)!=len(time):
                raise ValueError('Time and unc axes must be the same length.')
                
            self.unc = unc
            
        else:
            self.unc = np.array([np.nan]*len(self.time))
        
        # Set the time and flux axes
        self.time = time
        self.flux = flux
        
        # Set the units
        self.units = units
        
        # Store the orbital parameters
        if parameters is not None:
            self.parameters = parameters
        
        
    def fit(self, model, fitter='lmfit', **kwargs):
        """Fit the model to the lightcurve
        
        Parameters
        ----------
        model: ExoCTK.lightcurve_fitter.models.Model
            The model to fit to the data
        fitter: str
            The name of the fitter to use
        """
        # Interpolate model to data
        model.interp(self.time, self.units)
        
        if fitter=='lmfit':
            
            return 'lmfit'


    def plot(self):
        """Plot the light curve with all available fits"""
        plt.figure()
        
        plt.errorbar(self.time, self.flux, yerr=self.unc, marker='o', ls='none')
        
        plt.xlabel(self.units.long_names[0])
        plt.ylabel('Flux')


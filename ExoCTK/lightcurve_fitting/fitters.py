"""Functions used to fit models to light curve data

Author: Joe Filippazzo
Email: jfilippazzo@stsci.edu
"""
import numpy as np
import matplotlib.pyplot as plt
import batman
import lmfit

from .models import Model
    
def lmfitter(data, model, unc=None, method='leastsq', verbose=True):
    """Use lmfit
    
    Parameters
    ----------
    data: sequence
        The observational data
    model: ExoCTK.lightcurve_fitting.models.Model
        The model to fit
    unc: np.ndarray (optional)
        The uncertainty on the (same shape) data
    
    Returns
    -------
    lmfit.Model.fit.fit_report
        The results of the fit
    """
    # Initialize lmfit Params object
    initialParams = lmfit.Parameters()
    
    param_list = []
    indep_vars = {}
    for param in model.parameters.list:
        param = list(param)
        if param[2] == 'free':
            param[2] = True
            param_list.append(tuple(param))
        elif param[2] == 'fixed':
            param[2] = False
            param_list.append(tuple(param))
        else:
            indep_vars[param[0]] = param[1]

    # Get values from input ExoCTK.lightcurve_fitting.parameters.Parameters instance
    # Format: (key, value, vary?, min, max)
    initialParams.add_many(*param_list)
        
    # Create the lightcurve model
    # lcmodel = lmfit.Model(model, independent_vars=['times', 'ldtype', 'transitType'])
    lcmodel = lmfit.Model(model.eval)
    lcmodel.independent_vars = indep_vars.keys()
    
    # Set the unc
    if unc is None:
        unc = np.ones(len(data))
        
    # Fit light curve model to the simulated data
    result = lcmodel.fit(data, weights=1/unc, params=initialParams, method=method, **indep_vars)
    if verbose:
        print(result.fit_report())
    
    # Create new model with best fit parameters
    best_fit = Model(time=model.time, flux=result.best_fit, name='Best Fit')
    
    return best_fit
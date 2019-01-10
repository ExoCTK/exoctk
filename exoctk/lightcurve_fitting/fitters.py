"""Functions used to fit models to light curve data

Author: Joe Filippazzo
Email: jfilippazzo@stsci.edu
"""
import numpy as np
import lmfit
import copy

from .parameters import Parameters


def lmfitter(time, data, model, unc=None, method='powell', verbose=True, **kwargs):
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

    # Concatenate the lists of parameters
    all_params = [i for j in [model.components[n].parameters.list
                  for n in range(len(model.components))] for i in j]

    # Group the different variable types
    param_list = []
    indep_vars = {}
    for param in all_params:
        param = list(param)
        if param[2] == 'free':
            param[2] = True
            param_list.append(tuple(param))
        elif param[2] == 'fixed':
            param[2] = False
            param_list.append(tuple(param))
        else:
            indep_vars[param[0]] = param[1]

    # Add the time as an independent variable
    indep_vars['time'] = time

    # Get values from input parameters.Parameters instances
    initialParams.add_many(*param_list)

    # Create the lightcurve model
    lcmodel = lmfit.Model(model.eval)
    lcmodel.independent_vars = indep_vars.keys()

    # Set the unc
    if unc is None:
        unc = np.ones(len(data))

    # Fit light curve model to the simulated data
    result = lcmodel.fit(data, weights=1/unc, params=initialParams,
                         method=method, **indep_vars, **kwargs)
    if verbose:
        print(result.fit_report())

    # Get the best fit params
    fit_params = result.__dict__['params']
    new_params = [(fit_params.get(i).name, fit_params.get(i).value,
                   fit_params.get(i).vary, fit_params.get(i).min,
                   fit_params.get(i).max) for i in fit_params]

    # Create new model with best fit parameters
    params = Parameters()

    # Try to store each as an attribute
    for param in new_params:
        setattr(params, param[0], param[1:])

    # Make a new model instance
    best_model = copy.copy(model)
    best_model.name = 'Best Fit'
    best_model.parameters = params

    return best_model

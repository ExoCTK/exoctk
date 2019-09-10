#! /usr/bin/env python

"""This module provides a few examples of how to use the
``platon_wrapper`` module in ExoCTK's ``atmopsheric_retrieval``
subpackage.

Authors
-------

    - Matthew Bourque

Use
---

    To run all examples, one can execute this script via the command
    line as such:

        >>> python platon_example.py

    To run individual examples, one can import the ``example()``
    function and pass as a parameter which example to run.  Available
    examples include:

        >>> from platon_example import example
        example('emcee')
        example('multinest')

Dependencies
------------

    - numpy
    - platon
"""

import numpy as np
from platon.constants import R_sun, R_jup, M_jup

from platon_wrapper import PlatonWrapper


def example(method):
    """Performs an example run of the emcee and multinest retrievals

    Parameters
    ----------
    method : str
        The method to use to perform the atmopsheric retrieval; can
        either be ``multinest`` or ``emcee``
    """

    # Ensure that the method parameter is valid
    assert method in ['multinest', 'emcee'], \
        'Unrecognized method: {}'.format(method)

    # Define the fit parameters
    params = {
        'Rs': 1.19,  # Required
        'Mp': 0.73,  # Required
        'Rp': 1.4,  # Required
        'T': 1200.0,  # Required
        'logZ': 0,  # Optional
        'CO_ratio': 0.53,  # Optional
        'log_cloudtop_P': 4,  # Optional
        'log_scatt_factor': 0,  # Optional
        'scatt_slope': 4,  # Optional
        'error_multiple': 1,  # Optional
        'T_star': 6091}  # Optional

    # Initialize the object and set the parameters
    pw = PlatonWrapper()
    pw.set_parameters(params)

    # Fit for the stellar radius and planetary mass using Gaussian priors.  This
    # is a way to account for the uncertainties in the published values
    pw.fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
    pw.fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)

    # Fit for other parameters using uniform priors
    R_guess = 1.4 * R_jup
    T_guess = 1200
    pw.fit_info.add_uniform_fit_param('Rp', 0.9*R_guess, 1.1*R_guess)
    pw.fit_info.add_uniform_fit_param('T', 0.5*T_guess, 1.5*T_guess)
    pw.fit_info.add_uniform_fit_param("log_scatt_factor", 0, 1)
    pw.fit_info.add_uniform_fit_param("logZ", -1, 3)
    pw.fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 5)
    pw.fit_info.add_uniform_fit_param("error_multiple", 0.5, 5)

    # Define bins, depths, and errors
    pw.wavelengths = 1e-6*np.array([1.119, 1.1387])
    pw.bins = [[w-0.0095e-6, w+0.0095e-6] for w in pw.wavelengths]
    pw.depths = 1e-6 * np.array([14512.7, 14546.5])
    pw.errors = 1e-6 * np.array([50.6, 35.5])

    # Do some retrievals
    pw.retrieve(method)
    if method == 'multinest':
        pw.save_results()

    # Make corner plot of results
    pw.make_plot()

    return pw


if __name__ == '__main__':

    example('emcee')
    example('multinest')

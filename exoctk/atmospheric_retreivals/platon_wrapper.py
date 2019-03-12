"""A wrapper around the PLATON software.

This module serves as a wrapper around the atmospheric retreival
software for ``platon``.  For more information about ``platon``, please
see ``https://platon.readthedocs.io``.

Authors
-------
    - Matthew Bourque

Use
---

    Coming soon.

Dependencies
------------

    - ``corner``
    - ``numpy``
    - ``platon``
"""

import argparse
import os
import yaml

import corner
import numpy as np
from platon.fit_info import FitInfo
from platon.retriever import Retriever
from platon.constants import R_sun, R_jup, M_jup


def _apply_factors(params):
    """Apply appropriate multiplication factors to parameters.

    Parameters
    ----------
    params : dict
        A dictionary of parameters and their values for running the
        software.  See "Use" documentation for further details.
    """

    params['Rs'] = params['Rs'] * R_sun
    params['Mp'] = params['Mp'] * M_jup
    params['Rp'] = params['Rp'] * R_jup

    return params


def _parse_args():
    """Parse command line arguments.

    Returns
    -------
    args : obj
        An ``argparse`` object containing all of the added arguments.
    """

    param_file_help = 'The path to a file containing the parameters '
    param_file_help += 'needed for PLATON'

    parser = argparse.ArgumentParser()
    parser.add_argument('parameter_file', type=str, help=param_file_help)

    args = parser.parse_args()

    return args


def _parse_parameter_file(parameter_file):
    """Parses the supplied parameter file, ensures that the required
    parameters exist, and ensures that all supplied parameters are of a
    valid data type.  Also applies appropriate multiplication factors
    to the appropriate parameters.

    Parameters
    ----------
    parameter_file : str
        The path to the parameter file.

    Returns
    -------
    params : dict
        A dictionary containing parameter name/value pairs.
    """

    # Parse the parameter file
    with open(parameter_file, 'r') as f:
        params = yaml.load(f)

    return params


def _test_args(args):
    """Ensures that the command line arguments are of proper format and
    valid. If they are not, an assertion error is raised.

    Parameters
    ----------
    args : obj
        The ``argparse`` object containing the command line arguments.
    """

    assert os.path.exists(args.parameter_file), 'Parameter file does not exist.'


def _validate_parameters(params):
    """Ensure the supplied parameters are valid.  Throw assertion
    errors if they are not.

    Parameters
    ----------
    params : dict
        A dictionary of parameters and their values for running the
        software.  See "Use" documentation for further details.
    """

   # Make sure the parameter file contains the required keywords and they have values
    required_parameters = ['Rs', 'Mp', 'Rp', 'T']
    param_types = [float, float, float, int]
    for param, param_type in zip(required_parameters, param_types):
        assert param in params, '{} missing from parameter file'.format(param)
        assert type(params[param]) == param_type, '{} is not of type {}'.format(param, param_type)

    # Make sure the optional keywords are of proper type
    optional_parameters = ['logZ', 'CO_ratio', 'log_cloudtop_P', 'log_scatt_factor',
                           'scatt_slope', 'error_multiple', 'T_star']
    param_types = [int, float, int, int, int, int, int]
    for param, param_type in zip(optional_parameters, param_types):
        if param in params:
            assert type(params[param]) == param_type, '{} is not of type {}'.format(param, param_type)


class PlatonWrapper():
    """Class object for running the platon atmospheric retrieval
    software."""

    def __init__(self):
        """Initialize the class object."""

        self.retriever = Retriever()


    def retrieve(self):
        """Perform the atmopsheric retreival."""

        # Run nested sampling
        result = self.retriever.run_multinest(self.bins, self.depths, self.errors, self.fit_info, plot_best=True)
        print(result)

        # Do some plotting
        fig = corner.corner(result.samples, weights=result.weights,
                            range=[0.99] * result.samples.shape[1],
                            labels=self.fit_info.fit_param_names)
        fig.savefig("emcee_corner.png")

    def set_parameters(self, params):
        """Set necessary parameters to perform the retrieval.

        Required parameters include ``Rs``, ``Mp``, ``Rp``, and ``T``.
        Optional parameters include ``logZ``, ``CO_ratio``,
        ``log_cloudtop_P``, ``log_scatt_factor``, ``scatt_slope``,
        ``error_multiple``, and ``T_star``.

        Parameters
        ----------
        params : dict
            A dictionary of parameters and their values for running the
            software.  See "Use" documentation for further details.
        """

        _validate_parameters(params)
        _apply_factors(params)
        self.params = params
        self.fit_info = self.retriever.get_default_fit_info(**self.params)


if __name__ == '__main__':

    # Parse command line arguments
    args = _parse_args()

    # Test command line arguments
    _test_args(args)

    # Parse and test the parameter file
    params = _parse_parameter_file(args.parameter_file)

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
    # pw.fit_info.add_uniform_fit_param("log_scatt_factor", 0, 1)
    # pw.fit_info.add_uniform_fit_param("logZ", -1, 3)
    # pw.fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 5)
    # pw.fit_info.add_uniform_fit_param("error_multiple", 0.5, 5)

    # Define bins, depths, and errors
    pw.wavelengths = 1e-6*np.array([1.119, 1.138, 1.157, 1.175, 1.194, 1.213, 1.232, 1.251, 1.270, 1.288, 1.307, 1.326, 1.345, 1.364, 1.383, 1.401, 1.420, 1.439, 1.458, 1.477, 1.496, 1.515, 1.533, 1.552, 1.571, 1.590, 1.609, 1.628])
    pw.bins = [[w-0.0095e-6, w+0.0095e-6] for w in pw.wavelengths]
    pw.depths = 1e-6 * np.array([14512.7, 14546.5, 14566.3, 14523.1, 14528.7, 14549.9, 14571.8, 14538.6, 14522.2, 14538.4, 14535.9, 14604.5, 14685.0, 14779.0, 14752.1, 14788.8, 14705.2, 14701.7, 14677.7, 14695.1, 14722.3, 14641.4, 14676.8, 14666.2, 14642.5, 14594.1, 14530.1, 14642.1])
    pw.errors = 1e-6 * np.array([50.6, 35.5, 35.2, 34.6, 34.1, 33.7, 33.5, 33.6, 33.8, 33.7, 33.4, 33.4, 33.5, 33.9, 34.4, 34.5, 34.7, 35.0, 35.4, 35.9, 36.4, 36.6, 37.1, 37.8, 38.6, 39.2, 39.9, 40.8])

    pw.retrieve()

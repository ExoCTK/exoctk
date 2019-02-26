"""
"""

import argparse
import os
import yaml

import corner
import numpy as np
from platon.fit_info import FitInfo
from platon.retriever import Retriever
from platon.constants import R_sun, R_jup, M_jup


def parse_args():
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


def parse_parameter_file(parameter_file):
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
    param_dict : dict
        A dictionary containing parameter name/value pairs.
    """

    # Parse the parameter file
    with open(parameter_file, 'r') as f:
        param_dict = yaml.load(f)

    # Make sure the parameter file contains the required keywords and they have values
    required_parameters = ['Rs', 'Mp', 'Rp', 'T']
    param_types = [float, float, float, int]
    for param, param_type in zip(required_parameters, param_types):
        assert param in param_dict, '{} missing from parameter file'.format(param)
        assert type(param_dict[param]) == param_type, '{} is not of type {}'.format(param, param_type)

    # Make sure the optional keywords are of proper type
    optional_parameters = ['logZ', 'CO_ratio', 'log_cloudtop_P', 'log_scatt_factor',
                           'scatt_slope', 'error_multiple', 'T_star']
    param_types = [int, float, int, int, int, int, int]
    for param, param_type in zip(optional_parameters, param_types):
        if param in param_dict:
            assert type(param_dict[param]) == param_type, '{} is not of type {}'.format(param, param_type)

    # Apply proper multiplication factors
    param_dict['Rs'] = param_dict['Rs'] * R_sun
    param_dict['Mp'] = param_dict['Mp'] * M_jup
    param_dict['Rp'] = param_dict['Rp'] * R_jup

    return param_dict


def test_args(args):
    """Ensures that the command line arguments are of proper format and
    valid. If they are not, an assertion error is raised.

    Parameters
    ----------
    args : obj
        The ``argparse`` object containing the command line arguments.
    """

    assert os.path.exists(args.parameter_file), 'Parameter file does not exist.'


class PlatonWrapper():
    """
    """

    def __init__(self, kwargs):
        """
        """

        self.kwargs = kwargs


    def retrieve(self):
        """
        """

        # Construct fit_info object
        retriever = Retriever()
        fit_info = retriever.get_default_fit_info(**self.kwargs)

        # Fit for the stellar radius and planetary mass using Gaussian priors.  This
        # is a way to account for the uncertainties in the published values
        fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
        fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)

        # Fit for other parameters using uniform priors
        R_guess = 1.4 * R_jup
        T_guess = 1200
        fit_info.add_uniform_fit_param('Rp', 0.9*R_guess, 1.1*R_guess)
        fit_info.add_uniform_fit_param('T', 0.5*T_guess, 1.5*T_guess)
        fit_info.add_uniform_fit_param("log_scatt_factor", 0, 1)
        fit_info.add_uniform_fit_param("logZ", -1, 3)
        fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 5)
        fit_info.add_uniform_fit_param("error_multiple", 0.5, 5)

        # Define bins, depths, and errors
        wavelengths = 1e-6*np.array([1.119, 1.138, 1.157, 1.175, 1.194, 1.213, 1.232, 1.251, 1.270, 1.288, 1.307, 1.326, 1.345, 1.364, 1.383, 1.401, 1.420, 1.439, 1.458, 1.477, 1.496, 1.515, 1.533, 1.552, 1.571, 1.590, 1.609, 1.628])
        bins = [[w-0.0095e-6, w+0.0095e-6] for w in wavelengths]
        depths = 1e-6 * np.array([14512.7, 14546.5, 14566.3, 14523.1, 14528.7, 14549.9, 14571.8, 14538.6, 14522.2, 14538.4, 14535.9, 14604.5, 14685.0, 14779.0, 14752.1, 14788.8, 14705.2, 14701.7, 14677.7, 14695.1, 14722.3, 14641.4, 14676.8, 14666.2, 14642.5, 14594.1, 14530.1, 14642.1])
        errors = 1e-6 * np.array([50.6, 35.5, 35.2, 34.6, 34.1, 33.7, 33.5, 33.6, 33.8, 33.7, 33.4, 33.4, 33.5, 33.9, 34.4, 34.5, 34.7, 35.0, 35.4, 35.9, 36.4, 36.6, 37.1, 37.8, 38.6, 39.2, 39.9, 40.8])

        # Run nested sampling
        result = retriever.run_multinest(bins, depths, errors, fit_info, plot_best=True)
        print(result)

        # Do some plotting
        fig = corner.corner(result.samples, weights=result.weights,
                            range=[0.99] * result.samples.shape[1],
                            labels=fit_info.fit_param_names)
        fig.savefig("emcee_corner.png")


if __name__ == '__main__':

    # Parse command line arguments
    args = parse_args()

    # Test command line arguments
    test_args(args)

    # Parse and test the parameter file
    kwargs = parse_parameter_file(args.parameter_file)

    platon_wrapper = PlatonWrapper(kwargs)
    platon_wrapper.retrieve()

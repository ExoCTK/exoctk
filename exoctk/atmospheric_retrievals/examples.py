#! /usr/bin/env python

"""This module provides examples of how to use the ``platon_wrapper``
module in ExoCTK's ``atmopsheric_retrieval`` subpackage.  Each example
can be run using either the ``multinest`` method or the ``emcee``
method.  See the docstring of each example function for further
details.

Authors
-------

    - Matthew Bourque

Use
---

    To run all examples, one can execute this script via the command
    line as such:

        >>> python platon_example_aws.py

    To run individual examples, one can import the example functions
    and pass as a parameter which method to run.  Available
    examples include:

        from examples import example, example_aws, example_aws_long
        example('emcee')
        example('multinest')
        example_aws('emcee')
        example_aws('multinest')
        example_aws_long('emcee')
        example_aws_long('multinest')

Dependencies
------------

    Dependent libraries include:

    - numpy
    - platon

    Users must also have a ``aws_config.json`` file present within the
    ``atmospheric_retrievals`` subdirectory.  This file must be of
    a valid JSON format and contain two key/value pairs,
    ``ec2_id`` and ``ssh_file``, e.g.:

    {
    "ec2_id" : "lt-021de8b904bc2b728",
    "ssh_file" : "~/.ssh/my_ssh_key.pem"
    }

    where the ``ec2_id`` contains the ID for an EC2 launch template
    or an existing EC2 instance, and ``ssh_file`` points to the SSH
    public key used for logging into an AWS account.

    Note that if the ``ec2_id`` points to a launch template (i.e. the
    string starts with ``lt-``), a new EC2 instance will be created and
    launched.  However, if the ``ec2_id`` points to an existing EC2
    instance (i.e. the string starts with ``i-``), the existing EC2
    instance will be started and used.
"""

import logging
import time

import numpy as np
import pandas
from platon.constants import R_sun, R_jup, M_jup

from exoctk.atmospheric_retrievals.aws_tools import get_config
from exoctk.atmospheric_retrievals.platon_wrapper import PlatonWrapper


def example(method):
    """Performs a short example run of the retrievals.

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


def example_aws(method):
    """Performs an example run of the ``emcee`` or ``multinest``
    retrievals using AWS.

    Parameters
    ----------
    method : str
        The method to use to perform the atmopsheric retrieval; can
        either be ``multinest`` or ``emcee``
    """

    ssh_file = get_config()['ssh_file']
    ec2_id = get_config()['ec2_id']

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

    # Initialize the object, set parameters, and perform retreival
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

    # Set use for AWS and perform retreival
    pw.use_aws(ssh_file, ec2_id)
    pw.retrieve(method)


def example_aws_long(method):
    """Performs an longer example run of the ``emcee`` and
    `multinest`` retrievals using AWS.

    Parameters
    ----------
    method : str
        The method to use to perform the atmopsheric retrieval; can
        either be ``multinest`` or ``emcee``
    """

    ssh_file = get_config()['ssh_file']
    ec2_id = get_config()['ec2_id']

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

    # Initialize the object, set parameters, and perform retreival
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

    # Get bins, depths, and errors
    bins, depths, errors = get_example_data('wasp-19b')

    pw.bins = bins
    pw.depths = depths
    pw.errors = errors

    # Set use for AWS and perform retreival
    pw.use_aws(ssh_file, ec2_id)
    pw.retrieve(method)


def get_example_data(object_name):
    """Return ``bins``, ``depths``, and ``errors`` for the given
    ``object_name``.  Data is read in from a ``csv`` file with a
    corresponding filename.

    Parameters
    ----------
    object_name : str
        The object of interest (e.g. ``hd209458b``)

    Returns
    -------
    bins : np.array
        A 2xN ``numpy`` array of wavelength bins, of the form
        ``[[wavelength_bin_min, wavelength_bin_max], ...]``
    depths : np.array
        A 1D ``numpy`` array of depth values
    errors: np.array
        A 1D ``numpy`` array of depth error values.
    """

    logging.info('Using data for {}'.format(object_name))

    # Read in the data
    df = pandas.read_csv('data/{}.csv'.format(object_name), names=['wavelengths', 'bin_sizes', 'depths', 'errors'])

    # Remove and rows outside of wavelength range (3e-7 to 3e-5)
    df = df.loc[(1e-6*df['wavelengths'] - 1e-6*df['bin_sizes'] >= 3e-7) & (1e-6*df['wavelengths'] + 1e-6*df['bin_sizes'] <= 3e-5)]

    # Parse the data
    wavelengths = 1e-6*np.array(df['wavelengths'])
    bin_sizes = 1e-6*np.array(df['bin_sizes'])
    depths = np.array(df['depths'])
    errors = np.array(df['errors'])

    # Calculate bins
    bins = [[wavelength - bin_size, wavelength + bin_size] for wavelength, bin_size in zip(wavelengths, bin_sizes)]

    return bins, depths, errors


if __name__ == '__main__':

    example_aws_long('multinest')
    time.sleep(120)  # Wait a few minutes for the existing EC2 instance to completely stop
    example_aws_long('emcee')

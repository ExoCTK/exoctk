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

        >>> python examples.py

    To run individual examples, one can import the example functions
    and pass as a parameter which method to run.  Available
    examples include:

        from examples import example, example_aws_short, example_aws_long
        example('emcee')
        example('multinest')
        example_aws_short('emcee')
        example_aws_short('multinest')
        example_aws_long('emcee')
        example_aws_long('multinest')

Dependencies
------------

    Dependent libraries include:

    - ``exoctk``
    - ``numpy``
    - ``pandas``
    - ``platon``

    To run the examples that use AWS, users must also have a
    ``aws_config.json`` file present within the
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

References
----------

    Example data was pulled from the "Transmission Spectra" tab from
    corresponding ExoMAST pages, available at
    https://exo.mast.stsci.edu/
"""

import logging
import os
import time

import numpy as np
import pandas
from platon.constants import R_sun, R_jup, M_jup

from exoctk.atmospheric_retrievals.aws_tools import get_config
from exoctk.atmospheric_retrievals.platon_wrapper import PlatonWrapper


def example(method):
    """Performs a short example run of the retrievals using local
    machine.

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
    pw.fit_info.add_uniform_fit_param('Rp', 0.9*(1.4 * R_jup), 1.1*(1.4 * R_jup))
    pw.fit_info.add_uniform_fit_param('T', 0.5*1200, 1.5*1200)
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

    # Save the results
    pw.save_results()

    # Make corner plot of results
    pw.make_plot()

    return pw


def example_aws_short(method):
    """Performs an short example run of the retrievals using AWS.

    Parameters
    ----------
    method : str
        The method to use to perform the atmopsheric retrieval; can
        either be ``multinest`` or ``emcee``
    """

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

    # Initialize the object, set parameters, and perform retrieval
    pw = PlatonWrapper()
    pw.set_parameters(params)

    # Fit for the stellar radius and planetary mass using Gaussian priors.  This
    # is a way to account for the uncertainties in the published values
    pw.fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
    pw.fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)

    # Fit for other parameters using uniform priors
    pw.fit_info.add_uniform_fit_param('Rp', 0.9*(1.4 * R_jup), 1.1*(1.4 * R_jup))
    pw.fit_info.add_uniform_fit_param('T', 0.5*1200, 1.5*1200)
    pw.fit_info.add_uniform_fit_param("log_scatt_factor", 0, 1)
    pw.fit_info.add_uniform_fit_param("logZ", -1, 3)
    pw.fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 5)
    pw.fit_info.add_uniform_fit_param("error_multiple", 0.5, 5)

    # Define bins, depths, and errors
    pw.wavelengths = 1e-6*np.array([1.119, 1.1387])
    pw.bins = [[w-0.0095e-6, w+0.0095e-6] for w in pw.wavelengths]
    pw.depths = 1e-6 * np.array([14512.7, 14546.5])
    pw.errors = 1e-6 * np.array([50.6, 35.5])

    # Set use for AWS and perform retrieval
    ssh_file = get_config()['ssh_file']
    ec2_id = get_config()['ec2_id']
    pw.use_aws(ssh_file, ec2_id)
    pw.retrieve(method)


def example_aws_long(method):
    """Performs an longer example run of the retrievals using AWS.

    Parameters
    ----------
    method : str
        The method to use to perform the atmopsheric retrieval; can
        either be ``multinest`` or ``emcee``
    """

    params = {
        'Rs': 1.19,
        'Mp': 0.73,
        'Rp': 1.39,
        'T': 1476.81,
        'logZ': 0,
        'CO_ratio': 0.53,
        'log_cloudtop_P': 4,
        'log_scatt_factor': 0,
        'scatt_slope': 4,
        'error_multiple': 1,
        'log_cloudtop_P': 4}

    # Initialize the object, set parameters, and perform retrieval
    pw = PlatonWrapper()
    pw.set_parameters(params)

    if method == 'multinest':
        pw.fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
        pw.fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)
        pw.fit_info.add_uniform_fit_param('Rp', 0.9*(1.39 * R_jup), 1.1*(1.39 * R_jup))
        pw.fit_info.add_uniform_fit_param('T', 300, 3000)
        pw.fit_info.add_uniform_fit_param("log_scatt_factor", 0, 2)
        pw.fit_info.add_uniform_fit_param("logZ", -1, 3)
        pw.fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 7)
        pw.fit_info.add_uniform_fit_param("error_multiple", 0.5, 5)
    elif method == 'emcee':
        pw.fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
        pw.fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)
        pw.fit_info.add_uniform_fit_param('Rp', 0, np.inf, 0.9*(1.39 * R_jup), 1.1*(1.39 * R_jup))
        pw.fit_info.add_uniform_fit_param('T', 300, 3000, 0.5*1476.81, 1.5*1476.81)
        pw.fit_info.add_uniform_fit_param("log_scatt_factor", 0, 5, 0, 2)
        pw.fit_info.add_uniform_fit_param("logZ", -1, 3)
        pw.fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 7)
        pw.fit_info.add_uniform_fit_param("error_multiple", 0, np.inf, 0.5, 5)

    # Get bins, depths, and errors
    bins, depths, errors = get_example_data('hd209458b')
    pw.bins = bins
    pw.depths = depths
    pw.errors = errors

    # Set use for AWS and perform retrieval
    ssh_file = get_config()['ssh_file']
    ec2_id = get_config()['ec2_id']
    pw.use_aws(ssh_file, ec2_id)
    pw.retrieve(method)


def get_example_data(object_name):
    """Return ``bins``, ``depths``, and ``errors`` for the given
    ``object_name``.  Data is read in from a ``csv`` file with a
    filename corresponding to ``object_name``.

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
    data_file = os.path.join(os.path.dirname(__file__), 'example_data', '{}.csv'.format(object_name))
    df = pandas.read_csv(data_file, names=['wavelengths', 'bin_sizes', 'depths', 'errors'])

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

    # A short example using local machine
    example('multinest')
    example('emcee')

    # A short example using AWS
    example_aws_short('multinest')
    time.sleep(120)  # Allow time for the EC2 instance to restart
    example_aws_short('emcee')

    # A long example using AWS
    time.sleep(120)  # Allow time for the EC2 instance to restart
    example_aws_long('multinest')
    time.sleep(120)  # Allow time for the EC2 instance to restart
    example_aws_long('emcee')

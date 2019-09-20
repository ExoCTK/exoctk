#! /usr/bin/env python

"""This module provides an example of how to use the
``platon_wrapper`` module in ExoCTK's ``atmopsheric_retrieval``
subpackage, with focus on utilizing the feature for performing the
processing on AWS.

Authors
-------

    - Matthew Bourque

Use
---

    To run the example, one can execute this script via the command
    line as such:

        >>> python platon_example_aws.py

    To run individual examples, one can import the ``example_aws()``
    function and pass as a parameter which example to run.  Available
    examples include:

        from platon_example_aws import example_aws, example_aws_long
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

import time

import numpy as np
from platon.constants import R_sun, R_jup, M_jup

from exoctk.atmospheric_retrievals.aws_tools import get_config
from exoctk.atmospheric_retrievals.platon_wrapper import PlatonWrapper


def example_aws(method):
    """Performs an example run of the ``emcee`` and ``multinest``
    retrievals using AWS

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
    `multinest`` retrievals using AWS

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
    stis_bins, stis_depths, stis_errors = hd209458b_stis()
    wfc3_bins, wfc3_depths, wfc3_errors = hd209458b_wfc3()
    spitzer_bins, spitzer_depths, spitzer_errors = hd209458b_spitzer()

    pw.wavelengths = 1e-6*np.array([1.119, 1.138, 1.157, 1.175, 1.194, 1.213, 1.232, 1.251, 1.270, 1.288, 1.307, 1.326, 1.345, 1.364, 1.383, 1.401, 1.420, 1.439, 1.458, 1.477, 1.496, 1.515, 1.533, 1.552, 1.571, 1.590, 1.609, 1.628])
    pw.bins = np.concatenate([stis_bins, wfc3_bins, spitzer_bins])
    pw.depths = np.concatenate([stis_depths, wfc3_depths, spitzer_depths])
    pw.errors = np.concatenate([stis_errors, wfc3_errors, spitzer_errors])

    # Set use for AWS and perform retreival
    pw.use_aws(ssh_file, ec2_id)
    pw.retrieve(method)


def hd209458b_stis():
    """Returns a ``hd209458b`` example for STIS.  Taken from ``platon`` examples"""
    #http://iopscience.iop.org/article/10.1086/510111/pdf
    star_radius = 1.125 * R_sun
    wave_bins = [[300.1,347], [348,402], [403,457], [458,512], [512,567], [532,629], [629,726], [727,824], [825,922], [922,1019]]
    wave_bins = 1e-9 * np.array(wave_bins)

    planet_radii = [1.3263, 1.3254, 1.32, 1.3179, 1.3177, 1.3246, 1.3176, 1.3158, 1.32, 1.3268]
    radii_errors = [0.0018, 0.0010, 0.0006, 0.0006, 0.0010, 0.0006, 0.0005, 0.0006, 0.0006, 0.0013]
    transit_depths = (np.array(planet_radii)*R_jup/star_radius)**2 + 60e-6
    transit_errors = np.array(radii_errors)/np.array(planet_radii) * 2 * transit_depths
    return wave_bins, transit_depths, transit_errors

def hd209458b_wfc3():
    """Returns a ``hd209458b`` example for WFC3.  Taken from ``platon`` examples"""
    #https://arxiv.org/pdf/1302.1141.pdf
    wavelengths = 1e-6*np.array([1.119, 1.138, 1.157, 1.175, 1.194, 1.213, 1.232, 1.251, 1.270, 1.288, 1.307, 1.326, 1.345, 1.364, 1.383, 1.401, 1.420, 1.439, 1.458, 1.477, 1.496, 1.515, 1.533, 1.552, 1.571, 1.590, 1.609, 1.628])
    wavelength_bins = [[w-0.0095e-6, w+0.0095e-6] for w in wavelengths]
    depths = 1e-6 * np.array([14512.7, 14546.5, 14566.3, 14523.1, 14528.7, 14549.9, 14571.8, 14538.6, 14522.2, 14538.4, 14535.9, 14604.5, 14685.0, 14779.0, 14752.1, 14788.8, 14705.2, 14701.7, 14677.7, 14695.1, 14722.3, 14641.4, 14676.8, 14666.2, 14642.5, 14594.1, 14530.1, 14642.1])
    errors = 1e-6 * np.array([50.6, 35.5, 35.2, 34.6, 34.1, 33.7, 33.5, 33.6, 33.8, 33.7, 33.4, 33.4, 33.5, 33.9, 34.4, 34.5, 34.7, 35.0, 35.4, 35.9, 36.4, 36.6, 37.1, 37.8, 38.6, 39.2, 39.9, 40.8])
    return np.array(wavelength_bins), depths, errors

def hd209458b_spitzer():
    """Returns a ``hd209458b`` example for Spitzer.  Taken from ``platon`` examples"""
    #https://arxiv.org/pdf/1504.05942.pdf

    wave_bins = []
    depths = []
    errors = []

    wave_bins.append([3.2, 4.0])
    RpRs = np.average([0.12077, 0.1222, 0.11354, 0.11919], weights=1.0/np.array([0.00085, 0.00062, 0.00087, 0.00032]))
    depths.append(RpRs**2)
    errors.append(0.00032/RpRs * 2 * depths[-1])

    wave_bins.append([4.0, 5.0])
    RpRs = np.average([0.12199, 0.12099], weights=1.0/np.array([0.00094, 0.00029]))
    depths.append(RpRs**2)
    errors.append(0.00029/RpRs * 2 * depths[-1])

    wave_bins.append([5.1, 6.3])
    RpRs = np.average([0.12007, 0.11880], weights=1.0/np.array([0.00248, 0.00272]))
    depths.append(RpRs**2)
    errors.append(0.00248/RpRs * 2 * depths[-1])

    wave_bins.append([6.6, 9.0])
    RpRs = np.average([0.12007, 0.11991], weights=1.0/np.array([0.00114, 0.00073]))
    depths.append(RpRs**2)
    errors.append(0.00073/RpRs * 2 * depths[-1])

    return 1e-6*np.array(wave_bins), np.array(depths), np.array(errors)


if __name__ == '__main__':

    example_aws_long('multinest')
    time.sleep(120)  # Wait a few minutes for the existing EC2 instance to completely stop
    example_aws_long('emcee')

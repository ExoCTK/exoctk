#! /usr/bin/env python

"""Tests for the ``atmopshric_retrievals`` package.

Authors
-------

    Matthew Bourque

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_atmospheric_retrievals.py
"""

import glob
import shutil

import numpy as np
import os
from platon.constants import R_sun, R_jup, M_jup
import pytest

from ..atmospheric_retrievals.aws_tools import configure_logging
from ..atmospheric_retrievals.aws_tools import get_config
from ..atmospheric_retrievals.platon_wrapper import _apply_factors
from ..atmospheric_retrievals.platon_wrapper import PlatonWrapper

ON_TRAVIS = os.path.expanduser('~') == '/Users/travis'


def initialize_platon_wrapper_object():
    """Return a ``PlatonWrapper`` object for use by the tests within
    this module.

    The ``PlatonWrapper`` object contains basic attributes for a simple
    example to test with.

    Returns
    -------
    pw : obj
        The ``PlatonWrapper`` object
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
    R_guess = 1.4 * R_jup
    T_guess = 1200

    # Initialize the object and set the parameters
    pw = PlatonWrapper()
    pw.set_parameters(params)
    pw.fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
    pw.fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)
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

    return pw


def test_apply_factors():
    """Test the ``_apply_factors()`` function in ``platon_wrapper``
    module.
    """

    params = {'Rs': 1.19, 'Mp': 0.73, 'Rp': 1.4}
    params = _apply_factors(params)

    assert isinstance(params, dict)
    assert params['Rs'] == 827883000.0
    assert params['Mp'] == 1.3856787e+27
    assert params['Rp'] == 100088800.0


def test_configure_logging():
    """Tests the ``configure_logging`` function in ``aws_tools``
    module.
    """

    configure_logging()
    log_file = glob.glob('logs/aws_wrapper_????-??-??-??-??.log')[0]
    assert log_file
    os.remove(log_file)
    shutil.rmtree('/logs', ignore_errors=True)


def test_get_config():
    """Tests the ``get_config`` function in ``aws_tools`` module."""

    settings = get_config()

    assert isinstance(settings, dict)
    assert 'ec2_id' in settings
    assert 'ssh_file' in settings


@pytest.mark.skipif(ON_TRAVIS, reason='Test takes too long on Travis server.  Try testing locally.')
def test_retrieve_emcee():
    """Test that the ``emcee`` method of ``platon_wrapper``
    produces results for a small example.
    """

    pw = initialize_platon_wrapper_object()
    pw.retrieve('emcee')

    assert pw.result


@pytest.mark.skipif(ON_TRAVIS, reason='Test takes too long on Travis server.  Try testing locally.')
def test_retrieve_multinest():
    """Test that the ``multinest`` method of ``platon_wrapper``
    produces results for a small example.
    """

    pw = initialize_platon_wrapper_object()
    pw.retrieve('multinest')

    assert pw.result

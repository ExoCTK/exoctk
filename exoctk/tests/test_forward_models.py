#! /usr/bin/env python

"""Tests for the ``forward_models`` module

Authors
-------

    Joe Filippazzo

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_forward_models.py
"""

import unittest
import pytest
import os

from exoctk.forward_models.forward_models import fortney_grid, generic_grid

ON_GITHUB_ACTIONS = os.path.expanduser('~') in ['/home/runner', '/Users/runner']


class TestFortney(unittest.TestCase):
    """Tests for the fortney_grid function"""
    def setUp(self):
        """Setup"""
        temp = 1000
        chem = 'noTiO'
        cloud = '0'
        pmass = '1.5'
        m_unit = 'M_jup'
        reference_radius = '1'
        r_unit = 'R_jup'
        rstar = '1'
        rstar_unit = 'R_sun'

        self.args = {'temp': temp, 'chem': chem, 'cloud': cloud, 'pmass': pmass,
                      'm_unit': m_unit, 'reference_radius': reference_radius,
                      'r_unit': r_unit, 'rstar': rstar, 'rstar_unit': rstar_unit}

    @pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Grid data not available')
    def test_model(self):
        """Test that fortney_grid can be run"""
        fig, fh, temp_out = fortney_grid(self.args, write_plot=False, write_table=False)


class TestGeneric(unittest.TestCase):
    """Tests for the generic_grid function"""
    def setUp(self):
        """Setup"""
        self.args = {'r_star' : 1, 'r_planet' : 0.1, 'gravity': '10', 'temperature': 1000,
                     'condensation': 'local', 'metallicity': '+1.0', 'c_o': '0.56', 'haze': '0010',
                     'cloud': '0.00'}

    @pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Grid data not available')
    def test_model(self):
        """Test that fortney_grid can be run"""
        fig, fh, closest_match, error_message = generic_grid(self.args, write_plot=False, write_table=False)

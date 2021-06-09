#! /usr/bin/env python

"""Tests for the ``utils`` script.

Authors
-------

    Mees Fix

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_utils.py
"""

import numpy as np
import pytest

from exoctk.utils import get_canonical_name, get_target_data

args = ['planet_name', 'planet_data']
test_data = [('HD108236f', {'canonical_name':'HD 108236 f', 'Fe/H':-0.28, 'Teff':5660.0, 'stellar_gravity':4.49, 'transit_duration':0.13625, 'RA':186.5740627, 'DEC':-51.3630519}), 
             ('NGTS-14Ab', {'canonical_name':'NGTS-14 A b', 'Fe/H':0.1, 'Teff':5187.0, 'stellar_gravity':4.2, 'transit_duration':0.09333333333333334, 'RA':328.5174799, 'DEC':-38.3774193}), 
             ('2MASSJ10193800-0948225b', {'canonical_name':'WASP-43 b', 'Fe/H':-0.01, 'Teff':4400.0, 'stellar_gravity':3.707669652604325, 'transit_duration':0.0483, 'RA':154.9081869, 'DEC':-9.8064431})]


@pytest.mark.parametrize(args, test_data)
def test_get_canoical_name(planet_name, planet_data):
    'Test that the canonical name of the planet is returned by exomast'

    canonical_name = get_canonical_name(planet_name)
    assert canonical_name == planet_data['canonical_name']


@pytest.mark.parametrize(args, test_data)
def test_get_target_data(planet_name, planet_data):
    'Test that the canonical name of the planet is returned by exomast'

    data, _ = get_target_data(planet_name)

    # these are some params that the webapp uses for it's tools
    exomast_params = ['Fe/H', 'Teff', 'stellar_gravity', 'transit_duration', 'RA', 'DEC']

    for value in exomast_params:
        assert data[value] == planet_data[value]
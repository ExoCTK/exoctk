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

from astropy.table import Table
import numpy as np
import pytest

from exoctk.utils import filter_table, get_canonical_name, get_target_data, color_gen, medfilt

args = ['planet_name', 'planet_data']

# Include 3 actual planets with data and one fake planet to xfail.
test_data = [('HD108236f', {'canonical_name':'HD 108236 f', 'Fe/H':-0.28, 'Teff':5660.0, 'stellar_gravity':4.49, 'transit_duration':0.13625, 'RA':186.5740627, 'DEC':-51.3630519}), 
             ('NGTS-14Ab', {'canonical_name':'NGTS-14 A b', 'Fe/H':0.1, 'Teff':5187.0, 'stellar_gravity':4.2, 'transit_duration':0.09333333333333334, 'RA':328.5174799, 'DEC':-38.3774193}), 
             ('2MASSJ10193800-0948225b', {'canonical_name':'WASP-43 b', 'Fe/H':-0.01, 'Teff':4400.0, 'stellar_gravity':4.49, 'transit_duration':0.0483, 'RA':154.9081869, 'DEC':-9.8064431}),
 pytest.param('djgfjhsg', {'canonical_name':'sfghsfkjg', 'Fe/H':-999, 'Teff':-999, 'stellar_gravity':-999, 'transit_duration':-999, 'RA':-999, 'DEC':-999}, marks=pytest.mark.xfail)]


@pytest.mark.parametrize(args, test_data)
def test_get_canoical_name(planet_name, planet_data):
    '''Test that the canonical name of the planet is returned by exomast'''

    canonical_name = get_canonical_name(planet_name)
    assert canonical_name == planet_data['canonical_name']


@pytest.mark.parametrize(args, test_data)
def test_get_target_data(planet_name, planet_data):
    '''Test that the canonical name of the planet is returned by exomast'''

    data, _ = get_target_data(planet_name)

    # these are some params that the webapp uses for it's tools
    exomast_params = ['Fe/H', 'Teff', 'stellar_gravity', 'transit_duration', 'RA', 'DEC']

    for value in exomast_params:
        assert data[value] == planet_data[value]


@pytest.fixture
def build_table():
    """ Fixture to build table for tests requiring table """
    column_names = ['wavelength', 'flux']
    wavelength = np.arange(3000, 5000)
    flux = np.random.rand(wavelength.shape[0])*10e-13

    table = Table(data=[wavelength, flux], names=column_names)

    return table


@pytest.mark.parametrize("operator", ['>4021', '<3856', '>=4928', '<=3740', '==4000', '4*'])
def test_filter_table(build_table, operator):
    '''test table filter function with fake table'''
    # Test wavelength sort
    test_data = filter_table(build_table, wavelength=operator)
    assert all(test_data['wavelength'] == eval("build_table[np.where(build_table['wavelength'] {})]['wavelength']".format(operator)))


@pytest.mark.parametrize("colormap", ['viridis', pytest.param('hjdsgfdhsf', marks=pytest.mark.xfail)])
def test_color_gen(colormap):
    color_gen(colormap)


medfilt_data = [(np.array([ 1,  1,  8, 12,  2, 10,  5,  2,  5, 2]), 4),
                (np.array([ 1,  1,  8, 12,  2, 10,  5,  2,  5, 2]), 4)]
@pytest.mark.parametrize(['data', 'filter_window'], medfilt_data)
def test_medfilt(data,filter_window):
    filtered_data = medfilt(data, filter_window)
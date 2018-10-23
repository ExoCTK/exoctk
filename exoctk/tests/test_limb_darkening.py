#! /usr/bin/env python

"""Tests for the ``limb_darkening_fit`` module.

Authors
-------

    Joe Filippazzo

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_limb_darkening.py
"""

import os
from pkg_resources import resource_filename

from svo_filters import Filter

from exoctk import modelgrid as mg
from exoctk.limb_darkening import limb_darkening_fit as ldf


MODELGRID = mg.ModelGrid(resource_filename('ExoCTK', 'data/core/modelgrid/'))


def test_ldc_object():
    """Test to see that an LDC object can be loaded"""
    ld_session = ldf.LDC(MODELGRID)

    assert isinstance(ld_session, ldf.LDC)


def test_ldc_calculation_no_filter():
    """Test to see if a calculation can be performed with no filter and
    that they are appended to the results table"""
    # Make the session
    ld_session = ldf.LDC(MODELGRID)

    # Run the calculations
    ld_session.calculate(Teff=4000, logg=5, FeH=0, profile='quadratic')
    ld_session.calculate(Teff=4000, logg=5, FeH=0, profile='4-parameter')

    assert len(ld_session.results) == 2


def test_ldc_calculation_filter():
    """Test to see if a calculation can be performed with a filter and
    that they are appended to the results table"""
    # Make the session
    ld_session = ldf.LDC(MODELGRID)

    # Make a filter
    filt = Filter('2MASS.H')

    # Run the calculations
    ld_session.calculate(Teff=4000, logg=5, FeH=0, profile='quadratic',
                         bandpass=filt)
    ld_session.calculate(Teff=4000, logg=5, FeH=0, profile='4-parameter',
                         bandpass=filt)

    assert len(ld_session.results) == 2


def test_ldc_calculation_interpolation():
    """Test to see if a calculation can be performed with no filter and
    an interpolated grid point and that they are appended to the results
    table"""
    # Make the session
    ld_session = ldf.LDC(MODELGRID)

    # Run the calculations
    ld_session.calculate(Teff=4023, logg=5, FeH=0, profile='quadratic')
    ld_session.calculate(Teff=4000, logg=5.1, FeH=0, profile='quadratic')

    assert len(ld_session.results) == 2

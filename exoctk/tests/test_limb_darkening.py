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

from pkg_resources import resource_filename

from svo_filters import Filter

from exoctk import modelgrid as mg
from exoctk.limb_darkening import limb_darkening_fit as ldf


# Load local modelgrid
MODELGRID = mg.ModelGrid(resource_filename('exoctk', 'data/core/modelgrid/'))


def test_ldc_calculation_filter():
    """Test to see if a calculation can be performed with a filter and
    that they are appended to the results table"""
    print('Testing LDC calculation using a filter bandpass...')

    # Make the session
    ld_session = ldf.LDC(MODELGRID)

    # Make a filter
    filt = Filter('2MASS.H')

    # Run the calculations
    ld_session.calculate(Teff=4000, logg=4.5, FeH=0, profile='quadratic',
                         bandpass=filt)
    ld_session.calculate(Teff=4000, logg=4.5, FeH=0, profile='4-parameter',
                         bandpass=filt)

    assert len(ld_session.results) == 2


def test_ldc_calculation_grism():
    """Test to see if a calculation can be performed with a grism and
    that they are appended to the results table"""
    print('Testing LDC calculation using a filter grism...')

    # Make the session
    ld_session = ldf.LDC(MODELGRID)

    # Make a filter
    filt = Filter('2MASS.H', n_bins=10)

    # Run the calculations
    ld_session.calculate(Teff=4000, logg=4.5, FeH=0, profile='quadratic',
                         bandpass=filt)
    ld_session.calculate(Teff=4000, logg=4.5, FeH=0, profile='4-parameter',
                         bandpass=filt)

    # Two profiles split into 10 measurements = 20 calculations
    assert len(ld_session.results) == 20


def test_ldc_calculation_interpolation():
    """Test to see if a calculation can be performed with no filter and
    an interpolated grid point and that they are appended to the results
    table"""
    # Make the session
    ld_session = ldf.LDC(MODELGRID)

    # Run the calculations with one parameter off grid...
    ld_session.calculate(Teff=4023, logg=4.5, FeH=0, profile='quadratic')
    ld_session.calculate(Teff=4000, logg=4.1, FeH=0, profile='quadratic')
    ld_session.calculate(Teff=4000, logg=4.5, FeH=-0.1, profile='quadratic')

    # ...and two parameters off grid...
    ld_session.calculate(Teff=4023, logg=4.1, FeH=0, profile='quadratic')
    ld_session.calculate(Teff=4023, logg=4.5, FeH=-0.1, profile='quadratic')
    ld_session.calculate(Teff=4000, logg=4.1, FeH=-0.1, profile='quadratic')

    # ...and all three parameters off grid.
    ld_session.calculate(Teff=4023, logg=4.1, FeH=-0.1, profile='quadratic')

    assert len(ld_session.results) == 7


def test_ldc_calculation_no_filter():
    """Test to see if a calculation can be performed with no filter and
    that they are appended to the results table"""
    print('Testing LDC calculation with no filter bandpass...')

    # Make the session
    ld_session = ldf.LDC(MODELGRID)

    # Run the calculations
    ld_session.calculate(Teff=4000, logg=4.5, FeH=0, profile='quadratic')
    ld_session.calculate(Teff=4000, logg=4.5, FeH=0, profile='4-parameter')

    assert len(ld_session.results) == 2


def test_ldc_object():
    """Test to see that an LDC object can be loaded"""
    print('Testing LDC object creation...')

    ld_session = ldf.LDC(MODELGRID)

    assert isinstance(ld_session, ldf.LDC)


def test_ldc_plot():
    """Test that the LDC plots work"""
    print('Testing LDC object plotting...')

    ld_session = ldf.LDC(MODELGRID)

    # Run the calculations
    ld_session.calculate(Teff=4000, logg=4.5, FeH=0, profile='quadratic')

    # Regular plot
    fig = ld_session.plot()
    assert str(type(fig)) == "<class 'bokeh.plotting._figure.figure'>"

    # Tabbed plot
    fig = ld_session.plot_tabs()
    assert str(type(fig)) in ["<class 'bokeh.models.layouts.TabPanel'>", "<class 'bokeh.models.layouts.Tabs'>"]

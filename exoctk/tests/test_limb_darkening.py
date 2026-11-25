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

import pytest
from svo_filters import Filter

from exoctk import modelgrid as mg
from exoctk.limb_darkening import limb_darkening_fit as ldf
from exoctk.pkgdata import resource_filename


# Load local modelgrid
MODELGRID = mg.ModelGrid(resource_filename('exoctk', 'data/core/modelgrid/'))


@pytest.fixture
def limb_darkening_client(monkeypatch):
    web = pytest.importorskip(
        "exoctk.exoctk_app.app_exoctk",
        reason="Form tests require the optional web stack")
    monkeypatch.setitem(web.app_exoctk.config, "TESTING", True)
    monkeypatch.setitem(web.app_exoctk.config, "WTF_CSRF_ENABLED", False)
    monkeypatch.setattr(
        web, "render_template",
        lambda template, **context: {
            "template": template,
            "n_bins": context["form"].n_bins.data,
        })
    with web.app_exoctk.test_client() as client:
        yield client


@pytest.mark.parametrize("bandpass, expected_bins", [
    ("MIRI.F1500W.CLEARP.IMAGER", 1),
    ("2MASS.H", 1),
    ("TESS_TESS.Red", 1),
    ("tophat", 1),
    ("MIRI.CLEAR.P750L.LRSSLITLESS", 30),
    ("MIRI.CLEAR.SHORT.CH1", 30),
    ("WFC3_IR.G141", 30),
])
def test_limb_darkening_bandpass_channel_defaults(
        limb_darkening_client, bandpass, expected_bins):
    """Selecting a bandpass replaces the previous channel count."""
    response = limb_darkening_client.post("/limb_darkening", data={
        "filter_submit": "Filter Selected",
        "bandpass": bandpass,
        "n_bins": "7",
    })
    assert response.status_code == 200
    assert response.get_json()["n_bins"] == expected_bins


def test_limb_darkening_imaging_channel_override(limb_darkening_client):
    """An explicit channel count survives a calculation validation error."""
    selected = limb_darkening_client.post("/limb_darkening", data={
        "filter_submit": "Filter Selected",
        "bandpass": "MIRI.F1500W.CLEARP.IMAGER",
        "n_bins": "30",
    })
    assert selected.get_json()["n_bins"] == 1

    response = limb_darkening_client.post("/limb_darkening", data={
        "calculate_submit": "Calculate Coefficients",
        "bandpass": "MIRI.F1500W.CLEARP.IMAGER",
        "n_bins": "4",
    })
    assert response.status_code == 200
    assert response.get_json() == {
        "template": "limb_darkening.html",
        "n_bins": 4,
    }


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

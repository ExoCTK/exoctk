#! /usr/bin/env python

"""Tests for various modules within the ``contam_visibility``
subpackage.

Authors
-------

    Joe Filippazzo
    Matthew Bourque
    Ben Falk

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_contam_visibility.py
"""

import json
import os
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
import pysiaf

from pandas import DataFrame
from astropy.io import fits
from astropy.table import MaskedColumn, Table
import astropy.units as u

from exoctk.contam_visibility import field_simulator
from exoctk.contam_visibility import contamination_figure
from exoctk.contam_visibility import miri_lrs
from exoctk.contam_visibility import make_miri_lrs_traces
from exoctk.contam_visibility import precompute
from exoctk.contam_visibility import resolve
from exoctk.contam_visibility import new_vis_plot
from exoctk.contam_visibility import contamination_figure
from exoctk.contam_visibility.modes import CONTAM_VISIBILITY_MODES

# Determine if tests are being run on Github Actions
ON_GITHUB_ACTIONS = '/home/runner' in os.path.expanduser('~') or '/Users/runner' in os.path.expanduser('~')


def classification_row(mask_dsc=False, mask_columns=(), **values):
    """Build one synthetic Gaia row for source-classification tests."""

    defaults = {
        'classprob_dsc_combmod_star': np.nan,
        'classprob_dsc_combmod_galaxy': np.nan,
        'classprob_dsc_combmod_quasar': np.nan,
        'parallax': 1.,
        'astrometric_excess_noise': 0.,
        'phot_bp_rp_excess_factor': 1.,
        'bp_rp': 1.,
    }
    defaults.update(values)
    table = Table({name: [value] for name, value in defaults.items()},
                  masked=True)
    if mask_dsc:
        for name in (
                'classprob_dsc_combmod_star',
                'classprob_dsc_combmod_galaxy',
                'classprob_dsc_combmod_quasar'):
            table[name].mask[0] = True
    for name in mask_columns:
        table[name].mask[0] = True
    return table[0]


def gaia_source_table(source_ids, excess_noise=0.):
    """Build a minimal Gaia result table for find_sources tests."""

    count = len(source_ids)
    return Table({
        'source_id': source_ids,
        'ra': 10. + np.arange(count) * 0.001,
        'dec': np.full(count, 20.),
        'pmra': np.zeros(count),
        'pmdec': np.zeros(count),
        'ref_epoch': np.full(count, 2016.),
        'phot_g_mean_flux': np.full(count, 100.),
        'bp_rp': np.ones(count),
        'parallax': np.full(count, 0.1),
        'astrometric_excess_noise': np.full(count, excess_noise),
        'phot_bp_rp_excess_factor': np.ones(count),
        'classprob_dsc_combmod_star': np.full(count, np.nan),
        'classprob_dsc_combmod_galaxy': np.full(count, np.nan),
        'classprob_dsc_combmod_quasar': np.full(count, np.nan),
    })


def test_new_vis_plot():
    """Tests the `new_vis_plot.py` module"""
    ra, dec = '24.3544618', '-45.6777937' # WASP-18
    table = new_vis_plot.get_exoplanet_positions(ra, dec)

    assert isinstance(table, DataFrame)

    plt = new_vis_plot.build_visibility_plot('WASP-18b', 'NIRISS', ra, dec)

    assert str(type(plt)) == "<class 'bokeh.plotting._figure.figure'>"


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Need access to trace data FITS files.  Please try running locally')
def test_field_simulation():
    """Tests the ``field_simulation`` function in the ``field_simulator`` module"""

    ra = '04 25 29.0162'
    dec = '-30 36 01.603'
    aperture = 'NIS_SUBSTRIP256'

    targframe, starcube, results = field_simulator.field_simulation(ra, dec, aperture)

    assert isinstance(targframe, (np.ndarray, list)) and isinstance(starcube, (np.ndarray, list))


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Need access to trace data FITS files.  Please try running locally')
def test_precomputed_field_simulation():
    """Tests the ``field_simulation`` function in the ``field_simulator`` module using precomputed cache"""

    # Should pass
    targname = 'TRAPPIST-1'
    aperture = 'NIS_SUBSTRIP256'
    target_db = Path(__file__).parent / "test_data" / "NIS_SUBSTRIP256_test_db.h5"

    targframe, starcube, results = field_simulator.field_simulation(targname=targname, aperture=aperture, target_db=target_db)

    assert isinstance(targframe, (np.ndarray, list)) and isinstance(starcube, (np.ndarray, list))

    # Should still pass if target is not in the database
    bad_targname = 'WASP-18'
    targframe, starcube, results = field_simulator.field_simulation(targname=bad_targname, aperture=aperture, target_db=target_db)

    assert isinstance(targframe, (np.ndarray, list)) and isinstance(starcube, (np.ndarray, list))


def test_miri_position_angle_failure_is_raised(monkeypatch):
    """LRS propagates a PA calculation failure like every other mode."""

    aperture_name = miri_lrs.APERTURE
    full_name = field_simulator.APERTURES[aperture_name]['full']

    class FakeFullAperture:
        @staticmethod
        def corners(frame):
            assert frame == 'det'
            return np.array([0., 1.]), np.array([0., 1.])

    class FakeScienceAperture:
        XSciSize = 2
        YSciSize = 2

    class FakeSiaf:
        apertures = {
            full_name: FakeFullAperture(),
            aperture_name: FakeScienceAperture(),
        }

    positions = DataFrame({
        'V3PA_min_pa_angle': [0.],
        'V3PA_nominal_angle': [0.],
        'V3PA_max_pa_angle': [0.],
    })
    monkeypatch.delenv('EXOCTK_CONTAM_CACHE', raising=False)
    monkeypatch.setattr(field_simulator, 'check_for_data', lambda *args: None)
    monkeypatch.setattr(field_simulator.pysiaf, 'Siaf', lambda inst: FakeSiaf())
    monkeypatch.setattr(field_simulator, 'find_sources', lambda *args, **kwargs: Table())
    monkeypatch.setattr(
        field_simulator, 'get_exoplanet_positions',
        lambda *args, **kwargs: positions)
    monkeypatch.setattr(
        field_simulator, 'calculation_v3pas',
        lambda aperture, observable: np.array([0, 1]))

    def calculate(pa, **kwargs):
        if pa == 1:
            raise RuntimeError('synthetic PA failure')
        return {
            'pa': int(pa),
            'target_traces': [np.ones((2, 2))],
            'contaminants': np.zeros((2, 2)),
        }

    monkeypatch.setattr(field_simulator, 'calc_v3pa', calculate)

    with pytest.raises(RuntimeError, match='synthetic PA failure'):
        field_simulator.field_simulation(
            ra=10., dec=20., aperture=aperture_name)


def test_resolve_target():
    """Tests the ``resolve_target`` function in the ``resolve`` module"""

    ra, dec = resolve.resolve_target('Wasp-18 b')

    assert ra == 24.3544618
    assert dec == -45.6777937


@pytest.mark.parametrize(('probabilities', 'expected'), [
    ((0.9, 0.05, 0.05), 'STAR'),
    ((0.05, 0.9, 0.05), 'GALAXY'),
    ((0.05, 0.05, 0.9), 'STAR'),
    ((0.3, 0.5, 0.2), 'GALAXY'),
    ((0.5, 0.5, 0.), 'STAR'),
    ((0.4, 0.4, 0.2), 'STAR'),
])
def test_classify_source_uses_gaia_dsc(probabilities, expected):
    """Gaia DSC drives the binary contamination-rendering classification."""

    star, galaxy, quasar = probabilities
    row = classification_row(
        classprob_dsc_combmod_star=star,
        classprob_dsc_combmod_galaxy=galaxy,
        classprob_dsc_combmod_quasar=quasar,
        parallax=0.)

    assert field_simulator.classify_source(row) == expected


def test_classify_source_masked_dsc_uses_fallback():
    """Strong excess noise remains a legacy extension proxy."""

    row = classification_row(
        mask_dsc=True, parallax=0.1, astrometric_excess_noise=2.)

    assert field_simulator.classify_source(row) == 'GALAXY'


def test_classify_source_bp_rp_excess_uses_fallback_proxy():
    """Strong BP/RP excess remains a legacy extension proxy."""

    row = classification_row(
        mask_dsc=True, parallax=np.nan, bp_rp=1.,
        phot_bp_rp_excess_factor=2.)

    assert field_simulator.classify_source(row) == 'GALAXY'


def test_classify_source_one_masked_dsc_probability_defaults_to_star():
    """One missing DSC probability makes the full DSC result unusable."""

    row = classification_row(
        mask_columns=('classprob_dsc_combmod_quasar',),
        classprob_dsc_combmod_star=0.1,
        classprob_dsc_combmod_galaxy=0.9,
        classprob_dsc_combmod_quasar=0.,
        parallax=np.nan)

    assert field_simulator.classify_source(row) == 'STAR'


def test_classify_source_nan_and_masked_values_match():
    """NaN and masked DSC and parallax values have the same safe default."""

    nan_row = classification_row(parallax=np.nan)
    masked_row = classification_row(
        mask_dsc=True, mask_columns=('parallax',))

    assert field_simulator.classify_source(nan_row) == 'STAR'
    assert field_simulator.classify_source(masked_row) == 'STAR'


def test_classify_source_missing_parallax_is_not_a_galaxy():
    """Missing parallax alone must not suppress a dispersed contaminant."""

    row = classification_row(mask_dsc=True, parallax=np.nan)

    assert field_simulator.classify_source(row) == 'STAR'


def test_classify_source_low_parallax_is_not_a_galaxy():
    """Low finite parallax alone must not imply an extended source."""

    row = classification_row(mask_dsc=True, parallax=0.1)

    assert field_simulator.classify_source(row) == 'STAR'


def test_classify_source_regression_gaia_3910744542517589888():
    """The independently confirmed point source remains a STAR."""

    row = classification_row(
        classprob_dsc_combmod_star=0.997072,
        classprob_dsc_combmod_galaxy=0.000248,
        classprob_dsc_combmod_quasar=0.001099,
        parallax=np.nan)

    assert field_simulator.classify_source(row) == 'STAR'


@pytest.mark.parametrize('source_type', ['STAR', 'GALAXY'])
def test_classify_source_preserves_valid_upstream_type(source_type):
    """A valid upstream type survives when Gaia DSC is unavailable."""

    row = classification_row(type=source_type, parallax=10.)

    assert field_simulator.classify_source(row) == source_type


@pytest.mark.parametrize(('sdss_type', 'expected'), [
    ('STAR', 'STAR'),
    ('GALAXY', 'GALAXY'),
    ('QSO', 'STAR'),
    (' qso ', 'STAR'),
    ('star', 'STAR'),
    ('Galaxy', 'GALAXY'),
    ('', 'STAR'),
    ('<masked>', 'STAR'),
    (None, 'STAR'),
    ('UNKNOWN', 'STAR'),
])
def test_find_sources_normalizes_explicit_sdss_type(
        monkeypatch, sdss_type, expected):
    """Production flow normalizes usable SDSS classes by rendering type."""

    stars = gaia_source_table([1])
    if sdss_type is None:
        xmatch = Table(names=('source_id', 'spCl'), dtype=('i8', 'U10'))
    elif sdss_type == '<masked>':
        xmatch = Table({'source_id': [1], 'spCl': ['STAR']}, masked=True)
        xmatch['spCl'].mask[0] = True
    else:
        xmatch = Table({'source_id': [1], 'spCl': [sdss_type]})
    monkeypatch.setattr(
        field_simulator.GAIA_TAP, 'query_region',
        lambda *args, **kwargs: stars.copy())
    monkeypatch.setattr(
        field_simulator.XMatch, 'query',
        lambda *args, **kwargs: xmatch)

    result = field_simulator.find_sources(10., 20., pm_corr=False)

    assert result['type'][0] == expected


def test_find_sources_associates_sdss_types_by_source_id(monkeypatch):
    """SDSS types remain attached when Gaia IDs are out of key order."""

    stars = gaia_source_table([30, 10, 20])
    xmatch = Table({
        'source_id': [20, 10],
        'spCl': ['GALAXY', 'STAR'],
    })
    monkeypatch.setattr(
        field_simulator.GAIA_TAP, 'query_region',
        lambda *args, **kwargs: stars.copy())
    monkeypatch.setattr(
        field_simulator.XMatch, 'query',
        lambda *args, **kwargs: xmatch)

    result = field_simulator.find_sources(10., 20., pm_corr=False)

    assert list(result['source_id']) == [30, 10, 20]
    assert list(result['type']) == ['STAR', 'STAR', 'GALAXY']
    assert len(result) == len(stars)


@pytest.mark.parametrize(('sdss_types', 'expected'), [
    (['STAR', ' qso '], 'STAR'),
    (['STAR', 'GALAXY'], 'GALAXY'),
])
def test_find_sources_handles_duplicate_sdss_matches(
        monkeypatch, sdss_types, expected):
    """Duplicate matches are preserved only when usable classes agree."""

    stars = gaia_source_table([1], excess_noise=2.)
    xmatch = Table({
        'source_id': [1] * len(sdss_types),
        'spCl': sdss_types,
    })
    monkeypatch.setattr(
        field_simulator.GAIA_TAP, 'query_region',
        lambda *args, **kwargs: stars.copy())
    monkeypatch.setattr(
        field_simulator.XMatch, 'query',
        lambda *args, **kwargs: xmatch)

    result = field_simulator.find_sources(10., 20., pm_corr=False)

    assert result['type'][0] == expected
    assert len(result) == 1


def test_find_sources_identifies_high_proper_motion_target(monkeypatch):
    """Gaia query order must not replace a moving target with a field star."""

    # At the Gaia 2016 epoch the faint source is closer to the supplied J2000
    # coordinate than HD 189733. Back-propagating the proper motion identifies
    # the bright second row as the intended target.
    stars = Table({
        'source_id': [1827242816182176512, 1827242816201846144],
        'ra': [300.18318127620523, 300.1821218062407],
        'dec': [22.711090387243765, 22.709741105168273],
        'pmra': [0., -3.208339784864691],
        'pmdec': [0., -250.32333085817578],
        'ref_epoch': [2016., 2016.],
        'phot_g_mean_flux': [402.513622149291, 2.0118629422041267e7],
        'bp_rp': [1.5, 1.0959163],
        'parallax': [1., 50.],
        'astrometric_excess_noise': [0., 0.],
        'phot_bp_rp_excess_factor': [1., 1.],
    })

    monkeypatch.setattr(
        field_simulator.GAIA_TAP, 'query_region',
        lambda *args, **kwargs: stars.copy())

    def unavailable_xmatch(*args, **kwargs):
        raise RuntimeError('XMatch unavailable in unit test')

    monkeypatch.setattr(
        field_simulator.XMatch, 'query', unavailable_xmatch)
    result = field_simulator.find_sources(
        300.1821375, 22.7108528, target_date=2026)

    assert result['source_id'][0] == 1827242816201846144
    assert result['fluxscale'][0] == pytest.approx(1.)
    assert result['fluxscale'][1] < 1.e-4
    assert result['distance'][0] == pytest.approx(0.)


def test_source_projection_batches_pysiaf_detector_transform(monkeypatch):
    """Batch projection preserves the scalar placement equations."""

    stars = Table({
        'ra': [0., 1.25, -2.75],
        'dec': [0., 3.5, -4.5],
        'xtel': np.zeros(3), 'ytel': np.zeros(3),
        'xdet': np.zeros(3), 'ydet': np.zeros(3),
        'xsci': [100., 0., 0.], 'ysci': [200., 0., 0.],
        'xord0': [210, 0, 0], 'yord0': [410, 0, 0],
        'xord1': np.zeros(3), 'yord1': np.zeros(3),
    })
    aper = {
        'c0x0': 100., 'c0y0': 200., 'c1x0': -0.1, 'c1y0': 0.2,
        'c1x1': -0.03, 'c1y1': 0.12, 'c2y1': -0.011,
        'xord0to1': -500, 'yord0to1': 75,
    }

    def fake_sky_to_tel(attitude, ra, dec):
        return np.asarray(ra) * u.arcsec, np.asarray(dec) * u.arcsec

    class FakeAperture:
        def __init__(self):
            self.tel_to_det_calls = []

        def tel_to_det(self, xtel, ytel):
            self.tel_to_det_calls.append((np.asarray(xtel), np.asarray(ytel)))
            return np.asarray(xtel) + 10., np.asarray(ytel) + 20.

        def det_to_sci(self, xdet, ydet):
            return np.asarray(xdet) * 2., np.asarray(ydet) * 3.

    monkeypatch.setattr(
        field_simulator.pysiaf.utils.rotations, 'sky_to_tel', fake_sky_to_tel)
    aperture = FakeAperture()

    field_simulator._project_sources_to_detector(
        np.eye(3), stars, aperture, aper)

    assert len(aperture.tel_to_det_calls) == 1
    np.testing.assert_allclose(aperture.tel_to_det_calls[0][0], [1.25, -2.75])
    np.testing.assert_allclose(aperture.tel_to_det_calls[0][1], [3.5, -4.5])

    # These are the scalar formulas formerly evaluated once per source.
    for index in range(1, len(stars)):
        xsci = (stars['ra'][index] + 10.) * 2.
        ysci = (stars['dec'][index] + 20.) * 3.
        xord0 = int(xsci + aper['c0x0']
                    + aper['c1x0'] * (stars['xsci'][0] - xsci))
        yord0 = int(ysci + aper['c0y0']
                    + aper['c1y0'] * (stars['ysci'][0] - ysci))
        xord1 = xord0 + aper['xord0to1'] + int(
            aper['c1x1'] * (stars['xord0'][0] - xord0))
        yord1 = yord0 + aper['yord0to1'] + int(
            aper['c1y1'] * (stars['yord0'][0] - yord0)) + int(
            aper['c2y1'] * (stars['xord0'][0] - xord0))

        assert stars['xsci'][index] == xsci
        assert stars['ysci'][index] == ysci
        assert stars['xord0'][index] == xord0
        assert stars['yord0'][index] == yord0
        assert stars['xord1'][index] == xord1
        assert stars['yord1'][index] == yord1


def test_batched_source_projection_matches_scalar_pysiaf():
    """Real pySIAF vector inputs retain the previous scalar placements."""

    aperture = field_simulator.pysiaf.Siaf('NIRISS')['NIS_SUBSTRIP96']
    aper = field_simulator.APERTURES[aperture.AperName]
    stars = Table({
        'ra': [300.1821375, 300.1830, 300.1814],
        'dec': [22.7108528, 22.7113, 22.7100],
        'xtel': np.zeros(3), 'ytel': np.zeros(3),
        'xdet': np.zeros(3), 'ydet': np.zeros(3),
        'xsci': np.zeros(3), 'ysci': np.zeros(3),
        'xord0': np.zeros(3), 'yord0': np.zeros(3),
        'xord1': np.zeros(3), 'yord1': np.zeros(3),
    })

    stars['xdet'][0], stars['ydet'][0] = aperture.reference_point('det')
    stars['xtel'][0], stars['ytel'][0] = aperture.det_to_tel(
        stars['xdet'][0], stars['ydet'][0])
    stars['xsci'][0], stars['ysci'][0] = aperture.det_to_sci(
        stars['xdet'][0], stars['ydet'][0])
    stars['xord0'][0] = int(stars['xsci'][0] + aper['c0x0'])
    stars['yord0'][0] = int(stars['ysci'][0] + aper['c0y0'])

    v3pa = 240.
    apa = (v3pa + aperture.V3IdlYAngle) % 360
    attitude = field_simulator.pysiaf.utils.rotations.attitude_matrix(
        stars['xtel'][0], stars['ytel'][0], stars['ra'][0], stars['dec'][0],
        apa)

    expected = stars.copy()
    for index in range(1, len(expected)):
        v2, v3 = field_simulator.pysiaf.utils.rotations.sky_to_tel(
            attitude, expected['ra'][index], expected['dec'][index])
        expected['xtel'][index] = v2.to_value(u.arcsec)
        expected['ytel'][index] = v3.to_value(u.arcsec)
        expected['xdet'][index], expected['ydet'][index] = aperture.tel_to_det(
            expected['xtel'][index], expected['ytel'][index])
        expected['xsci'][index], expected['ysci'][index] = aperture.det_to_sci(
            expected['xdet'][index], expected['ydet'][index])
        expected['xord0'][index] = int(
            expected['xsci'][index] + aper['c0x0'] + aper['c1x0']
            * (expected['xsci'][0] - expected['xsci'][index]))
        expected['yord0'][index] = int(
            expected['ysci'][index] + aper['c0y0'] + aper['c1y0']
            * (expected['ysci'][0] - expected['ysci'][index]))
        expected['xord1'][index] = expected['xord0'][index] + aper['xord0to1'] + int(
            aper['c1x1'] * (expected['xord0'][0] - expected['xord0'][index]))
        expected['yord1'][index] = expected['yord0'][index] + aper['yord0to1'] + int(
            aper['c1y1'] * (expected['yord0'][0] - expected['yord0'][index])) + int(
            aper['c2y1'] * (expected['xord0'][0] - expected['xord0'][index]))

    field_simulator._project_sources_to_detector(attitude, stars, aperture, aper)

    for name in ('xtel', 'ytel', 'xdet', 'ydet', 'xsci', 'ysci'):
        np.testing.assert_allclose(stars[name], expected[name])
    for name in ('xord0', 'yord0', 'xord1', 'yord1'):
        np.testing.assert_array_equal(stars[name], expected[name])


def test_trace_templates_are_cached_across_position_angles(monkeypatch):
    """Repeated source templates avoid repeated trace-file reads."""

    calls = []
    monkeypatch.setenv('EXOCTK_DATA', '/synthetic-data')
    monkeypatch.setattr(
        field_simulator.glob, 'glob',
        lambda path: ['/synthetic-data/exoctk_contam/traces/'
                      'NIS_SUBSTRIP256/trace_5000.fits'])

    def synthetic_trace(filename, ext=0):
        calls.append((filename, ext))
        return np.full((2, 2), ext + 1., dtype=float)

    monkeypatch.setattr(field_simulator.fits, 'getdata', synthetic_trace)
    field_simulator._get_trace_cached.cache_clear()
    try:
        first = field_simulator.get_trace('NIS_SUBSTRIP256', 5000., 'STAR')
        second = field_simulator.get_trace('NIS_SUBSTRIP256', 5000., 'STAR')
    finally:
        field_simulator._get_trace_cached.cache_clear()

    assert len(calls) == 3
    assert all(np.array_equal(before, after)
               for before, after in zip(first, second))
    assert all(not trace.flags.writeable for trace in first)


def test_order_zero_templates_are_cached_across_position_angles(monkeypatch):
    """Order-zero templates are also reused when rendering crowded fields."""

    calls = []
    monkeypatch.setenv('EXOCTK_DATA', '/synthetic-data')
    monkeypatch.setattr(
        field_simulator.glob, 'glob',
        lambda path: ['/synthetic-data/exoctk_contam/order0/NIS_order0_5000.npy'])

    def synthetic_order_zero(filename):
        calls.append(filename)
        return np.ones((2, 2), dtype=float)

    monkeypatch.setattr(field_simulator.np, 'load', synthetic_order_zero)
    field_simulator._get_order0_cached.cache_clear()
    try:
        first = field_simulator.get_order0('NIS_SUBSTRIP256', 5000., 'STAR')
        second = field_simulator.get_order0('NIS_SUBSTRIP256', 5000., 'STAR')
    finally:
        field_simulator._get_order0_cached.cache_clear()

    assert calls == [
        '/synthetic-data/exoctk_contam/order0/NIS_order0_5000.npy']
    assert np.array_equal(first, second)
    assert not first.flags.writeable
def test_contamination_slider_uses_percent_and_common_display_cap():
    """All supported modes share a 0--10% percent-based display."""

    fractions = np.full((360, 3), 0.125)
    plot = contamination_figure.contam_slider_plot(
        [fractions], badPA_list=[0])
    spectrum_plot, slider_row, pa_plot = plot.children
    slider = slider_row.children[1]
    line_renderer = spectrum_plot.renderers[0]
    shading_renderer = pa_plot.renderers[0]
    threshold_renderer = pa_plot.renderers[-1]

    assert spectrum_plot.yaxis.axis_label == 'Contamination (%)'
    assert pa_plot.yaxis.axis_label == 'Mean Contamination (%)'
    assert spectrum_plot.y_range.end == 10
    assert pa_plot.y_range.end == 10
    assert slider.title == 'V3 Position Angle'
    assert np.all(line_renderer.data_source.data['contam1'] == 12.5)
    assert shading_renderer.data_source.data['top'][0] == 10
    assert threshold_renderer.data_source.data['top'][0] == 10


def test_contamination_slider_fill_stops_at_trace_boundary():
    """Filled order regions do not extend beyond their valid columns."""

    fractions = np.full((360, 5), 0.01)
    fractions[:, 3:] = np.nan
    plot = contamination_figure.contam_slider_plot(
        [fractions], badPA_list=[])
    spectrum_plot = plot.children[0]
    source = spectrum_plot.renderers[0].data_source
    area = spectrum_plot.renderers[1].glyph

    assert area.y1 == 'baseline1'
    np.testing.assert_array_equal(
        np.isfinite(source.data['baseline1']),
        np.isfinite(source.data['contam1']))
    assert np.all(source.data['baseline1'][:3] == 0)
    assert np.isnan(source.data['baseline1'][3:]).all()


def test_contamination_slider_breaks_mean_curve_at_bad_pas():
    """Unobservable PAs are gaps rather than zeroes in the mean curve."""

    fractions = np.full((360, 3), 0.01)
    plot = contamination_figure.contam_slider_plot(
        [fractions], badPA_list=[10, 11])
    pa_plot = plot.children[2]
    mean_renderer = pa_plot.renderers[1]
    mean_data = mean_renderer.data_source.data[mean_renderer.glyph.y]

    assert mean_data[9] == 1
    assert np.isnan(mean_data[10:12]).all()
    assert mean_data[12] == 1


def test_observable_pa_ranges_use_v3pa_not_instrument_angle():
    """Visibility selection must not apply the SIAF offset a second time."""

    v3pas = np.arange(87., 107.)
    position_table = DataFrame({
        'V3PA_min_pa_angle': v3pas,
        'V3PA_nominal_angle': v3pas,
        'V3PA_max_pa_angle': v3pas,
        'NIRISS_min_pa_angle': v3pas + 2.0,
        'NIRISS_nominal_angle': v3pas + 2.0,
        'NIRISS_max_pa_angle': v3pas + 2.0,
    })

    position_angles, bounds, sampled = (
        field_simulator.observable_v3pa_ranges(position_table))

    assert bounds == [(87, 106)]
    assert np.array_equal(position_angles, np.arange(87, 107))
    assert np.array_equal(sampled, np.arange(87, 107))


@pytest.mark.parametrize('instrument, aperture_name', [
    ('NIRISS', 'NIS_SUBSTRIP256'),
    ('NIRCAM', 'NRCA5_41STRIPE1_DHS_F444W'),
])
def test_siaf_aperture_pa_is_derived_once_from_v3pa(
        instrument, aperture_name):
    """SOSS and DHS geometry receives the SIAF aperture PA exactly once."""

    aperture = field_simulator.pysiaf.Siaf(instrument)[aperture_name]
    v3pa = 123.4

    assert field_simulator.aperture_pa_from_v3pa(v3pa, aperture) == (
        pytest.approx((v3pa + aperture.V3IdlYAngle) % 360))
def test_find_sources_ignores_invalid_gaia_fluxes(monkeypatch):
    """Invalid Gaia G fluxes cannot become normalized field sources."""

    stars = Table({
        'source_id': [1, 2, 3, 4, 5, 6],
        'ra': [10., 10.001, 10.002, 10.003, 10.004, 10.005],
        'dec': [20., 20., 20., 20., 20., 20.],
        'pmra': [0.] * 6,
        'pmdec': [0.] * 6,
        'ref_epoch': [2016.] * 6,
        'bp_rp': [1.] * 6,
        'parallax': [1.] * 6,
        'astrometric_excess_noise': [0.] * 6,
        'phot_bp_rp_excess_factor': [1.] * 6,
    })
    stars['phot_g_mean_flux'] = MaskedColumn(
        [1000., 1., 0., -1., np.nan, np.inf],
        mask=[False, True, False, False, False, False])
    monkeypatch.setattr(
        field_simulator.GAIA_TAP, 'query_region',
        lambda *args, **kwargs: stars.copy())
    monkeypatch.setattr(
        field_simulator.XMatch, 'query',
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError()))

    result = field_simulator.find_sources(10., 20., target_date=2026)

    assert list(result['source_id']) == [1]
    assert result['fluxscale'][0] == pytest.approx(1.)


def test_custom_source_flux_validation_rejects_invalid_target_and_sources():
    """The rendering boundary validates direct custom source tables."""

    sources = Table()
    sources['fluxscale'] = MaskedColumn(
        [1., 1., 0., -1., np.nan, np.inf],
        mask=[False, True, False, False, False, False])
    filtered = field_simulator._filter_valid_flux_sources(
        sources, 'fluxscale', 'flux scale')
    assert len(filtered) == 1
    assert filtered['fluxscale'][0] == pytest.approx(1.)

    sources['fluxscale'][0] = 0.
    with pytest.raises(ValueError, match='target does not have a valid flux scale'):
        field_simulator._filter_valid_flux_sources(
            sources, 'fluxscale', 'flux scale')


def test_fraction_contaminated_returns_nan_for_empty_channel_without_warning(
        monkeypatch):
    """Expected empty extraction channels do not globally hide warnings."""

    target = np.array([[1., 0.], [1., 0.]])
    contaminants = np.zeros((1, 2, 2))
    monkeypatch.setattr(
        field_simulator, 'NIRISS_SOSS_trace_mask',
        lambda aperture: [np.ones_like(target)])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result = field_simulator.fraction_contaminated(
            'NIS_SUBSTRIP256', [target], contaminants)[0]

    assert result.shape == (1, 2)
    assert np.allclose(result[:, 0], 0.)
    assert np.isnan(result[:, 1]).all()
    assert not any('Mean of empty slice' in str(warning.message)
                   for warning in caught)


@pytest.mark.parametrize('aperture, mask_function', [
    ('NIS_SUBSTRIP96', 'NIRISS_SOSS_trace_mask'),
    ('NRCA5_41STRIPE1_DHS_F444W', 'NIRCam_DHS_trace_mask'),
])
def test_fraction_contaminated_ignores_flux_outside_extraction(
        monkeypatch, aperture, mask_function):
    """Off-mask sources do not dilute SOSS or DHS contamination."""

    mask = np.array([[1.], [1.], [0.]])
    target = np.array([[1.], [1.], [0.]])
    contaminants = np.array([
        [[1.], [0.], [0.]],
        [[1.], [0.], [100.]],
    ])
    monkeypatch.setattr(
        field_simulator, mask_function, lambda aperture: [mask])

    result = field_simulator.fraction_contaminated(
        aperture, [target], contaminants)[0]

    assert np.allclose(result, 0.25)


@pytest.mark.parametrize('aperture', [
    'NIS_SUBSTRIP96',
    'NIS_SUBSTRIP256',
    'NRCA5_41STRIPE1_DHS_F322W2',
    'NRCA5_41STRIPE1_DHS_F444W',
    'MIRIM_SLITLESSPRISM_IP',
])
def test_contamination_supported(aperture):
    """The web interface enables contamination only for its supported modes."""

    assert field_simulator.contamination_supported(aperture)


@pytest.mark.parametrize('aperture', [
    'NIS_SOSSFULL',
    'NRCA5_GRISM256_F322W2',
    'NRCA5_GRISM256_F444W',
    'MIRIM_SLITLESSPRISM',
    'NIRSpec',
    None,
])
def test_contamination_not_supported(aperture):
    """Visibility-only modes must not run contamination calculations."""

    assert not field_simulator.contamination_supported(aperture)


@pytest.mark.parametrize(('aperture', 'expected'), [
    ('NIS_SUBSTRIP96', False),
    ('NIS_SUBSTRIP256', True),
    ('NRCA5_41STRIPE1_DHS_F322W2', True),
])
def test_order2_contamination_availability(aperture, expected):
    """Only SUBSTRIP96 omits its uncalculated SOSS Order 2 output."""

    assert contamination_figure.has_order2_contamination(aperture) is expected


def test_substrip96_contamination_plot_omits_order2(monkeypatch):
    """SUBSTRIP96 renders only its calculated Order 1 contamination panel."""

    class Cube:
        shape = (362, 2048, 96)

    monkeypatch.setattr(
        contamination_figure, 'nirissContam',
        lambda cube, lam_file: (np.zeros((2048, 360)), np.zeros((2048, 360))))

    plot = contamination_figure.contam(Cube(), 'NIS_SUBSTRIP96')

    assert len(plot.children) == 2

def test_substrip96_slider_omits_uncalculated_orders():
    """The web slider exposes only the calculated SUBSTRIP96 order."""

    fractions = [np.full((360, 3), value) for value in (0.01, 0.02, 0.03)]
    plot = contamination_figure.contam_slider_plot(
        fractions, badPA_list=[], instrument='NIS_SUBSTRIP96')
    spectrum_plot, _, pa_plot = plot.children

    # One order produces one line and one filled area in the upper plot, and
    # one mean curve plus one threshold region in the lower plot.
    assert len(spectrum_plot.renderers) == 2
    assert len(pa_plot.renderers) == 3
    assert 'contam1' in spectrum_plot.renderers[0].data_source.data
    assert 'contam2' not in spectrum_plot.renderers[0].data_source.data
    assert not any(key.startswith('contam2_')
                   for key in spectrum_plot.renderers[0].data_source.data)


def test_soss_layout_places_legacy_plot_above_slider(monkeypatch):
    """SOSS results retain the legacy view above the slider summary."""

    legacy_plot = contamination_figure.Spacer(width=1, height=1)
    slider_plot = contamination_figure.Spacer(width=2, height=2)
    captured = {}

    def fake_legacy(cube, instrument, targetName, badPAs):
        captured['cube'] = cube
        captured['instrument'] = instrument
        captured['target_name'] = targetName
        captured['bad_pas'] = badPAs
        return legacy_plot

    def fake_slider(pctlines, badPA_list, instrument):
        captured['pctlines'] = pctlines
        captured['slider_bad_pas'] = badPA_list
        captured['slider_instrument'] = instrument
        return slider_plot

    monkeypatch.setattr(contamination_figure, 'contam', fake_legacy)
    monkeypatch.setattr(contamination_figure, 'contam_slider_plot', fake_slider)
    targframes = [np.arange(6).reshape(2, 3), np.arange(6, 12).reshape(2, 3)]
    starcube = np.arange(24).reshape(4, 2, 3)
    pctlines = [np.zeros((4, 3)), np.ones((4, 3))]

    layout = contamination_figure.soss_contamination_plot_layout(
        targframes, starcube, pctlines, [7], 'NIS_SUBSTRIP96', 'Target')

    assert layout.children == [legacy_plot, slider_plot]
    assert captured['instrument'] == 'NIS_SUBSTRIP96'
    assert captured['slider_instrument'] == 'NIS_SUBSTRIP96'
    assert captured['target_name'] == 'Target'
    assert captured['bad_pas'] == [7]
    assert captured['slider_bad_pas'] == [7]
    np.testing.assert_array_equal(
        captured['cube'][0], targframes[0].T[::-1, ::-1])
    np.testing.assert_array_equal(
        captured['cube'][1], targframes[1].T[::-1, ::-1])
    np.testing.assert_array_equal(
        captured['cube'][2:], starcube.swapaxes(1, 2)[:, ::-1, ::-1])


def test_dhs_modes_available_in_web_form():
    """The web form exposes both supported NIRCam DHS filters."""

    choices = dict(CONTAM_VISIBILITY_MODES)
    assert choices['NRCA5_41STRIPE1_DHS_F322W2'] == 'NIRCam - DHS - F322W2'
    assert choices['NRCA5_41STRIPE1_DHS_F444W'] == 'NIRCam - DHS - F444W'


def _synthetic_miri_frames(contaminant_scale=0.):
    mask = np.zeros(miri_lrs.SHAPE)
    mask[:, 28:40] = 1
    target = mask.copy()
    contaminant = mask * contaminant_scale
    return target, contaminant, mask


def test_miri_isolated_and_overlapping_contamination():
    """MIRI preserves ExoCTK's pixel-level fraction-then-mean statistic."""

    target, isolated, mask = _synthetic_miri_frames(0.)
    result = miri_lrs.contamination_fraction(target, isolated, mask)
    assert np.allclose(result, 0.)

    target, overlapping, mask = _synthetic_miri_frames(1.)
    result = miri_lrs.contamination_fraction(target, overlapping, mask)
    assert np.allclose(result, 0.5)


def test_miri_contamination_empty_channels_do_not_warn():
    """Expected empty LRS wavelength channels remain NaN without warnings."""

    target = np.array([[1., 0.], [1., 0.]])
    contaminants = np.zeros_like(target)
    mask = np.ones_like(target)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result = miri_lrs.contamination_fraction(
            target, contaminants, mask, collapse_axis=0)

    assert result[0] == pytest.approx(0.)
    assert np.isnan(result[1])
    assert not any('Mean of empty slice' in str(warning.message)
                   for warning in caught)


def test_miri_fluxscale_changes_result_predictably():
    target, contaminant, mask = _synthetic_miri_frames(3.)
    result = miri_lrs.contamination_fraction(target, contaminant, mask)
    assert np.allclose(result, 0.75)


def test_miri_contamination_ignores_flux_outside_extraction():
    """An off-aperture trace must not dilute the in-aperture mean."""

    target, contaminant, mask = _synthetic_miri_frames(1.)
    expected = miri_lrs.contamination_fraction(target, contaminant, mask)
    contaminant[:, :10] = 100.

    result = miri_lrs.contamination_fraction(target, contaminant, mask)

    assert np.array_equal(result, expected)


def test_miri_all_pa_reducer_ignores_flux_outside_extraction():
    """The slider reducer must preserve the same MIRI mask semantics."""

    target, contaminant, mask = _synthetic_miri_frames(1.)
    contaminant_cube = np.stack([contaminant, contaminant.copy()])
    expected = field_simulator.fraction_contaminated(
        miri_lrs.APERTURE, [target], contaminant_cube,
        trace_masks=[mask], collapse_axis=miri_lrs.CROSS_DISPERSION_AXIS)[0]
    contaminant_cube[1, :, :10] = 100.

    result = field_simulator.fraction_contaminated(
        miri_lrs.APERTURE, [target], contaminant_cube,
        trace_masks=[mask], collapse_axis=miri_lrs.CROSS_DISPERSION_AXIS)[0]

    assert np.array_equal(result[0], expected[0])
    assert np.array_equal(result[1], expected[1])


def test_miri_ignores_contaminant_with_masked_fluxscale(monkeypatch):
    """A masked flux scale must never be interpreted as unity."""

    trace = np.ones(miri_lrs.SHAPE)
    asset = miri_lrs.MiriTraceAsset(
        trace=trace,
        wavelength=np.linspace(13.8, 5., miri_lrs.SHAPE[0]),
        extraction_mask=trace.copy(),
        valid_wavelength=np.ones(miri_lrs.SHAPE[0], dtype=bool),
        reference_pixel=(290, 34),
        metadata={},
    )
    monkeypatch.setattr(miri_lrs, 'load_trace', lambda teff: asset)
    stars = Table({
        'name': ['target', 'unknown-flux source'],
        'ra': [10., 10.],
        'dec': [20., 20.],
        'Teff': [5000., 4000.],
        'type': ['STAR', 'STAR'],
        **{name: np.zeros(2) for name in
           ('xdet', 'ydet', 'xtel', 'ytel', 'xsci', 'ysci')},
    })
    stars['fluxscale'] = MaskedColumn(
        [1., 1.], mask=[False, True])
    aperture = pysiaf.Siaf('MIRI')[miri_lrs.APERTURE]

    result = miri_lrs.calc_v3pa(0., stars, aperture)

    assert np.allclose(result['contaminants'], 0.)
    assert np.allclose(result['contamination'], 0.)


@pytest.mark.parametrize(
    'shift, expected',
    [((2, 3), (2, 3)), ((-2, -3), (0, 0))],
)
def test_miri_trace_shift_sign_and_clipping(shift, expected):
    trace = np.zeros((4, 4))
    trace[0, 0] = 1
    shifted = miri_lrs.shift_trace(trace, *shift, shape=(4, 4))
    if shift[0] < 0:
        assert not shifted.any()
    else:
        assert shifted[expected] == 1


def test_miri_expanded_source_rate_moves_onto_detector_after_shift():
    trace = np.zeros(miri_lrs.SHAPE)
    trace[-1, 34] = 10.
    trace[0, 34] = 20.
    wavelength = np.linspace(5.02, 13.863, miri_lrs.SHAPE[0])
    red_extension = np.zeros((2, miri_lrs.SHAPE[1]))
    red_extension[0, 34] = 5.
    blue_extension = np.zeros((2, miri_lrs.SHAPE[1]))
    blue_extension[-1, 34] = 4.
    asset = miri_lrs.MiriTraceAsset(
        trace=trace,
        wavelength=wavelength,
        extraction_mask=np.ones(miri_lrs.SHAPE),
        valid_wavelength=np.isfinite(wavelength),
        reference_pixel=(290, 34),
        metadata={},
        red_extension=red_extension,
        red_extension_y0=miri_lrs.SHAPE[0],
        blue_extension=blue_extension,
        blue_extension_y0=-2,
    )

    assert not np.any(miri_lrs.render_trace(asset, 0, 0) - trace)
    upstream = miri_lrs.render_trace(asset, -1, 0)
    assert upstream[-2, 34] == trace[-1, 34]
    assert upstream[-1, 34] == red_extension[0, 34]
    downstream = miri_lrs.render_trace(asset, 1, 0)
    assert downstream[0, 34] == blue_extension[-1, 34]
    assert downstream[1, 34] == trace[0, 34]


def test_miri_axes_mask_and_wavelength_contract():
    target, _, mask = _synthetic_miri_frames()
    wavelength = np.full(miri_lrs.SHAPE[0], np.nan)
    wavelength[6:-6] = np.linspace(13.8, 5.0, miri_lrs.SHAPE[0] - 12)
    valid = np.isfinite(wavelength)
    assert target.shape == mask.shape == (384, 68)
    assert miri_lrs.DISPERSION_AXIS == 0
    assert miri_lrs.CROSS_DISPERSION_AXIS == 1
    assert valid.sum() == 372
    assert np.all(np.diff(wavelength[valid]) < 0)
    assert np.all(target[mask.astype(bool)] > 0)

    # Native MIRI dispersion is detector Y, so all-PA reduction must collapse
    # detector X and retain 384 wavelength rows.
    cube_result = field_simulator.fraction_contaminated(
        miri_lrs.APERTURE, [target], np.zeros_like(target),
        trace_masks=[mask], collapse_axis=miri_lrs.CROSS_DISPERSION_AXIS)
    assert cube_result[0].shape == (1, miri_lrs.SHAPE[0])


def test_miri_result_pickle_round_trip():
    target, contaminants, mask = _synthetic_miri_frames(0.25)
    result = {
        'pa': 42., 'target': target, 'contaminants': contaminants,
        'extraction_mask': mask, 'wavelength': np.linspace(13., 5., 384),
    }
    restored = pickle.loads(pickle.dumps(result))
    assert restored['pa'] == 42.
    for key in ('target', 'contaminants', 'extraction_mask', 'wavelength'):
        assert np.array_equal(restored[key], result[key])
    pa_status = miri_lrs.PositionAngleResults([10, 11], [0, 1])
    restored_status = pickle.loads(pickle.dumps(pa_status))
    assert restored_status == [10, 11]
    assert restored_status.inaccessible == [0, 1]


def test_miri_precompute_can_add_target_to_existing_database(
        tmp_path, monkeypatch):
    """Living caches retain the LRS wavelength products."""

    trace = np.ones(miri_lrs.SHAPE)
    wavelength = np.linspace(5., 14., miri_lrs.SHAPE[0])
    mask = np.ones(miri_lrs.SHAPE)
    asset = miri_lrs.MiriTraceAsset(
        trace=trace, wavelength=wavelength, extraction_mask=mask,
        valid_wavelength=np.ones(miri_lrs.SHAPE[0], dtype=bool),
        reference_pixel=(0, 0), metadata={'source': 'synthetic'})
    monkeypatch.setattr(precompute.miri_lrs, 'load_trace', lambda teff: asset)

    filename = tmp_path / 'miri-cache.h5'
    contamination = np.zeros((360, *miri_lrs.SHAPE))
    contamination[240] = 0.5
    pa_results = miri_lrs.PositionAngleResults(
        successful=[240], inaccessible=[0, 1])
    precompute.save_exoplanet_data(
        filename, 'HD 189733 b', miri_lrs.APERTURE, 300.182, 22.711,
        trace[None], contamination, goodPA_list=pa_results)

    with precompute.h5py.File(filename, 'r') as handle:
        group = handle['HD 189733 b']
        assert group.attrs['filled']
        assert list(group['plane_index'][:]) == [240]
        assert list(group.attrs['inaccessiblePA_list']) == [0, 1]
        np.testing.assert_array_equal(group['wavelength'][:], wavelength)
        np.testing.assert_array_equal(group['extraction_mask'][:], mask)


def test_miri_pa_rotation_uses_siaf_coordinates():
    aperture = pysiaf.Siaf('MIRI')[miri_lrs.APERTURE]
    xdet, ydet = aperture.reference_point('det')
    xtel, ytel = aperture.det_to_tel(xdet, ydet)
    ra, dec = 10., 20.
    attitude0 = pysiaf.utils.rotations.attitude_matrix(
        xtel, ytel, ra, dec, 0.)
    source_ra, source_dec = pysiaf.utils.rotations.tel_to_sky(
        attitude0, xtel + 1., ytel)
    stars = Table({
        'ra': [ra, source_ra.to_value(u.deg)],
        'dec': [dec, source_dec.to_value(u.deg)],
        **{name: np.zeros(2) for name in
           ('xdet', 'ydet', 'xtel', 'ytel', 'xsci', 'ysci')},
    })
    offset0 = miri_lrs.source_offsets(stars.copy(), aperture, attitude0)[1]
    attitude90 = pysiaf.utils.rotations.attitude_matrix(
        xtel, ytel, ra, dec, 90.)
    offset90 = miri_lrs.source_offsets(stars.copy(), aperture, attitude90)[1]
    assert offset0 != offset90
    assert np.hypot(*offset0) == pytest.approx(np.hypot(*offset90), abs=1)


def test_hd189733b_moves_toward_lower_detector_rows_at_v3pa_238():
    """The upstream companion's +Y_SCI offset maps to a lower array row."""

    stars = Table({
        'ra': [300.182137, 300.1789791392316],
        'dec': [22.710853, 22.707659537029453],
        **{name: np.zeros(2) for name in
           ('xdet', 'ydet', 'xtel', 'ytel', 'xsci', 'ysci')},
    })
    aperture = pysiaf.Siaf('MIRI')[miri_lrs.APERTURE]
    xdet, ydet = aperture.reference_point('det')
    xtel, ytel = aperture.det_to_tel(xdet, ydet)
    attitude = pysiaf.utils.rotations.attitude_matrix(
        xtel, ytel, stars['ra'][0], stars['dec'][0], 238.)

    offsets = miri_lrs.source_offsets(stars, aperture, attitude)

    assert offsets == [(0, 0), (-132, 49)]


def test_miri_v3pa_has_no_additional_aperture_angle(monkeypatch):
    """MIRI source placement passes V3PA directly to the SIAF attitude."""

    aperture = pysiaf.Siaf('MIRI')[miri_lrs.APERTURE]
    xdet, ydet = aperture.reference_point('det')
    xtel, ytel = aperture.det_to_tel(xdet, ydet)
    ra, dec, v3pa = 10., 20., 123.4
    attitude = pysiaf.utils.rotations.attitude_matrix(
        xtel, ytel, ra, dec, v3pa)
    source_ra, source_dec = pysiaf.utils.rotations.tel_to_sky(
        attitude, xtel + 5., ytel - 2.)
    stars = Table({
        'name': ['target', 'offset source'],
        'ra': [ra, source_ra.to_value(u.deg)],
        'dec': [dec, source_dec.to_value(u.deg)],
        'Teff': [5000., 4000.],
        'fluxscale': [1., 0.1],
        'type': ['STAR', 'STAR'],
        **{name: np.zeros(2) for name in
           ('xdet', 'ydet', 'xtel', 'ytel', 'xsci', 'ysci')},
    })
    expected = miri_lrs.source_offsets(
        stars.copy(), aperture, attitude)[1]
    trace = np.ones(miri_lrs.SHAPE)
    asset = miri_lrs.MiriTraceAsset(
        trace=trace,
        wavelength=np.linspace(13.8, 5., miri_lrs.SHAPE[0]),
        extraction_mask=trace.copy(),
        valid_wavelength=np.ones(miri_lrs.SHAPE[0], dtype=bool),
        reference_pixel=(290, 34),
        metadata={},
    )
    monkeypatch.setattr(miri_lrs, 'load_trace', lambda teff: asset)

    result = miri_lrs.calc_v3pa(v3pa, stars, aperture)
    source = result['intersecting_sources'][0]

    assert (source['y_shift'], source['x_shift']) == expected


def test_miri_nearby_source_pruning_is_conservative():
    stars = Table({
        'ra': [10., 10.005, 10.1],
        'dec': [20., 20., 20.],
    })
    indices = miri_lrs.nearby_source_indices(stars, 60.)
    assert np.array_equal(indices, [0, 1])


def test_miri_source_query_covers_pruning_radius():
    width = field_simulator.source_query_width(miri_lrs.APERTURE)
    enclosing_radius = np.hypot(width, width) / 2.
    assert enclosing_radius.to_value(u.arcsec) == pytest.approx(60.)
    assert field_simulator.source_query_width(
        'NIS_SUBSTRIP256') == 7.5 * u.arcmin


def test_miri_trace_assets_are_cached(tmp_path):
    trace_dir = (tmp_path / 'exoctk_contam' / 'traces' /
                 miri_lrs.APERTURE)
    trace_dir.mkdir(parents=True)
    reference_teff = 4000
    path = trace_dir / f'{miri_lrs.APERTURE}_{reference_teff}.fits'
    header = fits.Header({
        'APERTURE': miri_lrs.APERTURE,
        'REFYPIX': 290,
        'REFXPIX': 34,
        'METAJSON': json.dumps({
            'blue_extension_y0': -2,
            'red_extension_y0': miri_lrs.SHAPE[0],
        }),
    })
    fits.HDUList([
        fits.PrimaryHDU(header=header),
        fits.ImageHDU(np.ones(miri_lrs.SHAPE), name='TRACE'),
        fits.ImageHDU(np.arange(miri_lrs.SHAPE[0]), name='WAVELENGTH'),
        fits.ImageHDU(np.ones(miri_lrs.SHAPE), name='EXTRACTION_MASK'),
        fits.ImageHDU(np.ones(miri_lrs.SHAPE[0]), name='VALID_WAVELENGTH'),
        fits.ImageHDU(np.ones((2, miri_lrs.SHAPE[1])),
                      name='BLUE_EXTENSION'),
        fits.ImageHDU(np.ones((6, miri_lrs.SHAPE[1])),
                      name='RED_EXTENSION'),
    ]).writeto(path)

    miri_lrs.load_trace.cache_clear()
    first = miri_lrs.load_trace(reference_teff, str(tmp_path))
    path.unlink()
    second = miri_lrs.load_trace(reference_teff, str(tmp_path))
    assert second is first
    assert miri_lrs.load_trace.cache_info().hits == 1
    miri_lrs.load_trace.cache_clear()


def test_miri_reference_trace_uses_named_calibration_asset(monkeypatch):
    """Shared calibration lookups use the documented reference asset."""

    loaded = []
    sentinel = object()

    def load_trace(teff, exoctk_data=None):
        loaded.append((teff, exoctk_data))
        return sentinel

    monkeypatch.setattr(miri_lrs, 'load_trace', load_trace)

    assert miri_lrs.load_reference_trace('/test/data') is sentinel
    assert loaded == [(4000, '/test/data')]


def test_miri_product_retains_detector_rows_outside_wavelength_grid():
    detector = np.zeros((427, 55))
    detector[0, :] = 7.
    detector[21, :] = 6.
    detector[22, :] = 1.
    detector[27, :] = 2.
    detector[28, :] = 3.
    detector[405, :] = 4.
    detector[406, :] = 5.
    report = {
        '1d': {'wave_pix': np.linspace(13.86, 5.02, 372)},
        'scalar': {'extraction_area': 12.},
    }

    (trace, wavelength, valid, mask, blue_extension, red_extension,
     registration) = make_miri_lrs_traces._miri_product(
        report, detector, 1.32)

    assert np.all(trace[0, 7:62] == 1.)
    assert np.all(trace[5, 7:62] == 2.)
    assert np.all(trace[6, 7:62] == 3.)
    assert np.all(trace[-1, 7:62] == 4.)
    assert np.all(blue_extension[0, 7:62] == 7.)
    assert np.all(blue_extension[-1, 7:62] == 6.)
    assert np.all(red_extension[0, 7:62] == 5.)
    assert registration['pandeia_detector_row0'] == 22
    assert registration['pandeia_detector_product'] == 'noiseless_source_rate'
    assert registration['blue_extension_y0'] == -22
    assert registration['red_extension_y0'] == miri_lrs.SHAPE[0]
    assert wavelength[6] == pytest.approx(5.02)
    assert wavelength[377] == pytest.approx(13.86)
    assert np.all(np.diff(wavelength[valid]) > 0.)
    assert np.count_nonzero(valid) == 372
    assert np.count_nonzero(mask) == 372 * 12


def test_miri_product_registers_extended_projection_to_native_grid():
    detector = np.zeros((439, 55))
    detector[22, :] = 1.
    detector[405, :] = 2.
    detector[406, :] = 3.
    report = {
        '1d': {'wave_pix': np.linspace(14.015, 5.02, 384)},
        'scalar': {'extraction_area': 12.},
    }
    calibrated = np.linspace(13.86, 5.02, 372)

    (trace, wavelength, valid, mask, blue_extension, red_extension,
     registration) = make_miri_lrs_traces._miri_product(
        report, detector, 1.32, calibrated_wave=calibrated)

    assert np.all(trace[0, 7:62] == 1.)
    assert np.all(trace[-1, 7:62] == 2.)
    assert np.all(red_extension[0, 7:62] == 3.)
    assert registration['pandeia_detector_row0'] == 22
    assert registration['pandeia_projected_wavelength_rows'] == 384
    assert registration['pandeia_calibrated_wavelength_rows'] == 372
    assert wavelength[6] == pytest.approx(5.02)
    assert wavelength[377] == pytest.approx(13.86)
    assert np.count_nonzero(valid) == 372
    assert np.count_nonzero(mask) == 372 * 12


def test_miri_extended_wave_grid_reaches_reference_response_zero():
    native = np.linspace(5., 13.86, 380)
    extended = make_miri_lrs_traces._extended_miri_wave_pix(
        native, 14.01584)

    assert np.array_equal(extended[:len(native)], native)
    assert extended[-1] == pytest.approx(14.01584)
    assert np.all(np.diff(extended) > 0.)
    assert np.max(np.diff(extended)[len(native) - 1:-1]) < 0.03
    assert len(extended) > len(native)


def test_miri_templates_use_gaia_g_normalization():
    """Template normalization must match field_simulator's Gaia flux ratios."""

    normalization = make_miri_lrs_traces._scene(4000, 10.)[
        'spectrum']['normalization']

    assert normalization == {
        'type': 'photsys',
        'bandpass': 'gaia,g',
        'norm_flux': 10.,
        'norm_fluxunit': 'vegamag',
    }


def test_miri_observable_pa_ranges_use_v3pa_not_instrument_angle():
    """The all-PA renderer must not apply the SIAF angle a second time."""

    v3pas = np.arange(87., 107.)
    table = DataFrame({
        'V3PA_min_pa_angle': v3pas,
        'V3PA_nominal_angle': v3pas,
        'V3PA_max_pa_angle': v3pas,
        'MIRI_min_pa_angle': v3pas + 4.83544897,
        'MIRI_nominal_angle': v3pas + 4.83544897,
        'MIRI_max_pa_angle': v3pas + 4.83544897,
    })
    position_angles, bounds, sampled = (
        field_simulator.observable_v3pa_ranges(table))
    assert bounds == [(87, 106)]
    assert np.array_equal(position_angles, np.arange(87, 107))
    assert np.array_equal(sampled, np.arange(87, 107))


def test_miri_calculates_every_pa_but_other_modes_do_not():
    observable = np.arange(87, 107)

    assert np.array_equal(
        field_simulator.calculation_v3pas(miri_lrs.APERTURE, observable),
        np.arange(360))
    assert np.array_equal(
        field_simulator.calculation_v3pas('NIS_SUBSTRIP256', observable),
        observable)


def test_miri_observability_is_independent_of_calculated_pas():
    results = miri_lrs.PositionAngleResults(
        successful=np.arange(360), inaccessible=[0, 1, 359])

    assert field_simulator.unobservable_v3pas(results) == [0, 1, 359]


def test_miri_apt_v3pa_has_extraction_overlapping_source(monkeypatch):
    """APT V3PA=53.82 places a nearby WASP-12 Gaia trace in the mask."""

    trace = np.zeros(miri_lrs.SHAPE)
    trace[:, 28:40] = 1.
    wavelength = np.linspace(13.8, 5., miri_lrs.SHAPE[0])
    mask = trace.copy()
    asset = miri_lrs.MiriTraceAsset(
        trace=trace, wavelength=wavelength, extraction_mask=mask,
        valid_wavelength=np.ones(miri_lrs.SHAPE[0], dtype=bool),
        reference_pixel=(290, 34), metadata={})
    loaded_temperatures = []

    def load_trace(teff):
        loaded_temperatures.append(float(teff))
        return asset

    monkeypatch.setattr(miri_lrs, 'load_trace', load_trace)

    # WASP-12 and VizieR/Gaia DR3 source 3435282866759615488.
    stars = Table({
        'name': ['WASP-12', '3435282866759615488', 'far source'],
        'ra': [97.63664477033, 97.63945094102, 98.],
        'dec': [29.67226537027, 29.67372250549, 30.],
        'Teff': [6000., 4000., 5000.],
        'fluxscale': [1., 0.1, 0.1],
        'type': ['STAR', 'STAR', 'STAR'],
        **{name: np.zeros(3) for name in
           ('xdet', 'ydet', 'xtel', 'ytel', 'xsci', 'ysci')},
    })
    aperture = pysiaf.Siaf('MIRI')[miri_lrs.APERTURE]
    result = miri_lrs.calc_v3pa(53.82, stars, aperture)

    assert result['intersecting_sources'][0]['name'] == (
        '3435282866759615488')
    assert result['intersecting_sources'][0]['x_shift'] == -1
    assert np.nanmax(result['contamination']) > 0
    assert loaded_temperatures == [6000., 4000.]


def test_miri_mode_available_in_web_form():
    choices = dict(CONTAM_VISIBILITY_MODES)
    assert choices['MIRIM_SLITLESSPRISM_IP'] == (
        'MIRI - LRS - SLITLESSPRISM_IP')
    assert 'MIRIM_SLITLESSPRISM' not in choices


def test_miri_all_pa_plot_accepts_native_shape_and_wavelength():
    fractions = [np.zeros((360, miri_lrs.SHAPE[0]))]
    fractions[0][175, 100] = 0.05
    wavelength = np.linspace(13.8, 5.0, miri_lrs.SHAPE[0])
    plot = contamination_figure.contam_slider_plot(
        fractions, [], wavelength=wavelength, trace_names=['MIRI LRS'],
        y_max=0.1)
    assert plot is not None
    assert plot.children[0].y_range.end == pytest.approx(10.)
    assert plot.children[-1].y_range.end == pytest.approx(10.)
    assert plot.children[0].yaxis.axis_label == 'Contamination (%)'
    assert plot.children[-1].yaxis.axis_label == 'Mean Contamination (%)'
    plotted = plot.children[0].renderers[0].data_source.data['contam1']
    assert np.nanmax(plotted) == pytest.approx(5.)
    threshold_legend = plot.children[-1].below[-1]
    labels = [item.label.value for item in threshold_legend.items]
    assert 'Target not observable' in labels
    assert 'Spectrum > 5.0% Contaminated' in labels


def test_miri_single_pa_plot_accepts_runtime_result():
    """The plotting module accepts the MIRI runtime result contract."""

    valid = np.ones(miri_lrs.SHAPE[0], dtype=bool)
    result = {
        'pa': 42.,
        'target': np.ones(miri_lrs.SHAPE),
        'contaminants': np.zeros(miri_lrs.SHAPE),
        'extraction_mask': np.ones(miri_lrs.SHAPE),
        'valid_wavelength': valid,
        'wavelength': np.linspace(5., 13.8, miri_lrs.SHAPE[0]),
        'contamination': np.zeros(miri_lrs.SHAPE[0]),
        'intersecting_sources': [],
    }

    plot = contamination_figure.miri_single_pa_plot(result)

    assert len(plot.children) == 2

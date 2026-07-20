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

import os
import sys
from pathlib import Path

import numpy as np
import pytest

from pandas import DataFrame
from astropy.table import Table

from exoctk.contam_visibility import field_simulator
from exoctk.contam_visibility import contamination_figure
from exoctk.contam_visibility import resolve
from exoctk.contam_visibility import new_vis_plot
from exoctk.contam_visibility import contamination_figure
from exoctk.contam_visibility.modes import CONTAM_VISIBILITY_MODES

# Determine if tests are being run on Github Actions
ON_GITHUB_ACTIONS = '/home/runner' in os.path.expanduser('~') or '/Users/runner' in os.path.expanduser('~')


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


def test_resolve_target():
    """Tests the ``resolve_target`` function in the ``resolve`` module"""

    ra, dec = resolve.resolve_target('Wasp-18 b')

    assert ra == 24.3544618
    assert dec == -45.6777937


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


@pytest.mark.parametrize('aperture', [
    'NIS_SUBSTRIP96',
    'NIS_SUBSTRIP256',
    'NRCA5_41STRIPE1_DHS_F322W2',
    'NRCA5_41STRIPE1_DHS_F444W',
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

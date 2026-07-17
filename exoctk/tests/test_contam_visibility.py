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
import warnings
from pathlib import Path

import numpy as np
import pytest

from pandas import DataFrame
from astropy.table import MaskedColumn, Table

from exoctk.contam_visibility import field_simulator
from exoctk.contam_visibility import resolve
from exoctk.contam_visibility import new_vis_plot
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


def test_dhs_modes_available_in_web_form():
    """The web form exposes both supported NIRCam DHS filters."""

    choices = dict(CONTAM_VISIBILITY_MODES)
    assert choices['NRCA5_41STRIPE1_DHS_F322W2'] == 'NIRCam - DHS - F322W2'
    assert choices['NRCA5_41STRIPE1_DHS_F444W'] == 'NIRCam - DHS - F444W'

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
import requests

from exoctk.utils import DATA_URLS, color_gen, download_exoctk_data, filter_table, get_canonical_name, get_target_data, medfilt

ARGS = ['planet_name', 'planet_data']

# Include 3 actual planets with data and one fake planet to xfail.
TEST_DATA = [('HD108236f', {'canonical_name': 'HD 108236 f', 'Fe/H': -0.28, 'Teff': 5660.0, 'stellar_gravity': 4.49, 'transit_duration': 0.13625, 'RA': 186.5740627, 'DEC': -51.3630519}),
             ('NGTS-14Ab', {'canonical_name': 'NGTS-14 A b', 'Fe/H': 0.1, 'Teff': 5187.0, 'stellar_gravity': 4.2, 'transit_duration': 0.09333333333333334, 'RA': 328.5174799, 'DEC': -38.3774193}),
             ('2MASSJ10193800-0948225b', {'canonical_name': 'WASP-43 b', 'Fe/H': -0.05, 'Teff': 4400.0, 'stellar_gravity': 4.49, 'transit_duration': 0.0513, 'RA': 154.9081869, 'DEC': -9.8064431}),
             pytest.param('djgfjhsg', {'canonical_name': 'sfghsfkjg', 'Fe/H': -999, 'Teff': -999, 'stellar_gravity': -999, 'transit_duration': -999, 'RA': -999, 'DEC': -999}, marks=pytest.mark.xfail)]

MEDFILT_DATA = [(np.array([1, 1, 8, 12, 2, 10, 5, 2, 5, 2]), 4),
                (np.array([1, 1, 8, 12, 2, 10, 5, 2, 5, 2]), 4)]


@pytest.mark.parametrize(ARGS, TEST_DATA)
def test_get_canoical_name(planet_name, planet_data):
    """Test that the canonical name of the planet is returned by exomast"""

    canonical_name = get_canonical_name(planet_name)
    assert canonical_name == planet_data['canonical_name']


@pytest.mark.parametrize(ARGS, TEST_DATA)
def test_get_target_data(planet_name, planet_data):
    """Test that the canonical name of the planet is returned by exomast"""

    data, _ = get_target_data(planet_name)

    # these are some params that the webapp uses for it's tools
    exomast_params = ['Fe/H', 'Teff', 'stellar_gravity', 'transit_duration', 'RA', 'DEC']

    for value in exomast_params:
        assert data[value] == pytest.approx(planet_data[value], 0.1)


@pytest.fixture
def build_table():
    """ Fixture to build table for tests requiring table """
    column_names = ['wavelength', 'flux']
    wavelength = np.arange(3000, 5000)
    flux = np.random.rand(wavelength.shape[0]) * 10e-13

    table = Table(data=[wavelength, flux], names=column_names)

    return table


@pytest.mark.parametrize("operator", ['>4021', '<3856', '>=4928', '<=3740', '==4000'])
def test_filter_table(build_table, operator):
    """test table filter function with fake table"""
    # Test wavelength sort
    TEST_DATA = filter_table(build_table, wavelength=operator)
    assert all(TEST_DATA['wavelength'] == eval("build_table[np.where(build_table['wavelength'] {})]['wavelength']".format(operator)))


@pytest.mark.parametrize("colormap", ['viridis', pytest.param('hjdsgfdhsf', marks=pytest.mark.xfail)])
def test_color_gen(colormap):
    color_gen(colormap)


@pytest.mark.parametrize(['data', 'filter_window'], MEDFILT_DATA)
def test_medfilt(data, filter_window):
    medfilt(data, filter_window)


def test_download_exoctk_data_uses_selected_package_and_log(tmp_path, monkeypatch):
    """Only the requested package and log archive should be downloaded."""

    requests_seen = []
    archives_seen = []

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def raise_for_status(self):
            return None

        def iter_content(self, chunk_size):
            yield b'archive data'

    def get(url, **kwargs):
        requests_seen.append((url, kwargs))
        return Response()

    def unpack_archive(path, destination):
        archives_seen.append((path, destination))

    monkeypatch.setattr(requests, 'get', get)
    monkeypatch.setattr('exoctk.utils.shutil.unpack_archive', unpack_archive)

    contam_urls = list(DATA_URLS['exoctk_contam'])
    log_urls = list(DATA_URLS['exoctk_log'])
    download_exoctk_data('exoctk_contam', str(tmp_path))

    assert [url for url, _ in requests_seen] == contam_urls + log_urls
    assert all(kwargs == {'stream': True, 'timeout': 60} for _, kwargs in requests_seen)
    assert len(archives_seen) == 2
    assert all(destination == str(tmp_path) for _, destination in archives_seen)
    assert DATA_URLS['exoctk_contam'] == contam_urls


def test_download_exoctk_data_reports_http_failure(tmp_path, monkeypatch):
    """An HTTP error should not be passed to the archive extractor."""

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def raise_for_status(self):
            raise requests.HTTPError('404 Not Found')

        def iter_content(self, chunk_size):
            yield b'not reached'

    monkeypatch.setattr('exoctk.utils.requests.get', lambda *args, **kwargs: Response())

    with pytest.raises(RuntimeError, match='Unable to download ExoCTK data'):
        download_exoctk_data('exoctk_contam', str(tmp_path))

    assert list(tmp_path.iterdir()) == []

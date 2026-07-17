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

import os

from astropy.table import Table
import numpy as np
import pytest
import requests

from exoctk import utils
from exoctk.utils import (DATA_URLS, ExoMASTError, color_gen,
                          download_exoctk_data, filter_table,
                          get_canonical_name, get_target_data, medfilt)

TEST_TARGET = {
    'canonical_name': 'HD 108236 f',
    'Fe/H': -0.28,
    'Teff': 5660.0,
    'stellar_gravity': 4.49,
    'transit_duration': 0.13625,
    'RA': 186.5740627,
    'DEC': -51.3630519,
}

MEDFILT_DATA = [(np.array([1, 1, 8, 12, 2, 10, 5, 2, 5, 2]), 4),
                (np.array([1, 1, 8, 12, 2, 10, 5, 2, 5, 2]), 4)]


class MockResponse:
    """Minimal requests response used to test ExoMAST handling."""

    def __init__(self, payload, error=None, json_error=None):
        self.payload = payload
        self.error = error
        self.json_error = json_error

    def raise_for_status(self):
        if self.error is not None:
            raise self.error

    def json(self):
        if self.json_error is not None:
            raise self.json_error
        return self.payload


def test_get_canonical_name_uses_valid_exomast_response(monkeypatch):
    """Canonical-name lookup should parse a successful ExoMAST response."""

    calls = []

    def mock_get(url, **kwargs):
        calls.append((url, kwargs))
        return MockResponse({'canonicalName': TEST_TARGET['canonical_name']})

    monkeypatch.setattr(utils.requests, 'get', mock_get)
    assert get_canonical_name('HD108236f') == TEST_TARGET['canonical_name']
    assert calls[0][1]['params'] == {'name': 'HD108236f'}
    assert calls[0][1]['timeout'] == utils.EXOMAST_TIMEOUT_SECONDS


def test_get_target_data_prefers_nexsci_catalog(monkeypatch):
    """Target lookup should use the preferred catalog from valid responses."""

    def mock_get(url, **kwargs):
        if 'identifiers' in url:
            return MockResponse(
                {'canonicalName': TEST_TARGET['canonical_name']})
        alternate = dict(TEST_TARGET, catalog_name='other', Teff=5000.)
        preferred = dict(TEST_TARGET, catalog_name='nexsci')
        return MockResponse([alternate, preferred])

    monkeypatch.setattr(utils.requests, 'get', mock_get)
    data, url = get_target_data('HD108236f')
    assert data['catalog_name'] == 'nexsci'
    assert data['Teff'] == pytest.approx(TEST_TARGET['Teff'])
    assert TEST_TARGET['canonical_name'].replace(' ', '') in url


@pytest.mark.parametrize('payload', [{}, None, []])
def test_get_canonical_name_rejects_malformed_response(monkeypatch, payload):
    monkeypatch.setattr(
        utils.requests, 'get', lambda *args, **kwargs: MockResponse(payload))
    with pytest.raises(ExoMASTError, match='no canonical name'):
        get_canonical_name('not-a-target')


def test_get_canonical_name_wraps_request_errors(monkeypatch):
    monkeypatch.setattr(
        utils.requests, 'get',
        lambda *args, **kwargs: MockResponse(
            None, error=requests.HTTPError('404 Client Error')))
    with pytest.raises(ExoMASTError, match='request failed'):
        get_canonical_name('not-a-target')


def test_get_canonical_name_wraps_invalid_json(monkeypatch):
    monkeypatch.setattr(
        utils.requests, 'get',
        lambda *args, **kwargs: MockResponse(None, json_error=ValueError()))
    with pytest.raises(ExoMASTError, match='invalid JSON'):
        get_canonical_name('not-a-target')


def test_get_target_data_rejects_empty_response(monkeypatch):
    def mock_get(url, **kwargs):
        if 'identifiers' in url:
            return MockResponse(
                {'canonicalName': TEST_TARGET['canonical_name']})
        return MockResponse([])

    monkeypatch.setattr(utils.requests, 'get', mock_get)
    with pytest.raises(ExoMASTError, match='no usable target data'):
        get_target_data('HD108236f')


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
    assert [os.path.basename(path) for path, _ in archives_seen] == [
        'exoctk_contam.tar.gz',
        'exoctk_log.tar.gz',
    ]
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

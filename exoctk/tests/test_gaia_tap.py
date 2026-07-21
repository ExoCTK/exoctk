"""Tests for endpoint-aware Gaia TAP querying."""

import numpy as np
import pytest

from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

from exoctk.contam_visibility.gaia_tap import (
    GAIA_TAP_ENDPOINTS,
    NORMALIZED_COLUMNS,
    GaiaFailoverTAP,
    GaiaTAPError,
)
from exoctk.contam_visibility.field_simulator import (
    calculate_current_coordinates,
)


CENTER = SkyCoord(10. * u.deg, 20. * u.deg)


def _native_table(endpoint, source_ids=(30, 10, 20), masked=False):
    """Return a small table using an endpoint's native column names."""

    count = len(source_ids)
    values = {
        'source_id': source_ids,
        'ra': [10.001, 10., 10.001][:count],
        'dec': np.full(count, 20.),
        'ref_epoch': np.full(count, 2016.),
        'pmra': np.arange(count, dtype=float),
        'pmdec': np.arange(count, dtype=float),
        'parallax': np.arange(count, dtype=float),
        'astrometric_excess_noise': np.zeros(count),
        'phot_g_mean_flux': np.full(count, 100.),
        'phot_g_mean_mag': np.full(count, 12.),
        'bp_rp': np.full(count, 1.),
        'phot_bp_rp_excess_factor': np.full(count, 1.),
        'classprob_dsc_combmod_star': np.asarray(source_ids) / 100.,
        'classprob_dsc_combmod_galaxy': np.asarray(source_ids) / 1000.,
        'classprob_dsc_combmod_quasar': np.asarray(source_ids) / 10000.,
    }
    internal_to_native = {
        internal: native.strip('"') for native, internal in endpoint.columns}
    columns = {}
    for internal, column in values.items():
        native = internal_to_native.get(internal)
        if native is not None:
            columns[native] = column
    table = Table(columns, masked=masked)
    if masked:
        table[internal_to_native['parallax']].mask[0] = True
    return table


@pytest.mark.parametrize(('index', 'fragments'), [
    (0, ('FROM gaiadr3.gaia_source', 'CONTAINS', 'ra AS ra')),
    (1, ('FROM gaia_dr3.gaia_source', 'q3c_radial_query', "'t'=")),
    (2, ('FROM "I/355/gaiadr3"', '"Source" AS source_id',
         '"RA_ICRS"', '"DE_ICRS"')),
])
def test_endpoint_native_query_construction(index, fragments):
    """Each service receives its own table, columns, and cone syntax."""

    query = GAIA_TAP_ENDPOINTS[index].build_query(10., 20., 0.1)

    assert 'SELECT *' not in query
    assert all(fragment in query for fragment in fragments)


@pytest.mark.parametrize('endpoint', GAIA_TAP_ENDPOINTS)
def test_endpoint_results_share_normalized_schema(endpoint):
    """Native endpoint tables normalize to one forward-compatible contract."""

    native = _native_table(endpoint)
    result = GaiaFailoverTAP._normalize(native, endpoint, CENTER)

    assert result.colnames == list(NORMALIZED_COLUMNS) + ['dist']
    assert all(name in result.colnames for name in (
        'classprob_dsc_combmod_star',
        'classprob_dsc_combmod_galaxy',
        'classprob_dsc_combmod_quasar'))


@pytest.mark.parametrize('endpoint', GAIA_TAP_ENDPOINTS)
def test_dsc_probabilities_remain_associated_with_source_ids(endpoint):
    """Renaming and deterministic sorting do not misalign DSC values."""

    result = GaiaFailoverTAP._normalize(
        _native_table(endpoint), endpoint, CENTER)
    probabilities = dict(zip(
        result['source_id'], result['classprob_dsc_combmod_star']))

    assert probabilities == {10: 0.1, 20: 0.2, 30: 0.3}


@pytest.mark.parametrize('endpoint', GAIA_TAP_ENDPOINTS)
def test_masked_values_remain_masked(endpoint):
    """Native nulls and unavailable optional fields stay masked."""

    result = GaiaFailoverTAP._normalize(
        _native_table(endpoint, masked=True), endpoint, CENTER)
    row = np.flatnonzero(result['source_id'] == 30)[0]

    assert np.ma.is_masked(result['parallax'][row])
    if endpoint.name == 'TAPVizieR':
        assert np.all(result['ref_epoch'].mask)


def test_unavailable_reference_epoch_skips_proper_motion_correction():
    """A masked optional epoch cannot break downstream source handling."""

    ra, dec = calculate_current_coordinates(
        10., 20., 1., 1., np.ma.masked, target_date=2026)

    assert np.ma.is_masked(ra)
    assert np.ma.is_masked(dec)


def test_rows_sort_by_distance_then_source_id():
    """Equal angular separations use source ID as a stable tie breaker."""

    endpoint = GAIA_TAP_ENDPOINTS[0]
    native = _native_table(endpoint, source_ids=(30, 10, 20))
    native['ra'] = [10.001, 10., 10.001]
    result = GaiaFailoverTAP._normalize(native, endpoint, CENTER)

    assert list(result['source_id']) == [10, 20, 30]


def test_missing_required_core_column_is_rejected():
    """A malformed result cannot pass adapter schema validation."""

    endpoint = GAIA_TAP_ENDPOINTS[0]
    native = _native_table(endpoint)
    native.remove_column('phot_g_mean_flux')

    with pytest.raises(GaiaTAPError, match='phot_g_mean_flux'):
        GaiaFailoverTAP._normalize(native, endpoint, CENTER)


@pytest.mark.parametrize('failure_count, expected_endpoint', [
    (1, 'NOIRLab'),
    (2, 'TAPVizieR'),
])
def test_query_region_fails_over_to_next_endpoint(
        monkeypatch, failure_count, expected_endpoint):
    """Transport/schema failures advance through endpoint priority order."""

    adapter = GaiaFailoverTAP(max_retries=1)
    attempts = []

    def run(endpoint, query):
        attempts.append(endpoint.name)
        if len(attempts) <= failure_count:
            raise GaiaTAPError('synthetic failure')
        return _native_table(endpoint)

    monkeypatch.setattr(adapter, '_run_query', run)
    result = adapter.query_region(CENTER, 1. * u.arcmin)

    assert attempts[-1] == expected_endpoint
    assert adapter.last_endpoint == next(
        item.url for item in GAIA_TAP_ENDPOINTS
        if item.name == expected_endpoint)
    assert len(result) == 3


def test_all_endpoint_failures_are_summarized(monkeypatch):
    """Total failure reports every service and retains exception chaining."""

    adapter = GaiaFailoverTAP(max_retries=1)

    def fail(endpoint, query):
        raise GaiaTAPError(f'{endpoint.name} unavailable')

    monkeypatch.setattr(adapter, '_run_query', fail)
    with pytest.raises(GaiaTAPError) as caught:
        adapter.query_region(CENTER, 1. * u.arcmin)

    assert all(endpoint.name in str(caught.value)
               for endpoint in GAIA_TAP_ENDPOINTS)
    assert caught.value.__cause__ is not None


class _Response:
    def __init__(self, text='', content=None, status=200, headers=None,
                 url='https://example.test/async'):
        self.text = text
        self.content = text.encode() if content is None else content
        self.status_code = status
        self.headers = headers or {}
        self.url = url
        self.ok = status < 400

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f'HTTP {self.status_code}')


class _Session:
    def __init__(self, responses=()):
        self.responses = iter(responses)
        self.calls = []

    def get(self, url, **kwargs):
        self.calls.append(('get', url, kwargs))
        return next(self.responses)

    def post(self, url, **kwargs):
        self.calls.append(('post', url, kwargs))
        return next(self.responses)


@pytest.mark.parametrize('phase', ('ERROR', 'ABORTED'))
def test_terminal_tap_failure_phases_are_rejected(phase):
    """Terminal unsuccessful UWS phases fail without fetching a result."""

    session = _Session((
        _Response(text=phase),
        _Response(text='service detail'),
    ))
    adapter = GaiaFailoverTAP(session=session, max_polls=1)

    with pytest.raises(GaiaTAPError, match=phase):
        adapter._poll_job('https://example.test/job/1')


def test_polling_timeout_is_finite():
    """A perpetually executing UWS job reaches the configured timeout."""

    session = _Session((_Response(text='EXECUTING'),) * 2)
    adapter = GaiaFailoverTAP(
        session=session, max_polls=2, poll_interval=0)

    with pytest.raises(TimeoutError, match='polling timeout'):
        adapter._poll_job('https://example.test/job/1')


def test_missing_result_resource_is_rejected():
    """A completed job must expose a retrievable result resource."""

    adapter = GaiaFailoverTAP(
        session=_Session((_Response(status=404),)))

    with pytest.raises(GaiaTAPError, match='no result resource'):
        adapter._fetch_result('https://example.test/job/1')


def test_malformed_table_is_rejected_before_find_sources(monkeypatch):
    """A parseable but unexpected table cannot become simulator input."""

    adapter = GaiaFailoverTAP(
        endpoints=(GAIA_TAP_ENDPOINTS[0],), max_retries=1)
    monkeypatch.setattr(
        adapter, '_run_query', lambda endpoint, query: Table({'not': [1]}))

    with pytest.raises(GaiaTAPError, match='missing required Gaia columns'):
        adapter.query_region(CENTER, 1. * u.arcmin)


def test_noirlab_uses_native_job_start_sequence():
    """NOIRLab receives its phase URL and explicit UWS start request."""

    endpoint = GAIA_TAP_ENDPOINTS[1]
    created = _Response(
        status=303, headers={'Location': '/tap/async/42'},
        url=f'{endpoint.url}/async?PHASE=RUN')
    session = _Session((created, _Response()))
    adapter = GaiaFailoverTAP(session=session)

    job_url = adapter._submit_async_job(endpoint, 'SELECT source_id')

    assert job_url == 'https://datalab.noirlab.edu/tap/async/42'
    assert session.calls[0][1].endswith('/async?PHASE=RUN')
    assert 'PHASE' not in session.calls[0][2]['data']
    assert session.calls[1][1].endswith('/tap/async/42/phase')
    assert session.calls[1][2]['data'] == {'PHASE': 'RUN'}

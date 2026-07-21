"""Endpoint-aware Gaia DR3 TAP queries for the contamination tool."""

from dataclasses import dataclass
import io
import logging
import time
from urllib.parse import urljoin

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import MaskedColumn, Table
import numpy as np
import requests


NORMALIZED_COLUMNS = (
    'source_id',
    'designation',
    'ra',
    'dec',
    'ref_epoch',
    'pmra',
    'pmdec',
    'parallax',
    'astrometric_excess_noise',
    'phot_g_mean_flux',
    'phot_g_mean_mag',
    'bp_rp',
    'phot_bp_rp_excess_factor',
    'classprob_dsc_combmod_star',
    'classprob_dsc_combmod_galaxy',
    'classprob_dsc_combmod_quasar',
)
REQUIRED_COLUMNS = frozenset((
    'source_id', 'designation', 'ra', 'dec', 'phot_g_mean_flux'))


@dataclass(frozen=True)
class GaiaTAPEndpoint:
    """Native query and UWS details for one Gaia TAP service."""

    name: str
    url: str
    table: str
    columns: tuple
    cone_predicate: str
    phase_in_url: bool = False
    start_after_submit: bool = False
    constant_columns: tuple = ()

    def build_query(self, ra, dec, radius):
        """Build an explicit native-column cone query."""

        select = ',\n    '.join(
            f'{native} AS {internal}' for native, internal in self.columns)
        predicate = self.cone_predicate.format(
            ra=ra, dec=dec, radius=radius)
        return (
            'SELECT\n    ' + select + f'\nFROM {self.table}\n'
            f'WHERE {predicate}')


_GAIA_COLUMNS = tuple((name, name) for name in NORMALIZED_COLUMNS)
_VIZIER_COLUMNS = (
    ('"Source"', 'source_id'),
    ('"DR3Name"', 'designation'),
    ('"RA_ICRS"', 'ra'),
    ('"DE_ICRS"', 'dec'),
    ('"pmRA"', 'pmra'),
    ('"pmDE"', 'pmdec'),
    ('"Plx"', 'parallax'),
    ('"epsi"', 'astrometric_excess_noise'),
    ('"FG"', 'phot_g_mean_flux'),
    ('"Gmag"', 'phot_g_mean_mag'),
    ('"BP-RP"', 'bp_rp'),
    ('"E(BP/RP)"', 'phot_bp_rp_excess_factor'),
    ('"PSS"', 'classprob_dsc_combmod_star'),
    ('"PGal"', 'classprob_dsc_combmod_galaxy'),
    ('"PQSO"', 'classprob_dsc_combmod_quasar'),
)

GAIA_TAP_ENDPOINTS = (
    GaiaTAPEndpoint(
        name='ESAC',
        url='https://gea.esac.esa.int/tap-server/tap',
        table='gaiadr3.gaia_source',
        columns=_GAIA_COLUMNS,
        cone_predicate=(
            "1=CONTAINS(POINT('ICRS', ra, dec), "
            "CIRCLE('ICRS', {ra}, {dec}, {radius}))"),
    ),
    GaiaTAPEndpoint(
        name='NOIRLab',
        url='https://datalab.noirlab.edu/tap',
        table='gaia_dr3.gaia_source',
        columns=_GAIA_COLUMNS,
        cone_predicate=(
            "'t'=q3c_radial_query(ra, dec, {ra}, {dec}, {radius})"),
        phase_in_url=True,
        start_after_submit=True,
    ),
    GaiaTAPEndpoint(
        name='TAPVizieR',
        url='https://tapvizier.cds.unistra.fr/TAPVizieR/tap',
        table='"I/355/gaiadr3"',
        columns=_VIZIER_COLUMNS,
        cone_predicate=(
            "1=CONTAINS(POINT('ICRS', \"RA_ICRS\", \"DE_ICRS\"), "
            "CIRCLE('ICRS', {ra}, {dec}, {radius}))"),
        # Gaia DR3 astrometry is at J2016.0; TAPVizieR omits the
        # catalog-level reference epoch from its native table.
        constant_columns=(('ref_epoch', 2016.0),),
    ),
)


class GaiaTAPError(RuntimeError):
    """A Gaia query failed at one or all configured TAP services."""


class GaiaFailoverTAP:
    """Query endpoint-native Gaia tables and normalize their results."""

    def __init__(
            self, timeout=30, poll_interval=1.0, max_polls=120,
            max_retries=3, retry_delay=5, endpoints=GAIA_TAP_ENDPOINTS,
            session=None):
        self.endpoints = tuple(endpoints)
        self.timeout = timeout
        self.poll_interval = poll_interval
        self.max_polls = max_polls
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.session = session or requests
        self.last_endpoint = None

    def query_region(self, coordinate, width, height=None):
        """Return a normalized table for the enclosing-circle query."""

        height = width if height is None else height
        ra = float(coordinate.ra.deg)
        dec = float(coordinate.dec.deg)
        width_deg = float(width.to_value(u.deg))
        height_deg = float(height.to_value(u.deg))
        radius = 0.5 * np.sqrt(width_deg**2 + height_deg**2)
        if not np.isfinite(radius) or radius <= 0:
            raise ValueError(f'Invalid radius: {radius}')

        errors = []
        last_error = None
        for retry in range(self.max_retries):
            for endpoint in self.endpoints:
                logging.info('[Gaia TAP] Attempting %s (%s)',
                             endpoint.name, endpoint.url)
                try:
                    query = endpoint.build_query(ra, dec, radius)
                    native = self._run_query(endpoint, query)
                    result = self._normalize(native, endpoint, coordinate)
                    self.last_endpoint = endpoint.url
                    logging.info('[Gaia TAP] %s succeeded with %d rows',
                                 endpoint.name, len(result))
                    return result
                except Exception as exc:
                    last_error = exc
                    reason = ' '.join(str(exc).split()) or type(exc).__name__
                    errors.append(f'{endpoint.name}: {reason}')
                    logging.warning('[Gaia TAP] %s failed: %s',
                                    endpoint.name, reason)
            if retry < self.max_retries - 1:
                time.sleep(self.retry_delay * (2**retry))

        attempted = '; '.join(errors)
        raise GaiaTAPError(
            f'All Gaia TAP endpoints failed: {attempted}') from last_error

    def _run_query(self, endpoint, query):
        job_url = self._submit_async_job(endpoint, query)
        self._poll_job(job_url)
        return self._fetch_result(job_url)

    def _submit_async_job(self, endpoint, query):
        submit_url = f'{endpoint.url}/async'
        payload = {
            'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv',
            'QUERY': query,
        }
        if endpoint.phase_in_url:
            submit_url += '?PHASE=RUN'
        else:
            payload['PHASE'] = 'RUN'

        response = self.session.post(
            submit_url, data=payload, timeout=self.timeout,
            allow_redirects=False)
        response.raise_for_status()
        location = response.headers.get('Location')
        if not location:
            detail = ' '.join(response.text.split())[:300]
            raise GaiaTAPError(
                f'no TAP job URL returned (HTTP {response.status_code}): '
                f'{detail}')
        job_url = urljoin(response.url, location).rstrip('/')

        if endpoint.start_after_submit:
            start = self.session.post(
                f'{job_url}/phase', data={'PHASE': 'RUN'},
                timeout=self.timeout)
            start.raise_for_status()
        return job_url

    def _poll_job(self, job_url):
        for _ in range(self.max_polls):
            response = self.session.get(
                f'{job_url}/phase', timeout=self.timeout)
            response.raise_for_status()
            phase = response.text.strip().upper()
            if phase == 'COMPLETED':
                return
            if phase in ('ERROR', 'ABORTED'):
                detail = ''
                try:
                    error = self.session.get(
                        f'{job_url}/error', timeout=self.timeout)
                    if error.ok:
                        detail = ' '.join(error.text.split())[:300]
                except requests.RequestException:
                    pass
                suffix = f': {detail}' if detail else ''
                raise GaiaTAPError(f'TAP job {phase}{suffix}')
            time.sleep(self.poll_interval)
        raise TimeoutError('Gaia TAP polling timeout')

    def _fetch_result(self, job_url):
        response = self.session.get(
            f'{job_url}/results/result', timeout=self.timeout)
        if response.status_code == 404:
            raise GaiaTAPError('completed TAP job has no result resource')
        response.raise_for_status()
        if not response.content.strip():
            raise GaiaTAPError('completed TAP job returned an empty result')
        upper = response.content[:2000].upper()
        if (b'QUERY_STATUS' in upper and b'ERROR' in upper
                and b'VOTABLE' in upper):
            raise GaiaTAPError('result resource contains a TAP error document')
        try:
            return Table.read(io.BytesIO(response.content), format='ascii.csv')
        except Exception as exc:
            raise GaiaTAPError(f'malformed TAP result table: {exc}') from exc

    @staticmethod
    def _normalize(table, endpoint, coordinate):
        """Map a native result table onto the contamination schema."""

        result = table.copy(copy_data=True)
        lower_names = {name.lower(): name for name in result.colnames}
        native_by_internal = {
            native.strip('"'): internal
            for native, internal in endpoint.columns}
        for native, internal in native_by_internal.items():
            if internal in result.colnames:
                continue
            actual = lower_names.get(native.lower())
            if actual is not None:
                result.rename_column(actual, internal)

        for name, value in endpoint.constant_columns:
            if name not in result.colnames:
                result.add_column(
                    np.full(len(result), value, dtype=float), name=name)

        missing_required = sorted(REQUIRED_COLUMNS - set(result.colnames))
        if missing_required:
            raise GaiaTAPError(
                'result missing required Gaia columns: '
                + ', '.join(missing_required))

        for name in NORMALIZED_COLUMNS:
            if name not in result.colnames:
                dtype = np.int64 if name == 'source_id' else np.float64
                result.add_column(MaskedColumn(
                    np.zeros(len(result), dtype=dtype), mask=True, name=name))

        result = result[list(NORMALIZED_COLUMNS)]
        try:
            positions = SkyCoord(
                np.asarray(result['ra'], dtype=float) * u.deg,
                np.asarray(result['dec'], dtype=float) * u.deg)
            distance = coordinate.separation(positions).to_value(u.deg)
        except Exception as exc:
            raise GaiaTAPError(
                f'invalid Gaia coordinates in result: {exc}') from exc
        result.add_column(distance, name='dist')
        result.sort(['dist', 'source_id'])
        return result

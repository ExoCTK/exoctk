import datetime
import os
import re
import warnings

from astropy.time import Time
from erfa import ErfaWarning

from bokeh.models import Band, ColumnDataSource, HoverTool
from bokeh.plotting import figure, show

import numpy as np
import requests

from jwst_gtvt.constants import URL
from jwst_gtvt.jwst_tvt import Ephemeris, LAUNCH_DATE
from jwst_gtvt.plotting import get_visibility_windows


HORIZONS_TIMEOUT = (5, 30)
VISIBILITY_RANGE_YEARS = 2
LOCAL_EPHEMERIS_PATTERN = re.compile(
    r'ephemeris_(\d{4}-\d{2}-\d{2})_(\d{4}-\d{2}-\d{2})\.txt$')


class EphemerisUnavailableError(RuntimeError):
    """Raised when neither Horizons nor a complete local ephemeris is available."""


class InvalidHorizonsResponseError(RuntimeError):
    """Raised when Horizons responds without ephemeris vector data."""


def visibility_date_range(reference_date=None):
    """Return a rolling visibility range starting on the current UTC date."""

    if reference_date is None:
        reference_date = datetime.datetime.now(datetime.timezone.utc).date()

    try:
        end_date = reference_date.replace(
            year=reference_date.year + VISIBILITY_RANGE_YEARS)
    except ValueError:
        # Map February 29 to February 28 when the end year is not a leap year.
        end_date = reference_date.replace(
            year=reference_date.year + VISIBILITY_RANGE_YEARS,
            month=2, day=28)

    return Time(reference_date.isoformat()), Time(end_date.isoformat())


class BoundedEphemeris(Ephemeris):
    """Retrieve Horizons data with bounded requests and a safe local fallback."""

    request_timeout = HORIZONS_TIMEOUT

    def __init__(self, start_date=None, end_date=None):
        if start_date is None or end_date is None:
            default_start, default_end = visibility_date_range()
            start_date = default_start if start_date is None else start_date
            end_date = default_end if end_date is None else end_date

        start_date = Time(start_date)
        end_date = Time(end_date)
        self._requested_end_date = end_date.strftime('%Y-%m-%d')
        super().__init__(start_date=start_date, end_date=end_date)

    def ephemeris_maximum_date(self):
        """Retrieve the last available date without waiting indefinitely."""

        request_url = URL.format(LAUNCH_DATE, '9999-01-01')
        try:
            response = requests.get(request_url, timeout=self.request_timeout)
            response.raise_for_status()
            match = re.search(
                r'after A\.D\. (\d{4}-[a-zA-Z]+-\d{1,2})',
                response.text)
            if match is None:
                raise ValueError('Unable to parse the Horizons maximum date')
            return datetime.datetime.strptime(
                match.group(1), '%Y-%b-%d').strftime('%Y-%m-%d')
        except (requests.RequestException, ValueError):
            # Let the bounded data request determine whether this date is
            # available instead of relying on a maximum date that goes stale.
            return self._requested_end_date

    def get_ephemeris_data(self, start_date=None, end_date=None):
        """Retrieve ephemeris data or use a local file with full coverage."""

        self.url = URL.format(start_date, end_date)
        try:
            self.eph_request = requests.get(
                self.url, timeout=self.request_timeout)
            self.eph_request.raise_for_status()
            ephemeris = np.asarray(self.eph_request.text.splitlines())
            if '$$SOE' not in ephemeris or '$$EOE' not in ephemeris:
                raise InvalidHorizonsResponseError(
                    'Horizons response did not contain ephemeris vector data')
            return ephemeris
        except (requests.RequestException,
                InvalidHorizonsResponseError) as exc:
            return self._read_local_ephemeris(start_date, end_date, exc)

    def _read_local_ephemeris(self, start_date, end_date, request_error):
        """Read the packaged ephemeris only when it covers the request."""

        filename = os.path.basename(self.ephemeris_filename)
        match = LOCAL_EPHEMERIS_PATTERN.fullmatch(filename)
        if match is None:
            raise EphemerisUnavailableError(
                'Horizons is unavailable and the packaged ephemeris '
                f'filename does not declare its coverage: {filename}'
            ) from request_error

        local_start, local_end = match.groups()
        if start_date < local_start or end_date > local_end:
            raise EphemerisUnavailableError(
                'Horizons is unavailable and the packaged ephemeris covers '
                f'{local_start} through {local_end}, not the requested '
                f'{start_date} through {end_date}.'
            ) from request_error

        warnings.warn(
            f'Horizons is unavailable; using packaged ephemeris {filename}.',
            RuntimeWarning)
        with open(self.ephemeris_filename) as handle:
            return np.asarray(handle.read().splitlines())


def get_exoplanet_positions(ra, dec, in_FOR=None):
    """Use the jwst_gtvt to obtain positions of exoplanet.
    """

    # ***** HACK HACK HACK HACK HACK *****
    # I don't know why the input RA and DEC values sometimes end in a colon (':'),
    # I only know that, when it does, the ephemeris errors out, and so does the
    # calculation.
    #
    # As such, the only thing I can do right now to fix it is to kill any trailing
    # characters that aren't numeric or decimals.
    #
    # When Joe is able to, this should be fixed at the source.
    # ***** HACK HACK HACK HACK HACK *****
    while ra[-1] not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.']:
        ra = ra[:-1]
    while dec[-1] not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.']:
        dec = dec[:-1]

    # jwst_gtvt validates its future maximum date during construction, which
    # can emit an ERFA warning about uncertain leap seconds. Keep this
    # suppression local so other ERFA warnings remain visible to callers.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message=(r'ERFA function "dtf2d" yielded .*'
                     r'"dubious year \(Note 6\)"'),
            category=ErfaWarning,
        )
        eph = BoundedEphemeris()
    exoplanet_data = eph.get_fixed_target_positions(ra, dec)

    if in_FOR is None:
        return exoplanet_data
    else:
        return exoplanet_data.loc[exoplanet_data['in_FOR']==in_FOR]


def build_visibility_plot(target_name, instrument, ra, dec,
                          exoplanet_df=None):
    """Build bokeh figure for visibility windows
    """

    instrument = instrument.upper()

    if instrument not in ['NIRCAM', 'NIRISS', 'MIRI', 'NIRSPEC']:
        raise ValueError(f'{instrument} not supported for this tool!')

    nominal_angle_column_name = instrument + '_nominal_angle'
    min_pa_column_name = instrument + '_min_pa_angle'
    max_pa_column_name = instrument + '_max_pa_angle'

    # Obtain exoplanet data when it was not already calculated for the table.
    if exoplanet_df is None:
        exoplanet_df = get_exoplanet_positions(ra, dec)
    exoplanet_df = exoplanet_df.loc[exoplanet_df['in_FOR']].copy()
    window_indices = get_visibility_windows(exoplanet_df.index.tolist())

    exoplanet_df['times'] = Time(exoplanet_df['MJD'], format='mjd').datetime

    # source = ColumnDataSource(exoplanet_df)

    # define bokeh figure
    TOOLTIPS = [
    ("Date", "@times{%F}"),
    ("Nominal Position Angle", "@{}".format(nominal_angle_column_name)),
    ("Min Position Angle", "@{}".format(min_pa_column_name)),
    ("Max Position Angle", "@{}".format(max_pa_column_name)),]

    p = figure(title=f"{target_name} Visibility with {instrument}",
               height=400, width=800, x_axis_type='datetime')

    p.xaxis.axis_label = 'Date'
    p.yaxis.axis_label = 'Available Aperture Position Angles (Degrees)'

    for start, end in window_indices:
        data_to_plot = exoplanet_df.loc[start:end]
        source = ColumnDataSource(data_to_plot)

        p.line("times", instrument + "_nominal_angle", line_dash=(10, 7), line_width=1, source=source, legend_label="Nominal Angle")

        band = Band(base='times', lower=instrument + '_min_pa_angle', upper=instrument + '_max_pa_angle', source=source, 
                    level='underlay', fill_alpha=1.0, line_width=4, line_color='green')

        p.add_layout(band)
        p.xaxis.major_label_orientation = 3.14/4

        p.add_tools(HoverTool(tooltips=TOOLTIPS, formatters={'@times': 'datetime'}))

    return p

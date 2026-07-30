import datetime
import os
import re
import warnings

from astropy.time import Time

from bokeh.models import Band, ColumnDataSource, HoverTool
from bokeh.plotting import figure, show

import numpy as np
import requests

from jwst_gtvt.constants import URL
from jwst_gtvt.jwst_tvt import Ephemeris, LAUNCH_DATE
from jwst_gtvt.plotting import get_visibility_windows


HORIZONS_TIMEOUT = (5, 30)
HORIZONS_KNOWN_MAX_DATE = '2030-03-16'
LOCAL_EPHEMERIS_PATTERN = re.compile(
    r'ephemeris_(\d{4}-\d{2}-\d{2})_(\d{4}-\d{2}-\d{2})\.txt$')


class EphemerisUnavailableError(RuntimeError):
    """Raised when neither Horizons nor a complete local ephemeris is available."""


class BoundedEphemeris(Ephemeris):
    """Retrieve Horizons data with bounded requests and a safe local fallback."""

    request_timeout = HORIZONS_TIMEOUT

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
            return HORIZONS_KNOWN_MAX_DATE

    def get_ephemeris_data(self, start_date=None, end_date=None):
        """Retrieve ephemeris data or use a local file with full coverage."""

        self.url = URL.format(start_date, end_date)
        try:
            self.eph_request = requests.get(
                self.url, timeout=self.request_timeout)
            self.eph_request.raise_for_status()
            return np.asarray(self.eph_request.text.splitlines())
        except requests.RequestException as exc:
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

    # Set ephemeris to go from Cycle 3 to Cycle 6:
    eph = BoundedEphemeris(
        start_date=Time('2024-07-30'), end_date=Time('2028-07-30'))
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

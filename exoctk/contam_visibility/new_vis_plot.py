from astropy.time import Time

from bokeh.models import Band, ColumnDataSource, HoverTool
from bokeh.plotting import figure, show

from jwst_gtvt.jwst_tvt import Ephemeris
from jwst_gtvt.plotting import get_visibility_windows


def get_exoplanet_positions(ra, dec, in_FOR=None):
    """Use the jwst_gtvt to obtain positions of exoplanet.
    """

    # Set ephemeris to go from Cycle 3 to Cycle 6:
    eph = Ephemeris(start_date=Time('2024-07-30'), end_date=Time('2028-07-30'))
    exoplanet_data = eph.get_fixed_target_positions(ra, dec)

    if in_FOR is None:
        return exoplanet_data
    else:
        return exoplanet_data.loc[exoplanet_data['in_FOR']==in_FOR]


def build_visibility_plot(target_name, instrument, ra, dec):
    """Build bokeh figure for visibility windows
    """

    instrument = instrument.upper()

    if instrument not in ['NIRCAM', 'NIRISS', 'MIRI', 'NIRSPEC']:
        raise ValueError(f'{instrument} not supported for this tool!')

    nominal_angle_column_name = instrument + '_nominal_angle'
    min_pa_column_name = instrument + '_min_pa_angle'
    max_pa_column_name = instrument + '_max_pa_angle'

    # obtain exoplanet data and filter visibility windows
    exoplanet_df = get_exoplanet_positions(ra, dec, in_FOR=True)
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

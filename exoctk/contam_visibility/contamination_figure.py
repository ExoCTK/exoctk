import os
import sys
from itertools import groupby, count

from astropy.io import fits
from bokeh.layouts import gridplot, column, row
from bokeh.models import Range1d, LinearColorMapper, LogColorMapper, Label, ColorBar, ColumnDataSource, HoverTool, Slider, CustomJS, VArea, CrosshairTool, TapTool, OpenURL, Span, Legend, Spacer
from bokeh.palettes import PuBu, Spectral6
from bokeh.plotting import figure
import numpy as np

from . import visibilityPA as vpa
from ..utils import fill_between
from exoctk.pkgdata import resource_filename


EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
if not EXOCTK_DATA:
    print(
        'WARNING: The $EXOCTK_DATA environment variable is not set. '
        'Contamination overlap will not work. Please set the '
        'value of this variable to point to the location of the exoctk_data '
        'download folder.  Users may retreive this folder by clicking the '
        '"ExoCTK Data Download" button on the ExoCTK website, or by using '
        'the exoctk.utils.download_exoctk_data() function.')
    TRACES_PATH = None

else:
    TRACES_PATH = os.path.join(EXOCTK_DATA, 'exoctk_contam', 'traces')

LAM_FILE = resource_filename('exoctk', 'data/contam_visibility/lambda_order1-2.txt')

disp_nircam = 0.001  # microns
lam0_nircam322w2 = 2.369
lam1_nircam322w2 = 4.417
lam0_nircam444w = 3.063
lam1_nircam444w = 5.111


def has_order2_contamination(instrument):
    """Return whether the contamination calculation includes SOSS Order 2."""

    # SUBSTRIP96 trace products contain only Order 1.  The legacy simulator
    # reuses SUBSTRIP256 templates for source placement, so plotting its Order
    # 2 panel would incorrectly imply a contamination calculation was done.
    return instrument != 'NIS_SUBSTRIP96'


def contam_slider_plot(pctlines, badPA_list, threshold=0.05, y_max=0.1,
                       instrument=None):
    """
    Make the contamination plot with a slider

    Parameters
    ----------
    pctlines: list
        The list of fractional contamination arrays
    badPA_list: array-like
        The position angles that are not observable
    threshold: float
        The threshold for contamination reporting in the plot
    y_max: float
        The fractional contamination at the top of both displayed axes. Values
        above this limit remain in the data and are clipped only visually.
    instrument: str, optional
        Instrument aperture name. ``NIS_SUBSTRIP96`` displays only the
        calculated SOSS Order 1 result.

    Returns
    -------
    bokeh.layouts.column
        The column of plots
    """
    # SUBSTRIP96 trace products provide only SOSS Order 1. The simulator
    # reuses SUBSTRIP256 templates for source placement, so its extra internal
    # trace arrays must not be presented as calculated contamination results.
    if instrument == 'NIS_SUBSTRIP96':
        pctlines = pctlines[:1]

    # Quantities
    pa_list = np.arange(360)
    orders = np.arange(1, len(pctlines)+1)
    n_channels = pctlines[0].shape[1]
    if not 0 < y_max <= 1:
        raise ValueError('y_max must be greater than zero and at most one')
    percent_lines = [np.asarray(line) * 100 for line in pctlines]
    display_y_max = y_max * 100

    # Store individual contamination fractions for each order+PA combination
    contam_dict = {}
    for order in orders:
        for pa, percent in enumerate(percent_lines[order - 1]):
            contam_dict[f'contam{order}_{pa}'] = percent

    # Choose initial values and build visible and available dicts
    pa_init = int(np.argmax(np.nansum(np.asarray(pctlines), axis=(0, 2)))) # Maximum contamination
    vis_dict = {f'contam{order}': contam_dict[f'contam{order}_{pa_init}'] for order in orders}
    vis_dict.update({'col': np.arange(n_channels), 'zeros': np.zeros(n_channels)})
    source_visible = ColumnDataSource(data=vis_dict)
    source_available = ColumnDataSource(data=contam_dict)

    # Contamination fraction plot
    plt = figure(width=900, height=300, tools=['reset', 'save'])
    colors = ['blue', 'red', 'green', 'cyan', 'dodgerblue', 'purple', 'orange', 'lime', 'yellow', 'magenta']
    for order in orders:
        plt.line('col', f'contam{order}', source=source_visible, color=colors[order - 1], line_width=2, line_alpha=0.6, legend_label=f'Order {order}')
        glyph = VArea(x="col", y1="zeros", y2=f"contam{order}", fill_color=colors[order - 1], fill_alpha=0.3)
        plt.add_glyph(source_visible, glyph)

    plt.y_range = Range1d(0, display_y_max)
    plt.x_range = Range1d(0, n_channels)
    plt.xaxis.axis_label = 'Column Index'
    plt.yaxis.axis_label = 'Contamination (%)'
    slider = Slider(title='V3 Position Angle',
                    value=pa_init,
                    start=min(pa_list),
                    end=max(pa_list),
                    step=1,
                    width=830)

    span = Span(line_width=2, location=slider.value, dimension='height')

    # Define CustomJS callback, which updates the plot based on selected function by updating the source_visible
    js_orders = "\n".join([f"data_visible['contam{order}'] = data_available['contam{order}_' + selected_pa];" for order in orders])
    callback = CustomJS(
        args=dict(source_visible=source_visible, source_available=source_available, span=span), code=f"""
            var selected_pa = (cb_obj.value).toString();
            var data_visible = source_visible.data;
            var data_available = source_available.data;
            {js_orders}
            span.location = cb_obj.value;
            source_visible.change.emit();
        """)

    # Make a guide that shows which PAs are unobservable
    viz_none = np.array([
        display_y_max if i in badPA_list else 0 for i in pa_list])
    viz_plt = figure(width=900, height=300, x_range=Range1d(0, 359))
    c0 = viz_plt.vbar(x=np.arange(360), top=viz_none, width=1.1, alpha=0.5, line_color=None, fill_color='#555555')

    # Add contamination for each order
    vis_ords = []
    for order in orders:

        # Plot the mean contamination across spectral channels.
        viz_plt.step(pa_list, np.nanmean(percent_lines[order - 1], axis=1), line_width=2, color=colors[order - 1], mode="center")

        # PLot columns that indicate >5% contamination
        viz_ord = np.array([
            display_y_max if value > threshold else 0
            for value in np.nanmax(pctlines[order - 1], axis=1)
        ])
        vis_ords.append(viz_plt.vbar(x=pa_list, top=viz_ord, width=1.1, alpha=0.2, line_color=None, fill_color=colors[order - 1]))

    # Formatting
    viz_plt.x_range = Range1d(0, 359)
    viz_plt.y_range = Range1d(0, display_y_max)
    viz_plt.xaxis.axis_label = 'V3 Position Angle'
    viz_plt.yaxis.axis_label = 'Mean Contamination (%)'
    viz_plt.add_layout(span)

    items = [("Target not observable", [c0])]
    for order, vis_ord in zip(orders, vis_ords):
        items.append((f'Ord {order} > {threshold * 100}% Contaminated', [vis_ord]))
    legend = Legend(items=items, location=(50, 0), orientation='horizontal', border_line_alpha=0)
    viz_plt.add_layout(legend, 'below')

    # Put plot together
    slider.js_on_change('value', callback)
    slider_row = row(Spacer(width=40), slider, Spacer(width=30))
    layout = column(plt, slider_row, viz_plt)

    return layout


def nirissContam(cube, paRange=[0, 360], lam_file=LAM_FILE):
    """ Generates the contamination figure that will be plotted on the website
    for NIRISS SOSS.
    """
    # Get data from FITS file
    if isinstance(cube, str):
        hdu = fits.open(cubeName)
        cube = hdu[0].data
        hdu.close()

    # Pull out the target trace and cube of neighbor traces
    trace1 = cube[0, :, :]
    trace2 = cube[1, :, :]
    cube = cube[2:, :, :]

    plotPAmin, plotPAmax = paRange

    # Start calculations
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    nPA = cube.shape[0]
    rows = cube.shape[1]
    cols = cube.shape[2]
    dPA = 360 // nPA
    PA = np.arange(nPA) * dPA

    contamO1 = np.zeros([rows, nPA])
    contamO2 = np.zeros([rows, nPA])

    low_lim_col = 20
    high_lim_col = 41

    for row in np.arange(rows):
        # Contamination for order 1 of target trace
        i = np.argmax(trace1[row, :])
        tr = trace1[row, i - low_lim_col:i + high_lim_col]
        w = tr / np.sum(tr**2)
        ww = np.tile(w, nPA).reshape([nPA, tr.size])
        contamO1[row, :] = np.sum(
            cube[:, row, i - low_lim_col:i + high_lim_col] * ww, axis=1)

        # Contamination for order 2 of target trace
        if lamO2[row] < 0.6:
            continue
        i = np.argmax(trace2[row, :])
        tr = trace2[row, i - 20:i + 41]
        w = tr / np.sum(tr**2)
        ww = np.tile(w, nPA).reshape([nPA, tr.size])
        contamO2[row, :] = np.sum(cube[:, row, i - 20:i + 41] * ww, axis=1)

    return contamO1, contamO2


def nircamContam(cube, paRange=[0, 360]):
    """ Generates the contamination figure that will be plotted on the website
    for NIRCam Grism Time Series mode.

    Parameters
    ----------
    cube : arr or str
        A 3D array of the simulated field at every Aperture Position Angle (APA).
        The shape of the cube is (361, subY, subX).
        or
        The name of an HDU .fits file sthat has the cube.

    Returns
    -------
    bokeh plot
    """
    # Get data from FITS file
    if isinstance(cube, str):
        hdu = fits.open(cubeName)
        cube = hdu[0].data
        hdu.close()

    # Pull out the target trace and cube of neighbor traces
    targ = cube[0, :, :]  # target star order 1 trace
    # neighbor star order 1 and 2 traces in all the angles
    cube = cube[1:, :, :]

    # Remove background values < 1 as it can blow up contamination
    targ = np.where(targ < 1, 0, targ)

    PAmin, PAmax = paRange[0], paRange[1]
    PArange = np.arange(PAmin, PAmax, 1)

    nPA, rows, cols = cube.shape[0], cube.shape[1], cube.shape[2]

    contamO1 = np.zeros([nPA, cols])

    # the width of the trace (in Y-direction for NIRCam GTS)
    peak = targ.max()
    low_lim_row = np.where(targ > 0.0001 * peak)[0].min()
    high_lim_row = np.where(targ > 0.0001 * peak)[0].max()

    # the length of the trace (in X-direction for NIRCam GTS)
    targ_trace_start = np.where(targ > 0.0001 * peak)[1].min()
    targ_trace_stop = np.where(targ > 0.0001 * peak)[1].max()

    # Begin contam calculation at each channel (column) X
    for X in np.arange(cols):
        if (X < targ_trace_start) or (X > targ_trace_stop):
            continue

        peakY = np.argmax(targ[:, X])
        TOP, BOT = peakY + high_lim_row, peakY - low_lim_row

        tr = targ[BOT:TOP, X]

        # calculate weights
        wt = tr / np.sum(tr**2)
        ww = np.tile(wt, nPA).reshape([nPA, tr.size])

        contamO1[:, X] = np.sum(cube[:, BOT:TOP, X] * ww, axis=1)

    contamO1 = contamO1[:, targ_trace_start:targ_trace_stop]
    return contamO1


def miriContam(cube, paRange=[0, 360]):
    """ Generates the contamination figure that will be plotted on the website
    for MIRI LRS.
    """
    # Get data from FITS file
    if isinstance(cube, str):
        hdu = fits.open(cubeName)
        cube = hdu[0].data
        hdu.close()

    # Pull out the target trace and cube of neighbor traces
    targ = cube[0, :, :]  # target star order 1 trace
    # neighbor star order 1 and 2 traces in all the angles
    cube = cube[1:, :, :]

    # Remove background values < 1 as it can blow up contamination
    targ = np.where(targ < 1, 0, targ)

    PAmin, PAmax = paRange[0], paRange[1]
    PArange = np.arange(PAmin, PAmax, 1)

    nPA, rows, cols = cube.shape[0], cube.shape[1], cube.shape[2]

    contamO1 = np.zeros([rows, nPA])

    # the width of the trace (in Y-direction for NIRCam GTS)
    peak = targ.max()

    low_lim_col = np.where(targ > 0.0001 * peak)[1].min()
    high_lim_col = np.where(targ > 0.0001 * peak)[1].max()

    # the length of the trace (in X-direction for NIRCam GTS)
    targ_trace_start = np.where(targ > 0.0001 * peak)[0].min()
    targ_trace_stop = np.where(targ > 0.0001 * peak)[0].max()
    # Begin contam calculation at each channel (row) Y
    for Y in np.arange(rows):
        if (Y < targ_trace_start) or (Y > targ_trace_stop):
            continue

        peakX = np.argmax(targ[Y, :])
        LEFT, RIGHT = peakX - low_lim_col, peakX + high_lim_col

        tr = targ[Y, LEFT:RIGHT]

        # calculate weights
        wt = tr / np.sum(tr**2)
        ww = np.tile(wt, nPA).reshape([nPA, tr.size])

        contamO1[Y, :] = np.sum(cube[:, Y, LEFT:RIGHT] * wt,
                                where=~np.isnan(cube[:, Y, LEFT:RIGHT] * wt),
                                axis=1)

        #target = np.sum(cube[0, Y, LEFT:RIGHT], axis=0)

        # contamO1[Y, :] = np.sum(cube[:, Y, LEFT:RIGHT]*ww,
        #                        where=~np.isnan(cube[:, Y, LEFT:RIGHT]),
        #                        axis=1)#/target
    contamO1 = contamO1[targ_trace_start:targ_trace_stop, :]
    return contamO1


def contam(cube, instrument, targetName='noName', paRange=[0, 360], badPAs=[]):
    print("Started contam")

    rows, cols = cube.shape[1], cube.shape[2]

    PAmin, PAmax = paRange[0], paRange[1]
    PA = np.arange(PAmin, PAmax, 1)

    # Generate the contam figure
    if instrument in ['NIS_SUBSTRIP256', 'NIS_SUBSTRIP96']:
        contamO1, contamO2 = nirissContam(cube, lam_file=LAM_FILE)
        ypix, lamO1, lamO2 = np.loadtxt(LAM_FILE, unpack=True)
        xlim0 = lamO1.min()
        xlim1 = lamO1.max()
    elif instrument == 'NRCA5_GRISM256_F444W':
        contamO1 = nircamContam(cube)
        xlim0 = lam0_nircam444w
        xlim1 = lam1_nircam444w
    elif instrument == 'NRCA5_GRISM256_F322W2':
        contamO1 = nircamContam(cube)
        xlim0 = lam0_nircam322w2
        xlim1 = lam1_nircam322w2
    elif instrument == 'MIRIM_SLITLESSPRISM':
        contamO1 = miriContam(cube)
        xlim0 = 5
        xlim1 = 12

    TOOLS = 'pan, box_zoom, reset, save'
    dPA = 1
    print("Contam figure setup done")

    # Order 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Contam plot
    ylim0 = PAmin - 0.5
    ylim1 = PAmax + 0.5
    color_mapper = LinearColorMapper(palette=PuBu[8][::-1][2:], low=-4, high=1)
    color_mapper.low_color = 'white'
    color_mapper.high_color = 'black'
    width = Span(dimension="width", line_width=1)
    height = Span(dimension="height", line_width=1)

    orders = 'Orders 1 & 2' if instrument.startswith('NRCA') else 'Order 1'
    s2 = figure(tools=TOOLS, width=800, height=600, title='{} {} Contamination with {}'.format(orders, targetName, instrument), x_range=Range1d(xlim0, xlim1), y_range=Range1d(ylim0, ylim1))
    o1_crosshair = CrosshairTool(overlay=[width, height])
    s2.add_tools(o1_crosshair)

    contamO1 = contamO1 if 'NRCA' in instrument else contamO1.T
    contamO1 = np.fliplr(contamO1) if (instrument == 'MIRIM_SLITLESSPRISM') or (instrument == 'NRCA5_GRISM256_F322W2') else contamO1
    fig_data = np.log10(np.clip(contamO1, 1.e-10, 1.))
    print("Contam order 1 started")

    # Begin plotting ~~~~~~~~~~~~~~~~~~~~~~~~

    s2.image([fig_data], x=xlim0, y=ylim0, dw=xlim1 - xlim0, dh=ylim1 - ylim0, color_mapper=color_mapper)
    s2.xaxis.axis_label = 'Wavelength (um)'
    s2.yaxis.axis_label = 'Aperture Position Angle (degrees)'

    print("Contam plot started")

    # Line plot
    #ax = 1 if 'NIRCam' in instrument else 0
    channels = cols if 'NRCA' in instrument else rows
    s3 = figure(tools=TOOLS, width=150, height=600, x_range=Range1d(0, 100), y_range=s2.y_range, title=None)
    s3.add_tools(o1_crosshair)

    try:
        s3.line(100 * np.sum(contamO1 >= 0.001, axis=1) / channels, PA - dPA / 2, line_color='blue', legend_label='> 0.001')
        s3.line(100 * np.sum(contamO1 >= 0.01, axis=1) / channels, PA - dPA / 2, line_color='green', legend_label='> 0.01')
    except AttributeError:
        s3.line(100 * np.sum(contamO1 >= 0.001, axis=1) / channels, PA - dPA / 2, line_color='blue', legend='> 0.001')
        s3.line(100 * np.sum(contamO1 >= 0.01, axis=1) / channels, PA - dPA / 2, line_color='green', legend='> 0.01')

    s3.xaxis.axis_label = '% channels contam.'
    s3.yaxis.major_label_text_font_size = '0pt'
    s3.ygrid.grid_line_color = None

    print("Contam done line plot")

    # Add shaded region for bad PAs
    bad_PA_color = '#555555'
    bad_PA_alpha = 0.6
    if len(badPAs) > 0:

        # Group bad PAs
        badPA_groups = [list(map(int, g)) for _, g in groupby(badPAs, lambda n, c=count(): n-next(c))]

        tops, bottoms, lefts, rights, lefts_line, rights_line = [], [], [], [], [], []
        for idx in range(0, len(badPA_groups)):
            PAgroup = badPA_groups[idx]
            top_idx = np.max(PAgroup)
            bot_idx = np.min(PAgroup)

            tops.append(top_idx)
            bottoms.append(bot_idx)
            lefts.append(xlim0)
            rights.append(xlim1)
            lefts_line.append(0)
            rights_line.append(100)

        s2.quad(top=tops, bottom=bottoms, left=lefts, right=rights, color=bad_PA_color, alpha=bad_PA_alpha)
        s3.quad(top=tops, bottom=bottoms, left=lefts_line, right=rights_line, color=bad_PA_color, alpha=bad_PA_alpha)

    print("Contam order 1 done")

    # ~~~~~~ Order 2 ~~~~~~

    # Contam plot
    if instrument.startswith('NIS') and has_order2_contamination(instrument):
        xlim0 = lamO2.min()
        xlim1 = lamO2.max()
        ylim0 = PA.min() - 0.5 * dPA
        ylim1 = PA.max() + 0.5 * dPA
        xlim0 = 0.614
        s5 = figure(tools=TOOLS, width=800, height=600, title='Order 2 {} Contamination with {}'.format(targetName, instrument), x_range=Range1d(xlim0, xlim1), y_range=s2.y_range)
        fig_data = np.log10(np.clip(contamO2.T, 1.e-10, 1.))[:, 300:]
        s5.image([fig_data], x=xlim0, y=ylim0, dw=xlim1 - xlim0, dh=ylim1 - ylim0, color_mapper=color_mapper)
        s5.xaxis.axis_label = 'Wavelength (um)'
        s5.yaxis.axis_label = 'Aperture Position Angle (degrees)'
        o2_crosshair = CrosshairTool(overlay=[width, height])
        s5.add_tools(o2_crosshair)

        # Line plot
        s6 = figure(tools=TOOLS, width=150, height=600, y_range=s2.y_range, x_range=Range1d(0, 100), title=None)
        s6.add_tools(o2_crosshair)

        try:
            s6.line(100 * np.sum(contamO2 >= 0.001, axis=0) / rows, PA - dPA / 2, line_color='blue', legend_label='> 0.001')
            s6.line(100 * np.sum(contamO2 >= 0.01, axis=0) / rows, PA - dPA / 2, line_color='green', legend_label='> 0.01')
        except AttributeError:
            s6.line(100 * np.sum(contamO2 >= 0.001, axis=0) / rows, PA - dPA / 2, line_color='blue', legend='> 0.001')
            s6.line(100 * np.sum(contamO2 >= 0.01, axis=0) / rows, PA - dPA / 2, line_color='green', legend='> 0.01')

        s6.xaxis.axis_label = '% channels contam.'
        s6.yaxis.major_label_text_font_size = '0pt'
        s6.ygrid.grid_line_color = None

        # Add shaded region for bad PAs
        if len(badPAs) > 0:

            # Group bad PAs
            badPA_groups = [list(map(int, g)) for _, g in groupby(badPAs, lambda n, c=count(): n - next(c))]

            tops, bottoms, lefts, rights, lefts_line, rights_line = [], [], [], [], [], []
            for idx in range(0, len(badPA_groups)):
                PAgroup = badPA_groups[idx]
                top_idx = np.max(PAgroup)
                bot_idx = np.min(PAgroup)

                tops.append(top_idx)
                bottoms.append(bot_idx)
                lefts.append(xlim0)
                rights.append(xlim1)
                lefts_line.append(0)
                rights_line.append(100)

            s5.quad(top=tops, bottom=bottoms, left=lefts, right=rights, color=bad_PA_color, alpha=bad_PA_alpha)
            s6.quad(top=tops, bottom=bottoms, left=lefts_line, right=rights_line, color=bad_PA_color,
                    alpha=bad_PA_alpha)

    # ~~~~~~ Plotting ~~~~~~

    if instrument.startswith('NIS') and has_order2_contamination(instrument):
        fig = gridplot(children=[[s2, s3], [s5, s6]])
    else:
        fig = gridplot(children=[[s2, s3]])

    return fig  # , contamO1


if __name__ == "__main__":
    # arguments RA & DEC, conversion to radians
    argv = sys.argv

    ra = argv[1]
    dec = argv[2]
    cubeNameSuf = argv[3]

    pamin = 0 if len(argv) < 5 else int(argv[4])
    pamax = 360 if len(argv) < 6 else int(argv[5])

    cubeName = argv[6]
    targetName = None if len(argv) < 8 else argv[7]
    save = False if len(argv) < 8 else True  # if name provided -> save
    tmpDir = "." if len(argv) < 9 else argv[8]
    os.makedirs(tmpDir, exist_ok=True)

    goodPA, badPA, _ = vpa.checkVisPA(ra, dec, targetName)

    contam(cubeName, targetName=targetName, paRange=[pamin, pamax],
           badPA=badPA, tmpDir=tmpDir)

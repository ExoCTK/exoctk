import datetime
from functools import wraps
import os
import json
from pkg_resources import resource_filename

import astropy.constants as constants
from astropy.extern.six.moves import StringIO
import astropy.table as at
import astropy.units as u
from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8
from bokeh.embed import components
from bokeh.models import Range1d
from bokeh.models.widgets import Panel, Tabs
from bokeh.plotting import figure
import flask
from flask import Flask, Response
from flask import request, send_file, make_response, render_template
import form_validation as fv
import h5py
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
from wtforms.validators import InputRequired, NumberRange
from wtforms import DecimalField

from exoctk.contam_visibility import resolve
from exoctk.contam_visibility import visibilityPA as vpa
from exoctk.contam_visibility import sossFieldSim as fs
from exoctk.contam_visibility import sossContamFig as cf
from exoctk.forward_models.forward_models import fortney_grid, generic_grid
from exoctk.groups_integrations.groups_integrations import perform_calculation
from exoctk.limb_darkening import limb_darkening_fit as lf
from exoctk.modelgrid import ModelGrid
from exoctk.utils import find_closest, filter_table, get_target_data, get_canonical_name, FORTGRID_DIR
import log_exoctk
from svo_filters import svo

# FLASK SET UP
app_exoctk = Flask(__name__)

# define the cache config keys, remember that it can be done in a settings file
app_exoctk.config['CACHE_TYPE'] = 'null'
app_exoctk.config['SECRET_KEY'] = 'Thisisasecret!'

EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
if EXOCTK_DATA == '':
    raise NameError("You need to have an exported 'EXOCTK_DATA' environment variable and data set up before we can continue.")

EXOCTKLOG_DIR = os.path.join(EXOCTK_DATA, 'exoctk_log')
GROUPS_INTEGRATIONS_DIR = os.path.join(EXOCTK_DATA, 'groups_integrations/groups_integrations.json')
MODELGRID_DIR = os.path.join(EXOCTK_DATA, 'modelgrid/default')

# # Load the database to log all form submissions
# if EXOCTKLOG_DIR is None:
#     dbpath = ':memory:'
# else:
#     dbpath = os.path.realpath(os.path.join(EXOCTKLOG_DIR, 'exoctk_log.db'))
#     if not os.path.isfile(dbpath):
#         log_exoctk.create_db(dbpath)
# try:
#     DB = log_exoctk.load_db(dbpath)
# except IOError:
#     DB = None


# Redirect to the index
@app_exoctk.route('/')
@app_exoctk.route('/index')
def index():
    """The Index page"""
    return render_template('index.html')


@app_exoctk.route('/limb_darkening', methods=['GET', 'POST'])
def limb_darkening():
    """The limb darkening form page"""
    # Load default form
    form = fv.LimbDarkeningForm()

    # Reload page with stellar data from ExoMAST
    if form.resolve_submit.data: 

        # Resolve the target in exoMAST
        data = get_target_data(form.targname.data)

        # Update the form data
        form.feh.data = float(data['Fe/H'])
        form.teff.data = float(data['Teff'])
        form.logg.data = float(data['stellar_gravity'])

        # Send it back to the main page
        return render_template('limb_darkening.html', form=form)

    # Reload page with appropriate filter data
    if form.filter_submit.data:

        kwargs = {}
        if form.bandpass.data == 'tophat':
            kwargs['n_bins'] = 1
            kwargs['pixels_per_bin'] = 100
            kwargs['wave_min'] = 1*u.um
            kwargs['wave_max'] = 2*u.um

        # Get the filter
        bandpass = svo.Filter(form.bandpass.data, **kwargs)

        # Update the form data
        form.wave_min.data = bandpass.wave_min.value
        form.wave_max.data = bandpass.wave_max.value

        # Send it back to the main page
        return render_template('limb_darkening.html', form=form)

    # Update validation values after a model grid is selected
    if form.modelgrid_submit.data:

        # Load the modelgrid
        mg = ModelGrid(form.modeldir.data, resolution=500)
        teff_rng = mg.Teff_vals.min(), mg.Teff_vals.max()
        logg_rng = mg.logg_vals.min(), mg.logg_vals.max()
        feh_rng = mg.FeH_vals.min(), mg.FeH_vals.max()

        # Update the validation parameters by setting validator attributes
        setattr(form.teff.validators[1], 'min', teff_rng[0])
        setattr(form.teff.validators[1], 'max', teff_rng[1])
        setattr(form.teff.validators[1], 'message', 'Effective temperature must be between {} and {}'.format(*teff_rng))
        setattr(form.logg.validators[1], 'min', logg_rng[0])
        setattr(form.logg.validators[1], 'max', logg_rng[1])
        setattr(form.logg.validators[1], 'message', 'Surface gravity must be between {} and {}'.format(*logg_rng))
        setattr(form.feh.validators[1], 'min', feh_rng[0])
        setattr(form.feh.validators[1], 'max', feh_rng[1])
        setattr(form.feh.validators[1], 'message', 'Metallicity must be between {} and {}'.format(*feh_rng))

        # Send it back to the main page
        return render_template('limb_darkening.html', form=form)

    # Validate form and submit for results
    if form.validate_on_submit() and form.calculate_submit.data:

        # Get the stellar parameters
        star_params = [form.teff.data, form.logg.data, form.feh.data]

        # Log the form inputs
        try:
            log_exoctk.log_form_input(request.form, 'limb_darkening', DB)
        except:
            pass

        # Load the model grid
        model_grid = ModelGrid(form.modeldir.data, resolution=500)
        form.modeldir.data = [j for i, j in form.modeldir.choices if i == form.modeldir.data][0]

        # Grism details
        if '.G' in form.bandpass.data.upper() and 'GAIA' not in form.bandpass.data.upper():
            kwargs = {'n_bins': form.n_bins.data, 'pixels_per_bin': form.n_pix.data,
                      'wl_min': form.wave_min.data*u.um, 'wl_max': form.wave_max.data*u.um}
        else:
            kwargs = {}

        # Make filter object and plot
        bandpass = svo.Filter(form.bandpass.data, **kwargs)
        bp_name = bandpass.name
        bk_plot = bandpass.plot(draw=False)
        bk_plot.plot_width = 580
        bk_plot.plot_height = 280
        js_resources = INLINE.render_js()
        css_resources = INLINE.render_css()
        filt_script, filt_plot = components(bk_plot)

        # Trim the grid to nearby grid points to speed up calculation
        full_rng = [model_grid.Teff_vals, model_grid.logg_vals, model_grid.FeH_vals]
        trim_rng = find_closest(full_rng, star_params, n=1, values=True)

        # Calculate the coefficients for each profile
        ld = lf.LDC(model_grid)
        for prof in form.profiles.data:
            ld.calculate(*star_params, prof, mu_min=form.mu_min.data, bandpass=bandpass)

        # Draw a figure for each wavelength bin
        tabs = []
        for wav in np.unique(ld.results['wave_eff']):

            # Plot it
            TOOLS = 'box_zoom, box_select, crosshair, reset, hover'
            fig = figure(tools=TOOLS, x_range=Range1d(0, 1), y_range=Range1d(0, 1),
                        plot_width=800, plot_height=400)
            ld.plot(wave_eff=wav, fig=fig)

            # Plot formatting
            fig.legend.location = 'bottom_right'
            fig.xaxis.axis_label = 'mu'
            fig.yaxis.axis_label = 'Intensity'

            tabs.append(Panel(child=fig, title=str(wav)))

        final = Tabs(tabs=tabs)

        # Get HTML
        script, div = components(final)

        # Store the tables as a string
        file_as_string = str(ld.results[[c for c in ld.results.dtype.names if
                                        ld.results.dtype[c] != object]])

        # Make a table for each profile with a row for each wavelength bin
        profile_tables = []
        for profile in form.profiles.data:

            # Make LaTeX for polynomials
            latex = lf.ld_profile(profile, latex=True)
            poly = '\({}\)'.format(latex).replace('*', '\cdot').replace('\e', 'e')

            # Make the table into LaTeX
            table = filter_table(ld.results, profile=profile)
            co_cols = [c for c in ld.results.colnames if (c.startswith('c') or
                    c.startswith('e')) and len(c) == 2 and not
                    np.all([np.isnan(i) for i in table[c]])]
            table = table[['wave_min', 'wave_max'] + co_cols]
            table.rename_column('wave_min', '\(\lambda_\mbox{min}\hspace{5px}(\mu m)\)')
            table.rename_column('wave_max', '\(\lambda_\mbox{max}\hspace{5px}(\mu m)\)')

            # Add the results to the lists
            html_table = '\n'.join(table.pformat(max_width=500, html=True))\
                        .replace('<table', '<table id="myTable" class="table table-striped table-hover"')

            # Add the table title
            header = '<br></br><strong>{}</strong><br><p>\(I(\mu)/I(\mu=1)\) = {}</p>'.format(profile, poly)
            html_table = header + html_table

            profile_tables.append(html_table)

        return render_template('limb_darkening_results.html', form=form,
                               table=profile_tables, script=script, plot=div,
                               file_as_string=repr(file_as_string),
                               filt_plot=filt_plot, filt_script=filt_script,
                               js=js_resources, css=css_resources)

    return render_template('limb_darkening.html', form=form)


@app_exoctk.route('/limb_darkening_error', methods=['GET', 'POST'])
def limb_darkening_error():
    """The limb darkening error page"""

    return render_template('limb_darkening_error.html')


@app_exoctk.route('/groups_integrations', methods=['GET', 'POST'])
def groups_integrations():
    """The groups and integrations calculator form page"""

    # Print out pandeia sat values
    with open(resource_filename('exoctk', 'data/groups_integrations/groups_integrations_input_data.json')) as f:
        sat_data = json.load(f)['fullwell']

    if request.method == 'POST':
        if request.form['submit'] == "Retrieve Parameters":
            target_name = request.form['targetname']
            canoncial_name = get_canonical_name(target_name)
            # Ping exomast api and get data
            data = get_target_data(target_name)
            Kmag = data['Kmag']
            obs_duration = data['transit_duration'] * 24. # Transit duration in exomast is in days, need it in hours
            
            groupsintegrationVars = {'targname':canoncial_name, 'Kmag':Kmag, 'obs_duration':obs_duration}

            return render_template('groups_integrations.html', sat_data=sat_data, groupsintegrationVars=groupsintegrationVars)

    return render_template('groups_integrations.html', sat_data=sat_data)


@app_exoctk.route('/groups_integrations_results', methods=['GET', 'POST'])
def groups_integrations_results():
    """The groups and integrations calculator results page"""

    # Read in parameters from form
    params = {}
    for key in dict(request.form).keys():
        params[key] = dict(request.form)[key][0]
    try:
        err = 0

        # Specific error catching
        if params['n_group'].isdigit():
            params['n_group'] = int(params['n_group'])
            if params['n_group'] <= 1:
                err = 'Please try again with at least one group.'
            else:
                if params['n_group'] != 'optimize':
                    err = "You need to double check your group input. Please put the number of groups per integration or 'optimize' and we can calculate it for you."

        if (False in [params['mag'].isdigit(), params['obs_time'].isdigit()]) and ('.' not in params['mag']) and ('.' not in params['obs_time']):
            err = 'Your magnitude or observation time is not a number, or you left the field blank.'

        else:
            if (4.5 > float(params['mag'])) or (12.5 < float(params['mag'])):
                err = 'Looks like we do not have useful approximations for your magnitude. Could you give us a number between 5.5-12.5?'
            if float(params['obs_time']) <= 0:
                err = 'You have a negative transit time -- I doubt that will be of much use to anyone.'

        if float(params['sat_max']) <= 0:
            err = 'You put in an underwhelming saturation level. There is something to be said for being too careful...'
        if (params['sat_mode'] == 'well') and float(params['sat_max']) > 1:
            err = 'You are saturating past the full well. Is that a good idea?'

        if type(err) == str:
            return render_template('groups_integrations_error.html', err=err)

        # Only create the dict if the form input looks okay
        # Make sure everything is the right type
        ins = params['ins']
        float_params = ['obs_time', 'mag', 'sat_max']
        str_params = ['mod', 'band', 'time_unit', '{}_filt'.format(ins),
                      '{}_ta_filt'.format(ins), 'ins',
                      '{}_subarray'.format(ins), '{}_subarray_ta'.format(ins),
                      'sat_mode']
        for key in params:
            if key in float_params:
                params[key] = float(params[key])
            if key in str_params:
                params[key] = str(params[key])

        # Also get the data path in there
        params['infile'] = resource_filename('exoctk', 'data/groups_integrations/groups_integrations_input_data.json')

        # Rename the ins-mode params to more general counterparts
        params['filt'] = params['{}_filt'.format(ins)]
        params['filt_ta'] = params['{}_filt_ta'.format(ins)]
        params['subarray'] = params['{}_subarray'.format(ins)]
        params['subarray_ta'] = params['{}_subarray_ta'.format(ins)]
        
        # Convert the obs_time to hours
        if params['time_unit'] != 'hours':
            params['obs_time'] = params['obs_time']*24
            params['time_unit'] = 'hours'

        results = perform_calculation(params)

        if type(results) == dict:
            results_dict = results
            one_group_error = ""
            zero_group_error = ""
            if results_dict['n_group'] == 1:
                one_group_error = 'Be careful! This only predicts one group, and you may be in danger of oversaturating!'
            if results_dict['max_ta_groups'] == 0:
                zero_group_error = 'Be careful! This oversaturated the TA in the minimum groups. Consider a different TA setup.'
            if results_dict['max_ta_groups'] == -1:
                zero_group_error = 'This object is too faint to reach the required TA SNR in this filter. Consider a different TA setup.'
                results_dict['min_sat_ta'] = 0
                results_dict['t_duration_ta_max'] = 0
                results_dict['max_sat_ta'] = 0
                results_dict['t_duration_ta_max'] = 0
            if results_dict['max_sat_prediction'] > results_dict['sat_max']:
                one_group_error = 'Hold up! You chose to input your own groups, and you have oversaturated the detector! Proceed with caution!'
            # Do some formatting for a prettier end product
            results_dict['filt'] = results_dict['filt'].upper()
            results_dict['filt_ta'] = results_dict['filt_ta'].upper()
            results_dict['band'] = results_dict['band'].upper()
            results_dict['mod'] = results_dict['mod'].upper()
            if results_dict['ins'] == 'niriss':
                if results_dict['subarray_ta'] == 'nrm':
                    results_dict['subarray_ta'] = 'SUBTASOSS -- BRIGHT'
                else:
                    results_dict['subarray_ta'] = 'SUBTASOSS -- FAINT'
            results_dict['subarray'] = results_dict['subarray'].upper()
            results_dict['subarray_ta'] = results_dict['subarray_ta'].upper()

            form_dict = {'miri': 'MIRI', 'nircam': 'NIRCam', 'nirspec': 'NIRSpec', 'niriss': 'NIRISS'}
            results_dict['ins'] = form_dict[results_dict['ins']]

            return render_template('groups_integrations_results.html',
                                   results_dict=results_dict,
                                   one_group_error=one_group_error,
                                   zero_group_error=zero_group_error)

        else:
            err = results
            return render_template('groups_integrations_error.html', err=err)

    except IOError:
        err = 'One of you numbers is NOT a number! Please try again!'
    except Exception as e:
        err = 'This is not an error we anticipated, but the error caught was : ' + str(e)
        return render_template('groups_integrations_error.html', err=err)


@app_exoctk.route('/contam_visibility', methods=['GET', 'POST'])
def contam_visibility():
    """The contamination and visibility form page"""

    # Log the form inputs
    try:
        log_exoctk.log_form_input(request.form, 'contam_visibility', DB)
    except:
        pass

    contamVars = {}
    if request.method == 'POST':
        tname = request.form['targetname']
        contamVars['tname'] = tname
        contamVars['ra'], contamVars['dec'] = request.form['ra'], request.form['dec']
        contamVars['PAmax'] = request.form['pamax']
        contamVars['PAmin'] = request.form['pamin']
        contamVars['inst'] = request.form['inst'].split()[0]

        if request.form['bininfo'] != '':
            contamVars['binComp'] = list(map(float, request.form['bininfo'].split(', ')))
        else:
            contamVars['binComp'] = request.form['bininfo']

        radec = ', '.join([contamVars['ra'], contamVars['dec']])

        if contamVars['PAmax'] == '':
            contamVars['PAmax'] = 359
        if contamVars['PAmin'] == '':
            contamVars['PAmin'] = 0

        if request.form['submit'] == 'Resolve Target':
            contamVars['ra'], contamVars['dec'] = resolve.resolve_target(tname)

            return render_template('contam_visibility.html', contamVars=contamVars)

        else:

            try:

                contamVars['visPA'] = True

                # Make plot
                TOOLS = 'crosshair, reset, hover, save'
                fig = figure(tools=TOOLS, plot_width=800, plot_height=400,
                             x_axis_type='datetime',
                             title=contamVars['tname'] or radec)
                pG, pB, dates, vis_plot, table = vpa.using_gtvt(contamVars['ra'],
                                                         contamVars['dec'],
                                                         contamVars['inst'],
                                                         )
                fh = StringIO()
                table.write(fh, format='ascii')
                visib_table = fh.getvalue()

                # Format x axis
                day0 = datetime.date(2019, 6, 1)
                dtm = datetime.timedelta(days=367)
                #vis_plot.x_range = Range1d(day0, day0 + dtm)

                # Get scripts
                vis_js = INLINE.render_js()
                vis_css = INLINE.render_css()
                vis_script, vis_div = components(vis_plot)

                if request.form['submit'] == 'Calculate Visibility and Contamination':

                    contamVars['contam'] = True

                    # Make field simulation
                    contam_cube = fs.sossFieldSim(contamVars['ra'],
                                                  contamVars['dec'],
                                                  binComp=contamVars['binComp'])
                    contam_plot = cf.contam(contam_cube,
                                            contamVars['tname'] or radec,
                                            paRange=[int(contamVars['PAmin']),
                                                     int(contamVars['PAmax'])],
                                            badPA=pB, fig='bokeh')

                    # Get scripts
                    contam_js = INLINE.render_js()
                    contam_css = INLINE.render_css()
                    contam_script, contam_div = components(contam_plot)

                else:

                    contamVars['contam'] = False
                    contam_script = contam_div = contam_js = contam_css = ''

                return render_template('contam_visibility_results.html',
                                       contamVars=contamVars, vis_plot=vis_div,
                                       vis_table=visib_table,
                                       vis_script=vis_script, vis_js=vis_js,
                                       vis_css=vis_css, contam_plot=contam_div,
                                       contam_script=contam_script,
                                       contam_js=contam_js,
                                       contam_css=contam_css)

            except Exception as e:
                err = 'The following error occurred: ' + str(e)
                return render_template('groups_integrations_error.html', err=err)

    return render_template('contam_visibility.html', contamVars=contamVars)

@app_exoctk.route('/visib_result', methods=['POST'])
def save_visib_result():
    """Save the results of the Visibility Only calculation"""

    visib_table = flask.request.form['data_file']
    return flask.Response(visib_table, mimetype="text/dat",
                          headers={"Content-disposition": "attachment; filename=visibility.txt"})

@app_exoctk.route('/download', methods=['POST'])
def exoctk_savefile():
    """Save results to file"""

    file_as_string = eval(request.form['file_as_string'])

    response = make_response(file_as_string)
    response.headers["Content-type"] = 'text; charset=utf-8'
    response.headers["Content-Disposition"] = "attachment; filename=ExoXTK_results.txt"
    return response


def _param_fort_validation(args):
    """Validates the input parameters for the forward models"""

    temp = args.get('ptemp', 1000)
    chem = args.get('pchem', 'noTiO')
    cloud = args.get('cloud', '0')
    pmass = args.get('pmass', '1.5')
    m_unit = args.get('m_unit', 'M_jup')
    reference_radius = args.get('refrad', 1)
    r_unit = args.get('r_unit', 'R_jup')
    rstar = args.get('rstar', 1)
    rstar_unit = args.get('rstar_unit', 'R_sun')
    
    input_args = {'temp': temp, 'chem': chem, 'cloud': cloud, 'pmass': pmass, 
                  'm_unit': m_unit, 'reference_radius': reference_radius, 
                  'r_unit': r_unit, 'rstar': rstar, 'rstar_unit': rstar_unit}

    return input_args 


@app_exoctk.route('/fortney', methods=['GET', 'POST'])
def fortney():
    """
    Pull up Forntey Grid plot the results and download
    """

    # Grab the inputs arguments from the URL
    args = flask.request.args

    input_args = _param_fort_validation(args)
    
    fig, fh = fortney_grid(input_args)
    table_string = fh.getvalue()

    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    script, div = components(fig)

    html = flask.render_template('fortney.html',
                                 plot_script=script,
                                 plot_div=div,
                                 js_resources=js_resources,
                                 css_resources=css_resources,
                                 temp=list(map(str, header.temp.unique())),
                                 table_string=table_string
                                 )
    return encode_utf8(html)

    
@app_exoctk.route('/generic', methods=['GET', 'POST'])
def generic():
    """
    Pull up Generic Grid plot the results and download
    """

    # Grab the inputs arguments from the URL
    args = dict(flask.request.args)
    for key in args:
        args[key] = args[key][0]
    
    fig, fh = generic_grid(args)

    # Write table string
    table_string = fh.getvalue()
    
    # Web-ify bokeh plot
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    script, div = components(fig)

    html = flask.render_template('generic.html',
                                 inputs= inputs,
                                 closest_match = closest_match,
                                 error_message=error_message,
                                 table_string=table_string,
                                 plot_script=script,
                                 plot_div=div,
                                 js_resources=js_resources,
                                 css_resources=css_resources,
                                 )
    return encode_utf8(html)


@app_exoctk.route('/fortney_result', methods=['POST'])
def save_fortney_result():
    """SAve the results of the Fortney grid"""

    table_string = flask.request.form['data_file']
    return flask.Response(table_string, mimetype="text/dat",
                          headers={"Content-disposition": "attachment; filename=fortney.dat"})


@app_exoctk.route('/generic_result', methods=['POST'])
def save_generic_result():
    """Save the results of the generic grid"""

    table_string = flask.request.form['data_file']
    return flask.Response(table_string, mimetype="text/dat",
                          headers={"Content-disposition": "attachment; filename=generic.dat"})


@app_exoctk.route('/zip_data_download')
def zip_data_download():
    """Download the exoctk data."""

    return send_file(resource_filename('exoctk', 'data/exoctk_data.zip'), mimetype="text/json",
                     attachment_filename='exoctk_data.zip',
                     as_attachment=True)


def check_auth(username, password):
    """This function is called to check if a username password combination is
    valid

    Parameters
    ----------
    username: str
        The username
    password: str
        The password
    """

    return username == 'admin' and password == 'secret'


def authenticate():
    """Sends a 401 response that enables basic auth"""

    return Response('Could not verify your access level for that URL.\n'
                    'You have to login with proper credentials', 401,
                    {'WWW-Authenticate': 'Basic realm="Login Required"'})


def requires_auth(page):
    """Requires authentication for a page before loading

    Parameters
    ----------
    page: function
        The function that sets a route

    Returns
    -------
    function
        The decorated route
    """

    @wraps(page)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if not auth or not check_auth(auth.username, auth.password):
            return authenticate()
        return page(*args, **kwargs)
    return decorated


@app_exoctk.route('/admin')
@requires_auth
def secret_page():
    """Shhhhh! This is a secret page of admin stuff"""

    tables = [i[0] for i in DB.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()]
    print(tables)

    log_tables = []
    for table in tables:

        try:
            data = log_exoctk.view_log(DB, table)

            # Add the results to the lists
            html_table = '\n'.join(data.pformat(max_width=500, html=True)).replace('<table', '<table id="myTable" class="table table-striped table-hover"')

        except:
            html_table = '<p>No data to display</p>'

        # Add the table title
        header = '<h3>{}</h3>'.format(table)
        html_table = header + html_table

        log_tables.append(html_table)

    return render_template('admin_page.html', tables=log_tables)


if __name__ == '__main__':
    # os.chmod('/internal/data1/app_data/.astropy/cache/', 777)
    port = int(os.environ.get('PORT', 5000))
    app_exoctk.run(host='0.0.0.0', port=port, debug=True)

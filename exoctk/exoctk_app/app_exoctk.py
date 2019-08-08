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
import urllib
from wtforms.validators import InputRequired, NumberRange
from wtforms import DecimalField

from exoctk.contam_visibility import resolve
from exoctk.contam_visibility import visibilityPA as vpa
from exoctk.contam_visibility import sossFieldSim as fs
from exoctk.contam_visibility import sossContamFig as cf
from exoctk.forward_models.forward_models import fortney_grid, generic_grid
from exoctk.groups_integrations.groups_integrations import perform_calculation
from exoctk.limb_darkening import limb_darkening_fit as lf
from exoctk.utils import find_closest, filter_table, get_target_data, get_canonical_name, FORTGRID_DIR, EXOCTKLOG_DIR, GENERICGRID_DIR
from exoctk.modelgrid import ModelGrid

import log_exoctk
from svo_filters import svo

# FLASK SET UP
app_exoctk = Flask(__name__)

# define the cache config keys, remember that it can be done in a settings file
app_exoctk.config['CACHE_TYPE'] = 'null'
app_exoctk.config['SECRET_KEY'] = 'Thisisasecret!'


# Load the database to log all form submissions
if EXOCTKLOG_DIR is None:
    dbpath = ':memory:'
else:
    dbpath = os.path.realpath(os.path.join(EXOCTKLOG_DIR, 'exoctk_log.db'))
    if not os.path.isfile(dbpath):
        log_exoctk.create_db(dbpath)
try:
    DB = log_exoctk.load_db(dbpath)
except IOError:
    DB = None


# Redirect to the index
@app_exoctk.route('/')
@app_exoctk.route('/index')
def index():
    """The Index page"""
    return render_template('index.html')

@app_exoctk.route('/limb_darkening', methods=['GET', 'POST'])
def limb_darkening():
    """The limb darkening form page. """
    # Load default form
    form = fv.LimbDarkeningForm()

    # Reload page with stellar data from ExoMAST
    if form.resolve_submit.data:

        if form.targname.data.strip() != '':

            try:

                # Resolve the target in exoMAST
                form.targname.data = get_canonical_name(form.targname.data)
                data, target_url = get_target_data(form.targname.data)

                # Update the form data
                form.feh.data = data.get('Fe/H')
                form.teff.data = data.get('Teff')
                form.logg.data = data.get('stellar_gravity')
                form.target_url.data = str(target_url)

            except:
                form.target_url.data = ''
                form.targname.errors = ["Sorry, could not resolve '{}' in exoMAST.".format(form.targname.data)]

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
        star_params = [float(form.teff.data), float(form.logg.data), float(form.feh.data)]

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
            ld.calculate(*star_params, prof, mu_min=float(form.mu_min.data), bandpass=bandpass)

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
    # Load default form
    form = fv.GroupsIntsForm()

    if request.method == 'GET':       

        # http://0.0.0.0:5000/groups_integrations?k_mag=8.131&transit_duration=0.09089&target=WASP-18+b
        target_name = request.args.get('target')
        form.targname.data = target_name

        k_mag = request.args.get('k_mag')
        form.kmag.data = k_mag

        # According to Kevin the obs_dur = 3*trans_dur+1 hours
        # transit_dur is in days from exomast, convert first.
        try:
            trans_dur = float(request.args.get('transit_duration'))
            trans_dur *= u.day.to(u.hour)
            obs_dur = 3*trans_dur + 1
            form.obs_duration.data = obs_dur
        except TypeError:
            trans_dur = request.args.get('transit_duration')
            if trans_dur == None:
                pass
            else:
                err = 'The Transit Duration from ExoMAST experienced some issues. Try a different spelling or source.'
                return render_template('groups_integrations_error.html', err=err)
        return render_template('groups_integrations.html', form=form, sat_data=sat_data)
        
    # Reload page with stellar data from ExoMAST
    if form.resolve_submit.data:

        if form.targname.data.strip() != '':

            # Resolve the target in exoMAST
            try:
                form.targname.data = get_canonical_name(form.targname.data)
                data, url = get_target_data(form.targname.data)

                # Update the Kmag
                kmag = data.get('Kmag')

                # Transit duration in exomast is in days, need it in hours
                if form.time_unit.data == 'day':
                    trans_dur = data.get('transit_duration')
                    obs_dur = 3*trans_dur + (1/24.) 
                else:
                    trans_dur = data.get('transit_duration')
                    trans_dur *= u.Unit('day').to('hour')
                    obs_dur = 3*trans_dur + 1
                    

                # Model guess
                logg_targ = data.get('stellar_gravity') or 4.5
                teff_targ = data.get('Teff') or 5500
                arr = np.array([tuple(i[1].split()) for i in form.mod.choices], dtype=[('spt', 'O'), ('teff', '>f4'), ('logg', '>f4')])
                mod_table = at.Table(arr)

                # If high logg, remove giants from guess list
                if logg_targ < 4:
                    mod_table = filter_table(mod_table, logg=">=4")
                teff = min(arr['teff'], key=lambda x:abs(x-teff_targ))
                mod_table.add_column(at.Column(np.array([i[0] for i in form.mod.choices]), name='value'))
                mod_table = filter_table(mod_table, teff="<={}".format(teff))
                mod_table.sort(['teff', 'logg'])

                # Set the form values
                form.mod.data = mod_table[-1]['value']
                form.kmag.data = kmag
                form.obs_duration.data = obs_dur
                form.target_url.data = url

            except:
                form.target_url.data = ''
                form.targname.errors = ["Sorry, could not resolve '{}' in exoMAST.".format(form.targname.data)]

        # Send it back to the main page
        return render_template('groups_integrations.html', form=form, sat_data=sat_data)

    if form.validate_on_submit() and form.calculate_submit.data:

        # Get the form data
        ins = form.ins.data
        params = {'ins': ins,
                  'mag': form.kmag.data,
                  'obs_time': form.obs_duration.data,
                  'sat_max': form.sat_max.data,
                  'sat_mode': form.sat_mode.data,
                  'time_unit': form.time_unit.data,
                  'band': 'K',
                  'mod': form.mod.data,
                  'filt'.format(ins): getattr(form, '{}_filt'.format(ins)).data,
                  'subarray'.format(ins): getattr(form, '{}_subarray'.format(ins)).data,
                  'filt_ta'.format(ins): getattr(form, '{}_filt_ta'.format(ins)).data,
                  'subarray_ta'.format(ins): getattr(form, '{}_subarray_ta'.format(ins)).data}

        # Get ngroups
        params['n_group'] = 'optimize' if form.n_group.data == 0 else int(form.n_group.data)

        # Also get the data path in there
        params['infile'] = resource_filename('exoctk', 'data/groups_integrations/groups_integrations_input_data.json')

        # Convert the obs_time to hours
        if params['time_unit'] == 'day':
            params['obs_time'] = params['obs_time']*24
            params['time_unit'] = 'hours'

        # Run the calculation
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
                one_group_error = 'This many groups will oversaturate the detector! Proceed with caution!'
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

    return render_template('groups_integrations.html', form=form, sat_data=sat_data)


@app_exoctk.route('/contam_visibility', methods=['GET', 'POST'])
def contam_visibility():
    """The contamination and visibility form page"""
    # Load default form
    form = fv.ContamVisForm()
    form.calculate_contam_submit.disabled = False

    if request.method == 'GET':       

        # http://0.0.0.0:5000/contam_visibility?ra=24.354208334287005&dec=-45.677930555343636&target=WASP-18%20b       
        target_name = request.args.get('target')
        form.targname.data = target_name

        ra = request.args.get('ra')
        form.ra.data = ra
        
        dec = request.args.get('dec')
        form.dec.data = dec
                
        return render_template('contam_visibility.html', form=form)

    # Reload page with stellar data from ExoMAST
    if form.resolve_submit.data:

        if form.targname.data.strip() != '':

            # Resolve the target in exoMAST
            try:
                form.targname.data = get_canonical_name(form.targname.data)
                data, url = get_target_data(form.targname.data)

                # Update the coordinates
                ra = data.get('RA')
                dec = data.get('DEC')

                # Set the form values
                form.ra.data = ra
                form.dec.data = dec
                form.target_url.data = url

            except:
                form.target_url.data = ''
                form.targname.errors = ["Sorry, could not resolve '{}' in exoMAST.".format(form.targname.data)]

        # Send it back to the main page
        return render_template('contam_visibility.html', form=form)

    # Reload page with appropriate mode data
    if form.mode_submit.data:

        # Update the button
        if form.inst.data != 'NIRISS':
            form.calculate_contam_submit.disabled = True
        else:
            form.calculate_contam_submit.disabled = False

        # Send it back to the main page
        return render_template('contam_visibility.html', form=form)

    if form.validate_on_submit() and (form.calculate_submit.data or form.calculate_contam_submit.data):

        try:

            # Log the form inputs
            try:
                log_exoctk.log_form_input(request.form, 'contam_visibility', DB)
            except:
                pass



            # Make plot
            title = form.targname.data or ', '.join([form.ra.data, form.dec.data])
            pG, pB, dates, vis_plot, table = vpa.using_gtvt(str(form.ra.data), str(form.dec.data), form.inst.data.split(' ')[0])

            # Make output table
            vers = '0.3'
            today = datetime.datetime.now()
            fh = StringIO()
            fh.write('# Hi! This is your Visibility output file for... \n')
            fh.write('# Target: {} \n'.format(form.targname.data))
            fh.write('# Instrument: {} \n'.format(form.inst.data))
            fh.write('# \n')
            fh.write('# This file was generated using ExoCTK v{} on {} \n'.format(vers, today))
            fh.write('# Visit our GitHub: https://github.com/ExoCTK/exoctk \n')
            fh.write('# \n')
            table.write(fh, format='csv', delimiter=',')
            visib_table = fh.getvalue()

            # Format x axis
            day0 = datetime.date(2019, 6, 1)
            dtm = datetime.timedelta(days=367)
            #vis_plot.x_range = Range1d(day0, day0 + dtm)


            # TODO: Fix this so bad PAs are included
            pB = []

            # Make plot
            TOOLS = 'crosshair, reset, hover, save'
            fig = figure(tools=TOOLS, plot_width=800, plot_height=400, x_axis_type='datetime', title=title)
            fh = StringIO()
            table.write(fh, format='ascii')
            visib_table = fh.getvalue()

            # Format x axis
            day0 = datetime.date(2019, 6, 1)
            dtm = datetime.timedelta(days=367)

            # Get scripts
            vis_js = INLINE.render_js()
            vis_css = INLINE.render_css()
            vis_script, vis_div = components(vis_plot)

            # Contamination plot too
            if form.calculate_contam_submit.data:

                # Make field simulation
                contam_cube = fs.sossFieldSim(form.ra.data, form.dec.data, binComp=form.companion.data)
                contam_plot = cf.contam(contam_cube, title, paRange=[int(form.pa_min.data), int(form.pa_max.data)], badPA=pB, fig='bokeh')

                # Get scripts
                contam_js = INLINE.render_js()
                contam_css = INLINE.render_css()
                contam_script, contam_div = components(contam_plot)

            else:

                contam_script = contam_div = contam_js = contam_css = ''

            return render_template('contam_visibility_results.html',
                                   form=form, vis_plot=vis_div,
                                   vis_table=visib_table,
                                   vis_script=vis_script, vis_js=vis_js,
                                   vis_css=vis_css, contam_plot=contam_div,
                                   contam_script=contam_script,
                                   contam_js=contam_js,
                                   contam_css=contam_css)

        except IOError:#Exception as e:
            err = 'The following error occurred: ' + str(e)
            return render_template('groups_integrations_error.html', err=err)

    return render_template('contam_visibility.html', form=form)


@app_exoctk.route('/visib_result', methods=['POST'])
def save_visib_result():
    """Save the results of the Visibility Only calculation"""

    visib_table = flask.request.form['data_file']
    targname = flask.request.form['targetname']
    instname = flask.request.form['instrumentname']

    return flask.Response(visib_table, mimetype="text/dat",
                          headers={"Content-disposition": "attachment; filename={}_{}_visibility.csv".format(targname, instname)})


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
    fig, fh, temp_out = fortney_grid(input_args)

    table_string = fh.getvalue()

    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    script, div = components(fig)

    html = flask.render_template('fortney.html',
                                 plot_script=script,
                                 plot_div=div,
                                 js_resources=js_resources,
                                 css_resources=css_resources,
                                 temp=temp_out,
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
    fig, fh, closest_match, error_message = generic_grid(args)

    # Write table string
    table_string = fh.getvalue()
    
    # Web-ify bokeh plot
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    script, div = components(fig)

    html = flask.render_template('generic.html',
                                 inputs= args,
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


@app_exoctk.route('/groups_integrations_download')
def groups_integrations_download():
    """Download the groups and integrations calculator data"""

    return send_file(resource_filename('exoctk', 'data/groups_integrations/groups_integrations_input_data.json'), mimetype="text/json",
                     attachment_filename='groups_integrations_input_data.json',
                     as_attachment=True)


@app_exoctk.route('/fortney_download')
def fortney_download():
    """Download the fortney grid data"""

    fortney_data = os.path.join(FORTGRID_DIR, 'fortney_grid.db')
    return send_file(fortney_data, attachment_filename='fortney_grid.db',
                     as_attachment=True)

@app_exoctk.route('/zip_data_download')
def zip_data_download():
    """Download the zipped ExoCTK data"""
    
    return send_file(resource_filename('exoctk', 'data/exoctk_data.zip'), mimetype='application/zip',
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

@app_exoctk.route('/lightcurve_fitting')
def lightcurve_fitting():
    """A landing page for the lightcurve_fitting tool"""

    return render_template('lightcurve_fitting.html')


@app_exoctk.route('/atmospheric_retrievals')
def atmospheric_retrievals():
    """A landing page for the atmospheric_retrievals tools"""

    return render_template('atmospheric_retrievals.html')


if __name__ == '__main__':
    # os.chmod('/internal/data1/app_data/.astropy/cache/', 777)
    port = int(os.environ.get('PORT', 5000))
    app_exoctk.run(host='0.0.0.0', port=port)

from functools import wraps
import io
import json
import os
from pkg_resources import resource_filename

from astropy.coordinates import SkyCoord
import astropy.table as at
from astropy.time import Time
import astropy.units as u
from bokeh.embed import components
from bokeh.resources import INLINE
import flask
from flask import Flask, make_response, render_template, Response, request, send_file
import form_validation as fv
import numpy as np

from exoctk import log_exoctk
from exoctk.contam_visibility.new_vis_plot import build_visibility_plot, get_exoplanet_positions
from exoctk.contam_visibility import field_simulator as fs
from exoctk.contam_visibility import contamination_figure as cf
from exoctk.contam_visibility.miniTools import contamVerify
from exoctk.forward_models.forward_models import fortney_grid, generic_grid
from exoctk.groups_integrations.groups_integrations import perform_calculation
from exoctk.limb_darkening import limb_darkening_fit as lf
from exoctk.limb_darkening import spam
from exoctk.modelgrid import ModelGrid
from exoctk.phase_constraint_overlap.phase_constraint_overlap import phase_overlap_constraint, calculate_pre_duration
from exoctk.throughputs import Throughput
from exoctk.utils import filter_table, get_env_variables, get_target_data, get_canonical_name

# FLASK SET UP
app_exoctk = Flask(__name__)

# define the cache config keys, remember that it can be done in a settings file
app_exoctk.config['CACHE_TYPE'] = 'null'
app_exoctk.config['SECRET_KEY'] = 'Thisisasecret!'

# Load the database to log all form submissions
if get_env_variables()['exoctklog_dir'] is None:
    DBPATH = ':memory:'
else:
    DBPATH = os.path.realpath(os.path.join(get_env_variables()['exoctklog_dir'], 'exoctk_log.db'))
    if not os.path.isfile(DBPATH):
        log_exoctk.create_db(DBPATH)
try:
    DB = log_exoctk.load_db(DBPATH)
except IOError:
    DB = None


def _param_fort_validation(args):
    """Validates the input parameters for the forward models

    Returns
    -------
    input_args : dict
        Dictionary with the input parameters for the forward models.
    """

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


def check_auth(username, password):
    """This function is called to check if a username password
    combination is valid

    Parameters
    ----------
    username: str
        The username
    password: str
        The password
    """

    return username == 'admin' and password == 'secret'


@app_exoctk.route('/download', methods=['POST'])
def exoctk_savefile():
    """Save results to file

    Returns
    -------
    ``flask.make_response`` obj
        Returns response including results in txt form.
    """

    file_as_string = eval(request.form['file_as_string'])

    response = make_response(file_as_string)
    response.headers["Content-type"] = 'text; charset=utf-8'
    response.headers["Content-Disposition"] = "attachment; filename=ExoCTK_results.txt"

    return response


@app_exoctk.route('/fortney', methods=['GET', 'POST'])
def fortney():
    """
    Pull up Forntey Grid plot the results and download

    Returns
    -------
    html : ``flask.render_template`` obj
        The rendered template for the Fortney results.
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
                                 temp=sorted(temp_out, key=float),
                                 table_string=table_string
                                 )

    # Log the form inputs
    if len(args) > 0:
        log_exoctk.log_form_input(args, 'fortney', DB)

    return html


@app_exoctk.route('/fortney_download')
def fortney_download():
    """Download the fortney grid data"""

    fortney_data = os.path.join(get_env_variables()['fortgrid_dir'], 'fortney_grid.db')
    return send_file(fortney_data, attachment_filename='fortney_grid.db', as_attachment=True)


@app_exoctk.route('/generic', methods=['GET', 'POST'])
def generic():
    """
    Pull up Generic Grid plot the results and download

    Returns
    -------
    html: ``flask.render_template`` obj
        The rendered template for the generic grid page.
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
                                 inputs=args,
                                 closest_match=closest_match,
                                 error_message=error_message,
                                 table_string=table_string,
                                 plot_script=script,
                                 plot_div=div,
                                 js_resources=js_resources,
                                 css_resources=css_resources,
                                 )

    # Log the form inputs
    if len(args) > 0:
        log_exoctk.log_form_input(args, 'generic', DB)

    return html


@app_exoctk.route('/groups_integrations', methods=['GET', 'POST'])
def groups_integrations():
    """The groups and integrations calculator form page

    Returns
    -------
    ``flask.render_template`` obj
        The rendered template for the Groups & Integrations calculator
        page.
    """

    # Print out pandeia sat values
    with open(resource_filename('exoctk', 'data/groups_integrations/groups_integrations_input_data.json')) as f:
        sat_data = json.load(f)['fullwell']

    # Load default form
    form = fv.GroupsIntsForm()

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
                    obs_dur = 3 * trans_dur + (1 / 24.)
                else:
                    trans_dur = data.get('transit_duration')
                    trans_dur *= u.Unit('day').to('hour')
                    obs_dur = 3 * trans_dur + 1

                # Model guess
                logg_targ = data.get('stellar_gravity') or 4.5
                teff_targ = data.get('Teff') or 5500
                arr = np.array([tuple(i[1].split()) for i in form.mod.choices], dtype=[('spt', 'O'), ('teff', '>f4'), ('logg', '>f4')])
                mod_table = at.Table(arr)

                # If high logg, remove giants from guess list
                if logg_targ < 4:
                    mod_table = filter_table(mod_table, logg=">=4")
                teff = min(arr['teff'], key=lambda x: abs(x - teff_targ))
                mod_table.add_column(at.Column(np.array([i[0] for i in form.mod.choices]), name='value'))
                mod_table = filter_table(mod_table, teff="<={}".format(teff))
                mod_table.sort(['teff', 'logg'])

                # Set the form values
                form.mod.data = mod_table[-1]['value']
                form.kmag.data = kmag
                form.obs_duration.data = obs_dur
                form.target_url.data = url

            except Exception:
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
            params['obs_time'] = params['obs_time'] * 24
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
                results_dict['min_saturation_ta'] = 0
                results_dict['duration_time_ta_max'] = 0
                results_dict['max_saturation_ta'] = 0
                results_dict['duration_time_ta_max'] = 0
            if results_dict['max_saturation_prediction'] > results_dict['sat_max']:
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

            # Log the successful form inputs
            params['kmag'] = form.kmag.data
            params['targname'] = form.targname.data
            params['num_groups'] = form.n_group.data
            log_exoctk.log_form_input(params, 'groups_integrations', DB)

            return render_template('groups_integrations_results.html',
                                   results_dict=results_dict,
                                   one_group_error=one_group_error,
                                   zero_group_error=zero_group_error)

        else:
            err = results
            return render_template('groups_integrations_error.html', err=err)

    return render_template('groups_integrations.html', form=form, sat_data=sat_data)


@app_exoctk.route('/pa_contam', methods=['GET', 'POST'])
def pa_contam():
    """The contamination and visibility form page

    Returns
    -------
    ``flask.render_template`` obj
        The rendered template for the contamination and visibility page.

    """
    return render_template('pa_contam.html')


@app_exoctk.route('/contam_visibility', methods=['GET', 'POST'])
def contam_visibility():
    """The contamination and visibility form page
    
    Returns
    -------
    ``flask.render_template`` obj
        The rendered template for the contamination and visibility page.

    """
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
                ra_deg = data.get('RA')
                dec_deg = data.get('DEC')

                # Set the form values
                form.ra.data = ra_deg
                form.dec.data = dec_deg
                form.target_url.data = url

            except Exception:
                form.target_url.data = ''
                form.targname.errors = ["Sorry, could not resolve '{}' in exoMAST.".format(form.targname.data)]

        # Send it back to the main page
        return render_template('contam_visibility.html', form=form)

    # Reload page with appropriate mode data
    if form.mode_submit.data:

        # Update the button
        if form.inst.data == 'NIRSpec':
            form.calculate_contam_submit.disabled = True
        else:
            form.calculate_contam_submit.disabled = False

        # Send it back to the main page
        return render_template('contam_visibility.html', form=form)

    if form.validate_on_submit() and (form.calculate_submit.data or form.calculate_contam_submit.data):

        if form.inst.data == "NIRSpec":
            instrument = form.inst.data
        else:
            instrument = fs.APERTURES[form.inst.data]['inst']

        try:
            # Log the form inputs
            log_exoctk.log_form_input(request.form, 'contam_visibility', DB)

            # Make plot
            title = form.targname.data or ', '.join([str(form.ra.data), str(form.dec.data)])
            vis_plot = build_visibility_plot(str(title), instrument, str(form.ra.data), str(form.dec.data))
            table = get_exoplanet_positions(str(form.ra.data), str(form.dec.data))

            # Make output table
            vis_table = table.to_csv()

            # Get scripts
            vis_js = INLINE.render_js()
            vis_css = INLINE.render_css()
            vis_script, vis_div = components(vis_plot)

            # Contamination plot too
            if form.calculate_contam_submit.data:

                # Get RA and Dec in degrees
                ra_deg, dec_deg = float(form.ra.data), float(form.dec.data)

                # Add companion
                try:
                    comp_teff = float(form.teff.data)
                except TypeError:
                    comp_teff = None
                try:
                    comp_mag = float(form.delta_mag.data)
                except TypeError:
                    comp_mag = None
                try:
                    comp_dist = float(form.dist.data)
                except TypeError:
                    comp_dist = None
                try:
                    comp_pa = float(form.pa.data)
                except TypeError:
                    comp_pa = None

                # Get PA value
                pa_val = float(form.v3pa.data)
                if pa_val == -1:

                    # Add a companion
                    companion = None
                    if comp_teff is not None and comp_mag is not None and comp_dist is not None and comp_pa is not None:
                        companion = {'name': 'Companion', 'ra': ra_deg, 'dec': dec_deg, 'teff': comp_teff, 'delta_mag': comp_mag, 'dist': comp_dist, 'pa': comp_pa}

                    # Make field simulation
                    targframe, starcube, results = fs.field_simulation(ra_deg, dec_deg, form.inst.data, target_date=form.epoch.data, binComp=companion, plot=False, multi=False)

                    # Make the plot
                    # contam_plot = fs.contam_slider_plot(results)

                    # Get bad PA list from missing angles between 0 and 360
                    badPAs = [j for j in np.arange(0, 360) if j not in [i['pa'] for i in results]]

                    # Make old contam plot
                    starCube = np.zeros((362, 2048, 96 if form.inst.data=='NIS_SUBSTRIP96' else 256))
                    starCube[0, :, :] = (targframe[0]).T[::-1, ::-1]
                    starCube[1, :, :] = (targframe[1]).T[::-1, ::-1]
                    starCube[2:, :, :] = starcube.swapaxes(1, 2)[:, ::-1, ::-1]
                    contam_plot = cf.contam(starCube, form.inst.data, targetName=form.targname.data, badPAs=badPAs)

                else:

                    # Get stars
                    stars = fs.find_sources(ra_deg, dec_deg, target_date=form.epoch.data, verbose=False)

                    # Add companion
                    if comp_teff is not None and comp_mag is not None and comp_dist is not None and comp_pa is not None:
                        stars = fs.add_source(stars, 'Companion', ra, dec, teff=comp_teff, delta_mag=comp_mag, dist=comp_dist, pa=comp_pa, type='STAR')

                    # Calculate contam
                    result, contam_plot = fs.calc_v3pa(pa_val, stars, form.inst.data, plot=True, verbose=False)

                # Get scripts
                contam_js = INLINE.render_js()
                contam_css = INLINE.render_css()
                contam_script, contam_div = components(contam_plot)

            else:

                contam_script = contam_div = contam_js = contam_css = pa_val = ''

            return render_template('contam_visibility_results.html',
                                   form=form, vis_plot=vis_div,
                                   vis_table=vis_table,
                                   vis_script=vis_script, vis_js=vis_js,
                                   vis_css=vis_css, contam_plot=contam_div,
                                   contam_script=contam_script,
                                   contam_js=contam_js,
                                   contam_css=contam_css, pa_val=pa_val, epoch=form.epoch.data)

        except Exception as e:
            err = 'The following error occurred: ' + str(e)
            return render_template('groups_integrations_error.html', err=err)

    return render_template('contam_visibility.html', form=form)


# Redirect to the index
@app_exoctk.route('/')
@app_exoctk.route('/index')
def index():
    """Returns the rendered index page

    Returns
    -------
    ``flask.render_template`` obj
        The rendered template for the index page.
    """
    return render_template('index.html')


@app_exoctk.route('/limb_darkening', methods=['GET', 'POST'])
def limb_darkening():
    """Returns the rendered limb darkening form page.

    Returns
    -------
    ``flask.render_template`` obj
        The rendered template for the limb-darkening page.
    """
    # Load default form
    form = fv.LimbDarkeningForm()

    # Planet properties
    planet_properties = ['transit_duration', 'orbital_period', 'rp_rs', 'a_rs', 'inclination', 'eccentricity', 'omega']

    def empty_fields(form):
        form.transit_duration.data = form.transit_duration.data or ''
        form.orbital_period.data = form.orbital_period.data or ''
        form.rp_rs.data = form.rp_rs.data or ''
        form.a_rs.data = form.a_rs.data or ''
        form.inclination.data = form.inclination.data or ''
        form.eccentricity.data = form.eccentricity.data or ''
        form.omega.data = form.omega.data or ''

        return form

    # Reload page with stellar data from ExoMAST
    if form.resolve_submit.data:

        if form.targname.data.strip() != '':

            try:

                # Resolve the target in exoMAST
                form.targname.data = get_canonical_name(form.targname.data)
                data, target_url = get_target_data(form.targname.data)

                # Update the star data
                form.feh.data = data.get('Fe/H')
                form.teff.data = data.get('Teff')
                form.logg.data = data.get('stellar_gravity')
                form.target_url.data = str(target_url)

                # Update the planet data
                form.transit_duration.data = data.get('transit_duration')
                form.orbital_period.data = data.get('orbital_period')
                form.rp_rs.data = data.get('Rp/Rs')
                form.a_rs.data = data.get('a/Rs')
                form.inclination.data = data.get('inclination')
                form.eccentricity.data = data.get('eccentricity')
                form.omega.data = data.get('omega')

            except Exception:
                form.target_url.data = ''
                form.targname.errors = ["Sorry, could not resolve '{}' in exoMAST.".format(form.targname.data)]

        # Ensure planet fields are not None
        form = empty_fields(form)

        # Send it back to the main page
        return render_template('limb_darkening.html', form=form)

    # Reload page with appropriate filter data
    if form.filter_submit.data:

        kwargs = {}
        if form.bandpass.data == 'tophat':
            kwargs['n_bins'] = 1
            kwargs['pixels_per_bin'] = 100
            kwargs['wave_min'] = 1 * u.um
            kwargs['wave_max'] = 2 * u.um

        # Get the filter
        bandpass = Throughput(form.bandpass.data, **kwargs)

        # Update the form data
        form.wave_min.data = bandpass.wave_min.value
        form.wave_max.data = bandpass.wave_max.value

        # Ensure planet fields are not None
        form = empty_fields(form)

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
        setattr(form.teff.validators[1], 'min', float(teff_rng[0]))
        setattr(form.teff.validators[1], 'max', float(teff_rng[1]))
        setattr(form.teff.validators[1], 'message', 'Effective temperature must be between {} and {}'.format(*teff_rng))
        setattr(form.logg.validators[1], 'min', float(logg_rng[0]))
        setattr(form.logg.validators[1], 'max', float(logg_rng[1]))
        setattr(form.logg.validators[1], 'message', 'Surface gravity must be between {} and {}'.format(*logg_rng))
        setattr(form.feh.validators[1], 'min', float(feh_rng[0]))
        setattr(form.feh.validators[1], 'max', float(feh_rng[1]))
        setattr(form.feh.validators[1], 'message', 'Metallicity must be between {} and {}'.format(*feh_rng))

        # Ensure planet fields are not None
        form = empty_fields(form)

        # Send it back to the main page
        return render_template('limb_darkening.html', form=form)

    # Validate form and submit for results
    if form.validate_on_submit() and form.calculate_submit.data:

        # Form inputs for logging
        form_input = dict(request.form)

        # Get the stellar parameters
        star_params = [float(form.teff.data), float(form.logg.data), float(form.feh.data)]

        # Load the model grid
        model_grid = ModelGrid(form.modeldir.data, resolution=500)
        form.modeldir.data = [j for i, j in form.modeldir.choices if i == form.modeldir.data][0]

        # Grism details
        kwargs = {'n_bins': form.n_bins.data, 'wave_min': form.wave_min.data * u.um, 'wave_max': form.wave_max.data * u.um}

        # Make filter object and plot
        bandpass = Throughput(form.bandpass.data, **kwargs)
        bk_plot = bandpass.plot(draw=False)
        bk_plot.width = 580
        bk_plot.height = 280
        js_resources = INLINE.render_js()
        css_resources = INLINE.render_css()
        filt_script, filt_plot = components(bk_plot)

        # Trim the grid to nearby grid points to speed up calculation
        # full_rng = [model_grid.Teff_vals, model_grid.logg_vals, model_grid.FeH_vals]
        # trim_rng = find_closest(full_rng, star_params, n=1, values=True)

        # Calculate the coefficients for each profile
        ld = lf.LDC(model_grid)
        for prof in form.profiles.data:
            ld.calculate(*star_params, prof, mu_min=float(form.mu_min.data), bandpass=bandpass)

        # Check if spam coefficients can be calculated
        planet_data = {param: getattr(form, param).data for param in planet_properties}
        planet_data['Rp/Rs'] = planet_data['rp_rs']
        planet_data['a/Rs'] = planet_data['a_rs']
        spam_calc = all([val is not None and val != '' for key, val in planet_data.items()])
        if spam_calc:

            # Make sure non-linear profile is included for spam calculation if all planet parameters are provided
            if '4-parameter' not in form.profiles.data:
                ld.calculate(*star_params, '4-parameter', mu_min=float(form.mu_min.data), bandpass=bandpass)

            # Calculate spam coeffs
            planet_data = {key: float(val) for key, val in planet_data.items()}
            ld.spam(planet_data=planet_data)

        # Draw tabbed figure
        final = ld.plot_tabs()

        # Get HTML
        script, div = components(final)

        # Store the tables as a string
        keep_cols = ['Teff', 'logg', 'FeH', 'profile', 'filter', 'wave_min', 'wave_eff', 'wave_max', 'c1', 'e1', 'c2', 'e2', 'c3', 'e3', 'c4', 'e4']
        print_table = ld.results[[col for col in keep_cols if col in ld.results.colnames]]
        file_as_string = '\n'.join(print_table.pformat(max_lines=-1, max_width=-1))

        # Make a table for each profile with a row for each wavelength bin
        profile_tables = []
        for profile in form.profiles.data:

            # Make LaTeX for polynomials
            latex = lf.ld_profile(profile, latex=True)
            poly = '\({}\)'.format(latex).replace('*', '\cdot').replace('\e', 'e')

            # Make the table into LaTeX
            table = filter_table(ld.results, profile=profile)
            co_cols = [c for c in ld.results.colnames if (c.startswith('c') or c.startswith('e')) and len(c) == 2 and not np.all([np.isnan(i) for i in table[c]])]
            table = table[['wave_eff', 'wave_min', 'wave_max'] + co_cols]
            table.rename_column('wave_eff', '\(\lambda_\mbox{eff}\hspace{5px}(\mu m)\)')
            table.rename_column('wave_min', '\(\lambda_\mbox{min}\hspace{5px}(\mu m)\)')
            table.rename_column('wave_max', '\(\lambda_\mbox{max}\hspace{5px}(\mu m)\)')

            # Add the results to the lists
            html_table = '\n'.join(table.pformat(max_width=-1, max_lines=-1, html=True)).replace('<table', '<table id="myTable" class="table table-striped table-hover"')

            # Add the table title
            header = '<br></br><strong>{}</strong><br><p>\(I(\mu)/I(\mu=1)\) = {}</p>'.format(profile, poly)
            html_table = header + html_table

            profile_tables.append(html_table)

            # Add the profile to the form inputs
            form_input[profile] = 'true'

        # Log the successful form inputs
        log_exoctk.log_form_input(form_input, 'limb_darkening', DB)

        # Make a table for each profile with a row for each wavelength bin
        profile_spam_tables = ''
        spam_file_as_string = ''
        if ld.spam_results is not None:

            # Store SPAM tables as string
            keep_cols = ['Teff', 'logg', 'FeH', 'profile', 'filter', 'wave_min', 'wave_eff', 'wave_max', 'c1', 'c2']
            print_spam_table = ld.spam_results[[col for col in keep_cols if col in ld.spam_results.colnames]]
            spam_file_as_string = '\n'.join(print_spam_table.pformat(max_lines=-1, max_width=-1))
            profile_spam_tables = []
            for profile in list(np.unique(ld.spam_results['profile'])):

                # Make LaTeX for polynomials
                latex = lf.ld_profile(profile, latex=True)
                poly = '\({}\)'.format(latex).replace('*', '\cdot').replace('\e', 'e')

                # Make the table into LaTeX
                table = filter_table(ld.spam_results, profile=profile)
                co_cols = [c for c in ld.spam_results.colnames if c.startswith('c') and c not in ['coeffs', 'color']]
                table = table[['wave_eff', 'wave_min', 'wave_max'] + co_cols]
                table.rename_column('wave_eff', '\(\lambda_\mbox{eff}\hspace{5px}(\mu m)\)')
                table.rename_column('wave_min', '\(\lambda_\mbox{min}\hspace{5px}(\mu m)\)')
                table.rename_column('wave_max', '\(\lambda_\mbox{max}\hspace{5px}(\mu m)\)')

                # Add the results to the lists
                html_table = '\n'.join(table.pformat(max_width=-1, max_lines=-1, html=True)).replace('<table', '<table id="myTable" class="table table-striped table-hover"')

                # Add the table title
                header = '<br></br><strong>{}</strong><br><p>\(I(\mu)/I(\mu=1)\) = {}</p>'.format(profile, poly)
                html_table = header + html_table
                profile_spam_tables.append(html_table)

        return render_template('limb_darkening_results.html', form=form,
                               table=profile_tables, spam_table=profile_spam_tables,
                               script=script, plot=div, spam_file_as_string=repr(spam_file_as_string),
                               file_as_string=repr(file_as_string),
                               filt_plot=filt_plot, filt_script=filt_script,
                               js=js_resources, css=css_resources)

    return render_template('limb_darkening.html', form=form)


@app_exoctk.route('/limb_darkening_error', methods=['GET', 'POST'])
def limb_darkening_error():
    """The limb darkening error page

    Returns
    -------
    ``flask.render_template`` obj
        The rendered template for the limb-darkening error page.
    """

    return render_template('limb_darkening_error.html')


@app_exoctk.route('/phase_constraint', methods=['GET', 'POST'])
def phase_constraint(transit_type='primary'):
    """Render page for the phase-constraint calculator

    Returns
    -------
    ``flask.render_template`` obj
        The rendered template for the phase-constraint calculator.
    """

    # Load default form
    form = fv.PhaseConstraint()

    # Reload page with stellar data from ExoMAST
    if form.resolve_submit.data:
        if form.targname.data.strip() != '':
            try:
                # Resolve the target in exoMAST
                form.targname.data = get_canonical_name(form.targname.data)
                data, target_url = get_target_data(form.targname.data)

                # Update the form data
                form.orbital_period.data = data.get('orbital_period')

                t_time = Time(data.get('transit_time'), format='mjd')
                form.transit_time.data = t_time.jd

                form.observation_duration.data = calculate_pre_duration(data.get('transit_duration') * 24.0)

                form.inclination.data = data.get('inclination')
                if form.inclination.data is None:
                    form.inclination.data = np.nan

                form.omega.data = data.get('omega')
                if form.omega.data is None:
                    form.omega.data = np.nan

                form.eccentricity.data = data.get('eccentricity')
                if form.eccentricity.data is None:
                    form.eccentricity.data = np.nan

                form.target_url.data = str(target_url)

                return render_template('phase_constraint.html', form=form)

            except Exception:
                form.target_url.data = ''
                form.targname.errors = ["Sorry, could not resolve '{}' in exoMAST.".format(form.targname.data)]

    # Extract transit type:
    transit_type = form.transit_type.data
    if form.validate_on_submit() and form.calculate_submit.data:

        # Log the form inputs
        log_exoctk.log_form_input(request.form, 'phase_constraint', DB)

        if transit_type == 'primary':
            minphase, maxphase = phase_overlap_constraint(target_name=form.targname.data,
                                                          period=form.orbital_period.data,
                                                          pretransit_duration=form.observation_duration.data,
                                                          window_size=form.window_size.data)
        elif transit_type == 'secondary':
            if (0. <= form.eccentricity.data < 1) and (-360. <= form.omega.data <= 360.) and (0 <= form.inclination.data <= 90.):
                # Use dummy time-of-transit as it doesn't matter for the phase-constraint calculation
                # (phase = 1 is always transit)
                minphase, maxphase = phase_overlap_constraint(target_name=form.targname.data,
                                                              period=form.orbital_period.data, t0=1.,
                                                              pretransit_duration=form.observation_duration.data,
                                                              window_size=form.window_size.data, secondary=True,
                                                              ecc=form.eccentricity.data, omega=form.omega.data,
                                                              inc=form.inclination.data)
            else:
                minphase, maxphase = np.nan, np.nan
        else:
            minphase, maxphase = np.nan, np.nan
        form.minimum_phase.data = minphase
        form.maximum_phase.data = maxphase

    # Send it back to the main page
    return render_template('phase_constraint.html', form=form)


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


@app_exoctk.route('/contam_verify', methods=['GET', 'POST'])
def save_contam_pdf():
    """Save the results of the Contamination Science FOV

    Returns
    -------
    ``flask.render_template`` obj
        The rendered template (and attachment) for the Contamination
        FOV.
    """

    RA, DEC = '19:50:50.2400', '+48:04:51.00'
    contam_pdf = contamVerify(RA, DEC, 'NIRISS', [1, 2], binComp=[], PDF='', web=True)
    filename = contam_pdf.split('/')[-1]

    return render_template(contam_pdf, filename, as_attachment=True)


@app_exoctk.route('/fortney_result', methods=['POST'])
def save_fortney_result():
    """Save the results of the Fortney grid"""

    table_string = flask.request.form['data_file']
    return flask.Response(table_string, mimetype="text/dat", headers={"Content-disposition": "attachment; filename=fortney.dat"})


@app_exoctk.route('/generic_result', methods=['POST'])
def save_generic_result():
    """Save the results of the generic grid"""

    table_string = flask.request.form['data_file']
    return flask.Response(table_string, mimetype="text/dat", headers={"Content-disposition": "attachment; filename=generic.dat"})


@app_exoctk.route('/groups_integrations_download')
def groups_integrations_download():
    """Download the groups and integrations calculator data"""

    return send_file(resource_filename('exoctk', 'data/groups_integrations/groups_integrations_input_data.json'), mimetype="text/json", attachment_filename='groups_integrations_input_data.json', as_attachment=True)


@app_exoctk.route('/visib_result', methods=['POST'])
def save_visib_result():
    """Save the results of the Visibility Only calculation

    Returns
    -------
    ``flask.Response`` obj
        flask.Response object with the results of the visibility only
        calculation.
    """
    visib_table = request.form['vis_table']
    targname = request.form['targetname']
    targname = targname.replace(' ', '_')  # no spaces
    instname = request.form['instrumentname']

    resp = make_response(visib_table)
    resp.headers["Content-Disposition"] = "attachment; filename={}_{}_visibility.csv".format(targname, instname)
    resp.headers["Content-Type"] = "text/csv"

    return resp


@app_exoctk.route('/admin')
@requires_auth
def secret_page():
    """Shhhhh! This is a secret page of admin stuff

    Returns
    -------
    ``flask.render_template`` obj
        The rendered template for the admin page.
    """
    # Reload the DB when this page loads so the data is current
    DB = log_exoctk.load_db(DBPATH)

    tables = [i[0] for i in DB.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()]

    log_tables = []
    for table in tables:

        try:
            data = log_exoctk.view_log(DB, table)

            # Add the results to the lists
            html_table = '\n'.join(data.pformat(max_width=500, html=True)).replace('<table', '<table id="myTable" class="table table-striped table-hover"')

        except Exception:
            html_table = '<p>No data to display</p>'

        # Add the table title
        header = '<h3>{}</h3>'.format(table)
        html_table = header + html_table

        log_tables.append(html_table)

    return render_template('admin_page.html', tables=log_tables)


if __name__ == '__main__':

    port = int(os.environ.get('PORT', 5000))
    app_exoctk.run(host='0.0.0.0', port=port)

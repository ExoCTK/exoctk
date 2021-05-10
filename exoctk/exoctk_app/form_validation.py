import numpy as np
import os

from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, DecimalField, RadioField, SelectField, SelectMultipleField, IntegerField, FloatField
from wtforms.validators import InputRequired, Length, NumberRange, AnyOf, ValidationError
from wtforms.widgets import ListWidget, CheckboxInput

from exoctk.modelgrid import ModelGrid
from exoctk.utils import get_env_variables, FILTERS_LIST, PROFILES
from svo_filters import svo


class MultiCheckboxField(SelectMultipleField):
    """Makes a list of checkbox inputs"""
    widget = ListWidget(prefix_label=False)
    option_widget = CheckboxInput()


class BaseForm(FlaskForm):
    """A generic form with target resolve built in"""
    # Target Resolve
    targname = StringField('targname', default='')
    target_url = StringField('target_url', default='')

    # Submit button
    resolve_submit = SubmitField('Resolve Target')


class FortneyModelForm(BaseForm):
    """Form validation for the forward model tools"""
    # Parameters
    planet_teff = SelectField('planet_teff', choices=[(500, '500'), (750, '750'), (1000, '1000'), (1250, '1250'), (1500, '1500'), (1750, '1750'), (2000, '2000'), (2250, '2250'), (2500, '2500')], validators=[InputRequired('An effective temperature is required')])
    planet_mass = DecimalField('planet_mass', default=1.5, validators=[InputRequired('A planet mass is required!'), NumberRange(min=0.0, message='Planet mass must be positive')])
    planet_mass_unit = SelectField('planet_mass_unit', choices=[('M_jup', 'Jupiter Mass'), ('kilogram', 'kilogram'), ('g', 'gram'), ('M_earth', 'Earth Mass'), ('M_sun', 'Solar Mass')], validators=[InputRequired('A mass unit is required')])
    planet_radius = DecimalField('planet_radius', default=1.25, validators=[InputRequired('A planet radius is required!'), NumberRange(min=0, message='Planet radius must be positive')])
    planet_radius_unit = SelectField('planet_radius_unit', choices=[('R_jup', 'Jupiter Radius'), ('kilometer', 'kilometer'), ('m', 'meter'), ('R_earth', 'Earth Radius'), ('R_sun', 'Solar Radius')], validators=[InputRequired('A planet radius unit is required')])
    stellar_radius = DecimalField('stellar_radius', default=1.0, validators=[InputRequired('A stellar radius is required!'), NumberRange(min=0, message='Stellar radius must be positive')])
    stellar_radius_unit = SelectField('stellar_radius_unit', choices=[('R_sun', 'Solar Radius'), ('R_jup', 'Jupiter Radius'), ('kilometer', 'kilometer'), ('m', 'meter'), ('R_earth', 'Earth Radius')], validators=[InputRequired('A stellar radius unit is required')])
    chemistry = SelectField('chemistry', choices=[('noTiO', '500'), ('eqchem', '750'), (1000, '1000'), (1250, '1250'), (1500, '1500'), (1750, '1750'), (2000, '2000'), (2250, '2250'), (2500, '2500')], validators=[InputRequired('A chemistry type is required')])
    clouds = SelectField('clouds', choices=[('0', 'Nothing'), ('ray10', 'Weak Rayleigh'), ('ray100', 'Medium Rayleigh'), ('ray1000', 'Strong Rayleigh'), ('flat10', 'Weak Cloud'), ('flat100', 'Medium Cloud'), ('flat1000', 'Strong Cloud')], validators=[InputRequired('A cloud model is required')])

    # Form submits
    calculate_submit = SubmitField('Calculate Forward Model')


class LimbDarkeningForm(BaseForm):
    """Form validation for the limb_darkening tool"""
    # Model grid
    modelgrid_dir = get_env_variables()['modelgrid_dir']
    default_modelgrid = os.path.join(modelgrid_dir, 'ATLAS9/')
    mg = ModelGrid(default_modelgrid, resolution=500)
    teff_rng = mg.Teff_vals.min(), mg.Teff_vals.max()
    logg_rng = mg.logg_vals.min(), mg.logg_vals.max()
    feh_rng = mg.FeH_vals.min(), mg.FeH_vals.max()
    modeldir = RadioField('modeldir', default=default_modelgrid, choices=[(os.path.join(modelgrid_dir, 'ATLAS9/'), 'Kurucz ATLAS9'), (os.path.join(modelgrid_dir, 'ACES/'), 'Phoenix ACES')], validators=[InputRequired('A model grid is required!')])

    # Stellar parameters
    teff = DecimalField('teff', default=3500, validators=[InputRequired('An effective temperature is required!'), NumberRange(min=float(teff_rng[0]), max=float(teff_rng[1]), message='Effective temperature must be between {} and {} for this model grid'.format(*teff_rng))])
    logg = DecimalField('logg', default=4.5, validators=[InputRequired('A surface gravity is required!'), NumberRange(min=float(logg_rng[0]), max=float(logg_rng[1]), message='Surface gravity must be between {} and {} for this model grid'.format(*logg_rng))])
    feh = DecimalField('feh', default=0.0, validators=[InputRequired('A surface gravity is required!'), NumberRange(min=float(feh_rng[0]), max=float(feh_rng[1]), message='Metallicity must be between {} and {} for this model grid'.format(*feh_rng))])
    mu_min = DecimalField('mu_min', default=0.1, validators=[InputRequired('A minimum mu value is required!'), NumberRange(min=0.0, max=1.0, message='Minimum mu must be between 0 and 1')])

    # LD profile
    profiles = MultiCheckboxField('profiles', choices=[(x, x) for x in PROFILES], validators=[InputRequired('At least one profile is required!')])

    # Bandpass
    default_filter = 'Kepler.K'
    defilt = svo.Filter(default_filter)
    bandpass = SelectField('bandpass', default=default_filter, choices=[('tophat', 'Top Hat')] + [(filt, filt) for filt in FILTERS_LIST], validators=[InputRequired('A filter is required!')])
    wave_min = DecimalField('wave_min', default=defilt.wave_min.value, validators=[NumberRange(min=0, max=30, message='Minimum wavelength must be between 0 and 30 microns!')])
    wave_max = DecimalField('wave_max', default=defilt.wave_max.value, validators=[NumberRange(min=0, max=30, message='Maximum wavelength must be between 0 and 30 microns!')])
    n_bins = IntegerField('n_bins', default=1)

    # Form submits
    calculate_submit = SubmitField('Calculate Coefficients')
    filter_submit = SubmitField('Filter Selected')
    modelgrid_submit = SubmitField('Model Grid Selected')


class GroupsIntsForm(BaseForm):
    """Form validation for the groups_integrations tool"""
    # Form submits
    calculate_submit = SubmitField('Calculate Groups and Integrations')

    # Stellar Parameters
    kmag = DecimalField('kmag', default=10.5, validators=[InputRequired('A K-band magnitude is required!'), NumberRange(min=5.1, max=11.9, message='K-band mag must be between 5-12, non-inclusive.')])
    obs_duration = DecimalField('obs_duration', default=3, validators=[InputRequired('An observation duration is required!'), NumberRange(min=0, message='Observation duration must be a positive number')])
    time_unit = SelectField('time_unit', default='hour', choices=[('hour', 'hours'), ('day', 'days')])
    models = [('a0i', 'A0I 9750 2.0'), ('aov', 'A0V 9500 2.0'), ('a1v', 'A1V 9250 4.0'), ('a5i', 'A5I 8500 2.0'), ('a3v', 'A3V 8250 4.0'), ('a5v', 'A5V 8250 4.0'), ('f0i', 'F0I 7750 2.0'), ('f0v', 'F0V 7250 1.5'), ('f5i', 'F5I 7000 4.0'), ('f2v', 'F2V 7000 4.0'), ('f5v', 'F5V 6500 4.0'), ('f8v', 'F8V 6250 4.5'), ('g0v', 'G0V 6000 4.5'), ('g0iii', 'G0III 5750 3.0'), ('g2v', 'G2V 5750 4.5'), ('g5v', 'G5V 5750 4.5'), ('g0i', 'G0I 5500 1.5'), ('g8v', 'G8V 5500 4.5'), ('g5iii', 'G5III 5250 2.5'), ('g5i', 'G5I 4740 1.0'), ('k0v', 'K0V 5250 4.5'), ('k0iii', 'K0III 4750 2.0'), ('k2v', 'K2V 4750 4.5'), ('k0i', 'K0I 4500 1.0'), ('k5v', 'K5V 4250 1.5'), ('k5iii', 'K5III 4000 1.5'), ('k7v', 'K7V 4000 4.5'), ('k5i', 'K5I 3750 0.5'), ('m0i', 'M0I 3750 0.0'), ('m0iii', 'M0III 3750 1.5'), ('m0v', 'M0V 3750 4.5'), ('m2i', 'M2I 3500 0.0'), ('m2v', 'M2V 3500 4.5'), ('m5v', 'M5V 3500 5.0')]
    mod = SelectField('mod', choices=models)
    n_group = IntegerField('n_group', default=0)
    ins = SelectField('ins', default='miri', choices=[('niriss', 'NIRISS'), ('nircam', 'NIRCam'), ('nirspec', 'NIRSpec'), ('miri', 'MIRI')])

    # Filter selects
    miri_filt = SelectField('miri_filt', choices=[('lrs', 'LRS')])
    nirspec_filt = SelectField('nirspec_filt', choices=[('f070lp_g140h', 'F070LP/G140H'), ('f100lp_g140h', 'F100LP/G140H'), ('f070lp_g140m', 'F070LP/G140M'), ('f100lp_g140m', 'F100LP/G140M'), ('f170lp_g235h', 'F170LP/G235H'), ('f170lp_g235m', 'F170LP/G235M'), ('f290lp_g395h', 'F290LP/G395H'), ('f290lp_g395m', 'F290LP/G395M')])
    niriss_filt = SelectField('niriss_filt', choices=[('soss', 'SOSS')])
    nircam_filt = SelectField('nircam_filt', choices=[('f322w2', 'F322W2'), ('f444w', 'F444W'), ('f277w', 'F277W')])

    # TA filter selects
    miri_filt_ta = SelectField('miri_filt_ta', choices=[('f560w', 'F560W'), ('f100w', 'F100W'), ('f1500w', 'F1500W')])
    nirspec_filt_ta = SelectField('nirspec_filt_ta', choices=[('f110w', 'F110W'), ('f140x', 'F140X'), ('clear', 'CLEAR')])
    niriss_filt_ta = SelectField('niriss_filt_ta', choices=[('f480m', 'F480M')])
    nircam_filt_ta = SelectField('nircam_filt_ta', choices=[('f335m', 'F335M')])

    # Subarray selects
    miri_subarray = SelectField('miri_subarray', choices=[('slitlessprism', 'SLITLESSPRISM')])
    nirspec_subarray = SelectField('nirspec_subarray', choices=[('sub2048', 'SUB2048'), ('sub1024a', 'SUB1024A'), ('sub1024b', 'SUB1024B'), ('sub512', 'SUB512')])
    niriss_subarray = SelectField('niriss_subarray', choices=[('substrip256', 'SUBSTRIP256'), ('substrip96', 'SUBSTRIP96')])
    nircam_subarray = SelectField('nircam_subarray', choices=[('full', 'FULL FRAME'), ('subgrism256', 'SUBGRISM256'), ('subgrism128', 'SUBGRISM128'), ('subgrism64', 'SUBGRISM64')])

    # TA subarray selects
    miri_subarray_ta = SelectField('miri_subarray_ta', choices=[('slitlessprism', 'SLITLESSPRISM')])
    nirspec_subarray_ta = SelectField('nirspec_subarray_ta', choices=[('full', 'FULL'), ('sub32', 'SUB32'), ('sub2048', 'SUB2048')])
    niriss_subarray_ta = SelectField('niriss_subarray_ta', choices=[('nrm', 'SUBTASOSS -- BRIGHT'), ('im', 'SUBTASOSS -- FAINT')])
    nircam_subarray_ta = SelectField('nircam_subarray_ta', choices=[('sub32tats', 'SUB32TATS')])

    # Saturation
    sat_mode = RadioField('sat_mode', default='well', choices=[('counts', 'Counts'), ('well', 'Full well fraction')])
    sat_max = DecimalField('sat_max', default=0.95, validators=[InputRequired('A saturation level is required!'), NumberRange(min=0.0, message='Saturation level must be positive.')])


class ContamVisForm(BaseForm):
    """Form validation for the contamination_visibility tool"""
    # Form submits
    calculate_submit = SubmitField('Calculate Visibility')
    calculate_contam_submit = SubmitField('Calculate Visibility and Contamination')
    mode_submit = SubmitField('Mode Selected')

    # Form inputs
    ra = DecimalField('ra', validators=[NumberRange(min=0, max=360, message='RA must be between 0 and 360 degrees')])
    dec = DecimalField('dec', validators=[NumberRange(min=-90, max=90, message='Declinaton must be between -90 and 90 degrees')])
    inst = SelectField('inst', choices=[('NIS_SUBSTRIP256', 'NIRISS - SOSS - SUBSTRIP256'), ('NIS_SUBSTRIP96', 'NIRISS - SOSS - SUBSTRIP96'), ('NRCA5_GRISM256_F322W2', 'NIRCam - Grism Time Series - F322W2'), ('NRCA5_GRISM256_F444W', 'NIRCam - Grism Time Series - F444W'), ('MIRI_SLITLESSPRISM', 'MIRI - LRS'), ('NIRSpec', 'NIRSpec (Visibility Only)')])
    companion = StringField('companion', default='')
    pa_min = DecimalField('pa_min', default=0, validators=[NumberRange(min=0, max=360, message='Minimum PA must be between 0 and 360 degrees')])
    pa_max = DecimalField('pa_max', default=360, validators=[NumberRange(min=0, max=360, message='Maximum PA must be between 0 and 360 degrees')])


class PhaseConstraint(BaseForm):
    """Form validation for the phase-constraint tool"""

    calculate_submit = SubmitField('Calculate Phase Constraint')

    orbital_period = FloatField('orbital_period', validators=[InputRequired('Orbital period is a required field')]) 
    eccentricity = FloatField('eccentricity', default=np.nan)
    transit_type = SelectField('transit_type', choices=[('primary', 'primary'), ('secondary', 'secondary')])
    omega = FloatField('omega', default=np.nan)
    inclination = FloatField('inclination', default=np.nan)
    transit_time = FloatField('transit_time', default=np.nan)
    window_size = FloatField('window_size', default=1.0)
    observation_duration = FloatField('observation_duration', default=2.0, validators=[InputRequired('Observation duration is a required field.')])
    minimum_phase = DecimalField('minimum_phase', default=0.0)
    maximum_phase = DecimalField('maximum_phase', default=0.0)
    minimum_phase_sec = DecimalField('minimum_phase_sec', default=0.0)
    maximum_phase_sec = DecimalField('maximum_phase_sec', default=0.0)

from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, DecimalField, RadioField, SelectField, SelectMultipleField, IntegerField
from wtforms.validators import InputRequired, Length, NumberRange, AnyOf
from wtforms.widgets import ListWidget, CheckboxInput

from exoctk.modelgrid import ModelGrid
from exoctk.utils import MODELGRID_DIR, FILTERS, PROFILES
from svo_filters import svo


class MultiCheckboxField(SelectMultipleField):
    """Makes a list of checkbox inputs"""
    widget = ListWidget(prefix_label=False)
    option_widget = CheckboxInput()


class LdcForm(FlaskForm):
    """Form validation for the LDC tool"""
    # Model grid
    default_modelgrid = MODELGRID_DIR
    mg = ModelGrid(default_modelgrid, resolution=500)
    teff_rng = mg.Teff_vals.min(), mg.Teff_vals.max()
    logg_rng = mg.logg_vals.min(), mg.logg_vals.max()
    feh_rng = mg.FeH_vals.min(), mg.FeH_vals.max()
    modeldir = RadioField('modeldir', default=default_modelgrid, choices=[('/user/jfilippazzo/Models/ACES/default', 'Phoenix ACES'), ('/user/jfilippazzo/Models/ATLAS9/default', 'Kurucz ATLAS9')], validators=[InputRequired('A model grid is required!')])

    # Target Resolve
    targname = StringField('targname', default='')

    # Stellar parameters
    teff = DecimalField('teff', default=3500, validators=[InputRequired('An effective temperature is required!'), NumberRange(min=teff_rng[0], max=teff_rng[1], message='Effective temperature must be between {} and {}'.format(*teff_rng))])
    logg = DecimalField('logg', default=4.5, validators=[InputRequired('A surface gravity is required!'), NumberRange(min=logg_rng[0], max=logg_rng[1], message='Surface gravity must be between {} and {}'.format(*logg_rng))])
    feh = DecimalField('feh', default=0.0, validators=[InputRequired('A surface gravity is required!'), NumberRange(min=feh_rng[0], max=feh_rng[1], message='Metallicity must be between {} and {}'.format(*feh_rng))])
    mu_min = DecimalField('mu_min', default=0.1, validators=[InputRequired('A minimum mu value is required!'), NumberRange(min=0.0, max=1.0, message='Minimum mu must be between 0 and 1')])

    # LD profile
    profiles = MultiCheckboxField('profiles', choices=[(x, x) for x in PROFILES], validators=[InputRequired('At least one profile is required!')])

    # Bandpass
    default_filter = 'Kepler.K'
    defilt = svo.Filter(default_filter)
    bandpass = SelectField('bandpass', default=default_filter, choices=[('tophat', 'Top Hat')]+[(filt, filt) for filt in sorted(FILTERS['Band'])], validators=[InputRequired('A filter is required!')])
    wave_min = DecimalField('wave_min', default=defilt.wave_min.value, validators=[NumberRange(min=0, max=30, message='Minimum wavelength must be between 0 and 30 microns!')])
    wave_max = DecimalField('wave_max', default=defilt.wave_max.value, validators=[NumberRange(min=0, max=30, message='Maximum wavelength must be between 0 and 30 microns!')])
    n_bins = IntegerField('n_bins', default=1)

    # Form submits
    resolve_submit = SubmitField('Resolve Target')
    calculate_submit = SubmitField('Calculate Coefficients')
    filter_submit = SubmitField('Filter Selected')
    modelgrid_submit = SubmitField('Model Grid Selected')

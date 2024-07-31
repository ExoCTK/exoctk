#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module to calculate limb darkening coefficients from a grid of model spectra
"""

from copy import copy
import inspect
import os
import warnings

from astropy.io import ascii as ii
import astropy.table as at
from astropy.utils.exceptions import AstropyWarning
import bokeh.plotting as bkp
from bokeh.models import Range1d, TabPanel, Tabs
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from svo_filters import svo

from .. import modelgrid
from .. import utils
from . import spam

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica'], 'size': 16})
rc('text', usetex=True)

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=FutureWarning)


def ld_profile(name='quadratic', latex=False):
    """
    Define the function to fit the limb darkening profile

    Reference:
        https://www.cfa.harvard.edu/~lkreidberg/batman/
        tutorial.html#limb-darkening-options

    Parameters
    ----------
    name: str
        The name of the limb darkening profile function to use,
        including 'linear', 'quadratic', 'square-root',
        'logarithmic', 'exponential', '3-parameter', and '4-parameter'
    latex: bool
        Return the function as a LaTeX formatted string

    Returns
    -------
    function, str
        The corresponding function for the given profile

    """
    # Supported profiles a la BATMAN
    names = ['linear', 'quadratic', 'square-root', 'logarithmic', 'exponential', '3-parameter', '4-parameter']

    # Check that the profile is supported
    if name in names:

        # Linear
        if name == 'linear':
            def profile(m, c1):
                return 1. - c1 * (1. - m)

        # Quadratic
        if name == 'quadratic':
            def profile(m, c1, c2):
                return 1. - c1 * (1. - m) - c2 * (1. - m)**2

        # Square-root
        if name == 'square-root':
            def profile(m, c1, c2):
                return 1. - c1 * (1. - m) - c2 * (1. - np.sqrt(m))

        # Logarithmic
        if name == 'logarithmic':
            def profile(m, c1, c2):
                return 1. - c1 * (1. - m) - c2 * m * np.log(m)

        # Exponential
        if name == 'exponential':
            def profile(m, c1, c2):
                return 1. - c1 * (1. - m) - c2 / (1. - np.e**m)

        # 3-parameter
        if name == '3-parameter':
            def profile(m, c1, c2, c3):
                return 1. - c1 * (1. - m) - c2 * (1. - m**1.5) - c3 * (1. - m**2)

        # 4-parameter
        if name == '4-parameter':
            def profile(m, c1, c2, c3, c4):
                return 1. - c1 * (1. - m**0.5) - c2 * (1. - m) - c3 * (1. - m**1.5) - c4 * (1. - m**2)

        if latex:
            profile = inspect.getsource(profile).replace('\n', '')
            profile = profile.replace('\\', '').split('return ')[1]

            for i, j in [('**', '^'), ('m', '\mu'), (' ', ''), ('np.', '\\'), ('0.5', '{0.5}'), ('1.5', '{1.5}')]:
                profile = profile.replace(i, j)

        return profile

    else:
        print("'{}' is not a supported profile. Try".format(name), names)
        return


class LDC:
    """A class to hold all the LDCs you want to run

    Example
    -------
    from exoctk.limb_darkening import limb_darkening_fit as lf
    from svo_filters import Filter
    ld = lf.LDC(model_grid='ACES')
    bp = Filter('WFC3_IR.G141', n_bins=5)
    ld.calculate(4000, 4.5, 0.0, 'quadratic', bandpass=bp)
    ld.calculate(4000, 4.5, 0.0, '4-parameter', bandpass=bp)
    ld.spam(planet_name='WASP-12b', profiles=['quadratic', 'logarithmic'])
    ld.plot_tabs(show=True)
    """
    def __init__(self, model_grid='ACES'):
        """Initialize an LDC object

        Parameters
        ----------
        model_grid: exoctk.modelgrid.ModelGrid
            The grid of synthetic spectra from which the coefficients
            will be calculated
        """
        # Try ACES or ATLAS
        if model_grid == 'ACES':
            model_grid = modelgrid.ACES(resolution=700)
        if model_grid == 'ATLAS9':
            model_grid = modelgrid.ATLAS9(resolution=700)

        # Make sure it's a ModelGrid object
        if not isinstance(model_grid, modelgrid.ModelGrid):
            raise TypeError("'model_grid' must be a exoctk.modelgrid.ModelGrid object.")

        # Load the model grid
        self.model_grid = model_grid
        self.model_cache = {}

        # Table for results
        columns = ['name', 'Teff', 'logg', 'FeH', 'profile', 'filter', 'models', 'coeffs', 'errors', 'wave', 'wave_min', 'wave_eff', 'wave_max', 'scaled_mu', 'raw_mu', 'mu_min', 'scaled_ld', 'raw_ld', 'ld_min', 'ldfunc', 'flux', 'bandpass', 'color']
        dtypes = ['|S20', float, float, float, '|S20', '|S20', '|S20', object, object, object, np.float16, np.float16, np.float16, object, object, np.float16, object, object, np.float16, object, object, object, '|S20']
        self.results = at.Table(names=columns, dtype=dtypes)

        # Table for spam results
        self.spam_results = None

        # Colors for plotting
        self.ld_color = {'quadratic': 'blue', '4-parameter': 'red', 'exponential': 'green', 'linear': 'orange', 'square-root': 'cyan', '3-parameter': 'magenta', 'logarithmic': 'pink'}

        self.count = 1

    @staticmethod
    def bootstrap_errors(mu_vals, func, coeffs, errors, n_samples=1000):
        """
        Bootstrapping LDC errors

        Parameters
        ----------
        mu_vals: sequence
            The mu values
        func: callable
            The LD profile function
        coeffs: sequence
            The coefficients
        errors: sequence
            The errors on each coeff
        n_samples: int
            The number of samples

        Returns
        -------
        tuple
            The lower and upper errors
        """
        # Generate n_samples
        vals = []
        for n in range(n_samples):
            co = np.random.normal(coeffs, errors)
            vals.append(func(mu_vals, *co))

        # r = np.array(list(zip(*vals)))
        dn_err = np.min(np.asarray(vals), axis=0)
        up_err = np.max(np.asarray(vals), axis=0)

        return dn_err, up_err

    def calculate(self, Teff, logg, FeH, profile, mu_min=0.05, ld_min=0.01,
                  bandpass=None, name=None, color=None, interp=False, **kwargs):
        """
        Calculates the limb darkening coefficients for a given synthetic
        spectrum. If the model grid does not contain a spectrum of the
        given parameters, the grid is interpolated to those parameters.

        Reference for limb-darkening laws:
        http://www.astro.ex.ac.uk/people/sing/David_Sing/Limb_Darkening.html

        Parameters
        ----------
        Teff: int
            The effective temperature of the model
        logg: float
            The logarithm of the surface gravity
        FeH: float
            The logarithm of the metallicity
        profile: str
            The name of the limb darkening profile function to use,
            including 'linear', 'quadratic', 'square-root',
            'logarithmic', 'exponential', and '4-parameter'
        mu_min: float
            The minimum mu value to consider
        ld_min: float
            The minimum limb darkening value to consider
        bandpass: svo_filters.svo.Filter() (optional)
            The photometric filter through which the limb darkening
            is to be calculated
        name: str (optional)
            A name for the calculation
        color: str (optional)
            A color for the plotted result
        interp: bool
            Interpolate spectra to given model parameters
        """
        # Define the limb darkening profile function
        ldfunc = ld_profile(profile)

        if not ldfunc:
            raise ValueError("No such LD profile:", profile)

        # Get the grid point (and update parameters if no interpolation)
        grid_point = self.get_model(Teff, logg, FeH, interp=interp)
        Teff = grid_point['Teff']
        logg = grid_point['logg']
        FeH = grid_point['FeH']

        # Retrieve the wavelength, flux, mu, and effective radius
        wave = grid_point.get('wave')
        flux = grid_point.get('flux')
        mu = grid_point.get('mu').squeeze()

        # Try to get bandpass if it is a string
        if isinstance(bandpass, str):
            try:
                bandpass = svo.Filter(bandpass, **kwargs)
            except Exception:
                raise ValueError("Could not find a bandpass named '{}'. Try passing a 'svo_filters.svo.Filter' object instead.".format(bandpass))

        # Use tophat if no bandpass
        if bandpass is None:
            units = self.model_grid.wave_units
            bandpass = svo.Filter('tophat', wave_min=np.min(wave) * units, wave_max=np.max(wave) * units)

        # Check if a bandpass is provided
        if not isinstance(bandpass, svo.Filter):
            raise TypeError("Invalid bandpass of type", type(bandpass))

        # # Make sure the bandpass has coverage
        # bp_min = bandpass.wave_min.value
        # bp_max = bandpass.wave_max.value
        # mg_min = self.model_grid.wave_rng[0].value
        # mg_max = self.model_grid.wave_rng[-1].value
        # if bp_min < mg_min or bp_max > mg_max:
        #     raise ValueError('Bandpass {} not covered by model grid of\
        #                       wavelength range {}'.format(bandpass.filterID,
        #                                                   self.model_grid
        #                                                       .wave_rng))

        # Apply the filter
        try:
            flux, _ = bandpass.apply([wave, flux])  # Sometimes this returns a tuple
        except ValueError:
            flux = bandpass.apply([wave, flux])  # Sometimes it returns one value

        # Make rsr curve 3 dimensions if there is only one
        # wavelength bin, then get wavelength only
        bp = bandpass.rsr
        if bp.ndim == 2:
            bp = bp[None, :]
        wave = bp[:, 0, :]

        # Calculate mean intensity vs. mu
        wave = wave[None, :] if wave.ndim == 1 else wave
        flux = flux[None, :] if flux.ndim == 2 else flux
        mean_i = np.nanmean(flux, axis=-1)
        # mean_i[mean_i == 0] = np.nan

        # Calculate limb darkening, I[mu]/I[1] vs. mu
        ld = mean_i / mean_i[:, np.where(mu == np.nanmax(mu))].squeeze(axis=-1)

        # Rescale mu values to make f(mu=0)=ld_min
        # for the case where spherical models extend beyond limb
        ld_avg = np.nanmean(ld, axis=0)
        muz = np.interp(ld_min, ld_avg, mu) if any(ld_avg < ld_min) else 0
        mu = (mu - muz) / (1 - muz)

        # Trim to useful mu range
        imu, = np.where(mu > mu_min)
        scaled_mu, scaled_ld = mu[imu], ld[:, imu]

        # Get effective wavelengths
        wave_effs = np.mean(bandpass.wave, axis=1)

        # Fit limb darkening coefficients for each wavelength bin
        for n, ldarr in enumerate(scaled_ld):

            # Get effective wavelength of bin
            wave_eff = wave_effs[n]

            try:

                # Fit polynomial to data
                coeffs, cov = curve_fit(ldfunc, scaled_mu, ldarr, method='lm')

                # Calculate errors from covariance matrix diagonal
                errs = np.sqrt(np.diag(cov))

                # Make a dictionary or the results
                result = {}

                # Check the count
                result['name'] = name or 'Calculation {}'.format(self.count)
                self.count += 1

                if bandpass.wave.shape[0] == len(scaled_ld) and name is None:
                    result['name'] = '{:.3f}'.format(wave_eff)

                # Set a color if possible
                result['color'] = color or self.ld_color[profile]

                # Add the results
                result['Teff'] = Teff
                result['logg'] = logg
                result['FeH'] = FeH
                result['filter'] = bandpass.filterID
                result['models'] = self.model_grid.path
                result['raw_mu'] = mu
                result['raw_ld'] = ld[n]
                result['scaled_mu'] = scaled_mu
                result['scaled_ld'] = ldarr
                result['flux'] = flux[n]
                result['wave'] = wave[n]
                result['mu_min'] = mu_min
                result['bandpass'] = bandpass
                result['ldfunc'] = ldfunc
                result['coeffs'] = coeffs
                result['errors'] = errs
                result['profile'] = profile
                result['n_bins'] = bandpass.n_bins
                result['pixels_per_bin'] = bandpass.pixels_per_bin
                result['wave_min'] = wave[n, 0].round(5)
                result['wave_eff'] = wave_eff
                result['wave_max'] = wave[n, -1].round(5)

                # Add the coeffs
                for n, (coeff, err) in enumerate(zip(coeffs, errs)):
                    cname = 'c{}'.format(n + 1)
                    ename = 'e{}'.format(n + 1)
                    result[cname] = coeff.round(3)
                    result[ename] = err.round(3)

                    # Add the coefficient column to the table if not present
                    if cname not in self.results.colnames:
                        self.results[cname] = [np.nan] * len(self.results)
                        self.results[ename] = [np.nan] * len(self.results)

                # Add the new row to the table
                result = {i: j for i, j in result.items() if i in self.results.colnames}
                self.results.add_row(result)

            except ValueError:
                print("Could not calculate coefficients at {}".format(wave_eff))

    def get_model(self, teff, logg, feh, interp=False):
        """
        Method to grab cached model or fetch new one

        Parameters
        ----------
        teff: float
            The effective temperature of the desired model
        logg: float
            The log surface gravity of the desired model
        feh: float
            The metallicity of the desired model
        interp: bool
            Interpolate the model spectra to the given parameters

        Returns
        -------
        dict
            The stellar intensity model
        """
        if not interp:

            teff_val, logg_val, feh_val = teff, logg, feh

            # Find the closest of each parameter
            teff = min(self.model_grid.Teff_vals, key=lambda x: abs(x - teff))
            logg = min(self.model_grid.logg_vals, key=lambda x: abs(x - logg))
            feh = min(self.model_grid.FeH_vals, key=lambda x: abs(x - feh))
            print('Closest model to [{}, {}, {}] => [{}, {}, {}]'.format(teff_val, logg_val, feh_val, teff, logg, feh))

        # Check if it is stored
        params = '{}_{}_{}_{}'.format(self.model_grid.path.split('/')[-2], teff, logg, feh)

        # Get cached or get new model
        if params in self.model_cache:
            model = self.model_cache[params]
        else:
            model = self.model_grid.get(teff, logg, feh, interp=interp)
            self.model_cache[params] = model
            print("Saving model '{}'".format(params))

        return model

    def plot(self, fig=None, show=False, **kwargs):
        """Plot the LDCs

        Parameters
        ----------
        fig: matplotlib.pyplot.figure, bokeh.plotting.figure (optional)
            An existing figure to plot on
        show: bool
            Show the figure
        """
        # Separate plotting kwargs from parameter kwargs
        pwargs = {i: j for i, j in kwargs.items() if i in self.results.columns}
        kwargs = {i: j for i, j in kwargs.items() if i not in pwargs.keys()}

        # Filter the table by given kwargs
        table = utils.filter_table(self.results, **pwargs)

        for row in table:

            # Set color and label for plot
            color = row['color']
            label = row['name']

            # Generate smooth curve
            ldfunc = row['ldfunc']
            mu_vals = np.linspace(0, 1, 1000)
            ld_vals = ldfunc(mu_vals, *row['coeffs'])

            # Generate smooth errors
            dn_err, up_err = self.bootstrap_errors(mu_vals, ldfunc, row['coeffs'], row['errors'])

            # Matplotlib fig by default
            if fig is None:
                fig = bkp.figure()

            # Add fits to matplotlib
            if isinstance(fig, matplotlib.figure.Figure):

                # Make axes
                ax = fig.add_subplot(111)

                # Plot the fitted points
                ax.errorbar(row['raw_mu'], row['raw_ld'], c='k', ls='None', marker='o', markeredgecolor='k', markerfacecolor='None')

                # Plot the mu cutoff
                ax.axvline(row['mu_min'], color='0.5', ls=':')

                # Draw the curve and error
                ax.plot(mu_vals, ld_vals, color=color, label=label, **kwargs)
                ax.fill_between(mu_vals, dn_err, up_err, color=color, alpha=0.1)
                ax.set_ylim(0, 1)
                ax.set_xlim(0, 1)

            # Or to bokeh!
            else:

                # Set the plot elements
                fig.x_range = Range1d(0, 1)
                fig.y_range = Range1d(0, 1)
                fig.xaxis.axis_label = 'mu'
                fig.yaxis.axis_label = 'Normalized Intensity'
                fig.legend.location = "bottom_right"

                # Plot the fitted points
                fig.circle(row['raw_mu'], row['raw_ld'], fill_color='black')

                # Plot the mu cutoff
                fig.line([row['mu_min']] * 2, [0, 1], legend_label='cutoff', line_color='#6b6ecf', line_dash='dotted')

                # Draw the curve and error
                fig.line(mu_vals, ld_vals, line_color=color, legend_label=label, **kwargs)
                vals = np.append(mu_vals, mu_vals[::-1])
                evals = np.append(dn_err, up_err[::-1])
                fig.patch(vals, evals, color=color, fill_alpha=0.2, line_alpha=0)

        if show:
            if isinstance(fig, matplotlib.figure.Figure):
                plt.xlabel('$\mu$')
                plt.ylabel('$I(\mu)/I(\mu = 1)$')
                plt.legend(loc=0, frameon=False)
                plt.show()
            else:
                bkp.show(fig)

        else:
            return fig

    def plot_tabs(self, show=False, **kwargs):
        """Plot the LDCs in a tabbed figure

        Parameters
        ----------
        fig: matplotlib.pyplot.figure, bokeh.plotting.figure (optional)
            An existing figure to plot on
        show: bool
            Show the figure
        """
        # Change names to reflect ld profile
        old_names = self.results['name']
        for n, row in enumerate(self.results):
            self.results[n]['name'] = row['profile']

        # Draw a figure for each wavelength bin
        tabs = []
        for wav in np.unique(self.results['wave_eff']):

            # Plot it
            TOOLS = 'box_zoom, box_select, crosshair, reset, hover'
            fig = bkp.figure(tools=TOOLS, x_range=Range1d(0, 1), y_range=Range1d(0, 1), width=800, height=400)
            self.plot(wave_eff=wav, fig=fig)

            # Plot formatting
            fig.legend.location = 'bottom_right'
            fig.xaxis.axis_label = 'mu'
            fig.yaxis.axis_label = 'Intensity'

            tabs.append(TabPanel(child=fig, title=str(wav)))

        # Make the final tabbed figure
        final = Tabs(tabs=tabs)

        # Put the names back
        self.results['name'] =  old_names

        if show:
            bkp.show(final)
        else:
            return final

    def save(self, filepath):
        """
        Save the LDC results to file

        Parameters
        ----------
        filepath: str
            The complete filepath to save the results to
        """
        # Make sure it is a string
        if not isinstance(filepath, str):
            raise TypeError("{}: Expecting a string for 'filepath' argument".format(type(filepath)))

        # Make sure it's a valid path
        if not os.path.exists(os.path.dirname(filepath)):
            raise ValueError("{}: Not a valid path")

        # Copy the results table
        keep_cols = ['Teff', 'logg', 'FeH', 'profile', 'filter', 'models', 'wave_min', 'wave_eff', 'wave_max', 'c1', 'e1', 'c2', 'e2', 'c3', 'e3', 'c4', 'e4']
        results = self.results[keep_cols]

        # Save to file
        ii.write(results, filepath, format='fixed_width_two_line')

    def spam(self, planet_name=None, planet_data=None, profiles=['quadratic'], **kwargs):
        """
        Calculates SPAM coefficients by transforming non-linear coefficients

        Parameters
        ----------
        planet_name: string
            The name of the input planet (e.g., 'WASP-19b'); this will be used to query the planet properties from MAST.
        planet_data: dict
            Dictionary containing the planet properties. Must include 'transit_duration', 'orbital_period' (days), 'Rp/Rs', 'a/Rs', 'inclination' (degrees), 'eccentricity' and 'omega' (degrees)
        profiles: sequence
            The profiles to calculate, ['quadratic', 'logarithmic', 'square-root']
        """
        # Profiles supported by SPAM transformation code
        spam_supported = ['quadratic', 'square-root', 'logarithmic']

        if not all([profile in spam_supported for profile in profiles]):
            raise ValueError("{}: Supported profiles include {}".format(profiles, spam_supported))

        # Require 4-parameter calculation
        if '4-parameter' in self.results['profile']:

            # For each desired profile...
            for profile in profiles:

                # Rows of non-linear calculations
                nonlin = copy(self.results[self.results['profile'] == '4-parameter'])

                # Update profile
                nonlin['profile'] = profile

                # ...and for each wavelength channel...
                for row in nonlin:

                    # Calculate SPAM coefficients
                    (c1, c2), properties = spam.transform_coefficients(row['c1'], row['c2'], row['c3'], row['c4'], ld_law=profile, planet_name=planet_name, planet_data=planet_data, **kwargs)
                    row['c1'], row['c2'] = c1.round(3), c2.round(3)

                # Remove unused columns
                omit_cols = ['c3', 'c4', 'e1', 'e2', 'e3', 'e4']
                nonlin = nonlin[[col for col in nonlin.columns if col not in omit_cols]]

                # Add planet properties to table
                planet_properties = ['transit_duration', 'orbital_period', 'Rp/Rs', 'a/Rs', 'inclination', 'eccentricity', 'omega']
                for prop in planet_properties:
                    nonlin.add_column([round(properties[prop], 4)] * len(nonlin), name=prop)

                # Add the results to the spam table
                if self.spam_results is None:
                    self.spam_results = nonlin
                else:
                    self.spam_results = at.vstack([self.spam_results, nonlin])

        else:

            print("Limb darkening coefficients must first be calculated using the 4-parameter profile to get SPAM coefficient transformations.")

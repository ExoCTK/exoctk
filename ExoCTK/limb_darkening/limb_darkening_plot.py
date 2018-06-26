#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module of plotting tools for the limb darkening subpackage.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import astropy.table as at
from matplotlib import rc

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

COLORS = ['blue', 'red', 'green', 'orange',
          'cyan', 'magenta', 'pink', 'purple']


def bootstrap_errors(mu_vals, func, coeffs, errors, n_samples=1000):
    """
    Bootstrap errors
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


def ld_plot(ldfuncs, grid_point, fig=None,
            colors='blue', bin_idx='', **kwargs):
    """
    Make a LD plot in Bokeh or Matplotlib

    Parameters
    ----------
    ldfuncs: function, list
        The limb darkening function to use
    grid_point: dict
        The model data for the grid point
    fig: matplotlib.figure, bokeh.plotting.figure
        The figure to plot in
    color: str, array-like
        The color(s) to use for the plot
    bin_idx: int (optional)
        The index of the wavelength bin to plot,
        otherwise plot all of them
    """
    # Get actual data points
    if isinstance(bin_idx, int):
        slc = slice(bin_idx, bin_idx+1)
    else:
        slc = slice(None)
    # flux = grid_point['flux'][slc]
    mu = grid_point['mu']
    mu_min = grid_point['mu_min']
    profiles = grid_point['profiles']
    mu_raw = grid_point['scaled_mu']
    ld_raw = grid_point['ld_raw'][slc]
    mu_vals = np.linspace(0, 1, 1000)

    # Make profile list and get colors
    if callable(ldfuncs):
        ldfuncs = [ldfuncs]
    if isinstance(colors, str):
        colors = [colors]
    if len(colors) != len(ldfuncs):
        colors = COLORS[:len(ldfuncs)]

    for color, profile, ldfunc in zip(colors, profiles, ldfuncs):
        # Get the coefficients for the given profile
        table = grid_point[profile]['coeffs']

        if bin_idx != '':
            table = at.Table(table[bin_idx])

        coeffs = table[[k for k in table.colnames if k.startswith('c')]]
        errs = table[[k for k in table.colnames if k.startswith('e')]]

        for n in range(len(table)):
            co = np.asarray(list(coeffs[n]))
            er = np.asarray(list(errs[n]))

            # Evaluate the limb darkening profile fit
            ld_vals = ldfunc(mu_vals, *co)

            # ==========================================
            # ==========================================
            # ==========================================
            # Bootstrap the results to get errors here!
            dn_err, up_err = bootstrap_errors(mu_vals, ldfunc, co, er)
            # dn_err = ldfunc(mu_vals, *co-er)
            # up_err = ldfunc(mu_vals, *co+er)
            # ==========================================
            # ==========================================
            # ==========================================

            if profile == 'uniform':
                ld_vals = [ld_vals]*len(mu_vals)

            if fig is None:
                fig = plt.gcf()

            # Add fits to matplotlib
            if isinstance(fig, matplotlib.figure.Figure):

                # Make axes
                ax = fig.add_subplot(111)

                # Plot the fitted points
                ax.errorbar(mu_raw, ld_raw[n], c='k', ls='None', marker='o',
                            markeredgecolor='k', markerfacecolor='None')

                # Plot the mu cutoff
                ax.axvline(mu_min, color='0.5', ls=':')

                # Draw the curve and error
                ax.plot(mu_vals, ld_vals, color=color, label=profile,
                        **kwargs)
                ax.fill_between(mu_vals, dn_err, up_err, color=color,
                                alpha=0.1)
                ax.set_ylim(0, 1)
                ax.set_xlim(0, 1)

            # Or to bokeh!
            else:

                # Plot the fitted points
                fig.circle(mu, ld_raw[n], fill_color='black')

                # Plot the mu cutoff
                fig.line([mu_min, mu_min], [0, 1], legend='cutoff',
                         line_color='#6b6ecf', line_dash='dotted')

                # Draw the curve and error
                fig.line(mu_vals, ld_vals, line_color=color, legend=profile,
                         **kwargs)
                vals = np.append(mu_vals, mu_vals[::-1])
                evals = np.append(dn_err, up_err[::-1])
                fig.patch(vals, evals, color=color, fill_alpha=0.2,
                          line_alpha=0)

# def ld_v_mu(model_grid, compare, profiles=('quadratic','nonlinear'),
#             **kwargs):
#     """
#     Produce a plot of mu vs. limb darkening values for a range of
#     model grid parameters compared to some interpolated value
#
#     Parameters
#     ----------
#     model_grid: core.ModelGrid
#         The model grid to plot. The ModelGrid.customize() method can be run
#         before hand or arguments can be passed via **kwargs argument
#     compare: list
#         The list of off-grid parameters to compare to the
#         model grid plots
#     profiles: list
#         The limb darkening profiles to include in the comparison
#
#     Example
#     -------
#     >>> model_comparison(model_grid, (2550, 5.22, 0),
#                          **{'Teff_rng':(2500,2600),
#                          'logg_rng':(5,5.5), 'FeH_rng':(0,0)})
#
#     """
#     # Make a copy of the ModelGrid so it doesn't change the
#     # input object in the Python session
#     grid = copy.copy(model_grid)
#     grid.customize(**kwargs)
#
#     # Plot the grids
#     for p,ls in zip(profiles,['-','--',':','-.']):
#         _ = ldcfit.ldc_grid(grid, p, plot=plt.gcf(), **{'ls':ls})
#
#     # Plot the interpolated comparison
#     for p,ls in zip(profiles,['-','--',':','-.']):
#         t, g, m = compare
#         _ = ldcfit.ldc(t, g, m, grid, p, plot=plt.gcf(),
#             **{'ls':ls, 'c':'k', 'lw':3})
#
#     # Delete the copy
#     del grid
#
#     # Plot labels
#     plt.xlabel(r'$\mu$')
#     plt.ylabel(r'$I(\mu)/I(\mu =0)$')

# def ldc_v_wavelength(model_grid, wave_ranges, profile, **kwargs):
#     """
#     Plot limb darkening coefficients vs. wavelength as resolution
#     of intensity spectra. Overplot mean value for specified
#     filters and/or wavelength intervals
#
#     Parameters
#     ----------
#     model_grid: core.ModelGrid
#         The model grid to use for the calculation
#     wave_ranges: list
#         The list of wavelength range tuples in [um],
#         e.g. [(0.6,0.9), (0.9,1.2), (1.2,1.5)]
#     profile: str
#         The limb darkening profile to use
#
#     """
#     # Get the number of coefficients for the limb darkening profile
#     rows = len(inspect.getargspec(ldcfit.ld_profile(profile)).args)-1
#     cols = len(wave_ranges)
#
#     # Initialize limb darkening coefficient, mu, and effecive radius grids
#     T = model_grid.Teff_vals
#     G = model_grid.logg_vals
#     M = model_grid.FeH_vals
#     coeff_grid = np.zeros((cols,rows,len(T),len(G),len(M)))
#     mu_grid = np.zeros((cols,len(T),len(G),len(M)))
#     r_grid = np.zeros((cols,len(T),len(G),len(M)))
#     w_mean, w_unc = [], []
#
#     # Calculate the grid in the given wavelength ranges
#     for n,wr in enumerate(wave_ranges):
#
#         # Get wave center and range
#         w = (wr[0]+wr[1])/2.
#         w_mean.append(w)
#         w_unc.append(wr[1]-w)
#
#         # Make a copy of the ModelGrid so it doesn't change the
#         # input object in the Python session
#         grid = copy.copy(model_grid)
#
#         # Apply wavelength segment to model grid
#         grid.customize(wave_rng=wr)
#
#         # Calculate the coefficient, mu, and radius grids
#         cg, mg, rg = ldcfit.ldc_grid(grid, profile)
#
#         # Add them to the arrays
#         coeff_grid[n] = cg
#         mu_grid[n] = mg
#         r_grid[n] = rg
#
#         del grid
#
#     # Draw plot
#     C = ['c{}'.format(i+1) for i in range(rows)]
#     fig = core.multiplot(rows, 1, xlabel='Wavelength', sharex=True,
#                          sharey=False,, ylabel=C, title=profile.title())
#
#     # Reshape the coeff grid so that the coefficients are
#     # the first dimension
#     coeff_grid = coeff_grid.reshape(rows,len(T),len(G),len(M),cols)
#
#     # For each coefficient, make a plot
#     for n,coeffs in enumerate(coeff_grid):
#         for nt,t in enumerate(T):
#             for ng,g in enumerate(G):
#                 for nm,m in enumerate(M):
#                     fig[n+1].plot(w_mean, coeff_grid[n,nt,ng,nm],
#                                   label=[t,g,m], marker='o')
#
#     # Plot a legend
#     fig[-1].legend(loc=0, frameon=False, fontsize=15)
#     diff = np.diff(w_mean)[0]/4.
#     plt.xlim(min(w_mean)-diff,max(w_mean)+diff)

"""Functions to simulate lightcurve data from ExoMAST parameters

Author: Joe Filippazzo
Email: jfilippazzo@stsci.edu
"""
import numpy as np
from bokeh.plotting import figure, show
try:
    import batman
except ImportError:
    print("Could not import batman. Functionality may be limited.")

from .. import utils


def simulate_lightcurve(target, radius=None, snr=1000., npts=1000, plot=False):
    """Simulate lightcurve data for the given target exoplanet

    Parameters
    ----------
    target: str
        The name of the target to simulate
    radius: array-like, float
        The radius or radii value(s) to use
    snr: float
        The signal to noise to use
    npts: int
        The number of points to plot
    nframes: int
        The number of frames to produce
    plot: bool
        Plot the figure

    Returns
    -------
    tuple
        The time, wavelengths, flux, and transit parameters
    """
    try:

        # Get the system parameters from ExoMAST
        targ, url = utils.get_target_data(target)
        name = targ.get('canonical_name') or target
        t0 = targ.get('transit_time', 0.)
        dt = targ.get('transit_duration', 1.)

        # Generate transit parameters with batman
        params = batman.TransitParams()
        params.t0 = t0
        params.rp = targ.get('Rp/Rs') or 0.1
        params.per = targ.get('orbital_period') or 0.5
        params.inc = targ.get('inclination') or 90.
        params.a = targ.get('a/Rs') or 15.
        params.ecc = targ.get('eccentricity') or 0.
        params.w = targ.get('omega') or 90.
        params.limb_dark = 'quadratic'
        params.transittype = 'primary'
        params.u = [0.1, 0.1]

        # Generate a time axis
        time = np.linspace(t0-dt, t0+dt, npts)

        # Make the transit model
        transit = batman.TransitModel(params, time, transittype='primary')

        # Generate the lightcurves
        flux = []
        if radius is None:
            radius = params.rp
        radii = [radius] if isinstance(radius, (int, float)) else radius
        for r in radii:
            params.rp = r
            flux.append(transit.light_curve(params))

        # Add noise
        ideal_flux = np.asarray(flux)
        flux = np.random.normal(loc=ideal_flux, scale=ideal_flux/snr)
        unc = flux - ideal_flux

        # Plot it
        if plot:
            fig = figure(title=name)
            fig.circle(time, flux[0])
            fig.xaxis.axis_label = targ.get('transit_time_unit') or 'MJD'
            fig.yaxis.axis_label = 'Relative Flux'
            show(fig)

        return time, flux, unc, params.__dict__

    except:
        raise ValueError('{}: Could not simulate light curve for this target'.format(target))

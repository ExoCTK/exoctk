"""
"""

import argparse

import corner
import numpy as np
from platon.fit_info import FitInfo
from platon.retriever import Retriever
from platon.constants import R_sun, R_jup, M_jup


class PlatonWrapper():
    """
    """

    def __init__(self, kwargs):
        """
        """

        self.kwargs = kwargs


    def retrieve(self):
        """
        """

        # Construct fit_info object
        retriever = Retriever()
        fit_info = retriever.get_default_fit_info(**self.kwargs)

        # Fit for the stellar radius and planetary mass using Gaussian priors.  This
        # is a way to account for the uncertainties in the published values
        fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
        fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)

        # Fit for other parameters using uniform priors
        R_guess = 1.4 * R_jup
        T_guess = 1200
        fit_info.add_uniform_fit_param('Rp', 0.9*R_guess, 1.1*R_guess)
        fit_info.add_uniform_fit_param('T', 0.5*T_guess, 1.5*T_guess)
        fit_info.add_uniform_fit_param("log_scatt_factor", 0, 1)
        fit_info.add_uniform_fit_param("logZ", -1, 3)
        fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 5)
        fit_info.add_uniform_fit_param("error_multiple", 0.5, 5)

        # Define bins, depths, and errors
        wavelengths = 1e-6*np.array([1.119, 1.138, 1.157, 1.175, 1.194, 1.213, 1.232, 1.251, 1.270, 1.288, 1.307, 1.326, 1.345, 1.364, 1.383, 1.401, 1.420, 1.439, 1.458, 1.477, 1.496, 1.515, 1.533, 1.552, 1.571, 1.590, 1.609, 1.628])
        bins = [[w-0.0095e-6, w+0.0095e-6] for w in wavelengths]
        depths = 1e-6 * np.array([14512.7, 14546.5, 14566.3, 14523.1, 14528.7, 14549.9, 14571.8, 14538.6, 14522.2, 14538.4, 14535.9, 14604.5, 14685.0, 14779.0, 14752.1, 14788.8, 14705.2, 14701.7, 14677.7, 14695.1, 14722.3, 14641.4, 14676.8, 14666.2, 14642.5, 14594.1, 14530.1, 14642.1])
        errors = 1e-6 * np.array([50.6, 35.5, 35.2, 34.6, 34.1, 33.7, 33.5, 33.6, 33.8, 33.7, 33.4, 33.4, 33.5, 33.9, 34.4, 34.5, 34.7, 35.0, 35.4, 35.9, 36.4, 36.6, 37.1, 37.8, 38.6, 39.2, 39.9, 40.8])

        # Run nested sampling
        result = retriever.run_multinest(bins, depths, errors, fit_info, plot_best=True)
        print(result)

        # Do some plotting
        fig = corner.corner(result.samples, weights=result.weights,
                            range=[0.99] * result.samples.shape[1],
                            labels=fit_info.fit_param_names)
        fig.savefig("emcee_corner.png")


if __name__ == '__main__':

    # Define various fitting parameters
    kwargs = {}
    kwargs['Rs'] = 1.19 * R_sun  # Required
    kwargs['Mp'] = 0.73 * M_jup  # Required
    kwargs['Rp'] = 1.4 * R_jup  # Required
    kwargs['T'] = 1200  # Required
    kwargs['logZ'] = 0  # Optional
    kwargs['CO_ratio'] = 0.53 # Optional
    kwargs['log_cloudtop_P'] = 4  # Optional
    kwargs['log_scatt_factor'] = 0  # Optional
    kwargs['scatt_slope'] = 4  # Optional
    kwargs['error_multiple'] = 1  # Optional
    kwargs['T_star'] = 6091  # Optional

    platon_wrapper = PlatonWrapper(kwargs)
    platon_wrapper.retrieve()

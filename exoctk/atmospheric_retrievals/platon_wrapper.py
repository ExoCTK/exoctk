"""A wrapper around the PLATON atmospheric retrieval tool.

This module serves as a wrapper around the atmospheric retrieval
software for ``platon``.  It provides methods for performing retreivals
through multinested sampling and EMCEE methods.  For more information
about ``platon``, please see ``https://platon.readthedocs.io``.

Authors
-------

    - Matthew Bourque

Use
---

    Users can perform the atmospheric retrieval by instantiating a
    ``PlatonWrapper`` object and passing fit parameters within the
    python environment.  An example of this is provided below
    ::

        import numpy as np
        from platon.constants import R_sun, R_jup, M_jup
        from exoctk.atmospheric_retrievals.platon_wrapper import PlatonWrapper

        # Build dictionary of parameters you wish to fit
        params = {
            'Rs': 1.19,  # Required
            'Mp': 0.73,  # Required
            'Rp': 1.4,  # Required
            'T': 1200,  # Required
            'logZ': 0,  # Optional
            'CO_ratio': 0.53,  # Optional
            'log_cloudtop_P': 4,  # Optional
            'log_scatt_factor': 0,  # Optional
            'scatt_slope': 4,  # Optional
            'error_multiple': 1,  # Optional
            'T_star': 6091}  # Optional

        # Initialize PlatonWrapper object and set the parameters
        pw = PlatonWrapper()
        pw.set_parameters(params)

        # Add any additional fit parameters
        R_guess = 1.4 * R_jup
        T_guess = 1200
        pw.fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
        pw.fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)
        pw.fit_info.add_uniform_fit_param('Rp', 0.9*R_guess, 1.1*R_guess)
        pw.fit_info.add_uniform_fit_param('T', 0.5*T_guess, 1.5*T_guess)
        pw.fit_info.add_uniform_fit_param("log_scatt_factor", 0, 1)
        pw.fit_info.add_uniform_fit_param("logZ", -1, 3)
        pw.fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 5)
        pw.fit_info.add_uniform_fit_param("error_multiple", 0.5, 5)

        # Define bins, depths, and errors
        pw.wavelengths = 1e-6*np.array([1.119, 1.138, 1.157, 1.175, 1.194, 1.213, 1.232, 1.251, 1.270, 1.288, 1.307, 1.326, 1.345, 1.364, 1.383, 1.401, 1.420, 1.439, 1.458, 1.477, 1.496, 1.515, 1.533, 1.552, 1.571, 1.590, 1.609, 1.628])
        pw.bins = [[w-0.0095e-6, w+0.0095e-6] for w in pw.wavelengths]
        pw.depths = 1e-6 * np.array([14512.7, 14546.5, 14566.3, 14523.1, 14528.7, 14549.9, 14571.8, 14538.6, 14522.2, 14538.4, 14535.9, 14604.5, 14685.0, 14779.0, 14752.1, 14788.8, 14705.2, 14701.7, 14677.7, 14695.1, 14722.3, 14641.4, 14676.8, 14666.2, 14642.5, 14594.1, 14530.1, 14642.1])
        pw.errors = 1e-6 * np.array([50.6, 35.5, 35.2, 34.6, 34.1, 33.7, 33.5, 33.6, 33.8, 33.7, 33.4, 33.4, 33.5, 33.9, 34.4, 34.5, 34.7, 35.0, 35.4, 35.9, 36.4, 36.6, 37.1, 37.8, 38.6, 39.2, 39.9, 40.8])

        # Perform the retrieval by your favorite method
        pw.retrieve_multinest()
        pw.retrieve_emcee()

        # Save the results to an output file
        pw.save_results()

        # Save a plot of the results
        pw.make_plot()

Dependencies
------------

    - ``corner``
    - ``numpy``
    - ``platon``
"""

import corner
import numpy as np
from platon.retriever import Retriever
from platon.constants import R_sun, R_jup, M_jup


def _apply_factors(params):
    """Apply appropriate multiplication factors to parameters.

    Parameters
    ----------
    params : dict
        A dictionary of parameters and their values for running the
        software.  See "Use" documentation for further details.
    """

    params['Rs'] = params['Rs'] * R_sun
    params['Mp'] = params['Mp'] * M_jup
    params['Rp'] = params['Rp'] * R_jup

    return params


def example(method):
    """Performs an example run of the emcee and multinest retrievals

    Parameters
    ----------
    method : str
        The method to use to perform the atmopsheric retrieval; can
        either be ``multinest`` or ``emcee``
    """

    # Ensure that the method parameter is valid
    assert method in ['multinest', 'emcee'], \
        'Unrecognized method: {}'.format(method)

    # Define the fit parameters
    params = {
        'Rs': 1.19,  # Required
        'Mp': 0.73,  # Required
        'Rp': 1.4,  # Required
        'T': 1200.0,  # Required
        'logZ': 0,  # Optional
        'CO_ratio': 0.53,  # Optional
        'log_cloudtop_P': 4,  # Optional
        'log_scatt_factor': 0,  # Optional
        'scatt_slope': 4,  # Optional
        'error_multiple': 1,  # Optional
        'T_star': 6091}  # Optional

    # Initialize the object and set the parameters
    pw = PlatonWrapper()
    pw.set_parameters(params)

    # Fit for the stellar radius and planetary mass using Gaussian priors.  This
    # is a way to account for the uncertainties in the published values
    pw.fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
    pw.fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)

    # Fit for other parameters using uniform priors
    R_guess = 1.4 * R_jup
    T_guess = 1200
    pw.fit_info.add_uniform_fit_param('Rp', 0.9*R_guess, 1.1*R_guess)
    pw.fit_info.add_uniform_fit_param('T', 0.5*T_guess, 1.5*T_guess)
    pw.fit_info.add_uniform_fit_param("log_scatt_factor", 0, 1)
    pw.fit_info.add_uniform_fit_param("logZ", -1, 3)
    pw.fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 5)
    pw.fit_info.add_uniform_fit_param("error_multiple", 0.5, 5)

    # Define bins, depths, and errors
    pw.wavelengths = 1e-6*np.array([1.119, 1.1387])
    pw.bins = [[w-0.0095e-6, w+0.0095e-6] for w in pw.wavelengths]
    pw.depths = 1e-6 * np.array([14512.7, 14546.5])
    pw.errors = 1e-6 * np.array([50.6, 35.5])

    # Do some retrievals
    if method == 'multinest':
        pw.retrieve_multinest()
        pw.save_results()
    elif method == 'emcee':
        pw.retrieve_emcee()

    # Make corner plot of results
    pw.make_plot()

    return pw


def _validate_parameters(supplied_params):
    """Ensure the supplied parameters are valid.  Throw assertion
    errors if they are not.

    Parameters
    ----------
    supplied_params : dict
        A dictionary of parameters and their values for running the
        software.  See "Use" documentation for further details.
    """

    # Define the parameters, their data types, and if they are required
    parameters = [('Rs', float, True),
                  ('Mp', float, True),
                  ('Rp', float, True),
                  ('T', float, True),
                  ('logZ', int, False),
                  ('CO_ratio', float, False),
                  ('log_cloudtop_P', int, False),
                  ('log_scatt_factor', int, False),
                  ('scatt_slope', int, False),
                  ('error_multiple', int, False),
                  ('T_star', int, False)]

    for parameter in parameters:
        name, data_type, required = parameter

        # Ensure that all required parameters are supplied
        if required:
            assert name in supplied_params, '{} missing from parameters'.format(parameter)

        # Ensure the supplied parameter is of a valid data type
        if name in supplied_params:
            assert type(supplied_params[name]) == data_type, '{} is not of type {}'.format(parameter, data_type)


class PlatonWrapper():
    """Class object for running the platon atmospheric retrieval
    software."""

    def __init__(self):
        """Initialize the class object."""

        self.retriever = Retriever()
        self.output_results = 'results.dat'
        self.output_plot = 'corner.png'

    def make_plot(self):
        """Create a corner plot that shows the results of the retrieval."""

        if self.method == 'emcee':
            fig = corner.corner(self.result.flatchain, range=[0.99] * self.result.flatchain.shape[1],
                        labels=self.fit_info.fit_param_names)

        elif self.method == 'multinest':
            fig = corner.corner(self.result.samples, weights=self.result.weights,
                                range=[0.99] * self.result.samples.shape[1],
                                labels=self.fit_info.fit_param_names)

        # Save the results
        self.output_plot = '{}_corner.png'.format(self.method)
        fig.savefig(self.output_plot)
        print('Corner plot saved to {}'.format(self.output_plot))

    def retrieve_emcee(self):
        """Perform the atmopsheric retrieval via emcee."""

        self.method = 'emcee'
        self.result = self.retriever.run_emcee(self.bins, self.depths, self.errors, self.fit_info)

    def retrieve_multinest(self):
        """Perform the atmopsheric retrieval via multinested sampling."""

        self.method = 'multinest'
        self.result = self.retriever.run_multinest(self.bins, self.depths, self.errors, self.fit_info, plot_best=False)

    def save_results(self):
        """Save the results of the retrieval to an output file."""

        # Save the results
        self.output_results = '{}_results.dat'.format(self.method)
        with open(self.output_results, 'w') as f:
            f.write(str(self.result))
        print('Results file saved to {}'.format(self.output_results))

    def set_parameters(self, params):
        """Set necessary parameters to perform the retrieval.

        Required parameters include ``Rs``, ``Mp``, ``Rp``, and ``T``.
        Optional parameters include ``logZ``, ``CO_ratio``,
        ``log_cloudtop_P``, ``log_scatt_factor``, ``scatt_slope``,
        ``error_multiple``, and ``T_star``.

        Parameters
        ----------
        params : dict
            A dictionary of parameters and their values for running the
            software.  See "Use" documentation for further details.
        """

        _validate_parameters(params)
        _apply_factors(params)
        self.params = params
        self.fit_info = self.retriever.get_default_fit_info(**self.params)


if __name__ == '__main__':

    example('emcee')
    example('multinest')

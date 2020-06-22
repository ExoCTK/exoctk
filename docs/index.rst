.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
    </style>

.. only:: html

   .. image:: ../exoctk/data/images/ExoCTK_logo.png

############################
ExoCTK Package Documentation
############################

ExoCTK is an open-source, modular data analysis package focused primarily on atmospheric characterization of exoplanets.

The subpackages currently included are:

Contamination and Visibility Calculators

Groups and Intrgrations Calculator

Transit Light-­Curve Fitting Tools

Limb-­Darkening Calculator

Phase Constraint Calculator
.. toctree::
  :maxdepth: 1

  tutorials/PhaseConstraintCalculator

Atmospheric Forward Modeling - Currently only available through the `website <https://exoctk.stsci.edu/fortney>`_.

All source code can be found on `GitHub <https://github.com/ExoCTK/exoctk>`_.

There is also a `website <https://exoctk.stsci.edu/>`_ where all current tools are available through interactive applications.


******************
User Documentation
******************


**Contamination and Visibility Calculator**

The Contamination Overlap tool is made up of two calculators: the visibility calculator, and the contamination calculator. The visibility calculator will determine at what position angles (PAs) your target is visible and when. The output results come in the form of interactive Bokeh plots, but users also have the option of downloading their data in an ascii file format. The contamination calculator will determine the percentage of spectral contamination that lands on your target at every PA. The visibility calculator is currently available for all JWST instruments, and the contamination calculator will be released for NIRISS (Mode: Single Object Slitless Spectroscopy), NIRCam (Mode: Grism Time Series), and MIRI (Mode: Low-Resolution Spectroscopy).

.. toctree::
  :maxdepth: 1

  source/exoctk.contam_visibility

**Groups and Integrations Calculator**

The groups_integrations tool is a JWST observation planning tool designed with exoplanet observations in mind. Given a potential observation (which requires transit time, and an estimate of model and magnitude for the host star, and specifics of instrument setup) it's simple to get an optimized groups and integrations plan for the observation. The example notebook also outlines cases for batch demos -- testing many transits/sources in a given instrument setup, or figuring out which instrument setup is best for a given transit.

The Groups and Integrations Calculator runs with pre-sampled pandeia data in the background -- so it can have the power of those carefully built instrument models, but still run 100 times faster.

.. toctree::
  :maxdepth: 1

  source/exoctk.groups_integrations

**Lightcurve Fitting Tool**

The light-curve fitting tool fits large numbers of spectroscopic light curves simultaneously while sharing model parameters across wavelengths and visits. It includes multiple uncertainty estimation algorithms and a comprehensive library of physical and systematic model components that are fully customizable.

.. toctree::
  :maxdepth: 1

  source/exoctk.lightcurve_fitting

**Limb Darkening Calculator**

The limb_darkening tool calculates limb-darkening coefficients for a specified stellar model, plotting results versus µ and wavelength. It uses high spectral resolution stellar atmospheric models, which are a necessity given JWST's expected precision.

.. toctree::
  :maxdepth: 1

  source/exoctk.limb_darkening

**Atmospheric Retrievals**

The ``atmospheric_retrievals`` subpackage within the ``exoctk`` package currently contains a module for performing retrievals via the `PLATON <https://platon.readthedocs.io/en/latest/>`_ package. `This Jupyter notebook <https://github.com/ExoCTK/exoctk/blob/master/exoctk/notebooks/atmospheric_retrievals_demo.ipynb>`_ contains a demo of how to use the `platon_wrapper <https://github.com/ExoCTK/exoctk/blob/master/exoctk/atmospheric_retrievals/platon_wrapper.py>`_ module.

Users who wish to use the ``atmospheric_retrievals`` tools may do so by installing the ``exoctk`` package.  Please see the `installation instructions <https://github.com/ExoCTK/exoctk#installation>`_ for further details.

.. toctree::
  :maxdepth: 1

  source/exoctk.atmospheric_retrievals

**Phase Constraint Calculator**

The Phase Constraint Calculator provides a simple interface for calculating JWST observation start windows in phase-space for both, transits and eclipse observations. This allows the user to quickly calculate 
minimum and maximum phase values that serve as inputs for the APT special requirements section when planning your observations.

.. toctree::
  :maxdepth: 1

  source/exoctk.phase_constraint_overlap

****************************************************
Installation Instructions and Notebook Availability
****************************************************

To install the ExoCTK package one can follow the instructions listed in our README available here on `GitHub <https://github.com/ExoCTK/exoctk#installation>`_.

There are also several Jupyter Notebooks available for users to aid in learning how to use the various tools included in the ExoCTK package which can be found within the 'Notebooks' folder following download of the ExoCTK repository or for viewing online `here <https://github.com/ExoCTK/exoctk/tree/master/exoctk/notebooks>`_.

***********************
Newsletter Subscription
***********************

If you'd like to stay up-to-date with our releases and updates we suggest subscribing to our newsletter. One can do so by following the instructions below:

**Subscribe by email - you are not require to log in to the ListServ interface.**

(1) Send an email from your mailbox to `exoctk-news-subscribe-request@maillist.stsci.edu`
(2) Subject and body should be blank. If you use an email signature, remove that as well and then send the message.
(3) You will receive a confirmation email - be sure to follow the instructions to ensure you are properly subscribed.

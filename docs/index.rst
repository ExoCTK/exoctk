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

Transit light-­curve fitting tools (TLC)

Limb-­darkening calculator (LDC)

The code can be found on `GitHub <https://github.com/ExoCTK/exoctk>`_ and there is also a  `website <https://exoctk.stsci.edu/>`_ the current tools are available through. 


******************
User Documentation
******************


**Contamination Visibility Tool**

.. toctree::
   :maxdepth: 1

   exoctk.contam_visibility

**Integrations Groups**

The groups_integrations tool is a JWST observation planning tool designed with exoplanet observations in mind. Given a potential observation (which requires transit time, and an estimate of model and magnitude for the host star, and specifics of instrument setup) it's simple to get an optimized groups and integrations plan for the observation. The example notebook also outlines cases for batch demos -- testing many transits/sources in a given instrument setup, or figuring out which instrument setup is best for a given transit.

The Groups and Integrations Calculator runs with pre-sampled pandeia data in the background -- so it can have the power of those carefully built instrument models, but still run 100 times faster.

.. toctree::
  :maxdepth: 1

  exoctk.integrations_groups

**Lightcurve Fitting Tool**

TLC fits large numbers of spectroscopic light curves simultaneously while sharing model parameters across wavelengths and visits. It includes multiple uncertainty estimation algorithms and a comprehensive library of physical and systematic model components that are fully customizable.

.. toctree::
  :maxdepth: 1

  exoctk.lightcurve_fitting

**Limb Darkening Calculator**

The limb_darkening tool calculates limb-darkening coefficients for a specified stellar model, plotting results versus µ and wavelength. It uses high spectral resolution stellar atmospheric models, which are a necessity given JWST's expected precision.

.. toctree::
   :maxdepth: 1

  exoctk.limb_darkening

**Nircam Coronagraphy**

.. toctree::
  :maxdepth: 1

  exoctk.nircam_coronagraphy
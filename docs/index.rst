.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
    </style>

.. only:: html

   .. image:: ../ExoCTK/data/images/ExoCTK_logo.png

############################
ExoCTK Package Documentation
############################

``ExoCTK`` is an open-source, modular data analysis package focused primarily on atmospheric characterization of exoplanets. 



Currently available subpackages:
-----------------------------------

Transit lightcurve fitting tools

Limb-­darkening calculator

Groups and Integrations Calculator

Contamination and Visibility Calculator

Atmoshperic Forward Modeling - Fortney Grid



Planned subpackages for future releases:
-----------------------------------

Atmoshperic Forward Modeling - Platon

Atmospheric Retrievals

Phase Constraint 



To sign up for update and new release emails see the newsletter subscription details at the bottom of the page. 



******************
User Documentation
******************


**Contamination and Visibility Calculator**

.. toctree::
   :maxdepth: 1

   source/ExoCTK.contam_visibility

**Groups and Integrations Calculator**

.. toctree::
  :maxdepth: 1

  source/ExoCTK.integrations_groups

**Lightcurve Fitting Tool**

TLC fits large numbers of spectroscopic light curves simultaneously while sharing model parameters across wavelengths and visits. It includes multiple uncertainty estimation algorithms and a comprehensive library of physical and systematic model components that are fully customizable.

.. toctree::
  :maxdepth: 1

  source/ExoCTK.lightcurve_fitting

**Limb Darkening Calculator**

LDC calculates limb-darkening coefficients for a specified stellar model, plotting results versus µ and wavelength. It uses high spectral resolution stellar atmospheric models, which are a necessity given JWST's expected precision.

.. toctree::
   :maxdepth: 1

  source/ExoCTK.limb_darkening
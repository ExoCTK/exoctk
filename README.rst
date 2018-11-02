.. image:: /exoctk/data/images/ExoCTK_logo.png
    :alt: ExoCTK Logo
    :scale: 10%

|build-status| |docs|

Introduction
------------
ExoCTK is an open-source, modular data analysis package focused primarily on atmospheric characterization of exoplanets. The subpackages included are:

* Transit light-­curve fitting tools
* Limb-­darkening calculator

Transit light-­curve fitting tools
---------------------------------
The ``lightcurve_fitting`` tool fits large numbers of spectroscopic light curves simultaneously while sharing model parameters across wavelengths and visits.  It includes multiple uncertainty estimation algorithms and a comprehensive library of physical and systematic model components that are fully customizable.

Limb Darkening Calculator
-------------------------
The ``limb_darkening`` tool calculates limb-darkening coefficients for a specified stellar model, plotting results versus µ and wavelength.  It uses high spectral resolution stellar atmospheric models, which are a necessity given JWST's expected precision.

.. figure:: /exoctk/data/images/limb_darkening.png
    :alt: LDC Demo
    :scale: 100%
    :align: center
    
    Coefficients of the quadratic and 4-parameter limb darkening profiles for the Phoenix ACES stellar atmosphere model [4000, 4.5, 0] through the WFC3_IR.G141 grism.



The Groups and Integrations Calculator
--------------------------------------
The ``groups_integrations`` tool is a JWST observation planning tool designed with
exoplanet observations in mind. Given a potential observation (which requires 
transit time, and an estimate of model and magnitude for the
host star, and specifics of instrument setup) it's simple to get an optimized
groups and integrations plan for the observation. The example notebook also
outlines cases for batch demos -- testing many transits/sources in a given instrument
setup, or figuring out which instrument setup is best for a given transit. 

The Groups and Integrations Calculator runs with pre-sampled `pandeia` data in
the background -- so it can have the power of those carefully built instrument
models, but still run 100 times faster. 


.. |build-status| image:: https://travis-ci.org/ExoCTK/exoctk.svg?branch=master
    :alt: build status
    :scale: 100%
    :target: https://travis-ci.org/ExoCTK/exoctk

.. |docs| image:: https://readthedocs.org/projects/docs/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: http://exoctk.readthedocs.io/en/latest/

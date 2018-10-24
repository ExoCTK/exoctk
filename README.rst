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
The `lightcurve_fitting` tool fits large numbers of spectroscopic light curves simultaneously while sharing model parameters across wavelengths and visits.  It includes multiple uncertainty estimation algorithms and a comprehensive library of physical and systematic model components that are fully customizable.

Limb Darkening Calculator
-------------------------
The `limb_darkening` tool calculates limb-darkening coefficients for a specified stellar model, plotting results versus µ and wavelength.  It uses high spectral resolution stellar atmospheric models, which are a necessity given JWST's expected precision.

.. figure:: /exoctk/data/images/LDC_demo.png
    :alt: LDC Demo
    :scale: 100%
    :align: center
    
    Limb darkening coefficients for the Phoenix ACES atmosphere models in the 1.5-1.7 micron range with (blue) and without (green) the 2MASS H band filter applied.

.. |build-status| image:: https://travis-ci.org/ExoCTK/exoctk.svg?branch=master
    :alt: build status
    :scale: 100%
    :target: https://travis-ci.org/ExoCTK/exoctk

.. |docs| image:: https://readthedocs.org/projects/docs/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: http://exoctk.readthedocs.io/en/latest/

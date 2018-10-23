.. image:: /exoctk/data/images/exoctk_logo.png
    :alt: ExoCTK Logo
    :scale: 50%

|build-status| |docs|

Introduction
------------
ExoCTK is an open-source, modular data analysis package focused primarily on atmospheric characterization of exoplanets. The subpackages included are:

* Transit light-­curve fitting tools (TLC)
* Limb-­darkening calculator (LDC)
* IFS exoplanet spectra extraction (IFS)
* Atmospheric forward models (AFM)
* Bayesian atmospheric retrieval (BAR)
* Planetary atmospheres libraries and tools (PAL)
* Transit Observation Tools (TOT)

Transit light-­curve fitting tools (TLC)
---------------------------------------
TLC fits large numbers of spectroscopic light curves simultaneously while sharing model parameters across wavelengths and visits.  It includes multiple uncertainty estimation algorithms and a comprehensive library of physical and systematic model components that are fully customizable.

Limb Darkening Calculator (LDC)
-------------------------------
LDC calculates limb-darkening coefficients for a specified stellar model, plotting results versus µ and wavelength.  It uses high spectral resolution stellar atmospheric models, which are a necessity given JWST's expected precision.

.. figure:: /exoctk/data/images/LDC_demo.png
    :alt: LDC Demo
    :scale: 100%
    :align: center
    
    Limb darkening coefficients for the Phoenix ACES atmosphere models in the 1.5-1.7 micron range with (blue) and without (green) the 2MASS H band filter applied.

IFS exoplanet spectra extraction (IFS)
--------------------------------------
IFS extracts the spectrum of a planet/brown-dwarf from GPI or JWST data that is compatible with retrieval codes. In particular uncertainties, along with their covariance, will be representative of the true statistic scatter in the data.

Atmospheric forward models (AFM)
--------------------------------
AFM provides a database of generic exoplanet atmospheric models similar to the ATLAS and Phoenix databases for stellar atmospheres.  The centrally located, well documented, and uniformly formatted grid of models will be used to plan exoplanet observations with JWST, WFIRST, and beyond.

Bayesian Atmospheric Retrievals (BAR)
-------------------------------------
BAR estimates parameter uncertainties, such as molecular abundances and thermal structures.

Planetary atmospheres libraries and tools (PAL)
-----------------------------------------------
PAL contains a robust set of molecular and atomic cross-section tables relevant to giant exoplanet atmospheres. It can also generate K-Coefficents on arbitrary wavelength grids and generate arbitrary exoplanet transmission spectra.

.. figure:: /exoctk/data/images/PAL_demo.png
    :alt: PAL Demo
    :scale: 100%
    :align: center

Transit Observation Tools (TOT)
-------------------------------
TOT is a transiting exoplanet noise simulator. The current implementation scales the measured flux, variance, and exposure time values from previously-observed systems, computes the expected rms per spectrophotometric channel, and estimates the transit/eclipse depth error based on the anticipated number of valid in- and out-of-transit data points. The uncertainty estimates depend on the orbital properties of the system, instrument configuration, and observation duration.

.. figure:: /exoctk/data/images/TOT_demo.png
    :alt: TOT Demo
    :scale: 100%
    :align: center

.. |build-status| image:: https://travis-ci.org/exoctk/exoctk.svg?branch=master
    :alt: build status
    :scale: 100%
    :target: https://travis-ci.org/exoctk/exoctk

.. |docs| image:: https://readthedocs.org/projects/docs/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: http://exoctk.readthedocs.io/en/latest/

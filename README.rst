.. image:: /exoctk/data/images/ExoCTK_logo.png
    :alt: ExoCTK Logo
    :scale: 10%

|build-status| |docs|


Introduction
------------
ExoCTK is an open-source, modular data analysis package focused primarily on atmospheric characterization of exoplanets. The subpackages included are:

* Transit Lightcurve Fitter
* Limb Darkening Calculator
* Groups and Integrations Calculator


Transit Lightcurve Fitter
-------------------------
The ``lightcurve_fitting`` tool fits large numbers of spectroscopic light curves simultaneously while sharing model parameters across wavelengths and visits.  It includes multiple uncertainty estimation algorithms and a comprehensive library of physical and systematic model components that are fully customizable.

.. figure:: /exoctk/data/images/lightcurve_fitting.png
    :alt: LCF Demo
    :scale: 100%
    :align: center


Limb Darkening Calculator
-------------------------
The ``limb_darkening`` tool calculates limb-darkening coefficients for a specified stellar model, plotting results versus µ and wavelength.  It uses high spectral resolution stellar atmospheric models, which are a necessity given JWST's expected precision.

.. figure:: /exoctk/data/images/limb_darkening.png
    :alt: LDC Demo
    :scale: 100%
    :align: center

    Coefficients of the quadratic and 4-parameter limb darkening profiles for the Phoenix ACES stellar atmosphere model [4000, 4.5, 0] through the WFC3_IR.G141 grism.


Groups and Integrations Calculator
----------------------------------
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


Installation
------------

The following are instructions on how to install the ``exoctk`` package for both users and contributors.  The ``exoctk`` repository provides a ``conda`` environment containing all of the dependencies needed to install and execute the ``exoctk`` software.

Download Anaconda or Miniconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You must first have a working installation of ``anaconda`` or ``miniconda`` for Python 3.  If you do not yet have this on your system, you can visit the following links for download and installation instructions:

- `Anaconda <https://www.anaconda.com/download/>`_
- `Miniconda <https://conda.io/en/latest/miniconda.html>`_

Package Installation
~~~~~~~~~~~~~~~~~~~~

The files needed for the environment installation are available within the `exoctk` package itself.  If you have not done so already, download and install the ``exoctk`` package, either by ``pip``:

::

  pip install exoctk


or by cloning the `exoctk` repository directly:

::

  git clone https://github.com/ExoCTK/exoctk.git
  cd exoctk
  python setup.py install


Environment Installation
~~~~~~~~~~~~~~~~~~~~~~~~
Following the installation of the ``exoctk`` package, you can install the ``conda`` environment via the ``env/environment-<PYTHON_VERSION>.yml`` files (relative to the parent directory of the repository).  Note that there is are separate environment files for each version of ``python`` that ``exoctk`` supports.  First, one should ensure that their version of ``conda`` is up to date:

::

  conda update conda


Next, one should activate the ``base`` environment:

::

  source activate base


Next, one can create the ``exoctk`` ``conda`` environment via the appropriate ``environment-<PYTHON_VERSION>.yml`` file:

::

  conda env create -f environment-<PYTHON_VERSION>.yml


where ``<PYTHON_VERSION>`` is the version of python you are using (e.g. ``environment-3.6.yml``)

Lastly, one can activate the newly-created environment with:

::

  source activate exoctk-<PYTHON_VERSION>
  
where again, ``<PYTHON_VERSION>`` is the version of python you are using (e.g. ``exoctk-3.6``)


Missing Dependencies?
~~~~~~~~~~~~~~~~~~~~~
If you find that the `exoctk` `conda` is missing a required dependency, please feel free to [submit a GitHub Issue](https://github.com/ExoCTK/exoctk/issues) detailing the problem.

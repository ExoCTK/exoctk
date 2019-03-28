.. image:: /exoctk/data/images/ExoCTK_logo.png
    :alt: ExoCTK Logo
    :scale: 10%

|build-status| |docs|


Introduction
------------
ExoCTK is an open-source, modular data analysis package focused primarily on atmospheric characterization of exoplanets. The subpackages included are:

* Contamination and Visibility Calculator
* Integrations and Groups Calculator
* Transit Light-Curve Fitter 
* Limb Darkening Calculator
* Atmospheric Forward Modeling - Currently only available through the _`website <https://exoctk.stsci.edu/fortney>`_. 

For more information on each package visit our documentation _`website <https://exoctk.readthedocs.io/en/latest/>`_. 

Most packages are also available through interactive tools at our _`website <https://exoctk.stsci.edu/>`_. 

Transit Light-Curve Fitter
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

Contamination Overlap
---------------------
The ``Contamination Overlap`` tool is made up of two calculators: the visibility
calculator, and the contamination calculator. The visibility calculator will
determine at what position angles (PAs) your target is visible and when. The
output results come in the form of an interactive Bokeh plot, but users also
have the option of downloading their data in an ascii file format. The
contamination calculator will determine the percentage of spectral contamination
that lands on your target at every PA. The visibility calculator is currently
available for all JWST instruments, and the contamination calculator will be
released for NIRISS (Mode: Single Object Slitless Spectroscopy), NIRCam
(Mode: Grism Time Series), and MIRI (Mode: Low-Resolution Spectroscopy).


.. figure:: /exoctk/data/images/visib_demo.png
:alt: VISdemo
:scale: 100%
:align: center

 The visibility is calculated for Kelt-8 with the NIRISS instrument. The
shaded region represents the PA range that a user can observe this target in.
The green line represents the nominal angle of the instrument for this target.

.. figure:: /exoctk/data/images/visib_table_demo.png
:alt: VISTdemo
:scale: 100%
:align: center

 Users also have the option to download their visibility data into an ascii
file for convenience. This is an example of an ascii file downloaded for the
Kelt-8 target using NIRISS. It lists the position angles (for the instrument
and JWST) with their corresponding dates.


Installation
------------

The following are instructions on how to install the ``exoctk`` package for both users and contributors.  The ``exoctk`` repository provides a ``conda`` environment containing all of the dependencies needed to install and execute the ``exoctk`` software.

Download Anaconda or Miniconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You must first have a working installation of ``anaconda`` or ``miniconda`` for Python 3.  If you do not yet have this on your system, you can visit the following links for download and installation instructions:

- `Anaconda <https://www.anaconda.com/download/>`_
- `Miniconda <https://conda.io/en/latest/miniconda.html>`_

Obtain the ``exoctk`` Package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To obtain the ``exoctk`` package with the necessary environment files, you can either install the package via ``pip``:

::

  pip install exoctk

or, clone the repository directly from GitHub:

::

  git clone https://github.com/ExoCTK/exoctk.git
  cd exoctk
  python setup.py [install|devlop]

Environment Installation
~~~~~~~~~~~~~~~~~~~~~~~~
You can install the ExoCTK ``conda`` environment via the ``env/environment-<PYTHON_VERSION>.yml`` files (relative to the parent directory of where the repository was installed).  Note that there are separate environment files for each version of ``python`` that ``exoctk`` supports.  First, one should ensure that their version of ``conda`` is up to date:

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

Package Installation
~~~~~~~~~~~~~~~~~~~~

In order to install the ``exoctk`` package within the newly-created ``conda`` environment, one must re-install the package, either via ``pip``:

::

  pip install exoctk


or by running the `exoctk` setup script:

::

  python setup.py [install|develop]



Missing Dependencies?
~~~~~~~~~~~~~~~~~~~~~
If you find that the `exoctk` `conda` is missing a required dependency, please feel free to `submit a GitHub Issue <https://github.com/ExoCTK/exoctk/issues>`_ detailing the problem.



Want to stay up-to-date with our releases and updates?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Subscribe to our newsletter by sending an email with a blank body and subject to `exoctk-news-subscribe-request@maillist.stsci.edu` from the email you want to enroll. You should then receive a confirmation email with instructions on how to confirm your subscription, please be sure to do so within 48 hours.

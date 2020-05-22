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
* Atmospheric Retrievals
* Phase Constraint Calculator
* Atmospheric Forward Modeling

For more information on each package visit our documentation at `readthedocs <https://exoctk.readthedocs.io/en/latest/>`_.

Most packages are also available through interactive tools at our `web portal <https://exoctk.stsci.edu/>`_.


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


Atmospheric Retrievals
----------------------

The ``atmospheric_retrievals`` subpackage within the ``exoctk`` package currently contains a module for performing retrievals via the `PLATON <https://platon.readthedocs.io/en/latest/>`_ package. `This Jupyter notebook <https://github.com/ExoCTK/exoctk/blob/master/exoctk/notebooks/atmospheric_retrievals_demo.ipynb>`_ contains a demo of how to use the `platon_wrapper <https://github.com/ExoCTK/exoctk/blob/master/exoctk/atmospheric_retrievals/platon_wrapper.py>`_ module.

Users who wish to use the ``atmospheric_retrievals`` tools may do so by installing the ``exoctk`` package.  Please see the `installation instructions <https://github.com/ExoCTK/exoctk#installation>`_ for further details.


Phase Constraint Calculator
-------------------------
The Phase Constraint Calculator provides a simple interface for calculating the JWST observation start window. The calculation currently only applies to transits, though one can subtract 0.5 from the phase values to compute the eclipse observation start window for planets on circular orbits. Enter the minimum and maximum phase values into the APT special requirements section when planning your observations.


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

To obtain the ``exoctk`` package with the necessary environment files, clone the repository directly from GitHub:

::

  git clone https://github.com/ExoCTK/exoctk.git
  cd exoctk


Environment Installation
~~~~~~~~~~~~~~~~~~~~~~~~
You can install the ExoCTK ``conda`` environment via the ``env/environment-<PYTHON_VERSION>.yml`` files (relative to the parent directory of where the repository was installed).  Note that there are separate environment files for each version of ``python`` that ``exoctk`` supports.  First, one should ensure that their version of ``conda`` is up to date:

::

  conda update conda


Next, one should activate the ``base`` environment:

::

  conda activate base


Next, one can create the ``exoctk`` ``conda`` environment via the appropriate ``environment-<PYTHON_VERSION>.yml`` file. One can find these files under the ``env`` directory and should run the following command in that directory:

::

  conda env create -f environment-<PYTHON_VERSION>.yml


where ``<PYTHON_VERSION>`` is the version of python you are using (e.g. ``environment-3.6.yml``)

Lastly, one can activate the newly-created environment with:

::

  conda activate exoctk-<PYTHON_VERSION>

where again, ``<PYTHON_VERSION>`` is the version of python you are using (e.g. ``exoctk-3.6``)


Package Installation
~~~~~~~~~~~~~~~~~~~~

In order to install the ``exoctk`` package within the newly-created ``conda``
environment, run the `exoctk` setup script:

::

  python setup.py [install|develop]


Obtain the ``exoctk`` Data
~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``exoctk`` data package will be available through the MAST portal soon!
Until then...

The suggested way to obtain the data is to execute the ``exoctk.utils.download_exoctk_data()`` function.  This function will download a series of compressed files from Box, extract the files, and organize them into a ``exoctk_data/`` directory.  Note that this can only be done once the ``exoctk`` package has been fully installed (see instructions above).

Lastly, export an environment variable for ``EXOCTK_DATA``.

- For Mac OS/Linux, add the line

::

    export EXOCTK_DATA='/path/to/your/unzipped/directory/exoctk_data/'

to your `.bashrc` or `.bash_profile`.

- For Windows, add an environment variable using System Utility.

Users may also download individual components of the ``exoctk`` data package directly through the `Box website <https://stsci.box.com/s/7ph64s6cfyusfcxjvih8ll5rn0ydzw86>`_.  Please note that materials must ultimately be placed within a ``exoctk_data/`` directory, and the ``EXOCTK_DATA`` environment variable be set in order for the ``exoctk`` package to work properly.


Missing Dependencies?
~~~~~~~~~~~~~~~~~~~~~
If you find that the `exoctk` `conda` is missing a required dependency, please feel free to `submit a GitHub Issue <https://github.com/ExoCTK/exoctk/issues>`_ detailing the problem.


Want to stay up-to-date with our releases and updates?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Subscribe to our newsletter by sending an email with a blank body and subject to ``exoctk-news-subscribe-request@maillist.stsci.edu`` from the email you want to enroll. You should then receive a confirmation email with instructions on how to confirm your subscription, please be sure to do so within 48 hours.

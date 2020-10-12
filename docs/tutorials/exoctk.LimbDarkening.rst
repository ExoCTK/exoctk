.. _LimbDarkening:

Calculate Limb Darkening Coefficients
=====================================

To calculate the limb darkening coefficients, we need a model grid.

In our first example, we use the Phoenix ACES models but any grid can be loaded into a modelgrid.ModelGrid() object if the spectra are stored as FITS files.

Two model grids are available in the EXOCTK_DATA directory and have corresponding child classes for convenience. The Phoenix ACES models and the Kurucz ATLAS9 models can be loaded with the modelgrid.ACES() and modelgrid.ATLAS9() classes respectively.

We can also use the resolution argument to resample the model spectra. This greatly speeds up the caluclations.

.. code-block:: python

	model_grid = modelgrid.ACES(resolution=700)
	print(model_grid.data)

``1056 models loaded from /Users/jfilippazzo/Documents/STScI/ExoCTK/EXOCTK_DATA/modelgrid/ACES/
 Teff  logg ...                          filename                         
------ ---- ... ----------------------------------------------------------
5800.0  3.0 ... lte05800-3.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
7600.0  5.0 ... lte07600-5.00+0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
4100.0  5.0 ... lte04100-5.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
6900.0  4.0 ... lte06900-4.00-0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
4400.0  3.0 ... lte04400-3.00+0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
2300.0  5.0 ... lte02300-5.00-0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
4200.0  4.0 ... lte04200-4.00-0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
4700.0  5.0 ... lte04700-5.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
6900.0  3.0 ... lte06900-3.00+0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
2600.0  4.0 ... lte02600-4.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
   ...  ... ...                                                        ...
5200.0  5.0 ... lte05200-5.00-0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
3500.0  3.0 ... lte03500-3.00+0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
3300.0  4.0 ... lte03300-4.00-0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
5100.0  4.0 ... lte05100-4.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
5400.0  5.0 ... lte05400-5.00-0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
3300.0  3.0 ... lte03300-3.00+0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
3500.0  4.0 ... lte03500-4.00-0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
3600.0  5.0 ... lte03600-5.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
5700.0  4.0 ... lte05700-4.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
6600.0  3.0 ... lte06600-3.00-0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
6000.0  4.0 ... lte06000-4.00+0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
Length = 1056 rows``

Now let's customize it to our desired effective temperature, surface gravity, metallicity, and wavelength ranges by running the customize() method on our grid.

.. code-block:: python

	model_grid.customize(Teff_rng=(2500,2600), logg_rng=(5,5.5), FeH_rng=(-0.5,0.5))

``12/1056 spectra in parameter range Teff:  (2500, 2600) , logg:  (5, 5.5) , FeH:  (-0.5, 0.5) , wavelength:  (<Quantity 0. um>, <Quantity 40. um>)
Loading flux into table...
100.00 percent complete!``

Now we can caluclate the limb darkening coefficients using the limb_darkening_fit.LDC() class.

.. code-block:: python 

	ld = lf.LDC(model_grid)


We just need to specify the desired effective temperature, surface gravity, metallicity, and the function(s) to fit to the limb darkening profile (including 'uniform', 'linear', 'quadratic', 'square-root', 'logarithmic', 'exponential', and 'nonlinear').

We can do this with for a single model on the grid:

..code-block:: python 

	teff, logg, FeH = 2500, 5, 0
	ld.calculate(teff, logg, FeH, 'quadratic', name='on-grid', color='blue')

``Bandpass trimmed to 0.05 um - 2.5999 um
1 bins of 100 pixels each.``

Or a single model off the grid, where the spectral intensity model is directly interpolated before the limb darkening coefficients are calculated. This takes a few seconds when calculating: 

..code-block:: python 

	teff, logg, FeH = 2511, 5.22, 0.04
	ld.calculate(teff, logg, FeH, 'quadratic', name='off-grid', color='red')

``
Interpolating grid point [2511/5.22/0.04]...
Run time in seconds:  11.166051149368286
Bandpass trimmed to 0.05 um - 2.5999 um
1 bins of 100 pixels each.``

Now we can see the results table using the following command:

..code-block:: python

	print(ld.results)

``  name    Teff  logg FeH   profile   filter ... color   c1    e1    c2    e2 
-------- ------ ---- ---- --------- ------- ... ----- ----- ----- ----- -----
 on-grid 2500.0  5.0  0.0 quadratic Top Hat ...  blue 0.218 0.024 0.391 0.033
off-grid 2511.0 5.22 0.04 quadratic Top Hat ...   red 0.224 0.025 0.398 0.033``

Using a Photometric Bandpass
----------------------------

Above we caluclated the limb darkening in a particular wavelength range set when we ran the ``customize()`` method on our ``core.ModelGrid()`` object.

Additionally, we can calculate the limb darkening through a particular photometric bandpass.

First we have to create a ``svo_filters.svo.Filter()`` object which we can then pass to the calculate method. Let's use 2MASS H-band for this example.

..code-block:: python 

	H_band = svo.Filter('2MASS.H')
	H_band.plot()


Now we can tell ``LDC.calculate()`` to apply the filter to the spectral intensity models before calculating the limb darkening coefficients using the bandpass argument. We'll compare the results of using the bandpass (purple line) to the results where we just used the wavelength window of 1.4-1.9 :math:`\mathcal $mu$ m` (green line).

..code-block:: python 

	ld = lf.LDC(model_grid)
	teff, logg, FeH = 2511, 5.22, 0.04
	ld.calculate(teff, logg, FeH, '4-parameter', name='Top Hat', color='green')
	ld.calculate(teff, logg, FeH, '4-parameter', bandpass=H_band, name='H band', color='purple')
	ld.plot(show=True)

``Interpolating grid point [2511/5.22/0.04]...
Run time in seconds:  12.76802396774292
Bandpass trimmed to 0.05 um - 2.5999 um
1 bins of 100 pixels each.
Interpolating grid point [2511/5.22/0.04]...
Run time in seconds:  12.711306095123291``

Using a Grism
-------------

Grisms are also supported. We can use the whole grism wavelength range (as if it was a bandpass) or truncate the grism to consider arbitrary wavelength ranges by setting the ``wave_min`` and ``wave_max`` arguments.

..code-block:: python 

	G141 = svo.Filter('WFC3_IR.G141', wave_min=1.61*q.um, wave_max=1.65*q.um)
	G141.plot()

``Bandpass trimmed to 1.11 um - 1.65 um
15 bins of 431 pixels each.``

Now we can caluclate the LDCs for each of the 15 wavelength bins of the G141 grism we just created, where the first column in the table is the bin center. This is not very useful to plot but... why not?

..code-block:: python 

	teff, logg, FeH = 2511, 5.22, 0.04
	ld.calculate(teff, logg, FeH, '4-parameter', bandpass=G141)
	print(ld.results)

``Interpolating grid point [2511/5.22/0.04]...
Run time in seconds:  12.591181993484497
  name   Teff  logg FeH    profile   ...   e2    c3     e3    c4     e4 
------- ------ ---- ---- ----------- ... ----- ------ ----- ------ -----
1.12 um 2511.0 5.22 0.04 4-parameter ... 0.011 -0.599 0.011  0.193 0.004
1.15 um 2511.0 5.22 0.04 4-parameter ... 0.016  0.454 0.017 -0.071 0.006
1.19 um 2511.0 5.22 0.04 4-parameter ...  0.01  0.458 0.011 -0.086 0.004
1.22 um 2511.0 5.22 0.04 4-parameter ...  0.01    0.7 0.011 -0.168 0.004
1.25 um 2511.0 5.22 0.04 4-parameter ... 0.008  0.321 0.009 -0.052 0.003
1.28 um 2511.0 5.22 0.04 4-parameter ... 0.019  0.832  0.02 -0.213 0.008
1.32 um 2511.0 5.22 0.04 4-parameter ... 0.006  0.766 0.007 -0.179 0.003
1.35 um 2511.0 5.22 0.04 4-parameter ... 0.024 -0.365 0.026  0.138  0.01
1.39 um 2511.0 5.22 0.04 4-parameter ... 0.048 -1.159 0.051  0.379 0.019
1.43 um 2511.0 5.22 0.04 4-parameter ... 0.012 -0.775 0.013  0.209 0.005
1.46 um 2511.0 5.22 0.04 4-parameter ... 0.033 -0.893 0.035  0.273 0.013
 1.5 um 2511.0 5.22 0.04 4-parameter ... 0.037 -0.776 0.039   0.26 0.015
1.54 um 2511.0 5.22 0.04 4-parameter ...  0.05 -0.623 0.053  0.235  0.02
1.59 um 2511.0 5.22 0.04 4-parameter ...  0.01  0.308 0.011 -0.049 0.004
1.63 um 2511.0 5.22 0.04 4-parameter ... 0.005   0.57 0.005 -0.131 0.002``





	
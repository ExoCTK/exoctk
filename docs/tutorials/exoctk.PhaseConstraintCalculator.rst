.. _PhaseConstraintCalculator:

Phase Constraint Calculator Tutorial
====================================

Due to uncertainties on the exact scheduling time of JWST observations (due, e.g., to previous observations, unforeseen scenarios, etc.), it is recommended that users of the observatory consider some flexibility in the start of their observations to give leeway of an hour to the observatory scheduling system (although observers can choose to narrow this window, there is a penalty in the charged time to the program). 

The time window can be defined in the time-domain in Astronomer's Proposal Tool (APT) under the `APT Special Requirements <https://jwst-docs.stsci.edu/jwst-astronomers-proposal-tool-overview/apt-workflow-articles/apt-special-requirements>`_, for periodic phenomena like transiting exoplanets this would be cumbersome to include as one would have to define a time-window for every possible transit/eclipse event on the current observing Cycle. Fortunately, APT also allows users to define this `window in phase-space <https://jwst-docs.stsci.edu/jppom/special-requirements/timing-special-requirements>`_ (i.e., in units of fractions of the orbital period), where the zero-phase can be arbitrarily defined. 

The **ExoCTK's phase-constraint package** was developed in order to perform the calculations of these windows in phase-space for any transiting exoplanet out there in a quick-and-easy way. Therefore this package greatly simplifies the work of observation planning when it comes to transiting exoplanet observations.

Quick Start
-----------
Let's suppose we want to obtain the 1-hour-window in phase space to schedule an observation of the primary transit of WASP-18b. To obtain this with **ExoCTK's phase-constraint** package, one would simply do:

.. code-block:: python

	import exoctk.phase_constraint_overlap.phase_constraint_overlap as pc
	min_phase, max_phase = pc.phase_overlap_constraint('WASP-18b', window_size = 1.)

Which would produce an output such as this: 

``Retrieved period is 0.94124. Retrieved t0 is 58374.669900000095.
Retrieved transit/eclipse duration is: 2.14368 hrs; implied pre mid-transit/eclipse on-target time: 2.89368 hrs.
Performing calculations with Period: 0.94124, t0: 58374.669900000095, ecc: None, omega: None degs, inc: None degs.
MINIMUM PHASE: 0.8276351762922669, MAXIMUM PHASE: 0.8719030215460457``

Two lines. That's all it took! In addition to the phase-constraints (the minimum and maximum phases), the phase_overlap_constraint call also returns the parameters it used to calculate the phase-constraint, along with some ephemerides of the planet, e.g., the period :math:`\mathcal P = 0.94124` days and time-of-transit center :math:`\mathcal {t}_{0} = 58374.6699` in Modified Julian Date (MJD, i.e., :math:`\mathcal JD - 2400000.5`). But how did this magic happen? What do these numbers actually mean? Keep reading to understand how the phase-constraint calculator actually works.

Primary Eclipses: Using the Phase-Constraint Calculator
-------------------------------------------------------
In the example above, the phase-constraint calculator returned the minimum and maximum phases for the exoplanet under study given only the planet name and the size of the window we were aiming to.

How Did It Do That?
~~~~~~~~~~~~~~~~~~~
In the background, the phase-constraint package automatically queries the exoplanet properties from exo.MAST given only the planet's name. Using this, it retrieves the properties of interest (period, :math:`\mathcal P`, and total transit duration, :math:`\mathcal {T}_{14}`, in this case) and, by default, assumes the observer wants to start the observations at the very least a time:

:math:`\mathcal {T}_{pre} = 0.75 + \textrm{MAX}(1, {T}_{14}/2) + {T}_{14}/2` hours

prior to mid-transit in this case. This time, by the way, is not arbitrary. Overall, the recommended (e.g., see this JWST observation planning step-by-step tutorial) time to spend on a target for a transit/eclipse observation is the above time :math:`\mathcal T_{pre}` prior to the mid-transit time, and :math:`\mathcal {T}_{post} = {T}_{14}/2 + MAX(1, {T}_{14}/2) + {T}_{W}` hours post mid-transit where :math:`\mathcal {T}_{W}` is the phase-constraint window (one hour in our example above). Using the retrieved properties for WASP-18b shown above, we can understand how the calculation was done in the background. The transit duration is :math:`\mathcal {T}_{14} = 2.14368` hours; the period is :math:`\mathcal P = 0.94124 = 22.58976` hours. The time :math:`\mathcal {T}_{pre}` is, thus, :math:`\mathcal {T}_{pre}\approx 2.89368`, which in phase-space units is

:math:`\mathcal {T}_{pre}/P \approx 0.128097`.

APT assumes the transit event is always located at phase 1 (or zero, whichever is more comfortable). Thus:

Maximum phase = :math:`\mathcal 1 - {T}_{pre}/P \approx 0.871903`,

which is exactly the maximum phase retrieved by the calculation. The minimum phase is simply one hour earlier in phase space. This gives:

Minimum phase = :math:`\mathcal 1 - ({T}_{pre}+1)/P \approx 0.827635`,

again, exactly the minimum phase quoted above.

Modifying Phase-Constraint Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The phase-constraint calculator allows to ingest a number of variables into the calculation in order to give control to the user in terms of the calculations they want to make. For instance, the pre-transit duration discussed above, :math:`\mathcal T_{pre}`, can be changed by the user. This is done using the ``pretransit_duration`` variable. Suppose we wanted to retrieve the phase-constraint that corresponds to a pre-transit duration of 4 hours instead. We can simply do:

.. code-block:: python

	minp, maxp = pc.phase_overlap_constraint('WASP-18b', window_size = 1., pretransit_duration = 4.)


``Retrieved period is 0.94124. Retrieved t0 is 58374.669900000095.
Performing calculations with Period: 0.94124, t0: 58374.669900000095, ecc: None, omega: None degs, inc: None degs.
MINIMUM PHASE: 0.7786607737311064, MAXIMUM PHASE: 0.8229286189848852``

Of course, that is not the only parameter we can change. In fact, every transit parameter of interest can be ingested to the phase_overlap_constraint function, in whose case the user-defined properties will override the exo.MAST ones. Let's use, for instance, the ephemerides found for WASP-18b by Shporer et al. (2019) - :math:`\mathcal P = 0.941452419`, :math:`\mathcal {t}_{0} = 2458002.354726`

.. code-block:: python

	minp, maxp = pc.phase_overlap_constraint('WASP-18b', window_size = 1., period = 0.941452419, t0 = 2458002.354726)

``Retrieved transit/eclipse duration is: 2.14368 hrs; implied pre mid-transit/eclipse on-target time: 2.89368 hrs.
Performing calculations with Period: 0.941452419, t0: 2458002.354726, ecc: None, omega: None degs, inc: None degs.
MINIMUM PHASE: 0.8276740668009621, MAXIMUM PHASE: 0.8719319239435721``

Note how they are only slightly differnt than the ones retrieved from exo.MAST! One important detail in the above calculation, is that the time-of-transit center is of no use in phase-space because, by definition, for APT this is at phase equals 1. This means one could put any place holder value for :math:`\mathcal {t}_{0}`, and the calculation would result in the exact same values:

.. code-block:: python 

	minp, maxp = pc.phase_overlap_constraint('WASP-18b', window_size = 1., period = 0.941452419, t0 = -1)

``Retrieved transit/eclipse duration is: 2.14368 hrs; implied pre mid-transit/eclipse on-target time: 2.89368 hrs.
Performing calculations with Period: 0.941452419, t0: -1, ecc: None, omega: None degs, inc: None degs.
MINIMUM PHASE: 0.8276740668009621, MAXIMUM PHASE: 0.8719319239435721``

Why does the phase-constraint overlap receives the time-of-transit center at all in the calculation? This will become clearer in the next section.

Secondary Eclipses: Using the Phase-Constraint Calculator
---------------------------------------------------------

Phase-Constraints for Secondary Eclipses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ExoCTK phase-constraint calculator can also obtain phase-constraints for secondary eclipses. This is indicated by the secondary flag in the ``phase_overlap_constraint`` function, which by default is ``False``. Setting it to ``True`` in the WASP-18b case gives:

.. code-block:: python

	minp, maxp = pc.phase_overlap_constraint('WASP-18b', window_size = 1., period = 0.941452419, secondary = True)

``Retrieved transit/eclipse duration is: 2.122865968966563 hrs; implied pre mid-transit/eclipse on-target time: 2.872865968966563 hrs.
Performing calculations with Period: 0.941452419, t0: None, ecc: 0.01, omega: 257.27 degs, inc: 85.68 degs.
MINIMUM PHASE: 0.3271883452721046, MAXIMUM PHASE: 0.3714462024147147``


Note that, given the small eccentricity and inclination of WASP-18b's orbit, in this case the maximum phase is almost equal to the value one would obtain assuming a circular orbit for this exoplanet, which would locate the maximum phase at :math:`\mathcal 0.5 - ({T}_{pre})/P \approx 0.3719` (i.e., with the secondary eclipse centered at phase :math:`\mathcal 0.5`). The difference is of seconds --- likely not critical for most JWST observations.

One important detail to remember before moving on: when ingesting the phase-constraints given above on APT, remember that we are still defining the zero-phase to be at the time of primary transit. This means that the phases given above only make sense to target eclipses in your observations if your "Zero Phase" in APT is set to the time of primary transit. This just makes it easier for the user: no need to compute times of secondary eclipses! (this is done in the background by the package). If you still want to know the time of secondary eclipse for some reason, keep reading. We got you covered!

Finding Secondary Eclipse Times
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To find the phase-constraints for secondary eclipses, in the background the ExoCTK phase-constraint package solves the proper minimization of the conjunction problem numerically (equation (5) in `Winn 2010 <https://arxiv.org/pdf/1001.2010v5.pdf>`_ ), and thus finds the time of secondary eclipse (in phase-space) to perform the calculation using the orbital elements retrieved from exo.MAST (for secondary eclipses, in addition to the period :math:`\mathcal P`, you need the inclination, :math:`\mathcal i`, the eccentricity, :math:`\mathcal e`, and the argument of periastron passage, :math:`\mathcal \omega` --- all of which can also be user-defined). This gives another functionality to the package: a secondary eclipse time calculator.

To retrieve the time of secondary eclipse, you can use the ``get_secondary_time`` flag in the ``phase_overlap_constraint`` function which, in addition to the minimum and maximum phases, returns the time of secondary eclipse just after the time of primary transit. Let's try this out for WASP-18b again:

.. code-block:: python

	minp, maxp, tsec = pc.phase_overlap_constraint('WASP-18b', window_size = 1., secondary = True, get_secondary_time = True)


``Retrieved period is 0.94124. Retrieved t0 is 58374.669900000095.
Retrieved transit/eclipse duration is: 2.122865968966563 hrs; implied pre mid-transit/eclipse on-target time: 2.872865968966563 hrs.
Performing calculations with Period: 0.94124, t0: 58374.669900000095, ecc: 0.01, omega: 257.27 degs, inc: 85.68 degs.
MINIMUM PHASE: 0.32714966265626544, MAXIMUM PHASE: 0.37141750791004413, TSEC: 58375.13919576395``

.. code-block:: python

	print('Secondary eclipse time:',tsec)

``Secondary eclipse time: 58375.13919576395``

As can be seen, the secondary eclipse time matches beautifully with our expectations for a non-eccentric orbit, which would give a secondary eclipse time of :math:`\mathcal {t}_{0} + P/2 \approx 58375.14052` MJD --- only a 5-second difference between the two results.

Exploring Challenges for Secondary-Eclipse Times: HD 80606b, GJ 436b and HAT-P-2b
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to showcase the power of the ExoCTK phase-constraint tool for secondary eclipse times and phase-constraints, we present here the results using our tool for more challenging systems in terms of predicting the location of their secondary eclipses. In order to compare with the literature values, however, we will be computing the phases at which secondary eclipses occur and not the times. This makes it easier to compare across datasets obtained at different epochs.

We start with HD 80606b, which is know to be very eccentric (:math:`\mathcal e=0.93`). A quick hack, if one is aiming at calculating the phase at which secondary eclipses occur is to let ``window_size = 0``. and ``pretransit_duration = 0`` (Of course, never input this in APT!). This will force the minimum and maximum phases to return the phase at which secondary eclipse occur (because one is forcing the window to be of zero width, and for the observations to start exactly at the time of secodary eclipse). Let's see how well our phase-constraint tool does in this challenguing system:

.. code-block:: python

	minp, maxp = pc.phase_overlap_constraint('HD80606 b', window_size = 0., pretransit_duration = 0., secondary = True)

``Retrieved period is 111.4367. Retrieved t0 is 55210.14280000003.
Performing calculations with Period: 111.4367, t0: 55210.14280000003, ecc: 0.93, omega: 301.03 degs, inc: 89.29 degs.
MINIMUM PHASE: 0.9455607255787186, MAXIMUM PHASE: 0.9455607255787186``

This matches pretty well with the phase at which secondary eclipse happens in the literature (0.947; `Laughlin, et al 2009 <https://www.nature.com/articles/nature07649>`_)! Note we are using more updated planetary parameters than the ones from Laughlin et al., 2009, which explains the slight discrepancy in phase-space.

Next, let's try GJ 436b --- a mildly eccentic system (:math:`\mathcal e = 0.138`):

.. code-block:: python

	minp, maxp = pc.phase_overlap_constraint('GJ 436b', window_size = 0., pretransit_duration = 0., secondary = True)

``Retrieved period is 2.64388312. Retrieved t0 is 54864.5839999998.
Performing calculations with Period: 2.64388312, t0: 54864.5839999998, ecc: 0.13827, omega: 351.0 degs, inc: 86.774 degs.
MINIMUM PHASE: 0.5868253469349103, MAXIMUM PHASE: 0.5868253469349103``

Woah! Excellent agreement with `Stevenson et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010Natur.464.1161S/abstract>`_, where the secondary eclipse phase is at :math:`\mathcal 0.5868 +/- 0.0003`. Finally, let's give the tool a shot with HAT-P-2b (:math:`\mathcal e=0.517`):

.. code-block:: python

	minp, maxp = pc.phase_overlap_constraint('HAT-P-2b', window_size = 0., pretransit_duration = 0., secondary = True)

``Retrieved period is 5.6335158. Retrieved t0 is 55288.349100000225.
Performing calculations with Period: 5.6335158, t0: 55288.349100000225, ecc: 0.5172, omega: 188.01 degs, inc: 86.16 degs.
MINIMUM PHASE: 0.1876234349401976, MAXIMUM PHASE: 0.1876234349401976``

Once again: beautiful agreement with de `Wit et al. (2017) <https://iopscience.iop.org/article/10.3847/2041-8213/836/2/L17/pdf>`_, where the secondary eclipse phase happens at 0.187.


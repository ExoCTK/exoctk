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

Of course, that is not the only parameter we can change. In fact, every transit parameter of interest can be ingested to the phase_overlap_constraint function, in whose case the user-defined properties will override the exo.MAST ones. Let's use, for instance, the ephemerides found for WASP-18b by Shporer et al. (2019) --- :math:`\mathcal P = 0.941452419`, :math:`\mathcal {t}_{0} = 2458002.354726`

.. code-block:: python

	minp, maxp = pc.phase_overlap_constraint('WASP-18b', window_size = 1., period = 0.941452419, t0 = 2458002.354726)

``Retrieved transit/eclipse duration is: 2.14368 hrs; implied pre mid-transit/eclipse on-target time: 2.89368 hrs.
Performing calculations with Period: 0.941452419, t0: 2458002.354726, ecc: None, omega: None degs, inc: None degs.
MINIMUM PHASE: 0.8276740668009621, MAXIMUM PHASE: 0.8719319239435721``

Note how they are only slightly differnt than the ones retrieved from exo.MAST! One important detail in the above calculation, is that the time-of-transit center is of no use in phase-space because, by definition, for APT this is at phase equals 1. This means one could put any place-holder value for :math: `\mathcal t0`, and the calculation would result in the exact same values:

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
words

Finding Secondary Eclipse Times
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
words

Exploring Challenges for Secondary-Eclipse Times: HD 80606b, GJ 436b and HAT-P-2b
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
words


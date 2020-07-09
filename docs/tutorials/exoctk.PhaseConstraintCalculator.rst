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

	Retrieved period is 0.94124. Retrieved t0 is 58374.669900000095.
	Retrieved transit/eclipse duration is: 2.14368 hrs; implied pre mid-transit/eclipse on-target time: 2.89368 hrs.
	Performing calculations with Period: 0.94124, t0: 58374.669900000095, ecc: None, omega: None degs, inc: None degs.
	MINIMUM PHASE: 0.8276351762922669, MAXIMUM PHASE: 0.8719030215460457

Primary Eclipses: Using the Phase-Constraint Calculator
-------------------------------------------------------
words

How Did It Do That?
~~~~~~~~~~~~~~~~~~~
words

Modifying Phase-Constraint Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
words

Secondary Eclipses: Using the Phase-Constraint Calculator
---------------------------------------------------------
words

Phase-Constraints for Secondary Eclipses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
words

Finding Secondary Eclipse Times
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
words

Exploring Challenges for Secondary-Eclipse Times: HD 80606b, GJ 436b and HAT-P-2b
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
words


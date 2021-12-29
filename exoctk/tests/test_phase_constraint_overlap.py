#! /usr/bin/env python

"""Tests for the ``phase_contraint_overlap`` module.

Authors
-------

    Matthew Bourque

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_phase_contraint_overlap.py
"""

from exoctk.phase_constraint_overlap import phase_constraint_overlap


def test_calculate_phase():
    """Tests the ``calculate_phase`` function"""

    minphase, maxphase = phase_constraint_overlap.calculate_phase(0.94124, 2.89368, 1)
    assert round(minphase, 3) == 0.828
    assert round(maxphase, 3) == 0.872


def test_calculate_pre_duration():
    """Tests the ``calculate_pre_duration`` function"""

    pretransit_duration = phase_constraint_overlap.calculate_pre_duration(3.5)
    assert pretransit_duration == 4.25


def test_calculate_tsec():
    """Tests the ``calculate_tsec`` function"""

    tsec = phase_constraint_overlap.calculate_tsec(0.94124, 0.01, 4.49020857, 1.4953981, 2458375.169883)
    assert round(tsec, 3) == 2458375.639


def test_drsky():
    """Tests the ``drsky`` function"""

    result = phase_constraint_overlap.drsky(0.0, 0.01, 4.49020857, 1.4953981)
    assert round(result, 3) == 0.216


def test_drsky_2prime():
    """Tests the ``drsky_2prime`` function"""

    result = phase_constraint_overlap.drsky_2prime(0.0, 0.01, 4.49020857, 1.4953981)
    assert round(result, 3) == -0.857


def test_drsky_prime():
    """Tests the ``drsky_prime`` function"""

    result = phase_constraint_overlap.drsky_prime(0.0, 0.01, 4.49020857, 1.4953981)
    assert round(result, 3) == -0.907


def test_getE():
    """Tests the ``getE`` function"""

    eccentric_anomaly = phase_constraint_overlap.getE(0.02, 0.01)
    assert round(eccentric_anomaly, 3) == 0.020


def test_getLTT():
    """Tests the ``getLTT`` function"""

    light_travel_time = phase_constraint_overlap.getLTT(1, 0.00200399, 0.01, 4.49020857, 1.4953981, 1.0)
    assert round(light_travel_time, 3) == -352.563


def test_getM():
    """Tests the ``getM`` function"""

    mean_anomaly = phase_constraint_overlap.getM(0.019801003143970763, 0.01)
    assert round(mean_anomaly, 3) == 0.020


def test_phase_overlap_constraint():
    """Tests the ``phase_overlap_constraint`` function"""

    minphase, maxphase = phase_constraint_overlap.phase_overlap_constraint('Wasp-18 b', 0.94124, None, 2.89368, 5.0, 1.0)
    assert round(minphase, 3) == 0.828
    assert round(maxphase, 3) == 0.872

#! /usr/bin/env python

"""Tests for the ``groups_integrations`` module.

Authors
-------

    Matthew Bourque

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_groups_integrations.py
"""

import os

from decimal import Decimal

from exoctk.groups_integrations import groups_integrations

INFILE = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'groups_integrations', 'groups_integrations_input_data.json')

TEST_DATA = {
    'ins': 'miri',
    'mag': Decimal('8.131'),
    'obs_time': Decimal('7.431039999999999'),
    'sat_max': Decimal('0.95'),
    'sat_mode': 'well',
    'time_unit': 'hour',
    'band': 'K',
    'mod': 'f5v',
    'filt': 'lrs',
    'subarray': 'slitlessprism',
    'filt_ta': 'f560w',
    'subarray_ta': 'slitlessprism',
    'n_group': 'optimize',
    'infile': INFILE}

EXPECTED_RESULTS = {
    'ins': 'miri',
    'mag': Decimal('8.131'),
    'obs_time': Decimal('7.431039999999999'),
    'sat_max': 183972.25,
    'sat_mode': 'well',
    'time_unit': 'hour',
    'band': 'K',
    'mod': 'f5v',
    'filt': 'lrs',
    'subarray': 'slitlessprism',
    'filt_ta': 'f560w',
    'subarray_ta': 'slitlessprism',
    'n_group': 36,
    'infile': INFILE,
    'duration_time': 7.432,
    'duration_time_ta_max': 0.63616,
    'duration_time_ta_min': 0.95424,
    'exposure_time': 7.225,
    'frames_per_group': 1,
    'frame_time': 0.159,
    'integration_time': 5.724,
    'max_saturation_prediction': 182002.902,
    'max_saturation_ta': 93529.062,
    'min_saturation_ta': 155881.77,
    'max_ta_groups': 3,
    'min_ta_groups': 5,
    'num_amps': 4,
    'num_columns': 72,
    'num_integrations': 4674,
    'num_reset_frames': 0,
    'num_rows': 416,
    'num_skips': 0,
    'observation_efficiency': 0.972,
    'ramp_time': 5.565,
    'ta_frame_time': 0.15904}


def test_calc_duration_time():
    """Tests the ``calc_duration_time`` function"""

    duration_time = groups_integrations.calc_duration_time(36, 4674, 0, 0.159, 1)
    assert round(duration_time, 2) == 26753.98


def test_calc_exposure_time():
    """Tests the ``calc_exposure_time`` function"""

    exposure_time = groups_integrations.calc_exposure_time(4674, 5.565)
    assert round(exposure_time, 2) == 26010.81


def test_calc_frame_time():
    """Tests the ``calc_frame_time`` function"""

    frame_time = groups_integrations.calc_frame_time(72, 416, 4, 'niriss')
    assert round(frame_time, 4) == 0.1251


def test_calc_groups_from_exp_time():
    """Tests the ``calc_groups_from_exp_time`` function"""

    groups = groups_integrations.calc_groups_from_exp_time(1000, 0.159)
    assert groups == 6289


def test_calc_integration_time():
    """Tests the ``calc_integration_time`` function"""

    integration_time = groups_integrations.calc_integration_time(36, 0.159, 1, 0)
    assert round(integration_time, 3) == 5.724


def test_calc_num_integrations():
    """Tests the ``calc_num_integrations`` function"""

    num_integrations = groups_integrations.calc_num_integrations(2.5, 36, 0, 0.159, 1)
    assert num_integrations == 1573


def test_calc_observation_efficiency():
    """Tests the ``calc_observation_efficiency`` function"""

    observation_efficiency = groups_integrations.calc_observation_efficiency(7.225, 7.432)
    assert round(observation_efficiency, 3) == 0.972


def test_calc_ramp_time():
    """Tests the ``calc_ramp_time`` function"""

    ramp_time = groups_integrations.calc_ramp_time(5.724, 0, 0.159)
    assert round(ramp_time, 3) == 5.565


def test_convert_saturation():
    """Tests the ``convert_saturation`` function"""

    max_saturation = groups_integrations.convert_saturation(100000, 'counts', 'miri', INFILE, True)
    assert round(max_saturation) == 193655


def test_interpolate_from_pandeia():
    """Tests the ``interpolate_from_pandeia`` function"""

    num_groups, max_saturation = groups_integrations.interpolate_from_pandeia(8.131, 'miri', 'lrs', 'slitlessprism', 'f5v', 'K', 0.159, 0.95, INFILE, False)
    assert num_groups == 1
    assert round(max_saturation, 3) == 31796.454


def test_map_to_ta_modes():
    """Tests the ``map_to_ta_modes`` function"""

    max_ta_groups, min_ta_groups = groups_integrations.map_to_ta_modes('miri', 10, 1)
    assert max_ta_groups == 9
    assert min_ta_groups == 3


def test_min_num_groups_for_sat():
    """Tests the ``min_num_groups_for_sat`` function"""

    minimum_num_groups_for_sat = groups_integrations.min_num_groups_for_sat(8.131, 'niriss', 'f480m', 'im', 'f5v', 'K', INFILE)
    assert minimum_num_groups_for_sat == 3


def test_perform_calculation():
    """Tests the ``perform_calculation`` function"""

    params = groups_integrations.perform_calculation(TEST_DATA)
    assert params == EXPECTED_RESULTS


def test_set_params_from_instrument():
    """Tests the ``set_params_from_instrument`` function"""

    rows, cols, amps, pixel_size, frame_time, num_reset_frames = groups_integrations.set_params_from_instrument('miri', 'slitlessprism')
    assert rows == 416
    assert cols == 72
    assert amps == 4
    assert pixel_size == 6.25e-06
    assert frame_time == 0.159
    assert num_reset_frames == 0


def test_set_frame_time():
    """Tests the ``set_frame_time`` function"""

    frame_time = groups_integrations.set_frame_time(INFILE, 'miri', 'slitlessprism', False)
    assert round(frame_time, 3) == 0.159

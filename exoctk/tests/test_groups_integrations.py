"""
"""

from decimal import Decimal

from exoctk.groups_integrations.groups_integrations import perform_calculation


def test_perform_calculation():
    """
    """

    infile = 'groups_integrations_1.4.json'

    params = {
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
        'infile': infile
    }

    expected_results = {
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
        'infile': infile,
        'n_col': 72,
        'n_row': 416,
        'n_amp': 4,
        'n_reset': 0,
        'n_frame': 1,
        'n_skip': 0,
        't_frame': 0.159,
        't_int': 5.725,
        't_ramp': 5.566,
        'n_int': 4673,
        't_exp': 7.225,
        't_duration': 7.432,
        'obs_eff': 0.972,
        'ta_t_frame': 0.15904,
        'min_ta_groups': 5,
        'max_ta_groups': 3,
        't_duration_ta_min': 0.95424,
        't_duration_ta_max': 0.63616,
        'max_sat_prediction': 182048.689,
        'max_sat_ta': 93529.062,
        'min_sat_ta': 155881.77
    }

    results = perform_calculation(params)

    # Ensure the number of keys in each dictionary are the same
    assert len(results) == len(expected_results)

    for item in results:

        print('Testing {}'.format(item))

        # Ensure the key names in each dictionary match
        assert item in expected_results

        # Ensure the key values in each dictionary match
        assert results[item] == expected_results[item]

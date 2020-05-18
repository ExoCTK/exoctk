"""
"""

from decimal import Decimal
import json
import pytest

from exoctk.groups_integrations.groups_integrations import perform_calculation

GROUPS_INTEGRATIONS_FILE = 'groups_integrations_1.4.json'
INSTRUMENTS =  ['miri', 'nircam', 'niriss', 'nirspec']
MAGNITUDES = [4.5, 6.5, 8.5, 10.5, 12.5]
PHOENIX_MODEL_KEYS = [
    'f0i', 'f0v', 'f5i', 'f2v', 'f5v', 'f8v',
    'g0v', 'g0iii', 'g2v', 'g5v', 'g8v', 'g0i', 'g5iii', 'g5i',
    'a0i', 'a0v', 'a1v', 'a5i', 'a3v', 'a5v',
    'm0i', 'm0iii', 'm0v', 'm2i', 'm2v', 'm5v',
    'k0v', 'k0iii', 'k2v', 'k0i', 'k5v', 'k5iii', 'k7v', 'k5i'
]
SCI_MODES_DICT = {'miri': {'filter': ['lrs'],
                           'subarray': ['slitlessprism']},
                  'niriss': {'filter': ['soss'],
                             'subarray': ['substrip256', 'substrip96']},
                  'nircam': {'filter': ['f322w2', 'f444w', 'f277w', 'f356w'],
                             'subarray': ['full', 'subgrism256', 'subgrism128', 'subgrism64']},
                  'nirspec': {'filter': [('g140h', 'f070lp'), ('g140h', 'f100lp'), ('g140m', 'f070lp'),
                                         ('g140m', 'f100lp'), ('g235h', 'f170lp'), ('g235m', 'f170lp'),
                                         ('g395h', 'f290lp'), ('g395m', 'f290lp'), ('prism', 'clear')],
                              'subarray': ['sub2048', 'sub1024a', 'sub1024b', 'sub512']}
}
TA_MODES_DICT = {'miri': {'filter': ['f560w', 'f100w', 'f1500w'],
                          'subarray': ['slitlessprism']},
                 'niriss': {'filter': ['f480m'],
                            'subarray': ['subtasoss']},
                 'nircam': {'filter': ['f335m'],
                            'subarray': ['sub32tats']},
                 'nirspec': {'filter': ['f110w', 'f140x', 'clear'],
                             'subarray': ['full', 'sub32', 'sub2048']}
}

# Generate data of expected values
EXPECTED_VALUES = []
for mode in ['sci_sat', 'ta_sat', 'ta_snr']:

    if mode == 'sci_sat':
        params_dict = SCI_MODES_DICT
    elif mode in ['ta_sat', 'ta_snr'] :
        params_dict = TA_MODES_DICT

    for instrument in INSTRUMENTS:
        for filt in params_dict[instrument]['filter']:

            if mode == 'sci_sat' and instrument == 'nirspec':
                filter_key = '{}_{}'.format(filt[1], filt[0])
            else:
                filter_key = filt

            for subarray in params_dict[instrument]['subarray']:
                for model in PHOENIX_MODEL_KEYS:
                    EXPECTED_VALUES.append((mode, instrument, filter_key, subarray, model))

# Read in data that will be tested
with open(GROUPS_INTEGRATIONS_FILE, 'r') as f:
    DATA = json.load(f)


def test_perform_calculation():
    """
    """

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
        'infile': GROUPS_INTEGRATIONS_FILE
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
        'infile': GROUPS_INTEGRATIONS_FILE,
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
        # print(results[item])
        # print(expected_results[item])
        # print(results[item] == expected_results[item])

@pytest.mark.parametrize('mode, instrument, filt, subarray, model', EXPECTED_VALUES)
def test_sat_modes(mode, instrument, filt, subarray, model):
    """Ensures the ``sci_sat`` key data is present and properly
    formatted in the ``groups_integrations.json`` file.

    Parameters
    ----------
    instrument : str
        The instrument of interest
    """

    assert instrument in DATA[mode]
    assert filt in DATA[mode][instrument]
    assert subarray in DATA[mode][instrument][filt]
    assert model in DATA[mode][instrument][filt][subarray]
    assert len(DATA[mode][instrument][filt][subarray][model]) == len(MAGNITUDES)

def test_groups_integrations_file():
    """Ensures that the necessary keys are in the
    ``groups_integrations.json`` file"""

    # Ensure first-order keys are present
    expected_keys = ['sci_sat', 'mags', 'fullwell', 'frame_time', 'ta_sat', 'ta_snr']
    for key in expected_keys:
        assert key in DATA

    # For mags key
    for magnitude in MAGNITUDES:
        assert magnitude in DATA['mags']

    # For fullwell key
    for instrument in INSTRUMENTS:
        assert instrument in DATA['fullwell']

    # For frame_time key
    for instrument in INSTRUMENTS:
        assert instrument in DATA['frame_time']
        for subarray in SCI_MODES_DICT[instrument]['subarray']:
            assert subarray in DATA['frame_time'][instrument]
        for subarray in TA_MODES_DICT[instrument]['subarray']:
            assert subarray in DATA['frame_time'][instrument]['ta']

    # For ta_snr key
    for instrument in INSTRUMENTS:
        assert instrument in DATA['ta_snr']
        for filt in TA_MODES_DICT[instrument]['filter']:
            assert filt in DATA['ta_snr'][instrument]
            for subarray in TA_MODES_DICT[instrument]['subarray']:
                assert subarray in DATA['ta_snr'][instrument][filt]
                for phoenix_model_key in PHOENIX_MODEL_KEYS:
                    assert phoenix_model_key in DATA['ta_snr'][instrument][filt][subarray]
                    assert len(DATA['ta_snr'][instrument][filt][subarray][phoenix_model_key]) == len(MAGNITUDES)

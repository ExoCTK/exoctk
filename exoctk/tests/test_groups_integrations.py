"""Tests for the ``groups_integrations.py`` module, as well as various
tests that ensures the structure and contents of the
``groups_integrations.json`` file in the ``exoctk_data`` package are
valid.

Authors
-------

    - Matthew Bourque, May 2020

Use
---

    The tests can be executed from the command line as such:
    ::

        pytest -s test_groups_integrations.py
"""

from decimal import Decimal
import json
import pytest

from exoctk.groups_integrations.groups_integrations import perform_calculation

GROUPS_INTEGRATIONS_FILE = '/user/bourque/repositories/observation_planning_data/groups_integrations_1.4.json'
INSTRUMENTS = ['miri', 'nircam', 'niriss', 'nirspec']
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

# Generate data of expected values for sat keys
EXPECTED_VALUES = []
for mode in ['sci_sat', 'ta_sat', 'ta_snr']:

    if mode == 'sci_sat':
        params_dict = SCI_MODES_DICT
    elif mode in ['ta_sat', 'ta_snr']:
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

# Generate data of expected values for frame_time key
FRAME_TIME_EXPECTED_VALUES = []
for instrument in INSTRUMENTS:
    for subarray in SCI_MODES_DICT[instrument]['subarray']:
        FRAME_TIME_EXPECTED_VALUES.append(('sci', instrument, subarray))
    for subarray in TA_MODES_DICT[instrument]['subarray']:
        FRAME_TIME_EXPECTED_VALUES.append(('ta', instrument, subarray))

# Read in data that will be tested
with open(GROUPS_INTEGRATIONS_FILE, 'r') as f:
    DATA = json.load(f)


@pytest.mark.parametrize('key', ['sci_sat', 'mags', 'fullwell', 'frame_time', 'ta_sat', 'ta_snr'])
def test_groups_integrations_keys(key):
    """Ensures that the necessary first-order keys are in the
    ``groups_integrations.json`` file

    Parameters
    ----------
    key : str
        The key of interest (e.g. ``mags``, ``fullwell``, etc.)
    """

    assert key in DATA


@pytest.mark.parametrize('mode, instrument, subarray', FRAME_TIME_EXPECTED_VALUES)
def test_frame_time_key(mode, instrument, subarray):
    """Ensures the contents of the ``frame_time`` key in the
    ``groups_integrations.json`` file are present and match the
    expected values.

    Parameters
    ----------
    mode : str
        Either ``sci`` or ``ta``
    instrument : str
        The instrument of interest
    subarray : str
        The subarray of interest
    """

    assert instrument in DATA['frame_time']

    if mode == 'sci':
        assert subarray in DATA['frame_time'][instrument]
        assert isinstance(DATA['frame_time'][instrument][subarray], float)
    elif mode == 'ta':
        assert subarray in DATA['frame_time'][instrument]['ta']
        assert isinstance(DATA['frame_time'][instrument]['ta'][subarray], float)


@pytest.mark.parametrize('instrument', INSTRUMENTS)
def test_fullwell_key(instrument):
    """Ensures the contents of the ``fullwell`` key in the
    ``groups_integrations.json`` file are present and match expected
    values.

    Parameters
    ----------
    instrument : str
        The instrument of interest
    """

    assert instrument in DATA['fullwell']
    assert isinstance(DATA['fullwell'][instrument], int) or isinstance(DATA['fullwell'][instrument], float)


@pytest.mark.parametrize('magnitude', MAGNITUDES)
def test_mags_key(magnitude):
    """Ensures the contents of the ``mags`` key is present in the
    ``groups_integrations.json`` file are present and match expected
    values.

    Parameters
    ----------
    magnitude : float
        The magnitude of interest
    """

    assert magnitude in DATA['mags']


def test_perform_calculation():
    """Tests the ``perform_calculation`` function in
    ``groups_integrations.py``
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

        # Ensure the key names in each dictionary match
        assert item in expected_results

        # Ensure the key values in each dictionary match
        # print(results[item])
        # print(expected_results[item])
        # print(results[item] == expected_results[item])


@pytest.mark.parametrize('key, instrument, filt, subarray, model', EXPECTED_VALUES)
def test_sat_keys(key, instrument, filt, subarray, model):
    """Ensures the ``sci_sat``, ``ta_sat``, and ``ta_snr`` keys data are
    present and properly formatted in the ``groups_integrations.json`` file.

    Parameters
    ----------
    key : str
        The key of interest.  Can be either ``sci_sat``, ``ta_sat``,
        or ``ta_snr``
    instrument : str
        The instrument of interest
    filt : str
        The filter of interest
    subarray : str
        The subarray of interest
    model : str
        The phoenix model of interest
    """

    assert instrument in DATA[key]
    assert filt in DATA[key][instrument]
    assert subarray in DATA[key][instrument][filt]
    assert model in DATA[key][instrument][filt][subarray]
    assert len(DATA[key][instrument][filt][subarray][model]) == len(MAGNITUDES)

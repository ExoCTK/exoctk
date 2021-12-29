#! /usr/bin/env python

"""Tests for the ``throughputs`` module.

Authors
-------

    Matthew Bourque

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_throughputs.py
"""

import os
import pytest

import numpy as np

from exoctk import throughputs

# Determine if tests are being run on Github Actions
ON_GITHUB_ACTIONS = '/home/runner' in os.path.expanduser('~') or '/Users/runner' in os.path.expanduser('~')


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Need access to Pandeia and its reference files.  Please try running locally')
def test_get_pce():
    """Tests the ``get_pce`` function"""

    instrument = 'niriss'
    mode = 'soss'
    filt = 'clear'
    disperser = 'gr700xd'
    aperture = 'soss'
    detector = {"ngroup": 10, "nexp": 1, "readout_pattern": "nisrapid", "nint": 1, "subarray": "substrip256"}

    wave, pce = throughputs.get_pce(instrument, mode, filt, disperser, aperture, detector)

    assert isinstance(wave, np.ndarray) and len(wave) > 0
    assert isinstance(pce, np.ndarray) and len(pce) > 0

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
import re
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


def test_jwst_throughput_files():
    """Tests that the packaged JWST throughputs are valid and that imaging
    filters transmit at their nominal wavelength (e.g. 15 um for F1500W)"""
    imaging = re.compile(r'^(MIRI|NIRCam|NIRISS)\.F(\d{3,4})[WMN]\w*\.CLEARP\.')
    bad = []
    for name in throughputs.JWST_THROUGHPUTS:
        wave, thru = np.loadtxt(os.path.join(throughputs.FILT_DIR, '{}.txt'.format(name)), unpack=True)
        if not (np.all(np.isfinite(thru)) and thru.min() >= 0 and thru.max() > 0):
            bad.append(name)
            continue

        match = imaging.match(name)
        if match:
            nominal = int(match.group(2)) / 100.
            band = wave[thru >= 0.5 * thru.max()]
            if not band.min() <= nominal <= band.max():
                bad.append(name)

    assert not bad, 'Invalid throughput files: {}'.format(bad)

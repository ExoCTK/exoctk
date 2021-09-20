#! /usr/bin/env python

"""Tests for various modules within the ``contam_visibility``
subpackage.

Authors
-------

    Matthew Bourque
    Ben Falk

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_contam_visibility.py
"""

import os
import sys

import numpy as np
import pytest

from exoctk.contam_visibility import field_simulator
from exoctk.contam_visibility import miniTools
from exoctk.contam_visibility import resolve
from exoctk.contam_visibility import visibilityPA

# Determine if tests are being run on Github Actions
ON_GITHUB_ACTIONS = '/home/runner' in os.path.expanduser('~') or '/Users/runner' in os.path.expanduser('~')


def test_checkVisPA():
    """Tests the ``checkVisPA`` function in the ``visibilityPA`` module"""

    ra = '24.3544618'
    dec = '-45.6777937'
    pa_good, pa_bad, gd, figure = visibilityPA.checkVisPA(ra, dec)

    assert isinstance(pa_good, list) and len(pa_good) > 0
    assert isinstance(pa_bad, list) and len(pa_bad) > 0
    assert isinstance(gd, list) and len(gd) > 0
    assert isinstance(figure, object)


# Cant get this to work because having latex installation issues
# def test_contamVerify():
#     """Tests the ``contamVerify`` function in the ``miniTools`` module"""

#     ra = 24.3544618
#     dec = -45.6777937
#     instrument = 'NIRISS'
#     APAlist = [1, 25, 181, 205]

#     pdf = miniTools.contamVerify(ra, dec, instrument, APAlist, web=True)

#     assert os.path.exists(pdf)


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Need access to trace data FITS files.  Please try running locally')
def test_fieldSim():
    """Tests the ``fieldSim`` function in the ``field_simulator`` module"""

    ra = '04 25 29.0162'
    dec = '-30 36 01.603'
    instrument = 'NIRISS'

    sim_cube = field_simulator.fieldSim(ra, dec, instrument)

    assert isinstance(sim_cube, np.ndarray) and len(sim_cube) > 0


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Test takes a long time to complete.  Please try running locally')
def test_gtsFieldSim():
    """Tests the ``gtsFieldSim`` function in the ``field_simulator`` module"""

    ra = '04 25 29.0162'
    dec = '-30 36 01.603'
    filt = 'F322W2'

    sim_cube = field_simulator.gtsFieldSim(ra, dec, filt)

    assert isinstance(sim_cube, np.ndarray) and len(sim_cube) > 0


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Test takes a long time to complete.  Please try running locally')
def test_lrsFieldSim():
    """Tests the ``lrsFieldSim`` function in the ``field_simulator`` module"""

    ra = '04 25 29.0162'
    dec = '-30 36 01.603'

    sim_cube = field_simulator.lrsFieldSim(ra, dec)

    assert isinstance(sim_cube, np.ndarray) and len(sim_cube) > 0


def test_resolve_target():
    """Tests the ``resolve_target`` function in the ``resolve`` module"""

    ra, dec = resolve.resolve_target('Wasp-18 b')

    assert ra == 24.3544618
    assert dec == -45.6777937


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Need access to trace data FITS files.  Please try running locally')
def test_sossFieldSim():
    """Tests the ``sossFieldSim`` function in the ``field_simulator`` module"""

    ra = '04 25 29.0162'
    dec = '-30 36 01.603'

    sim_cube = field_simulator.sossFieldSim(ra, dec)

    assert isinstance(sim_cube, np.ndarray) and len(sim_cube) > 0


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Need access to trace data FITS files.  Please try running locally')
def test_traceLength():
    """Tests the ``traceLenght`` function in the ``miniTools`` module"""

    instruments = 'NIRCam', 'NIRISS', 'MIRI'
    expected_results = {}
    expected_results['NIRCam'] = (0, 1342)
    expected_results['NIRISS'] = (1693, 5935)
    expected_results['MIRI'] = (41, 419)

    for instrument in instruments:
        target_trace_start, target_trace_stop = miniTools.traceLength(instrument)
        assert expected_results[instrument] == (target_trace_start, target_trace_stop)


@pytest.mark.skipif(sys.version_info > (3, 9), reason='jwst_gtvt does not currently support python>=3.9.')
def test_using_gtvt():
    """Tests the ``using_gtvt`` function in the ``visibilityPA`` module"""

    instrument = 'NIRISS'

    # this ra/dec has bad PAs
    ra = "-66"
    dec = "44"
    paMin, paMax, gd, fig, table, grouped_badPAs = visibilityPA.using_gtvt(ra, dec, instrument, targetName="Target", output="bokeh")
    assert grouped_badPAs is not None

    # this ra/dec has 100% coverage (no bad PAs)
    ra = '88'
    dec = '-64'
    output = visibilityPA.using_gtvt(ra, dec, instrument, targetName="Target", output="bokeh")

    assert output is not None
    assert len(output) == 6

    paMin, paMax, gd, fig, table, grouped_badPAs = output

    assert len(grouped_badPAs) == 0

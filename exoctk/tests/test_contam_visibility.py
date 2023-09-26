#! /usr/bin/env python

"""Tests for various modules within the ``contam_visibility``
subpackage.

Authors
-------

    Joe Filippazzo
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
from exoctk.contam_visibility import resolve
from exoctk.contam_visibility import new_vis_plot

# Determine if tests are being run on Github Actions
ON_GITHUB_ACTIONS = '/home/runner' in os.path.expanduser('~') or '/Users/runner' in os.path.expanduser('~')


def test_new_vis_plot():
    """Tests the `new_vis_plot.py` module"""
    ra, dec = '24.3544618', '-45.6777937' # WASP-18
    table = new_vis_plot.get_exoplanet_positions(ra, dec)

    assert str(type(table)) == "<class 'pandas.core.frame.DataFrame'>"

    plt = new_vis_plot.build_visibility_plot('WASP-18b', 'NIRISS', ra, dec)

    assert str(type(plt)) == "<class 'bokeh.plotting._figure.figure'>"


@pytest.mark.skipif(ON_GITHUB_ACTIONS, reason='Need access to trace data FITS files.  Please try running locally')
def test_field_simulation():
    """Tests the ``field_simulation`` function in the ``field_simulator`` module"""

    ra = '04 25 29.0162'
    dec = '-30 36 01.603'
    instrument = 'NIS_SUBSTRIP256'

    targframe, starcube, results = field_simulator.field_simulation(ra, dec, instrument)

    assert isinstance(targframe, (np.ndarray, list)) and isinstance(starcube, (np.ndarray, list))


def test_resolve_target():
    """Tests the ``resolve_target`` function in the ``resolve`` module"""

    ra, dec = resolve.resolve_target('Wasp-18 b')

    assert ra == 24.3544618
    assert dec == -45.6777937

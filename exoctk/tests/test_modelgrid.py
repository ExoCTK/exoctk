#! /usr/bin/env python

"""Tests for the ``modelgrid`` module.

Authors
-------

    Joe Filippazzo

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s modelgrid.py
"""

from pkg_resources import resource_filename

import numpy as np

from exoctk import modelgrid as mg


def test_modelgrid_object():
    """Test to see that ModelGrid object can be created"""
    print('Testing ModelGrid object creation...')

    # Load model grid
    mgrid = mg.ModelGrid(resource_filename('exoctk', 'data/core/modelgrid/'))

    assert isinstance(mgrid, mg.ModelGrid)


def test_model_getter_on_grid():
    """Test to see that an on-grid model can be pulled from a ModelGrid
    object"""
    print('Testing on-grid model getter from ModelGrid object...')

    # Load model grid
    mgrid = mg.ModelGrid(resource_filename('exoctk', 'data/core/modelgrid/'))

    # Fetch model
    model = mgrid.get(4000, 4.5, 0)

    assert isinstance(model.get('flux'), np.ndarray)


def test_model_getter_off_grid():
    """Test to see that an off-grid model can be pulled from a ModelGrid
    object"""
    print('Testing off-grid model getter from ModelGrid object...')

    # Load model grid
    mgrid = mg.ModelGrid(resource_filename('exoctk', 'data/core/modelgrid/'))

    # Fetch model
    model = mgrid.get(4023, 4.1, -0.1)

    assert isinstance(model.get('flux'), np.ndarray)

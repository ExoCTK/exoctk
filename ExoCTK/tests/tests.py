# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This packages contains package tests.
"""
from .. import core
from ..ldc import ldcfits

def test_load_ModelGrid():
    model_grid = core.ModelGrid('/user/jfilippazzo/Models/ACES/default/')
    assert len(model_grid.data)>0
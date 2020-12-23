# !/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for creating and managing grids of model spectra
"""
from glob import glob
import os
from pkg_resources import resource_filename
import warnings

from svo_filters.svo import Filter


FILT_DIR = resource_filename('exoctk', 'data/throughputs/')
JWST_THROUGHPUTS = [os.path.basename(file).replace('.txt', '') for file in glob(FILT_DIR + '*')]


class Throughput(Filter):
    def __init__(self, name, **kwargs):
        """
        Initialize the Throughput object as a child class of svo_filters.svo.Filter

        Parameters
        ----------
        name: str
            The [instrument.element_name] of the filter or disperser
        """
        # Check if JWST instrument
        if name in JWST_THROUGHPUTS:

            # Change throughput directory
            super().__init__(name, filter_directory=FILT_DIR, **kwargs)

        # Otherwise, just use svo_filters default
        else:
            super().__init__(name, **kwargs)


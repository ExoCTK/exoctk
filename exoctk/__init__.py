# !/usr/bin/python

"""
The Exoplanet Characterization Tool Kit is a collection of packages
used to reduce and analyze observations of transiting exoplanets.
"""

import os

from . import modelgrid
from . import references
from . import utils
from . import contam_visibility
from . import groups_integrations
from . import limb_darkening
from . import phase_constraint_overlap
from . import lightcurve_fitting
from . import log_exoctk

try:
    setup_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'setup.py')
    with open(setup_file, 'r') as f:
        data = f.readlines()
    __version__ = [line for line in data if 'version' in line][0].strip().split("version='")[-1].split("'")[0]

except (FileNotFoundError, NotADirectoryError):
    print('Could not determine exoctk version')
    __version__ = 'null'
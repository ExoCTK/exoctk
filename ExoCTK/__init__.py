# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The Exoplanet Characterization Tool Kit is a collection of packages used to reduce and analyze observations of transiting exoplanets
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from .core import *

from . import ldc
from . import pal
from . import tot
# from . import bar

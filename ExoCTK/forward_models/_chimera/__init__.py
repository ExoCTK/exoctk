"""
Package to generate spectroscopic forward models
"""
from . import fm
from . import ctran
from . import thermo
try:
    from . import _tran_module
except ImportError:
    if not _ASTROPY_SETUP_:
        raise
"""
Package to generate spectroscopic forward models
"""
try:
    from . import exotransmit
except ImportError:
    if not _ASTROPY_SETUP_:
        raise
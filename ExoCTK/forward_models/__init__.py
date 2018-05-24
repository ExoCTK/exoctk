"""
Package to generate spectroscopic forward models
"""
try:
    from . import exotransmit
    from . import forward_models
except ImportError:
    if not _ASTROPY_SETUP_:
        raise
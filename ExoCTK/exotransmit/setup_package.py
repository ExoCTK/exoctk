import io
from os.path import join
import os.path
import shutil
import sys

from distutils.core import Extension
from astropy_helpers import setup_helpers
import numpy as np

from astropy.extern import six


EXOTRANSMIT = os.path.relpath(os.path.dirname(__file__))

def get_extensions():
    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'].append(join(EXOTRANSMIT, "include"))
    cfg['sources'].extend(join(EXOTRANSMIT, 'exotransmit.pyx'))
    cfg = dict((str(key), val) for key, val in six.iteritems(cfg))
    return [Extension(str('ExoCTK.exotransmit'), 
        sources=[join(EXOTRANSMIT, 'exotransmit.pyx')],
        include_dirs=[join(EXOTRANSMIT, "include"), join(EXOTRANSMIT, np.get_include())])]
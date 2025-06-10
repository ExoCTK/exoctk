import numpy as np
from setuptools import setup
import os
import re

def get_version():
    with open(os.path.join("exoctk", "_version.py")) as f:
        match = re.search(r'__version__ = "([^"]+)"', f.read())
        return match.group(1)

setup(
    name='exoctk',
    include_dirs=[np.get_include()],
    version=get_version(),
)

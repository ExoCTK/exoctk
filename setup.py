import os
import re
from setuptools import setup
import numpy as np


def get_version():
    here = os.path.abspath(os.path.dirname(__file__))
    version_path = os.path.join(here, "exoctk", "_version.py")

    with open(version_path, encoding="utf-8") as f:
        content = f.read()

    match = re.search(r'^__version__ = ["\']([^"\']+)["\']', content, re.M)
    if not match:
        raise RuntimeError("Unable to find __version__ string in _version.py")

    return match.group(1)


setup(
    name="exoctk",
    version=get_version(),
    include_dirs=[np.get_include()],
    # You likely have more fields to include from pyproject.toml
)
